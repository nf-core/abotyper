/*
  SUBWORKFLOW: VARIANTS_QUANTIFICATION (haploscan + clair3)

  Uses pysam_haploscan.py to compute per-position allele frequencies AND
  per-read haplotypes in a single pass over the aligned BAM -- no pileup file.

  Single-reference mode: every sample is aligned against the same combined
  ABO reference, so no per-exon metadata matching is required -- the FASTA
  and FAI are broadcast as value channels to every sample.

  Clair3 runs alongside HAPLOSCAN on the same BAM for independent small-
  variant calling with phasing.
*/

include { HAPLOSCAN         } from '../../../modules/local/haploscan/main'
include { CLAIR3            } from '../../../modules/nf-core/clair3/main'
include { ABO_PANEL2BED     } from '../../../modules/local/abo/panel2bed/main'
include { ABO_CLAIR2METRICS } from '../../../modules/local/abo/clair2metrics/main'
include { BCFTOOLS_NORM     } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW     } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_QUERY    } from '../../../modules/nf-core/bcftools/query/main'

workflow VARIANTS_QUANTIFICATION {
    take:
    ch_bam     // channel: [ val(meta),  path(bam)   ]
    ch_bai     // channel: [ val(meta),  path(bai)   ]
    ch_fasta   // channel: [ val(meta1), path(fasta) ]
    ch_fai     // channel: [ val(meta1), path(fai)   ]
    model_url  // val: URL to the Clair3 model archive
    panel_file // path: abo_variant_panel.yaml (single source of truth for diagnostic positions)

    main:

    // Join BAM and BAI (metadata matches exactly)
    ch_bam_bai = ch_bam.join(ch_bai, by: 0)

    // .combine always returns a queue channel; .first() makes it a reusable
    // value channel so it can be broadcast to every sample.
    ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .map { fasta_meta, fasta, fai_meta, fai -> [fasta_meta, fasta, fai] }
        .first()

    //
    // MODULE: HAPLOSCAN
    //
    // HAPLOSCAN(
    //     ch_bam_bai,
    //     ch_fasta_fai,
    // )

    //
    // PROCESS: DOWNLOAD_CLAIR3_MODEL
    //
    DOWNLOAD_CLAIR3_MODEL(model_url)
    // .first() converts this to a reusable value channel -- without it,
    // DOWNLOAD_CLAIR3_MODEL only runs once and its output queue channel is
    // exhausted after the first sample, so CLAIR3 would silently only fire
    // for sample #1 once there is more than one sample.
    ch_clair3_model_dir = DOWNLOAD_CLAIR3_MODEL.out.model_dir.first()

    //
    // MODULE: CLAIR3
    //
    ch_clair3_input = ch_bam_bai
        .combine(ch_clair3_model_dir)
        .map { meta, bam, bai, model_dir -> [meta, bam, bai, "", "${model_dir}/", "ont"] }

    CLAIR3(
        ch_clair3_input,
        ch_fasta.map { meta1, fasta -> [meta1, fasta] }.first(), // convert to value channel
        ch_fai.map { meta1, fai -> [meta1, fai] }.first(), // convert to value channel
    )

    //
    // MODULE: ABO_PANEL2BED
    //
    // Runs once against the panel YAML (reference-wide, not per-sample).
    // .first() for the same reason as ch_clair3_model_dir above: this only
    // ever emits once and must be broadcast to every sample.
    ABO_PANEL2BED(panel_file)
    ch_panel_bed = ABO_PANEL2BED.out.bed.first()

    //
    // MODULE: ABO_CLAIR2METRICS
    //
    // Converts Clair3's gVCF into HAPLOSCAN-format metrics, replacing
    // HAPLOSCAN as the ABO_GETABOSNPS/predict_abo_phenotype.py input.
    ABO_CLAIR2METRICS(
        CLAIR3.out.gvcf,
        panel_file,
    )

    //
    // MODULE: BCFTOOLS_NORM
    //
    // Left-align indels and split multiallelics on Clair3's own
    // PASS-filtered calls (not the gVCF -- symbolic <NON_REF> ref-blocks
    // aren't what bcftools norm/view/query are built for; the actual
    // phenotype-prediction metrics come from ABO_CLAIR2METRICS reading the
    // gVCF directly). This chain produces a panel-filtered, QC/export-ready
    // VCF and flat TSV alongside it, for review/cross-checking.
    BCFTOOLS_NORM(
        CLAIR3.out.vcf.join(CLAIR3.out.tbi, by: 0),
        ch_fasta.map { meta1, fasta -> [meta1, fasta] }.first(),
    )

    //
    // MODULE: BCFTOOLS_VIEW
    //
    // Restrict normalized calls to the ABO panel positions.
    BCFTOOLS_VIEW(
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.index, by: 0),
        ch_panel_bed,
        [],
        [],
    )

    //
    // MODULE: BCFTOOLS_QUERY
    //
    // Flatten the panel-filtered calls to a per-position TSV.
    BCFTOOLS_QUERY(
        BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.index, by: 0),
        [],
        [],
        [],
    )

    emit:
    metrics    = ABO_CLAIR2METRICS.out.tsv
    // haplotypes disabled until a Clair3 phased-VCF -> Haplotypes.tsv
    // adapter is built; PREDICTABOPHENOTYPE tolerates an empty channel
    // here (phase confidence just reports "No haplotype data available").
    haplotypes = Channel.empty()

    clair3_vcf        = CLAIR3.out.vcf
    clair3_tbi        = CLAIR3.out.tbi
    clair3_phased_vcf = CLAIR3.out.phased_vcf
    clair3_phased_tbi = CLAIR3.out.phased_tbi
    clair3_gvcf       = CLAIR3.out.gvcf
    clair3_gtbi       = CLAIR3.out.gtbi
    panel_bed         = ch_panel_bed

    panel_vcf       = BCFTOOLS_VIEW.out.vcf
    panel_vcf_index = BCFTOOLS_VIEW.out.index
    panel_query_tsv = BCFTOOLS_QUERY.out.output
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: DOWNLOAD_CLAIR3_MODEL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads Clair3 PyTorch model files (.pt) from specified URL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process DOWNLOAD_CLAIR3_MODEL {
    tag { url.tokenize('/').last() }

    input:
    val url

    output:
    path("clair3_model_dir"), emit: model_dir

    script:
    """
    mkdir -p clair3_model_dir

    # Download only the two PyTorch model files
    wget ${url}/pileup.pt -O clair3_model_dir/pileup.pt
    wget ${url}/full_alignment.pt -O clair3_model_dir/full_alignment.pt
    """
}
