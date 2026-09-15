/*
  SUBWORKFLOW: VARIANTS_QUANTIFICATION
*/

//   Description: Calls small variants (with phasing) from each sample's aligned BAM using
//   Clair3 (--gvcf --enable_phasing), then converts its gVCF and phased VCF
//   into the per-position metrics / haplotype table PREDICTABOPHENOTYPE
//   expects (ABO_CLAIR2METRICS, ABO_CLAIR2HAPLOTYPES). Also regenerates the
//   ABO panel BED from the panel YAML once per run (ABO_PANEL2BED), broadcast
//   to every sample here and to the sibling VARIANTS_QC subworkflow.

include { CLAIR3                } from '../../../modules/nf-core/clair3/main'
include { ABO_PANEL2BED         } from '../../../modules/local/abo/panel2bed/main'
include { ABO_CLAIR2METRICS     } from '../../../modules/local/abo/clair2metrics/main'
include { ABO_CLAIR2HAPLOTYPES  } from '../../../modules/local/abo/clair2haplotypes/main'

workflow VARIANTS_QUANTIFICATION {
    take:
    ch_bam     // channel: [ val(meta),  path(bam)   ]
    ch_bai     // channel: [ val(meta),  path(bai)   ]
    ch_fasta   // channel: [ val(meta1), path(fasta) ]
    ch_fai     // channel: [ val(meta1), path(fai)   ]
    model_url  // val: URL to the Clair3 model archive
    panel_file // path: abo_variant_panel.yaml

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
    // PROCESS: DOWNLOAD_CLAIR3_MODEL
    //
    DOWNLOAD_CLAIR3_MODEL(model_url)
    ch_clair3_model_dir = DOWNLOAD_CLAIR3_MODEL.out.model_dir

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
    ABO_PANEL2BED(panel_file)

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
    // MODULE: ABO_CLAIR2HAPLOTYPES
    //
    // Converts Clair3's phased VCF into HAPLOSCAN-format Haplotypes.tsv,
    // driving PREDICTABOPHENOTYPE's PhaseConfidence scoring. Synthesized
    // from phase-set allele depths, not literal per-read observations --
    // see bin/clair2haplotypes.py docstring.
    ABO_CLAIR2HAPLOTYPES(
        CLAIR3.out.phased_vcf,
        panel_file,
    )

    emit:
    metrics    = ABO_CLAIR2METRICS.out.tsv
    haplotypes = ABO_CLAIR2HAPLOTYPES.out.tsv

    clair3_vcf        = CLAIR3.out.vcf
    clair3_tbi        = CLAIR3.out.tbi
    clair3_phased_vcf = CLAIR3.out.phased_vcf
    clair3_phased_tbi = CLAIR3.out.phased_tbi
    clair3_gvcf       = CLAIR3.out.gvcf
    clair3_gtbi       = CLAIR3.out.gtbi
    panel_bed         = ABO_PANEL2BED.out.bed
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: DOWNLOAD_CLAIR3_MODEL
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
