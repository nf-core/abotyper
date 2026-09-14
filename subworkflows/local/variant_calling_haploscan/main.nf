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

include { HAPLOSCAN } from '../../../modules/local/haploscan/main'
include { CLAIR3    } from '../../../modules/nf-core/clair3/main'

workflow VARIANTS_QUANTIFICATION {
    take:
    ch_bam     // channel: [ val(meta),  path(bam)   ]
    ch_bai     // channel: [ val(meta),  path(bai)   ]
    ch_fasta   // channel: [ val(meta1), path(fasta) ]
    ch_fai     // channel: [ val(meta1), path(fai)   ]
    model_url  // val: URL to the Clair3 model archive

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
    ch_clair3_model_dir = DOWNLOAD_CLAIR3_MODEL.out.model_dir

    //
    // MODULE: CLAIR3
    //
    ch_clair3_input = ch_bam_bai
        .combine(ch_clair3_model_dir)
        .map { meta, bam, bai, model_dir -> [meta, bam, bai, "", "${model_dir}/", "ont"] } // Already a value chanel

    CLAIR3(
        ch_clair3_input,
        ch_fasta.map { meta1, fasta -> [meta1, fasta] }.first(), // convert to value chanel
        ch_fai.map { meta1, fai -> [meta1, fai] }.first(), // convert to value chanel
    )

    emit:
    // metrics    = HAPLOSCAN.out.tsv
    // haplotypes = HAPLOSCAN.out.haplotypes

    clair3_vcf        = CLAIR3.out.vcf
    clair3_tbi        = CLAIR3.out.tbi
    clair3_phased_vcf = CLAIR3.out.phased_vcf
    clair3_phased_tbi = CLAIR3.out.phased_tbi
    // clair3_gvcf       = CLAIR3.out.gvcf
    // clair3_gtbi       = CLAIR3.out.gtbi
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
