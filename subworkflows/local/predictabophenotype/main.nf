/*
  SUBWORKFLOW: PREDICTABOPHENOTYPE
*/

//  Description: Predicts ABO phenotype by combining variant frequency data
//   with BAM coverage, extracting SNPs, and mapping them to phenotypes.
//   Uses metadata-driven joining and structured per-sample output.


include { ABO_GETABOSNPS } from '../../../modules/local/abo/getabosnps/main'
include { ABO_SNPS2PHENO } from '../../../modules/local/abo/snps2pheno/main'

workflow PREDICTABOPHENOTYPE {
    take:
    ch_variants_freq // channel: [ val(meta), [ freq ] ]
    ch_bam_coverage  // channel: [ val(meta), [ cov ] ]
    ch_haplotypes    // channel: [ val(meta), path(*.Haplotypes.tsv) ] - from HAPLOSCAN

    main:

    // JOIN: Variant frequency with BAM coverage
    ch_combined_input = ch_variants_freq
        .join(ch_bam_coverage)

    /*
    MODULE: ABO_GETABOSNPS
    */
    ABO_GETABOSNPS(
        ch_combined_input
    )

    // PREP: Organize SNP reports AND Haplotypes.tsv by sample, single
    // "combined" folder per sample (single-reference mode -- one report
    // covers every exon section).
    ch_snp_reports = ABO_GETABOSNPS.out.phenotype
        .map { meta, file ->
            [meta.id, [file: file, type: 'phenotype']]
        }
        .mix(
            ch_haplotypes.map { meta, file ->
                [meta.id, [file: file, type: 'haplotype']]
            }
        )
        .groupTuple()
        .map { id, files ->
            def sample_dir = file("${params.outdir}/per_sample_processing/${id}")
            sample_dir.mkdirs()
            def combined_dir = sample_dir.resolve('combined')
            combined_dir.mkdirs()
            files.each {
                it.file.copyTo(combined_dir.resolve(it.file.name))
            }
            return sample_dir
        }
        .collect()

    // STAGE: Existing per_sample_processing directory in results
    ch_per_sample_processing = channel.fromPath("${params.outdir}/per_sample_processing", type: 'dir')

    /*
    MODULE: ABO_SNPS2PHENO
    */
    ABO_SNPS2PHENO(
        ch_snp_reports,
        ch_per_sample_processing,
    )

    emit:
    abo_results = ABO_SNPS2PHENO.out.txt
}
