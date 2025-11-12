/*
 * Subworkflow: PREDICTABOPHENOTYPE
 * Description: Predicts ABO phenotype by combining variant frequency data
 * with BAM coverage, extracting SNPs, and mapping them to phenotypes.
 * Uses metadata-driven joining and structured per-sample output.
 */

include { GETABOSNPS    } from '../../../modules/local/abo/abosnps/main'
include { ABOSNPS2PHENO } from '../../../modules/local/abo/snps2pheno/main'

workflow PREDICTABOPHENOTYPE {

    take:
    ch_variants_freq    // channel: [ val(meta), [ freq ] ] - with exon metadata
    ch_bam_coverage     // channel: [ val(meta), [ cov ] ] - with exon metadata

    main:

    ch_versions = Channel.empty()

    // JOIN: Variant frequency with BAM coverage
    ch_combined_input = ch_variants_freq
        .join(ch_bam_coverage)
        .map { meta, freq, cov ->
            [meta, freq, cov, meta.exon]  // Pass exon as parameter
        }

    /*
    MODULE: GETABOSNPS
    */
    GETABOSNPS (
        ch_combined_input
    )
    ch_versions = ch_versions.mix(GETABOSNPS.out.versions)

    // PREP:Organize SNP reports by sample and exon
    ch_SNP_reports = GETABOSNPS.out.phenotype
        .map { meta, file ->
            [meta.id, [exon: meta.exon, file: file]]
        }
        .groupTuple()
        .map { id, files ->
            def sample_dir = file("${params.outdir}/per_sample_processing/${id}")
            sample_dir.mkdirs()
            files.each {
                def exon_dir = sample_dir.resolve(it.exon)
                exon_dir.mkdirs()
                it.file.copyTo(exon_dir.resolve(it.file.name))
            }
            return sample_dir
        }
        .collect()

    // STAGE: Existing per_sample_processing directory in results
    ch_per_sample_processing = Channel.fromPath("${params.outdir}/per_sample_processing", type: 'dir')

    /*
    MODULE: ABOSNPS2PHENO
    */
    ABOSNPS2PHENO (
        ch_SNP_reports,
        ch_per_sample_processing
    )
    ch_versions = ch_versions.mix(ABOSNPS2PHENO.out.versions.first())

    emit:
    versions = ch_versions    // channel: [ versions.yml ]
}
