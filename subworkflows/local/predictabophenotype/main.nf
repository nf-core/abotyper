include { GETABOSNPS as SNPS_EXON6          } from '../../../modules/local/abo/abosnps/main'
include { GETABOSNPS as SNPS_EXON7          } from '../../../modules/local/abo/abosnps/main'
include { ABOSNPS2PHENO                     } from '../../../modules/local/abo/snps2pheno/main'

workflow PREDICTABOPHENOTYPE {

    take:
    ch_variants_freq_e6 // channel: [ val(meta), [ freq ] ]
    ch_variants_freq_e7 // channel: [ val(meta), [ freq ] ]

    main:

    ch_versions = Channel.empty()

    // Check channels for sanity
    // ch_variants_freq_e6.view { meta, freq -> "E6: meta=$meta, freq=$freq" }
    // ch_variants_freq_e7.view { meta, freq -> "E7: meta=$meta, freq=$freq" }

    SNPS_EXON6 ( ch_variants_freq_e6 )
    SNPS_EXON7 ( ch_variants_freq_e7 )

    ch_versions = ch_versions.mix(SNPS_EXON6.out.versions.first())
    ch_versions = ch_versions.mix(SNPS_EXON7.out.versions.first())

    ch_SNP_reports = Channel.empty()
    ch_SNP_reports = ch_SNP_reports.mix(
        SNPS_EXON6.out.phenotype.collect().ifEmpty([]),
        SNPS_EXON7.out.phenotype.collect().ifEmpty([])
        )

    // ABOSNPS2PHENO ('${params.outdir}/per_sample_processing', ch_SNP_reports)
    // ch_versions = ch_versions.mix(ABOSNPS2PHENO.out.versions.first())

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}

