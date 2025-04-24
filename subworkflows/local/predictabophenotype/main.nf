include { GETABOSNPS as SNPS_EXON6          } from '../../../modules/local/abo/abosnps/main'
include { GETABOSNPS as SNPS_EXON7          } from '../../../modules/local/abo/abosnps/main'
include { ABOSNPS2PHENO                     } from '../../../modules/local/abo/snps2pheno/main'

workflow PREDICTABOPHENOTYPE {

    take:
    ch_variants_freq_e6 // channel: [ val(meta), [ freq ] ]
    ch_bam_coverage_e6  // channel: [ val(meta), [ cov ] ]
    ch_variants_freq_e7 // channel: [ val(meta), [ freq ] ]
    ch_bam_coverage_e7  // channel: [ val(meta), [ cov ] ]

    main:

    ch_versions = Channel.empty()

    // Check channels for sanity
    // ch_variants_freq_e6.view { meta, freq -> "E6: meta=$meta, freq=$freq" }
    // ch_variants_freq_e7.view { meta, freq -> "E7: meta=$meta, freq=$freq" }


    SNPS_EXON6 ( ch_variants_freq_e6.join(ch_bam_coverage_e6), "exon6")
    SNPS_EXON7 ( ch_variants_freq_e7.join(ch_bam_coverage_e7), "exon7")

    ch_versions = ch_versions.mix(SNPS_EXON6.out.versions.first())
    ch_versions = ch_versions.mix(SNPS_EXON7.out.versions.first())

    // Just some shenanigans to keep the process waiting until all processes are completed.
    ch_SNP_reports = SNPS_EXON6.out.phenotype.map { meta, file -> 
            [meta.id, [exon: 'exon6', file: file]]
        }
        .mix(
            SNPS_EXON7.out.phenotype.map { meta, file -> 
                [meta.id, [exon: 'exon7', file: file]]
            }
        )
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

    // Stage the existing per_sample_processing directory
    ch_per_sample_processing = Channel.fromPath("${params.outdir}/per_sample_processing", type: 'dir')

    // One of these input channels needs to be removed later !!
    ABOSNPS2PHENO (
        ch_SNP_reports, 
        ch_per_sample_processing
    )
    ch_versions = ch_versions.mix(ABOSNPS2PHENO.out.versions.first())

    emit:
    versions = ch_versions    // channel: [ versions.yml ]
}
