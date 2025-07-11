/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_abotyper_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN_READS        } from '../subworkflows/local/minimap_align_exons/main'
include { PREDICTABOPHENOTYPE         } from '../subworkflows/local/predictabophenotype/main'
include { VARIANTS_QUANTIFICATION     } from '../subworkflows/local/variant_calling_mpileup/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ABOTYPER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fai            // channel: fasta from params.fai (ABO database)
    fasta          // channel: fasta from params.fasta (ABO database)
    exon6fai       // channel: fasta from params.exon6fai
    exon6fasta     // channel: fasta from params.exon6fasta
    exon7fai       // channel: fasta from params.exon7fai
    exon7fasta     // channel: fasta from params.exon7fasta
    logo           // channel: png from params.logo (custom pathwest logo)
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    // 
    // MODULE: FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Collect fastqc reports for multiqc
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    
    //
    // SUBWORKFLOW: minimap2/align || samtools/flagstat || samtools/stats
    MINIMAP2_ALIGN_READS (
        ch_samplesheet,
        exon6fai,
        exon6fasta,
        exon7fai,
        exon7fasta
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_READS.out.versions.first())
    
    //
    // SUBWORKFLOW: Run pileup for variants
    //
    VARIANTS_QUANTIFICATION (
        MINIMAP2_ALIGN_READS.out.exon6bai,
        MINIMAP2_ALIGN_READS.out.exon6bam,
        MINIMAP2_ALIGN_READS.out.exon7bai,
        MINIMAP2_ALIGN_READS.out.exon7bam,
        exon6fai,
        exon6fasta,
        exon7fai,
        exon7fasta
    )
    ch_versions = ch_versions.mix(VARIANTS_QUANTIFICATION.out.versions.first())
    
    //
    // SUBWORKFLOW: Run ABO prediction
    //
    PREDICTABOPHENOTYPE (
        VARIANTS_QUANTIFICATION.out.exon6metrics,
        MINIMAP2_ALIGN_READS.out.exon6cov,
        VARIANTS_QUANTIFICATION.out.exon7metrics,
        MINIMAP2_ALIGN_READS.out.exon7cov
    )
    ch_versions = ch_versions.mix(PREDICTABOPHENOTYPE.out.versions.first())
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
