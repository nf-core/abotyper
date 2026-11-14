/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CHOPPER                 } from '../modules/nf-core/chopper/main'
include { MULTIQC                 } from '../modules/nf-core/multiqc'
include { NANOPLOT                } from '../modules/nf-core/nanoplot/main'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_abotyper_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ABOTYPERCORE_READTYPER  } from '../modules/local/abotypercore/readtyper/main'
include { MINIMAP2_ALIGN_READS    } from '../subworkflows/local/minimap2_align_reads'
include { PREDICTABOPHENOTYPE     } from '../subworkflows/local/predictabophenotype'
include { VARIANTS_QUANTIFICATION } from '../subworkflows/local/variants_quantification'
include { VARIANTS_QC             } from '../subworkflows/local/variants_qc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ABOTYPER {
    take:
    ch_samplesheet       // channel: samplesheet read in from --input
    abo_reference_fai    // channel: fai from params.abo_reference_fai
    abo_reference_fasta  // channel: fasta from params.abo_reference_fasta
    clair3_model         // map: url (download), path (channel of a local model folder) or packaged (model name)
    abo_panel            // channel: path from params.abo_panel, or [] for the panel bundled with abotyper-core
    clair3_region_bed    // channel: path from params.abo_clair3_region_bed, or [] for none

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: NANOPLOT
    //
    NANOPLOT (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.map { _meta, txt -> txt })

    //
    // MODULE: CHOPPER (--filter_reads): drop short and low-quality reads before alignment
    //
    ch_reads = ch_samplesheet
    if (params.filter_reads) {
        CHOPPER(ch_samplesheet, [])
        ch_reads = CHOPPER.out.fastq
    }

    //
    // SUBWORKFLOW: MINIMAP2_ALIGN_READS
    //
    MINIMAP2_ALIGN_READS(
        ch_reads,
        abo_reference_fasta,
        abo_reference_fai,
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        MINIMAP2_ALIGN_READS.out.stats.map { _meta, stats -> stats },
        MINIMAP2_ALIGN_READS.out.flagstat.map { _meta, flagstat -> flagstat },
    )

    ch_bam_bai = MINIMAP2_ALIGN_READS.out.bam.join(MINIMAP2_ALIGN_READS.out.bai)

    //
    // MODULE: ABOTYPERCORE_READTYPER - typing A, the primary call: every sample, from its reads
    //
    ABOTYPERCORE_READTYPER(
        ch_bam_bai,
        MINIMAP2_ALIGN_READS.out.fasta.first(),
    )

    //
    // Typing B, the per-position caller (Clair3 + BAM pileup), for the samples --clair3_mode selects:
    // always = every sample, review = samples the read-level typer did not PASS (REVIEW or NOCALL), off = none
    //
    if (params.clair3_mode == 'off') {
        ch_position_input = channel.empty()
    } else if (params.clair3_mode == 'review') {
        ch_position_input = ch_bam_bai
            .join(ABOTYPERCORE_READTYPER.out.tsv)
            .filter { _meta, _bam, _bai, tsv -> readTypingStatus(tsv) != 'PASS' }
            .map { meta, bam, bai, _tsv -> [meta, bam, bai] }
    } else {
        ch_position_input = ch_bam_bai
    }

    //
    // SUBWORKFLOW: VARIANTS_QUANTIFICATION
    //
    VARIANTS_QUANTIFICATION(
        ch_position_input.map { meta, bam, _bai -> [meta, bam] },
        ch_position_input.map { meta, _bam, bai -> [meta, bai] },
        MINIMAP2_ALIGN_READS.out.fasta,
        MINIMAP2_ALIGN_READS.out.fai,
        clair3_model,
        abo_panel,
        clair3_region_bed,
    )

    //
    // SUBWORKFLOW: VARIANTS_QC
    //
    VARIANTS_QC(
        VARIANTS_QUANTIFICATION.out.clair3_vcf,
        VARIANTS_QUANTIFICATION.out.clair3_tbi,
        MINIMAP2_ALIGN_READS.out.fasta,
        VARIANTS_QUANTIFICATION.out.panel_bed,
    )

    //
    // SUBWORKFLOW: PREDICTABOPHENOTYPE
    //

    PREDICTABOPHENOTYPE(
        VARIANTS_QUANTIFICATION.out.metrics,
        MINIMAP2_ALIGN_READS.out.coverage,
        VARIANTS_QUANTIFICATION.out.haplotypes,
        VARIANTS_QUANTIFICATION.out.evidence_status,
        VARIANTS_QUANTIFICATION.out.concordance,
        MINIMAP2_ALIGN_READS.out.stats,
        MINIMAP2_ALIGN_READS.out.flagstat,
        ABOTYPERCORE_READTYPER.out.tsv,
        ABOTYPERCORE_READTYPER.out.haplotypes,
        abo_panel,
    )

    // Per-sample evidence summary (Clair3 + BAM pileup) and read-level typing as MultiQC tables
    ch_multiqc_files = ch_multiqc_files.mix(
        PREDICTABOPHENOTYPE.out.evidence_summary,
        PREDICTABOPHENOTYPE.out.readtyping_mqc,
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'abotyper_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MULTIQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    def ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'abotyper'],
                files,
                params.multiqc_config
                    ? file(params.multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Status of a read-level call (PASS, REVIEW or NOCALL) from its one-row ReadTyping.tsv; NOCALL if unreadable
def readTypingStatus(tsv) {
    def lines = tsv.readLines()
    if (lines.size() < 2) {
        return 'NOCALL'
    }
    def header = lines[0].split('	') as List
    def values = lines[1].split('	', -1) as List
    def i = header.indexOf('status')
    return (i >= 0 && i < values.size()) ? values[i] : 'NOCALL'
}
