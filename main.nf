#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/abotyper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/abotyper
    Website: https://nf-co.re/abotyper
    Slack  : https://nfcore.slack.com/channels/abotyper
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ABOTYPER                } from './workflows/abotyper'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_abotyper_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_abotyper_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_abotyper_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.abo_reference_fai   = getGenomeAttribute('abo_reference_fai')   ?: params.abo_reference_fai
params.abo_reference_fasta = getGenomeAttribute('abo_reference_fasta') ?: params.abo_reference_fasta
params.clair3_model_url    = getGenomeAttribute('clair3_model_url')    ?: params.clair3_model_url
params.abo_panel           = getGenomeAttribute('abo_panel')           ?: params.abo_panel

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFCORE_ABOTYPER {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    abo_reference_fai   = params.abo_reference_fai ? Channel.fromPath(params.abo_reference_fai)
        .map { it -> [[id: it.baseName], it] } : Channel.empty()
    abo_reference_fasta = params.abo_reference_fasta ? Channel.fromPath(params.abo_reference_fasta)
        .map { it -> [[id: it.baseName], it] } : Channel.empty()
    // Clair3 model: a URL is downloaded, else a local model folder is staged, else the packaged model is used.
    // The nf-core clair3 module has no packaged models under conda, so there the packaged model is downloaded by name.
    def no_packaged_models = !(workflow.containerEngine in ['docker', 'singularity', 'podman'])
    def clair3_url = params.clair3_model_url ?:
        (!params.clair3_model_path && no_packaged_models ? "${params.clair3_models_base_url.replaceAll('/+$', '')}/${params.clair3_model}/" : null)
    if (!params.clair3_model_url && clair3_url) {
        log.info "Clair3: packaged models are not available with this profile; downloading '${params.clair3_model}' from ${clair3_url}"
    }
    clair3_model = [
        url:      clair3_url,
        path:     (!clair3_url && params.clair3_model_path) ? Channel.fromPath(params.clair3_model_path, type: 'dir', checkIfExists: true).first() : null,
        packaged: params.clair3_model,
    ]
    abo_panel = params.abo_panel ? Channel.fromPath(params.abo_panel).collect() : Channel.empty()
    // staged as a task input so Clair3 can read it on any executor (cloud, HPC), not just locally
    abo_clair3_region_bed = params.abo_clair3_region_bed
        ? Channel.fromPath(params.abo_clair3_region_bed, checkIfExists: true).collect()
        : Channel.value([])

    ABOTYPER(
        samplesheet,
        abo_reference_fai,
        abo_reference_fasta,
        clair3_model,
        abo_panel,
        abo_clair3_region_bed,
    )

    emit:
    multiqc_report = ABOTYPER.out.multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_ABOTYPER(
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_ABOTYPER.out.multiqc_report
    )
}
