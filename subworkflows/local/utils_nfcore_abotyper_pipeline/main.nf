//
// Subworkflow with functionality specific to the nf-core/abotyper pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/abotyper ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/abotyper/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Point conda users at the pinned environments
    //
    condaLockNote()

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1 ->
                // `control` is empty unless the samplesheet declares one
                def control = meta.control ?: ''
                // sample name and barcode reported in the result tables (see sampleNameAndBarcode)
                def names   = sampleNameAndBarcode(meta.id, meta.barcode)
                def extra   = [ control:control, sample_name:names.sample_name, barcode:names.barcode ]
                // Nanopore reads: one FASTQ per row, single-end
                return [ meta.id, meta + extra + [ single_end:true ], [ fastq_1 ] ]
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.isBound() ? multiqc_reports.getVal() : [],
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def condaLockNote() {
    def profiles = workflow.profile.tokenize(',')
    if (profiles.intersect(['conda', 'mamba']).size() > 0 && !profiles.contains('conda_lock')) {
        log.info("Conda environments are built from each module's environment.yml. For the exact tested " +
            "package builds, add conda_lock: -profile conda,conda_lock (https://nf-co.re/abotyper/usage#conda)")
    }
}

def validateInputParameters() {
    if (params.clair3_model_url && params.clair3_model_path) {
        error("Please give either --clair3_model_url (download a Clair3 model) or --clair3_model_path (a local model folder), not both.")
    }
}

//
// Sample name and barcode reported for a samplesheet row
//
// With a `barcode` column value, the sample is reported under its full `sample` name with that
// barcode. Without one, a `sample` name ending in _barcode<digits> is split into name and barcode
// (S1_barcode07 -> S1, barcode07); any other name is reported whole, and the barcode is left empty
// (reported as barcode00). The typing step takes both from here, not from folder names.
//
def sampleNameAndBarcode(String id, barcode) {
    if (barcode) {
        return [ sample_name: id, barcode: normaliseBarcode(barcode.toString()) ]
    }
    def m = id =~ /(?i)^(.+)_barcode(\d+)$/
    if (m.matches()) {
        return [ sample_name: m.group(1), barcode: normaliseBarcode(m.group(2)) ]
    }
    return [ sample_name: id, barcode: '' ]
}

// 7, 07 or barcode7 -> barcode07 (two digits at least, as the barcode suffix of a sample name)
def normaliseBarcode(String value) {
    def num = value.replaceFirst(/(?i)^barcode/, '')
    return (num as Integer) <= 99 ? "barcode${num.padLeft(2, '0')}" : "barcode${num}"
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Rows of the same sample are concatenated, so they must agree on whether it is a control
    def control_ok = metas.collect{ meta -> meta.control }.unique().size == 1
    if (!control_ok) {
        error("Please check input samplesheet -> rows of the same sample must have the same 'control' value: ${metas[0].id}")
    }

    // ... and on its barcode
    def barcode_ok = metas.collect{ meta -> meta.barcode }.unique().size == 1
    if (!barcode_ok) {
        error("Please check input samplesheet -> rows of the same sample must have the same 'barcode' value: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "NanoPlot (De Coster et al. 2018),",
            "minimap2 (Li 2018),",
            "SAMtools and BCFtools (Danecek et al. 2021),",
            "Clair3 (Zheng et al. 2022),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>De Coster W, D\'Hert S, Schultz DT, Cruts M, Van Broeckhoven C (2018) NanoPack: visualizing and processing long-read sequencing data, Bioinformatics, 34(15), 2666-2669. doi: 10.1093/bioinformatics/bty149</li>",
            "<li>Li H (2018) Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, 34(18), 3094-3100. doi: 10.1093/bioinformatics/bty191</li>",
            "<li>Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H (2021) Twelve years of SAMtools and BCFtools, GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>",
            "<li>Zheng Z, Li S, Su J, Leung AW-S, Lam T-W, Luo R (2022) Symphonizing pileup and full-alignment for deep learning-based long-read variant calling, Nature Computational Science, 2, 797-803. doi: 10.1038/s43588-022-00387-x</li>",
            "<li>Ewels P, Magnusson M, Lundin S, K\u00e4ller M (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report, Bioinformatics, 32(19), 3047-3048. doi: 10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
