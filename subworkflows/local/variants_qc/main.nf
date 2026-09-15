/*
  SUBWORKFLOW: VARIANTS_QC
*/

// Description: Review/export side-branch, independent of the phenotype-prediction path.
// Normalizes Clair3's PASS-filtered calls, restricts them to the ABO panel
// positions, and flattens them to a flat TSV -- a panel-scoped VCF a human
// can open in IGV, and a TSV for spot-checking, alongside (not feeding)
// ABO_CLAIR2METRICS's own gVCF-derived metrics used by PREDICTABOPHENOTYPE.

include { BCFTOOLS_NORM  } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW  } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_QUERY } from '../../../modules/nf-core/bcftools/query/main'

workflow VARIANTS_QC {
    take:
    ch_clair3_vcf // channel: [ val(meta),  path(vcf) ] -- CLAIR3.out.vcf
    ch_clair3_tbi // channel: [ val(meta),  path(tbi) ] -- CLAIR3.out.tbi
    ch_fasta      // channel: [ val(meta1), path(fasta) ]
    panel_bed     // path: ABO panel BED (broadcastable value channel)

    main:

    //
    // MODULE: BCFTOOLS_NORM
    //
    // Left-align indels and split multiallelics.
    BCFTOOLS_NORM(
        ch_clair3_vcf.join(ch_clair3_tbi, by: 0),
        ch_fasta.map { meta1, fasta -> [meta1, fasta] }.first(),
    )

    //
    // MODULE: BCFTOOLS_VIEW
    //
    // Restrict normalized calls to the ABO panel positions.
    BCFTOOLS_VIEW(
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.index, by: 0),
        panel_bed,
        [],
        [],
    )

    //
    // MODULE: BCFTOOLS_QUERY
    //
    // Flatten the panel-filtered calls to a per-position TSV.
    BCFTOOLS_QUERY(
        BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.index, by: 0),
        [],
        [],
        [],
    )

    emit:
    panel_vcf       = BCFTOOLS_VIEW.out.vcf
    panel_vcf_index = BCFTOOLS_VIEW.out.index
    panel_query_tsv = BCFTOOLS_QUERY.out.output
}
