# Intended use

**nf-core/abotyper is for research use only (RUO).**

- It has not been validated as an in-vitro diagnostic (IVD) and carries no regulatory
  clearance or approval. Its output must not be used on its own for clinical diagnosis,
  patient management, or transfusion or transplantation decisions.
- Any result used clinically must be confirmed by an accredited method (for example
  serology or an accredited genotyping assay) under the laboratory's own quality system.
- BAM-pileup-only reports, primary positions filled from BAM, and disagreements at
  primary typing markers are marked `REVIEW`. Clair3-only evidence and disagreements
  confined to non-primary markers do not automatically trigger this flag. Samples with
  no reads over the ABO region (likely amplicon dropout) are reported as `Unknown`.
  These flags support, and do not replace, expert review.

## Alleles this assay cannot resolve

The pipeline types the exon 2-7 amplicon. The variant panel scores 135 positions in exons 2-7
and the adjacent introns (NG_006669.2 18064-24544). Three ISBT alleles are defined by changes
outside what that amplicon and the panel can see, and are reported as another allele without a
QC flag:

| ISBT allele | Phenotype | Defining changes                                 | Reported as                                                                                  |
| ----------- | --------- | ------------------------------------------------ | -------------------------------------------------------------------------------------------- |
| ABO\*AW.13  | Aweak     | c.2T>C (exon 1 start codon)                      | ABO\*A1.01, phenotype A1                                                                     |
| ABO\*AW.16  | Aweak     | c.1A>G (exon 1 start codon), c.467C>T, c.1061del | ABO\*A2.01, phenotype A2                                                                     |
| ABO\*O.16   | O         | c.106G>T, c.188G>A, c.204-98\_\*3252del          | the other allele may appear homozygous (no exon 5-7 reads from O.16), e.g. A1/A1 for A1/O.16 |

For AW.13 and AW.16 the genotype class (A) is right but the A subgroup is not: the weak expression
comes from the exon 1 change. A weaker serological A reaction than the reported A1 or A2 suggests
one of these alleles. Typing them needs reads that cover exon 1 (a forward primer upstream of exon 1),
scored panel markers for c.1A>G and c.2T>C, and a Clair3 region (`--abo_clair3_region_bed`) that
includes exon 1. The large O.16 deletion cannot be called from the pileup.

## What an IVD pathway would need

If the pipeline is ever to be used diagnostically, the work goes beyond this repository:

- a defined intended purpose, population and claims;
- an analytical and clinical validation study (accuracy against a reference method,
  precision and reproducibility, limits of detection, interfering conditions, sample-type
  and instrument/chemistry coverage);
- software documentation and risk management under the applicable regulation
  (for example software life-cycle processes, a hazard analysis and a traceability matrix);
- controlled change management: locked, versioned container images, panel and
  nomenclature, with revalidation on any change that can affect a call;
- a quality management system, user training and post-market surveillance.

## Provenance recorded with every run

`ABO_result.xlsx` has an `About` sheet, and `ABO_run_qc.txt` starts with the same lines:
the research-use statement, the script version, the panel version and source, and the ISBT
allele-table version and the date it was checked. Allele names follow the ISBT ABO allele
table; positions marked "legacy" in the panel (for example the c.804 homopolymer row) should be
verified against the current table.
