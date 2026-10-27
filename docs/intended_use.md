# Intended use

**nf-core/abotyper is for research use only (RUO).**

- It has not been validated as an in-vitro diagnostic (IVD) and carries no regulatory
  clearance or approval. Its output must not be used on its own for clinical diagnosis,
  patient management, or transfusion or transplantation decisions.
- Any result used clinically must be confirmed by an accredited method (for example
  serology or an accredited genotyping assay) under the laboratory's own quality system.
- Reports that rest on a single line of evidence (Clair3 only, BAM pileup only) or on
  disagreeing lines are marked `REVIEW`, and samples with no reads over the ABO region
  (likely amplicon dropout) are reported as `Unknown`. These flags support, and do not
  replace, expert review.

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
