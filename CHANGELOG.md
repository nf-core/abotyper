# nf-core/abotyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> Versions below `v1.0.0` are retrospective labels for unreleased development milestones. They group the work of the merged pull requests in the git history (December 2024 – October 2026) and were never released; `v1.0.0` is the first release.

## v1.0.0 - [2026-10-16]

Initial release of nf-core/abotyper, created with the [nf-core](https://nf-co.re/) template.

The pipeline types the ABO blood group from Oxford Nanopore long-read amplicons of the ABO gene (exons 2–7, aligned to the RefSeqGene record NG_006669.2 / LRG_792), as described in Mobegi et al. 2025, IJMS 26(12):5443. Positions, thresholds and allele associations come from an external variant panel aligned to the [ISBT ABO allele table](https://blooddatabase.isbtweb.org/system/ABO), so new variants are added by editing the panel and not the code.

**Research use only.** The pipeline has not been validated for clinical diagnosis; see [docs/intended_use.md](docs/intended_use.md).

### `Added`

- Two lines of evidence for every panel position: Clair3 and an independent BAM pileup (`ABO_PILEUP`, `bin/pileup2metrics.py`), combined and compared by `ABO_EVIDENCE`. Clair3 is the primary line, the pileup fills positions Clair3 had no evidence for, and a per-position concordance table and an evidence status file are written for every sample. The pipeline warns at run time when Clair3 produces no gVCF for a sample.
- Clair3 is restricted to the ABO exon 2–7 region (`abo_clair3_region_bed`, `assets/refs/abo_clair3_region.bed`).
- `ABO_result.xlsx` with three sheets: `Typing` (one row per sample: reads, read length, ABO-mapped reads, on-target %, depths, evidence source, typing result, phenotype, genotype, variants detected, QC), `Variants` (one column per panel position with the call and depth, such as `A 50% (215x)`; alleles that reach the marker's threshold are highlighted and low depth is coloured; the full per-base table stays in `ABO_result.txt`) and `About` (script, panel and ISBT-table versions and the research-use statement).
- ISBT-style allele names (for example `ABO*A1.01/ABO*O.01.01`) that start from the reference allele of each family and are narrowed to the specific allele when non-primary panel markers fire. The panel records its own version and the ISBT allele-table version it was checked against (v9).
- Variant panel v2.3.0, checked against ISBT ABO allele table v9 and the v1.1 PDF: every DNA change in the tables inside the exon 2-7 amplicon is a panel marker (135 markers; the 82 rare ones are marked `report: detected`, so they only get columns in `ABO_result.xlsx` and `ABO_result.txt` when a sample carries them), except the silent c.930G>A. The two exon 1 start-codon changes and the c.204-98_*3252del deletion are documented as structural or assay gaps, and `assets/refs/abo_isbt_v9_coverage.tsv` lists every change and where it is handled. Reference bases were verified and the background at the new positions in the reference runs is below 3.5%. Indel rows carry their own `search_window`, so neighbouring indel markers never read each other's calls. New markers c.542G>A (`onull_542`, ABO*O.06 / O.01.11 / O.01.75) and c.768C>A (`o1_768`, ABO*O.01.26) were confirmed in two independent replicate pairs of reference samples. c.539 and c.1054 now carry both ISBT alt bases (G>C / G>A and C>T / C>G).
- QC and review flags: low-read A1 calls, A2 evidence gaps, disagreement between Clair3 and the BAM pileup at primary markers, primary positions that rest on the pileup alone, and a specific "No reads over ABO region (likely amplicon dropout)" flag in place of a reference-only call. Replicates of the same sample on different barcodes that disagree (genotype, A2 evidence or allele level) are marked in red, with a `ReplicateNote` recording each run's genotype and zygosity and possible allele dropout.
- Run-level outputs: `ABO_run_qc.txt` and an `ABO_evidence_summary_mqc.tsv` table that is added to the MultiQC report.
- `aggregate_abo_reports.py --threads` for parallel parsing of large runs, and `--nanoplot-dir` to take library read counts from NanoPlot.
- Intended-use statement ([docs/intended_use.md](docs/intended_use.md)), linked from the README and written into every output.
- nf-test suite for the local modules and subworkflows, with small committed fixtures in `tests/fixtures/` and `tests/run_nf_tests.sh` for running them on Linux.

### `Fixed`

- The BAM pileup row of c.1061delC left the reference (C) column at 0, so a position filled from the BAM instead of Clair3 read as heterozygous whatever the deletion fraction.
- The reference FASTA is stored as a regular file. It had been committed as a Git LFS pointer, which broke the pipeline for any checkout without LFS.
- All `bin/` scripts are executable in git and read and write files as UTF-8, so results are identical on Windows and Linux.
- `ABO_SNPS2PHENO` stub runs create all declared outputs.
- A1 and A2 subtypes follow the ISBT definitions: c.1061delC with a reference c.467 is named ABO*A2.06 and c.467C>T with a clearly absent c.1061delC (30 reads or more) is named ABO*A1.02, both flagged for review as rare alleles. c.907G>A and c.1032G>A are not part of any ISBT definition, so they are reported in the QC flags but no longer name or support an allele.
- Ael alleles: c.804dupG and c.804delG are named ABO*AEL.01, ABO*AEL.08 (with c.467C>T) and ABO*AEL.03, with a review flag because Ael is a very weak or silent A. The BAM pileup now searches the whole 7-G homopolymer run (NG 24283-24289) for that indel, because the aligner anchors it on the base before the run or on its first G, outside the old +-2 bp window; the marker threshold is 25% because the run gives up to 13% deletion signal in reference runs. A1.02 and A2.06 review flags now give the read percentages that caused them.
- Allele names are decided per haplotype when the Clair3 phasing is clean (two haplotypes covering at least 90% of the phased reads, every fired marker phased): each marker is given to the haplotype that carries it and counts as absent for the other. This stops the markers of an O allele renaming the A allele (A2.01 reported as A2) and names the two alleles of an O1/O1 sample separately. The sample falls back to pooled naming when the phasing is not clean, and `PhaseConfidence` says when the haplotypes decided the names.
- Depth check for every sample (`Depth check` column in the `Typing` sheet, coloured red for very low or no data, orange for low): fewer than `--min-depth` reads at a primary position sends the sample to review. Sample names containing `NEG` (`--negative-control-pattern`) are negative controls, where few reads is the expected result (green) and reads above `--min-depth` are flagged as possible contamination; names containing `POS` (`--positive-control-pattern`) are positive controls, which need normal depth and fail with too few reads. Calls are never changed by the check.
- Panel `associated_alleles` completed from ISBT v9 for every non-primary marker and six unsupported claims removed (c.297, c.722, c.742, c.907, c.1032, c.1054). Allele naming now also rules out any allele that is defined by a marker that was covered but did not fire, so a sample with the c.542G>A + c.646/681/771/829 cassette pattern is named ABO*O.01.11 and with c.106/188/189/220 added ABO*O.01.75, not the family name ABO*O.01.
- Path parameters shipped in `assets/` no longer carry a hardcoded `default` in `nextflow_schema.json`; their defaults are set in `nextflow.config` from `${projectDir}`.

### `Dependencies`

| Tool                      | Version / source                                                         |
| ------------------------- | ------------------------------------------------------------------------ |
| Nextflow                  | `>=25.10.4`                                                              |
| Clair3                    | container `hkubal/clair3:v2.0.0` (PyTorch models downloaded at run time) |
| samtools / htslib         | 1.24 (alignment statistics, coverage and `ABO_PILEUP`)                   |
| minimap2                  | via the nf-core `minimap2/align` module                                  |
| pyyaml                    | required by `bin/abo_panel.py` and all panel-driven scripts              |
| pandas, numpy, xlsxwriter | `ABO_SNPS2PHENO` and `ABO_GETABOSNPS`                                    |

### `Known gaps / not yet done`

- The large ABO deletion (c.204-98_*3252del) and the exon 1 start-codon variants (c.1A>G, c.2T>C) are documented as out-of-scope structural/assay gaps; no current engine scores them.
- Samples with very few reads can still receive a genotype (a negative control with 3 ABO reads was typed); a minimum-depth gate is planned.
- The A2 subtype markers c.907G>A and c.1032G>A are pipeline-defined background markers; ISBT table v9 does not list them as allele-defining.

## v0.10.0dev - [2026-10-15]

Pull requests #40, #43, #44, #45.

- A1/A2 aggregation logic streamlined and an A1/A2 zygosity fix.
- FastQC replaced by NanoPlot for long-read QC; `_NanoStats.txt` is published and the NanoPlot module patched.
- Modules config and Python scripts updated; reference files housekeeping; README and metro map refreshed.

## v0.9.0dev - [2026-09-15]

Pull requests #33–#35, #37–#39. The move to the full-length reference and the Clair3 call path.

- Full-length ABO reference (NG_006669.2) and the ISBT ABO allele table v9 adopted as the reference data; `schema_input.json` updated.
- Clair3 replaces HAPLOSCAN: Clair3 gVCF metrics and bcftools panel filtering wired into `VARIANTS_QUANTIFICATION`; new `ABO_PANEL2BED` module (with a fix for unpadded indel and homopolymer BED intervals); bcftools modules added; final audit of the change-over.
- Fix for an A2.01 subtype miscall; `bin/` documentation tidied; pipeline configuration gaps closed.
- The real ABO reference FASTA added (stored in Git LFS at this stage; reversed in v1.0.0).

## v0.8.0dev - [2026-07-30]

Pull request #31.

- nf-core template sync to tools 4.1.0, with post-lint fixes to the utils subworkflow.
- Biocontainers replaced by Seqera containers for the HAPLOSCAN module.

## v0.7.0dev - [2026-05-26]

Pull requests #25 and #27.

- nf-core template syncs to tools 4.0.0 and 4.0.2.
- Final fine-tuning of A1/A2 subtyping; introductory text added to the MultiQC configuration.

## v0.6.0dev - [2026-04-02]

Pull requests #20, #21 and #23.

- nf-core template syncs (tools 3.4.1 and 3.5.2); version reporting migrated to topic channels in the local modules.
- Pre-commit fixes and lint clean-up; unused positions removed from `aggregate_abo_reports.py`; changelog, README and RO-Crate metadata updated.

## v0.5.0dev - [2025-09-03]

Pull request #13.

- Schema update and release lint; test-data paths fixed; review comments addressed.
- MultiQC configuration consolidated; containers added to local modules; Python version updated for the mpileup statistics module.

## v0.4.0dev - [2025-08-14]

Pull requests #5, #7, #10 and #11. Release candidate hardening.

- Modules reworked to avoid aliasing; Python scripts upgraded and example data added.
- Docker and Singularity profiles finalised; MultiQC report customised and pinned to 1.28 for broader AVX CPU compatibility.
- Test and debug checks, lint fixes and pre-commit formatting before the first release candidate.

## v0.3.0dev - [2025-07-30]

Pull requests #3 and #4.

- Adopted the nf-core pipeline template (sync from tools 3.2.1 to 3.3.2) with merge fixes and pre-lint clean-up.

## v0.2.0dev - [2025-04-24]

The first working analysis pipeline, built as a dual exon 6 / exon 7 mini-amplicon design.

- Alignment and BAM-statistics modules, then mpileup, mpileup statistics, phenotype/genotype modules and an aggregation script.
- A1/A2/A3 subtyping, a fix for mpileup computations, accurate capture of OB blood types, and the aggregation process finalised in the workflow.
- First nf-core template update (tools 3.2.1).

## v0.1.0dev - [2024-12-05]

- Repository created with the README and first modules (December 2024); project logo and metro map followed (January–February 2025).
