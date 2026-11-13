# nf-core/abotyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [unreleased]

Initial release of nf-core/abotyper, created with the [nf-core](https://nf-co.re/) template.

The pipeline types the ABO blood group from Oxford Nanopore long-read amplicons of the ABO gene (exons 2–7, aligned to the RefSeqGene record NG_006669.2 / LRG_792), as described in Mobegi et al. 2025, IJMS 26(12):5443. Reads are checked with NanoPlot, aligned with minimap2 and called with Clair3, with an independent BAM pileup as a second line of evidence. The typing steps come from the [abotyper-core](https://github.com/fmobegi/abotyper-core) package. Positions, thresholds and allele associations come from a variant panel aligned to the [ISBT ABO allele table](https://blooddatabase.isbtweb.org/system/ABO), so new variants are added by editing the panel and not the code.

**Research use only.** The pipeline has not been validated for clinical diagnosis; see [docs/intended_use.md](docs/intended_use.md).

### `Added`

- The Python typing steps are the [abotyper-core](https://github.com/fmobegi/abotyper-core) package (0.3.1), released on its own to PyPI and Bioconda with its own tests (including golden-output tests on the pipeline's test samples), CLI and CI. The `abotypercore/*` local modules install it from Bioconda, use Seqera Containers (Wave) images built from that conda environment with `nf-core modules containers create`, call its `abotyper-*` commands and report the `abotyper-core` version. The package has one version for all its steps, which `ABO_run_qc.txt` and the About sheet report as `abotyper-core: <version>`; the pipeline test leaves that line out of its snapshot, so a new package release does not change the snapshot. Standalone tools (`abotyper-rename-samples`, `abotyper-calibrate-panel`) come with the package.
- Variant panel v2.3.0, bundled with abotyper-core and versioned with it (`--abo_panel` overrides it), checked against ISBT ABO allele table v9: every DNA change inside the exon 2–7 amplicon is a panel marker (135 markers; the 82 rare ones are marked `report: detected` and only get columns when a sample carries them), except the silent c.930G>A. abotyper-core's `abo_isbt_v9_coverage.tsv` lists every ISBT change and where it is handled.
- Two lines of evidence for every panel position: Clair3, restricted to the ABO exon 2–7 region (`--abo_clair3_region_bed`), and a BAM pileup (nf-core `samtools/mpileup`, `-a -B -d 0 -Q 7 -q 1`), combined and compared by `ABOTYPERCORE_PILEUP2METRICS`. Clair3 is the primary line, the pileup fills positions Clair3 had no evidence for, and a per-position concordance table and an evidence status file are written for every sample. The pipeline warns at run time when Clair3 produces no gVCF for a sample.
- Clair3 model choice: `--clair3_model_url` downloads a model, else `--clair3_model_path` uses a local model folder, else the model packaged in the Clair3 container is used (`--clair3_model`, default `r1041_e82_400bps_sup_v500`, one of the ONT models in the container). The nf-core clair3 module only provides packaged models under Docker, Singularity and Podman, so with other profiles that model is downloaded by name.
- ISBT-style allele names (for example `ABO*A1.01/ABO*O.01.01`), decided per haplotype when the Clair3 phasing is clean, with A1/A2, Ael and O subtype rules that follow the ISBT definitions.
- Result tables: `ABO_result.xlsx` (`Typing`, `Variants` and `About` sheets), `ABO_result.txt`, `final_export.csv` for laboratory information systems, a run-level `ABO_run_qc.txt` and an `ABO_evidence_summary_mqc.tsv` table in the MultiQC report.
- QC and review flags: depth check per sample, low-read A1 calls, A2 evidence gaps, disagreement between Clair3 and the BAM pileup at primary markers, positions that rest on the pileup alone, likely amplicon dropout, and disagreeing replicates of the same sample.
- Samplesheet columns `control` (`negative` or `positive`) and `barcode`. A negative control with no reads passes QC and one with reads above `--min-depth` is flagged as possible contamination; a positive control needs normal depth. Without a `barcode` value, a sample name ending in `_barcode<digits>` is split into name and barcode.
- `-profile conda,conda_lock` installs the exact package builds of the modules from conda `@EXPLICIT` lock files (`assets/conda-lock/`, `conf/conda_lock_explicit_{amd64,arm64}.config`). The pixi locks referenced by the configs nf-core tools generates cannot be read by `conda create --file`; see [Conda](docs/usage.md#conda).
- Intended-use statement ([docs/intended_use.md](docs/intended_use.md)), linked from the README and written into every output.
- nf-test suite for the pipeline and its local modules, subworkflows and functions, with small committed fixtures in `tests/fixtures/` and `tests/run_nf_tests.sh` for running them on Linux.

### `Dependencies`

| Tool              | Version / source                                                                         |
| ----------------- | ---------------------------------------------------------------------------------------- |
| Nextflow          | `>=25.10.4`                                                                              |
| Clair3            | container `hkubal/clair3:v2.0.0` (packaged PyTorch models; other profiles download them) |
| samtools / htslib | 1.25 (alignment statistics, coverage and `samtools/mpileup`)                             |
| minimap2          | via the nf-core `minimap2/align` module                                                  |
| abotyper-core     | 0.3.1 (Bioconda, Seqera Containers); the `abotypercore/*` typing modules                 |

### `Known gaps / not yet done`

- Only Oxford Nanopore reads are supported (`--minimap2_preset` `map-ont` or `lr:hq`): Clair3 runs with platform `ont` and an ONT model. Ion Torrent and Illumina reads from the pilot study typed correctly with the earlier pipeline; their Clair3 settings are still to be validated.
- The large ABO deletion (c.204-98_*3252del) and the exon 1 start-codon variants (c.1A>G, c.2T>C) are out-of-scope structural/assay gaps; no current engine scores them.
- The A2 subtype markers c.907G>A and c.1032G>A are pipeline-defined background markers; ISBT table v9 does not list them as allele-defining.
