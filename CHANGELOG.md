# nf-core/abotyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [unreleased]

The [ISBT ABO blood group allele database](https://blooddatabase.isbtweb.org/system/ABO) is now actively maintained, with new variants expected to be submitted on an ongoing basis. This release moves the pipeline's position-consuming scripts to a panel-driven variant coordinate system in order to align with this database and futureproof the pipeline against upcoming variant additions, and migrates the pipeline topology from the previous dual exon6/exon7 mini-amplicon design to a single combined ABO reference (NG_006669.2, RefSeqGene/LRG_792), as described in Mobegi et al. 2025 (IJMS 26(12):5443).

### `Added`

- Migrated the Nextflow pipeline topology to the single combined ABO reference: `nextflow_schema.json` and `nextflow.config` now expose `abo_reference_fasta`/`abo_reference_fai` in place of the old `exon6fasta`/`exon7fasta`/`exon6fai`/`exon7fai` parameters; `main.nf`, `workflows/abotyper.nf`, and the `minimap_align_exons`, `variant_calling_haploscan`, and `predictabophenotype` subworkflows were rewritten for single-pass alignment and variant calling per sample (no more `meta.exon`-based channel splitting/matching).
- Introduced `bin/abo_panel.py`: a shared variant-panel loader and scoring engine (`primary_biallelic` / `named_marker` rule types) consumed by all four position-aware scripts, replacing hardcoded position dictionaries with an external, editable table.
- Added the variant panel itself, `assets/refs/abo_variant_panel.yaml` (source of truth) and `abo_variant_panel.csv` (flat export for spreadsheet review) — 51 scorable diagnostic positions plus 2 documented structural variants, expanding coverage from the legacy 24 exon6/exon7 positions to include A2/A3/Am/Aweak/Ael/Bweak/Bel/B3/Bw/B(A)/cisAB/O-null markers spanning exons 2–7.
- Added `bin/calibrate_panel_positions.py` with three calibration modes for populating `amplicon_pos`: anchor mode (legacy mini-amplicon → combined reference, no dependencies), SAM/CIGAR mode (spliced alignment against a supplied transcript), and NG_006669.2 mode (direct offset from a verified genomic coordinate) — plus a `--verify` flag to spot-check panel `ref_base` values against a reference FASTA.
- Populated `ng006669_2_pos` for all 51 panel entries, derived from the annotated NG_006669.2 mRNA/CDS `join()` feature and independently verified against the actual NG_006669.2 sequence with zero ref_base mismatches (24 legacy + 16 new exon7 positions cross-checked first; the remaining 10 exon2–5/intron positions computed from the CDS block map and verified against the same sequence).
- Confirmed the combined reference will be the **untrimmed** NG_006669.2 record (42,144 bp); set `amplicon_pos = ng006669_2_pos` directly for all 51 panel entries on that basis (zero coordinate offset between the two systems). Positions upstream of exon 2 (~NG 1–18,000) will show zero coverage in aligned samples since the PCR amplicon itself only spans exon 2–7 — expected and harmless.
- Generalized `predict_abo_phenotype.py` report generation to be exon-agnostic (loops over whatever exons the panel + input data cover) rather than hardcoded to exon6/exon7 only.
- Generalized `aggregate_abo_reports.py`: column layout and per-position parsing are now built dynamically from the panel; `scan_bw_markers()`/`scan_a2_markers()` are replaced by one generic `scan_category_markers(category)` covering every panel category, not just Bw/B(A)/cisAB.
- Added `--legacy` mode to `pysam_haploscan.py` to preserve the previous dual-mini-amplicon behaviour (length-based exon detection, `legacy_exon6_pos`/`legacy_exon7_pos` resolution); the pipeline now defaults to the combined single-reference mode described above.

### `Fixed`

- Corrected a Bw subtype mislabelling inherited from v1.2.0: c.657C>T (legacy `pos283`) was flagged as an `ABO*BW.06` marker but is in fact part of the shared B-lineage 7-SNP backbone present in nearly all B/B3/Bweak/Bel/cisAB/BA alleles. The true Bw.06-private variant, c.1036A>G, is now tracked as its own panel entry.

### `Changed`

- Report and aggregate-table column headers now key off genomic coordinates (e.g. `Exon6_pos22694`) rather than the old small legacy mini-amplicon offsets (e.g. `Exon6_pos22`), since `resolved_position()` prefers `amplicon_pos` once populated. Verified end-to-end (predict → aggregate → Excel/CSV/LIS export) against the new coordinate scheme with no code changes required — the column-naming and report-parsing logic was already coordinate-agnostic by design. Any external tooling matching on the old literal column names will need updating.

### `Dependencies`

- Added `pyyaml` as a dependency of `bin/abo_panel.py`; needs adding to the conda environments of `modules/local/{haploscan,mpileupstats,abo/getabosnps,abo/snps2pheno}` and to the corresponding pinned containers before deployment (see integration notes for the container-rebuild caveat).

### `Known gaps / not yet done`

- `ABO*O.15` (large deletion) and the exon 1 start-codon variants (c.1A>G, c.2T>C) remain documented as out-of-scope structural/assay gaps, not scored by any current engine.
