# nf-core/abotyper v2.0.0 — panel-driven ABO genotyping

Every diagnostic variant position, its interpretation text, calling
thresholds, and associated ISBT subtype is defined in an external panel
(`assets/refs/abo_variant_panel.yaml`, or the equivalent `.csv`) instead of
being hardcoded in each script. Adding a new marker to the assay means
adding a row to the panel — no script changes, no redeploy.

The pipeline aligns every sample once against a single combined ABO
reference (default: `assets/refs/NG_006669v2.fasta`, the full ABO
RefSeqGene, LRG_792) spanning exons 2–7 in one continuous read, as
described in Mobegi et al. 2025 (*Int J Mol Sci* 26(12):5443,
doi:10.3390/ijms26125443). This enables read-level cis-phasing across the
whole gene and typing of exon 2–6 variants that earlier assay designs
could not reach.

## Scripts in this directory

| File | Role |
| --- | --- |
| `abo_panel.py` | Shared loader + scoring engine. Every other script imports this. Also a CLI: `--summary`, `--export-csv`, `--export-yaml`. |
| `pysam_haploscan.py` | BAM → per-position stats + per-read haplotypes, in a single pass using pysam. This is the live path (`HAPLOSCAN` module) for the combined-reference topology. |
| `predict_abo_phenotype.py` | Per-position stats → human-readable `*.ABOPhenotype.txt`. Auto-detects which exon(s) a report covers from panel position overlap, so one combined-reference run produces one report covering every exon in a single file. |
| `aggregate_abo_reports.py` | Aggregates every sample's report + haplotype data into `ABO_result.txt/.xlsx` and `final_export.csv`. Column layout, position parsing, and subtype-marker scanning are all panel-driven. Supports both the current single `combined/` per-sample layout and the legacy per-exon subfolder layout. |
| `calibrate_panel_positions.py` | Populates `amplicon_pos` for panel rows against a specific reference FASTA. Not required for the default `NG_006669v2.fasta` reference — see "Coordinate systems" below — but needed if switching to a different combined reference. |
| `stats_from_pileup.py` | samtools mpileup → per-position stats. Retained for the legacy dual-mini-amplicon topology; not used by the pipeline's current single-reference default. |
| `rename_samples.py` | Optional sample renaming using a tab-delimited `sequencingID`/`sampleName` file (`--renaming_file`). |

## Coordinate systems in the panel

Each panel row can carry up to four position fields, resolved in this
priority order by `VariantMarker.resolved_position()`:

1. **`amplicon_pos`** — position in whatever reference FASTA is currently
   in use, if calibrated via `calibrate_panel_positions.py`. Takes
   priority when set, since it's explicitly calibrated for the active run.
2. **`ng006669_2_pos`** — position in the full `NG_006669.2` RefSeqGene
   (42,144 bp), the pipeline's default `abo_reference_fasta`. Populated
   and verified (zero `ref_base` mismatches) for 41 of 51 panel rows.
   This is the coordinate actually used by default — no calibration step
   is required when running against `NG_006669v2.fasta`.
3. **`legacy_exon6_pos` / `legacy_exon7_pos`** — position in the old,
   retired exon6-only / exon7-only mini-amplicon references used by
   pipeline v1.x. Only valid when running the legacy per-exon topology
   against those old references directly — **not** valid against
   `NG_006669v2.fasta` or any other combined reference.
4. **`cdna_pos` / `cdna_change`** — stable HGVS cDNA numbering
   (`NM_020469.3`), matching ISBT ABO allele nomenclature. Never changes
   with assay redesign; use this as the durable identifier.

### Switching to a different combined reference

If a different combined reference FASTA is used instead of the shipped
`NG_006669v2.fasta`, `amplicon_pos` needs to be (re)calibrated:

```bash
# Anchor mode: derive amplicon_pos from ng006669_2_pos via a substring
# match of the new reference inside NG_006669.2 (fast, no spliced aligner):
python3 calibrate_panel_positions.py \
    --panel abo_variant_panel.yaml \
    --out abo_variant_panel.calibrated.yaml \
    --ng006669-reference NG_006669.2.fasta \
    --combined-reference new_reference.fasta

# Spliced-alignment mode: required for any position that still lacks
# ng006669_2_pos/legacy coordinates (currently the 10 exon2-5/intron/
# structural rows listed below):
minimap2 -ax splice NM_020469.3.fasta new_reference.fasta > abo_vs_ref.sam
python3 calibrate_panel_positions.py \
    --panel abo_variant_panel.calibrated.yaml \
    --out abo_variant_panel.calibrated.yaml \
    --from-sam abo_vs_ref.sam --force
```

Then point `--panel` at the calibrated file (or set `params.abo_panel`).
Check calibration status at any time with:

```bash
python3 abo_panel.py --panel abo_variant_panel.yaml --summary
```

## Adding a new variant

Add a row to `abo_variant_panel.yaml` (or the `.csv`). Minimum fields:
`id`, `cdna_pos`, `cdna_change`, `exon`, `ref_base`, `alt_base`,
`variant_type`, `category`, `call_rule` (`primary_biallelic` or
`named_marker`), plus whichever thresholds/labels the call rule needs —
see the schema comments at the top of the YAML file. Populate
`ng006669_2_pos` directly if the position in `NG_006669.2` is known, or
run `calibrate_panel_positions.py` to derive it.

## Known gaps and caveats

- **10 panel rows remain uncalibrated**: `onull_106`, `onull_188`,
  `onull_189`, `onull_220`, `o02_53`, `o04_88`, `aweak_intron6_374`,
  `aweak_intron2_98`, `b3_intron3_155`, `abantu_intron4_203` (exon 2–5 /
  intronic positions), plus the `o16` structural deletion and the exon 1
  start-codon variants. These need either the annotated `NG_006669.2`
  GenBank flatfile (for exact intron lengths) or a spliced alignment of
  `NM_020469.3` against the reference to resolve.
- **`o34_homopolymer_804`** (legacy `pos431`): sits in/near a G-homopolymer
  associated with c.804dupG/delG. ONT is error-prone at homopolymers —
  cross-check against the `ael_804_indel` marker and raw indel percentages
  before trusting this call in isolation.
- **`ABO*O.16`** (c.204-98_\*3252del, ~3.15 kb deletion) is documented under
  `structural_variants` in the panel but is **not scored** by the SNP/indel
  engines — detecting it needs a long-read gap/soft-clip detector (a
  reasonable future addition). Validate against a known O.16 sample before
  relying on this; until then it's a plausible false-homozygote source.
- **c.1A>G / c.2T>C** (exon 1 start-codon variants) sit upstream of the
  amplicon's forward primer and are a genuine assay coverage gap, not a
  panel or script limitation.
- **Corrected mislabelling inherited from v1.2.0**: c.657C>T (legacy
  `pos283`) was previously flagged as an `ABO*BW.06` marker. It isn't — it's
  part of the shared B-lineage backbone present in nearly every
  B/B3/Bweak/Bel/cisAB/BA allele. The real Bw.06-private variant is
  c.1036A>G, a separate panel entry (currently among the 10 uncalibrated
  rows above).
