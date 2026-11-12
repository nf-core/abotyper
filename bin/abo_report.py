#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Output writers of the aggregation: ABO_result.txt/.xlsx, final_export.csv, the run QC report and the MultiQC evidence table.

Part of aggregate_abo_reports.py (nf-core/abotyper); split out of that script so each part can be read and
tested on its own. This file is part of nf-core/abotyper and is released under the MIT License (see LICENSE).
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List

import pandas as pd
from xlsxwriter.utility import xl_col_to_name

from abo_naming import _cdna_tokens



# Reference (default) ISBT allele per internal token: used only when no panel marker
# points to a more specific allele of the same family. Checked against the ISBT ABO allele
# table v9 (blooddatabase.isbtweb.org/system/ABO): ABO*A1.01 (reference, no changes),
# ABO*B.01 (B backbone), ABO*O.01.01 (c.261delG only), ABO*O.02.01 (c.53G>T, c.220C>T,
# c.297A>G, c.526C>G, c.802G>A) and ABO*O.03 (no numeric suffix in the table).
RESEARCH_USE_STATEMENT = (
    "RESEARCH USE ONLY. Not validated for, and not to be used for, clinical diagnosis or "
    "transfusion decisions; confirm results with an accredited method before any clinical use."
)


class ReportMixin:
    """Output methods of ABOReportParser."""

    # -----------------------------------------------------------------
    # Compact per-position view (Variants sheet, "Variants detected" column)
    # -----------------------------------------------------------------
    @staticmethod
    def _num(v):
        try:
            f = float(v)
        except (TypeError, ValueError):
            return None
        return None if pd.isna(f) else f

    @staticmethod
    def _change_token(marker, alt: str) -> str:
        """The cDNA change of a marker that belongs to the alt that was seen (c.539G>C / c.539G>A -> c.539G>A)."""
        tokens = _cdna_tokens(marker.cdna_change)
        if not tokens:
            return f"c.{marker.cdna_pos}"
        for tk in tokens:
            if len(alt) == 1 and tk.endswith(">" + alt):
                return tk
            if alt == "del" and "del" in tk:
                return tk
            if alt == "ins" and ("dup" in tk or "ins" in tk):
                return tk
        return tokens[0]

    def _position_calls(self, final_df) -> List[Dict[str, dict]]:
        """One dict per row of the merged (MultiIndex) table: column label -> summary of that position.

        state is 'nodata' (no reads), 'ref' or 'alt'. A position reads 'alt' when the strongest alt allele
        reaches the marker's own calling threshold (variant_threshold_pct for named markers, 20% for the
        primary markers, which carry no such threshold)."""
        markers = {lab: self._marker_for_column(lab) for lab in self.position_columns}
        calls = []
        for i in range(len(final_df)):
            row_calls = {}
            for lab, mk in markers.items():
                reads = self._num(final_df.at[final_df.index[i], (lab, "#Reads")])
                type_ = final_df.at[final_df.index[i], (lab, "Type")]
                type_ = "" if type_ is None or (not isinstance(type_, str) and pd.isna(type_)) else str(type_)
                if mk is None or not reads:
                    row_calls[lab] = {"state": "nodata", "reads": int(reads or 0), "type": type_}
                    continue
                vt = (mk.variant_type or "").lower()
                indel = vt in ("deletion", "insertion", "indel") or "indel" in vt or (mk.ref_base or "").upper() == "REF"
                cand = []
                if not indel or "snp" in vt:
                    cand += [(b, self._num(final_df.at[final_df.index[i], (lab, b)]) or 0.0) for b in mk.alt_bases() if b in "ACGT"]
                if indel:
                    cand += [("del", self._num(final_df.at[final_df.index[i], (lab, "Del")]) or 0.0),
                             ("ins", self._num(final_df.at[final_df.index[i], (lab, "Ins")]) or 0.0)]
                alt, pct = max(cand, key=lambda x: x[1]) if cand else ("", 0.0)
                threshold = mk.variant_threshold_pct if mk.variant_threshold_pct is not None else 20
                info = {"reads": int(reads), "type": type_, "alt": alt, "pct": pct}
                if alt and pct >= threshold:
                    info.update(state="alt", change=self._change_token(mk, alt), cdna_pos=mk.cdna_pos)
                else:
                    info["state"] = "ref"
                row_calls[lab] = info
            calls.append(row_calls)
        return calls

    @staticmethod
    def _cell_text(info: dict) -> str:
        if info["state"] == "nodata":
            return "no data"
        if info["state"] == "alt":
            return f"{info['alt']} {info['pct']:.0f}% ({info['reads']}x)"
        return f"ref ({info['reads']}x)"

    @staticmethod
    def _variants_text(row_calls: dict) -> str:
        hits = sorted((v["cdna_pos"] or 0, v["change"], v["pct"], v["reads"]) for v in row_calls.values() if v["state"] == "alt")
        if not hits:
            return "" if all(v["state"] == "nodata" for v in row_calls.values()) else "None"
        # a call on 40 reads or fewer carries its depth, so noise in a near-empty sample is not mistaken for a variant
        return "; ".join(f"{c} {p:.0f}%" + (f" ({n}x)" if n <= 40 else "") for _, c, p, n in hits)

    def _undetected_rare(self, calls) -> set:
        """Column labels of rare panel positions (report: detected) that no sample in the run carries."""
        if not calls:
            return set()
        hidden = set()
        for lab in self.position_columns:
            mk = self._marker_for_column(lab)
            if mk is not None and (mk.report or "").lower() == "detected" and not any(c[lab]["state"] == "alt" for c in calls):
                hidden.add(lab)
        return hidden

    def _build_variants_view(self, final_df, calls, hidden=()):
        """Compact Variants sheet: one column per panel position holding the call, e.g. 'A 50% (215x)'.

        Returns (frame, position_labels) where position_labels maps frame column index -> panel column label."""
        data = {"Sample name": [], "Barcode": [], "Phenotype": [], "Genotype": [], "Typing Result": []}
        for i in range(len(final_df)):
            idx = final_df.index[i]
            try:
                bc = f"barcode{int(final_df.at[idx, ('', 'Barcode')]):02d}"
            except (TypeError, ValueError):
                bc = ""
            data["Sample name"].append(str(final_df.at[idx, ("", "Sequencing_ID")]))
            data["Barcode"].append(bc)
            data["Phenotype"].append(final_df.at[idx, ("Result", "Phenotype")])
            data["Genotype"].append(final_df.at[idx, ("Result", "Genotype")])
            data["Typing Result"].append(final_df.at[idx, ("Result", "ExtendedGenotype")])
        names = {}
        for lab in self.position_columns:
            if lab in hidden:
                continue
            mk = self._marker_for_column(lab)
            first = (_cdna_tokens(mk.cdna_change) or [lab])[0] if mk and mk.cdna_change else lab
            pos = mk.resolved_position() if mk else ""
            name = f"{first} ({pos})"
            while name in data:
                name += " "
            names[lab] = name
            data[name] = [self._cell_text(calls[i][lab]) for i in range(len(final_df))]
        for nc in self.notes_columns:
            data[nc] = [final_df.at[final_df.index[i], ("Notes", nc)] for i in range(len(final_df))]
        frame = pd.DataFrame(data).astype(object)
        cols = list(frame.columns)
        return frame, {cols.index(n): lab for lab, n in names.items()}

    def _build_typing_df(self, final_df) -> pd.DataFrame:
        """One NGSengine-style row per sample, same order as final_df (single-level columns)."""
        def num(v, nd=None):
            if v is None or pd.isna(v):
                return None
            return round(float(v), nd) if nd is not None else int(round(float(v)))

        variants_detected = getattr(self, "_variants_detected", [])
        rows = []
        for _, r in final_df.iterrows():
            sid = str(r["Sequencing_ID"])
            try:
                bc = int(r["Barcode"])
            except (TypeError, ValueError):
                bc = None
            m = self.sample_metrics.get((sid, bc), {})
            rows.append({
                "Sample name": sid,
                "Barcode": f"barcode{bc:02d}" if bc is not None else "",
                "Number of reads": num(m.get("reads")),
                "Read Length": num(m.get("read_length")),
                "Insert Size": num(m.get("insert_size")),
                "ABO reads": num(m.get("abo_reads")),
                "Mappability %": num(m.get("mappability"), 1),
                "Mean depth (covered bases)": num(m.get("covered_depth")),
                "Median Read depth": num(m.get("median_depth")),
                "Median analyzed Read depth": num(m.get("median_analyzed_depth")),
                "Number of question mark positions": m.get("question_marks"),
                "Evidence source": m.get("evidence_source") or "",
                "BAM-only positions": m.get("bam_only"),
                "Clair3/BAM discordant": m.get("discordant"),
                "Discordant primary positions": m.get("discordant_primary"),
                "Typing Result": r["ExtendedGenotype"],
                "Phenotype": r["Phenotype"],
                "Genotype": r["Genotype"],
                "Variants detected": variants_detected[len(rows)] if len(rows) < len(variants_detected) else "",
                "Depth check": r.get("DepthCheck"),
                "QCStatus": r["QCStatus"],
                "QCFlags": r["QCFlags"],
                "Replicate note": r.get("ReplicateNote"),
            })
        return pd.DataFrame(rows, columns=[
            "Sample name", "Barcode", "Number of reads", "Read Length", "Insert Size",
            "ABO reads", "Mappability %", "Mean depth (covered bases)", "Median Read depth",
            "Median analyzed Read depth",
            "Number of question mark positions", "Evidence source", "BAM-only positions",
            "Clair3/BAM discordant", "Discordant primary positions", "Typing Result", "Phenotype", "Genotype",
            "Variants detected", "Depth check", "QCStatus", "QCFlags", "Replicate note"]).astype(object)

    def save_results_to_file(self, final_df):
        """Save results to text and Excel files; column count/labels are dynamic."""
        try:
            calls = self._position_calls(final_df)
        except Exception as calls_err:
            print(f"ERROR: could not summarise the positions: {calls_err}")
            self.output_errors.append(f"position summary: {calls_err}")
            calls = None
        hidden = self._undetected_rare(calls)      # rare positions nobody carries stay out of both tables
        if hidden:
            print(f"{len(hidden)} rare panel positions were not detected in any sample and are left out of the tables")
        try:
            (final_df.drop(columns=sorted(hidden), level=0) if hidden else final_df).to_csv(
                "./ABO_result.txt", sep="\t", index=False, encoding="utf-8")
            print("Results saved successfully to text file.")
        except OSError as txt_err:
            print(f"ERROR: could not write ABO_result.txt: {txt_err}")
            self.output_errors.append(f"ABO_result.txt: {txt_err}")
            return

        try:
            if calls is None:
                calls = self._position_calls(final_df)
            self._variants_detected = [self._variants_text(c) for c in calls]
            variants_df, pos_labels = self._build_variants_view(final_df, calls, hidden)

            writer = pd.ExcelWriter("./ABO_result.xlsx", engine="xlsxwriter")

            if isinstance(final_df.columns, pd.MultiIndex):
                final_df.columns = final_df.columns.droplevel()

            workbook = writer.book
            data_format = workbook.add_format({"bg_color": "white", "font_color": "black", "border": 1})
            header_format = workbook.add_format({"bold": True, "fg_color": "#007399", "border": 1, "font_color": "white"})
            red_bg_format = workbook.add_format({"bg_color": "#e2725b", "font_color": "black"})
            orange_bg_format = workbook.add_format({"bg_color": "#ff9a00", "font_color": "black"})
            header_format.set_align("center")
            header_format.set_align("vcenter")

            # Sheet 1: Typing (one NGSengine-style row per sample). Sheet order = creation order.
            try:
                typing_df = self._build_typing_df(final_df)
                typing_df.to_excel(writer, sheet_name="Typing", index=False)
                ws_t = writer.sheets["Typing"]
                for c, name in enumerate(typing_df.columns):
                    ws_t.write(0, c, name, header_format)
                int_fmt = workbook.add_format({"border": 1, "num_format": "#,##0"})
                dec_fmt = workbook.add_format({"border": 1, "num_format": "0.0"})
                for r in range(len(typing_df)):
                    for c, name in enumerate(typing_df.columns):
                        v = typing_df.iat[r, c]
                        if v is None or (not isinstance(v, str) and pd.isna(v)):
                            ws_t.write_blank(r + 1, c, None, data_format)
                        elif name == "Mappability %":
                            ws_t.write_number(r + 1, c, float(v), dec_fmt)
                        elif name in ("Number of reads", "Read Length", "Insert Size", "Median Read depth",
                                      "ABO reads", "Mean depth (covered bases)",
                                      "Median analyzed Read depth", "Number of question mark positions",
                                      "BAM-only positions", "Clair3/BAM discordant",
                                      "Discordant primary positions"):
                            ws_t.write_number(r + 1, c, float(v), int_fmt)
                        else:
                            ws_t.write(r + 1, c, v, data_format)
                for c, name in enumerate(typing_df.columns):
                    longest = max([len(str(name))] + [len(str(x)) for x in typing_df[name] if x is not None])
                    ws_t.set_column(c, c, min(max(longest + 2, 10), 80 if name in ("Replicate note", "Variants detected") else 60 if name == "QCFlags" else 40))
                ws_t.set_row(0, 45, header_format)
                ws_t.freeze_panes(1, 0)
                last = len(typing_df)
                if last:
                    qc_t = xl_col_to_name(list(typing_df.columns).index("QCStatus"))
                    ws_t.conditional_format(
                        f"{qc_t}2:{qc_t}{last + 1}",
                        {"type": "cell", "criteria": "==", "value": '"REVIEW"', "format": orange_bg_format},
                    )
                    dc_t = xl_col_to_name(list(typing_df.columns).index("Depth check"))
                    green_bg_format = workbook.add_format({"bg_color": "#c6e0b4", "font_color": "black"})
                    for text_, fmt_ in (("Very low", red_bg_format), ("No data", red_bg_format),
                                        ("Negative control: reads present", red_bg_format),
                                        ("Positive control: depth too low", red_bg_format),
                                        ("Low", orange_bg_format),
                                        ("Negative control", green_bg_format),
                                        ("Positive control", green_bg_format)):
                        ws_t.conditional_format(
                            f"{dc_t}2:{dc_t}{last + 1}",
                            {"type": "text", "criteria": "begins with", "value": text_, "format": fmt_, "stop_if_true": True},
                        )
                    rn_t = xl_col_to_name(list(typing_df.columns).index("Replicate note"))
                    ws_t.conditional_format(
                        f"A2:{xl_col_to_name(len(typing_df.columns) - 1)}{last + 1}",
                        {"type": "formula", "criteria": f"=LEN(${rn_t}2)>0", "format": red_bg_format},
                    )
                    rd_t = xl_col_to_name(list(typing_df.columns).index("Number of reads"))
                    ws_t.conditional_format(
                        f"{rd_t}2:{rd_t}{last + 1}",
                        {"type": "cell", "criteria": "<=", "value": self.min_depth, "format": orange_bg_format},
                    )
            except Exception as typing_err:
                print(f"ERROR: could not write Typing sheet: {typing_err}")
                self.output_errors.append(f"Typing sheet: {typing_err}")

            # Sheet 2: Variants (one column per panel position)
            variants_df.to_excel(writer, sheet_name="Variants", header=False, index=False, startrow=2)
            worksheet = writer.sheets["Variants"]

            about = workbook.add_worksheet("About")
            bold = workbook.add_format({"bold": True})
            about.set_column(0, 0, 34)
            about.set_column(1, 1, 110)
            about.write(0, 0, "ABO typing results", bold)
            about.write(1, 0, RESEARCH_USE_STATEMENT, bold)
            for r_, (k_, v_) in enumerate(self._about_items(), start=3):
                about.write(r_, 0, k_, bold)
                about.write(r_, 1, v_)

            num_rows, num_cols = variants_df.shape
            _cols = list(variants_df.columns)
            reliability_col = xl_col_to_name(_cols.index("ReadReliability"))
            print(f"Data has {num_rows} rows, starting at row 3 with two header rows")

            red_cell = workbook.add_format({"bg_color": "#e2725b", "font_color": "black", "border": 1})
            orange_cell = workbook.add_format({"bg_color": "#ff9a00", "font_color": "black", "border": 1})
            alt_cell = workbook.add_format({"bold": True, "bg_color": "#fff2cc", "font_color": "black", "border": 1})
            grey_cell = workbook.add_format({"bg_color": "#d9d9d9", "font_color": "#595959", "border": 1})
            wrap_header = workbook.add_format({"bold": True, "fg_color": "#007399", "border": 1, "font_color": "white",
                                               "align": "center", "valign": "vcenter", "text_wrap": True})

            for row in range(num_rows):
                for col in range(num_cols):
                    cell_value = variants_df.iat[row, col]
                    if cell_value is None or (not isinstance(cell_value, str) and pd.isna(cell_value)):
                        continue
                    fmt = data_format
                    lab = pos_labels.get(col)
                    if lab is not None:
                        info = calls[row][lab]
                        if info["state"] == "nodata":
                            fmt = grey_cell
                        elif info["reads"] <= self.panel.threshold("reliability_very_low_reads"):
                            fmt = red_cell          # very low depth: the call is not reliable
                        elif info["reads"] <= self.panel.threshold("reliability_low_reads"):
                            fmt = orange_cell       # low depth
                        elif info["state"] == "alt":
                            fmt = alt_cell
                        worksheet.write(row + 2, col, cell_value, fmt)
                        if info["state"] == "alt" and info["type"]:
                            worksheet.write_comment(row + 2, col, info["type"], {"x_scale": 2.2, "y_scale": 0.8})
                    else:
                        worksheet.write(row + 2, col, cell_value, data_format)

            if num_rows:
                last_col = xl_col_to_name(num_cols - 1)
                try:
                    worksheet.conditional_format(
                        f"A3:{last_col}{num_rows + 2}",
                        {"type": "formula", "criteria": f'=${reliability_col}3="{self.reliability_label("very_low")}"', "format": red_bg_format},
                    )
                    worksheet.conditional_format(
                        f"A3:{last_col}{num_rows + 2}",
                        {"type": "formula", "criteria": f'=${reliability_col}3="{self.reliability_label("low")}"', "format": orange_bg_format},
                    )
                except Exception as format_err:
                    print(f"Warning: Could not apply row-level conditional formatting: {format_err}")
                qc_letter = xl_col_to_name(_cols.index("QCStatus"))
                worksheet.conditional_format(
                    f"{qc_letter}3:{qc_letter}{num_rows + 2}",
                    {"type": "cell", "criteria": "==", "value": '"REVIEW"', "format": orange_bg_format},
                )
                rn_letter = xl_col_to_name(_cols.index("ReplicateNote"))
                worksheet.conditional_format(
                    f"A3:{last_col}{num_rows + 2}",
                    {"type": "formula", "criteria": f"=LEN(${rn_letter}3)>0", "format": red_bg_format},
                )

            def _group(first, last, title):
                if first == last:
                    worksheet.write(0, first, title, header_format)
                else:
                    worksheet.merge_range(0, first, 0, last, title, header_format)

            _group(0, 1, "Sample")
            _group(2, 4, "Result")
            pos_idx = sorted(pos_labels)
            run_start = pos_idx[0] if pos_idx else None
            for k, ci in enumerate(pos_idx):
                grp = pos_labels[ci].split("_")[0]
                nxt = pos_labels[pos_idx[k + 1]].split("_")[0] if k + 1 < len(pos_idx) else None
                if nxt != grp:
                    _group(run_start, ci, grp.replace("Exon", "Exon "))
                    run_start = pos_idx[k + 1] if k + 1 < len(pos_idx) else None
            if self.notes_columns:
                first_note = _cols.index(self.notes_columns[0])
                _group(first_note, num_cols - 1, "Notes")
            for col, name in enumerate(_cols):
                worksheet.write(1, col, name.strip(), wrap_header)
            worksheet.set_row(1, 34)
            worksheet.set_column(0, 1, 16)
            worksheet.set_column(2, 4, 18)
            if pos_idx:
                worksheet.set_column(pos_idx[0], pos_idx[-1], 17)
            worksheet.set_column(_cols.index(self.notes_columns[0]), num_cols - 1, 24)
            worksheet.freeze_panes(2, 2)

            writer.close()
            print("Results saved successfully to Excel file.")
        except Exception as excel_err:
            print(f"ERROR: could not write ABO_result.xlsx: {excel_err}")
            self.output_errors.append(f"ABO_result.xlsx: {excel_err}")
            import traceback
            traceback.print_exc()

        # LIS export
        self.df_for_lis_soft = pd.DataFrame()
        self.df_for_lis_soft["Sample ID"] = final_df["Sequencing_ID"]
        self.df_for_lis_soft["Shipment Date"] = ""

        if ("Genotype" in final_df.columns and not final_df["Genotype"].isnull().all()
                and not (final_df["Genotype"] == "Unknown").all()):
            valid_genotype_mask = (final_df["Genotype"] != "Unknown") & final_df["Genotype"].notnull()
            self.df_for_lis_soft.loc[valid_genotype_mask, "ABO Geno Type1"] = final_df.loc[valid_genotype_mask, "Genotype"].str[0]
            self.df_for_lis_soft.loc[valid_genotype_mask, "ABO Geno Type2"] = final_df.loc[valid_genotype_mask, "Genotype"].str[1]
        else:
            self.df_for_lis_soft["ABO Geno Type1"] = ""
            self.df_for_lis_soft["ABO Geno Type2"] = ""

        plain_phenotype = final_df["Phenotype"].copy()

        if "BwSubtype" in final_df.columns:
            bw_mask = final_df["BwSubtype"].astype(str).str.len() > 0
            annotated_phenotype = final_df["Phenotype"].copy()
            annotated_phenotype.loc[bw_mask] = (
                final_df.loc[bw_mask, "Phenotype"] + " [" + final_df.loc[bw_mask, "BwSubtype"] + "]"
            )
            final_df["Phenotype"] = annotated_phenotype

        self.df_for_lis_soft["ABO Pheno Type"] = plain_phenotype
        self.df_for_lis_soft["RH"] = ""
        self.df_for_lis_soft["Blood Type"] = plain_phenotype
        self.df_for_lis_soft["ABORH Comments"] = ""

        if "ASubtype" in final_df.columns:
            self.df_for_lis_soft["A Subtype"] = final_df["ASubtype"]
        if "BwSubtype" in final_df.columns:
            self.df_for_lis_soft["Bw Subtype"] = final_df["BwSubtype"]
        if "PhaseConfidence" in final_df.columns:
            self.df_for_lis_soft["Phase Confidence"] = final_df["PhaseConfidence"]
        if "ReadReliability" in final_df.columns:
            self.df_for_lis_soft["Read Reliability"] = final_df["ReadReliability"]

        if isinstance(final_df.columns, pd.MultiIndex):
            reads_df = final_df.loc[:, (slice(None), "#Reads")]
            self.df_for_lis_soft["#Reads"] = reads_df.mean(axis=1)
        else:
            reads_columns = [col for col in final_df.columns if "#Reads" in str(col)]
            if reads_columns:
                self.df_for_lis_soft["#Reads"] = final_df[reads_columns].mean(axis=1)
            else:
                self.df_for_lis_soft["#Reads"] = 0

        for _src, _dst in (("QCStatus", "QC Status"), ("QCFlags", "QC Flags"),
                           ("A2Trace", "A2 Trace"), ("EvidenceBasis", "Evidence Basis"),
                           ("ReplicateNote", "Replicate Note")):
            if _src in final_df.columns:
                self.df_for_lis_soft[_dst] = final_df[_src]

        self.df_for_lis_soft.drop_duplicates(inplace=True)
        self.df_for_lis_soft.to_csv("./final_export.csv", index=False, encoding="utf-8")
        print(f"LIS export file created successfully with {len(self.df_for_lis_soft)} samples")

    def _about_items(self):
        """Version / provenance lines for the About sheet and the QC report."""
        meta = getattr(self.panel, "meta", {}) or {}
        items = [("aggregate_abo_reports.py", getattr(self, "script_version", "")),
                 ("Panel version", str(meta.get("panel_version", "unknown"))),
                 ("Panel source", str(meta.get("source", "unknown")))]
        if meta.get("isbt_allele_table_version"):
            items.append(("ISBT ABO allele table", str(meta["isbt_allele_table_version"])))
        if meta.get("isbt_allele_table_checked"):
            items.append(("ISBT table checked", str(meta["isbt_allele_table_checked"])))
        if meta.get("isbt_allele_table_release_date"):
            items.append(("ISBT table release date", str(meta["isbt_allele_table_release_date"])))
        return items

    def _evidence_summary(self, final_df):
        """Per-sample evidence rows (None when no sample carries an evidence status)."""
        rows = []
        for _, r in final_df.iterrows():
            try:
                key = (str(r["Sequencing_ID"]), int(r["Barcode"]))
            except (TypeError, ValueError):
                continue
            m = self.sample_metrics.get(key, {})
            if not m.get("evidence_source"):
                continue
            rows.append({
                "Sample": key[0],
                "Barcode": f"barcode{key[1]:02d}",
                "Evidence source": m["evidence_source"],
                "Evidence status": m.get("evidence_label") or "OK",
                "BAM-only positions": m.get("bam_only"),
                "Clair3/BAM discordant": m.get("discordant"),
                "Discordant primary positions": m.get("discordant_primary"),
                "QC status": r.get("QCStatus", ""),
                "Typing result": r.get("ExtendedGenotype", ""),
            })
        return rows or None

    def write_evidence_summary(self, final_df, path="ABO_evidence_summary_mqc.tsv"):
        """MultiQC custom-content table: where each sample's evidence came from, and which
        samples had failed or empty Clair3 output, no reads over the ABO region, or disagreement
        between Clair3 and the BAM pileup."""
        rows = self._evidence_summary(final_df)
        if not rows:
            return
        cols = list(rows[0].keys())
        try:
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("# id: 'abo_evidence_summary'\n")
                fh.write("# section_name: 'ABO evidence summary'\n")
                fh.write("# description: 'Agreement between the two lines of evidence (Clair3 and BAM pileup) "
                         "per sample. Samples with no reads over the ABO region, Clair3-only or BAM-only "
                         "evidence, or discordant primary positions need review. "
                         + RESEARCH_USE_STATEMENT.replace("'", "") + "'\n")
                fh.write("# plot_type: 'table'\n")
                fh.write("\t".join(cols) + "\n")
                for row in rows:
                    fh.write("\t".join("" if row[c] is None else str(row[c]) for c in cols) + "\n")
            print(f"Evidence summary written to {path}")
        except OSError as e:
            print(f"Could not write evidence summary: {e}")

    def report_qc(self, final_df):
        """Print (and optionally write) a run-level QC summary. Called after
        the results are saved, so final_df has single-level column names."""
        lines = ["=== QC SUMMARY ===", RESEARCH_USE_STATEMENT]
        lines += [f"{k}: {v}" for k, v in self._about_items()]
        n = len(final_df)
        if n == 0:
            return
        review = int((final_df["QCStatus"] == "REVIEW").sum()) if "QCStatus" in final_df.columns else 0
        lines.append(f"Samples: {n}   flagged REVIEW: {review} ({100.0 * review / n:.1f}%)")

        if "QCFlags" in final_df.columns:
            tally = Counter()
            for flags in final_df["QCFlags"].astype(str):
                for f in flags.split("; "):
                    f = f.strip()
                    if not f or f == "nan":
                        continue
                    f = re.sub(r"\(barcode.*?\)$", "", f).strip()     # drop per-group member lists
                    f = re.sub(r"\(\d+ < \d+,", "(<N,", f)           # collapse read counts
                    f = re.sub(r"\(O1 x\d, O2 x\d, B x\d\)", "", f).strip()
                    tally[f] += 1
            if tally:
                lines.append("Flags (a sample can carry more than one):")
                for f, c in tally.most_common():
                    lines.append(f"  {c:6d}  {f}")

        ev_rows = self._evidence_summary(final_df)
        if ev_rows:
            lines.append("Evidence (Clair3 + BAM pileup):")
            for k, v in Counter(r["Evidence source"] for r in ev_rows).most_common():
                lines.append(f"  {v:6d}  {k}")
            problems = [r for r in ev_rows if r["Evidence status"] != "OK" or r["Discordant primary positions"]]
            if problems:
                lines.append("Samples needing attention:")
                for r in problems:
                    extra = (f"; discordant at {r['Discordant primary positions']} primary position(s)"
                             if r["Discordant primary positions"] else "")
                    lines.append(f"  {r['Sample']} ({r['Barcode']}): {r['Evidence status']}{extra}")

        if "EvidenceBasis" in final_df.columns:
            lines.append("Evidence basis: " + ", ".join(
                f"{k}={v}" for k, v in final_df["EvidenceBasis"].value_counts().items()))
        if "A2Trace" in final_df.columns:
            lines.append("A2 trace:")
            for k, v in final_df["A2Trace"].value_counts().sort_index().items():
                lines.append(f"  {v:6d}  {k}")
        rs = self.replicate_stats
        if rs["groups"]:
            bad = rs["discordant_genotype"] + rs["discordant_a2"] + rs.get("discordant_allele", 0)
            lines.append(f"Replicate groups: {rs['groups']}   disagreeing: {bad} "
                         f"(genotype {rs['discordant_genotype']}, A2 evidence {rs['discordant_a2']}, "
                         f"allele level {rs.get('discordant_allele', 0)})")

        notes = []
        recs = self.qc_records
        if recs:
            uni = sum(1 for r in recs if r["reads_uniform"]) / len(recs)
            if uni >= 0.9:
                notes.append(
                    f"#Reads is identical at every primary position for {100 * uni:.0f}% of samples, i.e. it is a "
                    "sample-wide read count, not depth at the site. ReadReliability and the low-read flag therefore "
                    "cannot confirm site-level evidence. Emit per-position depth upstream (AlignmentStatistics 'Depth') "
                    "and write it into the report's 'Aligned Read Count'.")
            ref_only = sum(1 for r in recs if r["evidence"] == "reference-only") / len(recs)
            if ref_only >= 0.5:
                notes.append(
                    f"{100 * ref_only:.0f}% of samples are reference-only (no variant at any primary position). "
                    "These cannot be told apart from failed calls; confirm a subset with independent typing.")
            c = Counter(r["s261"] for r in recs if r["s261"] is not None)
            tot = sum(c.values())
            if tot >= 200:
                q = (c["het"] + 2 * c["hom"]) / (2 * tot)
                exp_het = 2 * q * (1 - q) * tot
                if exp_het > 0:
                    F = 1 - c["het"] / exp_het
                    if F > 0.30:
                        notes.append(
                            f"c.261delG shows a heterozygote deficit (observed {c['het']}, expected {exp_het:.0f} under "
                            f"Hardy-Weinberg; F={F:.2f}). Possible under selection or relatedness; also what missed "
                            "heterozygous calls look like.")
            aa = sum(1 for r in recs if r["genotype"] == "AA")
            ao = sum(1 for r in recs if r["genotype"] == "AO")
            if aa + ao >= 100 and aa / (aa + ao) > 0.5:
                notes.append(
                    f"A/A is {100 * aa / (aa + ao):.0f}% of group A (A/O {100 * ao / (aa + ao):.0f}%). In most populations "
                    "A/O carriers outnumber A/A; check against independently typed samples.")
        if notes:
            lines.append("Notes (informational):")
            for t in notes:
                lines.append("  - " + t)

        text = "\n".join(lines)
        print("\n" + text)
        if self.qc_report_path:
            try:
                with open(self.qc_report_path, "w", encoding="utf-8") as fh:
                    fh.write(text + "\n")
                print(f"QC report written to {self.qc_report_path}")
            except OSError as e:
                print(f"Could not write QC report: {e}")
