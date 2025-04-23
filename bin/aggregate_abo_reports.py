#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import glob
import pandas as pd
from xlsxwriter.utility import xl_col_to_name

__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2024, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = [
    "Fredrick Mobegi",
    "Benedict Matern",
    "Mathijs Groeneweg",
    "Claude 3.7 Sonnet Thinking (for rewrite to add ABO*A1/2/3 subtypes)",
]
__license__ = "GPL"
__version__ = "0.2.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Development"

"""
This file is part of the nf-core/abotyper pipeline "https://github.com/fmobegi/nf-core-abotyper".

nf-core/abotyper is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This pipeline is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with nf-core/abotyper. If not, see <http://www.gnu.org/licenses/>.


This class collates all ABO phenotype results from each sample into an a general table 
and generates an Excel worksheet and a CSV file for export to LIS soft or 
other general purpose lab management systems.
"""

print(
    "\033[92m\n ********* Started combining samples to single file ********* \033[0m\n"
)


class ABOReportParser:

    def __init__(self, input_dir):
        """
        Initialize the ABOReportParser.

        Args:
            input_dir (str): The input directory containing data files.
        """
        self.input_dir = input_dir
        self.results = []
        self.initialize_columns()
        self.failed_samples = []

    def initialize_columns(self):
        """
        Define the column headers for all exon positions.
        """
        # Primary ABO typing variant in exon 6
        exon6 = ["Exon6_pos22"] * 10  # c.261

        # ABO*A2 subtype positions in exon 6
        exon6_27 = ["Exon6_pos27"] * 10  # c.266
        exon6_29 = ["Exon6_pos29"] * 10  # c.268
        exon6_58 = ["Exon6_pos58"] * 10  # c.297

        # Primary variants in exon 7
        exon7_422 = ["Exon7_pos422"] * 10  # c.796C>A
        exon7_428 = ["Exon7_pos428"] * 10  # c.802G>A / c.802G>C
        exon7_429 = ["Exon7_pos429"] * 10  # c.803G>C / c.803G>T
        exon7_431 = ["Exon7_pos431"] * 10  # c.804dupG

        # A subtype positions in exon 7
        exon7_93 = ["Exon7_pos93"] * 10  # c.467C>T
        exon7_165 = ["Exon7_pos165"] * 10  # c.539G>C / c.539G>A
        exon7_272 = ["Exon7_pos272"] * 10  # c.646
        exon7_307 = ["Exon7_pos307"] * 10  # c.681
        exon7_371 = ["Exon7_pos371"] * 10  # c.745C>T
        exon7_446 = ["Exon7_pos446"] * 10  # c.820G>A
        exon7_680 = ["Exon7_pos680"] * 10  # c.1054C>T
        exon7_687 = ["Exon7_pos687"] * 10  # c.1061delC

        header_cols = (
            ["", ""]
            + exon6
            + exon6_27
            + exon6_29
            + exon6_58
            + exon7_422
            + exon7_428
            + exon7_429
            + exon7_431
            + exon7_93
            + exon7_165
            + exon7_272
            + exon7_307
            + exon7_371
            + exon7_446
            + exon7_680
            + exon7_687
            + ["", "", "", ""]
        )

        column_metrics = [
            "#Reads",
            "Mat",
            "Mis",
            "Ins",
            "Del",
            "A",
            "G",
            "C",
            "T",
            "Type",
        ]
        header_rows = (
            ["Barcode", "Sequencing_ID"]
            + column_metrics * 16
            + ["Phenotype", "Genotype", "ExtendedGenotype", "Reliability"]
        )

        self.columns = pd.MultiIndex.from_arrays([header_cols, header_rows])

    def parse_exon7(self, filename):
        """
        Open the file for reading and processing all exon 7 positions!
        """
        try:
            with open(filename, "r", encoding="utf-8") as f:
                lines = f.readlines()

            positions = []
            counts = []
            mat_values = []
            mis_values = []
            ins_values = []
            del_values = []
            a_values = []
            g_values = []
            c_values = []
            t_values = []

            i = 0
            while i < len(lines):
                line = lines[i].strip()

                if "Exon 7 position(1-based):" in line:
                    pos_match = re.search(r":\s*(\d+)", line)
                    if pos_match:
                        pos = int(pos_match.group(1))

                        for j in range(i, min(i + 10, len(lines))):
                            if "Aligned Read Count:" in lines[j]:
                                count_match = re.search(r":\s*(\d+)", lines[j])
                                if count_match:
                                    count = int(count_match.group(1))

                                    stats_idx = j + 2
                                    if (
                                        stats_idx < len(lines)
                                        and "Mat" in lines[stats_idx - 1]
                                    ):
                                        stats = lines[stats_idx].split()
                                        if len(stats) >= 8:
                                            mat, mis, ins, dele, a, g, c, t = [
                                                float(x) for x in stats
                                            ]

                                            positions.append(pos)
                                            counts.append(count)
                                            mat_values.append(mat)
                                            mis_values.append(mis)
                                            ins_values.append(ins)
                                            del_values.append(dele)
                                            a_values.append(a)
                                            g_values.append(g)
                                            c_values.append(c)
                                            t_values.append(t)

                                            break
                i += 1

            df = pd.DataFrame(
                {
                    "Exon": ["7"] * len(positions),
                    "Position": positions,
                    "#Reads": counts,
                    "Mat": mat_values,
                    "Mis": mis_values,
                    "Ins": ins_values,
                    "Del": del_values,
                    "A": a_values,
                    "G": g_values,
                    "C": c_values,
                    "T": t_values,
                }
            )

            # Make sure all needed positions are in the DataFrame
            # Primary positions + all subtype positions
            # pos. exon7= exonic(CDS) 422(796),428(802),429(803),431(805),93(467),165(539),687(1061)
            # Additional exonic(CDS) 272(646), 307(681), 371(745), 446(820), 680(1054)
            # [93, 165, 272, 307, 371, 446, 680, 687]
            all_positions = [
                422,
                428,
                429,
                431,  # Primary
                93,
                165,
                272,
                307,
                371,
                446,
                680,
                687,  # A1/A2/A3 subtypes
            ]

            for pos in all_positions:
                if pos not in df["Position"].values:
                    df = pd.concat(
                        [
                            df,
                            pd.DataFrame(
                                {
                                    "Exon": ["7"],
                                    "Position": [pos],
                                    "#Reads": [0],
                                    "Mat": [0],
                                    "Mis": [0],
                                    "Ins": [0],
                                    "Del": [0],
                                    "A": [0],
                                    "G": [0],
                                    "C": [0],
                                    "T": [0],
                                }
                            ),
                        ],
                        ignore_index=True,
                    )

            df = df.sort_values("Position").reset_index(drop=True)

            df["Type"] = df.apply(
                lambda row: self.get_type(
                    row["Position"], row["A"], row["G"], row["C"], row["T"], row["Del"]
                ),
                axis=1,
            )

            return df

        except Exception as e:
            print(f"Error parsing exon 7 file {filename}: {str(e)}")
            empty_df = pd.DataFrame(
                columns=[
                    "Exon",
                    "Position",
                    "#Reads",
                    "Mat",
                    "Mis",
                    "Ins",
                    "Del",
                    "A",
                    "G",
                    "C",
                    "T",
                    "Type",
                ]
            )

            for pos in [
                422,
                428,
                429,
                431,
                93,
                165,
                272,
                307,
                371,
                446,
                680,
                687,
            ]:
                empty_df = pd.concat(
                    [
                        empty_df,
                        pd.DataFrame(
                            {
                                "Exon": ["7"],
                                "Position": [pos],
                                "#Reads": [0],
                                "Mat": [0],
                                "Mis": [0],
                                "Ins": [0],
                                "Del": [0],
                                "A": [0],
                                "G": [0],
                                "C": [0],
                                "T": [0],
                                "Type": [""],
                            }
                        ),
                    ],
                    ignore_index=True,
                )

            return empty_df

    def parse_exon6(self, filename):
        """Parse exon 6 and extract data for all relevant positions (22, 27, 29, 58)."""
        try:
            with open(filename, "r", encoding="utf-8") as f:
                lines = f.readlines()

            positions = []
            counts = []
            mat_values = []
            mis_values = []
            ins_values = []
            del_values = []
            a_values = []
            g_values = []
            c_values = []
            t_values = []

            i = 0
            while i < len(lines):
                line = lines[i].strip()

                if "Exon 6 position(1-based):" in line:
                    pos_match = re.search(r":\s*(\d+)", line)
                    if pos_match:
                        pos = int(pos_match.group(1))

                        for j in range(i, min(i + 10, len(lines))):
                            if "Aligned Read Count:" in lines[j]:
                                count_match = re.search(r":\s*(\d+)", lines[j])
                                if count_match:
                                    count = int(count_match.group(1))

                                    stats_idx = j + 2
                                    if (
                                        stats_idx < len(lines)
                                        and "Mat" in lines[stats_idx - 1]
                                    ):
                                        stats = lines[stats_idx].split()
                                        if len(stats) >= 8:
                                            mat, mis, ins, dele, a, g, c, t = [
                                                float(x) for x in stats
                                            ]

                                            positions.append(pos)
                                            counts.append(count)
                                            mat_values.append(mat)
                                            mis_values.append(mis)
                                            ins_values.append(ins)
                                            del_values.append(dele)
                                            a_values.append(a)
                                            g_values.append(g)
                                            c_values.append(c)
                                            t_values.append(t)

                                            break
                i += 1

            df = pd.DataFrame(
                {
                    "Exon": ["6"] * len(positions),
                    "Position": positions,
                    "#Reads": counts,
                    "Mat": mat_values,
                    "Mis": mis_values,
                    "Ins": ins_values,
                    "Del": del_values,
                    "A": a_values,
                    "G": g_values,
                    "C": c_values,
                    "T": t_values,
                }
            )

            all_positions = [22, 27, 29, 58]

            for pos in all_positions:
                if pos not in df["Position"].values:
                    df = pd.concat(
                        [
                            df,
                            pd.DataFrame(
                                {
                                    "Exon": ["6"],
                                    "Position": [pos],
                                    "#Reads": [0],
                                    "Mat": [0],
                                    "Mis": [0],
                                    "Ins": [0],
                                    "Del": [0],
                                    "A": [0],
                                    "G": [0],
                                    "C": [0],
                                    "T": [0],
                                }
                            ),
                        ],
                        ignore_index=True,
                    )

            df = df.sort_values("Position").reset_index(drop=True)

            df["Type"] = df.apply(
                lambda row: self.get_type_exon6(
                    row["Position"], row["A"], row["G"], row["C"], row["T"], row["Del"]
                ),
                axis=1,
            )

            return df

        except Exception as e:
            print(f"Error parsing exon 6: {str(e)}")
            empty_df = pd.DataFrame(
                columns=[
                    "Exon",
                    "Position",
                    "#Reads",
                    "Mat",
                    "Mis",
                    "Ins",
                    "Del",
                    "A",
                    "G",
                    "C",
                    "T",
                    "Type",
                ]
            )

            for pos in [22, 27, 29, 58]:
                empty_df = pd.concat(
                    [
                        empty_df,
                        pd.DataFrame(
                            {
                                "Exon": ["6"],
                                "Position": [pos],
                                "#Reads": [0],
                                "Mat": [0],
                                "Mis": [0],
                                "Ins": [0],
                                "Del": [0],
                                "A": [0],
                                "G": [0],
                                "C": [0],
                                "T": [0],
                                "Type": [""],
                            }
                        ),
                    ],
                    ignore_index=True,
                )

            return empty_df

    def get_type_exon6(self, pos, a, g, c, t, dele):
        """
        Determine blood type or subtype for each exon 6 position based on nucleotide percentages.

        Args:
            pos: The position number
            a, g, c, t, dele: Nucleotide and deletion percentages
        """
        if pos == 22:  # c.261
            if g >= 80 and g > dele:
                return "A or B or O"
            elif dele >= 80 and dele > g:
                return "O1"
            elif abs(g + dele) >= 20:
                return "O1 and (A or B or O)"
            elif 20 < dele < 80 and 20 < g < 80:
                return "O1 and (A or B or O)"

        elif pos == 27:  # c.266
            if c >= 80:
                return "A1 or A3"
            elif t >= 80:
                return "A2"
            elif 20 < c < 80 and 20 < t < 80:
                return "A"

        elif pos == 29:  # c.268
            if t >= 80:
                return "A1 or A3"
            elif c >= 80:
                return "A2"
            elif 20 < t < 80 and 20 < c < 80:
                return "A"

        elif pos == 58:  # c.297
            if a >= 80:
                return "A1 or A3"
            elif g >= 80:
                return "A2"
            elif 20 < a < 80 and 20 < g < 80:
                return "A"

        return ""

    def get_type(self, pos, a, g, c, t, dele):
        """
        Determine the blood type or subtype for each position based on nucleotide percentages.
        This is used for first-pass analysis regardless of phenotype.

        Args:
            pos: The position number
            a, g, c, t: Nucleotide percentages
        """
        # Primary ABO variants
        if pos == 422:
            if c >= 80:
                return "A or O"
            elif a >= 80:
                return "B"
            elif abs(a - c) <= 20:
                return "(A or O) and B"
            elif 15 < a < 80 and 15 < c < 80:
                return "(A or O) and B"
        elif pos == 428:
            if g >= 70:
                return "O and (A or B)"
            elif a >= 70:
                return "O2"
            elif abs(g - a) <= 20:
                return "O2 and (O or A or B)"
            elif 15 < g < 70 and 15 < a < 70:
                return "O2 and (O or A or B)"
        elif pos == 429:
            if g >= 80:
                return "A or O"
            elif c >= 80:
                return "B"
            elif abs(g - c) <= 20:
                return "(A or O) and B"
            elif 15 < g < 80 and 20 < c < 80:
                return "(A or O) and B"
        elif pos == 431:
            if t >= 80:
                return "O and (A or B)"
            elif g >= 80:
                return "O3"
            elif abs(t - g) <= 20:
                return "O3 and (O or A or B)"
            elif 15 < t < 80 and 15 < g < 80:
                return "O3 and (O or A or B)"

        # A subtype positions
        elif pos == 93:  # genomic pos 467 /A3
            if c >= 80:
                return "A1"
            elif t >= 80:
                return "A2 or A3"
            # elif 20 < c < 80 and 20 < t < 80:
            #     return "A"
        elif pos == 165:  # genomic pos 539
            if c >= 80:
                return "A1 or A2"
            elif t >= 80:
                return "A3"
            # elif 20 < c < 80 and 20 < t < 80:
            #     return "A"
        elif pos == 272:  # genomic pos 646
            if t >= 80:
                return "A1"
            elif a >= 80:
                return "A2"
            # elif 20 < t < 80 and 20 < a < 80:
            #     return "A1 or A2"
        elif pos == 307:  # genomic pos 681
            if g >= 80:
                return "A1 or A2"
            elif a >= 80:
                return "A3"
            # elif 20 < g < 80 and 20 < a < 80:
            #     return "A"
        elif pos == 371:  # genomic pos 745
            if c >= 80:
                return "A1 or A2"
            elif t >= 80:
                return "A3"
            # elif 20 < c < 80 and 20 < t < 80:
            #     return "A"
        elif pos == 446:  # genomic pos 820
            if a >= 80:
                return "A1 or A2"
            elif c >= 80:
                return "A3"
            # elif 20 < a < 80 and 20 < c < 80:
            #     return "A"
        elif pos == 680:  # genomic pos 1054
            if g >= 80:
                return "A1 or A3"
            elif a >= 80:
                return "A2"
            elif 20 < g < 80 and 20 < a < 80:
                return "A"
        elif pos == 687:  # genomic pos 1061 /A3
            if c >= 80:
                return "A1"
            elif dele >= 80:  # Deletion indicates A2 or A3 subtypes
                return "A2 or A3"
            # elif 20 < c < 80 and 20 < dele < 80:
            #     return "A"

        return ""

    def assign_phenotype_genotype(self, df):
        """Assign the phenotype and genotype information"""
        try:
            # Primary positions
            type_exon6 = df.at[0, ("Exon6_pos22", "Type")]
            type_exon7_422 = df.at[0, ("Exon7_pos422", "Type")]
            type_exon7_428 = df.at[0, ("Exon7_pos428", "Type")]
            type_exon7_429 = df.at[0, ("Exon7_pos429", "Type")]
            type_exon7_431 = df.at[0, ("Exon7_pos431", "Type")]

            # Exon 6 A subtype positions
            type_exon6_27 = df.at[0, ("Exon6_pos27", "Type")]
            type_exon6_29 = df.at[0, ("Exon6_pos29", "Type")]
            type_exon6_58 = df.at[0, ("Exon6_pos58", "Type")]

            # Exon 7 A subtype positions
            type_exon7_93 = df.at[0, ("Exon7_pos93", "Type")]
            type_exon7_165 = df.at[0, ("Exon7_pos165", "Type")]
            type_exon7_272 = df.at[0, ("Exon7_pos272", "Type")]
            type_exon7_307 = df.at[0, ("Exon7_pos307", "Type")]
            type_exon7_371 = df.at[0, ("Exon7_pos371", "Type")]
            type_exon7_446 = df.at[0, ("Exon7_pos446", "Type")]
            type_exon7_680 = df.at[0, ("Exon7_pos680", "Type")]
            type_exon7_687 = df.at[0, ("Exon7_pos687", "Type")]

            nreads6 = df.at[0, ("Exon6_pos22", "#Reads")]
            nreads_exon7_p422 = df.at[0, ("Exon7_pos422", "#Reads")]
            nreads_exon7_p428 = df.at[0, ("Exon7_pos428", "#Reads")]
            nreads_exon7_p429 = df.at[0, ("Exon7_pos429", "#Reads")]
            nreads_exon7_p431 = df.at[0, ("Exon7_pos431", "#Reads")]

            Phenotype = "Unknown"
            Genotype = "Unknown"
            ExtendedGenotype = "Unknown"

            # PART 1: PRIMARY PHENOTYPING LOGIC

            ## OA COMBINATIONS ---------------------------------------------------------------------------
            ## combination 1 | AO1 --
            if (
                (type_exon6 == "O1 and (A or B or O)")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "A"
                Genotype = "AO"
                ExtendedGenotype = "AO1"

            ## combination 2 | AO2 --
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O2 and (O or A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "A"
                Genotype = "AO"
                ExtendedGenotype = "AO2"
                # Reliability = 'Enter-manually'

            ## combination 3 | AO3 --
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O3 and (O or A or B)")
            ):
                Phenotype = "A"
                Genotype = "AO"
                ExtendedGenotype = "AO3"

            ## OB COMBINATIONS ---------------------------------------------------------------------------
            ## combination 4 | BO1 --
            elif (
                (type_exon6 == "O1 and (A or B or O)")
                and (type_exon7_422 == "(A or O) and B")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "(A or O) and B")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "BO1"
                # Reliability = 'Enter-manually'

            ## combination 5 | O2B --
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "(A or O) and B")
                and (type_exon7_428 == "O2 and (O or A or B)")
                and (type_exon7_429 == "(A or O) and B")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "O2B"
                # Reliability = 'Enter-manually'

            ## combination 6 | AO3 --
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "(A or O) and B")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "(A or O) and B")
                and (type_exon7_431 == "O3 and (O or A or B)")
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "BO3"

            ## OO COMBINATIONS  ---------------------------------------------------------------------------
            ## combination 7 | O1O2 --
            elif (
                (type_exon6 == "O1 and (A or B or O)")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O2 and (O or A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1O2"

            ## combination 8 | O1O3 --
            elif (
                (type_exon6 == "O1 and (A or B or O)")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O3 and (O or A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1O3"
                # Reliability = 'Enter-manually'

            ## combination 9 | O2O3 --
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O2 and (O or A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O3 and (O or A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O2O3"

            ## combination 10 | O1O1 --
            elif (
                (type_exon6 == "O1")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1O1"
                # Reliability = 'Enter-manually'

            ## combination 11 | O2O2 --
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O2")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O2O2"
                # Reliability = 'Enter-manually'
                #
            ## combination 12 | O3O3 --
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O3")
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O3O3"

            ## combination 13 | AA ---------------------------------------------------------------------------
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "A or O")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "A or O")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "A"
                Genotype = "AA"
                ExtendedGenotype = "AA"

            ## combination 14 | BB ---------------------------------------------------------------------------
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "B")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "B")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "B"
                Genotype = "BB"
                ExtendedGenotype = "BB"

            ## combination 15 | AB ---------------------------------------------------------------------------
            elif (
                (type_exon6 == "A or B or O")
                and (type_exon7_422 == "(A or O) and B")
                and (type_exon7_428 == "O and (A or B)")
                and (type_exon7_429 == "(A or O) and B")
                and (type_exon7_431 == "O and (A or B)")
            ):
                Phenotype = "AB"
                Genotype = "AB"
                ExtendedGenotype = "AB"

            ## TODO extend to capture ABO*A subtypes (A1, A2, and A3)

            ## UNKNOWN None of the above --------------------------------------------------------------------
            else:
                Phenotype = "Unknown"
                Genotype = "Unknown"
                ExtendedGenotype = "Unknown"

            read_counts = [
                nreads6,
                nreads_exon7_p422,
                nreads_exon7_p428,
                nreads_exon7_p429,
                nreads_exon7_p431,
            ]
            read_counts = [x for x in read_counts if str(x) != "nan" and float(x) > 0]

            if read_counts:
                min_reads = min(read_counts)
                if min_reads <= 20:
                    Reliability = "Very Low(\u226420 reads)"
                elif 20 < min_reads <= 40:
                    Reliability = "Low (\u226440 reads)"
                elif min_reads >= 500:
                    Reliability = "Robust(\u2265500 reads)"
                else:
                    Reliability = "Normal"
            else:
                Reliability = "Unknown (no read data)"

            df[("", "Phenotype")] = Phenotype
            df[("", "Genotype")] = Genotype
            df[("", "ExtendedGenotype")] = ExtendedGenotype
            df[("", "Reliability")] = Reliability

            return df

        except Exception as e:
            print(f"Error in assign_phenotype_genotype: {str(e)}")
            import traceback

            traceback.print_exc()
            df[("", "Phenotype")] = "Error"
            df[("", "Genotype")] = "Error"
            df[("", "ExtendedGenotype")] = "Error"
            df[("", "Reliability")] = "Error processing"
            return df

    def process_file(self, filename):
        """Process a single file and extract all necessary data."""
        try:
            if "_" in filename:
                parts = filename.split("_")
                sample_name = "_".join(parts[:-1])
                barcode = parts[-1]
            else:
                sample_name = filename
                barcode = ""

            exon6_dir = os.path.join(self.input_dir, filename, "exon6")
            exon7_dir = os.path.join(self.input_dir, filename, "exon7")

            if not (os.path.exists(exon6_dir) and os.path.exists(exon7_dir)):
                print(f"Skipping file {filename}. Missing exon6 or exon7 directory.")
                return

            exon6_phenotypes = os.path.join(exon6_dir, "*.ABOPhenotype.txt")
            exon7_phenotypes = os.path.join(exon7_dir, "*.ABOPhenotype.txt")

            exon6_phenotype_files = glob.glob(exon6_phenotypes)
            exon7_phenotype_files = glob.glob(exon7_phenotypes)

            if not exon6_phenotype_files:
                print(f"Missing exon6 phenotype files for {filename}. Skipping.")
                self.failed_samples.append(
                    {"sample": filename, "reason": "Missing exon6 phenotype files"}
                )
                return

            if not exon7_phenotype_files:
                print(f"Missing exon7 phenotype files for {filename}. Skipping.")
                self.failed_samples.append(
                    {"sample": filename, "reason": "Missing exon7 phenotype files"}
                )
                return

            if os.path.getsize(exon6_phenotype_files[0]) == 0:
                print(f"Empty exon6 phenotype file (0 kb) for {filename}. Skipping.")
                self.failed_samples.append(
                    {"sample": filename, "reason": "Empty exon6 phenotype file (0 kb)"}
                )
                return

            if os.path.getsize(exon7_phenotype_files[0]) == 0:
                print(f"Empty exon7 phenotype file (0 kb) for {filename}. Skipping.")
                self.failed_samples.append(
                    {"sample": filename, "reason": "Empty exon7 phenotype file (0 kb)"}
                )
                return

            result_df = pd.DataFrame(columns=self.columns)
            result_df.loc[0, ("", "Barcode")] = barcode.replace("barcode", "")
            result_df.loc[0, ("", "Sequencing_ID")] = sample_name

            exon6_data = self.parse_exon6(exon6_phenotype_files[0])
            if not exon6_data.empty:
                for pos in [22, 27, 29, 58]:
                    pos_df = exon6_data[exon6_data["Position"] == pos]
                    if not pos_df.empty:
                        for col in [
                            "#Reads",
                            "Mat",
                            "Mis",
                            "Ins",
                            "Del",
                            "A",
                            "G",
                            "C",
                            "T",
                            "Type",
                        ]:
                            if col in pos_df.columns:
                                result_df.loc[0, (f"Exon6_pos{pos}", col)] = (
                                    pos_df.iloc[0][col]
                                )

            exon7_data = self.parse_exon7(exon7_phenotype_files[0])
            if not exon7_data.empty:
                all_positions = [
                    422,
                    428,
                    429,
                    431,
                    93,
                    165,
                    272,
                    307,
                    371,
                    446,
                    680,
                    687,
                ]

                for pos in all_positions:
                    pos_df = exon7_data[exon7_data["Position"] == pos]
                    if not pos_df.empty:
                        for col in [
                            "#Reads",
                            "Mat",
                            "Mis",
                            "Ins",
                            "Del",
                            "A",
                            "G",
                            "C",
                            "T",
                            "Type",
                        ]:
                            if col in pos_df.columns:
                                result_df.loc[0, (f"Exon7_pos{pos}", col)] = (
                                    pos_df.iloc[0][col]
                                )

            result_df = self.assign_phenotype_genotype(result_df)

            self.results.append(result_df)

            print(f"Successfully processed {filename}")

        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
            import traceback

            traceback.print_exc()

    def process_files(self):
        """Process all files in the input directory that match expected patterns."""
        for filename in os.listdir(self.input_dir):
            if os.path.isdir(os.path.join(self.input_dir, filename)):
                try:
                    pattern = r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[A-Z0-9]+)?(-[A-Z0-9]+)?_barcode\d+$"
                    match = re.match(pattern, filename)

                    if match:
                        print("\nProcessing file: " + filename)

                        parts = filename.split("_")

                        if len(parts) >= 2:
                            sample_name = parts[0]
                            barcode = parts[-1]
                            if barcode.startswith("barcode"):
                                print(
                                    f"Extracted Sample: {sample_name}, Barcode: {barcode}"
                                )
                                self.process_file(filename)
                                print(
                                    "Done adding Sample %s with barcode %s to merged data frame"
                                    % (sample_name, barcode)
                                )
                            else:
                                print(
                                    f"\nBarcode format incorrect in filename: {filename}"
                                )
                        else:
                            print(
                                f"Filename does not have the expected number of parts: {filename}"
                            )
                    else:
                        print(f"\nFile does not match expected patterns: {filename}")
                except Exception as e:
                    print(f"\nError processing file {filename}: {e}")
                finally:
                    print(f"Finished processing file: {filename}")

    def merge_dataframes(self):
        final_df = pd.concat(self.results)

        final_df[("", "Barcode")] = final_df[("", "Barcode")].astype(int)

        final_df = final_df.sort_values(
            by=[("", "Sequencing_ID"), ("", "Barcode")], ascending=True
        )

        return final_df

    def save_results_to_file(self, final_df):
        """Save results to text and Excel files with proper handling of headers and formatting."""
        try:
            final_df.to_csv("./ABO_result.txt", sep="\t", index=False)
            print("Results saved successfully to text file.")
        except Exception as txt_err:
            print(f"Error saving to text file: {txt_err}")
            return

        try:
            read_count_cols = []
            if isinstance(final_df.columns, pd.MultiIndex):
                for i, col in enumerate(final_df.columns):
                    if col[1] == "#Reads":
                        read_count_cols.append(i)
            else:
                for i, col in enumerate(final_df.columns):
                    if "#Reads" in str(col):
                        read_count_cols.append(i)

            writer = pd.ExcelWriter("./ABO_result.xlsx", engine="xlsxwriter")

            if isinstance(final_df.columns, pd.MultiIndex):
                final_df.columns = final_df.columns.droplevel()

            # Write the DataFrame to Excel
            final_df.to_excel(
                writer, sheet_name="ABO_Result", header=True, index=False, startrow=1
            )

            workbook = writer.book
            worksheet = writer.sheets["ABO_Result"]

            # Define Excel formats
            data_format = workbook.add_format(
                {"bg_color": "white", "font_color": "black", "border": 1}
            )
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "fg_color": "#007399",
                    "border": 1,
                    "font_color": "white",
                }
            )
            red_bg_format = workbook.add_format(
                {"bg_color": "#e2725b", "font_color": "black"}
            )
            orange_bg_format = workbook.add_format(
                {"bg_color": "#ff9a00", "font_color": "black"}
            )

            # Set header alignment
            header_format.set_align("center")
            header_format.set_align("vcenter")

            # Get dimensions
            num_rows, num_cols = final_df.shape

            # Find the Reliability column index (it's the last column)
            reliability_col = xl_col_to_name(num_cols - 1)

            # Sanity check in dev
            print(f"Data has {num_rows} rows, starting at row 3 with two header rows")

            # Apply conditional formatting to each read count column
            for col_idx in read_count_cols:
                col_letter = xl_col_to_name(col_idx)

                # Very low reads (≤20) - red background
                worksheet.conditional_format(
                    f"{col_letter}3:{col_letter}{num_rows + 2}",
                    {
                        "type": "cell",
                        "criteria": "<=",
                        "value": 20,
                        "format": red_bg_format,
                    },
                )

                # Low reads (21-49) - orange background
                worksheet.conditional_format(
                    f"{col_letter}3:{col_letter}{num_rows + 2}",
                    {
                        "type": "cell",
                        "criteria": "between",
                        "minimum": 21,
                        "maximum": 40,
                        "format": orange_bg_format,
                    },
                )

            # Print which columns are being formatted
            print(
                f"Applying read count conditional formatting to columns: {[xl_col_to_name(i) for i in read_count_cols]}"
            )

            try:
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",
                    {
                        "type": "formula",
                        "criteria": f'=${reliability_col}3="Very Low(\u226420 reads)"',
                        "format": red_bg_format,
                    },
                )
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",  # Changed to start at row 3
                    {
                        "type": "formula",
                        "criteria": f'=${reliability_col}3="Low (\u226440 reads)"',  # Changed to reference row 3
                        "format": orange_bg_format,
                    },
                )
            except Exception as format_err:
                print(
                    f"Warning: Could not apply row-level conditional formatting: {format_err}"
                )
            # Write data
            for row in range(num_rows):
                for col in range(num_cols):
                    cell_value = final_df.iat[row, col]
                    if not pd.isna(cell_value):
                        worksheet.write(row + 2, col, cell_value, data_format)

            header_columns = [
                "Exon6_pos22",
                "Exon6_pos27",
                "Exon6_pos29",
                "Exon6_pos58",
                "Exon7_pos422",
                "Exon7_pos428",
                "Exon7_pos429",
                "Exon7_pos431",
                "Exon7_pos93",
                "Exon7_pos165",
                "Exon7_pos272",
                "Exon7_pos307",
                "Exon7_pos371",
                "Exon7_pos446",
                "Exon7_pos680",
                "Exon7_pos687",
            ]

            column_start = 2
            merge_ranges = []

            for header in header_columns:
                start_col = column_start
                end_col = start_col + 9  # Each main header spans 10 columns

                start_letter = xl_col_to_name(start_col)
                end_letter = xl_col_to_name(end_col)

                merge_ranges.append((f"{start_letter}1:{end_letter}1", header))
                column_start = end_col + 1

            result_start = xl_col_to_name(column_start)
            result_end = xl_col_to_name(column_start + 3)

            # Merge header ranges
            worksheet.merge_range("A1:B1", "Sample", header_format)
            worksheet.merge_range(
                f"{result_start}1:{result_end}1", "Result", header_format
            )

            for merge_range in merge_ranges:
                worksheet.merge_range(merge_range[0], merge_range[1], header_format)

            for col in range(num_cols):
                cell_value = final_df.columns[col]
                if not pd.isna(cell_value):
                    worksheet.write(1, col, cell_value, header_format)

            writer.close()
            print("Results saved successfully to Excel file.")
        except Exception as excel_err:
            print(f"Error saving to Excel file: {excel_err}")
            import traceback

            traceback.print_exc()

        self.df_for_lis_soft = pd.DataFrame()
        self.df_for_lis_soft["Sample ID"] = final_df["Sequencing_ID"]
        self.df_for_lis_soft["Shipment Date"] = ""

        if (
            "Genotype" in final_df.columns
            and not final_df["Genotype"].isnull().all()
            and not (final_df["Genotype"] == "Unknown").all()
        ):
            valid_genotype_mask = (final_df["Genotype"] != "Unknown") & final_df[
                "Genotype"
            ].notnull()
            self.df_for_lis_soft.loc[valid_genotype_mask, "ABO Geno Type1"] = (
                final_df.loc[valid_genotype_mask, "Genotype"].str[0]
            )
            self.df_for_lis_soft.loc[valid_genotype_mask, "ABO Geno Type2"] = (
                final_df.loc[valid_genotype_mask, "Genotype"].str[1]
            )
        else:
            self.df_for_lis_soft["ABO Geno Type1"] = ""
            self.df_for_lis_soft["ABO Geno Type2"] = ""

        self.df_for_lis_soft["ABO Pheno Type"] = final_df["Phenotype"]
        self.df_for_lis_soft["RH"] = ""
        self.df_for_lis_soft["Blood Type"] = final_df["Phenotype"]
        self.df_for_lis_soft["ABORH Comments"] = ""

        if isinstance(final_df.columns, pd.MultiIndex):
            reads_df = final_df.loc[:, (slice(None), "#Reads")]
            self.df_for_lis_soft["#Reads"] = reads_df.mean(axis=1)
        else:
            reads_columns = [col for col in final_df.columns if "#Reads" in str(col)]
            if reads_columns:
                self.df_for_lis_soft["#Reads"] = final_df[reads_columns].mean(axis=1)
            else:
                self.df_for_lis_soft["#Reads"] = 0

        self.df_for_lis_soft.drop_duplicates(inplace=True)
        self.df_for_lis_soft.to_csv("./final_export.csv", index=False, encoding="utf-8")
        print(
            f"LIS export file created successfully with {len(self.df_for_lis_soft)} samples"
        )

    def run(self):
        """
        Run the ABOReportParser.
        This method processes files, merges dataframes, and saves results.
        """
        self.process_files()
        final_df = self.merge_dataframes()
        print("\n\nFinal Results:")
        print("-" * 336)
        print(final_df.to_string(index=False))
        print("-" * 336)
        self.save_results_to_file(final_df)

        if self.failed_samples:
            print("\n\nFailed Samples Summary:")
            print("-" * 80)
            for sample in self.failed_samples:
                print(f"Sample: {sample['sample']} - Reason: {sample['reason']}")
            print("-" * 80)
            print(f"Total failed samples: {len(self.failed_samples)}")
        else:
            print("\nAll samples processed successfully.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("\nUsage: python ABOReportParser.py <input_directory>\n")
        sys.exit(1)
    input_directory = sys.argv[1]
    parser = ABOReportParser(input_directory)
    parser.run()
    print("All done!\n")
