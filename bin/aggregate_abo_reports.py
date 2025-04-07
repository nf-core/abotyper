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
    "Claude 3.7 Sonnet Thinking (rewrite to add A1/A2/A3 subtypes)",
]
__license__ = "GPL"
__version__ = "0.2.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Development"

"""
A script to collate all ABO phenotype results from each sample into an 
Excel worksheet and a CSV for export to LIS soft.
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
        # Basic ABO typing positions
        exon6 = ["Exon6_pos22"] * 10  # c.261 position - del is ABO*O

        # A subtype positions in exon 6
        exon6_27 = ["Exon6_pos27"] * 10  # c.266 position - A1/A3 vs A2
        exon6_29 = ["Exon6_pos29"] * 10  # c.268 position - A1/A3 vs A2
        exon6_58 = ["Exon6_pos58"] * 10  # c.297 position - A1/A3 vs A2

        # Primary diagnostic positions in exon 7
        exon7_422 = ["Exon7_pos422"] * 10
        exon7_428 = ["Exon7_pos428"] * 10
        exon7_429 = ["Exon7_pos429"] * 10
        exon7_431 = ["Exon7_pos431"] * 10

        # A subtype positions in exon 7
        exon7_467 = ["Exon7_pos467"] * 10  # A1 vs A2/A3
        exon7_539 = ["Exon7_pos539"] * 10  # A1/A2 vs A3
        exon7_646 = ["Exon7_pos646"] * 10  # A1 vs A2
        exon7_681 = ["Exon7_pos681"] * 10  # A1/A2 vs A3
        exon7_745 = ["Exon7_pos745"] * 10  # A1/A2 vs A3
        exon7_820 = ["Exon7_pos820"] * 10  # A1/A2 vs A3
        exon7_1054 = ["Exon7_pos1054"] * 10  # A1/A3 vs A2
        exon7_1061 = ["Exon7_pos1061"] * 10  # A1 vs A2/A3

        # Ensure all arrays have the same length
        max_len = 10  # All arrays are initialized with 10 elements

        # Construct header columns
        header_cols = (
            ["", ""]  # Barcode and Sequencing_ID
            + exon6
            # Add exon 6 A subtype positions
            + exon6_27
            + exon6_29
            + exon6_58
            # Primary positions
            + exon7_422
            + exon7_428
            + exon7_429
            + exon7_431
            # All A subtype positions
            + exon7_467
            + exon7_539
            + exon7_646
            + exon7_681
            + exon7_745
            + exon7_820
            + exon7_1054
            + exon7_1061
            + ["", "", "", ""]  # Result columns
        )

        # Construct header rows
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
            + column_metrics * 16  # 16 positions total (4 exon6 + 12 exon7)
            + ["Phenotype", "Genotype", "ExtendedGenotype", "Reliability"]
        )

        # Create MultiIndex for DataFrame columns
        self.columns = pd.MultiIndex.from_arrays([header_cols, header_rows])

    def parse_exon7(self, filename):
        """
        Open the file for reading and processing all exon 7 positions!
        """
        try:
            with open(filename, "r", encoding="utf-8") as f:
                lines = f.readlines()

            # Lists to store data for each position
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

            # Process file line by line
            i = 0
            while i < len(lines):
                line = lines[i].strip()

                # Find position markers
                if "Exon 7 position(1-based):" in line:
                    pos_match = re.search(r":\s*(\d+)", line)
                    if pos_match:
                        pos = int(pos_match.group(1))

                        # Look for read count (typically 4 lines after position)
                        for j in range(i, min(i + 10, len(lines))):
                            if "Aligned Read Count:" in lines[j]:
                                count_match = re.search(r":\s*(\d+)", lines[j])
                                if count_match:
                                    count = int(count_match.group(1))

                                    # Look for stats line (typically 2 lines after read count)
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

                                            # Store all values
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

                                            # Skip to next position
                                            break
                i += 1

            # Create DataFrame
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
            # Primary positions + all A subtype positions
            all_positions = [
                422,
                428,
                429,
                431,  # Primary positions
                467,
                539,
                646,
                681,
                745,
                820,
                1054,
                1061,  # A subtypes
            ]

            for pos in all_positions:
                if pos not in df["Position"].values:
                    # Add empty row for missing position
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

            # Sort by position
            df = df.sort_values("Position").reset_index(drop=True)

            # AFTER adding all positions, apply get_type
            df["Type"] = df.apply(
                lambda row: self.get_type(
                    row["Position"], row["A"], row["G"], row["C"], row["T"]
                ),
                axis=1,
            )

            return df

        except Exception as e:
            print(f"Error parsing exon 7 file {filename}: {str(e)}")
            # Return empty DataFrame with all required positions
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

            # Use the same all_positions list as above
            for pos in [
                422,
                428,
                429,
                431,  # Primary positions
                467,
                539,
                646,
                681,
                745,
                820,
                1054,
                1061,  # A subtypes
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

            # Lists to store data for each position
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

            # Process file line by line
            i = 0
            while i < len(lines):
                line = lines[i].strip()

                # Find position markers
                if "Exon 6 position(1-based):" in line:
                    pos_match = re.search(r":\s*(\d+)", line)
                    if pos_match:
                        pos = int(pos_match.group(1))

                        # Look for read count (typically 4 lines after position)
                        for j in range(i, min(i + 10, len(lines))):
                            if "Aligned Read Count:" in lines[j]:
                                count_match = re.search(r":\s*(\d+)", lines[j])
                                if count_match:
                                    count = int(count_match.group(1))

                                    # Look for stats line (typically 2 lines after read count)
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

                                            # Store all values
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

                                            # Skip to next position
                                            break
                i += 1

            # Create DataFrame
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

            # Make sure all needed positions are in the DataFrame
            all_positions = [22, 27, 29, 58]  # All exon 6 positions

            for pos in all_positions:
                if pos not in df["Position"].values:
                    # Add empty row for missing position
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

            # Sort by position
            df = df.sort_values("Position").reset_index(drop=True)

            # Apply type determination for each position
            df["Type"] = df.apply(
                lambda row: self.get_type_exon6(
                    row["Position"], row["A"], row["G"], row["C"], row["T"], row["Del"]
                ),
                axis=1,
            )

            return df

        except Exception as e:
            print(f"Error parsing exon 6: {str(e)}")
            # Return empty DataFrame with all required positions
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

            # Add all required positions
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
        if pos == 22:
            # Original position 22 typing logic
            if g >= 80 and g > dele:
                return "A or B or O"
            elif dele >= 80 and dele > g:
                return "O1"
            elif abs(g + dele) >= 20:
                return "O1 and (A or B or O)"
            elif 20 < dele < 80 and 20 < g < 80:
                return "O1 and (A or B or O)"

        elif pos == 27:  # c.266 - A1/A3 vs A2
            if c >= 80:
                return "A1 or A3"
            elif t >= 80:
                return "A2"
            elif 20 < c < 80 and 20 < t < 80:
                return "A1/A3 and A2"

        elif pos == 29:  # c.268 - A1/A3 vs A2
            if t >= 80:
                return "A1 or A3"
            elif c >= 80:
                return "A2"
            elif 20 < t < 80 and 20 < c < 80:
                return "A1/A3 and A2"

        elif pos == 58:  # c.297 - A1/A3 vs A2
            if a >= 80:
                return "A1 or A3"
            elif g >= 80:
                return "A2"
            elif 20 < a < 80 and 20 < g < 80:
                return "A1/A3 and A2"

        return ""

    def add_phenotype_genotype(self, df):
        """
        Determine and add phenotype and genotype information based on all marker positions.

        Args:
            df: DataFrame containing all the position data

        Returns:
            Updated DataFrame with added phenotype and genotype columns
        """
        try:
            # Extract type information for all positions
            positions = {
                # Primary positions
                "exon6_22": (
                    df.loc[0, ("Exon6_pos22", "Type")]
                    if ("Exon6_pos22", "Type") in df.columns
                    else ""
                ),
                # A subtype positions in exon 6
                "exon6_27": (
                    df.loc[0, ("Exon6_pos27", "Type")]
                    if ("Exon6_pos27", "Type") in df.columns
                    else ""
                ),
                "exon6_29": (
                    df.loc[0, ("Exon6_pos29", "Type")]
                    if ("Exon6_pos29", "Type") in df.columns
                    else ""
                ),
                "exon6_58": (
                    df.loc[0, ("Exon6_pos58", "Type")]
                    if ("Exon6_pos58", "Type") in df.columns
                    else ""
                ),
                # Primary positions in exon 7
                "pos422": (
                    df.loc[0, ("Exon7_pos422", "Type")]
                    if ("Exon7_pos422", "Type") in df.columns
                    else ""
                ),
                "pos428": (
                    df.loc[0, ("Exon7_pos428", "Type")]
                    if ("Exon7_pos428", "Type") in df.columns
                    else ""
                ),
                "pos429": (
                    df.loc[0, ("Exon7_pos429", "Type")]
                    if ("Exon7_pos429", "Type") in df.columns
                    else ""
                ),
                "pos431": (
                    df.loc[0, ("Exon7_pos431", "Type")]
                    if ("Exon7_pos431", "Type") in df.columns
                    else ""
                ),
                # A subtype positions in exon 7
                "pos467": (
                    df.loc[0, ("Exon7_pos467", "Type")]
                    if ("Exon7_pos467", "Type") in df.columns
                    else ""
                ),
                "pos539": (
                    df.loc[0, ("Exon7_pos539", "Type")]
                    if ("Exon7_pos539", "Type") in df.columns
                    else ""
                ),
                "pos646": (
                    df.loc[0, ("Exon7_pos646", "Type")]
                    if ("Exon7_pos646", "Type") in df.columns
                    else ""
                ),
                "pos681": (
                    df.loc[0, ("Exon7_pos681", "Type")]
                    if ("Exon7_pos681", "Type") in df.columns
                    else ""
                ),
                "pos745": (
                    df.loc[0, ("Exon7_pos745", "Type")]
                    if ("Exon7_pos745", "Type") in df.columns
                    else ""
                ),
                "pos820": (
                    df.loc[0, ("Exon7_pos820", "Type")]
                    if ("Exon7_pos820", "Type") in df.columns
                    else ""
                ),
                "pos1054": (
                    df.loc[0, ("Exon7_pos1054", "Type")]
                    if ("Exon7_pos1054", "Type") in df.columns
                    else ""
                ),
                "pos1061": (
                    df.loc[0, ("Exon7_pos1061", "Type")]
                    if ("Exon7_pos1061", "Type") in df.columns
                    else ""
                ),
            }

            # Get read counts for reliability assessment (primary positions only)
            read_counts = {
                "exon6_22": (
                    df.loc[0, ("Exon6_pos22", "#Reads")]
                    if ("Exon6_pos22", "#Reads") in df.columns
                    else 0
                ),
                "pos422": (
                    df.loc[0, ("Exon7_pos422", "#Reads")]
                    if ("Exon7_pos422", "#Reads") in df.columns
                    else 0
                ),
                "pos428": (
                    df.loc[0, ("Exon7_pos428", "#Reads")]
                    if ("Exon7_pos428", "#Reads") in df.columns
                    else 0
                ),
                "pos429": (
                    df.loc[0, ("Exon7_pos429", "#Reads")]
                    if ("Exon7_pos429", "#Reads") in df.columns
                    else 0
                ),
                "pos431": (
                    df.loc[0, ("Exon7_pos431", "#Reads")]
                    if ("Exon7_pos431", "#Reads") in df.columns
                    else 0
                ),
            }

            # Initialize output data
            result = {
                "phenotype": "Unknown",
                "genotype": "Unknown",
                "extended_genotype": "Unknown",
                "reliability": "Unknown",
            }

            # Add more implementation specific to this method in the next step
            # Logic will be implemented to determine blood types

            # Return the updated dataframe with phenotype information
            df.loc[0, ("", "Phenotype")] = result["phenotype"]
            df.loc[0, ("", "Genotype")] = result["genotype"]
            df.loc[0, ("", "ExtendedGenotype")] = result["extended_genotype"]
            df.loc[0, ("", "Reliability")] = result["reliability"]

            return df

        except Exception as e:
            print(f"Error in add_phenotype_genotype: {str(e)}")
            # Set error values
            df.loc[0, ("", "Phenotype")] = "Error"
            df.loc[0, ("", "Genotype")] = "Error"
            df.loc[0, ("", "ExtendedGenotype")] = "Error"
            df.loc[0, ("", "Reliability")] = "Error processing"
            return df

    def get_type(self, pos, a, g, c, t):
        """
        Determine the blood type or subtype for each position based on nucleotide percentages.
        This is used for first-pass analysis regardless of phenotype.

        Args:
            pos: The position number
            a, g, c, t: Nucleotide percentages
        """
        # Primary diagnostic positions
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
            if g >= 80:
                return "O and (A or B)"
            elif a >= 80:
                return "O2"
            elif abs(g - a) <= 20:
                return "O2 and (O or A or B)"
            elif 15 < g < 80 and 15 < a < 80:
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
            elif 20 < t < 80 and 20 < g < 80:
                return "O3 and (O or A or B)"

        # A subtype positions
        elif pos == 467:
            if c >= 80:
                return "A1"
            elif t >= 80:
                return "A2 or A3"
            elif 20 < c < 80 and 20 < t < 80:
                return "A1 and (A2 or A3)"
        elif pos == 539:  # A1/A2 vs A3
            if c >= 80:
                return "A1 or A2"
            elif t >= 80:
                return "A3"
            elif 20 < c < 80 and 20 < t < 80:
                return "(A1 or A2) and A3"
        elif pos == 646:
            if t >= 80:
                return "A1"
            elif a >= 80:
                return "A2"
            elif 20 < t < 80 and 20 < a < 80:
                return "A1 and A2"
        elif pos == 681:
            if g >= 80:
                return "A1 or A2"
            elif a >= 80:
                return "A3"
            elif 20 < g < 80 and 20 < a < 80:
                return "(A1 or A2) and A3"
        elif pos == 745:  # A1/A2 vs A3
            if c >= 80:
                return "A1 or A2"
            elif t >= 80:
                return "A3"
            elif 20 < c < 80 and 20 < t < 80:
                return "(A1 or A2) and A3"
        elif pos == 820:  # A1/A2 vs A3
            if a >= 80:
                return "A1 or A2"
            elif c >= 80:
                return "A3"
            elif 20 < a < 80 and 20 < c < 80:
                return "(A1 or A2) and A3"
        elif pos == 1054:  # A1/A3 vs A2
            if g >= 80:
                return "A1 or A3"
            elif a >= 80:
                return "A2"
            elif 20 < g < 80 and 20 < a < 80:
                return "(A1 or A3) and A2"
        elif pos == 1061:  # A1 vs A2/A3
            if c >= 80:
                return "A1"
            elif t >= 80:
                return "A2 or A3"
            elif 20 < c < 80 and 20 < t < 80:
                return "A1 and (A2 or A3)"

        return ""

    def assign_phenotype_genotype(self, df):
        """Assign the phenotype and genotype information"""
        try:
            type_exon6 = df.at[0, ("Exon6_pos22", "Type")]
            type_exon7_422 = df.at[0, ("Exon7_pos422", "Type")]
            type_exon7_428 = df.at[0, ("Exon7_pos428", "Type")]
            type_exon7_429 = df.at[0, ("Exon7_pos429", "Type")]
            type_exon7_431 = df.at[0, ("Exon7_pos431", "Type")]

            nreads6 = df.at[0, ("Exon6_pos22", "#Reads")]
            nreads_exon7_p422 = df.at[0, ("Exon7_pos422", "#Reads")]
            nreads_exon7_p428 = df.at[0, ("Exon7_pos428", "#Reads")]
            nreads_exon7_p429 = df.at[0, ("Exon7_pos429", "#Reads")]
            nreads_exon7_p431 = df.at[0, ("Exon7_pos431", "#Reads")]

            # Default values
            Phenotype = "Unknown"
            Genotype = "Unknown"
            ExtendedGenotype = "Unknown"

            # PART 1: PRIMARY GENOTYPING LOGIC - Match exact patterns

            ## OA COMBINATIONS ---------------------------------------------------------------------------
            ## combination 1 | AO1 --
            if (
                type_exon6 == "O1 and (A or B or O)"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "A"
                Genotype = "AO"
                ExtendedGenotype = "AO1"

            ## combination 2 | AO2 --
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O2 and (O or A or B)"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "A"
                Genotype = "AO"
                ExtendedGenotype = "AO2"

            ## combination 3 | AO3 --
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O3 and (O or A or B)"
            ):
                Phenotype = "A"
                Genotype = "AO"
                ExtendedGenotype = "AO3"

            ## OB COMBINATIONS ---------------------------------------------------------------------------
            ## combination 4 | BO1 --
            elif (
                type_exon6 == "O1 and (A or B or O)"
                and type_exon7_422 == "B"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "B"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "BO1"

            ## combination 5 | BO2 --
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "B"
                and type_exon7_428 == "O2 and (O or A or B)"
                and type_exon7_429 == "B"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "BO2"

            ## combination 6 | BO3 --
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "B"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "B"
                and type_exon7_431 == "O3 and (O or A or B)"
            ):
                Phenotype = "B"
                Genotype = "BO"
                ExtendedGenotype = "BO3"

            ## OO COMBINATIONS  ---------------------------------------------------------------------------
            ## combination 7 | O1O2 --
            elif (
                type_exon6 == "O1"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O2"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1O2"

            ## combination 8 | O1O3 --
            elif (
                type_exon6 == "O1"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O3"
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1O3"

            ## combination 9 | O2O3 --
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O2"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O3"
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O2O3"

            ## combination 10 | O1O1 --
            elif (
                type_exon6 == "O1"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O1O1"

            ## combination 11 | O2O2 --
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O2"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O2O2"

            ## combination 12 | O3O3 --
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O3"
            ):
                Phenotype = "O"
                Genotype = "OO"
                ExtendedGenotype = "O3O3"

            ## combination 13 | AA ---------------------------------------------------------------------------
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "A or O"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "A or O"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "A"
                Genotype = "AA"
                ExtendedGenotype = "AA"

                # PART 2: A SUBTYPING - ONLY APPLY FOR AA GENOTYPE
                try:
                    # Get all A subtype positions
                    # From exon 6
                    type_exon6_27 = (
                        df.at[0, ("Exon6_pos27", "Type")]
                        if ("Exon6_pos27", "Type") in df.columns
                        else ""
                    )
                    type_exon6_29 = (
                        df.at[0, ("Exon6_pos29", "Type")]
                        if ("Exon6_pos29", "Type") in df.columns
                        else ""
                    )
                    type_exon6_58 = (
                        df.at[0, ("Exon6_pos58", "Type")]
                        if ("Exon6_pos58", "Type") in df.columns
                        else ""
                    )

                    # From exon 7
                    type_exon7_467 = (
                        df.at[0, ("Exon7_pos467", "Type")]
                        if ("Exon7_pos467", "Type") in df.columns
                        else ""
                    )
                    type_exon7_539 = (
                        df.at[0, ("Exon7_pos539", "Type")]
                        if ("Exon7_pos539", "Type") in df.columns
                        else ""
                    )
                    type_exon7_646 = (
                        df.at[0, ("Exon7_pos646", "Type")]
                        if ("Exon7_pos646", "Type") in df.columns
                        else ""
                    )
                    type_exon7_681 = (
                        df.at[0, ("Exon7_pos681", "Type")]
                        if ("Exon7_pos681", "Type") in df.columns
                        else ""
                    )
                    type_exon7_745 = (
                        df.at[0, ("Exon7_pos745", "Type")]
                        if ("Exon7_pos745", "Type") in df.columns
                        else ""
                    )
                    type_exon7_820 = (
                        df.at[0, ("Exon7_pos820", "Type")]
                        if ("Exon7_pos820", "Type") in df.columns
                        else ""
                    )
                    type_exon7_1054 = (
                        df.at[0, ("Exon7_pos1054", "Type")]
                        if ("Exon7_pos1054", "Type") in df.columns
                        else ""
                    )
                    type_exon7_1061 = (
                        df.at[0, ("Exon7_pos1061", "Type")]
                        if ("Exon7_pos1061", "Type") in df.columns
                        else ""
                    )

                    # Determine A subtypes based on markers
                    a1_markers = False
                    a2_markers = False
                    a3_markers = False

                    # Process exon 6 positions
                    if (
                        "A1 or A3" in str(type_exon6_27)
                        or "A1 or A3" in str(type_exon6_29)
                        or "A1 or A3" in str(type_exon6_58)
                    ):
                        if (
                            "A3" in str(type_exon7_539)
                            or "A3" in str(type_exon7_681)
                            or "A3" in str(type_exon7_745)
                            or "A3" in str(type_exon7_820)
                        ):
                            a3_markers = True
                        else:
                            a1_markers = True

                    if (
                        "A2" in str(type_exon6_27)
                        or "A2" in str(type_exon6_29)
                        or "A2" in str(type_exon6_58)
                    ):
                        a2_markers = True

                    # Position 467 rules
                    if type_exon7_467 == "A1":
                        a1_markers = True
                    elif type_exon7_467 == "A2 or A3":
                        if (
                            "A3" in str(type_exon7_681)
                            or "A3" in str(type_exon7_539)
                            or "A3" in str(type_exon7_745)
                            or "A3" in str(type_exon7_820)
                        ):
                            a3_markers = True
                        else:
                            a2_markers = True

                    # Position 539 rules
                    if "A1 or A2" in str(type_exon7_539):
                        if (
                            "A1" in str(type_exon7_467)
                            or "A1" in str(type_exon7_1061)
                            or "A1" in str(type_exon7_646)
                        ):
                            a1_markers = True
                        else:
                            a2_markers = True
                    elif "A3" in str(type_exon7_539):
                        a3_markers = True

                    # Position 646 rules
                    if type_exon7_646 == "A1":
                        a1_markers = True
                    elif type_exon7_646 == "A2":
                        a2_markers = True

                    # Position 681 rules
                    if "A1 or A2" in str(type_exon7_681):
                        if (
                            "A1" in str(type_exon7_467)
                            or "A1" in str(type_exon7_646)
                            or "A1" in str(type_exon7_1061)
                        ):
                            a1_markers = True
                        else:
                            a2_markers = True
                    elif "A3" in str(type_exon7_681):
                        a3_markers = True

                    # Position 745, 820, 1054, 1061 rules
                    if "A1" in str(type_exon7_1061):
                        a1_markers = True
                    elif "A2 or A3" in str(type_exon7_1061):
                        if (
                            "A3" in str(type_exon7_539)
                            or "A3" in str(type_exon7_681)
                            or "A3" in str(type_exon7_745)
                            or "A3" in str(type_exon7_820)
                        ):
                            a3_markers = True
                        else:
                            a2_markers = True

                    # Special mixed rules
                    if (
                        "A1 and" in str(type_exon7_467)
                        or "A1 and" in str(type_exon7_646)
                        or "A1 and" in str(type_exon7_1061)
                    ):
                        if (
                            "A3" in str(type_exon7_681)
                            or "A3" in str(type_exon7_539)
                            or "A3" in str(type_exon7_745)
                            or "A3" in str(type_exon7_820)
                        ):
                            a1_markers = True
                            a3_markers = True
                        else:
                            a1_markers = True
                            a2_markers = True

                    # Apply A extended genotype
                    if a1_markers and a2_markers and a3_markers:
                        # Just use the two clearest signals
                        if (
                            "A1" in str(type_exon7_467)
                            or "A1" in str(type_exon7_646)
                            or "A1" in str(type_exon7_1061)
                        ):
                            if (
                                "A3" in str(type_exon7_681)
                                or "A3" in str(type_exon7_539)
                                or "A3" in str(type_exon7_745)
                                or "A3" in str(type_exon7_820)
                            ):
                                ExtendedGenotype = "A1A3"
                            else:
                                ExtendedGenotype = "A1A2"
                        else:
                            ExtendedGenotype = "A2A3"
                    elif a1_markers and a2_markers:
                        ExtendedGenotype = "A1A2"
                    elif a1_markers and a3_markers:
                        ExtendedGenotype = "A1A3"
                    elif a2_markers and a3_markers:
                        ExtendedGenotype = "A2A3"
                    elif a1_markers:
                        ExtendedGenotype = "A1A1"
                    elif a2_markers:
                        ExtendedGenotype = "A2A2"
                    elif a3_markers:
                        ExtendedGenotype = "A3A3"
                except Exception as e:
                    print(f"Warning: A subtyping error: {str(e)}")
                    # Keep default AA

            ## combination 14 | BB ---------------------------------------------------------------------------
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "B"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "B"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "B"
                Genotype = "BB"
                ExtendedGenotype = "BB"
                # B subtyping removed as requested

            ## combination 15 | AB ---------------------------------------------------------------------------
            elif (
                type_exon6 == "A or B or O"
                and type_exon7_422 == "(A or O) and B"
                and type_exon7_428 == "O and (A or B)"
                and type_exon7_429 == "(A or O) and B"
                and type_exon7_431 == "O and (A or B)"
            ):
                Phenotype = "AB"
                Genotype = "AB"
                ExtendedGenotype = "AB"

            ## UNKNOWN None of the above --------------------------------------------------------------------
            else:
                Phenotype = "Unknown"
                Genotype = "Unknown"
                ExtendedGenotype = "Unknown"

            # Determine reliability based on read counts
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
                elif min_reads < 50:
                    Reliability = "Low (\u226450 reads)"
                elif min_reads >= 500:
                    Reliability = "Robust(\u2265500 reads)"
                else:
                    Reliability = "Normal"
            else:
                Reliability = "Unknown (no read data)"

            # Update DataFrame with results
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
            # Parse sample name and barcode from filename
            if "_" in filename:
                parts = filename.split("_")
                sample_name = "_".join(parts[:-1])
                barcode = parts[-1]
            else:
                sample_name = filename
                barcode = ""

            # Check for exon6 and exon7 directories
            exon6_dir = os.path.join(self.input_dir, filename, "exon6")
            exon7_dir = os.path.join(self.input_dir, filename, "exon7")

            # Check if both exon6 and exon7 directories exist
            if not (os.path.exists(exon6_dir) and os.path.exists(exon7_dir)):
                print(f"Skipping file {filename}. Missing exon6 or exon7 directory.")
                return

            exon6_phenotypes = os.path.join(exon6_dir, "*.ABOPhenotype.txt")
            exon7_phenotypes = os.path.join(exon7_dir, "*.ABOPhenotype.txt")

            # Find phenotype files
            exon6_phenotype_files = glob.glob(exon6_phenotypes)
            exon7_phenotype_files = glob.glob(exon7_phenotypes)

            # In process_file method where you check for missing files:
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

            # Check if files are empty (0 KB)
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

            # Create empty result DataFrame
            result_df = pd.DataFrame(columns=self.columns)
            result_df.loc[0, ("", "Barcode")] = barcode.replace("barcode", "")
            result_df.loc[0, ("", "Sequencing_ID")] = sample_name

            # Process exon6 data for all positions (22, 27, 29, 58)
            exon6_data = self.parse_exon6(exon6_phenotype_files[0])
            if not exon6_data.empty:
                # Process each exon6 position
                for pos in [22, 27, 29, 58]:  # All exon6 positions
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

            # Process exon7 data for all positions
            exon7_data = self.parse_exon7(exon7_phenotype_files[0])
            if not exon7_data.empty:
                # Process all positions - removed B subtype positions (526, 640, 657)
                all_positions = [
                    422,
                    428,
                    429,
                    431,  # Primary positions
                    467,
                    539,
                    646,
                    681,
                    745,
                    820,
                    1054,
                    1061,  # A subtypes
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

            # Determine phenotype and genotype
            result_df = self.assign_phenotype_genotype(result_df)

            # Add the processed result to the results list
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
                    # Updated regex pattern
                    # pattern = r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[0-9]+-[0-9]+)?_barcode\d+$"
                    pattern = r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[A-Z0-9]+)?(-[A-Z0-9]+)?_barcode\d+$"
                    match = re.match(pattern, filename)

                    if match:
                        print("\nProcessing file: " + filename)

                        # Extract barcode and sample_name from the filename
                        parts = filename.split("_")

                        if len(parts) >= 2:
                            sample_name = parts[0]
                            barcode = parts[
                                -1
                            ]  # Assuming barcode is always the last part
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
                    # Handle exceptions
                    print(f"\nError processing file {filename}: {e}")
                finally:
                    # Any cleanup code goes here
                    print(f"Finished processing file: {filename}")

    def merge_dataframes(self):
        final_df = pd.concat(self.results)

        # Force convert Barcode to int
        final_df[("", "Barcode")] = final_df[("", "Barcode")].astype(int)

        # Sort by two columns nested under "Sample"
        final_df = final_df.sort_values(
            by=[("", "Sequencing_ID"), ("", "Barcode")], ascending=True
        )

        return final_df

    def save_results_to_file(self, final_df):
        """Save results to text and Excel files with proper handling of headers and formatting."""
        # Save to text file
        try:
            final_df.to_csv("./ABO_result.txt", sep="\t", index=False)
            print("Results saved successfully to text file.")
        except Exception as txt_err:
            print(f"Error saving to text file: {txt_err}")
            return

        # Save to Excel file
        try:
            # IMPORTANT: Identify read count columns BEFORE dropping the MultiIndex level
            read_count_cols = []
            if isinstance(final_df.columns, pd.MultiIndex):
                # For MultiIndex columns, find each ExonX_posN/#Reads column
                for i, col in enumerate(final_df.columns):
                    if col[1] == "#Reads":  # Second level is "#Reads"
                        read_count_cols.append(i)
            else:
                # For regular columns, just find any with "#Reads"
                for i, col in enumerate(final_df.columns):
                    if "#Reads" in str(col):
                        read_count_cols.append(i)

            # Now create the Excel writer
            writer = pd.ExcelWriter("./ABO_result.xlsx", engine="xlsxwriter")

            # Check if columns are MultiIndex and drop a level if necessary
            if isinstance(final_df.columns, pd.MultiIndex):
                final_df.columns = final_df.columns.droplevel()

            # Write the DataFrame to Excel
            final_df.to_excel(
                writer, sheet_name="ABO_Result", header=True, index=False, startrow=1
            )

            workbook = writer.book
            worksheet = writer.sheets["ABO_Result"]

            # Define formats
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

            # Add this line to debug
            print(f"Data has {num_rows} rows, starting at row 3 with two header rows")

            # Apply conditional formatting to each read count column
            # We're using the indices we saved BEFORE dropping the level
            for col_idx in read_count_cols:
                col_letter = xl_col_to_name(col_idx)

                # Very low reads (≤20) - red background
                worksheet.conditional_format(
                    f"{col_letter}3:{col_letter}{num_rows + 2}",  # Changed to start at row 3
                    {
                        "type": "cell",
                        "criteria": "<=",
                        "value": 20,
                        "format": red_bg_format,
                    },
                )

                # Low reads (21-49) - orange background
                worksheet.conditional_format(
                    f"{col_letter}3:{col_letter}{num_rows + 2}",  # Changed to start at row 3
                    {
                        "type": "cell",
                        "criteria": "between",
                        "minimum": 21,
                        "maximum": 49,
                        "format": orange_bg_format,
                    },
                )

            # Print which columns are being formatted
            print(
                f"Applying read count conditional formatting to columns: {[xl_col_to_name(i) for i in read_count_cols]}"
            )

            # Apply row-level conditional formatting based on reliability
            try:
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",  # Changed to start at row 3
                    {
                        "type": "formula",
                        "criteria": f'=${reliability_col}3="Very Low(\u226420 reads)"',  # Changed to reference row 3
                        "format": red_bg_format,
                    },
                )
                worksheet.conditional_format(
                    f"A3:{xl_col_to_name(num_cols - 1)}{num_rows + 2}",  # Changed to start at row 3
                    {
                        "type": "formula",
                        "criteria": f'=${reliability_col}3="Low (\u226450 reads)"',  # Changed to reference row 3
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

            # Define all column headers and their ranges
            header_columns = [
                "Exon6_pos22",
                "Exon6_pos27",
                "Exon6_pos29",
                "Exon6_pos58",
                "Exon7_pos422",
                "Exon7_pos428",
                "Exon7_pos429",
                "Exon7_pos431",
                "Exon7_pos467",
                "Exon7_pos539",
                "Exon7_pos646",
                "Exon7_pos681",
                "Exon7_pos745",
                "Exon7_pos820",
                "Exon7_pos1054",
                "Exon7_pos1061",
            ]

            # Calculate merge ranges dynamically
            column_start = 2  # C is column 2 (0-indexed)
            merge_ranges = []

            for header in header_columns:
                start_col = column_start
                end_col = start_col + 9  # Each header spans 10 columns

                # Convert to Excel column letters
                start_letter = xl_col_to_name(start_col)
                end_letter = xl_col_to_name(end_col)

                merge_ranges.append((f"{start_letter}1:{end_letter}1", header))
                column_start = end_col + 1

            # Calculate where the Result columns start
            result_start = xl_col_to_name(column_start)
            result_end = xl_col_to_name(column_start + 3)  # 4 result columns

            # Merge header ranges
            worksheet.merge_range("A1:B1", "Sample", header_format)
            worksheet.merge_range(
                f"{result_start}1:{result_end}1", "Result", header_format
            )

            for merge_range in merge_ranges:
                worksheet.merge_range(merge_range[0], merge_range[1], header_format)

            # Write column headers
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

        ## Creating Final Export file
        # Initialize df_for_lis_soft DataFrame
        self.df_for_lis_soft = pd.DataFrame()
        self.df_for_lis_soft["Sample ID"] = final_df["Sequencing_ID"]
        self.df_for_lis_soft["Shipment Date"] = ""

        # Only process if "Genotype" column exists and all values are not 'Unknown' or blank
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

        # Handle columns based on index type
        if isinstance(final_df.columns, pd.MultiIndex):
            # Multi-level index case
            reads_df = final_df.loc[:, (slice(None), "#Reads")]
            # Compute the average of the '#Reads' columns
            self.df_for_lis_soft["#Reads"] = reads_df.mean(axis=1)
        else:
            # Single-level index case
            reads_columns = [col for col in final_df.columns if "#Reads" in str(col)]
            if reads_columns:
                self.df_for_lis_soft["#Reads"] = final_df[reads_columns].mean(axis=1)
            else:
                self.df_for_lis_soft["#Reads"] = 0

        self.df_for_lis_soft.drop_duplicates(inplace=True)
        self.df_for_lis_soft.to_csv("./final_export.csv", index=False, encoding="utf-8")
        print("Final export file created successfully.")

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

        # Print summary of failed samples
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
    # Check command-line arguments
    if len(sys.argv) != 2:
        print("\nUsage: python ABOReportParser.py <input_directory>\n")
        sys.exit(1)
    input_directory = sys.argv[1]
    parser = ABOReportParser(input_directory)
    parser.run()
    print("All done!\n")

sys.exit(0)
