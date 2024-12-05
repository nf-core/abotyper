#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import sys
import pandas as pd


__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2024, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = ["Fredrick Mobegi", "Benedict Matern", "Mathijs Groeneweg"]
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

    def initialize_columns(self):
        """
        Define the column headers.
        """
        exon6 = ["Exon6_pos22"] * 10
        exon7_422 = ["Exon7_pos422"] * 10
        exon7_428 = ["Exon7_pos428"] * 10
        exon7_429 = ["Exon7_pos429"] * 10
        exon7_431 = ["Exon7_pos431"] * 10
        max_len = max(
            len(exon6), len(exon7_422), len(exon7_428), len(exon7_429), len(exon7_431)
        )
        exon6 += [""] * (max_len - len(exon6))
        exon7_422 += [""] * (max_len - len(exon7_422))
        exon7_428 += [""] * (max_len - len(exon7_428))
        exon7_429 += [""] * (max_len - len(exon7_429))
        exon7_431 += [""] * (max_len - len(exon7_431))
        header_cols = (
            ["", ""]
            + exon6
            + exon7_422
            + exon7_428
            + exon7_429
            + exon7_431
            + ["", "", "", ""]
        )
        header_rows = (
            ["Barcode", "Sequencing_ID"]
            + ["#Reads", "Mat", "Mis", "Ins", "Del", "A", "G", "C", "T", "Type"] * 5
            + ["Phenotype", "Genotype", "ExtendedGenotype", "Reliability"]
        )
        self.columns = pd.MultiIndex.from_arrays([header_cols, header_rows])
        # print(self.columns)

    def parse_exon7(self, filename):
        """
        Open the file for reading and processing!
        """
        with open(filename, "r", encoding="utf-8") as f:
            lines = f.readlines()
            # Pos 422
            pos1 = int(lines[2].split(":")[1].strip())
            count1 = int(lines[6].split(":")[1].strip())
            mat1, mis1, ins1, del1, a1, g1, c1, t1 = [int(x) for x in lines[8].split()]
            # Pos 428
            pos2 = int(lines[10].split(":")[1].strip())
            count2 = int(lines[14].split(":")[1].strip())
            mat2, mis2, ins2, del2, a2, g2, c2, t2 = [int(x) for x in lines[16].split()]
            # Pos 429
            pos3 = int(lines[18].split(":")[1].strip())
            count3 = int(lines[22].split(":")[1].strip())
            mat3, mis3, ins3, del3, a3, g3, c3, t3 = [int(x) for x in lines[24].split()]
            # Pos 431
            pos4 = int(lines[26].split(":")[1].strip())
            count4 = int(lines[30].split(":")[1].strip())
            mat4, mis4, ins4, del4, a4, g4, c4, t4 = [int(x) for x in lines[32].split()]

        df = pd.DataFrame(
            {
                "Exon": ["7", "7", "7", "7"],
                "Position": [pos1, pos2, pos3, pos4],
                "#Reads": [count1, count2, count3, count4],
                "Mat": [mat1, mat2, mat3, mat4],
                "Mis": [mis1, mis2, mis3, mis4],
                "Ins": [ins1, ins2, ins3, ins4],
                "Del": [del1, del2, del3, del4],
                "A": [a1, a2, a3, a4],
                "G": [g1, g2, g3, g4],
                "C": [c1, c2, c3, c4],
                "T": [t1, t2, t3, t4],
            }
        )

        df["Position"] = df["Position"].apply(lambda x: x)

        def get_type(row):
            """Add the 'Type' column based on values in the [Ins Del A G C T] columns"""
            if row["Position"] == 422:
                if row["A"] >= 80:
                    return "B"
                elif row["C"] >= 80:
                    return "A or O"
                elif abs(row["A"] - row["C"]) <= 20:
                    return "(A or O) and B"
                elif 15 < row["A"] < 80 and 15 < row["C"] < 80:
                    return "(A or O) and B"

            if row["Position"] == 428:
                if row["G"] >= 80:
                    return "O and (A or B)"
                elif row["A"] >= 80:
                    return "O2"
                elif abs(row["G"] - row["A"]) <= 20:
                    return "O2 and (O or A or B)"
                elif 15 < row["G"] < 80 and 15 < row["A"] < 80:
                    return "O2 and (O or A or B)"

            if row["Position"] == 429:
                if row["G"] >= 80:
                    return "A or O"
                elif row["C"] >= 80:
                    return "B"
                elif abs(row["G"] - row["C"]) <= 20:
                    return "(A or O) and B"
                elif 15 < row["G"] < 80 and 20 < row["C"] < 80:
                    return "(A or O) and B"

            if row["Position"] == 431:
                if row["T"] >= 80:
                    return "O and (A or B)"
                elif row["G"] >= 80:
                    return "O3"
                elif abs(row["T"] - row["G"]) <= 20:
                    return "O3 and (O or A or B)"
                elif 20 < row["T"] < 80 and 20 < row["G"] < 80:
                    return "O3 and (O or A or B)"
            return ""

        df["Type"] = df.apply(get_type, axis=1)

        # Reorder the columns
        df = df[
            [
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
        ]
        # print(df)
        return df

    def parse_exon6(self, filename):
        """Parse exon 6 and extract relevant data"""
        with open(filename, "r", encoding="utf-8") as f:
            lines = f.readlines()
            # Extract the required lines from the first section
            pos1 = int(lines[2].split(":")[1].strip())
            count1 = int(lines[6].split(":")[1].strip())
            mat1, mis1, ins1, del1, a1, g1, c1, t1 = [int(x) for x in lines[8].split()]

        # Create a DataFrame from the extracted values
        df = pd.DataFrame(
            {
                "Exon": ["6"],
                "Position": [pos1],
                "#Reads": [count1],
                "Mat": [mat1],
                "Mis": [mis1],
                "Ins": [ins1],
                "Del": [del1],
                "A": [a1],
                "G": [g1],
                "C": [c1],
                "T": [t1],
            }
        )

        # Remove the string "Exon 6 pos" from the Position column
        df["Position"] = df["Position"].apply(lambda x: x)
        max_g = df["G"].max()
        max_del = df["Del"].max()

        if max_g >= 80 and max_g > max_del:
            df["Type"] = "A or B or O"
        elif max_del >= 80 and max_del > max_g:
            df["Type"] = "O1"
        elif abs(max_g + max_del) >= 20:
            df["Type"] = "O1 and (A or B or O)"
        elif 20 < df["Del"] < 80 and 20 < df["G"] < 80:
            df["Type"] = "O1 and (A or B or O)"
        else:
            df["Type"] = ""
        # Reorder the columns
        df = df[
            [
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
        ]
        return df

    def add_phenotype_genotype(self, df):
        """Find the column numbers with "Type" in the column name"""
        type_cols = [i for i, col in enumerate(df.columns) if "Type" in col]
        # Assign the column numbers to Type1, Type2, and Type3
        Type1 = df.columns[type_cols[0]]
        Type2 = df.columns[type_cols[1]]
        Type3 = df.columns[type_cols[2]]
        Type4 = df.columns[type_cols[3]]
        Type5 = df.columns[type_cols[4]]
        # Initialize the Reliability column
        df["Reliability"] = ""

        for i in range(len(df)):
            if (
                (df.loc[i, Type1] == "A or B").all()
                and (df.loc[i, Type2] == "A or O").all()
                and (df.loc[i, Type3] == "A or O").all()
                and (df.loc[i, Type4] == "A or O").all()
                and (df.loc[i, Type5] == "A or O")
            ).any():
                df.at[i, "Phenotype"] = "A"
                df.at[i, "Genotype"] = "AA"
                df.at[i, "ExtendedGenotype"] = "AA"
            else:
                df.at[i, "Phenotype"] = ""
                df.at[i, "Genotype"] = ""
                df.at[i, "ExtendedGenotype"] = ""
        return df

    def assign_phenotype_genotype(self, df):
        """Assign the phenotype and genotype information"""
        type_exon6 = df[("Exon6_pos22", "Type")]
        type_exon7_422 = df[("Exon7_pos422", "Type")]
        type_exon7_428 = df[("Exon7_pos428", "Type")]
        type_exon7_429 = df[("Exon7_pos429", "Type")]
        type_exon7_431 = df[("Exon7_pos431", "Type")]

        nreads6 = df[("Exon6_pos22", "#Reads")]
        nreads_exon7_p422 = df[("Exon7_pos422", "#Reads")]
        nreads_exon7_p428 = df[("Exon7_pos428", "#Reads")]
        nreads_exon7_p429 = df[("Exon7_pos429", "#Reads")]
        nreads_exon7_p431 = df[("Exon7_pos431", "#Reads")]

        ## OA COMBINATIONS ---------------------------------------------------------------------------
        ## combination 1 | AO1 --
        if (
            (type_exon6 == "O1 and (A or B or O)").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "A"
            Genotype = "AO"
            ExtendedGenotype = "AO1"
            # Reliability = 'Enter-manually'

        ## combination 2 | AO2 --
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O2 and (O or A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "A"
            Genotype = "AO"
            ExtendedGenotype = "AO2"
            # Reliability = 'Enter-manually'

        ## combination 3 | AO3 --
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O3 and (O or A or B)").all()
        ):
            Phenotype = "A"
            Genotype = "AO"
            ExtendedGenotype = "AO3"
            # Reliability = 'Enter-manually'

        ## OB COMBINATIONS ---------------------------------------------------------------------------
        ## combination 4 | BO1 --
        elif (
            (type_exon6 == "O1 and (A or B or O)").all()
            and (type_exon7_422 == "(A or O) and B").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "(A or O) and B").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "B"
            Genotype = "BO"
            ExtendedGenotype = "BO1"
            # Reliability = 'Enter-manually'

        ## combination 5 | O2B --
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "(A or O) and B").all()
            and (type_exon7_428 == "O2 and (O or A or B)").all()
            and (type_exon7_429 == "(A or O) and B").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "B"
            Genotype = "BO"
            ExtendedGenotype = "O2B"
            # Reliability = 'Enter-manually'

        ## combination 6 | AO3 --
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "(A or O) and B").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "(A or O) and B").all()
            and (type_exon7_431 == "O3 and (O or A or B)").all()
        ):
            Phenotype = "B"
            Genotype = "BO"
            ExtendedGenotype = "BO3"
            # Reliability = 'Enter-manually'

        ## OO COMBINATIONS  ---------------------------------------------------------------------------
        ## combination 7 | O1O2 --
        elif (
            (type_exon6 == "O1 and (A or B or O)").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O2 and (O or A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "O"
            Genotype = "OO"
            ExtendedGenotype = "O1O2"
            # Reliability = 'Enter-manually'

        ## combination 8 | O1O3 --
        elif (
            (type_exon6 == "O1 and (A or B or O)").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O3 and (O or A or B)").all()
        ):
            Phenotype = "O"
            Genotype = "OO"
            ExtendedGenotype = "O1O3"
            # Reliability = 'Enter-manually'

        ## combination 9 | O2O3 --
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O2 and (O or A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O3 and (O or A or B)").all()
        ):
            Phenotype = "O"
            Genotype = "OO"
            ExtendedGenotype = "O2O3"
            # Reliability = 'Enter-manually'

        ## combination 10 | O1O1 --
        elif (
            (type_exon6 == "O1").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "O"
            Genotype = "OO"
            ExtendedGenotype = "O1O1"
            # Reliability = 'Enter-manually'

        ## combination 11 | O2O2 --
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O2").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "O"
            Genotype = "OO"
            ExtendedGenotype = "O2O2"
            # Reliability = 'Enter-manually'
            #
        ## combination 12 | O3O3 --
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O3").all()
        ):
            Phenotype = "O"
            Genotype = "OO"
            ExtendedGenotype = "O3O3"
            # Reliability = 'Enter-manually'

        ## combination 13 | AA ---------------------------------------------------------------------------
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "A or O").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "A or O").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "A"
            Genotype = "AA"
            ExtendedGenotype = "AA"
            # Reliability = 'Enter-manually'

        ## combination 14 | BB ---------------------------------------------------------------------------
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "B").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "B").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "B"
            Genotype = "BB"
            ExtendedGenotype = "BB"
            # Reliability = 'Enter-manually'

        ## combination 15 | AB ---------------------------------------------------------------------------
        elif (
            (type_exon6 == "A or B or O").all()
            and (type_exon7_422 == "(A or O) and B").all()
            and (type_exon7_428 == "O and (A or B)").all()
            and (type_exon7_429 == "(A or O) and B").all()
            and (type_exon7_431 == "O and (A or B)").all()
        ):
            Phenotype = "AB"
            Genotype = "AB"
            ExtendedGenotype = "AB"
            # Reliability = 'Enter-manually'

        ## UNKNOWN None of the above --------------------------------------------------------------------
        else:
            Phenotype = "Unknown"
            Genotype = "Unknown"
            ExtendedGenotype = "Unknown"
            # Reliability = 'Enter-manually'

        if (
            (nreads6 <= 20).all()
            and (nreads_exon7_p422 <= 20).all()
            and (nreads_exon7_p428 <= 20).all()
            and (nreads_exon7_p429 <= 20).all()
            and (nreads_exon7_p431 <= 20).all()
        ):
            Reliability = "Very Low(\u226420 reads)"

        elif (
            ((nreads6 > 20) & (nreads6 < 50)).all()
            and ((nreads_exon7_p422 > 20) & (nreads_exon7_p422 < 50)).all()
            and ((nreads_exon7_p428 > 20) & (nreads_exon7_p428 < 50)).all()
            and ((nreads_exon7_p429 > 20) & (nreads_exon7_p429 < 50)).all()
            and ((nreads_exon7_p431 > 20) & (nreads_exon7_p431 < 50)).all()
        ):
            Reliability = "Low (\u226450 reads)"

        elif (
            (nreads6 >= 500).all()
            and (nreads_exon7_p422 >= 500).all()
            and (nreads_exon7_p428 >= 500).all()
            and (nreads_exon7_p429 >= 500).all()
            and (nreads_exon7_p431 >= 500).all()
        ):
            Reliability = "Robust(\u2265500 reads)"

        else:
            Reliability = ""

        df[("", "Phenotype")] = Phenotype
        df[("", "Genotype")] = Genotype
        df[("", "ExtendedGenotype")] = ExtendedGenotype
        df[("", "Reliability")] = Reliability
        return df

    def process_file(self, filename):
        sample_name, barcode = filename.split("_")
        exon6_dir = os.path.join(self.input_dir, filename, "exon6")
        exon7_dir = os.path.join(self.input_dir, filename, "exon7")

        # Check if both exon6 and exon7 directories exist
        if not (os.path.exists(exon6_dir) and os.path.exists(exon7_dir)):
            print(f"Skipping file {filename}. Missing exon6 or exon7 directory.")
            return

        exon6_phenotypes = os.path.join(exon6_dir, "ABOPhenotype.txt")
        exon7_phenotypes = os.path.join(exon7_dir, "ABOPhenotype.txt")
        exon6_polymorphisms = os.path.join(exon6_dir, "ABOReadPolymorphisms.txt")
        exon7_polymorphisms = os.path.join(exon7_dir, "ABOReadPolymorphisms.txt")

        # Check if both exon6 and exon7 ABOPhenotype.txt files exist
        if not (os.path.exists(exon6_phenotypes) and os.path.exists(exon7_phenotypes)):
            print(
                f"Skipping file {filename}. Missing ABOPhenotype.txt file for exon6 or exon7."
            )
            return
        if (
            os.path.getsize(exon6_polymorphisms) == 0
            or os.path.getsize(exon7_polymorphisms) == 0
        ):
            print(
                f"Skipping file {filename}. Expected polymorphisms not present in alignment file."
            )
            return

        # Define sample_df at the beginning of the method
        sample_df = pd.DataFrame({"Barcode": [barcode], "Sequencing_ID": [sample_name]})

        try:
            df_exon6 = self.parse_exon6(exon6_phenotypes)
        except Exception as e:
            print(f"Error processing exon6 for file {filename}: {str(e)}")
            df_exon6 = None

        try:
            df_exon7 = self.parse_exon7(exon7_phenotypes)
        except Exception as e:
            print(f"Error processing exon7 for file {filename}: {str(e)}")
            df_exon7 = None
        # print(df_exon7)

        if df_exon6 is not None and df_exon7 is not None:
            df_exon7_pos422 = df_exon7.iloc[[0]].reset_index(drop=True)
            df_exon7_pos428 = df_exon7.iloc[[1]].reset_index(drop=True)
            df_exon7_pos429 = df_exon7.iloc[[2]].reset_index(drop=True)
            df_exon7_pos431 = df_exon7.iloc[[3]].reset_index(drop=True)
            merged_df = pd.concat(
                [
                    sample_df,
                    df_exon6,
                    df_exon7_pos422,
                    df_exon7_pos428,
                    df_exon7_pos429,
                    df_exon7_pos431,
                ],
                axis=1,
                join="inner",
            ).drop(["Exon", "Position"], axis=1)
            merged_df["Barcode"] = merged_df["Barcode"].str.replace(
                "barcode", "", case=False
            )
            merged_df["Barcode"] = pd.to_numeric(merged_df["Barcode"], errors="coerce")
            merged_df = self.add_phenotype_genotype(merged_df)
            merged_df.columns = self.columns
            merged_df = self.assign_phenotype_genotype(merged_df)
            self.results.append(merged_df)

    # Module updated to capture the full spectrum of filenames possible from ONT setups
    def process_files(self):
        for filename in os.listdir(self.input_dir):
            if os.path.isdir(os.path.join(self.input_dir, filename)):
                try:
                    # Updated regex pattern
                    # pattern = r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[0-9]+-[0-9]+)?_barcode\d+$"
                    pattern = (
                        r"^(IMM|INGS|NGS|[A-Z0-9]+)(-[0-9]+)?(-[A-Z0-9]+)?_barcode\d+$"
                    )
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
        final_df = final_df.sort_values(by=[("", "Sequencing_ID")])
        return final_df

    def save_results_to_file(self, final_df):
        # Write to text file
        final_df.to_csv("./ABO_result.txt", sep="\t", index=False)
        writer = pd.ExcelWriter("./ABO_result.xlsx", engine="xlsxwriter")
        final_df.columns = final_df.columns.droplevel()
        final_df.to_excel(
            writer, sheet_name="ABO_Result", header=True, index=False, startrow=1
        )

        workbook = writer.book
        worksheet = writer.sheets["ABO_Result"]

        data_format = workbook.add_format(
            {"bg_color": "white", "font_color": "black", "border": 1}
        )
        header_format = workbook.add_format(
            {"bold": True, "fg_color": "#007399", "border": 1, "font_color": "white"}
        )
        red_bg_format = workbook.add_format(
            {"bg_color": "#e2725b", "font_color": "black"}
        )
        orange_bg_format = workbook.add_format(
            {"bg_color": "#ff9a00", "font_color": "black"}
        )
        # Remove the assignment statement for green_bg_format
        workbook.add_format({"bg_color": "#9caf88", "font_color": "black"})

        header_format.set_align("center")
        header_format.set_align("vcenter")

        num_rows, num_cols = final_df.shape
        f'A1:{chr(ord("A") + num_cols - 1)}{num_rows}'

        try:
            # Very low number of reads < 20
            worksheet.conditional_format(
                "A1:BD5000",
                {
                    "type": "formula",
                    "criteria": '=$BD1="Very Low(\u226420 reads)"',
                    "format": red_bg_format,
                },
            )
            ## Low >20 but < 50
            worksheet.conditional_format(
                "A1:BD5000",
                {
                    "type": "formula",
                    "criteria": '=$BD1="Low (\u226450 reads)"',
                    "format": orange_bg_format,
                },
            )
        except Exception as e:
            print(
                f"\nAn error occurred while applying conditional formatting: {str(e)}"
            )

        for row in range(2, num_rows + 2):  # Add 2 to account for the header row
            for col in range(num_cols):
                cell_value = final_df.iat[row - 2, col]
                if not pd.isna(cell_value):
                    worksheet.write(row, col, cell_value, data_format)

        header_columns = [
            "Exon6_pos22",
            "Exon7_pos422",
            "Exon7_pos428",
            "Exon7_pos429",
            "Exon7_pos431",
        ]
        merge_ranges = [
            ("C1:L1", header_columns[0]),
            ("M1:V1", header_columns[1]),
            ("W1:AF1", header_columns[2]),
            ("AG1:AP1", header_columns[3]),
            ("AQ1:AZ1", header_columns[4]),
        ]
        worksheet.merge_range("A1:B1", "Sample", header_format)
        worksheet.merge_range("BA1:BD1", "Result", header_format)

        for merge_range in merge_ranges:
            worksheet.merge_range(merge_range[0], merge_range[1], header_format)

        for col in range(num_cols):
            cell_value = final_df.columns[col]
            if not pd.isna(cell_value):
                worksheet.write(1, col, cell_value, header_format)

        writer.close()

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
        else:
            # Single-level index case
            reads_df = final_df.filter(like="#Reads")

        # Compute the average of the '#Reads' columns
        self.df_for_lis_soft["#Reads"] = reads_df.mean(axis=1)

        self.df_for_lis_soft.drop_duplicates(inplace=True)
        self.df_for_lis_soft.to_csv("./final_export.csv", index=False, encoding="utf-8")

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
