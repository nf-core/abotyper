#!/usr/bin/env python3

import argparse
import datetime
import io
import json
import os
import pprint
import re
import sys
import numpy as np
import pandas as pd


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
A script to extract variants relevant for determining ABO phenotypes from SAMtools pileup results.
"""


print("=" * 80)
print(f"ABO Blood Type Prediction Script - Started at {datetime.datetime.now()}")
print("=" * 80)

print("Initializing diagnostic variants in exon 6 and 7 ...\n")


def format_number(value):
    """
    Format number as integer if it's whole, otherwise as float with 4 decimals.
    """

    if abs(value - round(value)) < 1e-10:
        return f"{int(value)}"
    else:
        return f"{value:.2f}"


def read_nucleotide_frequencies(input_file):
    """
    Read nucleotide frequency data from input file and return positions and max position.
    """

    print(
        f"\n[{datetime.datetime.now()}] Reading nucleotide frequencies from: {input_file}"
    )
    positions = {}
    max_position = 0
    any_data = False

    try:
        print(f"Checking if input file exists and has content...")
        if os.stat(input_file).st_size == 0:
            print(f"Input file {input_file} is empty.")
            return {}, 0, True

        print(f"Opening input file for reading with pandas...")
        df = pd.read_csv(input_file, sep="\t")
        print(f"Data shape: {df.shape}")
        print(f"Columns found: {list(df.columns)}")

        if df.empty or len(df) < 134:  # Exon 6 minimum length is 135 bp
            print(f"Input file {input_file} contains insufficient data.")
            return {}, 0, True

        print(f"Processing {len(df)} data rows...")
        valid_row_count = 0

        for index, row in df.iterrows():
            try:
                pos = int(row["Ref_Position_1based"])
                ref_base = row["Ref_Base"]
                match_percent = float(row["Match_Percent"])
                mismatch_percent = float(row["Mismatch_Percent"])
                insertion_percent = float(row["Insertion_Percent"])
                deletion_percent = float(row["Deletion_Percent"])
                a_percent = float(row["A_Percent"])
                g_percent = float(row["G_Percent"])
                c_percent = float(row["C_Percent"])
                t_percent = float(row["T_Percent"])

                any_data = True

                positions[pos] = {
                    "ref_base": ref_base,
                    "match_percent": match_percent,
                    "mismatch_percent": mismatch_percent,
                    "insertion_percent": insertion_percent,
                    "deletion_percent": deletion_percent,
                    "A_percent": a_percent,
                    "G_percent": g_percent,
                    "C_percent": c_percent,
                    "T_percent": t_percent,
                    "coverage": 0,
                }

                if pos > max_position:
                    max_position = pos

                valid_row_count += 1

            except (ValueError, KeyError) as e:
                print(f"Warning: Could not parse row {index}: {row}. Error: {e}")
                continue

        print(
            f"Finished processing {len(df)} rows, {valid_row_count} valid positions extracted"
        )

        if not any_data:
            print("Alignment file contains no valid data.")
            return {}, 0, True

        print(
            f"Successfully read {len(positions)} positions, max position: {max_position}"
        )

    except Exception as e:
        print(f"Error reading nucleotide frequency file: {str(e)}", file=sys.stderr)
        print(
            f"Exception details: {type(e).__name__} at line {sys.exc_info()[2].tb_lineno}"
        )
        return {}, 0, True
    return (
        positions,
        max_position,
        False,
    )


def read_coverage_file(coverage_file):
    """
    Read coverage information from coverage statistics file.
    Returns a dictionary with numreads and covbases if found, otherwise empty dict.
    """

    print(
        f"\n[{datetime.datetime.now()}] Reading coverage information from: {coverage_file}"
    )

    try:
        if not os.path.exists(coverage_file):
            print(f"Coverage file does not exist: {coverage_file}")
            return {}

        print(f"Reading coverage file with pandas...")
        df = pd.read_csv(coverage_file, sep="\t", comment=None)

        if len(df) < 1:
            print("Coverage file has no data rows")
            return {}

        print(f"Found {len(df)} data rows with columns: {list(df.columns)}")

        numreads_col = None
        covbases_col = None

        for col in df.columns:
            if "numreads" in col.lower():
                numreads_col = col
            elif "covbases" in col.lower():
                covbases_col = col

        if numreads_col and covbases_col:
            numreads = int(df.iloc[0][numreads_col])
            covbases = int(df.iloc[0][covbases_col])

            coverage_info = {"numreads": numreads, "covbases": covbases}

            print(f"Extracted coverage info: numreads={numreads}, covbases={covbases}")
            return coverage_info
        else:
            print(f"Required columns not found. Available columns: {list(df.columns)}")
            return {}

    except Exception as e:
        print(f"Error reading coverage file: {str(e)}")
        return {}


def determine_exon_type(filename):
    """
    Determine if the input file is for exon6 or exon7 based on filename.
    """

    print(
        f"\n[{datetime.datetime.now()}] Determining exon type from filename: {filename}"
    )
    filename = filename.lower()
    if "exon6" in filename:
        print(f"✓ Exon type determined: exon6 (from filename)")
        return "exon6"
    elif "exon7" in filename:
        print(f"✓ Exon type determined: exon7 (from filename)")
        return "exon7"
    else:
        print(
            f"! Cannot determine exon type from filename. Will fall back to position-based method"
        )
        return None


def generate_exon6_report(positions, output_file, numreads=0, covbases=0):
    """
    Generate report specifically for exon 6.
    """

    print(f"\n[{datetime.datetime.now()}] Generating exon 6 report")
    print(f"Output file: {output_file}")
    print(f"Positions to process: {list(positions.keys())}")
    print(f"Coverage information: numreads={numreads}, covbases={covbases}")

    coverage_display = f" (reads={numreads}, cov={covbases})" if numreads > 0 else ""

    try:
        print(f"Opening output file for writing: {output_file}")
        with open(output_file, "w") as f:
            print("Writing header to file...")
            f.write("Exon 6:\n")

            # Primary ABO*O1 marker - position 22 (c.261)
            if 22 in positions:
                print("Position 22 data found in file")
                pos_data = positions[22]
                print(f"Position 22 data: {pos_data}")

                print("Writing position information to file...")
                f.write("\nExon 6 position(1-based): 22\n")

                print("Writing interpretation text...")
                # Simple interpretation text
                f.write("G nucleotide: A or B blood type.\n")
                f.write("Deletion    : O blood type(O1).\n")

                print("Writing raw data to file...")
                f.write(
                    f"(1-based) Position:22, Reference Base={pos_data['ref_base']}\n"
                )
                f.write(
                    f"Aligned Read Count:{pos_data['coverage']}{coverage_display}\n"
                )
                f.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                f.write(
                    f"{format_number(pos_data['match_percent'])}\t"
                    f"{format_number(pos_data['mismatch_percent'])}\t"
                    f"{format_number(pos_data['insertion_percent'])}\t"
                    f"{format_number(pos_data['deletion_percent'])}\t"
                    f"{format_number(pos_data['A_percent'])}\t"
                    f"{format_number(pos_data['G_percent'])}\t"
                    f"{format_number(pos_data['C_percent'])}\t"
                    f"{format_number(pos_data['T_percent'])}\n"
                )
                print("Finished writing position 22 data")
            else:
                print("! Position 22 data not found in positions dictionary")

            # Variants to determine ABO*A2 subtype in exon 6
            # Map CDS positions to exon positions: c.266 → 27, c.268 → 29, c.297 → 58
            a_subtype_positions = [27, 29, 58]
            available_subtype_positions = sorted(
                [p for p in a_subtype_positions if p in positions]
            )

            if available_subtype_positions:
                print("Writing A subtype header...")
                f.write("\n# -------- A subtypes variants in exon 6 --------\n")

                for pos in available_subtype_positions:
                    print(f"Processing A subtype position {pos}")
                    pos_data = positions[pos]

                    f.write(f"\nExon 6 position(1-based): {pos}\n")

                    if pos == 27:  # c.266C>T
                        f.write("C nucleotide: A1 or A3 subtype.\n")
                        f.write("T nucleotide: A2 subtype.\n")
                    elif pos == 29:  # c.268T>C
                        f.write("T nucleotide: A1 or A3 subtype.\n")
                        f.write("C nucleotide: A2 subtype.\n")
                    elif pos == 58:  # c.297A>G
                        f.write("A nucleotide: A1 or A3 subtype.\n")
                        f.write("G nucleotide: A2 subtype.\n")

                    f.write(
                        f"(1-based) Position:{pos}, Reference Base={pos_data['ref_base']}\n"
                    )
                    f.write(
                        f"Aligned Read Count:{pos_data['coverage']}{coverage_display}\n"
                    )
                    f.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                    f.write(
                        f"{format_number(pos_data['match_percent'])}\t"
                        f"{format_number(pos_data['mismatch_percent'])}\t"
                        f"{format_number(pos_data['insertion_percent'])}\t"
                        f"{format_number(pos_data['deletion_percent'])}\t"
                        f"{format_number(pos_data['A_percent'])}\t"
                        f"{format_number(pos_data['G_percent'])}\t"
                        f"{format_number(pos_data['C_percent'])}\t"
                        f"{format_number(pos_data['T_percent'])}\n"
                    )

        print(f"✓ Exon 6 report successfully written to {output_file}")

    except Exception as e:
        print(f"Error generating Exon 6 report: {str(e)}", file=sys.stderr)
        print(
            f"Exception details: {type(e).__name__} at line {sys.exc_info()[2].tb_lineno}"
        )
        print(f"Creating empty file due to error")
        open(output_file, "w").close()
        print(f"Empty file created due to error: {output_file}")


def generate_exon7_report(positions, output_file, numreads=0, covbases=0):
    """
    Generate report specifically for exon 7.
    """

    print(f"\n[{datetime.datetime.now()}] Generating exon 7 report")
    print(f"Output file: {output_file}")
    print(f"Positions to process: {list(positions.keys())}")
    print(f"Coverage information: numreads={numreads}, covbases={covbases}")

    coverage_display = f" (reads={numreads}, cov={covbases})" if numreads > 0 else ""
    print(f"Coverage display string: '{coverage_display}'")

    try:
        print(f"Opening output file for writing: {output_file}")
        with open(output_file, "w") as f:
            print("Writing header to file...")
            f.write("Exon 7:\n")

            primary_positions = [422, 428, 429, 431]
            print(f"Processing primary diagnostic positions: {primary_positions}")

            primary_exists = any(pos in positions for pos in primary_positions)
            if primary_exists:
                print("Found primary diagnostic positions in data")
            else:
                print("No primary diagnostic positions found in exon 7 data")

            for pos in primary_positions:
                if pos in positions:
                    print(f"Processing primary position {pos}")
                    pos_data = positions[pos]
                    print(f"Position {pos} data: {pos_data}")

                    f.write(f"\nExon 7 position(1-based): {pos}\n")

                    print(f"Writing interpretation text for position {pos}...")
                    if pos == 422:
                        f.write("A nucleotide: B blood type.\n")
                        f.write("C nucleotide: A or O blood type.\n")
                    elif pos == 428:
                        f.write("A nucleotide: O blood type (O2).\n")
                        f.write("G nucleotide: A or B or O blood type.\n")
                    elif pos == 429:
                        f.write("G nucleotide: A or O blood type.\n")
                        f.write("C nucleotide: B blood type.\n")
                    elif pos == 431:
                        f.write("G nucleotide: O blood type (O3).\n")
                        f.write("A nucleotide: O blood type (O4).\n")
                        f.write("T nucleotide: A or B or O blood type.\n")

                    print(f"Writing raw data for position {pos}...")
                    f.write(
                        f"(1-based) Position:{pos}, Reference Base={pos_data['ref_base']}\n"
                    )
                    f.write(
                        f"Aligned Read Count:{pos_data['coverage']}{coverage_display}\n"
                    )
                    f.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                    f.write(
                        f"{format_number(pos_data['match_percent'])}\t"
                        f"{format_number(pos_data['mismatch_percent'])}\t"
                        f"{format_number(pos_data['insertion_percent'])}\t"
                        f"{format_number(pos_data['deletion_percent'])}\t"
                        f"{format_number(pos_data['A_percent'])}\t"
                        f"{format_number(pos_data['G_percent'])}\t"
                        f"{format_number(pos_data['C_percent'])}\t"
                        f"{format_number(pos_data['T_percent'])}\n"
                    )
                else:
                    print(f"Primary position {pos} not found in data")

            # pos. exon7= exonic(CDS) 422(796),428(802),429(803),431(805),93(467),165(539),687(1061)
            # Additional exonic(CDS) 272(646), 307(681), 371(745), 446(820), 680(1054)
            a_subtype_positions = [93, 165, 272, 307, 371, 446, 680, 687]
            print(f"Processing A subtype positions: {a_subtype_positions}")

            available_subtype_positions = sorted(
                [p for p in a_subtype_positions if p in positions]
            )
            print(f"Available A subtype positions: {available_subtype_positions}")

            if available_subtype_positions:
                print("Writing subtype variant header...")
                f.write("\n# -------- A subtypes variants --------\n")

                for pos in available_subtype_positions:
                    print(f"Processing A subtype position {pos}")
                    pos_data = positions[pos]
                    print(f"Position {pos} data: {pos_data}")

                    f.write(f"\nExon 7 position(1-based): {pos}\n")

                    print(f"Writing interpretation text for position {pos}...")
                    # Add interpretation for A subtype positions only
                    # update for A subtype positions [93, 165, 272, 307, 371, 446, 680, 687]
                    if pos == 93:
                        f.write("C nucleotide: A1 subtype.\n")
                        f.write("T nucleotide: A2 or A3 subtype.\n")
                    elif pos == 165:
                        f.write("G nucleotide: A1 or A2 subtype.\n")
                        f.write("A nucleotide: A3 subtype.\n")
                    elif pos == 272:
                        f.write("T nucleotide: A1 subtype.\n")
                        f.write("A nucleotide: A2 subtype.\n")
                    elif pos == 307:
                        f.write("G nucleotide: A1 or A2 subtype.\n")
                        f.write("A nucleotide: A3 subtype.\n")
                    elif pos == 371:
                        f.write("C nucleotide: A1 or A2 subtype.\n")
                        f.write("T nucleotide: A3 subtype.\n")
                    elif pos == 446:
                        f.write("G nucleotide: A1 or A2 subtype.\n")
                        f.write("A nucleotide: A3 subtype.\n")
                    elif pos == 680:
                        f.write("C nucleotide: A1 or A3 subtype.\n")
                        f.write("T nucleotide: A2 subtype.\n")
                    elif pos == 687:
                        f.write("C nucleotide: A1 subtype.\n")
                        f.write("Deletion: A2 or A3 subtype (weaker expression).\n")

                    print(f"Writing raw data for position {pos}...")
                    f.write(
                        f"(1-based) Position:{pos}, Reference Base={pos_data['ref_base']}\n"
                    )
                    f.write(
                        f"Aligned Read Count:{pos_data['coverage']}{coverage_display}\n"
                    )
                    f.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                    f.write(
                        f"{format_number(pos_data['match_percent'])}\t"
                        f"{format_number(pos_data['mismatch_percent'])}\t"
                        f"{format_number(pos_data['insertion_percent'])}\t"
                        f"{format_number(pos_data['deletion_percent'])}\t"
                        f"{format_number(pos_data['A_percent'])}\t"
                        f"{format_number(pos_data['G_percent'])}\t"
                        f"{format_number(pos_data['C_percent'])}\t"
                        f"{format_number(pos_data['T_percent'])}\n"
                    )
            else:
                print("No A subtype positions found in data")

        print(f"✓ Exon 7 report successfully written to {output_file}")

    except Exception as e:
        print(f"Error generating Exon 7 report: {str(e)}", file=sys.stderr)
        print(
            f"Exception details: {type(e).__name__} at line {sys.exc_info()[2].tb_lineno}"
        )
        print(f"Creating empty file due to error")
        open(output_file, "w").close()
        print(f"Empty file created due to error: {output_file}")


def main():
    """
    Main function to process nucleotide frequencies into ABO phenotypes.
    """

    print(f"\n[{datetime.datetime.now()}] Starting main function")
    print("=" * 60)

    print("Parsing command line arguments...")
    parser = argparse.ArgumentParser(
        description="Extract metrics for ABO typing positions from nucleotide frequency data."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input samtools mpileup nucleotide variants metrics",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output filename for ABO phenotype report"
    )
    parser.add_argument(
        "-c", "--coverage", help="Input samtools coverage statistics file"
    )
    parser.add_argument(
        "-e",
        "--exon",
        choices=["exon6", "exon7"],
        help="Explicitly specify which exon data is being processed",
    )

    args = parser.parse_args()
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Coverage file: {args.coverage if args.coverage else 'Not provided'}")
    print(f"Exon type specified: {args.exon if args.exon else 'Not specified'}")

    try:
        print("\nStep 1: Reading nucleotide frequencies...")
        positions, max_position, empty_file = read_nucleotide_frequencies(args.input)
        print(f"Read {len(positions)} positions, max position: {max_position}")
        print(f"Empty file flag: {empty_file}")

        if empty_file or not positions:
            print("No data available. Creating empty file for tracking failures.")
            # Create a completely empty file
            open(args.output, "w").close()
            print(f"Empty file written to: {args.output}")
            return  # Exit function early

        numreads = 0
        covbases = 0

        if args.coverage:
            print("\nStep 2: Reading coverage information...")
            coverage_info = read_coverage_file(args.coverage)
            print(f"Coverage info contains {len(coverage_info)} entries")

            if coverage_info:
                print("Extracting coverage data")
                numreads = coverage_info.get("numreads", 0)
                covbases = coverage_info.get("covbases", 0)
                print(
                    f"Coverage data extracted: numreads={numreads}, covbases={covbases}"
                )

                print(
                    f"Updating coverage information for {len(positions)} positions..."
                )
                for pos in positions:
                    positions[pos]["coverage"] = numreads
                print("Coverage information updated for all positions")
            else:
                print("No coverage information found in coverage file")
        else:
            print("No coverage file provided, using default values")

        print("\nStep 3: Determining exon type...")
        exon_type = args.exon
        if exon_type:
            print(f"Using explicitly specified exon type: {exon_type}")
        else:
            print(
                "Exon type not explicitly specified, trying to determine from filename..."
            )
            exon_type = determine_exon_type(args.input)
            if not exon_type:
                print("Exon type could not be determined from filename")
                print("Falling back to position-based determination...")
                exon_type = "exon7" if max_position > 135 else "exon6"
                print(
                    f"Exon type determined by position: {exon_type} (max position = {max_position})"
                )

        print(f"\nProcessing {exon_type} data")

        print("\nStep 4: Filtering positions relevant to the exon type...")
        filtered_positions = {}

        if exon_type == "exon6":
            print("Processing as exon6 - position 22 and A subtype positions")
            exon6_positions = [22, 27, 29, 58]
            for pos in exon6_positions:
                if pos in positions:
                    filtered_positions[pos] = positions[pos]
                    print(f"Position {pos} found and kept for processing")
                else:
                    print(f"Position {pos} not found in data")

            print("\nStep 5: Generating exon 6 report...")
            generate_exon6_report(filtered_positions, args.output, numreads, covbases)

        else:  # exon7
            print("Processing as exon7")
            # pos. exon7= exonic(CDS) 422(796),428(802),429(803),431(805),93(467),165(539),687(1061)
            # Additional exonic(CDS) 272(646), 307(681), 371(745), 446(820), 680(1054)

            exon7_positions = [
                422,
                428,
                429,
                431,
                93,
                165,
                687,
                272,
                307,
                371,
                446,
                680,
            ]
            print(f"Relevant exon7 positions: {exon7_positions}")

            for pos in exon7_positions:
                if pos in positions:
                    filtered_positions[pos] = positions[pos]
                    print(f"Position {pos} found and kept for processing")
                else:
                    print(f"Position {pos} not found in data")

            print(f"Filtered to {len(filtered_positions)} relevant positions for exon7")

            print("\nStep 5: Generating exon 7 report...")
            generate_exon7_report(filtered_positions, args.output, numreads, covbases)

        print(f"\n✓ ABO phenotype metrics written to: {args.output}")

    except Exception as e:
        print(f"\n! ERROR: {str(e)}", file=sys.stderr)
        print(
            f"! Exception details: {type(e).__name__} at line {sys.exc_info()[2].tb_lineno}"
        )
        print("Creating empty output file due to error")
        open(args.output, "w").close()
        print(f"Empty file written to: {args.output}")

    print(f"\n[{datetime.datetime.now()}] Main function completed")
    print("=" * 60)


if __name__ == "__main__":
    try:
        print(f"Script execution started at: {datetime.datetime.now()}")
        main()
        print(f"Script execution completed successfully at: {datetime.datetime.now()}")
    except Exception as e:
        print(f"! CRITICAL ERROR: Unhandled exception: {str(e)}")
        print(
            f"! Exception details: {type(e).__name__} at line {sys.exc_info()[2].tb_lineno}"
        )
        sys.exit(1)
