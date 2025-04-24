#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import re
import argparse
import sys
import os

__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2024, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = ["Fredrick Mobegi", "Benedict Matern", "Mathijs Groeneweg"]
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

This class uses SAMtools mpileup file to calculates alignment frequencies for all
nucleotides and indels per reference position.
"""


def create_zero_stats(pos, ref_base):
    """Create a stats dictionary with zero values."""
    return {
        "pos": pos,
        "ref_base": ref_base,
        "match_percent": 0,
        "mismatch_percent": 0,
        "insertion_percent": 0,
        "deletion_percent": 0,
        "A_percent": 0,
        "G_percent": 0,
        "C_percent": 0,
        "T_percent": 0,
        "depth": 0,
    }


def write_stats(output_file, stats):
    """Write statistics to the output file."""
    output_file.write(
        f"{stats['pos']}\t{stats['ref_base']}\t{stats['match_percent']}\t"
        f"{stats['mismatch_percent']}\t{stats['insertion_percent']}\t{stats['deletion_percent']}\t"
        f"{stats['A_percent']}\t{stats['G_percent']}\t{stats['C_percent']}\t{stats['T_percent']}\t"
        f"{stats['depth']}\n"
    )


def parse_mpileup_line(line, exon_info=None, total_rows=0):
    """
    Except for 2 positions (exon6pos22 or c.261, and exon7pos687 or c.1061),
    indels in other ABO associated SVN positions are most likely sequencing errors
    and must be handled correctly to avoid mistyping.
    Parse a single line from a mpileup file and calculate nucleotide statistics.
    """
    try:
        fields = line.strip().split("\t")

        if len(fields) < 6:
            return None

        pos = int(fields[1])
        ref_base = fields[2].upper()
        coverage = int(fields[3])
        read_bases = fields[4]

        if coverage == 0:
            return create_zero_stats(pos, ref_base)

        read_bases = re.sub(r"\^.", "", read_bases)
        read_bases = read_bases.replace("$", "")

        matches = 0
        mismatches = 0
        insertions = 0
        deletions = 0
        base_counts = {"A": 0, "G": 0, "C": 0, "T": 0}

        i = 0
        while i < len(read_bases):
            if read_bases[i] in ".,":
                matches += 1
                base_counts[ref_base] += 1
                i += 1
            elif read_bases[i].upper() in "ACGT":
                mismatches += 1
                base = read_bases[i].upper()
                base_counts[base] += 1
                i += 1
            elif read_bases[i] == "+":
                i += 1
                ins_len_str = ""
                while i < len(read_bases) and read_bases[i].isdigit():
                    ins_len_str += read_bases[i]
                    i += 1
                ins_len = int(ins_len_str) if ins_len_str else 0
                i += ins_len
                insertions += 1
            elif read_bases[i] == "-":
                i += 1
                del_len_str = ""
                while i < len(read_bases) and read_bases[i].isdigit():
                    del_len_str += read_bases[i]
                    i += 1
                del_len = int(del_len_str) if del_len_str else 0
                i += del_len
                deletions += 1
            elif read_bases[i] == "*":
                deletions += 1
                i += 1
            else:
                i += 1

        total_events = matches + mismatches + insertions + deletions
        total_nucleotides = matches + mismatches  # Sum of just ATGC counts

        if total_events == 0:
            return create_zero_stats(pos, ref_base)

        # Determine if we should include insertions and deletions in the calculation
        include_indels = True
        is_exon6 = exon_info == "exon6"

        if is_exon6 and coverage < 200 and pos != 22:
            include_indels = False
        elif total_rows > 140 and coverage < 200:
            # Only include indels at key diagnostic positions
            if pos not in [431, 687]:
                include_indels = False

        # When we're ignoring indels, we need to ensure A+T+G+C = 100%
        if include_indels:
            # Using all events as denominator
            denominator = total_events

            # Calculate percentages for bases and indels
            base_percentages = {
                base: int((count / denominator) * 100)
                for base, count in base_counts.items()
            }
            insertion_percent = int((insertions / denominator) * 100)
            deletion_percent = int((deletions / denominator) * 100)
        else:
            # When ignoring indels, use only ATGC counts as denominator
            denominator = total_nucleotides

            if denominator == 0:
                return create_zero_stats(pos, ref_base)

            # Calculate percentages for bases only - should sum to 100%
            base_percentages = {
                base: int((count / denominator) * 100)
                for base, count in base_counts.items()
            }

            # Force sum of ATGC to be 100% by adjusting the reference base
            # This handles any rounding issues
            atgc_sum = sum(base_percentages.values())
            if atgc_sum != 100 and atgc_sum > 0:
                # Adjust the reference base percentage to make sum exactly 100%
                diff = 100 - atgc_sum
                base_percentages[ref_base] += diff

            # Set indel percentages to zero when ignoring them
            insertion_percent = 0
            deletion_percent = 0

        match_percent = base_percentages[ref_base]
        mismatch_percent = sum(
            base_percentages[base] for base in "ACGT" if base != ref_base
        )

        return {
            "pos": pos,
            "ref_base": ref_base,
            "match_percent": match_percent,
            "mismatch_percent": mismatch_percent,
            "insertion_percent": insertion_percent,
            "deletion_percent": deletion_percent,
            "A_percent": base_percentages["A"],
            "G_percent": base_percentages["G"],
            "C_percent": base_percentages["C"],
            "T_percent": base_percentages["T"],
            "depth": coverage,
        }
    except Exception as e:
        print(f"Error processing line: {line.strip()}", file=sys.stderr)
        print(f"Exception: {str(e)}", file=sys.stderr)
        return None


def process_mpileup_file(input_file, output_file, summary_file=None):
    """Process an mpileup file and output nucleotide statistics."""
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist", file=sys.stderr)
        return False

    # Check if output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError as e:
            print(f"Error creating output directory: {str(e)}", file=sys.stderr)
            return False

    # Process the file
    try:
        # First pass: determine reference length and total rows
        unique_positions = set()
        total_rows = 0

        if input_file.endswith(".gz"):
            with gzip.open(input_file, "rt") as f:
                for line in f:
                    total_rows += 1
                    fields = line.strip().split("\t")
                    if len(fields) >= 2:
                        try:
                            unique_positions.add(int(fields[1]))
                        except ValueError:
                            print(
                                f"Warning: Invalid position value in line: {line.strip()}",
                                file=sys.stderr,
                            )
        else:
            with open(input_file, "r") as f:
                for line in f:
                    total_rows += 1
                    fields = line.strip().split("\t")
                    if len(fields) >= 2:
                        try:
                            unique_positions.add(int(fields[1]))
                        except ValueError:
                            print(
                                f"Warning: Invalid position value in line: {line.strip()}",
                                file=sys.stderr,
                            )

        # Determine exon type based on reference length
        ref_length = len(unique_positions)
        exon_info = None
        if 130 <= ref_length <= 140:
            exon_info = "exon6"
        elif 800 <= ref_length <= 830:
            exon_info = "exon7"

        # Second pass: process the file with exon info
        with open(output_file, "w") as out:
            out.write(
                "Ref_Position_1based\tRef_Base\tMatch_Percent\tMismatch_Percent\tInsertion_Percent\tDeletion_Percent\tA_Percent\tG_Percent\tC_Percent\tT_Percent\tDepth\n"
            )

            if input_file.endswith(".gz"):
                with gzip.open(input_file, "rt") as f:
                    for line in f:
                        stats = parse_mpileup_line(line, exon_info, total_rows)
                        if stats:
                            write_stats(out, stats)
            else:
                with open(input_file, "r") as f:
                    for line in f:
                        stats = parse_mpileup_line(line, exon_info, total_rows)
                        if stats:
                            write_stats(out, stats)

        if summary_file:
            generate_summary(output_file, summary_file)

        return True
    except IOError as e:
        print(f"I/O error processing file: {str(e)}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error processing file: {str(e)}", file=sys.stderr)
        return False


def generate_summary(stats_file, summary_file, threshold=10):
    """Generate a summary of polymorphic positions."""
    try:
        with open(stats_file, "r") as stats, open(summary_file, "w") as summary:
            try:
                next(stats)  # Skip header
            except StopIteration:
                print(
                    f"Warning: Stats file {stats_file} appears to be empty",
                    file=sys.stderr,
                )
                return False

            for line in stats:
                try:
                    fields = line.strip().split("\t")
                    pos, ref = fields[0], fields[1]
                    match_percent = int(fields[2])
                    mismatch_percent = int(fields[3])
                    ins_percent, del_percent = int(fields[4]), int(fields[5])
                    a_percent, g_percent, c_percent, t_percent = map(int, fields[6:10])
                    depth = int(fields[10])

                    if (
                        mismatch_percent >= threshold
                        or ins_percent >= threshold
                        or del_percent >= threshold
                    ):
                        summary.write(
                            f"(1-based) Position:{pos}, Reference Base={ref}\n"
                        )
                        summary.write(f"Aligned Read Count:{depth}\n")
                        summary.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                        summary.write(
                            f"{match_percent}\t{mismatch_percent}\t{ins_percent}\t{del_percent}\t"
                            f"{a_percent}\t{g_percent}\t{c_percent}\t{t_percent}\n\n"
                        )
                except (IndexError, ValueError) as e:
                    print(
                        f"Error processing line in stats file: {line.strip()}",
                        file=sys.stderr,
                    )
                    print(f"Exception: {str(e)}", file=sys.stderr)
                    continue
        return True
    except IOError as e:
        print(f"I/O error generating summary: {str(e)}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error generating summary: {str(e)}", file=sys.stderr)
        return False


def main():
    """Parse command line arguments and process mpileup file."""
    parser = argparse.ArgumentParser(
        description="Calculate nucleotide statistics from mpileup format"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Mpileup file to analyze (gzipped or uncompressed)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file for nucleotide stats",
    )
    parser.add_argument(
        "-s",
        "--summary",
        help="Output file for SNP summary (optional)",
        default=None,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=10,
        help="Threshold percentage for considering a position polymorphic (default: 10)",
    )

    args = parser.parse_args()

    success = process_mpileup_file(args.input, args.output, args.summary)

    if success:
        print(f"Nucleotide stats written to: {args.output}")

        if args.summary:
            print(f"SNP summary written to: {args.summary}")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
