#!/usr/bin/env python3
import gzip
import re
import argparse
import sys
import os


def parse_mpileup_line(line):
    """Parse a single line from a mpileup file and calculate nucleotide statistics."""
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

        if total_events == 0:
            return create_zero_stats(pos, ref_base)

        def calc_percent(count):
            return int((count / total_events) * 100)

        base_percentages = {
            base: calc_percent(count) for base, count in base_counts.items()
        }
        insertion_percent = calc_percent(insertions)
        deletion_percent = calc_percent(deletions)

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


def process_mpileup(mpileup_file, output_file):
    """Process a mpileup file and output nucleotide statistics."""
    try:
        with open(output_file, "w") as out_f:
            out_f.write(
                "Ref_Position_1based\tRef_Base\tMatch_Percent\tMismatch_Percent\t"
                "Insertion_Percent\tDeletion_Percent\tA_Percent\tG_Percent\tC_Percent\tT_Percent\tDepth\n"
            )

            open_func = gzip.open if mpileup_file.endswith(".gz") else open
            mode = "rt" if mpileup_file.endswith(".gz") else "r"

            with open_func(mpileup_file, mode) as f:
                for line_num, line in enumerate(f, 1):
                    try:
                        stats = parse_mpileup_line(line)
                        if stats:
                            out_f.write(
                                f"{stats['pos']}\t{stats['ref_base']}\t"
                                f"{stats['match_percent']}\t{stats['mismatch_percent']}\t"
                                f"{stats['insertion_percent']}\t{stats['deletion_percent']}\t"
                                f"{stats['A_percent']}\t{stats['G_percent']}\t{stats['C_percent']}\t{stats['T_percent']}\t"
                                f"{stats['depth']}\n"
                            )
                    except Exception as e:
                        print(f"Error on line {line_num}: {str(e)}", file=sys.stderr)
                        continue
    except Exception as e:
        print(f"Error processing file: {str(e)}", file=sys.stderr)
        sys.exit(1)


def generate_summary(stats_file, summary_file, threshold=10):
    try:
        with open(stats_file, "r") as stats, open(summary_file, "w") as summary:
            next(stats)

            for line in stats:
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
                    summary.write(f"(1-based) Position:{pos}, Reference Base={ref}\n")
                    summary.write(f"Aligned Read Count:{depth}\n")
                    summary.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                    summary.write(
                        f"{match_percent}\t{mismatch_percent}\t{ins_percent}\t{del_percent}\t"
                        f"{a_percent}\t{g_percent}\t{c_percent}\t{t_percent}\t\n"
                    )

                    context = get_sequence_context(stats_file, int(pos))
                    summary.write(f"{context}\n\n")

    except Exception as e:
        print(f"Error generating summary: {str(e)}", file=sys.stderr)
        sys.exit(1)


def get_sequence_context(stats_file, pos, window=10):
    context = ["N"] * (2 * window + 1)
    with open(stats_file, "r") as stats:
        next(stats)
        for line in stats:
            fields = line.strip().split("\t")
            current_pos = int(fields[0])
            if abs(current_pos - pos) <= window:
                context[current_pos - pos + window] = fields[1]

    center_index = window
    before = "".join(context[:center_index])
    center = context[center_index]
    after = "".join(context[center_index + 1 :])
    return f"{before} {center} {after}"


def main():
    parser = argparse.ArgumentParser(
        description="Analyze mpileup file for nucleotide statistics."
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
        "-s", "--summary", help="Output file for SNP summary (optional)", default=None
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=10,
        help="Threshold percentage for considering a position polymorphic (default: 10)",
    )

    args = parser.parse_args()

    process_mpileup(args.input, args.output)

    print(f"Nucleotide stats written to: {args.output}")

    if args.summary:
        generate_summary(args.output, args.summary, args.threshold)
        print(f"SNP summary written to: {args.summary}")


if __name__ == "__main__":
    main()
