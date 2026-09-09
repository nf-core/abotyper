#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import logging
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Set, Tuple, Union

from abo_panel import load_panel, Panel

__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2024-2025, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = ["Fredrick Mobegi", "Benedict Matern", "Mathijs Groeneweg",
               "Claude Sonnet 5 (v2.0.0 rewrite for panel-driven indel handling)"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Production"


"""
SAMtools Pileup Statistics Calculator — v2.0.0

This file is part of the nf-core/abotyper pipeline "https://github.com/fmobegi/nf-core-abotyper".

CHANGES IN v2.0.0
------------------
  - The set of positions where indels are the DIAGNOSTIC variant (previously
    the hardcoded KEY_DIAGNOSTIC_POSITIONS = {431, 687} / EXON6_INDEL_POSITION
    = 22) is now loaded from the external variant panel (--panel), via
    Panel.indel_diagnostic_positions(). Add a new indel-diagnostic marker by
    editing the panel file, not this script.
  - Exon-type classification by reference length is now OPTIONAL context
    used only to decide the low-coverage indel-inclusion heuristic; it is
    no longer required for indel-diagnostic-position lookup. This means the
    script degrades gracefully when run against the new combined exon2-7
    reference (where the old length-range heuristic doesn't apply) --
    pass --exon-label explicitly (e.g. "combined") or let it default to
    "combined", in which case the more permissive high-coverage rule is
    always used for non-diagnostic positions.

This pipeline is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

This script processes SAMtools mpileup files to calculate alignment
frequencies for all nucleotides and indels per reference position.
"""


@dataclass
class BaseCount:
    """Data class for nucleotide base counts."""
    A: int = 0
    G: int = 0
    C: int = 0
    T: int = 0

    def total(self) -> int:
        return self.A + self.G + self.C + self.T

    def get_count(self, base: str) -> int:
        return getattr(self, base.upper(), 0)

    def set_count(self, base: str, count: int) -> None:
        if hasattr(self, base.upper()):
            setattr(self, base.upper(), count)


@dataclass
class PositionStats:
    """Data class for position statistics."""
    pos: int
    ref_base: str
    match_percent: int = 0
    mismatch_percent: int = 0
    insertion_percent: int = 0
    deletion_percent: int = 0
    A_percent: int = 0
    G_percent: int = 0
    C_percent: int = 0
    T_percent: int = 0
    depth: int = 0

    @classmethod
    def create_zero_stats(cls, pos: int, ref_base: str) -> 'PositionStats':
        return cls(pos=pos, ref_base=ref_base)


@dataclass
class ParsedBases:
    """Data class for parsed base information from mpileup."""
    matches: int = 0
    mismatches: int = 0
    insertions: int = 0
    deletions: int = 0
    base_counts: BaseCount = field(default_factory=BaseCount)


class PileupProcessor:
    """Main class for processing mpileup files. v2.0.0: indel-diagnostic
    positions are panel-driven rather than hardcoded."""

    # Retained only as a fallback hint for --exon-label auto-detection
    # against v1.x-style short mini-amplicon files. Meaningless for the
    # new combined exon2-7 reference -- pass --exon-label "combined"
    # (the default) in that case.
    LEGACY_EXON6_LENGTH_RANGE = (130, 140)
    LEGACY_EXON7_LENGTH_RANGE = (800, 830)
    LOW_COVERAGE_THRESHOLD = 200

    NUCLEOTIDES = {"A", "G", "C", "T"}

    def __init__(self, panel: Panel, exon_label: str = "combined", loglevel: str = "INFO"):
        """
        Args:
            panel: loaded variant panel (see abo_panel.load_panel)
            exon_label: "combined" (default -- treat as the new single
                exon2-7 amplicon; indel-diagnostic positions are the union
                across all exons, resolved via amplicon_pos), or a specific
                legacy label like "Exon 6"/"Exon 7" (resolved via
                legacy_exon6_pos/legacy_exon7_pos) for v1.x-style
                mini-amplicon files.
        """
        self.panel = panel
        self.exon_label = exon_label
        self._setup_logging(loglevel)
        self.logger = logging.getLogger(__name__)
        self.indel_diagnostic_positions = self._resolve_indel_positions()
        self.logger.info(
            f"Loaded {len(self.indel_diagnostic_positions)} indel-diagnostic "
            f"positions from panel for exon_label='{exon_label}'"
        )

    def _setup_logging(self, loglevel: str) -> None:
        logging.basicConfig(
            level=getattr(logging, loglevel.upper()),
            format='[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

    def _resolve_indel_positions(self) -> Set[int]:
        if self.exon_label == "combined":
            positions: Set[int] = set()
            for exon in self.panel.exons():
                positions.update(self.panel.indel_diagnostic_positions(exon))
            return positions
        return set(self.panel.indel_diagnostic_positions(self.exon_label))

    def determine_legacy_exon_hint(self, ref_length: int) -> str:
        """Best-effort legacy hint only; not required for correctness."""
        if self.LEGACY_EXON6_LENGTH_RANGE[0] <= ref_length <= self.LEGACY_EXON6_LENGTH_RANGE[1]:
            return "Exon 6"
        if self.LEGACY_EXON7_LENGTH_RANGE[0] <= ref_length <= self.LEGACY_EXON7_LENGTH_RANGE[1]:
            return "Exon 7"
        return "combined"

    def _parse_bases_string(self, read_bases: str, ref_base: str) -> ParsedBases:
        """Parse the read bases string from mpileup format."""
        read_bases = re.sub(r'\^.', '', read_bases)
        read_bases = read_bases.replace('$', '')

        parsed = ParsedBases()
        i = 0

        while i < len(read_bases):
            char = read_bases[i]

            if char in '.,':
                parsed.matches += 1
                parsed.base_counts.set_count(ref_base, parsed.base_counts.get_count(ref_base) + 1)
                i += 1

            elif char.upper() in self.NUCLEOTIDES:
                parsed.mismatches += 1
                base = char.upper()
                parsed.base_counts.set_count(base, parsed.base_counts.get_count(base) + 1)
                i += 1

            elif char == '+':
                i += 1
                ins_len_str = ""
                while i < len(read_bases) and read_bases[i].isdigit():
                    ins_len_str += read_bases[i]
                    i += 1
                ins_len = int(ins_len_str) if ins_len_str else 0
                i += ins_len
                parsed.insertions += 1

            elif char == '-':
                i += 1
                del_len_str = ""
                while i < len(read_bases) and read_bases[i].isdigit():
                    del_len_str += read_bases[i]
                    i += 1
                del_len = int(del_len_str) if del_len_str else 0
                i += del_len
                parsed.deletions += 1

            elif char == '*':
                parsed.deletions += 1
                i += 1

            else:
                logging.getLogger(__name__).debug(f"Unknown character in read bases: {char}")
                i += 1

        return parsed

    def _should_include_indels(self, pos: int, coverage: int, exon_hint: str) -> bool:
        """
        Determine if indels should be included in calculations for this
        position. Panel-driven: always include at any position flagged as
        indel-diagnostic in the variant panel (was: hardcoded {431, 687, 22}).
        """
        if pos in self.indel_diagnostic_positions:
            return True

        if coverage < self.LOW_COVERAGE_THRESHOLD:
            if exon_hint == "Exon 6":
                return False
            # For "Exon 7" / "combined": exclude indels at non-diagnostic
            # positions in long low-coverage stretches, matching v1.x
            # behaviour for exon7 (total_rows > 140).
            return False if exon_hint in ("Exon 7", "combined") else True

        return True

    def _calculate_percentages(self, parsed: ParsedBases, ref_base: str, include_indels: bool) -> Dict[str, int]:
        """Calculate percentages for each nucleotide and indel type."""
        total_events = parsed.matches + parsed.mismatches + parsed.insertions + parsed.deletions
        total_nucleotides = parsed.matches + parsed.mismatches

        if total_events == 0:
            return {
                'match_percent': 0, 'mismatch_percent': 0,
                'insertion_percent': 0, 'deletion_percent': 0,
                'A_percent': 0, 'G_percent': 0, 'C_percent': 0, 'T_percent': 0
            }

        if include_indels:
            denominator = total_events
            base_percentages = {
                f"{base}_percent": int((parsed.base_counts.get_count(base) / denominator) * 100)
                for base in self.NUCLEOTIDES
            }
            insertion_percent = int((parsed.insertions / denominator) * 100)
            deletion_percent = int((parsed.deletions / denominator) * 100)
        else:
            denominator = total_nucleotides
            if denominator == 0:
                return {
                    'match_percent': 0, 'mismatch_percent': 0,
                    'insertion_percent': 0, 'deletion_percent': 0,
                    'A_percent': 0, 'G_percent': 0, 'C_percent': 0, 'T_percent': 0
                }
            base_percentages = {
                f"{base}_percent": int((parsed.base_counts.get_count(base) / denominator) * 100)
                for base in self.NUCLEOTIDES
            }
            atgc_sum = sum(base_percentages.values())
            if atgc_sum != 100 and atgc_sum > 0:
                diff = 100 - atgc_sum
                ref_key = f"{ref_base}_percent"
                base_percentages[ref_key] += diff
            insertion_percent = 0
            deletion_percent = 0

        match_percent = base_percentages[f"{ref_base}_percent"]
        mismatch_percent = sum(
            base_percentages[f"{base}_percent"] for base in self.NUCLEOTIDES if base != ref_base
        )

        return {
            'match_percent': match_percent,
            'mismatch_percent': mismatch_percent,
            'insertion_percent': insertion_percent,
            'deletion_percent': deletion_percent,
            **base_percentages
        }

    def parse_mpileup_line(self, line: str, exon_hint: str = "combined") -> Optional[PositionStats]:
        """Parse a single line from a mpileup file and calculate nucleotide statistics."""
        try:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                self.logger.warning(f"Insufficient fields in line: {line.strip()}")
                return None

            pos = int(fields[1])
            ref_base = fields[2].upper()
            coverage = int(fields[3])
            read_bases = fields[4]

            if coverage == 0:
                return PositionStats.create_zero_stats(pos, ref_base)

            parsed = self._parse_bases_string(read_bases, ref_base)
            include_indels = self._should_include_indels(pos, coverage, exon_hint)
            percentages = self._calculate_percentages(parsed, ref_base, include_indels)

            return PositionStats(pos=pos, ref_base=ref_base, depth=coverage, **percentages)

        except (ValueError, IndexError) as e:
            self.logger.error(f"Error parsing line: {line.strip()[:100]}... - {e}")
            return None
        except Exception as e:
            self.logger.error(f"Unexpected error parsing line: {e}")
            return None

    def _get_unique_positions(self, input_file: Path) -> Tuple[Set[int], int]:
        """Get unique positions and total row count from input file."""
        unique_positions = set()
        total_rows = 0

        self.logger.info(f"Analyzing file structure: {input_file}")

        try:
            open_func = gzip.open if input_file.suffix == '.gz' else open
            mode = 'rt' if input_file.suffix == '.gz' else 'r'

            with open_func(input_file, mode) as f:
                for line_num, line in enumerate(f, 1):
                    total_rows += 1
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        try:
                            unique_positions.add(int(fields[1]))
                        except ValueError:
                            self.logger.warning(f"Invalid position at line {line_num}: {fields[1]}")

        except Exception as e:
            self.logger.error(f"Error analyzing file structure: {e}")
            raise

        self.logger.info(f"Found {len(unique_positions)} unique positions in {total_rows} rows")
        return unique_positions, total_rows

    def write_stats_header(self, output_file) -> None:
        header = [
            "Ref_Position_1based", "Ref_Base", "Match_Percent", "Mismatch_Percent",
            "Insertion_Percent", "Deletion_Percent", "A_Percent", "G_Percent",
            "C_Percent", "T_Percent", "Depth"
        ]
        output_file.write('\t'.join(header) + '\n')

    def write_stats_line(self, output_file, stats: PositionStats) -> None:
        line = [
            str(stats.pos), stats.ref_base, str(stats.match_percent),
            str(stats.mismatch_percent), str(stats.insertion_percent),
            str(stats.deletion_percent), str(stats.A_percent),
            str(stats.G_percent), str(stats.C_percent),
            str(stats.T_percent), str(stats.depth)
        ]
        output_file.write('\t'.join(line) + '\n')

    def process_mpileup_file(
        self,
        input_file: Union[str, Path],
        output_file: Union[str, Path],
        summary_file: Optional[Union[str, Path]] = None,
        exon_label_override: Optional[str] = None,
    ) -> bool:
        """Process an mpileup file and output nucleotide statistics."""
        input_path = Path(input_file)
        output_path = Path(output_file)

        if not input_path.exists():
            self.logger.error(f"Input file does not exist: {input_path}")
            return False

        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            unique_positions, total_rows = self._get_unique_positions(input_path)

            if not unique_positions:
                self.logger.error("No valid positions found in input file")
                return False

            ref_length = len(unique_positions)
            exon_hint = exon_label_override or self.exon_label
            if exon_hint == "combined":
                # Provide an informational best-guess only; does not affect
                # indel-diagnostic lookups (panel-driven), only the
                # low-coverage indel-inclusion heuristic.
                guessed = self.determine_legacy_exon_hint(ref_length)
                self.logger.info(
                    f"Reference length {ref_length} -> legacy-style hint "
                    f"'{guessed}' (informational only; indel-diagnostic "
                    f"positions come from the panel regardless)"
                )
                exon_hint = guessed if guessed != "combined" else "combined"

            self.logger.info(f"Processing with exon_hint='{exon_hint}' "
                              f"(reference length: {ref_length})")

            stats_processed = 0
            open_func = gzip.open if input_path.suffix == '.gz' else open
            mode = 'rt' if input_path.suffix == '.gz' else 'r'

            with open_func(input_path, mode) as infile, open(output_path, 'w') as outfile:
                self.write_stats_header(outfile)

                for line_num, line in enumerate(infile, 1):
                    if line_num % 10000 == 0:
                        self.logger.debug(f"Processed {line_num} lines...")

                    stats = self.parse_mpileup_line(line, exon_hint)
                    if stats:
                        self.write_stats_line(outfile, stats)
                        stats_processed += 1

            self.logger.info(f"Successfully processed {stats_processed} positions")

            if summary_file:
                self.generate_summary(output_path, Path(summary_file))

            return True

        except Exception as e:
            self.logger.error(f"Error processing mpileup file: {e}")
            return False

    def generate_summary(self, stats_file: Path, summary_file: Path, threshold: int = 10) -> bool:
        """Generate a summary of polymorphic positions."""
        self.logger.info(f"Generating summary with threshold {threshold}%")

        try:
            with open(stats_file, 'r') as stats, open(summary_file, 'w') as summary:
                next(stats)
                polymorphic_count = 0

                for line_num, line in enumerate(stats, 1):
                    try:
                        fields = line.strip().split('\t')
                        if len(fields) < 11:
                            self.logger.warning(f"Insufficient fields in stats line {line_num}")
                            continue

                        pos, ref = fields[0], fields[1]
                        match_percent = int(fields[2])
                        mismatch_percent = int(fields[3])
                        ins_percent, del_percent = int(fields[4]), int(fields[5])
                        a_percent, g_percent, c_percent, t_percent = map(int, fields[6:10])
                        depth = int(fields[10])

                        if (mismatch_percent >= threshold or ins_percent >= threshold
                                or del_percent >= threshold):
                            summary.write(f"(1-based) Position:{pos}, Reference Base={ref}\n")
                            summary.write(f"Aligned Read Count:{depth}\n")
                            summary.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                            summary.write(
                                f"{match_percent}\t{mismatch_percent}\t{ins_percent}\t{del_percent}\t"
                                f"{a_percent}\t{g_percent}\t{c_percent}\t{t_percent}\n\n"
                            )
                            polymorphic_count += 1

                    except (IndexError, ValueError) as e:
                        self.logger.error(f"Error processing stats line {line_num}: {e}")
                        continue

                self.logger.info(f"Found {polymorphic_count} polymorphic positions")
                return True

        except Exception as e:
            self.logger.error(f"Error generating summary: {e}")
            return False


def main():
    """Main function to parse command line arguments and process mpileup file."""
    parser = argparse.ArgumentParser(
        description="Calculate nucleotide statistics from mpileup format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s -i input.mpileup -o output.stats
    %(prog)s -i input.mpileup.gz -o output.stats -s summary.txt
    %(prog)s -i input.mpileup -o output.stats --exon-label "Exon 6" --panel abo_variant_panel.yaml
        """
    )

    parser.add_argument("-i", "--input", required=True,
                         help="Mpileup file to analyze (gzipped or uncompressed)")
    parser.add_argument("-o", "--output", required=True,
                         help="Output file for nucleotide statistics")
    parser.add_argument("-s", "--summary",
                         help="Output file for polymorphic position summary (optional)")
    parser.add_argument("-t", "--threshold", type=int, default=10,
                         help="Threshold percentage for considering a position polymorphic (default: %(default)s)")
    parser.add_argument("--panel", default="abo_variant_panel.yaml",
                         help="Path to the variant panel file (.yaml/.json/.csv/.tsv). Default: %(default)s")
    parser.add_argument("--exon-label", default="combined",
                         help='Which panel exon this mpileup file covers, e.g. "Exon 6", '
                              '"Exon 7", or "combined" (default) for the new exon2-7 amplicon '
                              "spanning multiple exons in one file.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()
    loglevel = "DEBUG" if args.verbose else "INFO"

    try:
        panel = load_panel(args.panel)
    except Exception as exc:
        print(f"! CRITICAL: could not load variant panel '{args.panel}': {exc}")
        sys.exit(1)

    processor = PileupProcessor(panel=panel, exon_label=args.exon_label, loglevel=loglevel)

    success = processor.process_mpileup_file(
        input_file=args.input,
        output_file=args.output,
        summary_file=args.summary,
    )

    if success:
        processor.logger.info(f"[OK] Nucleotide statistics written to: {args.output}")
        if args.summary:
            processor.logger.info(f"[OK] Polymorphic position summary written to: {args.summary}")
        sys.exit(0)
    else:
        processor.logger.error("Processing failed")
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n! Script interrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"! CRITICAL ERROR: Unhandled exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
