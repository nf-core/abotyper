#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import logging
import os
import re
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, Set, Tuple, Union

__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2024, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = ["Fredrick Mobegi", "Benedict Matern", "Mathijs Groeneweg"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Production"


"""
SAMtools Pileup Statistics Calculator - Improved Version

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

This script processes SAMtools mpileup files to calculate alignment frequencies
for all nucleotides and indels per reference position.
"""

class ExonType(Enum):
    """Enumeration for exon types based on reference length."""
    EXON6 = "exon6"
    EXON7 = "exon7"
    UNKNOWN = "unknown"


@dataclass
class BaseCount:
    """Data class for nucleotide base counts."""
    A: int = 0
    G: int = 0
    C: int = 0
    T: int = 0

    def total(self) -> int:
        """Return total count of all bases."""
        return self.A + self.G + self.C + self.T

    def get_count(self, base: str) -> int:
        """Get count for a specific base."""
        return getattr(self, base.upper(), 0)

    def set_count(self, base: str, count: int) -> None:
        """Set count for a specific base."""
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
        """Create a stats object with zero values."""
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
    """Main class for processing mpileup files."""

    # Class constants
    EXON6_LENGTH_RANGE = (130, 140)
    EXON7_LENGTH_RANGE = (800, 830)
    LOW_COVERAGE_THRESHOLD = 200
    KEY_DIAGNOSTIC_POSITIONS = {431, 687}
    EXON6_INDEL_POSITION = 22  # Position where indels are expected in exon6

    NUCLEOTIDES = {"A", "G", "C", "T"}

    def __init__(self, loglevel: str = "INFO"):
        """Initialize the processor with logging configuration."""
        self._setup_logging(loglevel)
        self.logger = logging.getLogger(__name__)

    def _setup_logging(self, loglevel: str) -> None:
        """Setup logging configuration."""
        logging.basicConfig(
            level=getattr(logging, loglevel.upper()),
            format='[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

    def determine_exon_type(self, ref_length: int) -> ExonType:
        """Determine exon type based on reference length."""
        if self.EXON6_LENGTH_RANGE[0] <= ref_length <= self.EXON6_LENGTH_RANGE[1]:
            return ExonType.EXON6
        elif self.EXON7_LENGTH_RANGE[0] <= ref_length <= self.EXON7_LENGTH_RANGE[1]:
            return ExonType.EXON7
        else:
            return ExonType.UNKNOWN

    def _parse_bases_string(self, read_bases: str, ref_base: str) -> ParsedBases:
        """Parse the read bases string from mpileup format."""
        # Clean up read bases string
        read_bases = re.sub(r'\^.', '', read_bases)  # Remove start markers
        read_bases = read_bases.replace('$', '')  # Remove end markers

        parsed = ParsedBases()
        i = 0

        while i < len(read_bases):
            char = read_bases[i]

            if char in '.,':
                # Match to reference
                parsed.matches += 1
                parsed.base_counts.set_count(ref_base, parsed.base_counts.get_count(ref_base) + 1)
                i += 1

            elif char.upper() in self.NUCLEOTIDES:
                # Mismatch
                parsed.mismatches += 1
                base = char.upper()
                parsed.base_counts.set_count(base, parsed.base_counts.get_count(base) + 1)
                i += 1

            elif char == '+':
                # Insertion
                i += 1
                ins_len_str = ""
                while i < len(read_bases) and read_bases[i].isdigit():
                    ins_len_str += read_bases[i]
                    i += 1
                ins_len = int(ins_len_str) if ins_len_str else 0
                i += ins_len  # Skip the inserted sequence
                parsed.insertions += 1

            elif char == '-':
                # Deletion
                i += 1
                del_len_str = ""
                while i < len(read_bases) and read_bases[i].isdigit():
                    del_len_str += read_bases[i]
                    i += 1
                del_len = int(del_len_str) if del_len_str else 0
                i += del_len  # Skip the deleted sequence
                parsed.deletions += 1

            elif char == '*':
                # Deletion (alternative format)
                parsed.deletions += 1
                i += 1

            else:
                # Unknown character, skip
                self.logger.debug(f"Unknown character in read bases: {char}")
                i += 1

        return parsed

    def _should_include_indels(
        self,
        exon_type: ExonType,
        pos: int,
        coverage: int,
        total_rows: int
    ) -> bool:
        """
        Determine if indels should be included in calculations for this position.

        Except for 2 positions (exon6pos22 or c.261, and exon7pos687 or c.1061),
        indels in other ABO associated SNV positions are most likely sequencing errors
        and must be handled correctly to avoid mistyping.
        """
        # Always include indels at key diagnostic positions
        if pos in self.KEY_DIAGNOSTIC_POSITIONS or pos == self.EXON6_INDEL_POSITION:
            return True

        # For low coverage, be more conservative about including indels
        if coverage < self.LOW_COVERAGE_THRESHOLD:
            if exon_type == ExonType.EXON6:
                return False  # Exclude indels for most exon6 positions with low coverage
            elif total_rows > 140:
                return False  # Exclude indels for non-diagnostic positions in long sequences

        return True

    def _calculate_percentages(
        self,
        parsed: ParsedBases,
        ref_base: str,
        include_indels: bool
    ) -> Dict[str, int]:
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
            # When ignoring indels, use only ATGC counts as denominator
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

            # Ensure ATGC percentages sum to 100% by adjusting reference base
            atgc_sum = sum(base_percentages.values())
            if atgc_sum != 100 and atgc_sum > 0:
                diff = 100 - atgc_sum
                ref_key = f"{ref_base}_percent"
                base_percentages[ref_key] += diff

            insertion_percent = 0
            deletion_percent = 0

        # Calculate match and mismatch percentages
        match_percent = base_percentages[f"{ref_base}_percent"]
        mismatch_percent = sum(
            base_percentages[f"{base}_percent"]
            for base in self.NUCLEOTIDES
            if base != ref_base
        )

        return {
            'match_percent': match_percent,
            'mismatch_percent': mismatch_percent,
            'insertion_percent': insertion_percent,
            'deletion_percent': deletion_percent,
            **base_percentages
        }

    def parse_mpileup_line(
        self,
        line: str,
        exon_type: ExonType = ExonType.UNKNOWN,
        total_rows: int = 0
    ) -> Optional[PositionStats]:
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

            # Parse the read bases
            parsed = self._parse_bases_string(read_bases, ref_base)

            # Determine if we should include indels
            include_indels = self._should_include_indels(exon_type, pos, coverage, total_rows)

            # Calculate percentages
            percentages = self._calculate_percentages(parsed, ref_base, include_indels)

            return PositionStats(
                pos=pos,
                ref_base=ref_base,
                depth=coverage,
                **percentages
            )

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
        """Write the header line to the output file."""
        header = [
            "Ref_Position_1based", "Ref_Base", "Match_Percent", "Mismatch_Percent",
            "Insertion_Percent", "Deletion_Percent", "A_Percent", "G_Percent",
            "C_Percent", "T_Percent", "Depth"
        ]
        output_file.write('\t'.join(header) + '\n')

    def write_stats_line(self, output_file, stats: PositionStats) -> None:
        """Write a statistics line to the output file."""
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
        summary_file: Optional[Union[str, Path]] = None
    ) -> bool:
        """Process an mpileup file and output nucleotide statistics."""
        input_path = Path(input_file)
        output_path = Path(output_file)

        # Validate input file
        if not input_path.exists():
            self.logger.error(f"Input file does not exist: {input_path}")
            return False

        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # First pass: analyze file structure
            unique_positions, total_rows = self._get_unique_positions(input_path)

            if not unique_positions:
                self.logger.error("No valid positions found in input file")
                return False

            # Determine exon type
            ref_length = len(unique_positions)
            exon_type = self.determine_exon_type(ref_length)
            self.logger.info(f"Detected exon type: {exon_type.value} (reference length: {ref_length})")

            # Second pass: process the file
            self.logger.info("Processing mpileup data...")
            stats_processed = 0

            open_func = gzip.open if input_path.suffix == '.gz' else open
            mode = 'rt' if input_path.suffix == '.gz' else 'r'

            with open_func(input_path, mode) as infile, open(output_path, 'w') as outfile:
                self.write_stats_header(outfile)

                for line_num, line in enumerate(infile, 1):
                    if line_num % 10000 == 0:
                        self.logger.debug(f"Processed {line_num} lines...")

                    stats = self.parse_mpileup_line(line, exon_type, total_rows)
                    if stats:
                        self.write_stats_line(outfile, stats)
                        stats_processed += 1

            self.logger.info(f"Successfully processed {stats_processed} positions")

            # Generate summary if requested
            if summary_file:
                self.generate_summary(output_path, Path(summary_file))

            return True

        except Exception as e:
            self.logger.error(f"Error processing mpileup file: {e}")
            return False

    def generate_summary(
        self,
        stats_file: Path,
        summary_file: Path,
        threshold: int = 10
    ) -> bool:
        """Generate a summary of polymorphic positions."""
        self.logger.info(f"Generating summary with threshold {threshold}%")

        try:
            with open(stats_file, 'r') as stats, open(summary_file, 'w') as summary:
                next(stats)  # Skip header

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

                        # Check if position is polymorphic
                        if (mismatch_percent >= threshold or
                            ins_percent >= threshold or
                            del_percent >= threshold):

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
    %(prog)s -i input.mpileup -o output.stats -t 5 --verbose
        """
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Mpileup file to analyze (gzipped or uncompressed)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file for nucleotide statistics"
    )
    parser.add_argument(
        "-s", "--summary",
        help="Output file for polymorphic position summary (optional)"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=int,
        default=10,
        help="Threshold percentage for considering a position polymorphic (default: 10)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )

    args = parser.parse_args()

    # Set logging level
    loglevel = "DEBUG" if args.verbose else "INFO"

    # Initialize processor
    processor = PileupProcessor(loglevel=loglevel)

    # Process the file
    success = processor.process_mpileup_file(
        input_file=args.input,
        output_file=args.output,
        summary_file=args.summary
    )

    if success:
        processor.logger.info(f"✓ Nucleotide statistics written to: {args.output}")
        if args.summary:
            processor.logger.info(f"✓ Polymorphic position summary written to: {args.summary}")
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
