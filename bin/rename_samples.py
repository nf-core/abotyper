#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sample Renaming Utility

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

This script processes ABO pipeline export files and laboratory deobfuscation files
to rename samples for MatchPoint export compatibility.
"""

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd

__author__ = "Fredrick Mobegi"
__copyright__ = "Copyright 2024, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = ["Fredrick Mobegi", "Benedict Matern", "Mathijs Groeneweg"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Production"


class SampleRenamer:
    """Main class for sample renaming operations."""

    # Constants
    DATE_FORMAT = "%Y_%m_%d"
    REQUIRED_COLUMNS = ["Acc#", "Patient Name"]
    MIN_READ_THRESHOLD = 20
    SAMPLE_ID_PATTERN = r"(.+)_barcode\d+$"

    def __init__(self, loglevel: str = "INFO"):
        """Initialize the renamer with logging configuration."""
        self._setup_logging(loglevel)
        self.logger = logging.getLogger(__name__)

    def _setup_logging(self, loglevel: str) -> None:
        """Setup logging configuration."""
        logging.basicConfig(
            level=getattr(logging, loglevel.upper()),
            format='[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

    def read_renaming_file(self, deobfuscation_file: Path) -> pd.DataFrame:
        """
        Read the Excel file from soft PCR HLA exported file and select specific columns by name.

        Args:
            deobfuscation_file: Path to the Excel deobfuscation file

        Returns:
            DataFrame with selected columns
        """
        self.logger.info(f"Reading deobfuscation file: {deobfuscation_file}")

        if not deobfuscation_file.exists():
            raise FileNotFoundError(f"Deobfuscation file not found: {deobfuscation_file}")

        try:
            df = pd.read_excel(
                deobfuscation_file,
                index_col=None,
                na_values=["NA"],
                usecols=self.REQUIRED_COLUMNS
            )

            self.logger.info(f"Successfully read {len(df)} rows from deobfuscation file")
            self.logger.debug(f"Columns found: {list(df.columns)}")

            # Validate required columns are present
            missing_cols = [col for col in self.REQUIRED_COLUMNS if col not in df.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")

            return df

        except Exception as e:
            self.logger.error(f"Error reading deobfuscation file: {e}")
            raise

    def preprocess_renaming_file(self, renaming_file: pd.DataFrame) -> pd.DataFrame:
        """
        Preprocess the renaming file by renaming columns, extracting the first value
        from comma-separated Grid_number column, and filtering invalid entries.

        Args:
            renaming_file: Raw DataFrame from deobfuscation file

        Returns:
            Processed DataFrame ready for merging
        """
        self.logger.info("Preprocessing renaming file...")

        # Create a copy to avoid modifying original
        df = renaming_file.copy()

        # Rename columns
        df = df.rename(columns={
            "Acc#": "Sample ID",
            "Patient Name": "Grid_number"
        })

        if "Grid_number" not in df.columns:
            self.logger.warning("Grid_number column not found after renaming")
            return df

        # Process Grid_number column
        original_count = len(df)

        # Split comma-separated values, keeping only the first part
        df["Grid_number"] = df["Grid_number"].apply(
            lambda x: str(x).split(",")[0] if pd.notna(x) and "," in str(x) else str(x)
        )

        # Ensure Grid_number is string type
        df["Grid_number"] = df["Grid_number"].astype(str)

        # Filter out rows where Grid_number starts with a letter (invalid samples)
        df = df[~df["Grid_number"].str.contains("^[a-zA-Z]", regex=True, na=False)]

        # Remove rows with 'nan' string values
        df = df[df["Grid_number"] != "nan"]

        filtered_count = len(df)
        removed_count = original_count - filtered_count

        self.logger.info(f"Preprocessed renaming file: {filtered_count} valid entries, {removed_count} removed")

        return df

    def apply_regex_pattern(self, series: pd.Series, pattern: str) -> pd.Series:
        """
        Apply regex pattern to extract base sample names from Sample ID column.

        Args:
            series: Pandas Series containing sample IDs
            pattern: Regex pattern to apply

        Returns:
            Series with pattern applied
        """
        import re

        def extract_match(x):
            if pd.isna(x):
                return x
            match = re.match(pattern, str(x))
            return match.group(1) if match else x

        return series.apply(extract_match)

    def read_final_export_file(self, final_export_file: Path) -> pd.DataFrame:
        """
        Read the final export file from ABO pipeline into a DataFrame.

        Args:
            final_export_file: Path to the CSV export file

        Returns:
            DataFrame with export data
        """
        self.logger.info(f"Reading final export file: {final_export_file}")

        if not final_export_file.exists():
            raise FileNotFoundError(f"Final export file not found: {final_export_file}")

        try:
            df = pd.read_csv(final_export_file, sep=",")
            self.logger.info(f"Successfully read {len(df)} rows from export file")
            self.logger.debug(f"Columns found: {list(df.columns)}")

            if "Sample ID" not in df.columns:
                raise ValueError("Sample ID column not found in export file")

            return df

        except Exception as e:
            self.logger.error(f"Error reading final export file: {e}")
            raise

    def merge_dataframes(self, final_export: pd.DataFrame, renaming_file: pd.DataFrame) -> pd.DataFrame:
        """
        Left join final_export with renaming data using Sample ID as the key.

        Args:
            final_export: DataFrame from ABO pipeline export
            renaming_file: Processed DataFrame from deobfuscation file

        Returns:
            Merged DataFrame
        """
        self.logger.info("Merging export and renaming dataframes...")

        before_count = len(final_export)
        merged_df = pd.merge(final_export, renaming_file, on="Sample ID", how="left")
        after_count = len(merged_df)

        # Count successful matches
        matched_count = merged_df["Grid_number"].notna().sum()
        unmatched_count = before_count - matched_count

        self.logger.info(f"Merge complete: {matched_count} samples matched, {unmatched_count} unmatched")

        if unmatched_count > 0:
            self.logger.warning(f"{unmatched_count} samples could not be matched with patient IDs")

        return merged_df

    def reorder_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Reorder columns to put Grid_number first.

        Args:
            df: Input DataFrame

        Returns:
            DataFrame with reordered columns
        """
        if "Grid_number" not in df.columns:
            self.logger.warning("Grid_number column not found for reordering")
            return df

        col_order = ["Grid_number"] + [col for col in df.columns if col != "Grid_number"]
        return df[col_order]

    def rename_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Rename columns for final output format.

        Args:
            df: Input DataFrame

        Returns:
            DataFrame with renamed columns
        """
        return df.rename(columns={
            "Sample ID": "SequencingAcc#",
            "Grid_number": "Sample ID"
        })

    def create_copy_without_sequencing_acc(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Create a copy without the sequencing accession number.

        Args:
            df: Input DataFrame with SequencingAcc# column

        Returns:
            DataFrame without SequencingAcc# column
        """
        df_copy = df.copy()
        if "SequencingAcc#" in df_copy.columns:
            df_copy.drop("SequencingAcc#", axis=1, inplace=True)
        return df_copy

    def filter_and_clean_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter data based on read count threshold and clean invalid entries.

        Args:
            df: Input DataFrame

        Returns:
            Cleaned and filtered DataFrame
        """
        self.logger.info("Filtering and cleaning data...")

        # Convert literal 'nan' strings to actual NaN
        df = df.replace("nan", np.nan)

        original_count = len(df)

        # Filter by read count if column exists
        if "#Reads" in df.columns:
            before_filter = len(df)
            df = df[df["#Reads"] >= self.MIN_READ_THRESHOLD]
            after_filter = len(df)

            self.logger.info(f"Filtered by read count (>={self.MIN_READ_THRESHOLD}): "
                           f"{after_filter} samples kept, {before_filter - after_filter} removed")

            # Remove the #Reads column as it's no longer needed
            df = df.drop(columns=["#Reads"])

        # Remove rows where Sample ID is NaN
        df = df[df["Sample ID"].notna()]

        final_count = len(df)
        total_removed = original_count - final_count

        self.logger.info(f"Data cleaning complete: {final_count} samples remaining, "
                        f"{total_removed} total samples removed")

        return df

    def write_output_files(
        self,
        df_with_acc: pd.DataFrame,
        df_without_acc: pd.DataFrame,
        output_dir: Path = Path(".")
    ) -> Tuple[Path, Path]:
        """
        Write DataFrames to output files with date suffix.

        Args:
            df_with_acc: DataFrame including sequencing accession numbers
            df_without_acc: DataFrame without sequencing accession numbers
            output_dir: Directory to write files to

        Returns:
            Tuple of (path_with_acc, path_without_acc)
        """
        self.logger.info("Writing output files...")

        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate date suffix
        date_suffix = datetime.now().strftime(self.DATE_FORMAT)

        # Define output paths
        path_with_acc = output_dir / f"MatchPointExport_with_sequencingAcc_{date_suffix}.txt"
        path_without_acc = output_dir / f"MatchPointExport_{date_suffix}.txt"

        try:
            # Write files
            df_with_acc.to_csv(
                path_with_acc,
                index=False,
                encoding="utf-8",
                lineterminator="\n"
            )

            df_without_acc.to_csv(
                path_without_acc,
                index=False,
                encoding="utf-8",
                lineterminator="\n"
            )

            self.logger.info(f"✓ Files written successfully:")
            self.logger.info(f"  - With accession: {path_with_acc} ({len(df_with_acc)} rows)")
            self.logger.info(f"  - Without accession: {path_without_acc} ({len(df_without_acc)} rows)")

            return path_with_acc, path_without_acc

        except Exception as e:
            self.logger.error(f"Error writing output files: {e}")
            raise

    def process_samples(
        self,
        final_export_file: Path,
        deobfuscation_file: Path,
        output_dir: Optional[Path] = None
    ) -> Tuple[Path, Path]:
        """
        Main processing function to rename samples and generate output files.

        Args:
            final_export_file: Path to ABO pipeline export file
            deobfuscation_file: Path to laboratory deobfuscation Excel file
            output_dir: Optional output directory (defaults to current directory)

        Returns:
            Tuple of output file paths
        """
        if output_dir is None:
            output_dir = Path(".")

        self.logger.info("Starting sample renaming process...")

        try:
            # Step 1: Read deobfuscation file
            renaming_df = self.read_renaming_file(deobfuscation_file)

            # Step 2: Preprocess renaming file
            renaming_df = self.preprocess_renaming_file(renaming_df)

            # Step 3: Apply regex pattern to extract base sample names
            renaming_df["Sample ID"] = self.apply_regex_pattern(
                renaming_df["Sample ID"], self.SAMPLE_ID_PATTERN
            )

            # Step 4: Read final export file
            export_df = self.read_final_export_file(final_export_file)

            # Step 5: Apply same regex pattern to export file
            export_df["Sample ID"] = self.apply_regex_pattern(
                export_df["Sample ID"], self.SAMPLE_ID_PATTERN
            )

            # Step 6: Merge dataframes
            merged_df = self.merge_dataframes(export_df, renaming_df)

            # Step 7: Reorder columns
            merged_df = self.reorder_columns(merged_df)

            # Step 8: Convert Grid_number to string
            if "Grid_number" in merged_df.columns:
                merged_df["Grid_number"] = merged_df["Grid_number"].astype(str)

            # Step 9: Rename columns
            merged_df = self.rename_columns(merged_df)

            # Step 10: Create copy without sequencing accession
            df_without_acc = self.create_copy_without_sequencing_acc(merged_df)

            # Step 11: Filter and clean the data (only for the version without accession)
            df_without_acc = self.filter_and_clean_data(df_without_acc)

            # Step 12: Write output files
            output_paths = self.write_output_files(merged_df, df_without_acc, output_dir)

            self.logger.info("✓ Sample renaming process completed successfully")
            return output_paths

        except Exception as e:
            self.logger.error(f"Error in sample processing: {e}")
            raise


def main():
    """Main function to parse command line arguments and process samples."""
    parser = argparse.ArgumentParser(
        description="Rename samples from ABO pipeline export for MatchPoint compatibility.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s -f final_export.csv -d deobfuscation.xlsx
    %(prog)s -f final_export.csv -d deobfuscation.xlsx -o output/
    %(prog)s -f final_export.csv -d deobfuscation.xlsx --verbose
        """
    )

    parser.add_argument(
        "-f", "--final-export",
        required=True,
        help="Path to the final export CSV file from ABO pipeline"
    )
    parser.add_argument(
        "-d", "--deobfuscation",
        required=True,
        help="Path to the Excel deobfuscation file from laboratory"
    )
    parser.add_argument(
        "-o", "--output",
        default=".",
        help="Output directory for renamed files (default: current directory)"
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

    # Initialize renamer
    renamer = SampleRenamer(loglevel=loglevel)

    try:
        # Convert string paths to Path objects
        final_export_path = Path(args.final_export)
        deobfuscation_path = Path(args.deobfuscation)
        output_path = Path(args.output)

        # Process samples
        output_files = renamer.process_samples(
            final_export_file=final_export_path,
            deobfuscation_file=deobfuscation_path,
            output_dir=output_path
        )

        print(f"\n✓ Successfully renamed samples!")
        print(f"Output files:")
        print(f"  - {output_files[0]}")
        print(f"  - {output_files[1]}")

        sys.exit(0)

    except Exception as e:
        renamer.logger.error(f"Processing failed: {e}")
        print(f"\n✗ Error: {e}")
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
