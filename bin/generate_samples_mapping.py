#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Benedict Matern, Fredrick Mobegi"
__copyright__ = "Copyright 2024, ABO blood group typing using third-generation sequencing (TGS) technology"
__credits__ = ["Fredrick Mobegi", "Benedict Matern", "Mathijs Groeneweg"]
__license__ = "GPL"
__version__ = "0.2.0"
__maintainer__ = "Fredrick Mobegi"
__email__ = "fredrick.mobegi@health.wa.gov.au"
__status__ = "Development"


"""
This script takes a filepath and identifies the _barcode* and exon* parts.
It uses these to create a new filename.
The alignment.samtools* files are copied into this new name.
This enables their identification by multiQC.
"""
import sys
import os
import argparse
import shutil


def generate_config(directory):
    """
    Generate MultiQC samples mapping and copy files.

    Args:
        directory (str): Directory containing the files for MultiQC mapping.
    """
    file_mapping = {}
    for root, _, files in os.walk(directory):
        for file in files:
            if file.startswith("alignment.samtools.") and (
                file.endswith(".stats") or file.endswith(".flagstat")
            ):
                filepath = os.path.join(root, file)
                parts = filepath.split("/")
                barcode_part = None
                exon_part = None
                for part in parts:
                    if "_barcode" in part:
                        barcode_part = part
                    if "exon" in part:
                        exon_part = part
                if barcode_part and exon_part:
                    sample_name = f"{barcode_part}_{exon_part}.{file.split('.')[1]}"
                    file_mapping[filepath] = sample_name
                    new_filename = os.path.basename(filepath).replace(
                        "alignment.samtools.", f"{sample_name}."
                    )
                    print(
                        "\033[1;32mCopying\033[0m",
                        filepath,
                        "\033[1;32minto\033[0m",
                        new_filename,
                    )
                    shutil.copy(filepath, os.path.join(os.getcwd(), new_filename))

    with open("samtools_files_mapping.yaml", "w", encoding="utf-8") as file:
        for filepath, sample_name in file_mapping.items():
            file.write(f"  {filepath}: {sample_name}\n")


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Generate MultiQC samples mapping and copy files."
    )
    parser.add_argument(
        "directory", help="Directory containing the files for MultiQC mapping"
    )
    args = parser.parse_args()
    generate_config(args.directory)
    print("Done renaming the SAMtools metrics files")


if __name__ == "__main__":
    main()
    sys.exit(0)
