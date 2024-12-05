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

# This file is part of abo-analysis.
#
# abo-analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# abo-analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with abo-analysis. If not, see <http://www.gnu.org/licenses/>.


import os

# from Bio.Seq import Seq
from Bio import SeqIO
from os.path import split


# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
    resultsOutput = open(outputfileName, "w")
    return resultsOutput


def loadInputRecords(recordFileName):
    if ".fasta" == recordFileName[-6:] or ".fa" == recordFileName[-3:]:
        FileOutputFormat = "fasta"
    elif ".fastq" == recordFileName[-6:] or ".fq" == recordFileName[-3:]:
        FileOutputFormat = "fastq"
    else:
        FileOutputFormat = "UnknownFormat"

    parsedInputReads = SeqIO.parse(recordFileName, FileOutputFormat)
    alignedSequences = enumerate(parsedInputReads)

    return alignedSequences
