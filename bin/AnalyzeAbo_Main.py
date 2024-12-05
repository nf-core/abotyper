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

# # This file is part of abo-analysis.
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

from AnalyzeAbo import *
import getopt
import sys
import os

SoftwareVersion = "abo-analysis Version 1.0"


def usage():
    print(
        "usage:\n"
        + "\tThis script is written for python 3.11\n"
        + "\tI haven't written the usage tutorial yet.  Oops.  Do this now please."
    )


# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.

    # Todo: I should pass in an Exon 6 and Exon 7 reference.
    # I don't need the A,B,O alleles, since I'm using logic to parse the phenotype from the allele name

    global referenceFileName
    global readsFileName
    global allelesFileName
    # global allelesAFileName
    # global allelesBFileName
    # global allelesOFileName
    global outputDirectoryName
    global analysisType

    referenceFileName = None
    readsFileName = None
    allelesFileName = None
    # allelesAFileName        = None
    # allelesBFileName        = None
    # allelesOFileName        = None
    outputDirectoryName = None
    analysisType = None

    if len(sys.argv) < 3:
        print("I don't think you have enough arguments.\n")
        usage()
        raise ()
        return False

    # getopt.getopt(..) is a function for parsing args the way smart people do it
    # For More info, use google or
    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],  # TODO Fix this, i don't accept abo alleles anymore
            "hvr:R:a:o:A:B:O:t:e6:e7:",
            [
                "help",
                "version",
                "reference=",
                "reads=",
                "alleles=",
                "output=",
                "alleles-a=",
                "alleles-b=",
                "alleles-o=",
                "analysis-type=",
                "exon-6=",
                "exon-7",
            ],
        )

        for opt, arg in opts:

            if opt in ("-h", "--help"):
                print(SoftwareVersion)
                usage()
                return False

            elif opt in ("-v", "--version"):
                print(SoftwareVersion)
                return False

            elif opt in ("-i", "--alleles"):
                allelesFileName = arg

            elif opt in ("-r", "--reference"):
                referenceFileName = arg

            elif opt in ("-R", "--reads"):
                readsFileName = arg

            elif opt in ("-o", "--output"):
                outputDirectoryName = arg

            elif opt in ("e6", "--exon-6"):
                print("YOU GAVE ME EXON 6:" + arg)

            # elif opt in ("-A", "--alleles-a"):
            #    allelesAFileName = arg

            # elif opt in ("-B", "--alleles-b"):
            #    allelesBFileName = arg

            # elif opt in ("-O", "--alleles-o"):
            #    allelesOFileName = arg

            elif opt in ("-t", "--analysis-type"):
                analysisType = arg

            else:
                print("Unknown commandline option: " + opt)
                raise ()

    except getopt.GetoptError as errorMessage:
        print("Something seems wrong with your commandline parameters.")
        print(errorMessage)
        usage()
        return False

    print("Reference File:" + str(referenceFileName))
    print("Reads File:" + str(readsFileName))
    print("Alleles File:" + str(allelesFileName))
    print("Output Directory:" + str(outputDirectoryName))

    # Quick sanity check.
    if len(referenceFileName) < 4:
        print("referenceFileName is too short:" + str(referenceFileName))
        return False
    # TODO: I don't feel like writing sanity checks rignt now.

    # If we're doing allele analysis we need:
    # fasta with alleles
    # exon 6 and7 referencx
    # If read analysis
    # need ex 6 and 7.

    # if(len(allelesFileName) < 4):
    #    print('allelesFileName is too short:' + str(allelesFileName))
    #    return False
    if len(outputDirectoryName) < 4:
        print("Output directory is too short:" + str(outputDirectoryName))
        return False

    if not os.path.isdir(outputDirectoryName):
        os.makedirs(outputDirectoryName)

    return True


if __name__ == "__main__":

    try:
        if readArgs():

            print("Starting to Analyze ABO")

            if analysisType == "READS":

                findPolymorphismsInReads(
                    referenceFileName, readsFileName, outputDirectoryName
                )

            elif analysisType == "ALLELES":
                # compareIndividualABOAllelesToReference(referenceFileName, allelesFileName, outputDirectoryName)
                batchABOAllelesAgainstReference(
                    referenceFileName, allelesFileName, outputDirectoryName
                )
            # findPolymorphismsInReads(referenceFileName, readsFileName, outputDirectoryName)

            else:
                print(
                    'Please provide an analysis type. Use the "-t" or "--analysis-type" with a type of "READS" or "ALLELES" '
                )
                print("I give up.")

            print("Done analyzing ABO.")
        else:
            print(
                "\nI'm giving up because I was not satisfied with your commandline arguments."
            )

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print("Fatal problem during read extraction:")
        print(sys.exc_info())
        raise ()
