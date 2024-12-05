#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
from AnalyzeAbo_Common import *
from PairwiseAlignmentResult import *
from BloodGroupStats import *
from os.path import join, isdir
from os import makedirs

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


# Parameters.
# If this proportion of reads match the reference, we assume the query sequence matches the reference.
ReadCutoffAgreementPercent = 80
AlleleCutoffAgreementPercent = 95


def findPolymorphismsInReads(referenceFileName, readsFileName, outputDirectoryName):
    """
    Find polymorphisms in the provided reads.

    Parameters:
    - reference_file_name (str): The name of the reference file.
    - reads_file_name (str): The name of the reads file.
    - output_directory (str): The directory where output files will be saved.
    """
    print(("Looking for polymorphisms in these reads:" + str(readsFileName)))
    alignSequencesAgainstReference(
        referenceFileName, readsFileName, outputDirectoryName
    )
    sequenceStats = analyzeReadAlignment(outputDirectoryName)
    writeReadAlignmentSpreadsheet(outputDirectoryName, sequenceStats)
    findReadPolymorphisms(outputDirectoryName, sequenceStats)
    analyzeChosenPolymorphicPositions(
        referenceFileName, outputDirectoryName, sequenceStats
    )


# Perform BW Alignment.  Align all reads against the Reference.


def batchABOAllelesAgainstReference(
    referenceFileName, allelesFileName, outputDirectoryName
):
    """
    Align ABO alleles against the reference and analyze the alignment.

    Parameters:
    - reference_file_name (str): The name of the reference file.
    - alleles_file_name (str): The name of the alleles file.
    - output_directory (str): The directory where output files will be saved.
    """
    alignSequencesAgainstReference(
        referenceFileName, allelesFileName, outputDirectoryName
    )
    bloodGroupStats = analyzeAlleleAlignment(outputDirectoryName)
    writeAlleleAlignmentSpreadsheet(outputDirectoryName, bloodGroupStats)
    findInterestingAllelePolymorphisms(outputDirectoryName, bloodGroupStats)


def alignSequencesAgainstReference(
    referenceLocation, alleleFileLocation, outputDirectory
):
    """
    Align sequences against a reference and analyze the alignment.

    Parameters:
    - reference_location (str): The location of the reference file.
    - allele_file_location (str): The location of the allele file.
    - output_directory (str): The directory where output files will be saved.
    """
    print("Aligning ABO allele sequences against the reference.")

    alignmentSubdir = join(outputDirectory, "alignment")

    # Align Alleles against Reference
    if not isdir(alignmentSubdir):
        makedirs(alignmentSubdir)

    # Part 1 Index the Reference
    try:
        # Copy the reference sequence to the alignment directory. This is a complicated way to do it.
        newReferenceLocation = join(alignmentSubdir, "AlignmentReference.fasta")
        refSequence = list(SeqIO.parse(referenceLocation, "fasta"))[0]
        refSequence.id = "AlignmentReference"
        sequenceWriter = createOutputFile(newReferenceLocation)
        SeqIO.write([refSequence], sequenceWriter, "fasta")
        sequenceWriter.close()

        # Index The Reference || BWA deprecated !!
        # cmd = ("bwa index " + newReferenceLocation)
        # os.system(cmd)

    except Exception:
        print(
            "Exception indexing alignment reference. Is minimap2 installed? Folder writing permission issue?"
        )
        raise

    # Part 2 Align
    try:
        # align | sam->bam | sort
        tempAlignmentName = join(alignmentSubdir, "alignment")
        alignmentOutputName = tempAlignmentName + ".bam"
        statOutputName = tempAlignmentName + ".samtools.flagstat"
        flagstatOutputName = tempAlignmentName + ".samtools.stats"
        minimap2Args = "-ax map-ont -t 8"

        ## BWA replaced with minimap2 for long-read alignment
        # bwaMemArgs = "-t 4 -x ont2d"
        # cmd = ("bwa mem " +
        #        bwaMemArgs + " " +
        #        newReferenceLocation + " " +
        #        alleleFileLocation +
        #        " | samtools view  -Sb - | samtools sort -o "
        #        + alignmentOutputName)
        cmd = (
            "minimap2 "
            + minimap2Args
            + " "
            + newReferenceLocation
            + " "
            + alleleFileLocation
            + " | samtools view  -Sb - | samtools sort -o "
            + alignmentOutputName
        )
        # print ('alignment command:\n' + cmd)

        os.system(cmd)

    except Exception:
        print(
            "Exception aligning reads against reference. Are minimap2 and samtools installed?"
        )
        raise

    # Part 3 Index Alignment and flagstat
    try:
        cmd = "samtools index " + alignmentOutputName
        # print ('alignment index command:\n' + cmd)
        os.system(cmd)
        # print ('index command:\n' + cmd)

        ## Samtools stats
        cmd = "samtools flagstat " + alignmentOutputName + ">" + flagstatOutputName
        os.system(cmd)

        cmd = "samtools stats " + alignmentOutputName + ">" + statOutputName
        os.system(cmd)

    except Exception:
        print("Exception indexing alignment reference. Is samtools installed?")
        raise


def analyzeReadAlignment(outputDirectory):
    """
    Analyze the read alignment and calculate statistics.

    Parameters:
    - output_directory (str): The directory where output files will be saved.

    Returns:
    - list: A list of sequence statistics.
    """
    print("Parse the alignment, finding allele patterns...")

    alignmentSubdir = join(outputDirectory, "alignment")

    # Load up the Alignment Reference file, we'll need it.
    alignmentReferenceFileName = join(alignmentSubdir, "AlignmentReference.fasta")
    alignmentRef = list(SeqIO.parse(alignmentReferenceFileName, "fasta"))[0]

    # Open the bam file
    bamfile = pysam.AlignmentFile(join(alignmentSubdir, "alignment.bam"), "rb")

    # I can use an array of BloodGroupColumn objects for this purpose.
    # This helps me keep track of matches and mismatches and indels
    sequenceStats = []

    # Iterate the reference sequence column by column.
    pileupIterator = bamfile.pileup(alignmentRef.id)

    for pileupColumnIndex, pileupColumn in enumerate(pileupIterator):
        # referencePosition = pileupColumn.reference_pos
        referenceBase = alignmentRef[pileupColumn.reference_pos].upper()
        # alignedCount = pileupColumn.nsegments

        currentColumnStats = BloodGroupColumn(referenceBase)

        currentColumnStats.alignedReadCount = pileupColumn.nsegments

        # Iterate the Reads at this position
        for alignedReadRowIndex, pileupRead in enumerate(pileupColumn.pileups):
            queryPosition = pileupRead.query_position

            # This piled up read is a deletion.
            if pileupRead.is_del == 1:
                currentColumnStats.deleteBase()

                # myBloodGroupStats.processDeletion(bloodGroup, pileupColumnIndex)
            # This one is an insertion.
            elif pileupRead.indel > 0:
                # Inertions can be multiple bases, I think.  I'm only grabbing one base here.
                # Seems to me that one base is probably enough.

                currentBase = pileupRead.alignment.query_sequence[queryPosition].upper()
                currentColumnStats.insertBases(currentBase)
                # myBloodGroupStats.processInsertion(bloodGroup, pileupColumnIndex, currentBase)

            elif queryPosition is not None:
                currentBase = pileupRead.alignment.query_sequence[queryPosition].upper()

                # I guess I should check for this but I don't really know what this means. Maybe it never happens.
                if pileupRead.is_refskip:
                    print(
                        (
                            "This read is a refskip, i dont know what that means:"
                            + pileupRead.alignment.query_name
                        )
                    )
                    raise Exception(
                        "This read is a refskip, i dont know what that means:"
                        + pileupRead.alignment.query_name
                    )
                # else this means we have a base aligned at this position for this read.

                else:
                    if currentBase == referenceBase:
                        currentColumnStats.matchBase(currentBase)
                        # myBloodGroupStats.processMatch(bloodGroup, pileupColumnIndex, currentBase)

                    else:
                        currentColumnStats.mismatchBase(currentBase)
                        # myBloodGroupStats.processMismatch(bloodGroup, pileupColumnIndex,currentBase)

            else:
                print("I did not get an aligned position.")
                raise Exception("I don't know what to do with this aligned Sequence.")

        sequenceStats.append(currentColumnStats)

    return sequenceStats


def writeReadAlignmentSpreadsheet(outputDirectory, sequenceStats):
    """
    Analyze the read alignment and calculate statistics.

    Parameters:
    - output_directory (str): The directory where output files will be saved.

    Returns:
    - list: A list of sequence statistics.
    """
    print("Writing Read Alignment Stats to output file.")

    alignmentSummaryFile = createOutputFile(
        join(outputDirectory, "ReadAlignmentSpreadsheet.csv")
    )
    # headers:
    alignmentSummaryFile.write(
        "Ref_Position_1based,Ref_Base"
        + ",Match_Percent,Mismatch_Percent,Insertion_Percent,Deletion_Percent"
        + ",A_Percent,G_Percent,C_Percent,T_Percent"
        + "\n"
    )

    for referencePosition, currentColumnStats in enumerate(sequenceStats):
        alignmentSummaryFile.write(
            str(referencePosition + 1)
            + ","
            + str(currentColumnStats.referenceBase)
            + ","
            + str(currentColumnStats.getMatchPercent())
            + ","
            + str(currentColumnStats.getMismatchPercent())
            + ","
            + str(currentColumnStats.getInsertionPercent())
            + ","
            + str(currentColumnStats.getDeletionPercent())
            + ","
            + str(currentColumnStats.getNucleotidePercentage("A"))
            + ","
            + str(currentColumnStats.getNucleotidePercentage("G"))
            + ","
            + str(currentColumnStats.getNucleotidePercentage("C"))
            + ","
            + str(currentColumnStats.getNucleotidePercentage("T"))
            + "\n"
        )

    alignmentSummaryFile.close()


def printBasePolymorphisms(
    alignmentSummaryFile, referencePosition0Based, currentSequenceStats
):
    alignmentSummaryFile.write(
        "\n(1-based) Position:"
        + str(referencePosition0Based + 1)
        + ", Reference Base="
        + currentSequenceStats.referenceBase
        + "\n"
    )

    alignmentSummaryFile.write(
        "Aligned Read Count:" + str(currentSequenceStats.alignedReadCount) + "\n"
    )

    # Header, blood group is in the columns.
    alignmentSummaryFile.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
    alignmentSummaryFile.write(str(int(currentSequenceStats.getMatchPercent())) + "\t")
    alignmentSummaryFile.write(
        str(int(currentSequenceStats.getMismatchPercent())) + "\t"
    )
    alignmentSummaryFile.write(
        str(int(currentSequenceStats.getInsertionPercent())) + "\t"
    )
    alignmentSummaryFile.write(
        str(int(currentSequenceStats.getDeletionPercent())) + "\t"
    )
    # loop Nucleotides
    for nucleotide in "A", "G", "C", "T":
        nucleotidePercentage = currentSequenceStats.getNucleotidePercentage(nucleotide)
        percentageInt = int(nucleotidePercentage)
        alignmentSummaryFile.write(str(percentageInt) + "\t")
    alignmentSummaryFile.write("\n")


def findReadPolymorphisms(outputDirectory, sequenceStats):
    print("Searching for the Important Polymorphisms in our ABO reads.")

    alignmentSummaryFile = createOutputFile(
        join(outputDirectory, "ABOReadPolymorphisms.txt")
    )

    # TODO: Print a read count. We want to know how many reads are aligned in here.

    for referencePosition0Based, currentSequenceStats in enumerate(sequenceStats):
        matchPercent = currentSequenceStats.getMatchPercent()

        if matchPercent > ReadCutoffAgreementPercent:
            # The sequences think we match. I don't care about this position.
            pass

        else:
            printBasePolymorphisms(
                alignmentSummaryFile, referencePosition0Based, currentSequenceStats
            )

            # Print the surrounding sequence
            # I'm using try catch for programming logic, there is certainly a smarter way to do it.
            try:
                baseCount = 10
                leftSequence = ""
                rightSequence = ""

                for i in range(0, baseCount):
                    leftSequence += sequenceStats[
                        referencePosition0Based - baseCount + i
                    ].referenceBase
                    rightSequence += sequenceStats[
                        referencePosition0Based + i
                    ].referenceBase

                pass

                alignmentSummaryFile.write(
                    leftSequence
                    + " "
                    + currentSequenceStats.referenceBase
                    + " "
                    + rightSequence
                    + "\n"
                )

            except Exception:
                print(
                    "I had an exception when trying to print the sequence surrounding an interesting ABO polymorphism."
                )
                print("I bet the reason is out of bounds, sequence is unknown.")
                print("I choose to just not respond right now.")

            # alignmentSummaryFile.write('\n')

    alignmentSummaryFile.close()


def analyzeChosenPolymorphicPositions(
    referenceFileName, outputDirectory, sequenceStats
):
    """
    Analyze chosen polymorphic positions.

    Parameters:
    - reference_file_name (str): The name of the reference file.
    - output_directory (str): The directory where output files will be saved.
    - sequence_stats (list): A list of sequence statistics.
    """
    print("Analyzing a few chosen polymorphic positions...")
    # Many hard coded values in here.  I think that's okay, these positions are set in stone for modern humans.

    phenotypeOutputFile = createOutputFile(join(outputDirectory, "ABOPhenotype.txt"))

    alignmentRef = list(SeqIO.parse(referenceFileName, "fasta"))[0]

    # Exon 6 starts at genomic nucleotide 240
    # Exon 7 starts at genomic nucelotide 375

    # Exon 6.
    if len(alignmentRef) == 135:
        phenotypeOutputFile.write("Exon 6:\n")

        # phenotypeOutputFile.write('Genomic position: 261')
        phenotypeOutputFile.write("\nExon 6 position(1-based): 22\n")
        phenotypeOutputFile.write("G nucleotide: A or B blood type.\n")
        phenotypeOutputFile.write("Deletion    : O blood type(O1).")
        printBasePolymorphisms(
            phenotypeOutputFile, 21, sequenceStats[21]
        )  # pass the 0-based

        # TODO: Make up some cutoffs.  If the numbers are close to 50/50, we can guess 2 alleles.

    # Exon 7
    elif len(alignmentRef) == 691:
        phenotypeOutputFile.write("Exon 7:\n")

        # phenotypeOutputFile.write('Genomic position: 796')
        phenotypeOutputFile.write("\nExon 7 position(1-based): 422\n")
        phenotypeOutputFile.write("A nucleotide: B blood type.\n")
        phenotypeOutputFile.write("C nucelotide: A or O blood type.")
        printBasePolymorphisms(phenotypeOutputFile, 421, sequenceStats[421])

        # phenotypeOutputFile.write('Genomic position: 802')
        phenotypeOutputFile.write("\nExon 7 position(1-based): 428\n")
        phenotypeOutputFile.write("A nucleotide: O blood type (O2).\n")
        phenotypeOutputFile.write("G nucelotide: A or B or O blood type.")
        printBasePolymorphisms(phenotypeOutputFile, 427, sequenceStats[427])

        # phenotypeOutputFile.write('Genomic position: 803')
        phenotypeOutputFile.write("\nExon 7 position(1-based): 429\n")
        phenotypeOutputFile.write("G nucleotide: A or O blood type.\n")
        phenotypeOutputFile.write("C nucelotide: B blood type.")
        printBasePolymorphisms(phenotypeOutputFile, 428, sequenceStats[428])

        # phenotypeOutputFile.write('Genomic position: 805')
        phenotypeOutputFile.write("\nExon 7 position(1-based): 431\n")
        phenotypeOutputFile.write("G nucleotide: O blood type (O3).\n")
        phenotypeOutputFile.write("T nucelotide: A or B or O blood type.")
        printBasePolymorphisms(phenotypeOutputFile, 430, sequenceStats[430])
    else:
        phenotypeOutputFile.write("This exon has length:" + str(len(alignmentRef)))
    phenotypeOutputFile.close()

    # Exon 7

    # TODO: Not sure how to implement this right now.
    # I think the important positions will just be hardcoded here.
    # I will have to check if we're working with Ex 6 or 7. Somehow.


def analyzeAlleleAlignment(outputDirectory):
    """
    Analyze the alignment of ABO alleles.

    Parameters:
    - output_directory (str): The directory where output files will be saved.

    Returns:
    - list: A list of blood group statistics.
    """

    print("Parse the alignment, finding allele patterns...")

    alignmentSubdir = join(outputDirectory, "alignment")

    # Load up the Alignment Reference file, we'll need it.
    alignmentReferenceFileName = join(alignmentSubdir, "AlignmentReference.fasta")
    alignmentRef = list(SeqIO.parse(alignmentReferenceFileName, "fasta"))[0]

    # Open the bam file
    bamfile = pysam.AlignmentFile(join(alignmentSubdir, "alignment.bam"), "rb")

    myBloodGroupStats = BloodGroupStats(alignmentRef.seq)

    # Iterate the reference sequence column by column.
    pileupIterator = bamfile.pileup(alignmentRef.id)

    for pileupColumnIndex, pileupColumn in enumerate(pileupIterator):
        # referencePosition = pileupColumn.reference_pos
        referenceBase = alignmentRef[pileupColumn.reference_pos].upper()
        # alignedCount = pileupColumn.nsegments

        # Iterate the Reads at this position
        for alignedReadRowIndex, pileupRead in enumerate(pileupColumn.pileups):
            # wait, what is the sequene name here?
            queryAlleleName = pileupRead.alignment.query_name
            bloodGroup = getPhenotype(queryAlleleName)

            queryPosition = pileupRead.query_position

            # This piled up read is a deletion.
            if pileupRead.is_del == 1:
                myBloodGroupStats.processDeletion(bloodGroup, pileupColumnIndex)
            # This one is an insertion.
            elif pileupRead.indel > 0:
                # Inertions can be multiple bases, I think.  I'm only grabbing one base here.
                # Seems to me that one base is probably enough.
                currentBase = pileupRead.alignment.query_sequence[queryPosition].upper()
                myBloodGroupStats.processInsertion(
                    bloodGroup, pileupColumnIndex, currentBase
                )

            elif queryPosition is not None:
                currentBase = pileupRead.alignment.query_sequence[queryPosition].upper()

                if pileupRead.is_refskip:
                    print(
                        (
                            "This read is a refskip, i dont know what that means:"
                            + pileupRead.alignment.query_name
                        )
                    )
                    raise Exception(
                        "This read is a refskip, i dont know what that means:"
                        + pileupRead.alignment.query_name
                    )
                # else this means we have a base aligned at this position for this read.

                else:
                    if currentBase == referenceBase:
                        myBloodGroupStats.processMatch(
                            bloodGroup, pileupColumnIndex, currentBase
                        )

                    else:
                        myBloodGroupStats.processMismatch(
                            bloodGroup, pileupColumnIndex, currentBase
                        )

            else:
                print(
                    (
                        "I did not get an alignemd position for this allele:"
                        + str(queryAlleleName)
                    )
                )
                raise Exception("I don't know what to do with this aligned abo allele.")

    return myBloodGroupStats


def writeAlleleAlignmentSpreadsheet(outputDirectory, myBloodGroupStats):
    """
    Write the allele alignment spreadsheet.

    Parameters:
    - output_directory (str): The directory where output files will be saved.
    - blood_group_stats (list): A list of blood group statistics.
    """
    print("Writing alignment Stats to output file.")

    alignmentSummaryFile = createOutputFile(
        join(outputDirectory, "AlignmentSpreadsheet.csv")
    )
    # headers:
    alignmentSummaryFile.write(
        "Ref_Position_1based,Ref_Base"
        # +_',Aligned_Sequence_Count'
        + ",Most_Common_Base_A,Most_Common_Base_B,Most_Common_Base_O"
        + ",A_Match_Percent,A_Mismatch_Percent,A_Insertion_Percent,A_Deletion_Percent"
        + ",B_Match_Percent,B_Mismatch_Percent,B_Insertion_Percent,B_Deletion_Percent"
        + ",O_Match_Percent,O_Mismatch_Percent,O_Insertion_Percent,O_Deletion_Percent"
        + "\n"
    )

    for referencePosition, referenceBase in enumerate(
        myBloodGroupStats.referenceSequence
    ):
        alignmentSummaryFile.write(
            str(referencePosition + 1)
            + ","
            + str(referenceBase)
            # + ',' + str(alignedCount)
            + ","
            + str(myBloodGroupStats.getMostCommonBase("A", referencePosition))
            + ","
            + str(myBloodGroupStats.getMostCommonBase("B", referencePosition))
            + ","
            + str(myBloodGroupStats.getMostCommonBase("O", referencePosition))
            + ","
            + str(myBloodGroupStats.getMatchPercent("A", referencePosition))
            + ","
            + str(myBloodGroupStats.getMismatchPercent("A", referencePosition))
            + ","
            + str(myBloodGroupStats.getInsertionPercent("A", referencePosition))
            + ","
            + str(myBloodGroupStats.getDeletionPercent("A", referencePosition))
            + ","
            + str(myBloodGroupStats.getMatchPercent("B", referencePosition))
            + ","
            + str(myBloodGroupStats.getMismatchPercent("B", referencePosition))
            + ","
            + str(myBloodGroupStats.getInsertionPercent("B", referencePosition))
            + ","
            + str(myBloodGroupStats.getDeletionPercent("B", referencePosition))
            + ","
            + str(myBloodGroupStats.getMatchPercent("O", referencePosition))
            + ","
            + str(myBloodGroupStats.getMismatchPercent("O", referencePosition))
            + ","
            + str(myBloodGroupStats.getInsertionPercent("O", referencePosition))
            + ","
            + str(myBloodGroupStats.getDeletionPercent("O", referencePosition))
            + "\n"
        )

    alignmentSummaryFile.close()


def findInterestingAllelePolymorphisms(outputDirectory, myBloodGroupStats):
    """
    Find interesting polymorphisms in the allele sequences.

    Parameters:
    - output_directory (str): The directory where output files will be saved.
    - blood_group_stats (list): A list of blood group statistics.
    """
    print("Searching for the Interesting Polymorphisms.")

    # If we match 95% of the time, we can consider it 100%
    # matchPercentageCutoff = 95

    alignmentSummaryFile = createOutputFile(
        join(outputDirectory, "Interesting_Polymorphisms.txt")
    )

    for referencePosition0Based, referenceBase in enumerate(
        myBloodGroupStats.referenceSequence
    ):
        aMatchPercent = myBloodGroupStats.getMatchPercent("A", referencePosition0Based)
        bMatchPercent = myBloodGroupStats.getMatchPercent("B", referencePosition0Based)
        oMatchPercent = myBloodGroupStats.getMatchPercent("O", referencePosition0Based)

        # If all blood groups agree 100%, this is not an interesting base.
        # if(aMatchPercent == 100 and bMatchPercent == 100 and oMatchPercent == 100):
        if (
            aMatchPercent > AlleleCutoffAgreementPercent
            and bMatchPercent > AlleleCutoffAgreementPercent
            and oMatchPercent > AlleleCutoffAgreementPercent
        ):
            pass

        else:
            alignmentSummaryFile.write(
                "\nPosition (1-based):"
                + str(referencePosition0Based + 1)
                + ", ReferenceBase="
                + referenceBase
                + "\n"
            )

            # Header, blood group is in the columns.
            alignmentSummaryFile.write("\tMat\tMis\tIns\tDel\tA\tG\tC\tT\n")

            # Loop Blood Group / phenotypes
            for bloodGroup in ("A", "B", "O"):
                alignmentSummaryFile.write(bloodGroup + "\t")

                alignmentSummaryFile.write(
                    str(
                        int(
                            myBloodGroupStats.getMatchPercent(
                                bloodGroup, referencePosition0Based
                            )
                        )
                    )
                    + "\t"
                )
                alignmentSummaryFile.write(
                    str(
                        int(
                            myBloodGroupStats.getMismatchPercent(
                                bloodGroup, referencePosition0Based
                            )
                        )
                    )
                    + "\t"
                )
                alignmentSummaryFile.write(
                    str(
                        int(
                            myBloodGroupStats.getInsertionPercent(
                                bloodGroup, referencePosition0Based
                            )
                        )
                    )
                    + "\t"
                )
                alignmentSummaryFile.write(
                    str(
                        int(
                            myBloodGroupStats.getDeletionPercent(
                                bloodGroup, referencePosition0Based
                            )
                        )
                    )
                    + "\t"
                )

                # loop Nucleotides
                for nucleotide in ("A", "G", "C", "T"):
                    nucleotidePercentage = myBloodGroupStats.getNucleotidePercentage(
                        bloodGroup, nucleotide, referencePosition0Based
                    )
                    percentageInt = int(nucleotidePercentage)

                    alignmentSummaryFile.write(str(percentageInt) + "\t")

                alignmentSummaryFile.write("\n")

            alignmentSummaryFile.write("\n")

    alignmentSummaryFile.close()


def getPhenotype(currentSeqID):
    # The character at [4] is the blood group. This method is not complicated.

    if currentSeqID.startswith("ABO"):
        # print('This read starts with ABO.')
        pass
    else:
        print("doesnt start with abo.")
        raise Exception("I want the sequence IDs to start with ABO")

    currentPhenotype = "??"

    # Strange logic to decide the ABO types
    if currentSeqID[4] == "O":
        currentPhenotype = "O"

    elif currentSeqID[4] == "A":
        currentPhenotype = "A"

    elif currentSeqID[4] == "B":
        currentPhenotype = "B"

    else:
        print(("I don't know what phenotype is this allele:" + str(currentSeqID)))
        raise Exception("I want the sequence IDs to start with ABO")

    return currentPhenotype


def scaleValues(inputValues):
    # scale from 0:1
    scaledData = []

    oldMinimum = min(inputValues)
    oldMaximum = max(inputValues)

    # Calculated based on a formula from here:
    # http://stats.stackexchange.com/questions/25894/changing-the-scale-of-a-variable-to-0-100

    if oldMaximum == oldMinimum:
        return inputValues
    else:
        for currentValue in inputValues:
            scaledData.append(
                ((1 - 0) / (oldMaximum - oldMinimum)) * (currentValue - oldMaximum) + 1
            )

    return scaledData
