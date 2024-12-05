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

# This class maintains statistics about allele polymorphisms relative to the standard ABO*A1.01.01.1 reference
# At each position, the aligned alleles show some polymorphims
# Polymorphisms are different depending on blood group.
# I want to cluster these polymorphisms by Blood Group.
# This means I can compare the expected polymorphisms, between A, B, and O alleles.


class BloodGroupStats:

    def __init__(self, referenceSequence):
        # An array. Each entry in this array is a position in the A1.01 reference.
        # There are statistics about the polymorphisms of A alles, B alleles, etc.
        self.aStats = []
        self.bStats = []
        self.oStats = []
        self.referenceSequence = referenceSequence

        for index, referenceBase in enumerate(referenceSequence):
            self.aStats.append(BloodGroupColumn(referenceBase))
            self.bStats.append(BloodGroupColumn(referenceBase))
            self.oStats.append(BloodGroupColumn(referenceBase))

    def getStatsForPhenotype(self, phenotype):
        if phenotype == "A":
            return self.aStats
        elif phenotype == "B":
            return self.bStats
        elif phenotype == "O":
            return self.oStats
        else:
            raise Exception("Unknown phenotype:" + str(phenotype))

    def getNucleotidePercentage(self, phenotype, nucleotide, position):
        return self.getStatsForPhenotype(phenotype)[position].getNucleotidePercentage(
            nucleotide
        )

    def getMostCommonBase(self, phenotype, position):
        return self.getStatsForPhenotype(phenotype)[position].getMostCommonBase()

    def getMatchPercent(self, phenotype, position):
        return self.getStatsForPhenotype(phenotype)[position].getMatchPercent()

    def getMismatchPercent(self, phenotype, position):
        return self.getStatsForPhenotype(phenotype)[position].getMismatchPercent()

    def getInsertionPercent(self, phenotype, position):
        return self.getStatsForPhenotype(phenotype)[position].getInsertionPercent()

    def getDeletionPercent(self, phenotype, position):
        return self.getStatsForPhenotype(phenotype)[position].getDeletionPercent()

    def processMismatch(self, phenotype, position, newBases):
        self.getStatsForPhenotype(phenotype)[position].mismatchBase(newBases)

    def processMatch(self, phenotype, position, newBases):
        self.getStatsForPhenotype(phenotype)[position].matchBase(newBases)

    def processDeletion(self, phenotype, position):
        self.getStatsForPhenotype(phenotype)[position].deleteBase()

    def processInsertion(self, phenotype, position, newBases):
        self.getStatsForPhenotype(phenotype)[position].mismatchBase(newBases)

        # print('I handled the insertion at position:' + str(position))


# This class represents a single column in an alignment, specific to a blood group.


class BloodGroupColumn:

    def __init__(self, refBase):
        self.referenceBase = refBase

        self.alignedReadCount = 0

        self.matchCount = 0
        self.mismatchCount = 0
        self.inCount = 0
        self.delCount = 0

        self.aCount = 0
        self.gCount = 0
        self.cCount = 0
        self.tCount = 0

    def getNucleotidePercentage(self, nucleotide):
        # I will try adding the delete count in here as well. The deletes aren't stored in the nucleotide counts
        # totalNucleotideCount = self.aCount + self.gCount + self.cCount + self.tCount
        totalNucleotideCount = (
            self.aCount + self.gCount + self.cCount + self.tCount + self.delCount
        )

        if totalNucleotideCount == 0:
            print(
                "I am returning a value of -1 because apparently all my nuclotide positions are equal to zero."
            )
            # raise Exception('what is going on here?')
            return -1
        else:

            # TODO: I'm not sure if these calculations are correct.
            # The way I store Insertions and Deletions is kind of confusing for this data.

            if nucleotide == "A":
                return (100.0 * self.aCount) / totalNucleotideCount
            elif nucleotide == "G":
                return (100.0 * self.gCount) / totalNucleotideCount
            elif nucleotide == "C":
                return (100.0 * self.cCount) / totalNucleotideCount
            elif nucleotide == "T":
                return (100.0 * self.tCount) / totalNucleotideCount
            else:
                raise Exception(
                    "I do not know what this nucleotide is:" + str(nucleotide)
                )

    def getMostCommonBase(self):
        maxPolymorphismTypeCount = max(
            (self.matchCount, self.mismatchCount, self.inCount, self.delCount)
        )
        maxNucleotideTypeCount = max(
            (self.aCount, self.gCount, self.cCount, self.tCount)
        )

        if maxPolymorphismTypeCount == self.matchCount:
            return self.referenceBase
        elif maxPolymorphismTypeCount == self.delCount:
            return "-"
        elif (
            maxPolymorphismTypeCount == self.inCount
            or maxPolymorphismTypeCount == self.mismatchCount
        ):

            if maxNucleotideTypeCount == self.aCount:
                return "A"
            elif maxNucleotideTypeCount == self.gCount:
                return "G"
            elif maxNucleotideTypeCount == self.cCount:
                return "C"
            elif maxNucleotideTypeCount == self.tCount:
                return "T"
            else:
                raise Exception(
                    "I do not know what to do with this max value, something is wrong."
                )

        else:
            raise Exception(
                "I do not know what to do with this max value, something is wrong."
            )

    # Actually kind of hard to calculate a total.
    # I think the nucleotide counts + deletion counts should do it.
    def getMappedReadCount(self):
        # return self.aCount + self.cCount + self.gCount + self.tCount + self.delCount
        return self.matchCount + self.mismatchCount + self.inCount + self.delCount

    # def getMatchPercent(self):
    #     #totalDataCount = self.matchCount + self.mismatchCount + self.inCount + self.delCount
    #     return ((100.0 * self.matchCount)/self.getMappedReadCount())

    # def getMismatchPercent(self):
    #     #totalDataCount = self.matchCount + self.mismatchCount + self.inCount + self.delCount
    #     return ((100.0 * self.mismatchCount)/self.getMappedReadCount())

    # def getInsertionPercent(self):
    #     #totalDataCount = self.matchCount + self.mismatchCount + self.inCount + self.delCount
    #     return ((100.0 * self.inCount)/self.getMappedReadCount())

    # def getDeletionPercent(self):
    #     #totalDataCount = self.matchCount + self.mismatchCount + self.inCount + self.delCount
    #     return ((100.0 * self.delCount)/self.getMappedReadCount())
    def getMatchPercent(self):
        try:
            mapped_count = self.getMappedReadCount()
            if mapped_count > 0:
                return (100.0 * self.matchCount) / mapped_count
            else:
                return 0.0  # Or any other appropriate value when count is zero
        except ZeroDivisionError:
            return 0.0

    def getMismatchPercent(self):
        try:
            mapped_count = self.getMappedReadCount()
            if mapped_count > 0:
                return (100.0 * self.mismatchCount) / mapped_count
            else:
                return 0.0
        except ZeroDivisionError:
            return 0.0

    def getInsertionPercent(self):
        try:
            mapped_count = self.getMappedReadCount()
            if mapped_count > 0:
                return (100.0 * self.inCount) / mapped_count
            else:
                return 0.0
        except ZeroDivisionError:
            return 0.0

    def getDeletionPercent(self):
        try:
            mapped_count = self.getMappedReadCount()
            if mapped_count > 0:
                return (100.0 * self.delCount) / mapped_count
            else:
                return 0.0
        except ZeroDivisionError:
            return 0.0

    def deleteBase(self):
        self.delCount += 1

    def addNewBase(self, newBases):
        if newBases == "A":
            self.aCount += 1
        elif newBases == "G":
            self.gCount += 1
        elif newBases == "C":
            self.cCount += 1
        elif newBases == "T":
            self.tCount += 1
        else:
            raise Exception("Unknown nucleotide base:" + newBases)

    def matchBase(self, newBases):
        self.matchCount += 1
        self.addNewBase(newBases)

    def insertBases(self, newBases):
        self.inCount += 1
        self.addNewBase(newBases)

    def mismatchBase(self, newBases):
        self.mismatchCount += 1
        self.addNewBase(newBases)
