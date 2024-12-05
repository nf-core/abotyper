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


# from Bio import pairwise2
from Bio.Align import substitution_matrices, PairwiseAligner


class PairwiseAlignmentResult:

    def __init__(self, sequenceID, referenceSequence, querySequence):
        self.sequence1Aligned = ""
        self.sequence2Aligned = ""
        self.alignmentScore = 0
        self.indexBegin = 0
        self.indexEnd = 0
        self.sequence1AlignedShort = ""
        self.sequence2AlignedShort = ""
        self.alignedSectionLength = 0

        # Instead of storing a count of these events
        # I just want to store an array of the regions.
        # It's a tuple (begin, end)
        # self.matchTuples=[]
        # self.mismatchTuples=[]
        # self.insertionTuples=[]
        # self.deletionTuples=[]

        # Except for all-tuples. Lets just store these tuples:
        # (alignment state, begin, end)
        self.allTuples = []

        # self.matchCount            = 0
        # self.mismatchCount         = 0
        # self.insertionCount        = 0
        # self.deletionCount         = 0
        # self.indelCount            = 0

        # SequenceIDs stores the ids of the alleles in this alignment.
        self.sequenceID = sequenceID
        self.sequence1 = referenceSequence
        self.sequence2 = querySequence

        self.alignPairwise()

    # TODO: Pairwise alignment takes a while.
    # I can try some of the other alignment algorithms.
    # But then, why not just do a batch processing?

    def alignPairwise(self):
        """
        Alignment Parameters.
        A match score is 0.  I chose this so a perfect match is always a score of 0.
        I chose to punish gaps more than mismatches.
        Likely, there is fine-tuning to be done here.
        I feel like match_score should be zero but I don't get any alignments when that happens.
        Maybe the gap penalty should be the same as a mismatch penalty, because errors with repeats are somewhat common.

        I want a local alignment, so it allows sequences of different lengths.
        The method to call depends on how i want to treat mismatches and indels.
        I found some hints in the man pages for biopython's pairwise2 package
        The match parameters are:

        CODE  DESCRIPTION
        x     No parameters. Identical characters have score of 1, otherwise 0.
        m     A match score is the score of identical chars, otherwise mismatch
             score.
        d     A dictionary returns the score of any pair of characters.
        c     A callback function returns scores.
        The gap penalty parameters are:

        CODE  DESCRIPTION
        x     No gap penalties.
        s     Same open and extend gap penalties for both sequences.
        d     The sequences have different open and extend gap penalties.
        c     A callback function returns the gap penalties.
        """

        match_score = 1
        mismatch_score = -1
        gap_open = -10
        gap_extend = -1

        ## ----------------------- OLD ALIGNER -------------------------------------------------
        # alignments = pairwise2.align.globalms(
        #     self.sequence1,
        #     self.sequence2,
        #     match_score,
        #     mismatch_score,
        #     gap_open,
        #     gap_extend)

        ## ----------------------- START CHANGES TO ALIGNER -------------------------------------------------
        # NEW:  TODO We have to make changes to replace Bio.pairwise2 above with PairwiseAligner!!
        # It is deprecated and likely to be removed from future releases of Biopython.
        aligner = PairwiseAligner(
            mode="global",
            # mode = 'local',
            substitution_matrix=substitution_matrices.load("BLOSUM62"),
            open_gap_score=gap_open,
            extend_gap_score=gap_extend,
            match_score=match_score,
            mismatch_score=mismatch_score,
        )

        # Let's try a global alignment with the new aligner :)
        alignments = aligner.align(self.sequence1, self.sequence2)

        ## ----------------------- END CHANGES TO ALIGNER -------------------------------------------------

        topAlignment = alignments[0]

        (
            self.sequence1Aligned,
            self.sequence2Aligned,
            self.alignmentScore,
            self.indexBegin,
            self.indexEnd,
        ) = topAlignment

        self.sequence1AlignedShort = self.sequence1Aligned[
            self.indexBegin : self.indexEnd
        ]
        self.sequence2AlignedShort = self.sequence2Aligned[
            self.indexBegin : self.indexEnd
        ]

        # Count the number of mismatches etc.
        # An "insertion" is an extra base in sequence 2, not accounted in sequence 1.
        # A "deletion" is a base in sequence 1, missing in sequence 2
        # An indel is both of those.
        self.alignedSectionLength = len(self.sequence1AlignedShort)

        # Possible States = 'match' 'mismatch' 'insertion' 'deletion
        # Start at 0, because python enumerate method starts at 0.
        currentState = "match"
        beginIndex = 0

        # Iterate the bases of the alignment, report regions that are match, mismatch, del, etc.
        for baseIndex, seq1Base in enumerate(self.sequence1AlignedShort):
            seq2Base = self.sequence2AlignedShort[baseIndex]
            if seq1Base == seq2Base:
                if currentState == "match":
                    # We're STILL In a match state.  Do nothing.
                    pass
                else:
                    self.storeAlignmentRegion(currentState, beginIndex, baseIndex - 1)
                    beginIndex = baseIndex
                    currentState = "match"

            elif seq1Base == "-":
                if currentState == "insertion":
                    pass
                else:
                    self.storeAlignmentRegion(currentState, beginIndex, baseIndex - 1)
                    beginIndex = baseIndex
                    currentState = "insertion"

            elif seq2Base == "-":
                if currentState == "deletion":
                    pass
                else:
                    self.storeAlignmentRegion(currentState, beginIndex, baseIndex - 1)
                    beginIndex = baseIndex
                    currentState = "deletion"

            else:
                if currentState == "mismatch":
                    pass
                else:
                    self.storeAlignmentRegion(currentState, beginIndex, baseIndex - 1)
                    beginIndex = baseIndex
                    currentState = "mismatch"

        # Guess I need to store the "last" aligned region.
        self.storeAlignmentRegion(currentState, beginIndex, baseIndex - 1)

    def storeAlignmentRegion(self, state, beginIndex, endIndex):
        self.allTuples.append((state, beginIndex, endIndex))
