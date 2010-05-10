/*==========================================================================
   SeqAn - The Library for Sequence Analysis
   http://www.seqan.de 
  ==========================================================================
   Copyright (C) 2010
  
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.
  
  ==========================================================================
   Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ==========================================================================
   Usage: compare_sam_wit [options] <contig.fasta> <result.sam> <golden.wit>

   Call as "compare_sam_wit --help" for a full list of options.
  
   This program is used to compare the read mapper results in a SAM file
   with the golden standard in a WIT file.  Log messages are printed to
   stderr, the result is written to the first line of stdout.

   The first line of stdout consists of a JSON encoded dictionary that has
   the following entries:

   * total_intervals       Total number of intervals with the given error
                           rate in the WIT file.
   * found_intervals       Number of intervals in the WIT file that were
                           found by the read mapper.
   * superflous_intervals  Number of alignments in the SAM file that do not
                           have their end position in an interval from the
                           WIT file and the alignment at this position has
                           a higher error rate than the specified one.
   * additional_intervals  Same as superflous_intervals, but these alignments
                           have an error rate that is as low as the specified
                           one or lower.  If this happens in the non-weighted
                           case then this is a bug in the benchmark tools.
                           For the weighted case, this happens if RazerS does
                           not find the alignment with low weighted but too
                           high unweighted error rate and the read mapper
                           generating the SAM file finds it.  Read the paper
                           and/or manual for more details.
  ==========================================================================*/

#include <map>

#include <seqan/find.h>                 // Finding infrastructure.
#include <seqan/misc/misc_cmdparser.h>  // Command parser.
#include <seqan/store.h>                // Fragment store et al.

#include "wit_store.h"
#include "find_hamming_simple_ext.h"
#include "find_myers_ukkonen_reads.h"
#include "find_myers_ukkonen_ext.h"
#include "find_approx_dp_quality.h"
#include "find_hamming_simple_quality.h"

using namespace seqan;  // Remove some syntatic noise.


// ======================================================================
// Options
// ======================================================================

// Container for the program options.
struct Options {
    // Maximum number or errors per read length in percent.
    int maxError;

    // Print the missed intervals to stderr for debugging purposes.
    bool showMissedIntervals;

    // Print the hit intervals to stderr for debugging purposes.
    bool showHitIntervals;

    // Print each end position that we try to match agains the interval.
    bool showTryHitIntervals;

    // If true, N matches as a wildcard.  Otherwise it matches none.
    bool matchN;

    // If true, use weighted distances instead of unit ones.
    bool weightedDistances;

    // Distance function to use, also see validDistanceFunction.
    String<char> distanceFunction;

    // Name of reference sequence file.
    String<char> seqFileName;

    // Name of WIT file.
    String<char> witFileName;

    // Name of SAM file.
    String<char> samFileName;

    // Return true iff distanceFunction is a valid distance function.
    // Valid distances are one of {"hamming", "edit"}.
    bool validDistanceFunction() const
    {
        if (distanceFunction == "hamming") return true;
        if (distanceFunction == "edit") return true;
        return false;
    }
};


// ======================================================================
// Comparison Result
// ======================================================================


// Counters for the comparison result.
struct ComparisonResult {
    // Total number of intervals in golden standard.
    size_t totalIntervalCount;

    // Number of intervals in golden standard that were found by the
    // read mapper.
    size_t foundIntervalCount;

    // Number of intervals with a too high distance that were found by
    // the read mappers, i.e. "junk output".
    size_t superflousIntervalCount;

    // Number of intervals with a good score that the read mapper
    // found which were not in our golden standard.
    size_t additionalIntervalCount;
    
    ComparisonResult()
            : totalIntervalCount(0), foundIntervalCount(0),
              superflousIntervalCount(0), additionalIntervalCount(0) {}
};


// Output-to-stream operator for ComparisonResult.
template <typename TStream>
TStream & operator<<(TStream & stream, ComparisonResult const & result) {
    stream << "{\"total_intervals\": " << result.totalIntervalCount
           << ", \"found_intervals\": " << result.foundIntervalCount
           << ", \"superflous_intervals\": " << result.superflousIntervalCount
           << ", \"additional_intervals\": " << result.additionalIntervalCount
           << "}";
    return stream;
}


// Resetting all counters to 0 for ComparisonResult.
void clear(ComparisonResult & result) {
    result.totalIntervalCount = 0;
    result.foundIntervalCount = 0;
    result.superflousIntervalCount = 0;
    result.additionalIntervalCount = 0;
}


// ======================================================================
// FlaggedInterval
// ======================================================================

struct FlaggedInterval {
    size_t firstPos;
    size_t lastPos;
    bool flag;

    FlaggedInterval() {}

    FlaggedInterval(size_t _firstPos, size_t _lastPos, bool _flag = false)
            : firstPos(_firstPos), lastPos(_lastPos), flag(_flag) {}
};


template <typename TStream>
TStream & operator<<(TStream & stream, FlaggedInterval const & interval) {
    return stream << "FlaggedInterval(" << interval.firstPos << ", " << interval.lastPos << ", " << interval.flag << ")";
}


// Copy-and-paste from reweight_wit.h
//
// Compute quality-based alignment score.  The read has to be given
// since we do not have qualities in the alignment object.
template <typename TAlign>
int computeQualityAlignmentScore(TAlign const & align, Score<int, ScoreMatrix<Dna5> > const & scoreMatrix, String<Dna5Q> const & read) {
    // TODO(holtgrew): Maybe convert to iterators for performance?
    typedef typename Row<TAlign const>::Type TRow;
    typedef typename Value<TRow>::Type TAlignChar;
    typedef typename Size<TRow>::Type TSize;

    TRow & rowContig = row(align, 0);
    TRow & rowRead = row(align, 1);

    int result = 0;
    for (TSize i = 0; i < length(rowContig); ++i) {
        TAlignChar contigChar = rowContig[i];
        TAlignChar readChar = rowRead[i];
        if (isGap(rowContig, i)) {
            result -= getQualityValue(read[toSourcePosition(rowRead, i)]);
        } else if (isGap(rowRead, i)) {
            if (toSourcePosition(rowRead, i) == 0) {
                result -= getQualityValue(read[0]);
            } else if (toSourcePosition(rowRead, i) == length(read)) {
                result -= getQualityValue(read[length(read) - 1]);
            } else {
                int x = 0;
                x += getQualityValue(read[toSourcePosition(rowRead, i) - 1]);
                x += getQualityValue(read[toSourcePosition(rowRead, i)]);
                result -= ceil(1.0 * x / 2);
            }
        } else {
            result += score(scoreMatrix, readChar, contigChar) * getQualityValue(read[toSourcePosition(rowRead, i)]);
        }
    }

    return result;
}


// Returns the best score for the alignment of the aligned read from
// the given fragment store.  The maximum error is given, to be able
// to limit the interval in the contig we are looking for.
template <typename TFragmentStore, typename TContigSeq2, typename TAlignedRead, typename TScore, typename TPatternSpec>
int bestScoreForAligned(TFragmentStore & fragments,
                        TContigSeq2 & contig2,
                        bool const & isForward,
                        TAlignedRead const & alignedRead,
                        int /*maxError*/,
                        Options const & options,
                        TScore const & scoringScheme,
                        TPatternSpec const &) {
    typedef size_t TContigId;  // TODO(holtgrew): Better type.
    typedef size_t TAlignedReadPos;  // TODO(holtgrew): Better type.
    typedef typename TFragmentStore::TContigStore      TContigStore;
    typedef typename TFragmentStore::TContigSeq        TContigSeq;
    typedef typename Value<TContigStore>::Type         TContig;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
    typedef typename TFragmentStore::TReadSeqStore     TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq        TReadSeq;
    typedef typename Size<TReadSeq>::Type TSize;
    TContigStore & contigs = fragments.contigStore;
    TContigId contigId = alignedRead.contigId;
    TContigSeq & contig = contigs[contigId].seq;
    TContigGaps contigGaps(contigs[contigId].seq, contigs[contigId].gaps);
    TAlignedReadPos endPos = positionGapToSeq(contigGaps, alignedRead.endPos);
    TAlignedReadPos beginPos = positionGapToSeq(contigGaps, alignedRead.beginPos);
    TReadSeqStore & readSeqs = fragments.readSeqStore;
    TReadSeq read = readSeqs[alignedRead.readId];
            
    // Maybe compute reverse of read, then beginPos and endPos
    // must be exchanged, too.
//     std::cout << "endPos = " << endPos << ", beginPos = " << beginPos << std::endl;
    if (!isForward) {
        SEQAN_ASSERT_GT(beginPos, endPos);
        beginPos = length(contig) - beginPos;
        endPos = length(contig) - endPos;
    }

    Finder<TContigSeq2> finder(contig2);
    Pattern<TReadSeq, TPatternSpec> pattern(read, -length(read) * 1000);
    _patternMatchNOfPattern(pattern, options.matchN);
    _patternMatchNOfFinder(pattern, options.matchN);
    bool ret = setEndPosition(finder, pattern, endPos);
//     std::cerr << __FILE__ << ":" << __LINE__ << " -- endPos = " << endPos << ", endPosition(finder) == " << endPosition(finder) << std::endl;
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(endPos, endPosition(finder));
    
    // No explicit alignment is required if distances are not to be weighted.
    if (!options.weightedDistances)
        return getScore(pattern);

    // Otherwise, we need to build an alignment and compute the score from it.
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT_TRUE(ret);

    // Prepare alignment datastructures.
    Align<String<Dna5>, ArrayGaps> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), infix(finder));
    assignSource(row(align, 1), read);

    // Perform banded Needleman-Wunsch alignment.
    // TODO(holtgrew): Use wrapper for trace pumping, once it is in place.
    StringSet<String<Dna5> > stringSet;
    appendValue(stringSet, infix(finder));
    appendValue(stringSet, read);
    _Align_Traceback<TSize> trace;
    int alignmentScore = globalAlignment(trace, stringSet, scoringScheme, getScore(pattern), -getScore(pattern), BandedNeedlemanWunsch());
    _pump_trace_2_Align(align, trace);
    SEQAN_ASSERT_EQ(alignmentScore, getScore(pattern));

    // Compute quality-based score of alignment.  We pass the
    // score matrix to allow for N-is-wildcard mode.
    int qualityValue = computeQualityAlignmentScore(align, scoringScheme, read);
    return qualityValue;
}

// Return maximum error count for maximum error rate
inline int maxErrorRateToMaxErrors(int maxErrorRate, size_t len) {
    return floor(maxErrorRate / 100.0 * len);
}


template <typename TFragmentStore, typename TAlignedReadIter, typename TWitRecordIter, typename TContigSeq, typename TContigGaps, typename TPatternSpec>
void
compareAlignedReadsToReferenceOnContigForOneRead(Options const & options,
                                                 TFragmentStore & fragments,
                                                 size_t const & contigId,
                                                 TContigSeq /*const*/ & contig,
                                                 TContigGaps /*const*/ & contigGaps,
                                                 bool const & isForward,
                                                 size_t const & readId,
                                                 TAlignedReadIter const & alignedReadsBegin,
                                                 TAlignedReadIter const & alignedReadsEnd,
                                                 TWitRecordIter const & witRecordsBegin,
                                                 TWitRecordIter const & witRecordsEnd,
                                                 ComparisonResult & result,
                                                 TPatternSpec const &) {
    // TODO(holtgrew): Enable n-matches-wildcard mode and use qualities for weights.
    typedef size_t TPos;
    typedef std::map<size_t, FlaggedInterval> TMap;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename TMap::iterator TMapIterator;

    // TODO(holtgrew): Only necessary for weighted variant.
    // Build scoring matrix that allows N to match with all.
    int gapExtensionScore = -1;
    int gapOpenScore = -1;
    if (TYPECMP<TPatternSpec, HammingSimple>::VALUE) {
        // No gaps for hamming distance.
        gapOpenScore = -length(fragments.readSeqStore[readId]);
        gapExtensionScore = -length(fragments.readSeqStore[readId]);
    }
    Score<int, ScoreMatrix<Dna5> > matrixScore(gapExtensionScore, gapOpenScore);
    for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x) {
        for (int y = 0; y < ValueSize<Dna5>::VALUE; ++y) {
            setScore(matrixScore, Dna5(x), Dna5(y), -1);
        }
        setScore(matrixScore, Dna5(x), Dna5('N'), 0);
        setScore(matrixScore, Dna5('N'), Dna5(x), 0);
        setScore(matrixScore, Dna5(x), Dna5(x), 0);
    }

//     std::cerr << "compare aligned reads to reference on contig (" << contigId << ") for one read (" << readId << ")" << std::endl;
//     std::cerr << ".--- readId == " << readId << std::endl;
//     for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
//         std::cerr << "| " << *it << std::endl;
//     }
//     std::cerr << "`---" << std::endl;

    // Build map of intervals to find.
    //
    // Map is "end position -> (found flag, interval first, interval last)".
    TMap intervalMap;
    for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
//         std::cerr << *it << std::endl;
        if (it->isForward != isForward) continue;  // Skip aligned reads on other strand.
        if (static_cast<int>(it->distance) != options.maxError) continue;  // Skip intervals with wrong distance.

//         std::cerr << __FILE__ << ":" << __LINE__ << " -- add(" << it->firstPos << ", " << it->lastPos << ", is forward = " << it->isForward << ")" << "  " << *it << std::endl;
        intervalMap[it->lastPos] = FlaggedInterval(it->firstPos, it->lastPos);
        result.totalIntervalCount += 1;
    }

    // Now, try to hit all entries in intervals with the aligned reads on this strand.
//     std::cerr << __FILE__ << ":" << __LINE__ << " -- aligned read # = " << alignedReadsEnd - alignedReadsBegin << std::endl;
    for (TAlignedReadIter it = alignedReadsBegin; it != alignedReadsEnd; ++it) {
        // Skip aligned reads on other strand.
        if (isForward && it->beginPos > it->endPos) {
//             std::cerr << __FILE__ << ":" << __LINE__ << " -- skipping, isForward = " << isForward << ", it->beginPos = " << it->beginPos << ", it->endPos = " << it->endPos << std::endl;
            continue;
        }
        if (!isForward && it->beginPos <= it->endPos) {
//             std::cerr << __FILE__ << ":" << __LINE__ << " -- skipping, isForward = " << isForward << ", it->beginPos = " << it->beginPos << ", it->endPos = " << it->endPos << std::endl;
            continue;
        }

        // Convert from gap space to sequence space and maybe into reverse strand position.
        TPos endPos = positionGapToSeq(contigGaps, it->endPos);
        TPos beginPos = positionGapToSeq(contigGaps, it->beginPos);
        if (!isForward) {
//             std::cerr << "beginPos (before conversion) "  << beginPos << std::endl;
//             std::cerr << "endPos (before conversion) "  << endPos << std::endl;
            endPos = length(contig) - endPos;
            beginPos = length(contig) - beginPos;
//             std::cerr << "beginPos (after conversion) "  << beginPos << std::endl;
//             std::cerr << "endPos (after conversion) "  << endPos << std::endl;
        }

        // Skip if aligning too far to the left.
        if (endPos < length(fragments.readSeqStore[it->readId])) 
            continue;

        // Compute last position of alignment and search for interval this
        // position falls into.
        TPos lastPos = endPos - 1;
        SEQAN_ASSERT_LEQ(beginPos, lastPos);
        if (options.showTryHitIntervals)
            std::cerr << "Searching for read " << fragments.readNameStore[it->readId]
                      << " on contig " << fragments.contigNameStore[contigId]
                      << " forward? " << isForward
                      << " last pos " << lastPos
                      << " (beginPos == " << beginPos
                      << ", endPos == " << endPos << ")" << std::endl;
        TMapIterator iter = intervalMap.lower_bound(lastPos);

        // Skip reads that aligned with a too bad score.
        int bestScore = bestScoreForAligned(fragments, contig, isForward, *it, maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId])), options, matrixScore, TPatternSpec());
        if (bestScore < -maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId]))) {
            std::cerr << "INFO: Skipping read " << fragments.readNameStore[it->readId] << " (" << fragments.readSeqStore[it->readId] << ")" << std::endl;
            std::cerr << "  contigId = " << contigId << std::endl;
            std::cerr << "  forward strand? " << isForward << std::endl;
            std::cerr << "  begin pos = " << beginPos << ", endPos == " << endPos << std::endl;
            std::cerr << "  infix(contig, beginPos, endPos) == " << infix(contig, beginPos, endPos) << std::endl;
            std::cerr << "  score is " << bestScore << std::endl;
            std::cerr << "  max error rate is " << options.maxError << std::endl;
            std::cerr << "  read length is " << length(fragments.readSeqStore[it->readId]) << std::endl;
            std::cerr << "  read is " << fragments.readSeqStore[it->readId] << std::endl;
            std::cerr << "  read qualities are ";
            for (unsigned i = 0; i < length(fragments.readSeqStore[it->readId]); ++i) {
                std::cerr << getQualityValue(fragments.readSeqStore[it->readId][i]) << ", ";
            }
            std::cerr << std::endl;
            std::cerr << "  max errors is " <<  maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId])) << std::endl;
            result.superflousIntervalCount += 1;
            continue;
        }

        // Handle alignment out of target intervals.
        if (iter == intervalMap.end() || (iter->second.firstPos > lastPos) || (iter->second.lastPos < lastPos)) {
            if (options.weightedDistances)
                std::cerr << "WARNING: ";
            else
                std::cerr << "PANIC: ";
            std::cerr << "A read in the SAM file aligns out of all target intervals for this read in the WIT file." << std::endl;
            if (iter == intervalMap.end())
                std::cerr << "not found in map" << std::endl;
            else
                std::cerr << "found was " << iter->second << std::endl;
            std::cerr << "bestScore = " << bestScore << std::endl;
            std::cerr << "read name = " << fragments.readNameStore[it->readId] << std::endl;
            std::cerr << "read is = " << fragments.readSeqStore[it->readId] << std::endl;
            Dna5String rcRead(fragments.readSeqStore[it->readId]);
            reverseComplementInPlace(rcRead);
            std::cerr << "          " << rcRead << std::endl;
            std::cerr << "on forward strand? " << isForward << std::endl;
            std::cerr << "original begin pos = " << it->beginPos << std::endl;
            std::cerr << "original end pos = " << it->endPos << std::endl;
            std::cerr << "begin pos = " << beginPos << std::endl;
            std::cerr << "end pos = " << endPos << std::endl;
            std::cerr << "last pos = " << lastPos << std::endl;
            std::cerr << "contigId = " << it->contigId << std::endl;
            std::cerr << "max error rate is " << options.maxError << std::endl;
            std::cerr << "read length is " << length(fragments.readSeqStore[it->readId]) << std::endl;
            std::cerr << "max errors is " <<  maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId])) << std::endl;
            if (options.weightedDistances) {
                result.additionalIntervalCount += 1;
                continue;
            } else {
                exit(1);
            }
        }

        SEQAN_ASSERT_LEQ(iter->second.firstPos, lastPos);
        SEQAN_ASSERT_GEQ(iter->second.lastPos, lastPos);

        // Count interval if hit for the this time.
        if (!iter->second.flag) {
            iter->second.flag = true;
            result.foundIntervalCount += 1;
        }
    }

    // If configured so, show hit and missed intervals.
    if (options.showHitIntervals || options.showMissedIntervals) {
        for (TMapIterator it = intervalMap.begin(); it != intervalMap.end(); ++it) {
            if (!it->second.flag && options.showMissedIntervals) {
                std::cerr << "Missed interval for read " << fragments.readNameStore[readId]
                          << " on contig " << fragments.contigNameStore[contigId]
                          << " forward strand? " << isForward << " -- ["
                          << it->second.firstPos << ", " << it->second.lastPos << "]" << std::endl;
            } else if (it->second.flag && options.showHitIntervals) {
                std::cerr << "Hit interval for read " << fragments.readNameStore[readId]
                          << " on contig " << fragments.contigNameStore[contigId]
                          << " forward strand? " << isForward << " -- ["
                          << it->second.firstPos << ", " << it->second.lastPos << "]" << std::endl;
            }
        }
    }
}


// Compare the aligned reads in [alignedReadsBegin, alignedReadsEnd)
// to the intervals in [witRecordsBegin, witRecordsEnd) on the given
// contig on the forward strand iff isForward.
//
// foundIntervalCount is incremented by the number of hit intervals,
// the number of relevant (on the selected strand) intervals is added
// to relevantIntervalCount.
template <typename TFragmentStore, typename TAlignedReadsIter, typename TWitRecordsIter, typename TContigSeq, typename TPatternSpec>
void
compareAlignedReadsToReferenceOnContig(Options const & options,
                                       TFragmentStore & fragments,
                                       size_t const & contigId,
                                       TContigSeq /*const*/ & contig,
                                       bool const & isForward,
                                       TAlignedReadsIter const & alignedReadsBegin,
                                       TAlignedReadsIter const & alignedReadsEnd,
                                       TWitRecordsIter const & witRecordsBegin,
                                       TWitRecordsIter const & witRecordsEnd,
                                       ComparisonResult & result,
                                       TPatternSpec const &) {
    typedef size_t TPos;
    typedef std::map<size_t, FlaggedInterval> TMap;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContigStoreElement::TGapAnchors> > TContigGaps;
    typedef typename TMap::iterator TMapIterator;

//     std::cerr << "contig id = " << contigId << " forward = " << isForward << std::endl;
//     std::cerr << "# of aligned reads " << (alignedReadsEnd - alignedReadsBegin) << ", # of wit records " << (witRecordsEnd - witRecordsBegin) << std::endl;

    // Build contig gaps datastructure for gap space to sequence space conversion.
    TContigGaps contigGaps(fragments.contigStore[contigId].seq, fragments.contigStore[contigId].gaps);

    // The aligned reads and wit records are filtered to this contig and
    // sorted by read id.

    TAlignedReadsIter alignedReadsBeginForRead = alignedReadsBegin;
    TAlignedReadsIter alignedReadsEndForRead = alignedReadsBeginForRead;
    TWitRecordsIter witRecordsBeginForRead = witRecordsBegin;
    TWitRecordsIter witRecordsEndForRead = witRecordsBegin;

    while (alignedReadsBeginForRead != alignedReadsEnd && witRecordsBeginForRead != witRecordsEnd) {
        // Get current contigId.
        size_t readId = _min(alignedReadsBeginForRead->readId, witRecordsBeginForRead->readId);
//         std::cerr << "alignedReadsBeginForRead->readId == " << alignedReadsBeginForRead->readId << ", witRecordsBeginForRead->readId == " <<  witRecordsBeginForRead->readId << std::endl;
//         std::cerr << ".--- readId == " << readId << std::endl;
        // Get aligned reads iterator for the next contigId.
        while (alignedReadsEndForRead != alignedReadsEnd && alignedReadsEndForRead->readId <= readId) {
//             std::cerr << "| read id = " << alignedReadsEndForRead->readId << ", begin pos = " << positionGapToSeq(contigGaps, alignedReadsEndForRead->beginPos) << ", end pos = " << positionGapToSeq(contigGaps, alignedReadsEndForRead->endPos) << std::endl;
            ++alignedReadsEndForRead;
        }
//         std::cerr << "+---" << std::endl;
        // Get wit records iterator for the next contigId.
        while (witRecordsEndForRead != witRecordsEnd && witRecordsEndForRead->readId <= readId) {
//             if (witRecordsEndForRead->distance == options.maxError)
//                 std::cerr << "| " << *witRecordsEndForRead << std::endl;
            ++witRecordsEndForRead;
        }
//         std::cerr << "`---" << std::endl;

        // Actually compare the aligned reads for this contig on forward and backwards strand.
        compareAlignedReadsToReferenceOnContigForOneRead(options, fragments, contigId, contig, contigGaps, isForward, readId, alignedReadsBeginForRead, alignedReadsEndForRead, witRecordsBeginForRead, witRecordsEndForRead, result, TPatternSpec());
        // This iteration's end iterators are the next iteration's begin iterators.
        alignedReadsBeginForRead = alignedReadsEndForRead;
        witRecordsBeginForRead = witRecordsEndForRead;
    }
}


// Compare aligned reads in fragment store to the intervals specified
// in witStore, results is used for the counting statistics.
template <typename TFragmentStore, typename TPatternSpec>
void
compareAlignedReadsToReference(ComparisonResult & result,
                               TFragmentStore & fragments,  // non-const so reads can be sorted
                               WitStore & witStore,  // non-const so it is sortable
                               Options const & options,
                               TPatternSpec const &) {
    // Type aliases.
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore>::Type TAlignedReadsIter;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef typename TContigStoreElement::TContigSeq TContigSeq;
    typedef typename WitStore::TIntervalStore TIntervalStore;
    typedef typename Iterator<TIntervalStore>::Type TIntervalIter;
    
    // Initialization.
    clear(result);

    // Sort aligned reads and wit records by contig id.
    sortAlignedReads(fragments.alignedReadStore, SortEndPos());
    sortAlignedReads(fragments.alignedReadStore, SortReadId());
    sortAlignedReads(fragments.alignedReadStore, SortContigId());
    sortWitRecords(witStore, SortLastPos());
    sortWitRecords(witStore, SortReadId());
    sortWitRecords(witStore, SortContigId());

    TAlignedReadsIter alignedReadsBegin = begin(fragments.alignedReadStore, Standard());
    TAlignedReadsIter alignedReadsEnd = alignedReadsBegin;
//     std::cerr << "# of alignd reads " << length(fragments.alignedReadStore) << std::endl;
    TIntervalIter witRecordsBegin = begin(witStore.intervals, Standard());
    TIntervalIter witRecordsEnd = witRecordsBegin;
//     std::cerr << "# of wit records " << length(witRecords) << std::endl;

//     std::cerr << ".--- wit records" << std::endl;
//     for (TWitRecordsIter it = begin(witRecords, Standard()); it != end(witRecords, Standard()); ++it)
//         std::cerr << "| " << *it << " --- " << it->contigId << std::endl;
//     std::cerr << "`---" << std::endl;

    while (alignedReadsBegin != end(fragments.alignedReadStore, Standard()) &&
           witRecordsBegin != end(witStore.intervals, Standard())) {
        // Get current contigId.
        size_t contigId = _min(alignedReadsBegin->contigId, witRecordsBegin->contigId);
//         std::cerr << "alignedReadsBegin->contigId == " << alignedReadsBegin->contigId << ", witRecordsBegin->contigId == " <<  witRecordsBegin->contigId << std::endl;
//         std::cerr << ".--- contigId == " << contigId << std::endl;
        // Get aligned reads iterator for the next contigId.
        while (alignedReadsEnd != end(fragments.alignedReadStore, Standard()) && alignedReadsEnd->contigId <= contigId) {
//             std::cerr << "| read id = " << alignedReadsEnd->readId << ", begin pos = " << positionGapToSeq(contigGaps, alignedReadsEnd->beginPos) << ", end pos = " << positionGapToSeq(contigGaps, alignedReadsEnd->endPos) << std::endl;
            ++alignedReadsEnd;
        }
//         std::cerr << "+---" << std::endl;
        // Get wit records iterator for the next contigId.
        while (witRecordsEnd != end(witStore.intervals, Standard()) && witRecordsEnd->contigId <= contigId) {
//             if (witRecordsEnd->distance == options.maxError)
//                 std::cerr << "| " << *witRecordsEnd << std::endl;
            ++witRecordsEnd;
        }
//         std::cerr << "`---" << std::endl;

        // Actually compare the aligned reads for this contig on forward and backwards strand.
        TContigSeq & contig = fragments.contigStore[contigId].seq;
        compareAlignedReadsToReferenceOnContig(options, fragments, contigId, contig, true, alignedReadsBegin, alignedReadsEnd, witRecordsBegin, witRecordsEnd, result, TPatternSpec());
        TContigSeq rcContig(contig);
        reverseComplementInPlace(rcContig);
        compareAlignedReadsToReferenceOnContig(options, fragments, contigId, rcContig, false, alignedReadsBegin, alignedReadsEnd, witRecordsBegin, witRecordsEnd, result, TPatternSpec());
        
        // This iteration's end iterators are the next iteration's begin iterators.
        alignedReadsBegin = alignedReadsEnd;
        witRecordsBegin = witRecordsEnd;
    }
}
