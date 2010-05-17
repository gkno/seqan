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
   This file contains the flesh of the program compare_sam_wit.
  ==========================================================================*/

#include <algorithm>
#include <map>

#include <seqan/find.h>                 // Finding infrastructure.
#include <seqan/misc/misc_cmdparser.h>  // Command parser.
#include <seqan/store.h>                // Fragment store et al.
#include <seqan/align.h>
#include <seqan/graph_align.h>

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

    // Print superflous intervals (intervals found in SAM file but have too bad score).
    bool showSuperflousIntervals;

    // Print additional intervals (intervals found in SAM with good score that are not in WIT file).
    bool showAdditionalIntervals;

    // Print the hit intervals to stderr for debugging purposes.
    bool showHitIntervals;

    // Print each end position that we try to match agains the interval.
    bool showTryHitIntervals;

    // If true, N matches as a wildcard.  Otherwise it matches none.
    bool matchN;

    // If true, use weighted distances instead of unit ones.
    bool weightedDistances;

    // The benchmark category, one of {"all", "any-best", "all-best"}.
    CharString benchmarkCategory;

    // Distance function to use, also see validDistanceFunction.
    CharString distanceFunction;

    // Name of reference sequence file.
    CharString seqFileName;

    // Name of WIT file.
    CharString witFileName;

    // Name of SAM file.
    CharString samFileName;

    // Return true iff distanceFunction is a valid distance function.
    // Valid distances are one of {"hamming", "edit"}.
    bool validDistanceFunction() const
    {
        if (distanceFunction == "hamming") return true;
        if (distanceFunction == "edit") return true;
        return false;
    }

    // Return true iff benchmarkCategory is a valid benchmark
    // category, i.e. one of {"all", "any-best", "all-best"}.
    bool validBenchmarkCategory() const
    {
        if (benchmarkCategory == "all") return true;
        if (benchmarkCategory == "any-best") return true;
        if (benchmarkCategory == "all-best") return true;
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
    size_t intervalId;
    size_t firstPos;
    size_t lastPos;
    bool flag;

    FlaggedInterval() {}

    FlaggedInterval(size_t _intervalId, size_t _firstPos, size_t _lastPos, bool _flag = false)
            : intervalId(_intervalId), firstPos(_firstPos), lastPos(_lastPos), flag(_flag) {}
};


template <typename TStream>
TStream & operator<<(TStream & stream, FlaggedInterval const & interval) {
    return stream << "FlaggedInterval(" << interval.intervalId << ", " << interval.firstPos << ", " << interval.lastPos << ", " << interval.flag << ")";
}


// Copy-and-paste from reweight_wit.h
//
// Compute quality-based alignment score.  The read has to be given
// since we do not have qualities in the alignment object.
template <typename TAlign>
int computeQualityAlignmentScore(TAlign const & align,
                                 Score<int, ScoreMatrix<Dna5> > const & scoreMatrix,
                                 String<Dna5Q> const & read) {
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

    // If we are aligning on the reverse strand then we have to compute the
    // begin and end position on this strand.
    if (!isForward) {
        SEQAN_ASSERT_GT(beginPos, endPos);
        beginPos = length(contig) - beginPos;
        endPos = length(contig) - endPos;
    }

    // Initialize finder and pattern, configure to match N with none or all,
    // depending on configuration.
    Finder<TContigSeq2> finder(contig2);
    Pattern<TReadSeq, TPatternSpec> pattern(read, -length(read) * 1000);
    _patternMatchNOfPattern(pattern, options.matchN);
    _patternMatchNOfFinder(pattern, options.matchN);
    bool ret = setEndPosition(finder, pattern, endPos);
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
    int alignmentScore = globalAlignment(align, stringSet, scoringScheme, getScore(pattern), -getScore(pattern), BandedNeedlemanWunsch());
    (void)alignmentScore; // Supress warning in non-debug mode.
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
                                                 String<size_t> & result,
                                                 TPatternSpec const &) {
    typedef size_t TPos;
    typedef std::map<size_t, FlaggedInterval> TMap;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename TMap::iterator TMapIterator;

    // Build scoring matrix that allows N to match with all.
    int gapExtensionScore = -1;
    int gapOpenScore = -1;
    if (TYPECMP<TPatternSpec, HammingSimple>::VALUE) {
        // No gaps for hamming distance.
        gapOpenScore = -length(fragments.readSeqStore[readId]);
        gapExtensionScore = -length(fragments.readSeqStore[readId]);
    }
    // Build scoring matrix.
    Score<int, ScoreMatrix<Dna5> > matrixScore(gapExtensionScore, gapOpenScore);
    for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x) {
        for (int y = 0; y < ValueSize<Dna5>::VALUE; ++y) {
            setScore(matrixScore, Dna5(x), Dna5(y), -1);
        }
        setScore(matrixScore, Dna5(x), Dna5(x), 0);
    }
    // Update score matrix if N works as a wildcard character.
    if (options.matchN) {
        for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x) {
            setScore(matrixScore, Dna5(x), Dna5('N'), 0);
            setScore(matrixScore, Dna5('N'), Dna5(x), 0);
        }
    }

    // Build map of intervals to find.
    //
    // Map is "end position -> (found flag, interval first, interval last)".
    TMap intervalMap;
    for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
        if (it->isForward != isForward) continue;  // Skip aligned reads on other strand.
        if (static_cast<int>(it->distance) != options.maxError) continue;  // Skip intervals with wrong distance.

        intervalMap[it->lastPos] = FlaggedInterval(it->id, it->firstPos, it->lastPos);
    }

    // Now, try to hit all entries in intervals with the aligned reads on this strand.
    for (TAlignedReadIter it = alignedReadsBegin; it != alignedReadsEnd; ++it) {
        // Skip aligned reads on other strand.
        if (isForward && it->beginPos > it->endPos)
            continue;
        if (!isForward && it->beginPos <= it->endPos)
            continue;

        // Convert from gap space to sequence space and maybe into reverse strand position.
        TPos endPos = positionGapToSeq(contigGaps, it->endPos);
        TPos beginPos = positionGapToSeq(contigGaps, it->beginPos);
        if (!isForward) {
            endPos = length(contig) - endPos;
            beginPos = length(contig) - beginPos;
        }

        // Skip if aligning too far to the left.
        if (endPos < length(fragments.readSeqStore[it->readId])) 
            continue;

        // Compute last position of alignment and search for interval this
        // position falls into.
        TPos lastPos = endPos - 1;
        SEQAN_ASSERT_LEQ(beginPos, lastPos);
        if (options.showTryHitIntervals) {
            std::cerr << "log> {\"type\": \"log.try_hit"
                      << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                      << "\", \"read_id\": \"" << fragments.readNameStore[it->readId]
                      << ", \"strand\": \"" << (isForward ? "forward" : "reverse")
                      << "\", \"alignment_first\": " << beginPos
                      << ", \"alignment_end\": " << endPos
                      << ", \"read_seq\": \"" << fragments.readSeqStore[it->readId]
                      << "\", \"contig_infix_seq\": \"" << infix(contig, beginPos, endPos);
            std::cerr << "\", \"qualities\": [";
            for (unsigned i = 0; i < length(fragments.readSeqStore[it->readId]); ++i) {
                if (i > 0)
                    std::cerr << ", ";
                std::cerr << getQualityValue(fragments.readSeqStore[it->readId][i]);
            }
            std::cerr << "]}" << std::endl;
        }
        TMapIterator iter = intervalMap.lower_bound(lastPos);

        // Skip reads that aligned with a too bad score.
        int bestScore = bestScoreForAligned(fragments, contig, isForward, *it, maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId])), options, matrixScore, TPatternSpec());
        if (bestScore < -maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId]))) {
            if (options.showSuperflousIntervals) {
                std::cerr << "log> {\"type\": \"log.superflous_hit"
                          << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                          << "\", \"read_id\": \"" << fragments.readNameStore[it->readId]
                          << "\", \"distance\": " << -bestScore
                          << ", \"strand\": \"" << (isForward ? "forward" : "reverse")
                          << "\", \"alignment_begin\": " << beginPos
                          << ", \"alignment_end\": " << endPos
                          << ", \"read_seq\": \"" << fragments.readSeqStore[it->readId]
                          << "\", \"contig_infix_seq\": \"" << infix(contig, beginPos, endPos);
                std::cerr << "\", \"qualities\": [";
                for (unsigned i = 0; i < length(fragments.readSeqStore[it->readId]); ++i) {
                    if (i > 0)
                        std::cerr << ", ";
                    std::cerr << getQualityValue(fragments.readSeqStore[it->readId][i]);
                }
                std::cerr << "]}" << std::endl;
            }
            appendValue(result, IntervalOfReadOnContig::superflousId());
            continue;
        }

        // Handle alignment out of target intervals.
        if (iter == intervalMap.end() || (iter->second.firstPos > lastPos) || (iter->second.lastPos < lastPos)) {
            if (options.showAdditionalIntervals) {
                if (options.weightedDistances) {
                    std::cerr << "log> {\"type\": \"log.additional_hit"
                              << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                              << "\", \"read_id\": \"" << fragments.readNameStore[it->readId]
                              << "\", \"distance\": " << -bestScore
                              << ", \"strand\": \"" << (isForward ? "forward" : "reverse")
                              << "\", \"begin_pos\": " << beginPos
                              << ", \"end_pos\": " << endPos
                              << ", \"read_seq\": \"" << fragments.readSeqStore[it->readId]
                              << "\", \"contig_infix_seq\": \"" << infix(contig, beginPos, endPos);
                    std::cerr << "\", \"qualities\": [";
                    for (unsigned i = 0; i < length(fragments.readSeqStore[it->readId]); ++i) {
                        if (i > 0)
                            std::cerr << ", ";
                        std::cerr << getQualityValue(fragments.readSeqStore[it->readId][i]);
                    }
                    std::cerr << "]}" << std::endl;
                } else {
                    std::cerr << "PANIC: A read in the SAM file aligns out of all target intervals for this read in the WIT file." << std::endl;
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
                }
            }
            if (options.weightedDistances) {
                appendValue(result, IntervalOfReadOnContig::additionalId());
                continue;
            } else {
                exit(1);
            }
        }

        SEQAN_ASSERT_LEQ(iter->second.firstPos, lastPos);
        SEQAN_ASSERT_GEQ(iter->second.lastPos, lastPos);

        // Count interval if hit for the first time.
        if (!iter->second.flag) {
            appendValue(result, iter->second.intervalId);
            iter->second.flag = true;
        }
    }

    // If configured so, show hit and missed intervals.
    if (options.showHitIntervals || options.showMissedIntervals) {
        for (TMapIterator it = intervalMap.begin(); it != intervalMap.end(); ++it) {
            if (!it->second.flag && options.showMissedIntervals) {
                std::cout << "log> {\"type\": \"log.missed_interval"
                          << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                          << "\", \"strand\": \"" << (isForward ? "forward" : "reverse")
                          << "\", \"read_id\": \"" << fragments.readNameStore[readId]
                          << "\", \"interval_first\": " << it->second.firstPos
                          << ", \"interval_last\": " << it->second.lastPos << "}" << std::endl;
            } else if (it->second.flag && options.showHitIntervals) {
                std::cout << "log> {\"type\": \"log.hit_interval"
                          << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                          << "\", \"strand\": \"" << (isForward ? "forward" : "reverse")
                          << "\", \"read_id\": \"" << fragments.readNameStore[readId]
                          << "\", \"interval_first\": " << it->second.firstPos
                          << ", \"interval_last\": " << it->second.lastPos << "}" << std::endl;
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
                                       String<size_t> & result,
                                       TPatternSpec const &) {
    typedef size_t TPos;
    typedef std::map<size_t, FlaggedInterval> TMap;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContigStoreElement::TGapAnchors> > TContigGaps;
    typedef typename TMap::iterator TMapIterator;

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
        // Get aligned reads iterator for the next contigId.
        while (alignedReadsEndForRead != alignedReadsEnd && alignedReadsEndForRead->readId <= readId)
            ++alignedReadsEndForRead;
        // Get wit records iterator for the next contigId.
        while (witRecordsEndForRead != witRecordsEnd && witRecordsEndForRead->readId <= readId)
            ++witRecordsEndForRead;

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
compareAlignedReadsToReference(String<size_t> & result,
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
    TIntervalIter witRecordsBegin = begin(witStore.intervals, Standard());
    TIntervalIter witRecordsEnd = witRecordsBegin;

    while (alignedReadsBegin != end(fragments.alignedReadStore, Standard()) &&
           witRecordsBegin != end(witStore.intervals, Standard())) {
        // Get current contigId.
        size_t contigId = _min(alignedReadsBegin->contigId, witRecordsBegin->contigId);
        // Get aligned reads iterator for the next contigId.
        while (alignedReadsEnd != end(fragments.alignedReadStore, Standard()) && alignedReadsEnd->contigId <= contigId)
            ++alignedReadsEnd;
        // Get wit records iterator for the next contigId.
        while (witRecordsEnd != end(witStore.intervals, Standard()) && witRecordsEnd->contigId <= contigId)
            ++witRecordsEnd;

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


// Tags for selecting an evaluateFoundIntervals_compareToIntervals specialization.
struct _CategoryAll;
typedef Tag<_CategoryAll> CategoryAll;
struct _CategoryAnyBest;
typedef Tag<_CategoryAnyBest> CategoryAnyBest;
struct _CategoryAllBest;
typedef Tag<_CategoryAllBest> CategoryAllBest;


// result must contain all ids of the intervals in
// [beginIntervalsForRead, endIntervalsForRead) that have the error
// rate options.maxError.  The intervals are sorted by (read id,
// error rate).
void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               Options const & options,
                                               CategoryAll const &) {
    if (beginIntervalsForRead == endIntervalsForRead)
        return;  // Guard: Skip if interval range empty.

    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // Now, get begin and end iterators to the intervals for the error rate
    // configured in options.
    IntervalOfReadOnContig refInterval;
    refInterval.distance = options.maxError;
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        comparisonResult.totalIntervalCount += 1;
        comparisonResult.foundIntervalCount += std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id);
    }
}


void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               Options const & /*options*/,
                                               CategoryAnyBest const &) {
    if (beginIntervalsForRead == endIntervalsForRead)
        return;  // Guard: Skip if interval range empty.

    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // The read mapper has to find one of the intervals with the smallest
    // error rate in [beginIntervalsForRead, endIntervalsForRead).
    //
    // The intervals in this range are sorted by error rate.  We get the
    // smallest error rate first, then get the subrange with the smallest
    // range.
    IntervalOfReadOnContig refInterval;
    refInterval.distance = value(beginIntervalsForRead).distance;
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    // One interval is to be found.
    comparisonResult.totalIntervalCount += 1;

    // Now, try to find one of the intervals.
    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        if (std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id)) {
            comparisonResult.foundIntervalCount += 1;
            break;
        }
    }
}


void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               Options const & /*options*/,
                                               CategoryAllBest const &) {
    if (beginIntervalsForRead == endIntervalsForRead)
        return;  // Guard: Skip if interval range empty.

    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // The read mapper has to find one of the intervals with the smallest
    // error rate in [beginIntervalsForRead, endIntervalsForRead).
    //
    // The intervals in this range are sorted by error rate.  We get the
    // smallest error rate first, then get the subrange with the smallest
    // range.
    IntervalOfReadOnContig refInterval;
    refInterval.distance = value(beginIntervalsForRead).distance;
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    // One interval is to be found.

    // Now, try to find one of the intervals.
    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        comparisonResult.totalIntervalCount += 1;
        comparisonResult.foundIntervalCount += std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id);
    }
}


// Evaluate the ids in result pointing to intervals in witStore.
// Result is written to comparisonResult.
void evaluateFoundIntervals(ComparisonResult & comparisonResult,
                            WitStore & witStore,
                            String<size_t> & result,
                            Options const & options)
{
    typedef Iterator<String<size_t>, Standard>::Type TIdIterator;
    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // Initialization.
    clear(comparisonResult);
    sortWitRecords(witStore, SortDistance());
    sortWitRecords(witStore, SortReadId());
    std::sort(begin(result, Standard()), end(result, Standard()));

    // Count ids of superflous and additional intervals (at end of result).
    std::pair<TIdIterator, TIdIterator> superflousBounds = std::equal_range(begin(result, Standard()), end(result, Standard()), IntervalOfReadOnContig::superflousId());
    comparisonResult.superflousIntervalCount = superflousBounds.second - superflousBounds.first;
    std::pair<TIdIterator, TIdIterator> additionalBounds = std::equal_range(begin(result, Standard()), end(result, Standard()), IntervalOfReadOnContig::additionalId());
    comparisonResult.additionalIntervalCount = additionalBounds.second - additionalBounds.first;

    // For each read: Get list of intervals that are to be found,
    // depending on options.  Then, look whether they are in result.
    for (size_t readId = 0; readId < length(value(witStore.readNames)); ++readId) {
        // We need a read id to search for in a sequence of
        // IntervalOfReadOnContig objects.  To do this, we create a
        // reference interval with the current read id.
        IntervalOfReadOnContig refInterval;
        refInterval.readId = readId;
        std::pair<TIntervalIterator, TIntervalIterator> boundsForRead;
        boundsForRead = std::equal_range(begin(witStore.intervals, Standard()), end(witStore.intervals, Standard()), refInterval, WitStoreLess<SortReadId>(witStore));
        if (boundsForRead.first == boundsForRead.second)
            continue;  // Skip if there are no intervals to be found for read.

        // Now, depending on the configured benchmark category, perform the
        // evaluation.
        if (options.benchmarkCategory == "all") {
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, options, CategoryAll());
        } else if (options.benchmarkCategory == "any-best") {
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, options, CategoryAnyBest());
        } else if (options.benchmarkCategory == "all-best") {
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, options, CategoryAllBest());
        } else {
            SEQAN_ASSERT_FAIL("Invalid benchmark category '%s'.", toCString(options.benchmarkCategory));
        }
    }
}
