// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): I think most time is spent reading the GSI file. We should be able to speed this up greatly with a more compact file format.
// TODO(holtgrew): Clean up, remove cruft.

#ifndef APPS_RABEMA_EVALUATION_H_
#define APPS_RABEMA_EVALUATION_H_

#include <seqan/basic.h>
#include <seqan/store.h>
#include <seqan/find.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/misc/misc_interval_tree.h>           // For interval trees.
#include <seqan/stream.h>
#include <seqan/bam_io.h>

#include "fai_index.h"

#include "rabema.h"

#include "wit_store.h"
#include "find_hamming_simple_ext.h"
#include "find_myers_ukkonen_reads.h"
#include "find_myers_ukkonen_ext.h"
#include "find_approx_dp_quality.h"
#include "find_hamming_simple_quality.h"
#include "evaluation_options.h"

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

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

    size_t totalReadCount;
    double normalizedIntervals;
    
    ComparisonResult()
            : totalIntervalCount(0), foundIntervalCount(0),
              superflousIntervalCount(0), additionalIntervalCount(0),
              totalReadCount(0), normalizedIntervals(0) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Resetting all counters to 0 for ComparisonResult.
void clear(ComparisonResult & result) {
    result.totalIntervalCount = 0;
    result.foundIntervalCount = 0;
    result.superflousIntervalCount = 0;
    result.additionalIntervalCount = 0;
    result.totalReadCount = 0;
    result.normalizedIntervals = 0;
}

// Output-to-stream operator for ComparisonResult.
template <typename TStream>
TStream & operator<<(TStream & stream, ComparisonResult const & result) {
    stream << "{\"total_intervals\": " << result.totalIntervalCount
           << ", \"found_intervals\": " << result.foundIntervalCount
           << ", \"superflous_intervals\": " << result.superflousIntervalCount
           << ", \"additional_intervals\": " << result.additionalIntervalCount
           << ", \"total_reads\": " << result.totalReadCount
           << ", \"normalized_intervals\": " << result.normalizedIntervals
           << "}";
    return stream;
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
                result -= static_cast<int>(ceil(static_cast<double>(x) / 2.0));
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
                        Options<EvaluateResults> const & options,
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
    Pattern<TReadSeq, TPatternSpec> pattern(read, -static_cast<int>(length(read)) * 1000);
    _patternMatchNOfPattern(pattern, options.matchN);
    _patternMatchNOfFinder(pattern, options.matchN);
    bool ret = setEndPosition(finder, pattern, endPos);
    (void)ret;  // If run without assertions.
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(endPos, endPosition(finder));
    
    // No explicit alignment is required if distances are not to be weighted.
    if (!options.weightedDistances)
        return getScore(pattern);

    // Otherwise, we need to build an alignment and compute the score from it.
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);

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
    return (int)floor(maxErrorRate / 100.0 * len);
}


template <typename TFragmentStore, typename TAlignedReadIter, typename TWitRecordIter, typename TContigSeq, typename TContigGaps, typename TPatternSpec>
void
compareAlignedReadsToReferenceOnContigForOneRead(Options<EvaluateResults> const & options,
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
    typedef typename TFragmentStore::TContigStore TContigStore;

//     std::cout << "wit records are" << std::endl;
//     for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
//         std::cout << value(it) << std::endl;
//     }
//     std::cout << "XXX" << std::endl;

    // Build scoring matrix that allows N to match with all.
    int gapExtensionScore = -1;
    int gapOpenScore = -1;
    if (IsSameType<TPatternSpec, HammingSimple>::VALUE) {
        // No gaps for hamming distance.
        gapOpenScore = -static_cast<int>(length(fragments.readSeqStore[readId]));
        gapExtensionScore = -static_cast<int>(length(fragments.readSeqStore[readId]));
    }
    // Build scoring matrix.
    Score<int, ScoreMatrix<Dna5> > matrixScore(gapExtensionScore, gapOpenScore);
    for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x) {
        for (int y = 0; y < ValueSize<Dna5>::VALUE; ++y)
            setScore(matrixScore, Dna5(x), Dna5(y), -1);
        setScore(matrixScore, Dna5(x), Dna5(x), 0);
    }
    // Update score matrix if N works as a wildcard character.
    if (options.matchN) {
        for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x) {
            setScore(matrixScore, Dna5(x), Dna5('N'), 0);
            setScore(matrixScore, Dna5('N'), Dna5(x), 0);
        }
    }
        
    //if (readId == 71)
    //    std::cout << "stop here" << std::endl;

    // Build interval tree.
    //std::cerr << ">>> readId == " << readId << ", contigId == " << contigId << std::endl;
    //std::cerr << "-----------" << std::endl;
    typedef IntervalAndCargo<size_t, size_t> TInterval;
    String<TInterval> intervals;
    for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
        // Skip aligned reads on other strand.
        //if (it->isForward != isForward) continue;
        // Skip intervals with too high distance, ignore distance in oracle wit mode.
        //std::cerr << "it->distance == " << it->distance << std::endl;
        if (!options.oracleWitMode && static_cast<int>(it->distance) > options.maxError) continue;
        // std::cerr << "insert(intervals, TInterval(" << value(it).firstPos << ", " << value(it).lastPos + 1 << ", " << value(it).id << "))" << std::endl;

        SEQAN_ASSERT_LEQ(value(it).firstPos, value(it).lastPos);
        appendValue(intervals, TInterval(value(it).firstPos, value(it).lastPos + 1, value(it).id));
    }
    IntervalTree<size_t, size_t> intervalTree(intervals, ComputeCenter());
    //std::cerr << "-----------" << std::endl;

    // Now, try to hit all entries in interval tree with the aligned reads on this strand.
    std::set<size_t> intervalsInResult;
    for (TAlignedReadIter it = alignedReadsBegin; it != alignedReadsEnd; ++it) {
        SEQAN_ASSERT_EQ(it->readId, readId);
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
//        if (endPos < length(fragments.readSeqStore[it->readId])) 
//            continue;

        // Skip reads that aligned with a too bad score.  Ignore score if in wit-oracle mode, i.e. comparing against simulated data.
        int bestScore = 1;  // Marker for "not computed, oracle wit mode."
        if (!options.oracleWitMode) {
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
        }

        // Compute last position of alignment and search for intervals this
        // position falls into.
        TPos lastPos = endPos - 1;
        SEQAN_ASSERT_LEQ(beginPos, lastPos);
        if (options.showTryHitIntervals) {
            std::cerr << "log> {\"type\": \"log.try_hit"
                      << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                      << "\", \"read_id\": \"" << fragments.readNameStore[it->readId]
                      << "\", \"strand\": \"" << (isForward ? "forward" : "reverse")
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

        // Query interval tree with lastPos.
        String<size_t> foundIntervalIds;
        //if (lastPos == 470875) {
        //  std::cout << "stop here" << std::endl;
        //  std::cout << "readId == " << readId << ", getMateNo() == " << getMateNo(fragments, readId) << std::endl;
        //}
        if (length(intervals) > 0u) {
            // std::cerr << "findIntervals(intervalTree, " << lastPos << ", foundIntervalIds);" << std::endl;
            findIntervals(intervalTree, lastPos, foundIntervalIds);
        }

        // Handle alignment out of target intervals.
        if (length(foundIntervalIds) == 0) {
            if (options.weightedDistances) {
                if (options.showAdditionalIntervals) {
                    std::cerr << "log> {\"type\": \"log.additional_hit"
                              << "\", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                              << "\", \"read_id\": \"" << fragments.readNameStore[it->readId]
                              << "\", \"distance\": " << -bestScore
                              << ", \"strand\": \"" << (isForward ? "forward" : "reverse")
                              << "\", \"begin_pos\": " << beginPos
                              << ", \"end_pos\": " << endPos
                              << ", \"mate_no\": " << getMateNo(fragments, readId)
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
                appendValue(result, IntervalOfReadOnContig::additionalId());
                continue;
            } else { // if (!options.weightedDistances) {
                if (options.dontPanic)
                    std::cerr << "WARNING: ";
                else
                    std::cerr << "PANIC: ";
                std::cerr << "A read in the Sam file aligns out of all target intervals for this read in the WIT file." << std::endl;
                std::cerr << "bestScore = " << bestScore << std::endl;
                std::cerr << "read name = " << fragments.readNameStore[it->readId] << std::endl;
                std::cerr << "read is = " << fragments.readSeqStore[it->readId] << std::endl;
                Dna5String rcRead(fragments.readSeqStore[it->readId]);
                reverseComplement(rcRead);
                std::cerr << "          " << rcRead << std::endl;
                std::cerr << "mate no is = " << static_cast<int>(getMateNo(fragments, it->readId)) << std::endl;
                std::cerr << "on forward strand? " << isForward << std::endl;
                std::cerr << "original begin pos = " << it->beginPos << std::endl;
                std::cerr << "original end pos = " << it->endPos << std::endl;
                std::cerr << "begin pos = " << beginPos << std::endl;
                std::cerr << "end pos = " << endPos << std::endl;
                std::cerr << "last pos = " << lastPos << std::endl;
                std::cerr << "contigId = " << it->contigId << std::endl;
                std::cerr << "mateNo " << getMateNo(fragments, it->readId) << std::endl;
                std::cerr << "max error rate is " << options.maxError << std::endl;
                std::cerr << "contig_infix_seq = " << infix(contig, beginPos, endPos) << std::endl;
                std::cerr << "read length is " << length(fragments.readSeqStore[it->readId]) << std::endl;
                std::cerr << "max errors is " <<  maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId])) << std::endl;
                if (options.oracleWitMode)
                    std::cerr << "NOTE: Oracle WIT mode is enabled!" << std::endl;
                else if (!options.dontPanic)
                    exit(1);
                appendValue(result, IntervalOfReadOnContig::additionalId());
            }
        }

        // Record found intervals.
        typedef typename Iterator<String<size_t>, Standard>::Type TFoundIdsIterator;
        for (TFoundIdsIterator iter = begin(foundIntervalIds, Standard()); iter != end(foundIntervalIds, Standard()); ++iter) {
            // Count interval if hit for the first time.
            if (intervalsInResult.find(value(iter)) == intervalsInResult.end()) {
                appendValue(result, value(iter));
                intervalsInResult.insert(value(iter));
            }
        }
    }

    // If configured so, show hit and missed intervals.
    if (options.showHitIntervals) {
        for (TWitRecordIter it = witRecordsBegin; it != witRecordsEnd; ++it) {
            bool found = intervalsInResult.find(value(it).id) != intervalsInResult.end();
			if (found && options.showHitIntervals) {
                std::cerr << "log> {\"type\": \"log.hit_interval"
                          << "\", \"interval_id\": " << value(it).id
                          << ", \"contig_id\": \"" << fragments.contigNameStore[contigId]
                          << "\", \"strand\": \"" << (isForward ? "forward" : "reverse")
                          << "\", \"read_id\": \"" << fragments.readNameStore[readId]
                          << "\", \"interval_first\": " << value(it).firstPos
                          << ", \"interval_last\": " << value(it).lastPos << "}" << std::endl;
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
compareAlignedReadsToReferenceOnContig(Options<EvaluateResults> const & options,
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
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContigStoreElement::TGapAnchors> > TContigGaps;

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
//         std::cout << "read id " << fragments.readNameStore[readId] << std::endl;
        compareAlignedReadsToReferenceOnContigForOneRead(options, fragments, contigId, contig, contigGaps, isForward, readId, alignedReadsBeginForRead, alignedReadsEndForRead, witRecordsBeginForRead, witRecordsEndForRead, result, TPatternSpec());
        // This iteration's end iterators are the next iteration's begin iterators.
        alignedReadsBeginForRead = alignedReadsEndForRead;
        witRecordsBeginForRead = witRecordsEndForRead;
    }
}

// TODO(holtgrew): Move this into its own header?

class RefIdMapping
{
public:
    String<unsigned> map;
    
    RefIdMapping() {}
};

inline unsigned length(RefIdMapping const & mapping)
{
    return length(mapping.map);
}

template <typename TTargetNameStore, typename TTargetNameStoreCache, typename TSourceNameStore>
void rebuildMapping(RefIdMapping & mapping,
                    TTargetNameStore const & targetNameStore,
                    TTargetNameStoreCache const & targetNameStoreCache,
                    TSourceNameStore const & sourceNameStore)
{
    clear(mapping.map);
    resize(mapping.map, length(sourceNameStore), maxValue<unsigned>());

    for (unsigned i = 0; i < length(sourceNameStore); ++i)
    {
        unsigned idx = 0;
        if (getIdByName(targetNameStore, sourceNameStore[i], idx, targetNameStoreCache))
            mapping.map[i] = idx;
    }
}

struct RabemaStats
{
    // Number of intervals that were to find.
    __uint64 intervalsToFind;
    // Number of intervals that were actually found.
    __uint64 intervalsFound;
    // SAM records that did not correspond to alignments below the configured maximal error rate.
    __uint64 invalidAlignments;

    // Total number of reads that we have GSI records for, equals number of normalized intervals to be found.
    __uint64 totalReads;
    // Normalized number of found intervals.
    double normalizedIntervals;

    // The following arrays are indexed by the integer value of the error rate.

    // Number of intervals that were to find for each error rate.
    String<unsigned> intervalsToFindForErrorRate;
    // Number of found intervals for each error rate.
    String<unsigned> intervalsFoundForErrorRate;

    // The following values are normalized towards all intervals.

    // Normalized number of intervals to find for each error rate.
    String<double> normalizedIntervalsToFindForErrorRate;
    // Normalized number of intervals found for each error rate.
    String<double> normalizedIntervalsFoundForErrorRate;

    RabemaStats() : intervalsToFind(0), intervalsFound(0), invalidAlignments(0), totalReads(0), normalizedIntervals(0)
    {}

    RabemaStats(unsigned maxErrorRate) :
            intervalsToFind(0), intervalsFound(0), invalidAlignments(0), totalReads(0), normalizedIntervals(0)
    {
        resize(intervalsToFindForErrorRate, maxErrorRate + 1, 0);
        resize(intervalsFoundForErrorRate, maxErrorRate + 1, 0);
        resize(normalizedIntervalsToFindForErrorRate, maxErrorRate + 1, 0.0);
        resize(normalizedIntervalsFoundForErrorRate, maxErrorRate + 1, 0.0);
    }
};

void updateMaximalErrorRate(RabemaStats & stats, unsigned maxErrorRate)
{
    if (length(stats.intervalsToFindForErrorRate) <= maxErrorRate + 1)
        resize(stats.intervalsToFindForErrorRate, maxErrorRate + 1, 0);
    if (length(stats.intervalsFoundForErrorRate) <= maxErrorRate + 1)
        resize(stats.intervalsFoundForErrorRate, maxErrorRate + 1, 0);
    if (length(stats.normalizedIntervalsToFindForErrorRate) <= maxErrorRate + 1)
        resize(stats.normalizedIntervalsToFindForErrorRate, maxErrorRate + 1, 0.0);
    if (length(stats.normalizedIntervalsFoundForErrorRate) <= maxErrorRate + 1)
        resize(stats.normalizedIntervalsFoundForErrorRate, maxErrorRate + 1, 0.0);
}

template <typename TStream>
void write(TStream & stream, RabemaStats const & stats, Options<EvaluateResults> const & options)
{
    stream << "Intervals to find:              " << stats.intervalsToFind << '\n'
           << "Intervals found:                " << stats.intervalsFound << '\n'
           << "Intervals found [%]             " << (100.0 * stats.intervalsFound / stats.intervalsToFind) << '\n'
           << "Invalid alignments:             " << stats.invalidAlignments << '\n'
           << '\n'
           << "Number of reads:                " << stats.totalReads << '\n'
           << "Normalized intervals found:     " << stats.normalizedIntervals << '\n'
           << "Normalized intervals found [%]: " << (100.0 * stats.normalizedIntervals / stats.totalReads) << '\n'
           << '\n';
    char buffer[1000];
    sprintf(buffer, "  ERR\t%8s\t%8s\t%8s\t%8s\t%10s\t%10s\n", "#max", "#found", "%found", "norm max", "norm found", "norm found [%]");
    stream << buffer;
    stream << "------------------------------------------------------------------------------------------------------\n";
    for (unsigned i = 0; i < length(stats.intervalsToFindForErrorRate); ++i)
    {
        if (!options.oracleWitMode && options.benchmarkCategory == "all" && (int)i != options.maxError)
            continue;
        sprintf(buffer, "%5u\t%8d\t%8d\t%8.2f\t%8.2f\t%10.2f\t%10.2f\n", i, stats.intervalsToFindForErrorRate[i], stats.intervalsFoundForErrorRate[i],
                100.0 * stats.intervalsFoundForErrorRate[i] / stats.intervalsToFindForErrorRate[i],
                stats.normalizedIntervalsToFindForErrorRate[i], stats.normalizedIntervalsFoundForErrorRate[i],
                100.0 * stats.normalizedIntervalsFoundForErrorRate[i] / stats.normalizedIntervalsToFindForErrorRate[i]);
        stream << buffer;
    }
    stream << '\n';
}

struct CmpWitRecordLowering
{
    // Lexicographically sort by (readId, contigId, first pos).
    bool operator()(WitRecord const & lhs, WitRecord const & rhs) const
    {
        return (lhs.readId < rhs.readId) || (lhs.readId == rhs.readId && lhs.contigId < rhs.contigId) ||
                (lhs.readId == rhs.readId && lhs.contigId == rhs.contigId && lhs.firstPos < rhs.firstPos);
    }
};

void performIntervalLowering(String<WitRecord> & gsiRecords, int maxError)
{
    if (empty(gsiRecords))
        return;

    typedef Iterator<String<WitRecord> >::Type TIterator;
    typedef IntervalAndCargo<unsigned, unsigned> TInterval;

    // Step 1: Adjust distances.
    std::sort(begin(gsiRecords, Standard()), end(gsiRecords, Standard()), CmpWitRecordLowering());

    // Add sentinel interval.
    WitRecord sentinel(back(gsiRecords));
    sentinel.firstPos = MaxValue<size_t>::VALUE;
    sentinel.lastPos = MaxValue<size_t>::VALUE;
    // sentinel.id = MaxValue<size_t>::VALUE;
    appendValue(gsiRecords, sentinel);

    String<TInterval> openIntervals;
    unsigned i = 0;
    for (TIterator it = begin(gsiRecords, Standard()), itEnd = end(gsiRecords, Standard()); it != itEnd; ++it, ++i)
    {
        unsigned count = 0;
        for (unsigned j = 0; j < length(openIntervals); ++j)
        {
            unsigned idx = length(openIntervals) - 1 - j;
            WitRecord const thisIntervalRecord = gsiRecords[cargo(openIntervals[idx])];
            SEQAN_ASSERT_EQ(thisIntervalRecord.readName, it->readName);
            if (thisIntervalRecord.contigId != it->contigId || thisIntervalRecord.lastPos < it->firstPos)
                count += 1;
        }
        resize(openIntervals, length(openIntervals) - count);

        // Perform distance lowering for containing intervals.
        for (unsigned j = 0; j < length(openIntervals); ++j)
        {
            unsigned idx = length(openIntervals) - 1 - j;
            unsigned id = cargo(openIntervals[idx]);
            if (gsiRecords[id].distance <= maxError)
                gsiRecords[id].distance = _min(gsiRecords[id].distance, it->distance);
            else
                break;  // All containing intervals must have a greater distance.
        }

        appendValue(openIntervals, TInterval(it->firstPos, it->lastPos + 1, i));
    }

    // Step 2: Filter out intervals that are contained in intervals of lesser/equal distance.
    String<WitRecord> filteredGsiRecords;
    clear(openIntervals);
    i = 0;
    for (TIterator it = begin(gsiRecords, Standard()), itend = end(gsiRecords, Standard()); it != itend; ++it, ++i)
    {
        // Remove non-overlapping intervals on top of openInterval stack, appending to filtered intervals
        unsigned count = 0;
        for (unsigned j = 0; j < length(openIntervals); ++j) {
            unsigned idx = length(openIntervals) - 1 - j;
            WitRecord const & thisIntervalRecord = gsiRecords[cargo(openIntervals[idx])];
            SEQAN_ASSERT_EQ(thisIntervalRecord.readName, it->readName);
            if (thisIntervalRecord.contigId != it->contigId || thisIntervalRecord.lastPos < it->firstPos)
            {
                count += 1;
                unsigned startDistance = gsiRecords[cargo(openIntervals[idx])].distance;
                if (!empty(filteredGsiRecords)) {
                    if (back(filteredGsiRecords).lastPos >= leftBoundary(openIntervals[idx]))
                    {
                        if (back(filteredGsiRecords).contigId == gsiRecords[cargo(openIntervals[idx])].contigId)
                        {
                            // Assert current containing already written out.
                            SEQAN_ASSERT_GEQ(back(filteredGsiRecords).firstPos, leftBoundary(openIntervals[idx]));
                            SEQAN_ASSERT_LEQ(back(filteredGsiRecords).lastPos + 1, rightBoundary(openIntervals[idx]));
                            // Get start distance.
                            startDistance = back(filteredGsiRecords).distance + 1;
                        }
                    }
                }
                unsigned upperLimit = maxError;
                if ((unsigned)maxError < startDistance)
                    upperLimit = startDistance;
                for (unsigned i = startDistance; i <= upperLimit; ++i)
                {
                    appendValue(filteredGsiRecords, gsiRecords[cargo(openIntervals[idx])]);
                    // std::cerr << "APPENDING " << cargo(openIntervals[idx]) << "\t" << idx << "\n"
                    //           << gsiRecords[cargo(openIntervals[idx])] << "\n"
                    //           << "now is\n"
                    //           << back(filteredGsiRecords) << "\t" << back(filteredGsiRecords).originalDistance << "\n";
                    back(filteredGsiRecords).distance = i;
                }
            }
        }
        resize(openIntervals, length(openIntervals) - count);

        // Add interval to the stack of intervals.
        if (empty(openIntervals) || gsiRecords[cargo(back(openIntervals))].distance > it->distance)
            appendValue(openIntervals, TInterval(it->firstPos, it->lastPos + 1, i));
    }
    move(gsiRecords, filteredGsiRecords);
}


template <typename TPatternSpec>
int benchmarkReadResult(RabemaStats & result,
                        String<BamAlignmentRecord> const & samRecords,
                        BamIOContext<StringSet<CharString > > const & bamIOContext,
                        String<WitRecord> const & gsiRecords,
                        StringSet<CharString> const & refSeqNames,
                        NameStoreCache<StringSet<CharString> > const & refSeqNamesCache,
                        StringSet<Dna5String> const & refSeqs,
                        RefIdMapping const & refIdMapping,
                        Options<EvaluateResults> const & options,
                        TPatternSpec const & /*tagPattern*/,
                        bool pairedEnd = false,
                        bool second = false)
{
    typedef IntervalAndCargo<unsigned, unsigned> TInterval;

    // TODO(holtgrew): While much clearer and simpler than the old version, I'm not quite sure that this one is bug free yet.
    // TODO(holtgrew): What about additional interval logging and counting, those that SHOULD be in the standard but are not there?

    // Select gold standard intervals (GSI) records.
    //
    // We only select intervals that match the specification of pairedEnd, second.  Also, the intervals must be for an
    // error rate less than or equal to the maximal configured one.  In case of any-best or all-best mode, we only
    // select those of the eligible intervals with the lowest error rate.
    //
    // In case of oracle mode, we ignore the distance of the intervals in the GSI file here but use it later on.
    //
    // Start with picking the smallest distance if *-best mode.
    int smallestDistance = options.oracleWitMode ? maxValue<int>() : options.maxError;
    if (options.oracleWitMode || options.benchmarkCategory == "any-best" || options.benchmarkCategory == "all-best")
        for (unsigned i = 0; i < length(gsiRecords); ++i)
            smallestDistance = std::min(smallestDistance, gsiRecords[i].distance);
    int largestDistance = options.maxError;
    if (options.oracleWitMode && smallestDistance != maxValue<int>())
        largestDistance = smallestDistance;
    // if (options.oracleWitMode)
    //     smallestDistance = 0;
    String<WitRecord> pickedGsiRecords;
    for (unsigned i = 0; i < length(gsiRecords); ++i)
    {
        // Note: In case of oracle mode, we ignore the distance.
        if (!options.oracleWitMode && gsiRecords[i].distance > smallestDistance)
            continue;  // Skip with wrong distance.
        if (!options.oracleWitMode && gsiRecords[i].distance > options.maxError)
            continue;  // Ignore intervals with too high error rate.
        if (!pairedEnd && (gsiRecords[i].flags & WitRecord::FLAG_PAIRED))
            continue;  // Skip paired if non-paired selected.
        if (pairedEnd && !(gsiRecords[i].flags & WitRecord::FLAG_PAIRED))
            continue;  // Skip non-paired if paired selected.
        if (pairedEnd && second && !(gsiRecords[i].flags & WitRecord::FLAG_SECOND_MATE))
            continue;  // Skip if second selected but this interval is not for second.

        appendValue(pickedGsiRecords, gsiRecords[i]);

        // Get index of the sequence from GSI record contig name.
        if (!getIdByName(refSeqNames, back(pickedGsiRecords).contigName, back(pickedGsiRecords).contigId, refSeqNamesCache))
        {
            std::cerr << "ERROR: Could not find reference sequence for name " << back(pickedGsiRecords).contigName << '\n';
            return 1;
        }
    }

    // On these selected GSI records, we now perform interval lowering in case of "any-best" and "all-best".  This means
    // that an interval I with distance k_i containing an interval J with distance k_j < k_i is re-labeled with a
    // distance of k_j for the smallest k_j of all contained intervals J.
    if (!options.oracleWitMode && (options.benchmarkCategory == "any-best" || options.benchmarkCategory == "all-best"))
    {
        // std::sort(begin(pickedGsiRecords, Standard()), end(pickedGsiRecords, Standard()), CmpWitRecordLowering());
        // std::cerr << ",-- pre-lowering\n";
        // for (unsigned i = 0; i < length(pickedGsiRecords); ++i)
        //     std::cerr << "| " << pickedGsiRecords[i] << "\n";
        // std::cerr << "`--\n";
        performIntervalLowering(pickedGsiRecords, options.maxError);
        // std::sort(begin(pickedGsiRecords, Standard()), end(pickedGsiRecords, Standard()), CmpWitRecordLowering());
        // std::cerr << ",-- post-lowering\n";
        // for (unsigned i = 0; i < length(pickedGsiRecords); ++i)
        //     std::cerr << "| " << pickedGsiRecords[i] << "\n";
        // std::cerr << "`--\n";
    }

    // Build string of intervals for each reference sequence from these filtered and lowered records.
    String<String<TInterval> > intervals;
    resize(intervals, length(refSeqs));
    String<unsigned> numIntervalsForErrorRate;
    resize(numIntervalsForErrorRate, options.maxError + 1, 0);
    String<int> intervalDistances;  // Distance of interval i.
    for (unsigned i = 0; i < length(pickedGsiRecords); ++i)
    {
        int distance = pickedGsiRecords[i].distance;
        if (!options.oracleWitMode && options.benchmarkCategory != "all" && distance != options.maxError)
            continue;  // Only search if interval distance matches -- in *-best.
        int originalDistance = pickedGsiRecords[i].originalDistance;

        appendValue(intervals[pickedGsiRecords[i].contigId], TInterval(pickedGsiRecords[i].firstPos, pickedGsiRecords[i].lastPos + 1, length(intervalDistances)));
        appendValue(intervalDistances, originalDistance);
        if (!options.oracleWitMode && options.benchmarkCategory != "any-best")
            numIntervalsForErrorRate[originalDistance] += 1;
    }
    // Marker array that states whether an interval was hit.
    String<bool> intervalHit;
    resize(intervalHit, length(intervalDistances), false);

    // Build interval trees.
    String<IntervalTree<unsigned> > intervalTrees;
    resize(intervalTrees, length(refSeqs));
    for (unsigned i = 0; i < length(intervalTrees); ++i)
        createIntervalTree(intervalTrees[i], intervals[i]);

    // One of the SAM records must have a non-empty SEQ field, extract the read seq.
    Dna5String readSeqL, readSeqR;
    bool seenL = false, seenR = false;
    for (unsigned i = 0; i < length(samRecords); ++i)
    {
        seenL |= hasFlagFirst(samRecords[i]) || (!hasFlagFirst(samRecords[i]) && !hasFlagLast(samRecords[i]));
        seenR |= hasFlagLast(samRecords[i]);
        if ((hasFlagFirst(samRecords[i]) || (!hasFlagFirst(samRecords[i]) && !hasFlagLast(samRecords[i]))) && !empty(samRecords[i].seq))
        {
            readSeqL = samRecords[i].seq;
            if (hasFlagRC(samRecords[i]))
                reverseComplement(readSeqL);
        }
        else if (hasFlagLast(samRecords[i]) && !empty(samRecords[i].seq))
        {
            readSeqR = samRecords[i].seq;
            if (hasFlagRC(samRecords[i]))
                reverseComplement(readSeqR);
        }
        if (!empty(readSeqL) && !empty(readSeqR))
            break;  // Short-circuit and break.
    }
    if (seenL && empty(readSeqL))
    {
        std::cerr << "ERROR: No alignment for query " << front(samRecords).qName << " (left-end)\n";
        return 1;
    }
    if (seenR && empty(readSeqR))
    {
        std::cerr << "ERROR: No alignment for query " << front(samRecords).qName << " (right-end)\n";
        return 1;
    }

    // Try to hit intervals.
    // std::cerr << "NUM SAM RECORDS\t" << length(samRecords) << '\n';
    for (unsigned i = 0; i < length(samRecords); ++i)
    {
        BamAlignmentRecord const & samRecord = samRecords[i];
        int seqId = refIdMapping.map[samRecord.rId];

        // Compute actual alignment score to rule out invalidly reported alignments.
        //
        // If we run in oracle mode then we ignore the actual alignment score and use the best alignment.  We only care
        // about the alignment's end position in this case.
        if (!options.oracleWitMode)
        {
            int bestDistance = minValue<int>();  // Marker for "not set yet".
            // Get best distance from NM tag if set and we are to trust it.
            if (options.trustNM)
            {
                // TODO(holtgrew): Remove const cast once we have const holders!
                BamTagsDict bamTags(const_cast<CharString &>(samRecord.tags));
                unsigned idx = 0;
                if (findTagKey(idx, bamTags, "NM"))
                {
                    if (!extractTagValue(bestDistance, bamTags, idx))
                        bestDistance = minValue<int>();  // Reset to sentinel.
                }
            }
            // Otherwise, perform a realignment.
            Dna5String readSeq;
            Dna5String contigSeq;
            if (bestDistance == minValue<int>())
            {
                if (hasFlagFirst(samRecord) || (!hasFlagFirst(samRecord) && !hasFlagLast(samRecord)))
                    readSeq = readSeqL;
                if (hasFlagLast(samRecord))
                    readSeq = readSeqR;
                int bandwidth = static_cast<int>(ceil(0.01 * options.maxError * length(readSeq)));
                __int32 beginPos = samRecord.pos - bandwidth;
                __int32 endPos = beginPos + getAlignmentLengthInRef(samRecord) - countPaddings(samRecord.cigar) + 2 * bandwidth;
                // if (hasFlagRC(samRecord))
                // {
                //     endPos = length(refSeqs[seqId]) - endPos;
                //     beginPos = length(refSeqs[seqId]) - beginPos;
                //     std::swap(endPos, beginPos);
                // }
                contigSeq = infix(refSeqs[seqId], beginPos, endPos);
                if (hasFlagRC(samRecord))
                    reverseComplement(contigSeq);
                Finder<Dna5String> finder(contigSeq);
                Pattern<Dna5String, TPatternSpec> pattern(readSeq, -static_cast<int>(length(readSeq)) * 1000);
                _patternMatchNOfPattern(pattern, options.matchN);
                _patternMatchNOfFinder(pattern, options.matchN);
                bool ret = setEndPosition(finder, pattern, length(contigSeq) - bandwidth);
                (void) ret;  // When run without assertions.
                SEQAN_CHECK(ret, "setEndPosition() must not fail!");
                bestDistance = -getScore(pattern);
            }

            // Skip invalid alignments.
            if (bestDistance > static_cast<int>(0.01 * options.maxError * length(readSeq)))
            {
                if (options.showSuperflousIntervals)
                {
                    std::cerr << "SUPERFLOUS/INVALID\t";
                    write2(std::cerr, samRecord, bamIOContext, Sam());
                    std::cerr << "  READ:  \t" << readSeq << '\n'
                              << "  CONTIG:\t" << contigSeq << '\n'
                              << "  DISTANCE:        \t" << bestDistance << '\n'
                              << "  ALLOWED DISTANCE:\t" << static_cast<int>(0.01 * options.maxError * length(readSeq)) << '\n';
                }
                result.invalidAlignments += 1;
                continue;
            }
        }

        // std::cerr << "SAM RECORD\t";
        // write2(std::cerr, samRecord, bamIOContext, Sam());
        // std::cerr << '\n';

        // Get sequence id and last position of alignment.  We try to hit the interval with the last position (not
        // C-style end) of the read.
        unsigned lastPos = 0;
        if (!hasFlagRC(samRecord))
            lastPos = samRecord.pos + getAlignmentLengthInRef(samRecord) - countPaddings(samRecord.cigar) - 1;
        else
            lastPos = length(refSeqs[seqId]) - samRecord.pos - 1;
        
        if (options.showTryHitIntervals)
            std::cerr << "TRY HIT\tchr=" << refSeqNames[seqId] << "\tlastPos=" << lastPos << "\tqName=" << samRecord.qName << "\n";
        
        // Try to hit any interval.
        String<unsigned> result;
        findIntervals(intervalTrees[seqId], lastPos, result);
        for (unsigned i = 0; i < length(result); ++i)
            intervalHit[result[i]] = true;
    }

    // Compute number of found intervals.
    unsigned numFound = 0;
    String<unsigned> foundIntervalsForErrorRate;
    if ((int)length(foundIntervalsForErrorRate) <= largestDistance + 1)
        resize(foundIntervalsForErrorRate, largestDistance + 1, 0);
    if (options.oracleWitMode || options.benchmarkCategory == "any-best")
    {
        int bestDistance = maxValue<int>();
        for (unsigned i = 0; i < length(intervalDistances); ++i)
            if (intervalHit[i])
            {
                if (options.showHitIntervals)
                    std::cerr << "HIT\t" << pickedGsiRecords[i] << "\n";
                bestDistance = std::min(bestDistance, intervalDistances[i]);
            }
        if (bestDistance != maxValue<int>())
        {
            numFound += 1;
            foundIntervalsForErrorRate[bestDistance] += 1;
        }
    }
    else  // !options.oracleWitMode && options.benchmarkCategory in ["all-best", "all"]
    {
        for (unsigned i = 0; i < length(intervalDistances); ++i)
        {
            if (options.benchmarkCategory == "all" && intervalDistances[i] != options.maxError)
                continue;  // Only count intervals on our maximal error rate in "all" mode.
            if (intervalHit[i])
            {
                if (options.showHitIntervals)
                    std::cerr << "HIT\t" << pickedGsiRecords[i] << "\n";
                numFound += 1;
                foundIntervalsForErrorRate[intervalDistances[i]] += 1;
            }
            else
            {
                if (options.showMissedIntervals)  // inside braces for consistency with above
                    std::cerr  << "MISSED\t" << pickedGsiRecords[i] << "\n";
            }
        }
        SEQAN_ASSERT_LEQ(numFound, length(intervalDistances));
    }
    
    // Update the resulting RabemaStats.
    updateMaximalErrorRate(result, largestDistance);
    result.totalReads += 1;
    if (options.oracleWitMode || options.benchmarkCategory == "any-best")
    {
        bool found = (numFound > 0u);
        result.intervalsToFind += 1;
        result.intervalsFound += found;
        result.normalizedIntervals += found;
        int d = (smallestDistance == maxValue<int>()) ? 0 : smallestDistance;
        result.intervalsToFindForErrorRate[d] += 1;
        result.intervalsFoundForErrorRate[d] += found;
        result.normalizedIntervalsToFindForErrorRate[d] += 1;
        result.normalizedIntervalsFoundForErrorRate[d] += found;
    }
    else  // all-best or all was selected
    {
        unsigned intervalsToFind = 0;
        unsigned intervalsFound = 0;
        for (unsigned d = 0; d < length(numIntervalsForErrorRate); ++d)
        {
            // In case of "all", we only count the intervals from with maximal error rate.
            if (options.benchmarkCategory != "all" || (int)d == options.maxError)
            {
                intervalsToFind += numIntervalsForErrorRate[d];
                intervalsFound += foundIntervalsForErrorRate[d];;
                result.intervalsToFindForErrorRate[d] += numIntervalsForErrorRate[d];
                result.intervalsFoundForErrorRate[d] += foundIntervalsForErrorRate[d];
            }
        }
        result.intervalsToFind += intervalsToFind;
        result.intervalsFound += intervalsFound;
        result.normalizedIntervals += 1.0 * intervalsFound / intervalsToFind;
        for (unsigned d = 0; d < length(numIntervalsForErrorRate); ++d)
        {
            // In case of "all", we only count the intervals from with maximal error rate.
            if (options.benchmarkCategory != "all" || (int)d == options.maxError)
            {
                if (intervalsToFind > 0u)
                {
                    result.normalizedIntervalsToFindForErrorRate[d] += 1.0 * numIntervalsForErrorRate[d] / intervalsToFind;
                    result.normalizedIntervalsFoundForErrorRate[d] += 1.0 * foundIntervalsForErrorRate[d] / intervalsToFind;
                }
            }
        }
    }

    return 0;
}

// Stream over both the SAM and GSI file and compare the hits in the SAM file against the intervals in the GSI file.
//
// Both the SAM file and the GSI file have to be sorted by queryname for this to work.

template <typename TSamStream, typename TSamReaderSpec, typename TGsiStream, typename TGsiStreamSpec, typename TPatternSpec>
int
compareAlignedReadsToReference(RabemaStats & result,
                               RecordReader<TSamStream, TSamReaderSpec> & samReader,
                               BamIOContext<StringSet<CharString> > & bamIOContext,
                               StringSet<CharString> const & refNameStore,
                               NameStoreCache<StringSet<CharString> > const & refNameStoreCache,
                               StringSet<Dna5String> const & refSeqs,
                               RecordReader<TGsiStream, TGsiStreamSpec> & gsiReader,
                               Options<EvaluateResults> const & options,
                               TPatternSpec const & tagPattern)
{
    // Read in initial SAM/GSI records.
    BamAlignmentRecord samRecord;
    if (atEnd(samReader) || readRecord(samRecord, bamIOContext, samReader, Sam()) != 0)
    {
        std::cerr << "ERROR: Could not read first SAM record.\n";
        return 1;
    }
    WitRecord gsiRecord;
    if (atEnd(gsiReader) || readRecord(gsiRecord, gsiReader, Gsi()) != 0)
    {
        std::cerr << "ERROR: Could not read first GSI record.\n";
        return 1;
    }

    // Mapping between ref IDs from SAM file and reference sequence (from SAM file to reference sequences).
    RefIdMapping refIdMapping;
    rebuildMapping(refIdMapping, refNameStore, refNameStoreCache, nameStore(bamIOContext));

    // Current SAM and GSI records are stored in these arrays.
    String<BamAlignmentRecord> currentSamRecords;
    String<WitRecord> currentGsiRecords;

    // These flags store whether we processed the last SAM/GSI record.
    bool samDone = false, gsiDone = false;

    // The main loop: We walk over both the SAM and GSI records.
    // unsigned chunkI = 0;
    std::cerr << "Each dot corresponds to 10k processed reads.\n"
              << "\n"
              << "Progress: ";
    unsigned i = 0;
    while (!samDone || !gsiDone)
    {
        if (i > 0u && i % (100*1000) == 0u)
            std::cerr << i / 100 / 1000 << "00k";
        else if (i > 0 && i % (10*1000) == 0u)
            std::cerr << '.';
        ++i;
        
        // We process the record for the next query/read.  Since records for this next query/read might be missing in
        // both files, we need to determine which is the next one.
        CharString currentReadName;
        if (gsiDone)
            currentReadName = samRecord.qName;
        else if (samDone)
            currentReadName = gsiRecord.readName;
        else
            currentReadName = (gsiRecord.readName < samRecord.qName) ? gsiRecord.readName : samRecord.qName;
        // std::cerr << "CURRENT\t" << currentReadName << '\n';

        // These flags determine whether evaluation is run for single-end and/or paired-end reads.
        bool seenSingleEnd = false, seenPairedEnd = false;

        // Read all SAM records with the same query name.
        clear(currentSamRecords);
        // std::cerr << "SAM\t" << samRecord.qName << '\n';
        while (!samDone && samRecord.qName == currentReadName)
        {
            if (!hasFlagUnmapped(samRecord))  // Ignore records with non-aligned reads.
            {
                seenSingleEnd |= !hasFlagMultiple(samRecord);
                seenPairedEnd |= hasFlagMultiple(samRecord);
                appendValue(currentSamRecords, samRecord);
            }
            if (atEnd(samReader))
            {
                // At end of SAM File, do not read next one.
                samDone = true;
                continue;
            }
            if (readRecord(samRecord, bamIOContext, samReader, Sam()) != 0)
            {
                std::cerr << "ERROR: Could not read SAM record.\n";
                return 1;
            }
            if (samRecord.qName < currentReadName)
            {
                std::cerr << "ERROR: Wrong order in SAM file: " << samRecord.qName << " succeeds " << currentReadName << " in file.\n";
                return 1;
            }
            // Rebuild ref ID mapping if we discovered a new reference sequence.
            if (length(nameStore(bamIOContext)) != length(refIdMapping))
                rebuildMapping(refIdMapping, refNameStore, refNameStoreCache, nameStore(bamIOContext));
        }

        // Read in the next block of GSI records.
        // std::cerr << "GSI\t" << gsiRecord.readName << '\n';
        clear(currentGsiRecords);
        while (!gsiDone && gsiRecord.readName == currentReadName)
        {
            seenSingleEnd |= !(gsiRecord.flags & WitRecord::FLAG_PAIRED);
            seenPairedEnd |= (gsiRecord.flags & WitRecord::FLAG_PAIRED);
            appendValue(currentGsiRecords, gsiRecord);
            if (atEnd(gsiReader))
            {
                // At end of GSI File, do not read next one.
                gsiDone = true;
                continue;
            }
            if (readRecord(gsiRecord, gsiReader, Gsi()) != 0)
            {
                std::cerr << "ERROR: Could not read GSI record.\n";
                return 1;
            }
            if (gsiRecord.readName < currentReadName)
            {
                std::cerr << "ERROR: Wrong order in GSI file: " << gsiRecord.readName << " succeeds " << currentReadName << " in file.\n";
                return 1;
            }
        }

        // Now, compare the SAM records against the intervals stored in the GSI records.
        //
        // We collected the records for all queries.  Here, we differentiate between the different cases.
        if (seenSingleEnd)
        {
            benchmarkReadResult(result, currentSamRecords, bamIOContext, currentGsiRecords, refNameStore, refNameStoreCache, refSeqs, refIdMapping, options, tagPattern, /*pairedEnd=*/false);
        }
        if (seenPairedEnd)
        {
            benchmarkReadResult(result, currentSamRecords, bamIOContext, currentGsiRecords, refNameStore, refNameStoreCache, refSeqs, refIdMapping, options, tagPattern, /*pairedEnd=*/true, /*second=*/false);
            benchmarkReadResult(result, currentSamRecords, bamIOContext, currentGsiRecords, refNameStore, refNameStoreCache, refSeqs, refIdMapping, options, tagPattern, /*pairedEnd=*/true, /*second=*/true);
        }
    }
    std::cerr << " DONE\n";

    return 0;
}


// Tags for selecting an evaluateFoundIntervals_compareToIntervals specialization.
struct CategoryAll_;
typedef Tag<CategoryAll_> CategoryAll;
struct CategoryAnyBest_;
typedef Tag<CategoryAnyBest_> CategoryAnyBest;
struct CategoryAllBest_;
typedef Tag<CategoryAllBest_> CategoryAllBest;


// result must contain all ids of the intervals in
// [beginIntervalsForRead, endIntervalsForRead) that have the error
// rate options.maxError.  The intervals are sorted by (read id,
// error rate).
template <typename TFragmentStore>
void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               TFragmentStore const & /*fragments*/,
                                               Options<EvaluateResults> const & options,
                                               CategoryAll const &) {
    if (beginIntervalsForRead == endIntervalsForRead)
        return;  // Guard: Skip if interval range empty.

    typedef Iterator<WitStore::TIntervalStore, Standard>::Type TIntervalIterator;

    // Now, get begin and end iterators to the intervals for the error rate
    // configured in options, 0 in oracle wit mode.
    IntervalOfReadOnContig refInterval;
    refInterval.distance = options.oracleWitMode ? 0 : options.maxError;
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    // Intervals to be found for this read and counter for found intervals.
    size_t intervalsForRead = boundsForDistance.second - boundsForDistance.first;
    size_t intervalsFound = 0;

    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        SEQAN_ASSERT_EQ(static_cast<int>(value(it).distance), options.maxError);
        bool found = std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id);
        intervalsFound += found;
        if (!found && options.showMissedIntervals) {
            std::cerr << "log> {\"type\": \"log.missed_interval"
                      << "\", \"interval_id\": " << value(it).id
                      << ", \"contig_id\": \"" << value(witStore.contigNames)[value(it).contigId]
                      << "\", \"strand\": \"" << (value(it).isForward ? "forward" : "reverse")
                      << "\", \"read_id\": " << value(it).readId
                      << ", \"read_name\": \"" << value(witStore.readNames)[value(it).readId]
                      << "\", \"interval_first\": " << value(it).firstPos
                      << ", \"interval_last\": " << value(it).lastPos << "}" << std::endl;
        }
    }
    if (intervalsFound > intervalsForRead)
      exit(-1);

    // Update comparison results.
    comparisonResult.totalIntervalCount += intervalsForRead;
    comparisonResult.foundIntervalCount += intervalsFound;
    comparisonResult.totalReadCount += 1;
    comparisonResult.normalizedIntervals += 1.0 * intervalsFound / intervalsForRead;
}


template <typename TFragmentStore>
void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               TFragmentStore const & /*fragments*/,
                                               Options<EvaluateResults> const & options,
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
    if (static_cast<int>(refInterval.distance) > options.maxError)
        return;  // Guard: Skip if best has too large distance.
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    // One interval is to be found.
    comparisonResult.totalIntervalCount += 1;
    comparisonResult.totalReadCount += 1;

    // Now, try to find one of the intervals.
    bool found = false;
    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        if (std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id)) {
            found = true;
            comparisonResult.foundIntervalCount += 1;
            comparisonResult.normalizedIntervals += 1;
            break;
        }
    }
    if (!found && options.showMissedIntervals) {
      for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
          std::cerr << "log> {\"type\": \"log.missed_interval"
                    << "\", \"interval_id\": " << value(it).id
                    << ", \"contig_id\": \"" << value(witStore.contigNames)[value(it).contigId]
                    << "\", \"strand\": \"" << (value(it).isForward ? "forward" : "reverse")
                    << "\", \"read_id\": \"" << value(witStore.readNames)[value(it).readId]
                    << "\", \"distance\": " << value(beginIntervalsForRead).distance
                    << ", \"interval_first\": " << value(it).firstPos
                    << ", \"interval_last\": " << value(it).lastPos << "}" << std::endl;
      }
    }
}


template <typename TFragmentStore>
void evaluateFoundIntervals_compareToIntervals(ComparisonResult & comparisonResult,
                                               WitStore & witStore,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & beginIntervalsForRead,
                                               Iterator<WitStore::TIntervalStore, Standard>::Type & endIntervalsForRead,
                                               String<size_t> & result,
                                               TFragmentStore const & /*fragments*/,
                                               Options<EvaluateResults> const & options,
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
    refInterval.distance = value(beginIntervalsForRead).distance;  // Difference to all: Distance of smallest here not options.maxError.
    if (static_cast<int>(refInterval.distance) > options.maxError)
        return;  // Guard: Skip if best has too large distance.
    std::pair<TIntervalIterator, TIntervalIterator> boundsForDistance;
    boundsForDistance = std::equal_range(beginIntervalsForRead, endIntervalsForRead, refInterval, WitStoreLess<SortDistance>(witStore));
    if (boundsForDistance.first == boundsForDistance.second)
        return;  // Guard:  Skip if there are no required intervals.

    comparisonResult.totalReadCount += 1;
    size_t foundReads = 0;
    size_t bestReadCount = boundsForDistance.second - boundsForDistance.first;

    // Now, try to find one of the intervals.
    for (TIntervalIterator it = boundsForDistance.first; it != boundsForDistance.second; ++it) {
        comparisonResult.totalIntervalCount += 1;
        bool found = std::binary_search(begin(result, Standard()), end(result, Standard()), value(it).id);
        comparisonResult.foundIntervalCount += found;
        foundReads += found;
        if (!found && options.showMissedIntervals) {
            std::cerr << "log> {\"type\": \"log.missed_interval"
                      << "\", \"interval_id\": " << value(it).id
                      << ", \"contig_id\": \"" << value(witStore.contigNames)[value(it).contigId]
                      << "\", \"strand\": \"" << (value(it).isForward ? "forward" : "reverse")
                      << "\", \"read_id\": \"" << value(witStore.readNames)[value(it).readId]
                      << "\", \"interval_first\": " << value(it).firstPos
                      << ", \"interval_last\": " << value(it).lastPos << "}" << std::endl;
        }
    }

    comparisonResult.normalizedIntervals += 1.0 * foundReads / bestReadCount;
}


// Evaluate the ids in result pointing to intervals in witStore.
// Result is written to comparisonResult.
template <typename TFragmentStore>
void evaluateFoundIntervals(ComparisonResult & comparisonResult,
                            WitStore & witStore,
                            String<size_t> & result,
                            TFragmentStore const & fragments,
                            Options<EvaluateResults> const & options)
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
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, fragments, options, CategoryAll());
        } else if (options.benchmarkCategory == "any-best") {
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, fragments, options, CategoryAnyBest());
        } else if (options.benchmarkCategory == "all-best") {
            evaluateFoundIntervals_compareToIntervals(comparisonResult, witStore, boundsForRead.first, boundsForRead.second, result, fragments, options, CategoryAllBest());
        } else {
            SEQAN_ASSERT_FAIL("Invalid benchmark category '%s'.", toCString(options.benchmarkCategory));
        }
    }
}

// Entry point for the read mapper evaluation subprogram.
int evaluateReadMapperResult(Options<EvaluateResults> const & options)
{
    double startTime = 0;  // For measuring time below.

    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;

    std::cerr << "==============================================================================\n"
              << "                RABEMA - Read Alignment BEnchMArk\n"
              << "==============================================================================\n"
              << "                        Result Comparison\n"
              << "==============================================================================\n"
              << "\n";
    std::cerr << "____OPTIONS___________________________________________________________________\n\n";

    std::cerr << "Max error rate [%]    " << options.maxError << "\n"
              << "Oracle mode           " << (options.oracleWitMode ? (char const *) "yes" : (char const *) "no") << "\n"
              << "Benchmark category    " << options.benchmarkCategory << "\n"
              << "Distance measure      " << options.distanceFunction << "\n"
              << "Match Ns              " << (options.matchN ? (char const *) "yes" : (char const *) "no") << '\n'
              << "GSI File              " << options.witFileName << '\n'
              << "SAM File              " << options.samFileName << '\n'
              << "Reference File        " << options.seqFileName << '\n'
              << "Show\n"
              << "    additional        " << (options.showAdditionalIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "    hit               " << (options.showHitIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "    missed            " << (options.showMissedIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "    superflous        " << (options.showSuperflousIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "    try hit           " << (options.showTryHitIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "\n";

    std::cerr << "____LOADING FILES_____________________________________________________________\n\n";

    // =================================================================
    // Prepare File I/O.
    // =================================================================

    startTime = sysTime();
    // Open reference FAI index.
    std::cerr << "Reference Index           " << options.seqFileName << ".fai ...";
    FaiIndex faiIndex;
    if (load(faiIndex, toCString(options.seqFileName)) != 0)
    {
        std::cerr << " FAILED (not fatal, we can just build it)\n";
        std::cerr << "Building Index            " << options.seqFileName << ".fai ...";
        if (buildIndex(toCString(options.seqFileName), Fai()) != 0)
        {
            std::cerr << "Could not build FAI index.\n";
            return 1;
        }
        std::cerr << " OK\n";
        std::cerr << "Reference Index           " << options.seqFileName << ".fai ...";
        if (load(faiIndex, toCString(options.seqFileName)) != 0)
        {
            std::cerr << "Could not load FAI index we just build.\n";
            return 1;
        }
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }
    else
    {
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }

    std::cerr << "Reference Sequences       " << options.seqFileName << " ...";
    StringSet<Dna5String> refSeqs;
    resize(refSeqs, length(faiIndex.refNameStore));
    for (unsigned i = 0; i < length(faiIndex.refNameStore); ++i)
    {
        reserve(refSeqs[i], sequenceLength(faiIndex, i), Exact());
        if (getSequence(refSeqs[i], faiIndex, i) != 0)
        {
            std::cerr << "ERROR: Could not read sequence " << faiIndex.refNameStore[i] << ".\n";
            return 0;
        }
    }
    std::cerr << " OK\n";

    // Open gold standard intervals (GSI) file and read in header.
    std::cerr << "Gold Standard Intervals   " << options.witFileName << " (header) ...";
    std::ifstream inGsi(toCString(options.witFileName), std::ios::in | std::ios::binary);
    if (!inGsi.is_open())
    {
        std::cerr << "Could not open GSI file.\n";
        return 1;
    }
    RecordReader<std::ifstream, SinglePass<> > gsiReader(inGsi);
    WitHeader gsiHeader;
    if (readRecord(gsiHeader, gsiReader, Gsi()) != 0)
    {
        std::cerr << "Could not read GSI header.\n";
        return 1;
    }
    std::cerr << " OK\n";
    // TODO(holtgrew): Do anything with GSI header?

    // Open SAM file and read in header.
    std::cerr << "Alignments                " << options.samFileName << " (header) ...";
    std::ifstream inSam(toCString(options.samFileName), std::ios_base::in | std::ios_base::binary);
    if (!inSam.is_open())
    {
        std::cerr << "Could not open SAM file." << std::endl;
        return 1;
    }
    RecordReader<std::ifstream, SinglePass<> > samReader(inSam);
    TNameStore refNameStore;
    TNameStoreCache refNameStoreCache(refNameStore);
    BamIOContext<TNameStore> samIOContext(refNameStore, refNameStoreCache);
    BamHeader samHeader;
    if (readRecord(samHeader, samIOContext, samReader, Sam()) != 0)
    {
        std::cerr << "Could not read SAM header.\n";
        return 1;
    }
    // Check that the SAM file is sorted by query name.
    // if (getSortOrder(samHeader) != BAM_SORT_QUERYNAME)
    // {
    //     std::cerr << "SAM file not sorted by 'queryname'!\n";
    //     return 1;
    // }
    std::cerr << " OK\n";
    
    std::cerr << "\nTook " << sysTime() - startTime << "s\n";

    // =================================================================
    // Stream the SAM hits against gold standard intervals and check.
    // =================================================================

    std::cerr << "\n____COMPARING SAM HITS WITH INTERVALS_________________________________________\n\n";

    startTime = sysTime();
    typedef Position<WitStore::TIntervalStore>::Type TPos;
    // The result will be a list of ids to entries in witStore.
    int res = 0;
    RabemaStats result(options.maxError);
    if (options.distanceFunction == "edit")
        res = compareAlignedReadsToReference(result, samReader, samIOContext, refNameStore, refNameStoreCache, refSeqs, gsiReader, options, MyersUkkonenReads());
    else  // options.distanceFunction == "hamming"
        res = compareAlignedReadsToReference(result, samReader, samIOContext, refNameStore, refNameStoreCache, refSeqs, gsiReader, options, HammingSimple());
    if (res != 0)
        return 1;

    std::cerr << "\nTook " << sysTime() - startTime << " s\n";

    std::cerr << "\n____RESULTING STATISTICS______________________________________________________\n\n";

    std::cerr << "Note that in all-best and any-best mode, we differentiate the intervals we\n"
              << "found and those we have to find by their distance.  This is not possible in\n"
              << "all mode since a multiple lower-error intervals might be contained in an\n"
              << "higher-error interval.\n"
              << '\n'
              << "Alignments will be marked as \"invalid\" if they have a higher error rate than\n"
              << "allowed, i.e. they might become valid when increasing the allowed error rate.\n"
              << '\n'
              << "In \"all\" mode, only intervals from the maximal level are required.  In \"all-best\"\n"
              << "and \"any-best\" mode, the intervals are relabeled with the smallest distance that\n"
              << "a containing interval has.  Contained intervals are then removed.\n\n\n";
    
    write(std::cerr, result, options);

    /*

    typedef FragmentStore<> TFragmentStore;
    TFragmentStore fragments;

    // Load Contigs.
    double startTime = sysTime();
    std::cerr << "Reading FASTA contigs sequence file " << options.seqFileName << " ..." << std::endl;
    if (!loadContigs(fragments, options.seqFileName)) {
        std::cerr << "Could not read contigs." << std::endl;
        return 1;
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // Load Sam File.
    std::cerr << "Reading Sam file file " << options.samFileName << " ..." << std::endl;
    startTime = sysTime();
    {
        std::fstream fstrm(toCString(options.samFileName),
                           std::ios_base::in | std::ios_base::binary);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open Sam file." << std::endl;
            return 1;
        }
        read(fstrm, fragments, Sam());
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;
    //for (unsigned i = 0; i < length(fragments.readNameStore); ++i) {
      //std::cerr << ">>>" << fragments.readNameStore[i] << " " << i << std::endl;
    //}

    // =================================================================
    // Load WIT file.
    // =================================================================
    std::cerr << "Loading intervals from " << options.witFileName << std::endl;
    startTime = sysTime();
    WitStore witStore;
    loadWitFile(witStore, fragments, options.witFileName);
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Compare The Sam Hits Against WIT Intervals.
    // =================================================================
    std::cerr << "Compare reader hits from Sam file against WIT file." << std::endl;
    startTime = sysTime();
    typedef Position<WitStore::TIntervalStore>::Type TPos;
    // The result will be a list of ids to entries in witStore.
    String<size_t> result;
    if (options.distanceFunction == "edit")
        compareAlignedReadsToReference(result, fragments, witStore, options, MyersUkkonenReads());
    else  // options.distanceFunction == "hamming"
        compareAlignedReadsToReference(result, fragments, witStore, options, HammingSimple());
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Perform Counting On Result Indices, Yield ComparisonResult.
    // =================================================================
    ComparisonResult comparisonResult;
    evaluateFoundIntervals(comparisonResult, witStore, result, fragments, options);

    // =================================================================
    // Write Output.
    // =================================================================
    startTime = sysTime();
    // The output consists of one line that describes the total and
    // found intervals as a JSON record with the entries
    // "total_intervals", "found_itervals", "superflous_intervals",
    // "additional_intervals".
    if (options.outFileName == "-") {
        // Print to stdout.
        std::cout << comparisonResult << std::endl;
    } else {
        // Write output to file.
        std::fstream fstrm(toCString(options.outFileName), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open output JSON file." << std::endl;
            return 1;
        }
        fstrm << comparisonResult << std::endl;
    }
    */

    return 0;
}

#endif  // #ifndef APPS_RABEMA_EVALUATION_H_
