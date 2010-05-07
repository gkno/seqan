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
  ==========================================================================*/

#include <map>

#include <seqan/find.h>                 // Finding infrastructure.
#include <seqan/misc/misc_cmdparser.h>  // Command parser.
#include <seqan/store.h>                // Fragment store et al.

#include "witio.h"
#include "intervals.h"
#include "find_hamming_simple_ext.h"
#include "find_myers_ukkonen_reads.h"
#include "find_myers_ukkonen_ext.h"
#include "find_approx_dp_quality.h"
#include "find_hamming_simple_quality.h"

using namespace seqan;  // Remove some syntatic noise.


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

    // Distance function to use, also see validDistanceFunction.
    String<char> distanceFunction;

    // Name of reference sequence file.
    String<char> seqFileName;

    // Name of WIT file.
    String<char> witFileName;

    // Name of SAM file.
    String<char> samFileName;

    // Return true iff distanceFunction is a valid distance function.
    // Can be one of {"hamming", "edit"}.
    bool validDistanceFunction() const
    {
        if (distanceFunction == "hamming") return true;
        if (distanceFunction == "hamming-weighted") return true;
        if (distanceFunction == "edit") return true;
        if (distanceFunction == "edit-weighted") return true;
        return false;
    }
};

// Returns the best score for the alignment of the aligned read from
// the given fragment store.  The maximum error is given, to be able
// to limit the interval in the contig we are looking for.
template <typename TFragmentStore, typename TContigSeq2, typename TAlignedRead, typename TPatternSpec>
int bestScoreForAligned(TFragmentStore & fragments,
                        TContigSeq2 & contig2,
                        bool const & isForward,
                        TAlignedRead const & alignedRead,
                        int /*maxError*/,
                        TPatternSpec const &) {
    typedef size_t TContigId;  // TODO(holtgrew): Better type.
    typedef size_t TAlignedReadPos;  // TODO(holtgrew): Better type.
    typedef typename TFragmentStore::TContigStore      TContigStore;
    typedef typename TFragmentStore::TContigSeq        TContigSeq;
    typedef typename Value<TContigStore>::Type         TContig;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;
    typedef typename TFragmentStore::TReadSeqStore     TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq        TReadSeq;
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
    bool ret = setEndPosition(finder, pattern, endPos);
//     std::cerr << __FILE__ << ":" << __LINE__ << " -- endPos = " << endPos << ", endPosition(finder) == " << endPosition(finder) << std::endl;
    SEQAN_ASSERT_TRUE(ret);
    SEQAN_ASSERT_EQ(endPos, endPosition(finder));
    return getScore(pattern);
    /*
//     std::cout << "endPos = " << endPos << ", beginPos = " << beginPos << std::endl;
    TAlignedReadPos contigOffset;
    if (endPos < length(read) + maxError) {
        contigOffset = 0;
    } else {
        contigOffset = endPos - length(read) - maxError;
    }
    // TODO(holtgrew): We actually want to use a Suffix Segment here, but the Finder breaks then.
    TAlignedReadPos contigLength = length(contig) - contigOffset;
    typedef Segment<TContigSeq> TContigSegment;
    TContigSegment contigSegment(contig, contigOffset,
                                 contigOffset + contigLength);
    // Search the best alignment ending at the end-align position.
//     std::cout << "Aligning read " << read << std::endl;
//     std::cout << "Aligning contig " << prefix(contigSegment, length(read) + maxError) << std::endl;
    Finder<TContigSegment> finder(contigSegment);
    const int kErrorValueInfty = 100000; // TODO(holtgrew): Should depend on scoring.
    Pattern<TReadSeq, TPatternSpec> pattern(read, -kErrorValueInfty);
    while (find(finder, pattern) and endPosition(finder) < length(read) + maxError)
        continue;  // Skip.
//     std::cerr << "end position: " << endPosition(finder) << std::endl;
//     std::cerr << "length(read) == " << length(read) << std::endl;
//     std::cerr << "length(contig) == " << length(contig) << std::endl;
//     std::cerr << "max error = " << maxError << std::endl;
//     std::cerr << "end position = " << endPosition(finder) << std::endl;
//     std::cerr << "contig offset = " << contigOffset << std::endl;
//     std::cerr << "sum offset = " << endPosition(finder) + contigOffset << std::endl;
//     std::cerr << "read is " << read << std::endl;
//     std::cerr << "contig segment[:length(read) + maxError + 2] is " << prefix(contigSegment, length(read) + maxError + 2) << std::endl;
//     std::cerr << "contig is " << infix(fragments.contigStore[alignedRead.contigId].seq, endPos - length(readSeqs[alignedRead.readId]), endPos) << std::endl;
//     std::cerr << "          " << fragments.contigNameStore[alignedRead.contigId] << std::endl;
    SEQAN_ASSERT_EQ(endPosition(finder), length(read) + maxError);

    return getScore(pattern);
    */
}


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


// Return maximum error count for maximum error rate
inline int maxErrorRateToMaxErrors(int maxErrorRate, size_t len) {
    return floor(maxErrorRate / 100.0 * len);
}


template <typename TFragmentStore, typename TAlignedReadIter, typename TWitRecordIter, typename TContigSeq, typename TContigGaps, typename TPatternSpec>
void
compareAlignedReadsToReferenceOnContigForOneRead(Options const & options,
                                                 TFragmentStore & fragments,
                                                 size_t const & contigId,
                                                 TContigSeq const & contig,
                                                 TContigGaps const & contigGaps,
                                                 bool const & isForward,
                                                 size_t const & readId,
                                                 TAlignedReadIter const & alignedReadsBegin,
                                                 TAlignedReadIter const & alignedReadsEnd,
                                                 TWitRecordIter const & witRecordsBegin,
                                                 TWitRecordIter const & witRecordsEnd,
                                                 size_t & foundIntervalCount,
                                                 size_t & relevantIntervalCount,
                                                 TPatternSpec const &) {
    // TODO(holtgrew): Enable n-matches-wildcard mode and use qualities for weights.
    typedef size_t TPos;
    typedef std::map<size_t, FlaggedInterval> TMap;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename TMap::iterator TMapIterator;

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
        if (it->distance != options.maxError) continue;  // Skip intervals with wrong distance.

//         std::cerr << __FILE__ << ":" << __LINE__ << " -- add(" << it->firstPos << ", " << it->lastPos << ", is forward = " << it->isForward << ")" << "  " << *it << std::endl;
        intervalMap[it->lastPos] = FlaggedInterval(it->firstPos, it->lastPos);
        relevantIntervalCount += 1;
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
        int bestScore = bestScoreForAligned(fragments, contig, isForward, *it, maxErrorRateToMaxErrors(options.maxError, length(fragments.readSeqStore[it->readId])), TPatternSpec());
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
            continue;
        }

        // Handle alignment out of target intervals.
        if (iter == intervalMap.end() || (iter->second.firstPos > lastPos) || (iter->second.lastPos < lastPos)) {
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
            exit(1);
        }

        SEQAN_ASSERT_LEQ(iter->second.firstPos, lastPos);
        SEQAN_ASSERT_GEQ(iter->second.lastPos, lastPos);

        // Count interval if hit for the this time.
        if (!iter->second.flag) {
            iter->second.flag = true;
            foundIntervalCount += 1;
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
                                       TContigSeq const & contig,
                                       bool const & isForward,
                                       TAlignedReadsIter const & alignedReadsBegin,
                                       TAlignedReadsIter const & alignedReadsEnd,
                                       TWitRecordsIter const & witRecordsBegin,
                                       TWitRecordsIter const & witRecordsEnd,
                                       size_t & foundIntervalCount,
                                       size_t & relevantIntervalCount,
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
        compareAlignedReadsToReferenceOnContigForOneRead(options, fragments, contigId, contig, contigGaps, isForward, readId, alignedReadsBeginForRead, alignedReadsEndForRead, witRecordsBeginForRead, witRecordsEndForRead, foundIntervalCount, relevantIntervalCount, TPatternSpec());
        // This iteration's end iterators are the next iteration's begin iterators.
        alignedReadsBeginForRead = alignedReadsEndForRead;
        witRecordsBeginForRead = witRecordsEndForRead;
    }
}


// Compare aligned reads in fragment store to the intervals specified
// in witRecords, the number of hits is written to foundIntervalCount,
// the number of relevant intervals is written to
// relevantIntervalCount.
template <typename TFragmentStore, typename TPatternSpec>
void
compareAlignedReadsToReference(Options const & options,
                               TFragmentStore & fragments,  // non-const so reads can be sorted
                               String<WitRecord> & witRecords,  // non-const so it is sortable
                               size_t & foundIntervalCount,
                               size_t & relevantIntervalCount,
                               TPatternSpec const &) {
    // Type aliases.
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore>::Type TAlignedReadsIter;
    typedef typename Iterator<String<WitRecord> >::Type TWitRecordsIter;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef typename TContigStoreElement::TContigSeq TContigSeq;
    
    // Initialization.
    foundIntervalCount = 0;
    relevantIntervalCount = 0;

    // Sort aligned reads and wit records by contig id.
    sortAlignedReads(fragments.alignedReadStore, SortEndPos());
    sortAlignedReads(fragments.alignedReadStore, SortReadId());
    sortAlignedReads(fragments.alignedReadStore, SortContigId());
    std::sort(begin(witRecords, Standard()), end(witRecords, Standard()), WitRecord_Lt_ContigIdReadIdLastPos());

    TAlignedReadsIter alignedReadsBegin = begin(fragments.alignedReadStore, Standard());
    TAlignedReadsIter alignedReadsEnd = alignedReadsBegin;
//     std::cerr << "# of alignd reads " << length(fragments.alignedReadStore) << std::endl;
    TWitRecordsIter witRecordsBegin = begin(witRecords, Standard());
    TWitRecordsIter witRecordsEnd = witRecordsBegin;
//     std::cerr << "# of wit records " << length(witRecords) << std::endl;

//     std::cerr << ".--- wit records" << std::endl;
//     for (TWitRecordsIter it = begin(witRecords, Standard()); it != end(witRecords, Standard()); ++it)
//         std::cerr << "| " << *it << " --- " << it->contigId << std::endl;
//     std::cerr << "`---" << std::endl;

    while (alignedReadsBegin != end(fragments.alignedReadStore, Standard()) &&
           witRecordsBegin != end(witRecords, Standard())) {
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
        while (witRecordsEnd != end(witRecords, Standard()) && witRecordsEnd->contigId <= contigId) {
//             if (witRecordsEnd->distance == options.maxError)
//                 std::cerr << "| " << *witRecordsEnd << std::endl;
            ++witRecordsEnd;
        }
//         std::cerr << "`---" << std::endl;

        // Actually compare the aligned reads for this contig on forward and backwards strand.
        TContigSeq const & contig = fragments.contigStore[contigId].seq;
        compareAlignedReadsToReferenceOnContig(options, fragments, contigId, contig, true, alignedReadsBegin, alignedReadsEnd, witRecordsBegin, witRecordsEnd, foundIntervalCount, relevantIntervalCount, TPatternSpec());
        TContigSeq rcContig(contig);
        reverseComplementInPlace(rcContig);
        compareAlignedReadsToReferenceOnContig(options, fragments, contigId, rcContig, false, alignedReadsBegin, alignedReadsEnd, witRecordsBegin, witRecordsEnd, foundIntervalCount, relevantIntervalCount, TPatternSpec());
        
        // This iteration's end iterators are the next iteration's begin iterators.
        alignedReadsBegin = alignedReadsEnd;
        witRecordsBegin = witRecordsEnd;
    }
}
