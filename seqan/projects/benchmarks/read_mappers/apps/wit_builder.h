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
   FUBerMarkRM -- WIT builder
  ==========================================================================*/

#include <algorithm>
#include <cmath>      // ceil(), floor()
#include <fstream>
#include <iostream>
#include <map>

#include <seqan/basic/basic_testing.h>
#include <seqan/find.h>                 // Finding infrastructure.
#include <seqan/misc/misc_cmdparser.h>  // Command parser.
#include <seqan/modifier.h>             // reverseAndComplement()
#include <seqan/store.h>                // Fragment store et al.

#include "curve_smoothing.h"
#include "find_myers_ukkonen_reads.h"

//#include "witio.h"
#include "verification.h"
#include "wit_builder_options.h"
#include "find_hamming_simple_ext.h"
#include "find_hamming_simple_quality.h"
#include "find_approx_dp_quality.h"

using namespace seqan;  // Remove some syntatic noise.

// Ceil away from Zero.
//
// ceilAwayFromZero(-1.5) == -2
// ceilAwayFromZero(1.5) == 2
template <typename T>
inline T ceilAwayFromZero(const T &x) {
    if (x < 0)
        return floor(x);
    return ceil(x);
}


// Build the error curve points around the end position of the given contig.
//
// This is pretty involved and this function easily is the most complex one.
//
// Returns rightmost border of the added points.
template <typename TContigSeq, typename TReadSeq, typename TPatternSpec, typename TReadNames>
size_t buildErrorCurvePoints(String<WeightedMatch> & errorCurve,
                             TContigSeq /*const*/ & contig,
                             size_t contigId,
                             bool isForward,
                             TReadSeq /*const*/ &read,
                             size_t readId,
                             size_t endPos,
                             int maxError,
                             size_t previousReadId,
                             size_t previousContigId,
                             size_t previousRightBorder,
                             TReadNames const & readNames,
                             bool matchN,
                             TPatternSpec const &) {
    typedef typename Position<TContigSeq>::Type TPosition;

//     std::cerr << __FILE__ << ":" << __LINE__ << " readId = " << readId << ", name = " << readNames[readId] << ", contigId = " << contigId << ", endPos = " << endPos << std::endl;
//     std::cerr << __FILE__ << ":" << __LINE__ << " previousRightBorder = " << previousRightBorder << std::endl;

    // Debug-adjustments.
    #define ENABLE 0
    #define ALL 0
    #define READID 4

//     if (readNames[readId] == CharString("SRR027007.862.1")) {
//         std::cerr << "**************** read id = " << readId << " readname = " << readNames[readId] << " read = " << read << std::endl;
//     }

    // The read maps with less than the given number of errors to [left, right]
    // to the contig sequence, is forward strand iff is Forwar is true.
    TPosition left = endPos, right = endPos;

    // Skip this alignment position if not right of previously
    // right border of the interval.
    // TODO(holtgrew): U-oh, what about a +-1 error here?
    if (readId == previousReadId && contigId == previousContigId && left <= previousRightBorder) {
//         std::cerr << __FILE__  << ":" << __LINE__ << " Skipping because " << left << " <= " << previousRightBorder << std::endl;
        return previousRightBorder;
    }

    bool ret;  // Flag used for assertions below.
    int relativeMinScore = ceilAwayFromZero(100.0 * -maxError / length(read));

    // Setup the finder and pattern.
    Finder<TContigSeq> finder(contig);
    Pattern<TReadSeq, TPatternSpec> pattern(read, -length(read) * 40);
    // If configured so, match N against all other values, otherwise match
    // against none.
    _patternMatchNOfPattern(pattern, matchN);
    _patternMatchNOfFinder(pattern, matchN);
    TPosition hitBeginPosition;

    if (ENABLE && (ALL || readId == READID)) {
        std::cerr << "**************** read id = " << readId << " readname = " << readNames[readId] << " read = " << read << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " relative min score = " << relativeMinScore << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " extending to the right." << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " contig id == " << contigId << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " is forward strand? " << isForward << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " contig length = " << length(contig) << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " endPos = " << endPos << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " infix(contig, endPos - length(read), endPos) == " << infix(contig, endPos - length(read), endPos) << std::endl;
    }

    // First, extend the interval to the right.
    {
        // Skip to original hit.
        ret = setEndPosition(finder, pattern, endPos);
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_EQ(endPos, endPosition(finder));
        ret = findBegin(finder, pattern, getScore(pattern));
        SEQAN_ASSERT_TRUE(ret);
        SEQAN_ASSERT_GT(endPos, beginPosition(finder));

        // Add original hit to the error curve points.
        int relativeScore = ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
        appendValue(errorCurve, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
        hitBeginPosition = beginPosition(finder);
        if (ENABLE && (ALL || readId == READID)) {
            std::cerr << __FILE__ << ":" << __LINE__ << " -- getScore(pattern) == " << getScore(pattern) << std::endl;
            std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(errorCurve) << " for read id " << readId << " (FIRST HIT)" << std::endl;
            std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << " endPos = " << endPos << std::endl;
        }

        // Now extend to the first hit with too low score.
        bool foundWithTooLowScore = false;
        while (find(finder, pattern)) {
            ret = findBegin(finder, pattern, getScore(pattern));
            SEQAN_ASSERT_TRUE(ret);
            if (getScore(pattern) < -maxError) {
                foundWithTooLowScore = true;
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- Found too low score." << std::endl;
                }
                break;
            }
            int relativeScore = ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
            appendValue(errorCurve, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
            if (ENABLE && (ALL || readId == READID)) {
                std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(errorCurve) << " for read id " << readId << std::endl;
                std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << std::endl;
            }
            right += 1;
        }
        // If we broke because of the score limit then collect the last not
        // yet added hit and the ones right of it until the beginPosition
        // changes.
        if (foundWithTooLowScore) {
            if (beginPosition(finder) == hitBeginPosition) {
                relativeScore = ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                appendValue(errorCurve, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
                TPosition currentBeginPosition = beginPosition(finder);
                // Add the rest until we hit one with a different begin position.
                //
                // Loop at most length(read) times (this limit is here
                // for the quality based DP algorithm which can suffer
                // from "infinite inserts").
                for (unsigned i = 0; find(finder, pattern) && i < length(read); ++i) {
                    ret = findBegin(finder, pattern, getScore(pattern));
                    SEQAN_ASSERT_TRUE(ret);
                    if (beginPosition(finder) != currentBeginPosition)
                        break;
                    relativeScore = ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                    appendValue(errorCurve, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
                    if (ENABLE && (ALL || readId == READID)) {
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(errorCurve) << " for read id " << readId << std::endl;
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << std::endl;
                    }
                    right += 1;
                }
            }
        }
        
        if (ENABLE && (ALL || readId == READID)) {
            std::cerr << __FILE__ << ":" << __LINE__ << " extending to the left." << std::endl;
        }
        // Then, extend the interval to the left.
        {
            // Length by which we extend the interval.  Note that this must be
            // at least the length of the read + max error count so no special
            // treatment of islands on the left side must be done.
            //
            // TODO(holtgrew): David said something about a parallelogram width of 16... Ask him again about this.
            TPosition kIntervalLen = length(read) + maxError;
            // Tentatively extend the interval to the left.  We will
            // extend until we hit a position with too low score.
            TPosition tentativeLeft = endPos;
            String<WeightedMatch> tempMatches;
            // Flag that indicates whether we found an entry with too low score.
            bool foundTooLowScore = false;
            // Flag that indicates whether we performed a loop iteration after having found a too low score.
            bool loopedAfterFoundTooLowScore = false;
            // Flag for breaking out of loop if we hit the right border of the last interval.
            bool hitLastRight = false;
            while (tentativeLeft > 0 && !loopedAfterFoundTooLowScore && !hitLastRight) {
                // Stop if went went to the rightmost position of the
                // last interval with the previous loop iteration.
                if (readId == previousReadId && contigId == previousContigId && tentativeLeft == previousRightBorder)
                    break;
                loopedAfterFoundTooLowScore = foundTooLowScore;
                TPosition oldTentativeLeft = tentativeLeft;
                // Move the tentative left position left by kIntervalLen but not
                // further left than 0.
                if (ENABLE && (ALL || readId == READID))
                    std::cerr << "tentativeLeft before = " << tentativeLeft << std::endl;
                tentativeLeft -= _min(kIntervalLen, tentativeLeft);
                if (ENABLE && (ALL || readId == READID))
                    std::cerr << "tentativeLeft after = " << tentativeLeft << std::endl;
                // Do not go further than the previous right position.
                if (readId == previousReadId && contigId == previousContigId && tentativeLeft <= previousRightBorder) {
                    tentativeLeft = previousRightBorder;
                    hitLastRight = true;
                }
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- tentative left = " << tentativeLeft << std::endl;
                }
                // Search from tentative left position to the previous tentative left position.
                ret = setEndPosition(finder, pattern, tentativeLeft);
                SEQAN_ASSERT_TRUE(ret);
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- endPosition(finder) = " << endPosition(finder) << std::endl;
                }
                if (endPosition(finder) > oldTentativeLeft)
                    break;  // Could not set position of the finder left of old tentative left.
                ret = findBegin(finder, pattern, getScore(pattern));
                SEQAN_ASSERT_TRUE(ret);
                int relativeScore = ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches) << " for read id " << readId << std::endl;
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << std::endl;
                }
                foundTooLowScore = foundTooLowScore || (relativeScore < relativeMinScore);
                while (find(finder, pattern) and endPosition(finder) != oldTentativeLeft) {
                    ret = findBegin(finder, pattern, getScore(pattern));
                    SEQAN_ASSERT_TRUE(ret);
                    relativeScore = ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                    appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
                    if (ENABLE && (ALL || readId == READID)) {
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches) << " for read id " << readId << std::endl;
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << std::endl;
                    }
                    foundTooLowScore = foundTooLowScore || (relativeScore < relativeMinScore);
                }
            }
            if (ENABLE && (ALL || readId == READID)) {
                std::cerr << __FILE__ << ":" << __LINE__ << " -- after loop" << std::endl;
            }
            // Now we can be sure that the temporary matches contain an entry
            // left of the first one with a too low score that has a different
            // begin position than hitBeginPosition.
            std::sort(begin(tempMatches, Standard()), end(tempMatches, Standard()));
            ModifiedString<String<WeightedMatch>, ModReverse> rTempMatches(tempMatches);
            TPosition prefixLength = 0;
            TPosition i;
            // Search up to the first too low score.
            for (i = 0; i < length(rTempMatches); ++i) {
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << "Considering " << rTempMatches[i] << std::endl;
                }
                if (rTempMatches[i].distance < relativeMinScore) {
                    if (rTempMatches[i].beginPos == hitBeginPosition) {
                        prefixLength += 1;
                        left -= 1;
                        i += 1;
                    }
                    break;
                }
                prefixLength += 1;
                left -= 1;
            }
            if (ENABLE && (ALL || readId == READID)) {
                std::cerr << "prefixLength == " << prefixLength << std::endl;
                std::cerr << "i == " << i << std::endl;
            }
            // Then, append until we find a position with a different begin position.
            for (; i < length(rTempMatches); ++i) {
                if (rTempMatches[i].beginPos == hitBeginPosition) {
                    prefixLength += 1;
                    left -= 1;
                } else {
                    break;
                }
            }
            // Finally, append the prefix of the given length.
            if (ENABLE && (ALL || readId == READID)) {
                std::cerr << "appending prefix(rTempMatches, " << prefixLength << ")" << std::endl;
            }
            append(errorCurve, prefix(rTempMatches, prefixLength));
        }
    }

//     std::cerr << __FILE__ << ":" << __LINE__ << " return " << right << std::endl;

    return right;
}



// Compute error curve from the reads and reference sequences in the
// fragment score.
template <typename TFragmentStore, typename TPatternSpec>
void matchesToErrorFunction(TFragmentStore /*const*/ & fragments,
                            TErrorCurves & errorCurves,
                            Options const & options,
                            TPatternSpec const &)
{
    typedef typename TFragmentStore::TAlignedReadStore                   TAlignedReadStore;
    typedef typename Value<TAlignedReadStore>::Type                      TAlignedRead;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type         TAlignedReadIterator;
    typedef typename TFragmentStore::TContigStore                        TContigStore;
    typedef typename TFragmentStore::TReadStore                          TReadStore;
    typedef typename TFragmentStore::TReadSeqStore                       TReadSeqStore;
    typedef typename Value<TReadStore>::Type                             TRead;
    typedef typename TRead::TId                                          TReadId;
    typedef typename Value<TContigStore>::Type                           TContig;
    typedef typename TContig::TId                                        TContigId;
    typedef typename TFragmentStore::TContigSeq                          TContigSeq;
    typedef typename TFragmentStore::TReadSeq                            TReadSeq;
    typedef typename Value<TAlignedReadStore>::Type                      TAlignedRead;
    typedef typename TAlignedRead::TPos                                  TAlignedReadPos;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;

    if (length(fragments.alignedReadStore) == 0)
        return;  // Do nothing if the aligned read store is empty.

    // TODO(holtgrew): Already sorted by contig id, iteration below much simpler!
    // Sort aligned reads by (contig index, read index, end position).
    sortAlignedReads(fragments.alignedReadStore, SortEndPos());
    sortAlignedReads(fragments.alignedReadStore, SortReadId());
    sortAlignedReads(fragments.alignedReadStore, SortContigId());

//     for (TAlignedReadIterator it = begin(fragments.alignedReadStore, Standard()); it != end(fragments.alignedReadStore, Standard()); ++it) {
//         fprintf(stderr, "%3u\t%3u\t%8lu\t%3s\n", it->contigId, it->readId, it->endPos, (it->endPos < it->beginPos ? "R" : "F"));
//     }

    // For each contig:
    //   For each read alignments:
    //     Skip if not aligned on read.
    //     If aligned on forward strand:
    //       Build error curve for this read alignment on forward strand.
    //     Else:
    //       Build error curve for this read alignment on backward strand.

    std::cerr << "Each dot represents roughly 5% of work for the given contig in the given orientation." << std::endl;
    TContigStore /*const*/ & contigStore = fragments.contigStore;
    TAlignedReadStore /*const*/ & alignedReadStore = fragments.alignedReadStore;
    for (size_t contigId = 0; contigId < length(contigStore); ++contigId) {
        TContigGaps contigGaps(contigStore[contigId].seq, contigStore[contigId].gaps);
        TContigSeq /*const*/ & contig = contigStore[contigId].seq;
        // Get reverse-complement of the contig.
        TContigSeq rcContig(contig);
        reverseComplementInPlace(rcContig);

        for (int isForward = 0; isForward <= 1; ++isForward) {
            std::cerr << "[" << fragments.contigNameStore[contigId] << "] (" << (contigId + 1) << "/" << length(contigStore) << ") " << (isForward ? "F" : "R") << " ";
            // Previous contig id and right border for both forward and backward
            // strand.  The number of contigs suffices as a sentinel value.
            TReadId previousReadId = length(fragments.readStore);
            TContigId previousContigId = length(fragments.contigStore);
            TAlignedReadPos previousRightBorder = 0;

            for (size_t i = 0; i < length(alignedReadStore); ++i) {
                if (length(alignedReadStore) / 20 > 0 && i % (length(alignedReadStore) / 20) == 0)
                    std::cerr << ".";
                size_t idx = isForward ? i : length(alignedReadStore) - i - 1;
                TAlignedRead const & alignedRead = alignedReadStore[idx];
                if (alignedRead.contigId != contigId)
                    continue;  // Skip alignments on other contig.
                if ((isForward && (alignedRead.beginPos > alignedRead.endPos)) ||
                    (!isForward && (alignedRead.beginPos < alignedRead.endPos)))
                    continue;  // Skip alignments on backwards strand when processing on forward strand.
                TReadId readId = alignedRead.readId;
                TReadSeq /*const*/ read = fragments.readSeqStore[readId];

                // Convert mapped begin and end positions from gap space
                // (with gaps in the alignment) to sequence space.
//                 std::cerr << "read name == " << fragments.readNameStore[readId] << std::endl;
//                 std::cerr << "alignedRead.contigId = " << alignedRead.contigId << std::endl;
//                 std::cerr << "alignedRead.readId = " << alignedRead.readId << std::endl;
//                 std::cerr << "alignedRead.beginPos = " << alignedRead.beginPos << std::endl;
//                 std::cerr << "alignedRead.endPos = " << alignedRead.endPos << std::endl;
                TAlignedReadPos endPos = positionGapToSeq(contigGaps, alignedRead.endPos);
//                 std::cerr << "endPos = " << endPos << std::endl;

                // Build error curve fragment around the aligned read's
                // position, depending on which strand the read was aligned.
                size_t right;
//                 std::cerr << __FILE__ << ":" << __LINE__ << " -- isForward = " << isForward << ", endPos = " << endPos << ", beginPos = " << positionGapToSeq(contigGaps, alignedRead.beginPos) << std::endl;

                // Convert from error rate from options to error count.
                int maxError = floor(0.01 * options.maxError * length(read));
                
                if (isForward) {
                    right = buildErrorCurvePoints(errorCurves[readId], contig, contigId, isForward, read, readId, endPos, maxError, previousReadId, previousContigId, previousRightBorder, fragments.readNameStore, options.matchN, TPatternSpec());
                } else {
                    right = buildErrorCurvePoints(errorCurves[readId], rcContig, contigId, isForward, read, readId, length(contig) - endPos, maxError, previousReadId, previousContigId, previousRightBorder, fragments.readNameStore, options.matchN, TPatternSpec());
                }
                previousReadId = readId;
                previousContigId = contigId;
                previousRightBorder = right;
            }
            std::cerr << std::endl;
        }
    }

    // For all reads:
    //   Sort all error curve points.
    //   Fill gaps.
    //   Smooth them.
    //   Filter out low scoring ones.

    TReadSeqStore const & readSeqs = fragments.readSeqStore;
    for (TReadId i = 0; i < length(readSeqs); ++i) {
        std::sort(begin(errorCurves[i], Standard()), end(errorCurves[i], Standard()));
        fillGaps(errorCurves[i]);
        smoothErrorCurve(errorCurves[i]);

        // Compute relative min score for the read.
        String<WeightedMatch> filtered;
        TReadSeq read = readSeqs[i];
        int maxError = floor(options.maxError / 100.0 * length(read));
        int relativeMinScore = ceilAwayFromZero(100.0 * -maxError / length(read));

        // Filter out low scoring ones.
        typedef typename Iterator<String<WeightedMatch> >::Type TIterator;
        for (TIterator it = begin(errorCurves[i]); it != end(errorCurves[i]); ++it) {
            if (value(it).distance >= relativeMinScore)
                appendValue(filtered, value(it));
        }
        move(errorCurves[i], filtered);
    }
}
