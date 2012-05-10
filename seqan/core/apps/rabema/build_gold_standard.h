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

#ifndef APPS_RABEMA_BUILD_GOLD_STANDARD_H_
#define APPS_RABEMA_BUILD_GOLD_STANDARD_H_

// TODO(holtgrew): Another possible way to save memory is not to store the first, always identical part of the read id.
// TODO(holtgrew): Yet another way to save memory is to load the contigs on demand, i.e. using FAI.

#include <fstream>
#include <iostream>
#include <map>

#include <seqan/store.h>
#include <seqan/bam_io.h>

#include "rabema.h"
#include "build_gold_standard_options.h"

#include "find_myers_ukkonen_reads.h"

#include "verification.h"
#include "find_hamming_simple_ext.h"
#include "find_hamming_simple_quality.h"
#include "find_approx_dp_quality.h"
#include "intervals.h"
#include "verification.h"

#include "fai_index.h"

// ============================================================================
// Enums, Tags, Classes, Typedefs.
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

void trimSeqHeaderToId(seqan::CharString & header)
{
    unsigned i = 0;
    for (; i < length(header); ++i)
        if (isspace(header[i]))
          break;
    resize(header, i);
}

// Build intervals from the error curves.
void intervalizeErrorCurves(String<WitRecord> & result,
                            TErrorCurves const & errorCurves,
                            StringSet<CharString> const & readNameStore,
                            StringSet<CharString> const & contigNameStore,
                            Options<BuildGoldStandard> const & options) {

    std::cerr << "\n____POINT TO INTERVAL CONVERSION______________________________________________\n\n"
              << "Progress: ";
    unsigned tenPercent = errorCurves.size() / 10;
    unsigned progressI = 0;
    typedef TErrorCurves::const_iterator TErrorCurvesIter;
    for (TErrorCurvesIter it = errorCurves.begin(); it != errorCurves.end(); ++it, ++progressI)
    {
        if (tenPercent > 0u && progressI % tenPercent == 0u)
            std::cerr << progressI / tenPercent * 10 << '%';
        else if (tenPercent > 5u && progressI % (tenPercent / 5) == 0u)
            std::cerr << '.';

        size_t readId = it->first;
        TWeightedMatches const & matches = it->second;

        // Sort the matches.  Matches with high scores (negative score and
        // low absolute value) come first.
        TWeightedMatches sortedMatches(matches);
        std::sort(begin(sortedMatches, Standard()), end(sortedMatches, Standard()));
        //std::cerr << "\n -- " << readNameStore[readId] << "\n";
        //for (unsigned i = 0; i < length(sortedMatches); ++i)
        //    std::cerr << sortedMatches[i] << "\n";
        //std::cerr << "\n";

        // intervals[e] holds the intervals for error e of the current read.
        String<String<ContigInterval> > intervals;
        int maxError = options.oracleSamMode ? 0 : (int)options.maxError;
        resize(intervals, maxError + 1);

        // Join the intervals stored in sortedMatches.
        //
        // The position of the previous match, so we can consider only the
        // ones with the smallest error.
        //
        // The following two vars should be != first pos and contigId.
        size_t previousPos = maxValue<size_t>();
        size_t previousContigId = maxValue<size_t>();
        typedef Iterator<TWeightedMatches>::Type TWeightedMatchesIter;
        for (TWeightedMatchesIter it = begin(sortedMatches);
             it != end(sortedMatches); ++it) {
            // Skip it if (it - 1) pointed to same pos (and must point to
            // one with smaller absolute distance.
            if (it->pos == previousPos && it->contigId == previousContigId)
                continue;
            // Consider all currently open intervals with a greater error
            // than the error in *it and extend them or create a new one.
            int error = options.oracleSamMode ? 0 : abs(it->distance);
            SEQAN_ASSERT_LEQ(error, maxError);
            for (int e = error; e <= maxError; ++e) {
                // Handle base case of no open interval:  Create new one.
                if (length(intervals[e]) == 0) {
                    appendValue(intervals[e], ContigInterval(it->contigId, it->isForward, it->pos, it->pos));
                    continue;
                }
                ContigInterval &interval = back(intervals[e]);
                // Either extend the interval or create a new one.
                if (interval.contigId == it->contigId && interval.isForward == it->isForward)
                    SEQAN_ASSERT_LEQ(interval.last, it->pos);
                if (interval.contigId == it->contigId && interval.isForward == it->isForward && interval.last + 1 == it->pos)
                    back(intervals[e]).last += 1;
                else
                    appendValue(intervals[e], ContigInterval(it->contigId, it->isForward, it->pos, it->pos));
            }
            // Book-keeping.
            previousPos = it->pos;
            previousContigId = it->contigId;
        }

        // Print the resulting intervals.
        typedef Iterator<String<String<ContigInterval> > >::Type TIntervalContainerIter;
        int distance = 0;
        for (TIntervalContainerIter it = begin(intervals);
             it != end(intervals); ++it, ++distance) {
            typedef Iterator<String<ContigInterval> >::Type TIntervalIter;
            for (TIntervalIter it2 = begin(*it); it2 != end(*it); ++it2) {
                int flags = 0;
                // We appended custom prefixes to the read ids.  "/S" means single-end, "/0" means left mate, "/1" means
                // right mate.
                int mateNo = -1;
                if (back(readNameStore[readId]) == '0')
                    mateNo = 0;
                if (back(readNameStore[readId]) == '1')
                    mateNo = 1;
                SEQAN_ASSERT_EQ(readNameStore[readId][length(readNameStore[readId]) - 2], '/');
                if (mateNo == 0)
                  flags = WitRecord::FLAG_PAIRED | WitRecord::FLAG_FIRST_MATE;
                else if (mateNo == 1)
                  flags = WitRecord::FLAG_PAIRED | WitRecord::FLAG_SECOND_MATE;

                appendValue(result,
                            WitRecord(prefix(readNameStore[readId], length(readNameStore[readId]) - 2), flags,
                                      distance, contigNameStore[it2->contigId],
                                      it2->isForward, it2->first, it2->last));
            }
        }
    }
    std::cerr << "100% DONE\n";
}


// Build the error curve points around the end position of the given contig.
//
// This is pretty involved and this function easily is the most complex one.
//
// Returns rightmost border of the added points.
template <typename TContigSeq, typename TReadSeq, typename TPatternSpec, typename TReadNames>
size_t buildErrorCurvePoints(String<WeightedMatch> & errorCurve,
                             int & maxError,
                             TContigSeq /*const*/ & contig,
                             size_t contigId,
                             bool isForward,
                             TReadSeq /*const*/ &read,
                             size_t readId,
                             size_t endPos,
                             size_t previousReadId,
                             size_t previousContigId,
                             size_t previousRightBorder,
                             TReadNames const & readNames,
                             bool matchN,
                             TPatternSpec const &) {
    typedef typename Position<TContigSeq>::Type TPosition;

//     std::cerr << __FILE__ << ":" << __LINE__ << " readId = " << readId << ", name = " << readNames[readId] << ", contigId = " << contigId << ", endPos = " << endPos << std::endl;
//     std::cerr << __FILE__ << ":" << __LINE__ << " previousRightBorder = " << previousRightBorder << std::endl;

    // In oracle Sam mode, the maximum error is the error at the position given in the Sam alignment.
    bool oracleSamMode = false;
    if (maxError == -1) {
        oracleSamMode = true;
        Finder<TContigSeq> finder(contig);
        Pattern<TReadSeq, TPatternSpec> pattern(read, -(int)length(read) * 40);
        bool ret = setEndPosition(finder, pattern, endPos);
        (void) ret; // If compiled without assertions.
        SEQAN_ASSERT(ret);
        maxError = -getScore(pattern);
    }
    
    // Debug-adjustments.
    #define ENABLE 0
    #define ALL 0
    #define READID 2

//     if (readNames[readId] == CharString("SRR027007.862.1")) {
//           std::cerr << "**************** read id = " << readId << " readname = " << readNames[readId] << " read = " << read << std::endl;
//     }

    // The read maps with less than the given number of errors to [left, right]
    // to the contig sequence, is forward strand iff is Forwar is true.
    TPosition /*left = endPos,*/ right = endPos;

    // Skip this alignment position if not right of previously
    // right border of the interval.
    // TODO(holtgrew): Remove this piece, it's unused right now. We could/should also fix this to work for streaming over SAM files.
    // TODO(holtgrew): U-oh, what about a +-1 error here?
    //if (readId == previousReadId && contigId == previousContigId && left <= previousRightBorder) {
    //    //std::cerr << __FILE__  << ":" << __LINE__ << " Skipping because " << left << " <= " << previousRightBorder << std::endl;
    //    return previousRightBorder;
    //}

    bool ret;  // Flag used for assertions below.
    (void) ret;  // If run without assertions.
    int relativeMinScore = (int)ceilAwayFromZero(100.0 * -maxError / length(read));

    // Setup the finder and pattern.
    Finder<TContigSeq> finder(contig);
    Pattern<TReadSeq, TPatternSpec> pattern(read, -(int)length(read) * 40);
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

    // We will first gather all results in tempMatches.  Below, we
    // will smooth the curve in these points and throw out too bad
    // entries.  Then, we will append tempMatches to errorCurve.
    String<WeightedMatch> tempMatches;

    // First, extend the interval to the right.
    {
        // Skip to original hit.
        ret = setEndPosition(finder, pattern, endPos);
        SEQAN_ASSERT(ret);
        SEQAN_ASSERT_EQ(endPos, endPosition(finder));
        SEQAN_ASSERT_GEQ(getScore(pattern), -maxError);
        while (findBegin(finder, pattern, getScore(pattern)))
            continue;  // Find leftmost begin position.
        SEQAN_ASSERT_GT(endPos, beginPosition(finder));
        SEQAN_ASSERT_EQ(getScore(pattern), getBeginScore(pattern));

        // Add original hit to the error curve points.
        int relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
        appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
        hitBeginPosition = beginPosition(finder);
        if (ENABLE && (ALL || readId == READID)) {
            std::cerr << __FILE__ << ":" << __LINE__ << " -- getScore(pattern) == " << getScore(pattern) << std::endl;
            std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches) << " for read id " << readId << " (FIRST HIT)" << std::endl;
            std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << " endPos = " << endPos << std::endl;
        }

        // Now extend to the first hit with too low score.
        bool foundWithTooLowScore = false;
        while (find(finder, pattern)) {
	        while (findBegin(finder, pattern, getScore(pattern)))
    	        continue;  // Find leftmost begin position.
            if (getScore(pattern) < -maxError && beginPosition(finder) != back(tempMatches).beginPos) {
                foundWithTooLowScore = true;
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- Found too low score." << std::endl;
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- low scoring match was " << WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)) << std::endl;
                }
                break;
            }
            int relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
            appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
            if (ENABLE && (ALL || readId == READID)) {
                std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches) << " for read id " << readId << std::endl;
                std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << std::endl;
            }
            right += 1;
        }
        // If we broke because of the score limit then collect the last not
        // yet added hit and the ones right of it until the beginPosition
        // changes.
        if (foundWithTooLowScore) {
            if (beginPosition(finder) == hitBeginPosition) {
                relativeScore = static_cast<int>(ceilAwayFromZero(100.0 * static_cast<double>(getScore(pattern)) / static_cast<double>(length(read))));
                appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
                TPosition currentBeginPosition = beginPosition(finder);
                // Add the rest until we hit one with a different begin position.
                //
                // Loop at most length(read) times (this limit is here
                // for the quality based DP algorithm which can suffer
                // from "infinite inserts").
                for (unsigned i = 0; find(finder, pattern) && i < length(read); ++i) {
			        while (findBegin(finder, pattern, getScore(pattern)))
			            continue;  // Find leftmost begin position.
                    SEQAN_ASSERT(ret);
                    if (beginPosition(finder) != currentBeginPosition)
                        break;
                    relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                    appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
                    if (ENABLE && (ALL || readId == READID)) {
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches) << " for read id " << readId << std::endl;
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
            // Flag that indicates whether we found an entry with too low score.
            bool foundTooLowScore = false;
            // Flag that indicates whether we performed a loop iteration after having found a too low score.
            bool loopedAfterFoundTooLowScore = false;
            // Flag for breaking out of loop if we hit the right border of the last interval.
            bool hitLastRight = false;
            while (tentativeLeft > 0 && !loopedAfterFoundTooLowScore && !hitLastRight) {
                // Stop if went went to the rightmost position of the
                // last interval with the previous loop iteration.
                //if (readId == previousReadId && contigId == previousContigId && tentativeLeft == previousRightBorder)
                //    break;
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
                SEQAN_ASSERT(ret);
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- endPosition(finder) = " << endPosition(finder) << std::endl;
                }
                if (endPosition(finder) > oldTentativeLeft)
                    break;  // Could not set position of the finder left of old tentative left.
			    while (findBegin(finder, pattern, getScore(pattern)))
			        continue;  // Find leftmost begin position.
                int relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches) << " for read id " << readId << std::endl;
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << std::endl;
                }
                foundTooLowScore = foundTooLowScore || (relativeScore < relativeMinScore);
                while (find(finder, pattern) && endPosition(finder) != oldTentativeLeft) {
			        while (findBegin(finder, pattern, getScore(pattern)))
            			continue;  // Find leftmost begin position.
                    SEQAN_ASSERT(ret);
                    relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                    appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore, beginPosition(finder)));
                    if (ENABLE && (ALL || readId == READID)) {
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- raw score is " << getScore(pattern) << std::endl;
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches) << " for read id " << readId << std::endl;
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read << std::endl;
                    }
                    foundTooLowScore = foundTooLowScore || (relativeScore < relativeMinScore);
                }
            }
            if (ENABLE && (ALL || readId == READID)) {
                std::cerr << __FILE__ << ":" << __LINE__ << " -- after loop" << std::endl;
            }
            /*
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
        */
        }
    }
    // Postprocessing: Sorting, smoothing, filtering.
    std::sort(begin(tempMatches, Standard()), end(tempMatches, Standard()));
    smoothErrorCurve(tempMatches);
    appendValue(tempMatches, WeightedMatch(0, 0, 0, relativeMinScore - 1, 0));  // Sentinel.
    if (oracleSamMode) {
        // In oracle Sam mode, we only want the lake with last pos endPos-1.
        String<WeightedMatch> buffer;
        bool flag = false;
        for (size_t i = 0; i < length(tempMatches); ++i) {
            if (tempMatches[i].distance < relativeMinScore) {
                if (flag) {
                    append(errorCurve, buffer);
                    break;
                } else {
                    clear(buffer);
                }
            } else {
                appendValue(buffer, tempMatches[i]);
                if (tempMatches[i].pos == endPos-1)
                    flag = true;
            }
        }
    } else {
        for (size_t i = 0; i < length(tempMatches); ++i) {
            if (tempMatches[i].distance >= relativeMinScore)
                appendValue(errorCurve, tempMatches[i]);
        }
    }

//     if (readId == READID) {
//         std::cerr << ",-- errorCurve is (read id = " << readId << " readname = " << readNames[readId] << ")" << std::endl;
//         for (unsigned i = 0; i < length(errorCurve); ++i) {
//             std::cerr << "| " << errorCurve[i] << std::endl;
//         }
//         std::cerr << "`--" << std::endl;
//     }
    
//     std::cerr << __FILE__ << ":" << __LINE__ << " return " << right << std::endl;

    return right;
}



// Compute error curve from the reads and reference sequences in the
// fragment score.
template <typename TPatternSpec>
int matchesToErrorFunction(TErrorCurves & errorCurves,
                           RecordReader<std::ifstream, SinglePass<> > & samReader,
                           BamIOContext<StringSet<CharString> > & samIOContext,
                           StringSet<CharString> & readNameStore,
                           StringSet<CharString> const & contigNameStore,
                           FaiIndex const & faiIndex,
                           Options<BuildGoldStandard> const & options,
                           TPatternSpec const &)
{
    // typedef typename TFragmentStore::TAlignedReadStore                   TAlignedReadStore;
    // typedef typename Value<TAlignedReadStore>::Type                      TAlignedRead;
    // typedef typename Iterator<TAlignedReadStore, Standard>::Type         TAlignedReadIterator;
    // typedef typename TFragmentStore::TContigStore                        TContigStore;
    // typedef typename TFragmentStore::TReadStore                          TReadStore;
    // typedef typename TFragmentStore::TReadSeqStore                       TReadSeqStore;
    // typedef typename Value<TReadStore>::Type                             TRead;
    // typedef typename TRead::TId                                          TReadId;
    // typedef typename Value<TContigStore>::Type                           TContig;
    // typedef typename TContig::TId                                        TContigId;
    // typedef typename TFragmentStore::TContigSeq                          TContigSeq;
    // typedef typename TFragmentStore::TReadSeq                            TReadSeq;
    // typedef typename Value<TAlignedReadStore>::Type                      TAlignedRead;
    // typedef typename TAlignedRead::TPos                                  TAlignedReadPos;
    // typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;

    double startTime = 0;

    if (atEnd(samReader))
        return 0;  // Do nothing if there are no more alignments in the file.

    // TODO(holtgrew): Can we circumvent using a store here? Requires lots of memory.

    // We store the read names and append "/0" and "/1", depending on their flag value ("/0" for first, "/1" for last
    // flag).  Sadly, there is no easy way out here.  Single-end reads are stored as "/S".
    NameStoreCache<StringSet<CharString> > readNameStoreCache(readNameStore);
    String<unsigned> readLengthStore;

    // In oracle SAM mode, we store the distance of the alignment from the SAM file for each read.  Otherwise, this
    // variable remains unused.
    String<int> readAlignmentDistances;

//     for (TAlignedReadIterator it = begin(fragments.alignedReadStore, Standard()); it != end(fragments.alignedReadStore, Standard()); ++it) {
//         fprintf(stderr, "%3u\t%3u\t%8lu\t%3s\n", it->contigId, it->readId, it->endPos, (it->endPos < it->beginPos ? "R" : "F"));
//     }

    std::cerr << "\n____EXTENDING INTERVALS_______________________________________________________\n\n"
              << "Each dot represents work for 100k nucleotides in the genome.\n";

    // Stream over SAM file.  The main assumption here is that the reads are sorted by coordinate, which is check when
    // loading the SAM header.
    //
    // For each read alignments:
    //   If aligned on forward strand:
    //     Build error curve for this read alignment on forward strand.
    //   Else:
    //     Build error curve for this read alignment on backward strand.

    startTime = sysTime();      // Time at beginning for total time display at end.
    int prevRefId = -1;      // Previous contig id.
    unsigned prevReadId = maxValue<unsigned>();
    size_t prevRightBorder = maxValue<size_t>();
    int posOnContig = 0;
    BamAlignmentRecord record;  // Current read record.
    Dna5String contig;
    Dna5String rcContig;
    Dna5String readSeq;
    CharString readName;
    while (!atEnd(samReader))
    {
        // -------------------------------------------------------------------
        // Read SAM Record and retrieve name and sequence.
        // -------------------------------------------------------------------

        // Read SAM record.
        if (readRecord(record, samIOContext, samReader, Sam()) != 0)
        {
            std::cerr << "ERROR reading SAM record!\n";
            return 1;
        }
        // Break if we have an unaligned SAM record.
        if (record.rId == -1)
            break;  // Done!
        // Get read name and sequence from record.
        readSeq = record.seq;  // Convert read sequence to Dna5.
        // Compute reverse complement since we align against reverse strand, SAM has aligned sequence against forward
        // strand.
        if (hasFlagRC(record))
            reverseComplement(readSeq);
        trimSeqHeaderToId(record.qName);  // Remove everything after the first whitespace.
        if (!hasFlagMultiple(record))
            append(record.qName, "/S");
        if (hasFlagMultiple(record) && hasFlagFirst(record))
            append(record.qName, "/0");
        if (hasFlagMultiple(record) && hasFlagLast(record))
            append(record.qName, "/1");
        // Translate read to read id.
        // TODO(holtgrew): Can we get rid of this translation to numeric ids?
        unsigned readId = 0;
        if (!getIdByName(readNameStore, record.qName, readId, readNameStoreCache))
        {
            readId = length(readNameStore);
            appendName(readNameStore, record.qName, readNameStoreCache);
            appendValue(readLengthStore, length(record.seq));
            if (options.oracleSamMode)
                appendValue(readAlignmentDistances, -1);
        }

        // -------------------------------------------------------------------
        // Handle Progress
        // -------------------------------------------------------------------

        // Handle progress: Change of contig.
        SEQAN_ASSERT_LEQ(prevRefId, record.rId);
        if (prevRefId != record.rId)
        {
            for (int i = prevRefId + 1; i <= record.rId; ++i)
            {
                if (i != prevRefId)
                    std::cerr << "\n";
                std::cerr << contigNameStore[record.rId] << " (" << record.rId + 1 << "/" << length(contigNameStore) << ") ";
            }
            posOnContig = 0;

            // Load reference sequence on the fly using FAI index to save memory.
            unsigned faiRefId = 0;
            if (!getIdByName(faiIndex, contigNameStore[record.rId], faiRefId))
            {
                std::cerr << "Reference sequence " << contigNameStore[record.rId] << " not known in FAI file.\n";
                return 1;
            }
            if (getSequence(contig, faiIndex, faiRefId) != 0)
            {
                std::cerr << "Could not load sequence " << contigNameStore[record.rId] << " from FASTA file with FAI.\n";
                return 1;
            }
            // std::cerr << "PREFIX " << prefix(contig, 100) << "\n";
            // std::cerr << "SUFFIX " << suffix(contig, length(contig) - 100) << "\n";
            rcContig = contig;
            reverseComplement(rcContig);
        }
        // Handle progress: Position on contig.
        for (; posOnContig < record.pos; posOnContig += 100*1000)
        {
            if (posOnContig % (1000*1000) == 0 && posOnContig > 0)
                std::cerr << posOnContig / 1000 / 1000 << "M";
            else
                std::cerr << ".";
        }

        // -------------------------------------------------------------------
        // Extend error curve points for read.
        // -------------------------------------------------------------------

        // In oracle SAM mode, set max error to -1, buildErrorCurvePoints() will use the error at the alignment position
        // from the SAM file.  In normal mode, convert from error rate from options to error count.
        int maxError = options.oracleSamMode ? -1 : static_cast<int>(floor(0.01 * options.maxError * length(record.seq)));

        // Compute end position of alignment.
        int endPos = record.pos + getAlignmentLengthInRef(record) - countPaddings(record.cigar);

        size_t right = 0;
        if (!hasFlagRC(record))
            right = buildErrorCurvePoints(errorCurves[readId], maxError, contig, record.rId, !hasFlagRC(record), readSeq, readId, endPos, prevReadId, prevRefId, prevRightBorder, readNameStore, options.matchN, TPatternSpec());
        else
            right = buildErrorCurvePoints(errorCurves[readId], maxError, rcContig, record.rId, !hasFlagRC(record), readSeq, readId, length(rcContig) - record.pos, prevReadId, prevRefId, prevRightBorder, readNameStore, options.matchN, TPatternSpec());
        if (options.oracleSamMode)
            readAlignmentDistances[readId] = maxError;

        // Update variables storing the previous read/contig id and position.
        prevReadId = readId;
        prevRefId = record.rId;
        prevRightBorder = right;
    }
    std::cerr << "\n\nTook " << sysTime() - startTime << " s\n";

#if 0
    TContigStore /*const*/ & contigStore = fragments.contigStore;
    TAlignedReadStore /*const*/ & alignedReadStore = fragments.alignedReadStore;
    for (size_t contigId = 0; contigId < length(contigStore); ++contigId) {
        TContigGaps contigGaps(contigStore[contigId].seq, contigStore[contigId].gaps);
        TContigSeq /*const*/ & contig = contigStore[contigId].seq;
        // Get reverse-complement of the contig.
        TContigSeq rcContig(contig);
        reverseComplement(rcContig);

        for (int isForward = 0; isForward <= 1; ++isForward) {
            std::cerr << "[" << fragments.contigNameStore[contigId] << "] (" << (contigId + 1) << "/" << length(contigStore) << ") " << (isForward ? "F" : "R") << " ";

            // Sort aligned reads by (contig index, read index, end
            // position) when iterating forward and by begin position
            // instead when iterating backwards.
            //
            // TODO(holtgrew): Only re-sort with a given contig id? We are sorting $contig-count times too often.
            if (isForward) {
                sortAlignedReads(fragments.alignedReadStore, SortEndPos());
                sortAlignedReads(fragments.alignedReadStore, SortReadId());
                sortAlignedReads(fragments.alignedReadStore, SortContigId());
            } else {
                sortAlignedReads(fragments.alignedReadStore, SortBeginPos());
                sortAlignedReads(fragments.alignedReadStore, SortReadId());
                sortAlignedReads(fragments.alignedReadStore, SortContigId());
            }
            
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
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- beginPos = " << positionGapToSeq(contigGaps, alignedRead.beginPos) << ", endPos = " << endPos << std::endl;
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- rev beginPos = " << length(contig) - positionGapToSeq(contigGaps, alignedRead.beginPos) << ", rev endPos = " << length(contig) - endPos << std::endl;
                }

                // Build error curve fragment around the aligned read's
                // position, depending on which strand the read was aligned.
                size_t right;
//                std::cerr << __FILE__ << ":" << __LINE__ << " -- read name == " << fragments.readNameStore[readId] << ", isForward = " << isForward << ", endPos = " << endPos << ", beginPos = " << positionGapToSeq(contigGaps, alignedRead.beginPos) << std::endl;

                int maxError;
                if (options.oracleSamMode) {
                    // In oracle Sam mode, set max error to -1, buildErrorCurvePoints() will use the error at the alignment position from the Sam file.
                    maxError = -1;
                } else {
                    // In normal mode, convert from error rate from options to error count.
                    maxError = (int)floor(0.01 * options.maxError * length(read));
                }
                
                if (isForward) {
                    right = buildErrorCurvePoints(errorCurves[readId], maxError, contig, contigId, isForward, read, readId, endPos, previousReadId, previousContigId, previousRightBorder, fragments.readNameStore, options.matchN, TPatternSpec());
                } else {
                    right = buildErrorCurvePoints(errorCurves[readId], maxError, rcContig, contigId, isForward, read, readId, length(contig) - endPos, previousReadId, previousContigId, previousRightBorder, fragments.readNameStore, options.matchN, TPatternSpec());
                }
                if (options.oracleSamMode)
                    readAlignmentDistances[readId] = maxError;
                previousReadId = readId;
                previousContigId = contigId;
                previousRightBorder = right;
            }
            std::cerr << std::endl;
        }
    }
    if (options.verbosity >= 2)
        std::cerr << "[timer] sorting/gap filling/smoothing/filtering: " << sysTime() - startTime << " s" << std::endl;
#endif  // #if 0

    // For all reads:
    //   Sort all error curve points.
    //   Fill gaps.
    //   Smooth them.
    //   Filter out low scoring ones.

    std::cerr << "\n____SMOOTHING ERROR CURVES____________________________________________________\n\n";
    startTime = sysTime();
    unsigned tenPercent = length(readLengthStore) / 10;
    std::cerr << "Progress: ";
    for (unsigned readId = 0; readId < length(readLengthStore); ++readId) {
        if (tenPercent > 0u && readId % tenPercent == 0u)
            std::cerr << readId / tenPercent * 10 << '%';
        else if (tenPercent > 5u && readId % (tenPercent / 5) == 0u)
            std::cerr << '.';

        std::sort(begin(errorCurves[readId], Standard()), end(errorCurves[readId], Standard()));
        fillGaps(errorCurves[readId]);
        smoothErrorCurve(errorCurves[readId]);

        // Compute relative min score for the read.
        String<WeightedMatch> filtered;
        int maxError = (int)floor(options.maxError / 100.0 * readLengthStore[readId]);
        if (options.oracleSamMode) {
            SEQAN_ASSERT_NEQ(readAlignmentDistances[readId], -1);
            maxError = readAlignmentDistances[readId];
        }
        int relativeMinScore = (int)ceilAwayFromZero(100.0 * -maxError / readLengthStore[readId]);

        // Filter out low scoring ones.
        typedef typename Iterator<String<WeightedMatch> >::Type TIterator;
        for (TIterator it = begin(errorCurves[readId]); it != end(errorCurves[readId]); ++it) {
            if (value(it).distance >= relativeMinScore) {
                #if ENABLED
                std::cerr << fragments.readNameStore[readId] << " accepting  " << value(it) << std::endl;
                #endif
                appendValue(filtered, value(it));
            } else {
                #if ENABLED
                std::cerr << fragments.readNameStore[readId] << " discarding " << value(it) << std::endl;
                #endif
            }
        }
        move(errorCurves[readId], filtered);
    }
    std::cerr << "100% DONE\n"
              << "\nTook: " << sysTime() - startTime << " s\n";

    return 0;
}

// Entry point for the gold standard building subprogram.
int buildGoldStandard(Options<BuildGoldStandard> const & options)
{
    double startTime = 0;  // For measuring time below.

    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache; // TODO(holtgrew): Remove?
    
    std::cerr << "==============================================================================\n"
              << "                RABEMA - Read Alignment BEnchMArk\n"
              << "==============================================================================\n"
              << "                      Building Gold Standard\n"
              << "==============================================================================\n"
              << "\n"
              << "____LOADING FILES_____________________________________________________________\n\n";

    // =================================================================
    // Load contigs.
    // =================================================================
    startTime = sysTime();
    std::cerr << "Reference Index      " << options.referenceSeqFilename << ".fai ...";
    FaiIndex faiIndex;
    if (load(faiIndex, toCString(options.referenceSeqFilename)) != 0)
    {
        std::cerr << "Could not read reference FAI index.\n";
        return 1;
    }
    std::cerr << " OK\n";
    
    // Open SAM file and read in header.
    std::cerr << "Alignments           " << options.perfectMapFilename << " (header) ...";
    std::ifstream inSam(toCString(options.perfectMapFilename), std::ios_base::in | std::ios_base::binary);
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
    // Check that the SAM file is sorted by COORDINATE.
    if (getSortOrder(samHeader) != BAM_SORT_COORDINATE)
    {
        std::cerr << "SAM file not sorted by 'coordinate'!\n";
        return 1;
    }
    std::cerr << " OK\n";
    std::cerr << "\nTook " << sysTime() - startTime << "s\n";

    // =================================================================
    // Build point-wise error curve.
    // =================================================================
    startTime = sysTime();
    TErrorCurves errorCurves;
    int res = 0;
    StringSet<CharString> readNameStore;
    if (options.distanceFunction == "edit")
        res = matchesToErrorFunction(errorCurves, samReader, samIOContext, readNameStore, refNameStore, faiIndex, options, MyersUkkonenReads());
    else // options.distanceFunction == "hamming"
        res = matchesToErrorFunction(errorCurves, samReader, samIOContext, readNameStore, refNameStore, faiIndex, options, HammingSimple());
    if (res != 0)
        return 1;
    if (options.verbosity >= 2)
        std::cerr << "[timer] building error curve: " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Verify result when enabled.
    // =================================================================
    if (options.verify) {
        // TODO(holtgrew): Re-enable?
        std::cerr << "VERIFICATION DISABLED AT THE MOMENT!\n";
#if 0
        startTime = sysTime();
        if (options.oracleSamMode) {
            std::cerr << "WARNING: Not verifying since we run in oracle Sam mode!" << std::endl;
        } else {
            bool valid;
            if (options.distanceFunction == "edit")
                valid = verifyMatchesToErrorFunctionResults(fragmentStore, errorCurves, options, MyersUkkonenReads());
            else // options.distanceFunction == "hamming"
                valid = verifyMatchesToErrorFunctionResults(fragmentStore, errorCurves, options, HammingSimple());
            
            if (!valid) {
                std::cerr << "ERROR: Result does not validate!" << std::endl;
                return 1;
            }
            // else...
            std::cerr << "Result DID validate!" << std::endl;
        }
        if (options.verbosity >= 2)
            std::cerr << "[timer] verification: " << sysTime() - startTime << " s" << std::endl;
#endif  // #if 0
    }

    // =================================================================
    // Convert points in error curves to intervals and write them to
    // stdout or a file.
    // =================================================================
    startTime = sysTime();
    String<WitRecord> witRecords;
    typedef Iterator<String<WitRecord>, Standard>::Type TWitRecordIterator;
    intervalizeErrorCurves(witRecords, errorCurves, readNameStore, refNameStore, options);
    std::cerr << "Took " << sysTime() - startTime << " s\n";
    // The two alternatives are equivalent after opening the file.
    startTime = sysTime();
    std::cerr << "\n____WRITING OUTPUT____________________________________________________________\n\n";
    if (options.outFileName == "-") {
        std::cerr << "Writing to stdout ...\n";
        writeWitHeader(std::cout);
        writeWitComment(std::cout, WIT_COLUMN_NAMES);
        for (TWitRecordIterator it = begin(witRecords, Standard());
             it != end(witRecords, Standard()); ++it)
            std::cout << *it << std::endl;
        std::cerr << "DONE\n";
    } else {
        std::cerr << "Writing to " << options.outFileName << " ...";
        std::fstream fstrm(toCString(options.outFileName), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open out file \"" << options.outFileName << "\""
                      << std::endl;
            return 1;
        }
        writeWitHeader(fstrm);
        writeWitComment(fstrm, WIT_COLUMN_NAMES);
        for (TWitRecordIterator it = begin(witRecords, Standard());
             it != end(witRecords, Standard()); ++it)
            fstrm << *it << std::endl;
        std::cerr << " DONE\n";
    }
    std::cerr << "\n Took " << sysTime() - startTime << " s\n";

    return 0;
}

#endif  // #ifndef APPS_RABEMA_BUILD_GOLD_STANDARD_H_
