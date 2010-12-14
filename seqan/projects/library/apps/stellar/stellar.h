 /*==========================================================================
                     STELLAR - Fast Local Alignment

 ============================================================================
  Copyright (C) 2010 by Birte Kehr

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#ifndef SEQAN_HEADER_STELLAR_H
#define SEQAN_HEADER_STELLAR_H

#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds2.h>
#include "stellar_types.h"
#include "stellar_extension.h"

using namespace seqan;


struct _VerifyAllLocal;
typedef Tag<_VerifyAllLocal> const AllLocal;

struct _VerifyBestLocal;
typedef Tag<_VerifyBestLocal> const BestLocal;

struct _VerifyBandedGlobal;
typedef Tag<_VerifyBandedGlobal> const BandedGlobal;

struct _VerifyBandedGlobalExtend;
typedef Tag<_VerifyBandedGlobalExtend> const BandedGlobalExtend;

///////////////////////////////////////////////////////////////////////////////

inline bool _verifyFast(BestLocal) {
	return true;
}

template<typename TTag>
inline bool _verifyFast(TTag) {
	return false;
}

///////////////////////////////////////////////////////////////////////////////
// Appends a segment of only error positions from align to queue.
template<typename TAlign, typename TPos, typename TScoreValue>
inline void
_appendNegativeSegment(TAlign const & align, 
                       TPos & pos, TPos len, 
                       Score<TScoreValue> const & scoreMatrix, 
                       String<Triple<TPos, TPos, TScoreValue> > & queue) {
SEQAN_CHECKPOINT
    typedef Triple<TPos, TPos, TScoreValue> TMerger;
    TPos beginPos = pos;

    TScoreValue score = 0;
    while (pos < len) {
        if (isGap(row(align, 0), pos) || isGap(row(align, 1), pos)) {
            score += scoreGap(scoreMatrix);
        } else if (value(row(align, 0), pos) != value(row(align, 1), pos)) {
            score += scoreMismatch(scoreMatrix);
        } else {
            break;
        }
        ++pos;
    }
    if (pos == len) appendValue(queue, TMerger(beginPos, pos, infimumValue<TScoreValue>()+1));
    else appendValue(queue, TMerger(beginPos, pos, score));
}

///////////////////////////////////////////////////////////////////////////////
// Appends a segment of only matching positions from align to queue.
template<typename TAlign, typename TPos, typename TScoreValue>
inline void
_appendPositiveSegment(TAlign const & align, 
                       TPos & pos, TPos len, 
                       Score<TScoreValue> const & scoreMatrix, 
                       String<Triple<TPos, TPos, TScoreValue> > & queue) {
SEQAN_CHECKPOINT
    if (pos == len) return;
    typedef Triple<TPos, TPos, TScoreValue> TMerger;
    TPos beginPos = pos;

    TScoreValue score = 0;
    while ((pos < len) && 
           (!isGap(row(align, 0), pos) &&
            !isGap(row(align, 1), pos) &&
            (value(row(align, 0), pos) == value(row(align, 1), pos)))) {
        score += scoreMatch(scoreMatrix);
        ++pos;
    }
    appendValue(queue, TMerger(beginPos, pos, score));
}

///////////////////////////////////////////////////////////////////////////////
// See Lemma 5 in Zhang et al., 1999.
template<typename TMerger>
inline bool
_negativeMerge(String<TMerger> & queue) {
SEQAN_CHECKPOINT
    typedef typename TMerger::T1 TPos;
    TPos len = length(queue);
    if (len < 3) return false;

    TMerger cd = value(queue, len-1);
    TMerger bc = value(queue, len-2);
    TMerger ab = value(queue, len-3);

    if ((bc.i3 < 0) || (bc.i3 >= abs(_max(ab.i3, cd.i3)))) {
        return false;
    } else {
        String<TMerger> newMerger;
        appendValue(newMerger, TMerger(ab.i1, cd.i2, ab.i3 + bc.i3 + cd.i3)); 
        replace(queue, len-3, len, newMerger);

        return true;
    }
}

///////////////////////////////////////////////////////////////////////////////
// See Lemma 6 in Zhang et al., 1999.
template<typename TMerger>
inline bool
_positiveMerge(String<TMerger> & queue) {
SEQAN_CHECKPOINT
    typedef typename TMerger::T1 TPos;
    TPos len = length(queue);
    if (len < 5) return false;

    TMerger ef = value(queue, len-1);
    TMerger de = value(queue, len-2);
    TMerger cd = value(queue, len-3);
    TMerger bc = value(queue, len-4);
    TMerger ab = value(queue, len-5);

    if ((cd.i3 >= 0) || (cd.i3 < _max(ab.i3, ef.i3))) {
        return false;
    } else {
        String<TMerger> newMerger;
        appendValue(newMerger, TMerger(bc.i1, de.i2, bc.i3 + cd.i3 + de.i3)); 
        replace(queue, len-4, len-1, newMerger);

        return true;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Implements the algorithm from Zhang et al. in Bioinformatics, 1999: "Post-processing long pairwise alignments".
// Splits an alignment into sub-alignments that contain no x-Drop.
template<typename TAlign, typename TScoreValue, typename TScoreValue1, typename TScoreValue2>
void
_splitAtXDrops(TAlign const & align,
              Score<TScoreValue> & scoreMatrix,
              TScoreValue1 scoreDropOff,
              TScoreValue2 minScore,
              String<TAlign> & alignmentString) {
SEQAN_CHECKPOINT
    typedef typename Position<Row<TAlign> >::Type TPos;
    typedef Triple<TPos, TPos, TScoreValue> TMerger;
    
    // initialization
    String<TMerger> queue;
    TPos pos = _min(toViewPosition(row(align, 0), clippedBeginPosition(row(align, 0))),
                    toViewPosition(row(align, 1), clippedBeginPosition(row(align, 1))));
    appendValue(queue, TMerger(pos, pos, infimumValue<TScoreValue1>()+1));

    TPos aliLength = _max(toViewPosition(row(align, 0), clippedEndPosition(row(align, 0))),
                          toViewPosition(row(align, 1), clippedEndPosition(row(align, 1))));
    TPos len;
    while ((pos < aliLength) || (length(queue) > 1)) {
        // construct useful tree
        if (!_negativeMerge(queue)) {
            if (!_positiveMerge(queue)) {
                _appendPositiveSegment(align, pos, aliLength, scoreMatrix, queue);
                _appendNegativeSegment(align, pos, aliLength, scoreMatrix, queue);
            }
        }

        // check for x-Drop
        len = length(queue);
        if ((len == 3) && (value(queue, 2).i3 < scoreDropOff * (-1))) {
            if (value(queue, 1).i3 >= minScore) {
				// create new sub-alignment
                TAlign ali(align);
                TPos begin0 = toSourcePosition(row(ali, 0), value(queue, 1).i1);
                TPos begin1 = toSourcePosition(row(ali, 1), value(queue, 1).i1);
                TPos end0 = toSourcePosition(row(ali, 0), value(queue, 1).i2);
                TPos end1 = toSourcePosition(row(ali, 1), value(queue, 1).i2);
                setClippedBeginPosition(row(ali, 0), begin0);
                setClippedBeginPosition(row(ali, 1), begin1);
				setBeginPosition(row(ali, 0), 0);
				setBeginPosition(row(ali, 1), 0);
                setClippedEndPosition(row(ali, 0), end0);
                setClippedEndPosition(row(ali, 1), end1);

                // append sub-alignment
                appendValue(alignmentString, ali);
            }
            replace(queue, 0, 2, String<TMerger>());
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Checks whether two matches overlap in seq2 and 
//  whether the non-overlaping parts are shorter than minLength.
template<typename TMatch, typename TSize>
bool
checkOverlap(TMatch & matchA, TMatch & matchB, TSize minLength) {
SEQAN_CHECKPOINT
	// check overlap in seq2
	if (matchA.begin2 >= matchB.begin2) {
		if (matchA.end2 >= matchB.end2) {
			// check length of non-overlapping parts of both matches
			if ((TSize)matchA.begin2 - (TSize)matchB.begin2 >= minLength &&
				(TSize)matchA.end2 - (TSize)matchB.end2 >= minLength) {
				return false;
			}
		}
		// check whether offset is the same in both sequences
		if (toViewPosition(matchA.row2, matchA.begin2) - toViewPosition(matchB.row2, matchB.begin2) != 
			toViewPosition(matchA.row1, matchA.begin1) - toViewPosition(matchB.row1, matchB.begin1)) {
			return false;
		}
	} else {
		if (matchA.end2 < matchB.end2) {
			// check length of non-overlapping parts of both matches
			if ((TSize)matchB.begin2 - (TSize)matchA.begin2 >= minLength &&
				(TSize)matchB.end2 - (TSize)matchA.end2 >= minLength) {
				return false;
			}
		}
		// check whether offset is the same in both sequences
		if (toViewPosition(matchB.row2, matchB.begin2) - toViewPosition(matchA.row2, matchA.begin2) != 
			toViewPosition(matchB.row1, matchB.begin1) - toViewPosition(matchA.row1, matchA.begin1)) {
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// Marks matches that overlap in both sequences with a longer match as invalid.
template<typename TSequence, typename TId, typename TSize>
void
maskOverlaps(String<StellarMatch<TSequence, TId> > & matches,
			 TSize minLength) {
SEQAN_CHECKPOINT
	typedef StellarMatch<TSequence, TId>					TMatch;
	typedef typename TMatch::TPos							TPos;
	typedef typename Iterator<String<TMatch>, Rooted>::Type	TIter;
	typedef typename Iterator<String<TSize>, Rooted>::Type	TOverlapIter;

	// sort matches by begin position in row0
	sortMatches(matches, LessPos<TMatch>());

	// list of "open" matches in row0
	String<TSize> overlaps;

	TIter it = begin(matches);
	TIter itEnd = end(matches);

	for (; it != itEnd; ++it) {
		if ((*it).id == TMatch::INVALID_ID) continue;

		TOverlapIter overlapIt = begin(overlaps);
		TOverlapIter overlapEnd = end(overlaps);
		
		while (overlapIt != overlapEnd && matches[*overlapIt].end1 >= (*it).end1) {
			if (matches[*overlapIt].id == TMatch::INVALID_ID) {
				++overlapIt;
				continue;
			}
			if (checkOverlap(*it,  matches[*overlapIt], minLength)) {
				// set shorter match invalid (*it is here shorter than *overlapIt)
				(*it).id = TMatch::INVALID_ID;
			}
			++overlapIt;
		}

		TPos insertPos = position(overlapIt);
		while (overlapIt != overlapEnd &&  matches[*overlapIt].end1 >= (*it).begin1) {
			if (matches[*overlapIt].id == TMatch::INVALID_ID) {
				++overlapIt;
				continue;
			}
			if (checkOverlap(*it,  matches[*overlapIt], minLength)) {
				// set shorter match invalid
				if (length((*it)) > length(matches[*overlapIt])) {
					 matches[*overlapIt].id = TMatch::INVALID_ID;
				} else {
					(*it).id = TMatch::INVALID_ID;
				}
			}
			++overlapIt;
		}
		resize(overlaps, position(overlapIt));

		// insert match into list
		if ((*it).id != TMatch::INVALID_ID)
			insertValue(overlaps, insertPos, position(it));
	}
}

///////////////////////////////////////////////////////////////////////////////
// Removes matches that are marked as invalid, and then keeps only the numMatches best matches.
template<typename TSequence, typename TId, typename TSize>
void
compactMatches(String<StellarMatch<TSequence, TId> > & matches, TSize numMatches) {
	typedef StellarMatch<TSequence, TId>						TMatch;
	typedef typename Iterator<String<TMatch>, Standard>::Type	TIterator;
	
	// sort matches by length (and validity)
	sortMatches(matches, LessLength<TMatch>());
	
	// count valid matches
	TSize num = 0;

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());

	for(; it != itEnd; ++it) {
		if ((*it).id != TMatch::INVALID_ID)
			++num;
	}

	// keep only valid and longest matches
	resize(matches, _min(num, numMatches));
}

///////////////////////////////////////////////////////////////////////////////
// Appends a match to matches container and removes overlapping matches if threshold is reached.
template<typename TSource, typename TId, typename TSize, typename TSize1>
inline bool
_insertMatch(QueryMatches<StellarMatch<TSource, TId> > & queryMatches,
			 StellarMatch<TSource, TId> match,
			 TSize minLength,
			 TSize1 disableThresh,
			 TSize1 compactThresh,
			 TSize1 numMatches) {
SEQAN_CHECKPOINT

	appendValue(queryMatches.matches, match);

	if (length(queryMatches.matches) > disableThresh) {
		queryMatches.disabled = true;
		clear(queryMatches.matches);
		return false;
	}
	if (length(queryMatches.matches) > compactThresh) {
		maskOverlaps(queryMatches.matches, minLength);		// remove overlaps and duplicates
		compactMatches(queryMatches.matches, numMatches);	// keep only the <numMatches> longest matches

		// raise compact threshold if many matches are kept
		if ((length(queryMatches.matches) << 1) > compactThresh)
			compactThresh += (compactThresh >> 1);
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Conducts banded alignment on swift hit and extracts longest contained eps-match.
template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TDelta, 
         typename TDrop, typename TSize1, typename TId, typename TSequence>
void
verifySwiftHit(TInfixA const & a,
			   TInfixB const & b,
			   TEpsilon eps,
			   TSize minLength,
			   TDrop /*xDrop*/,
			   TDelta delta,
			   TSize1 & disableThresh,
			   TSize1 & compactThresh,
			   TSize1 numMatches,
			   TId & queryId,
			   QueryMatches<StellarMatch<TSequence, TId> > & matches,
			   BandedGlobal) {
SEQAN_CHECKPOINT
	typedef typename StellarMatch<TSequence, TId>::TAlign TAlign;

    // define a scoring scheme
    typedef int TScore;
    TScore match = 1;
    TScore mismatchIndel = (TScore)_max((TScore) ceil(-1/eps) + 1, -(TScore)length(host(a)));
    Score<TScore> scoreMatrix(match, mismatchIndel, mismatchIndel);

    // diagonals for banded alignment
    __int64 upperDiag = 0;
    __int64 lowerDiag = endPosition(a) - (__int64)endPosition(b) - beginPosition(a) + beginPosition(b);
    if (beginPosition(b) == 0) upperDiag = lowerDiag + delta;
    if (endPosition(b) == endPosition(host(b))) lowerDiag = -(__int64)delta;

	// banded alignment on parallelogram
	Align<TInfixB> bandedAlign;
    resize(rows(bandedAlign), 2);
    assignSource(row(bandedAlign, 0), a);
    assignSource(row(bandedAlign, 1), b);
	StringSet<TInfixB> str;
	appendValue(str, a);
	appendValue(str, b);
	globalAlignment(bandedAlign, str, scoreMatrix, lowerDiag, upperDiag, BandedNeedlemanWunsch());

	longestEpsMatch(bandedAlign, minLength, eps);

	// integrate alignment in object of type TAlign
	TAlign align;
	resize(rows(align), 2);
	setSource(row(align, 0), host(a));
	setSource(row(align, 1), host(b));
	integrateAlign(align, bandedAlign);

	// set begin and end positions of align
	setClippedBeginPosition(row(align, 0), beginPosition(a) + clippedBeginPosition(row(bandedAlign, 0)));
	setClippedBeginPosition(row(align, 1), beginPosition(b) + clippedBeginPosition(row(bandedAlign, 1)));
	setBeginPosition(row(align, 0), 0);
	setBeginPosition(row(align, 1), 0);
	setClippedEndPosition(row(align, 0), beginPosition(a) + clippedEndPosition(row(bandedAlign, 0)));
	setClippedEndPosition(row(align, 1), beginPosition(b) + clippedEndPosition(row(bandedAlign, 1)));

	if ((TSize)length(row(align, 0)) < minLength)
		return;

	// insert eps-match in matches string
	StellarMatch<TSequence, TId> m(align, queryId);
	_insertMatch(matches, m, minLength, disableThresh, compactThresh, numMatches);
}

///////////////////////////////////////////////////////////////////////////////
// Conducts banded alignment on swift hit, extends alignment, and extracts longest contained eps-match.
template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TDelta, 
         typename TDrop, typename TSize1, typename TId, typename TSequence>
void
verifySwiftHit(TInfixA const & a,
			   TInfixB const & b,
			   TEpsilon eps,
			   TSize minLength,
			   TDrop xDrop,
			   TDelta delta,
			   TSize1 & disableThresh,
			   TSize1 & compactThresh,
			   TSize1 numMatches,
			   TId & queryId,
			   QueryMatches<StellarMatch<TSequence, TId> > & matches,
			   BandedGlobalExtend) {
SEQAN_CHECKPOINT
	typedef typename StellarMatch<TSequence, TId>::TAlign TAlign;

    // define a scoring scheme
    typedef int TScore;
    TScore match = 1;
    TScore mismatchIndel = (TScore)_max((TScore) ceil(-1/eps) + 1, -(TScore)length(host(a)));
    Score<TScore> scoreMatrix(match, mismatchIndel, mismatchIndel);
    TScore scoreDropOff = (TScore) _max((TScore) xDrop * (-mismatchIndel), infimumValue<TScore>()+1);

    // diagonals for banded alignment
    __int64 upperDiag = 0;
    __int64 lowerDiag = endPosition(a) - (__int64)endPosition(b) - beginPosition(a) + beginPosition(b);
    if (beginPosition(b) == 0) upperDiag = lowerDiag + delta;
    if (endPosition(b) == endPosition(host(b))) lowerDiag = -(__int64)delta;

	// banded alignment on parallelogram
	Align<TInfixB> bandedAlign;
    resize(rows(bandedAlign), 2);
    assignSource(row(bandedAlign, 0), a);
    assignSource(row(bandedAlign, 1), b);
	StringSet<TInfixB> str;
	appendValue(str, a);
	appendValue(str, b);
	globalAlignment(bandedAlign, str, scoreMatrix, lowerDiag, upperDiag, BandedNeedlemanWunsch());

	// create alignment object for the complete sequences
	TAlign align;
	resize(rows(align), 2);
	setSource(row(align, 0), host(a));
	setSource(row(align, 1), host(b));

	// extend alignment and obtain longest contained eps-match
	// TODO: something is wrong here, e.g. extract around seed, but also something else
	if (!_extendAndExtract(bandedAlign, scoreDropOff, scoreMatrix, a, b, EXTEND_BOTH, minLength, eps, align))
		return;

	// insert eps-match in matches string
	StellarMatch<TSequence, TId> m(align, queryId);
	_insertMatch(matches, m, minLength, disableThresh, compactThresh, numMatches);
}

///////////////////////////////////////////////////////////////////////////////
// Conducts banded local alignment on swift hit (= computes eps-cores),
//  splits eps-cores at X-drops, and calls _extendAndExtract for extension of eps-cores
template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TDelta, typename TDrop,
         typename TSize1, typename TId, typename TSequence, typename TTag>
void
verifySwiftHit(TInfixA const & a,
               TInfixB const & b,
               TEpsilon eps,
               TSize minLength,
               TDrop xDrop,
               TDelta delta,
			   TSize1 & disableThresh,
			   TSize1 & compactThresh,
			   TSize1 numMatches,
			   TId & queryId,
			   QueryMatches<StellarMatch<TSequence, TId> > & matches,
			   TTag tag) {
SEQAN_CHECKPOINT
	typedef typename StellarMatch<TSequence, TId>::TAlign TAlign;

    TSize maxLength = 1000000000;
    if ((TSize)length(a) > maxLength) {
        std::cerr << "Warning: SWIFT hit <" << beginPosition(a) << "," << endPosition(a);
        std::cerr << "> , <" << beginPosition(b) << "," << endPosition(b);
		std::cerr << "> too long. Verification skipped.\n" << std::flush;
        return;
    }

    // define a scoring scheme
    typedef int TScore;
    TScore match = 1;
    TScore mismatchIndel = (TScore)_max((TScore) ceil(-1/eps) + 1, -(TScore)length(host(a)));
    Score<TScore> scoreMatrix(match, mismatchIndel, mismatchIndel);
    TScore scoreDropOff = (TScore) _max((TScore) xDrop * (-mismatchIndel), infimumValue<TScore>()+1);

    // calculate minimal score for local alignments
    TEpsilon e = floor(eps*minLength);
    TSize minLength1 = _max(0, (TSize)ceil((e+1) / eps));
    TEpsilon e1 = floor(eps*minLength1);
    TSize minScore = _min((TSize)ceil((minLength-e) / (e+1)), (TSize)ceil((minLength1-e1) / (e1+1)));

    // diagonals for banded local alignment
    __int64 upperDiag = 0;
    __int64 lowerDiag = endPosition(a) - (__int64)endPosition(b) - beginPosition(a) + beginPosition(b);
	if (beginPosition(b) == 0) {
		if (endPosition(b) == endPosition(host(b))) {
			// TODO: is it possible to get a smaller band in this case?
			upperDiag = delta;
			lowerDiag = -(__int64)delta;
		} else
			upperDiag = lowerDiag + delta;
	} else if (endPosition(b) == endPosition(host(b)))
		lowerDiag = -(__int64)delta;

	// banded local alignment
    LocalAlignmentFinder<> finder = LocalAlignmentFinder<>();
	Align<TInfixB> localAlign;
    resize(rows(localAlign), 2);
    assignSource(row(localAlign, 0), a);
    assignSource(row(localAlign, 1), b);

	while (localAlignment(localAlign, finder, scoreMatrix, minScore, lowerDiag, upperDiag, BandedWatermanEggert())) {

        // split local alignments containing an X-drop
        String<Align<TInfixB> > seedAlignments;
        _splitAtXDrops(localAlign, scoreMatrix, scoreDropOff, minScore, seedAlignments);

        typename Iterator<String<Align<TInfixB> > >::Type aliIt = begin(seedAlignments);
        while (aliIt != end(seedAlignments)) {
			// create alignment object for the complete sequences
			TAlign align;
			resize(rows(align), 2);
			setSource(row(align, 0), host(a));
			setSource(row(align, 1), host(b));

            // determine extension direction
            ExtensionDirection direction;
            if (length(seedAlignments) == 1) direction = EXTEND_BOTH;
            else if (aliIt == begin(seedAlignments)) direction = EXTEND_RIGHT;
            else if (aliIt == end(seedAlignments)-1) direction = EXTEND_LEFT;
            else direction = EXTEND_NONE;

			// extend alignment and obtain longest contained eps-match
			if (!_extendAndExtract(*aliIt, scoreDropOff, scoreMatrix, a, b, direction, minLength, eps, align)) {
				aliIt++;
				continue;
			}

            // insert eps-match in matches string
			StellarMatch<TSequence, TId> m(align, queryId);
            if(!_insertMatch(matches, m, minLength, disableThresh, compactThresh, numMatches)) return;
            ++aliIt;
        }
		if (_verifyFast(tag)) break;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Calls swift filter and verifies swift hits. = Computes eps-matches.
template<typename TText, typename TIndex, typename TSize, typename TDrop, typename TSize1,
         typename TSource, typename TId, typename TTag>
int stellar(Finder<TText, Swift<SwiftLocal> > & finder,
               Pattern<TIndex, Swift<SwiftLocal> > & pattern,
               double epsilon,
               TSize minLength,
               TDrop xDrop,
			   TSize1 & disableThresh,
			   TSize1 & compactThresh,
			   TSize1 numMatches,
			   StringSet<TId> & queryIDs,
               StringSet<QueryMatches<StellarMatch<TSource, TId> > > & matches,
			   TTag tag) {
SEQAN_CHECKPOINT
	typedef StellarMatch<TSource, TId> TMatch;
	typedef typename Infix<TText>::Type TInfix;

    resize(matches, countSequences(needle(pattern)));
    TSize numSwiftHits = 0;

    TSize maxLength = 0;
    TSize totalLength = 0;

	while (find(finder, pattern, epsilon, minLength)) {
		TInfix finderInfix = infix(finder);

        ++numSwiftHits;
        totalLength += length(finderInfix);
        if ((TSize)length(finderInfix) > maxLength) maxLength = length(finderInfix);

		if (value(matches, pattern.curSeqNo).disabled) continue;

		TInfix patternInfix = infix(pattern, value(host(needle(pattern)), pattern.curSeqNo));
		////Debug output:
		//std::cout << beginPosition(finderInfix) << ",";
		//std::cout << endPosition(finderInfix) << "  ";
		//std::cout << beginPosition(patternInfix) << ",";
		//std::cout << endPosition(patternInfix) << std::endl;

        // verification
		verifySwiftHit(finderInfix, patternInfix, epsilon, minLength, xDrop,
					   pattern.bucketParams[0].delta + pattern.bucketParams[0].overlap, disableThresh, compactThresh,
					   numMatches, value(queryIDs, pattern.curSeqNo), value(matches, pattern.curSeqNo), tag);

	}

	//if (numSwiftHits > 0) {
	//	std::cout << "    # hits           : " << numSwiftHits << std::endl;
	//	std::cout << "    Longest hit      : " << maxLength << std::endl;
	//	std::cout << "    Avg hit length   : " << totalLength/numSwiftHits << std::endl;
	//}
	
	typedef typename Iterator<StringSet<QueryMatches<TMatch> >, Standard>::Type TIterator;
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());

	for(; it < itEnd; ++it) {
		if (length(*it) > 0 && !(*it).disabled) {
			maskOverlaps((*it).matches, minLength);		// remove overlaps and duplicates
			compactMatches((*it).matches, numMatches);	// keep only the <numMatches> longest matches
		}
	}
    
    return numSwiftHits;
}

#endif
