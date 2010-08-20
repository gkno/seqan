 /*==========================================================================
                     SwiftLocal - Fast Local Alignment

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

#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include "swift_local_types.h"

using namespace seqan;


struct _VerifyAllLocal;
typedef Tag<_VerifyAllLocal> const AllLocal;

struct _VerifyBestLocal;
typedef Tag<_VerifyBestLocal> const BestLocal;

struct _VerifyBandedGlobal;
typedef Tag<_VerifyBandedGlobal> const BandedGlobal;

struct _VerifyBandedGlobalExtend;
typedef Tag<_VerifyBandedGlobalExtend> const BandedGlobalExtend;


inline bool _verifyFast(BestLocal) {
	return true;
}

template<typename TTag>
inline bool _verifyFast(TTag) {
	return false;
}

// returns true if align has a match at pos, otherwise false
template<typename TSource, typename TSize>
inline bool
isMatch(Align<TSource> const & align, TSize pos) {
SEQAN_CHECKPOINT
    if(isGap(row(align, 0), pos)) {
        return false;
    } else if(isGap(row(align, 1), pos)) {
        return false;
    } else if(row(align, 0)[pos] != row(align, 1)[pos]) {
        return false;
    } else {
        return true;
    }
}

///////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////

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
    TPos pos = _min(toViewPosition(row(align, 0), sourceBeginPosition(row(align, 0))),
                    toViewPosition(row(align, 1), sourceBeginPosition(row(align, 1))));
    appendValue(queue, TMerger(pos, pos, infimumValue<TScoreValue1>()+1));

    TPos aliLength = _max(toViewPosition(row(align, 0), sourceEndPosition(row(align, 0))),
                          toViewPosition(row(align, 1), sourceEndPosition(row(align, 1))));
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
                // append new sub-alignment
                TPos begin = value(queue, 1).i1;
                TPos end = value(queue, 1).i2;

                TAlign ali(align);
                setSourceBeginPosition(row(ali, 0), toSourcePosition(row(ali, 0), begin));
                setSourceBeginPosition(row(ali, 1), toSourcePosition(row(ali, 1), begin));
                setSourceEndPosition(row(ali, 0), toSourcePosition(row(ali, 0), end));
                setSourceEndPosition(row(ali, 1), toSourcePosition(row(ali, 1), end));
                appendValue(alignmentString, ali);
            }
            replace(queue, 0, 2, String<TMerger>());
        }
    }
}

template<typename TString, typename TSeed, typename TScore, typename TDiag, typename TAlign>
void
_bandedInfixAlignment(StringSet<TString> & str,
					 TSeed & seed,
					 TScore & scoreMatrix,
					 TDiag startDiag,
					 TAlign & align) {
SEQAN_CHECKPOINT
	if (length(str[0]) != 0 && length(str[1]) != 0) {
		_Align_Traceback<unsigned> trace;
		globalAlignment(trace, str, scoreMatrix, startDiag - leftDiagonal(seed),
			startDiag - rightDiagonal(seed), BandedNeedlemanWunsch());
		
		Align<TString> infixAlign;
		resize(rows(infixAlign), 2);
		assignSource(row(infixAlign, 0), str[0]);
		assignSource(row(infixAlign, 1), str[1]);

		_pump_trace_2_Align(infixAlign, trace);
		integrateAlign(align, infixAlign);
	}
}

template<typename TSource, typename TPos>
void
_fillGapsString(Align<TSource> const & align,
                String<Triple<TPos, TPos, TPos> > & gaps) {
SEQAN_CHECKPOINT
    typedef Triple<TPos, TPos, TPos> TGapInfo;
    TPos totalErrors = 0;
	typename Row<Align<TSource> >::Type row0 = row(align, 0);
    TPos i = 0;//toViewPosition(row0, beginPosition(row0));
	TPos endPos = endPosition(row0);
    TPos gapBegin = i;

    // append gap starting at beginPosition (also if its length is 0!)
    while(i < endPos && !isMatch(align, i)) {
        ++i;
        ++totalErrors;
    }
    appendValue(gaps, TGapInfo(gapBegin, i, totalErrors));

    // iterate over alignment and append gaps
    while (i < endPos) {
        // skip matches
        while(i < endPos && isMatch(align, i)) {
            ++i;
        }
        gapBegin = i;
        // skip and count mismatches/indels
        while(i < endPos && !isMatch(align, i)) {
            ++i;
            ++totalErrors;
        }
        appendValue(gaps, TGapInfo(gapBegin, i, totalErrors));
    }
    /*for(unsigned l = 0; l < length(gaps); ++l) {
        std::cout << gaps[l].i1 << "  " << gaps[l].i2 << "  " << gaps[l].i3 << std::endl;
    }*/
}

// checks the error rate of the fragment between end of left and start of right
template<typename TPos, typename TFloat>
inline bool
_isEpsMatch(Triple<TPos, TPos, TPos> const & left,
           Triple<TPos, TPos, TPos> const & right,
           TFloat eps) {
SEQAN_CHECKPOINT
    // compute mismatches/indels and length
    TPos errors = right.i3 - left.i3 - (right.i2 - right.i1);
    TPos length = right.i1 - left.i2;

    // check error rate
    return errors/(TFloat)(length) <= eps;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align that spans seedBegin and seedEnd and sets the view positions of
// align to start and end position of the longest epsilon match
template<typename TSource, typename TPos, typename TSize, typename TFloat>
void 
longestEpsMatch(Align<TSource> & align,
						  TPos seedBegin,
						  TPos seedEnd,
						  TSize matchMinLength,
						  TFloat epsilon) {
SEQAN_CHECKPOINT
    // Preprocessing: compute and store gaps and lengths
    // A gap is a triple of gap begin position, gap end position, and total number of errors in sequence from begin
    //   to end position of this gap.
    typedef typename Position<Align<TSource> >::Type TPosition;
    typedef String<Triple<TPosition, TPosition, TPosition> > TGapsString;
    TGapsString gaps;
    _fillGapsString(align, gaps);

    // Identify longest eps match by iterating over combinations of left and right positions
    typename Iterator<TGapsString >::Type rightIt = end(gaps) - 1;
    typename Iterator<TGapsString >::Type leftIt = begin(gaps);

    TPosition beginPos = 0;
    TPosition endPos = 0;
    TSize minLength = matchMinLength - 1;

	typename Row<Align<TSource > >::Type row0 = row(align, 0);

	TPosition seedBeginView = (TPosition)toViewPosition(row0, seedBegin + sourceBeginPosition(row0));
	TPosition seedEndView = (TPosition)toViewPosition(row0, seedEnd + sourceBeginPosition(row0));
	if (sourceBeginPosition(row0) > 0) {
		seedBeginView -= (TPosition)toViewPosition(row0, sourceBeginPosition(row0));
		seedEndView -= (TPosition)toViewPosition(row0, sourceBeginPosition(row0));
	}
    
    while ((*leftIt).i2 + minLength < (*rightIt).i1 && (*leftIt).i2 <= seedBeginView) {
        while ((*leftIt).i2 + minLength < (*rightIt).i1 && (*rightIt).i2 >= seedEndView) {
            if(_isEpsMatch(*leftIt, *rightIt, epsilon)) {
                beginPos = (*leftIt).i2;
                endPos = (*rightIt).i1;
                minLength = endPos - beginPos;
                break;
            }
            --rightIt;
        }
        rightIt = end(gaps) - 1;
        ++leftIt;
    }

    // Set view positions to the eps-match
	setSourceBeginPosition(row(align, 0), toSourcePosition(row(align, 0), beginPos));
	setSourceBeginPosition(row(align, 1), toSourcePosition(row(align, 1), beginPos));
	setBeginPosition(row(align, 0), beginPos);
	setBeginPosition(row(align, 1), beginPos);
	setSourceEndPosition(row(align, 0), toSourcePosition(row(align, 0), endPos));
	setSourceEndPosition(row(align, 1), toSourcePosition(row(align, 1), endPos));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align that spans seedBegin and seedEnd and sets the view positions of
// align to start and end position of the longest epsilon match
template<typename TSource, typename TSize, typename TFloat>
void 
longestEpsMatch(Align<TSource> & align,
                TSize matchMinLength,
                TFloat epsilon) {
SEQAN_CHECKPOINT
    // Preprocessing: compute and store gaps and lengths
    // A gap is a triple of gap begin position, gap end position, and total number of errors in sequence from begin
    //   to end position of this gap.
    typedef typename Position<Align<TSource> >::Type TPosition;
    typedef String<Triple<TPosition, TPosition, TPosition> > TGapsString;
    TGapsString gaps;
    _fillGapsString(align, gaps);

    // Identify longest eps match by iterating over combinations of left and right positions
    typename Iterator<TGapsString >::Type rightIt = end(gaps) - 1;
    typename Iterator<TGapsString >::Type leftIt = begin(gaps);

    TPosition beginPos = 0;
    TPosition endPos = 0;
    TSize minLength = matchMinLength - 1;
    
    while ((*leftIt).i2 + minLength < (*rightIt).i1) {
        while ((*leftIt).i2 + minLength < (*rightIt).i1) {
            if(_isEpsMatch(*leftIt, *rightIt, epsilon)) {
                beginPos = (*leftIt).i2;
                endPos = (*rightIt).i1;
                minLength = endPos - beginPos;
                break;
            }
            --rightIt;
        }
        rightIt = end(gaps) - 1;
        ++leftIt;
    }

    // Set view positions to the eps-match
	setSourceBeginPosition(row(align, 0), toSourcePosition(row(align, 0), beginPos));
	setSourceBeginPosition(row(align, 1), toSourcePosition(row(align, 1), beginPos));
	setBeginPosition(row(align, 0), beginPos);
	setBeginPosition(row(align, 1), beginPos);
	setSourceEndPosition(row(align, 0), toSourcePosition(row(align, 0), endPos));
	setSourceEndPosition(row(align, 1), toSourcePosition(row(align, 1), endPos));
}

template<typename TScoreValue, typename TScore, typename TInfixA, typename TInfixB, typename TDir, typename TSize, typename TEps, typename TAlign>
bool
_extendAndExtract(Align<TInfixB> const & localAlign,
				  TScoreValue scoreDropOff,
				  TScore const & scoreMatrix,
				  TInfixA const & a,
				  TInfixB const & b,
				  TDir direction,
				  TSize minLength,
				  TEps eps,
				  TAlign & align) {
SEQAN_CHECKPOINT
    typedef Seed<int, SimpleSeed> TSeed;
    typedef typename Position<TInfixB>::Type TPos;

	integrateAlign(align, localAlign);

	// begin and end position of local alignment (seed)
	TPos seedBeginA = sourceBeginPosition(row(localAlign, 0)) + beginPosition(a);
	TPos seedBeginB = sourceBeginPosition(row(localAlign, 1)) + beginPosition(b);
	TPos seedEndA = sourceEndPosition(row(localAlign, 0)) + beginPosition(a);
	TPos seedEndB = sourceEndPosition(row(localAlign, 1)) + beginPosition(b);

	if (direction == 3) {
		// set begin and end positions of align
		setSourceBeginPosition(row(align, 0), seedBeginA);
		setSourceBeginPosition(row(align, 1), seedBeginB);
		setBeginPosition(row(align, 0), 0);
		setBeginPosition(row(align, 1), 0);
		setSourceEndPosition(row(align, 0), seedEndA);
		setSourceEndPosition(row(align, 1), seedEndB);
	} else {
		// gapped X-drop extension of seed alignments
		TSeed seed(seedBeginA, seedBeginB, seedEndA - 1, seedEndB - 1);
		extendSeed(seed, scoreDropOff, scoreMatrix, host(a), host(b), direction, GappedXDrop());

		if (length(seed) < minLength - (int)floor(minLength*eps))
			return false;

		// set extended begin and end positions of align
		setSourceBeginPosition(row(align, 0), leftPosition(seed, 0));
		setSourceBeginPosition(row(align, 1), leftPosition(seed, 1));
		setBeginPosition(row(align, 0), 0);
		setBeginPosition(row(align, 1), 0);
		setSourceEndPosition(row(align, 0), rightPosition(seed, 0)+1);
		setSourceEndPosition(row(align, 1), rightPosition(seed, 1)+1);

		// banded alignment on ...
		__int64 startDiag = leftPosition(seed, 1) - leftPosition(seed, 0);
		if (direction != 1) { // ... extension to the left
			StringSet<TInfixB> str;
			appendValue(str, infix(host(a), leftPosition(seed, 0), seedBeginA));
			appendValue(str, infix(host(b), leftPosition(seed, 1), seedBeginB));
			_bandedInfixAlignment(str, seed, scoreMatrix, startDiag, align);
		}
		if (direction != 0) { // ... extension to the right
			StringSet<TInfixB> str;
			appendValue(str, infix(host(a), seedEndA, rightPosition(seed, 0)+1));
			appendValue(str, infix(host(b), seedEndB, rightPosition(seed, 1)+1));
			_bandedInfixAlignment(str, seed, scoreMatrix, startDiag, align);
		}
	}

	if ((TSize)length(row(align, 0)) < minLength)
		return false;

	// cut ends to obtain longest epsilon-match that contains the seed alignment
	TPos extBegin = beginPosition(source(row(align, 0))) + sourceBeginPosition(row(align, 0));
	longestEpsMatch(align, seedBeginA - extBegin, seedEndA - extBegin, minLength, eps);


    if ((TSize)length(row(align, 0)) < minLength)
		return false;

	return true;
}

//template<typename TRow, typename TSize>
//inline bool
//_isUpstream(TRow & row1, TRow & row2, TSize minLength) {
//SEQAN_CHECKPOINT
//    typedef typename Position<TRow>::Type TPos;
//    
//    TPos end1 = endPosition(sourceSegment(row1));
//    TPos begin2 = beginPosition(sourceSegment(row2));
//    if (end1 <= begin2) return true;
//    
//    TPos begin1 = beginPosition(sourceSegment(row1));
//    if (begin1 < begin2 && (begin2 - begin1 >= (TPos)minLength)) {
//        TPos end2 = endPosition(sourceSegment(row2));
//        if ((end1 < end2) && (end2 - end1 >= (TPos)minLength)) return true;
//    }
//    
//    return false;
//}

//template<TRow>
//inline bool
//_isDownstream(TRow & row1, TRow & row2) {
//QAN_CHECKPOINT
//    _isUpstream(row2, row1);
//}

//template<typename TSource, typename TId, typename TSize>
//inline void
//_insertMatch(String<SwiftLocalMatch<TSource, TId> >/*String<Align<TSource> >*/ & matches, SwiftLocalMatch<TSource, TId> match, /*Align<TSource> & align,*/ TSize minLength) {
//SEQAN_CHECKPOINT
//	typedef SwiftLocalMatch<TSource, TId> TMatch;
//    typedef typename Position<String<typename TMatch::TAlign/*Align<TSource> */> >::Type TPosString;
//    //typedef typename Row<Align<TSource> >::Type TRow;
//
//    TPosString maxPos = length(matches);
//    TPosString minPos = 0;
//
//    if(maxPos == 0) {
//        // matches string is empty
//        appendValue(matches, match);
//        return;
//    }
//
//    //TRow matchRow0 = row(match.align, 0);
//    //TRow matchRow1 = row(match.align, 1);
//
//    // determine insertion position
//    TPosString pos = maxPos / 2;
//    while (pos < maxPos) {
//        //TRow row0 = row(value(matches, pos).align, 0);
//        if (_isUpstream(value(matches, pos), match, 0, minLength)) {
//        //if (_isUpstream(row0, matchRow0, minLength)) {
//            // match is downstream of matches[pos] in first row
//            if (minPos == pos) {
//                ++pos;
//                break;
//            }
//            minPos = pos;
//        } else if (_isUpstream(match, value(matches, pos), 0, minLength)) {
//        //} else if (_isUpstream(matchRow0, row0, minLength)) {
//            // match is upstream of matches[pos] in first row
//            maxPos = pos;
//        } else {
//            // match overlaps in first row with another match
//            //TRow row1 = row(value(matches, pos).align, 1);
//            if (_isUpstream(value(matches, pos), match, 1, minLength)) {
//            //if (_isUpstream(row1, matchRow1, minLength)) {
//                // match is downstream of matches[pos] in second row
//                if (minPos == pos) {
//                    ++pos;
//                    break;
//                }
//                minPos = pos;
//            } else if (_isUpstream(match, value(matches, pos), 1, minLength)) {
//            //} else if (_isUpstream(matchRow1, row1, minLength)) {
//                // match is upstream of matches[pos] in second row
//                maxPos = pos;
//            } else {
//                // match overlaps in both rows with another match -> keep longer match
//                if (length(row(match.align, 0)/*matchRow0*/) > length(row(value(matches, pos).align, 0)/*row0*/)) {
//                    // new match is longer
//                    replace(matches, pos, pos+1, String<TMatch/*Align<TSource>*/ >());
//                    --maxPos;
//                } else {
//                    return;
//                }
//            }
//        }
//        pos = minPos + ((maxPos-minPos) / 2);
//    }
//    // insert match in matches at position pos
//    String<TMatch/*Align<TSource>*/ > str;
//    appendValue(str, match);
//    replace(matches, pos, pos, str);
//}

// checks if two matches overlap also in seq2 and 
// whether the non-overlaping parts are shorter than minLength
template<typename TMatch, typename TSize>
bool
checkOverlap(TMatch & matchA, TMatch & matchB, TSize minLength) {
SEQAN_CHECKPOINT
	// check overlap in seq2
	if (((TSize)matchA.begin2 - (TSize)matchB.begin2 > minLength &&
		 (TSize)matchA.end2 - (TSize)matchB.end2 > minLength) ||
		((TSize)matchB.begin2 - (TSize)matchA.begin2 > minLength &&
		 (TSize)matchB.end2 - (TSize)matchA.end2 > minLength)) {
		// non-overlapping parts of both matches are longer than minLength
		return false;
	}

	// same offset in both sequences?
	if ((TSize)matchA.begin1 - (TSize)matchB.begin1 != (TSize)matchA.begin2 - (TSize)matchB.begin2) {
		return false;
	}

	// check length of non-overlapping parts in seq1
	if (((TSize)matchA.begin1 - (TSize)matchB.begin1 > minLength &&
		 (TSize)matchA.end1 - (TSize)matchB.end1 > minLength) ||
		((TSize)matchB.begin1 - (TSize)matchA.begin1 > minLength &&
		 (TSize)matchB.end1 - (TSize)matchA.end1 > minLength)) {
		// non-overlapping parts of both matches are longer than minLength
		return false;
	}

	return true;
}

template<typename TSequence, typename TId, typename TSize>
void
maskOverlaps(String<SwiftLocalMatch<TSequence, TId> > & matches,
			 TSize minLength) {
SEQAN_CHECKPOINT
	typedef SwiftLocalMatch<TSequence, TId>						TMatch;
	typedef typename TMatch::TPos								TPos;
	
	typedef unsigned											TCargo;
	typedef IntervalAndCargo<TPos, TCargo>						TInterval;
	typedef IntervalTree<TPos, TCargo>							TIntervalTree;

	typedef typename Iterator<String<TCargo>, Standard>::Type	TIntervalIterator;
	typedef typename Iterator<String<TMatch>, Rooted>::Type		TIterator;

	// put matches into interval tree
	String<TInterval> intervals;
	for (unsigned i = 0; i < length(matches); ++i) {
		TMatch m = value(matches, i);
		appendValue(intervals, TInterval(m.begin1, m.end1, i));
	}
	TIntervalTree tree(intervals);

	TIterator it = begin(matches, Rooted());
	TIterator itEnd = end(matches, Rooted());

	for (; it != itEnd; ++it) {
		if ((*it).id == TMatch::INVALID_ID) continue;
		
		// find overlapping matches (overlap in seq1)
		String<TCargo> overlaps;
		findIntervals(tree, (*it).begin1, (*it).end1, overlaps);

		TIntervalIterator overlapIt = begin(overlaps, Standard());
		TIntervalIterator overlapItEnd = end(overlaps, Standard());

		for (; overlapIt != overlapItEnd; ++overlapIt) {
			TMatch overlap = value(matches, *overlapIt);
			if (overlap.id == TMatch::INVALID_ID || position(it) == *overlapIt) continue;

			// check for overlap in both sequences and small non-overlapping parts
			if (checkOverlap(*it, overlap, minLength)) {
				// set shorter match invalid
				if (length((*it)) > length(overlap)) {
					overlap.id = TMatch::INVALID_ID;
				} else {
					(*it).id = TMatch::INVALID_ID;
					break;
				}
			}
		}
	}
}

template<typename TSequence, typename TId, typename TSize>
void
compactMatches(String<SwiftLocalMatch<TSequence, TId> > & matches, TSize numMatches) {
	typedef SwiftLocalMatch<TSequence, TId>						TMatch;
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

template<typename TSource, typename TId, typename TSize, typename TSize1>
inline void
_insertMatch(String<SwiftLocalMatch<TSource, TId> > & matches,
			 SwiftLocalMatch<TSource, TId> match,
			 TSize minLength,
			 TSize1 & compactThresh,
			 TSize1 numMatches) {
SEQAN_CHECKPOINT

	appendValue(matches, match);

	if (length(matches) > compactThresh) {
		maskOverlaps(matches, minLength);		// remove overlaps and duplicates
		compactMatches(matches, numMatches);	// keep only the <numMatches> longest matches

		// raise compact threshold if many matches are kept
		if ((length(matches) << 1) > compactThresh)
			compactThresh += (compactThresh >> 1);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// banded alignment on swift hit and extraction of longest contained eps-match
template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TDelta, 
         typename TDrop, typename TSize1, typename TId, typename TSequence/*, typename TAlign*/>
void
verifySwiftHit(TInfixA const & a,
			   TInfixB const & b,
			   TEpsilon eps,
			   TSize minLength,
			   TDelta delta,
			   TDrop /*xDrop*/,
			   TSize1 & compactThresh,
			   TSize1 numMatches,
			   TId & queryId,
			   //String<TAlign> & matches,
			   String<SwiftLocalMatch<TSequence, TId> > & matches,
			   BandedGlobal) {
SEQAN_CHECKPOINT
	typedef typename SwiftLocalMatch<TSequence, TId>::TAlign TAlign;

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
	setSourceBeginPosition(row(align, 0), beginPosition(a) + sourceBeginPosition(row(bandedAlign, 0)));
	setSourceBeginPosition(row(align, 1), beginPosition(b) + sourceBeginPosition(row(bandedAlign, 1)));
	setBeginPosition(row(align, 0), 0);
	setBeginPosition(row(align, 1), 0);
	setSourceEndPosition(row(align, 0), beginPosition(a) + sourceEndPosition(row(bandedAlign, 0)));
	setSourceEndPosition(row(align, 1), beginPosition(b) + sourceEndPosition(row(bandedAlign, 1)));

	if ((TSize)length(row(align, 0)) < minLength)
		return;

	// insert eps-match in matches string
	SwiftLocalMatch<TSequence, TId> m(align, queryId);
	_insertMatch(matches, m, minLength, compactThresh, numMatches);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// banded alignment on swift hit and extraction of longest contained eps-match
template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TDelta, 
         typename TDrop, typename TSize1, typename TId, typename TSequence/*, typename TAlign*/>
void
verifySwiftHit(TInfixA const & a,
			   TInfixB const & b,
			   TEpsilon eps,
			   TSize minLength,
			   TDelta delta,
			   TDrop xDrop,
			   TSize1 & compactThresh,
			   TSize1 numMatches,
			   TId & queryId,
			   //String<TAlign> & matches,
			   String<SwiftLocalMatch<TSequence, TId> > & matches,
			   BandedGlobalExtend) {
SEQAN_CHECKPOINT
	typedef typename SwiftLocalMatch<TSequence, TId>::TAlign TAlign;

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
	if (!_extendAndExtract(bandedAlign, scoreDropOff, scoreMatrix, a, b, 2, minLength, eps, align))
		return;

	// insert eps-match in matches string
	SwiftLocalMatch<TSequence, TId> m(align, queryId);
	_insertMatch(matches, m, minLength, compactThresh, numMatches);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TDelta, typename TDrop,
         typename TSize1, typename TId, typename TSequence/*, typename TAlign*/, typename TTag>
void
verifySwiftHit(TInfixA const & a,
               TInfixB const & b,
               TEpsilon eps,
               TSize minLength,
               TDelta delta,
               TDrop xDrop,
			   TSize1 & compactThresh,
			   TSize1 numMatches,
			   TId & queryId,
			   //String<TAlign> & matches,
			   String<SwiftLocalMatch<TSequence, TId> > & matches,
			   TTag tag) {
SEQAN_CHECKPOINT
	typedef typename SwiftLocalMatch<TSequence, TId>::TAlign TAlign;

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
            // determine extension direction
            char direction;
            if (length(seedAlignments) == 1) direction = 2;
            else if (aliIt == begin(seedAlignments)) direction = 0;
            else if (aliIt == end(seedAlignments)-1) direction = 1;
            else direction = 3;

			// create alignment object for the complete sequences
			TAlign align;
			resize(rows(align), 2);
			setSource(row(align, 0), host(a));
			setSource(row(align, 1), host(b));

			// extend alignment and obtain longest contained eps-match
			if (!_extendAndExtract(*aliIt, scoreDropOff, scoreMatrix, a, b, direction, minLength, eps, align)) {
				aliIt++;
				continue;
			}

            // insert eps-match in matches string
			SwiftLocalMatch<TSequence, TId> m(align, queryId);
            _insertMatch(matches, m, minLength, compactThresh, numMatches);
            ++aliIt;
        }
		if (_verifyFast(tag)) break;
    }
}

// calls swift, verifies swift hits, returns eps-matches
template<typename TText, typename TIndex, typename TSize, typename TDrop, typename TSize1,
         typename TSource, typename TId, typename TTag>
int localSwift(Finder<TText, Swift<SwiftLocal> > & finder,
                Pattern<TIndex, Swift<SwiftLocal> > & pattern,
                double epsilon,
                TSize minLength,
                TDrop xDrop,
			    TSize1 & compactThresh,
			    TSize1 numMatches,
				StringSet<TId> & queryIDs,
                StringSet<String<SwiftLocalMatch<TSource, TId> > > & matches,
                //StringSet<String<Align<TSource> > > & matches,
				TTag tag) {
SEQAN_CHECKPOINT
	typedef SwiftLocalMatch<TSource, TId> TMatch;
    resize(matches, countSequences(needle(pattern)));
    TSize numSwiftHits = 0;

    TSize maxLength = 0;
    TSize totalLength = 0;

	while (find(finder, pattern, epsilon, minLength)) {
        ++numSwiftHits;

        // verification
		verifySwiftHit(infix(finder), infix(pattern, value(host(needle(pattern)), pattern.curSeqNo)), epsilon,
					   minLength, pattern.bucketParams[0].delta + pattern.bucketParams[0].overlap, xDrop, compactThresh,
					   numMatches, value(queryIDs, pattern.curSeqNo), value(matches, pattern.curSeqNo), tag);

        totalLength += length(infix(finder));
        if ((TSize)length(infix(finder)) > maxLength) maxLength = length(infix(finder));
	}

	//if (numSwiftHits > 0) {
	//	std::cout << "    Longest hit      : " << maxLength << std::endl;
	//	std::cout << "    Avg hit length   : " << totalLength/numSwiftHits << std::endl;
	//}
	
	typedef typename Iterator<StringSet<String<TMatch> >, Standard>::Type TIterator;
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());

	for(; it < itEnd; ++it) {
		maskOverlaps(*it, minLength);		// remove overlaps and duplicates
		compactMatches(*it, numMatches);	// keep only the <numMatches> longest matches
	}
    
    return numSwiftHits;
}
