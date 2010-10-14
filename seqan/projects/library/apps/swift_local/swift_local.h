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
#include <seqan/seeds2.h>
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
                // append new sub-alignment
                TPos begin = value(queue, 1).i1;
                TPos end = value(queue, 1).i2;

                TAlign ali(align);
                setClippedBeginPosition(row(ali, 0), toSourcePosition(row(ali, 0), begin));
                setClippedBeginPosition(row(ali, 1), toSourcePosition(row(ali, 1), begin));
				setBeginPosition(row(ali, 0), 0);
				setBeginPosition(row(ali, 1), 0);
                setClippedEndPosition(row(ali, 0), toSourcePosition(row(ali, 0), end));
                setClippedEndPosition(row(ali, 1), toSourcePosition(row(ali, 1), end));
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
    TPos len = right.i1 - left.i2;

    // check error rate
    return errors/(TFloat)(len) <= eps;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align from possEndsLeft and possEndsRight and sets the view positions of
// align to start and end position of the longest epsilon match
template<typename TPos, typename TCoord, typename TSize, typename TEps>
Pair<typename Iterator<String<Triple<TPos, TPos, TCoord> > >::Type>
longestEpsMatch(String<Triple<TPos, TPos, TCoord> > const & possEndsLeft,
				String<Triple<TPos, TPos, TCoord> > const & possEndsRight,
				TPos const alignLen,
				TPos const alignErr,
				TSize const matchMinLength,
				TEps const epsilon) {
SEQAN_CHECKPOINT
	typedef Triple<TPos, TPos, TCoord> TPossibleEnd;
	typedef typename Iterator<String<TPossibleEnd> >::Type	TIterator;

    // Identify longest eps match by iterating over combinations of left and right positions
    TIterator rightIt = end(possEndsRight) - 1;
    TIterator leftIt = end(possEndsLeft) - 1;
	TIterator right, left;

    TSize minLength = matchMinLength - 1;
	bool found = false;

	TSize totalLen = (*leftIt).i1 + alignLen + (*rightIt).i1;
	while (leftIt >= begin(possEndsLeft)) {
		totalLen = (*leftIt).i1 + alignLen + (*rightIt).i1;
		if (totalLen < minLength) break;
		while (rightIt >= begin(possEndsRight)) {
			totalLen = (*leftIt).i1 + alignLen + (*rightIt).i1;
			if (totalLen < minLength) break;
			TSize totalErr = (*leftIt).i2 + alignErr + (*rightIt).i2;
			if ((TEps)totalErr/(TEps)totalLen <= epsilon) {
				right = rightIt;
				left = leftIt;
				//std::cout << totalLen << std::endl;
				minLength = totalLen;
				found = true;
				break;
			}
			--rightIt;
		}
		rightIt = end(possEndsRight) - 1;
		--leftIt;
	}

	if (found) return Pair<TIterator>(left, right);
	else return Pair<TIterator>(0,0);
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

	TPosition seedBeginView = (TPosition)toViewPosition(row0, seedBegin + clippedBeginPosition(row0));
	TPosition seedEndView = (TPosition)toViewPosition(row0, seedEnd + clippedBeginPosition(row0));
	if (clippedBeginPosition(row0) > 0) {
		seedBeginView -= (TPosition)toViewPosition(row0, clippedBeginPosition(row0));
		seedEndView -= (TPosition)toViewPosition(row0, clippedBeginPosition(row0));
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
	setClippedBeginPosition(row(align, 0), toSourcePosition(row(align, 0), beginPos));
	setClippedBeginPosition(row(align, 1), toSourcePosition(row(align, 1), beginPos));
	setBeginPosition(row(align, 0), beginPos);
	setBeginPosition(row(align, 1), beginPos);
	setClippedEndPosition(row(align, 0), toSourcePosition(row(align, 0), endPos));
	setClippedEndPosition(row(align, 1), toSourcePosition(row(align, 1), endPos));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align and sets the view positions of
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
	setClippedBeginPosition(row(align, 0), toSourcePosition(row(align, 0), beginPos));
	setClippedBeginPosition(row(align, 1), toSourcePosition(row(align, 1), beginPos));
	setBeginPosition(row(align, 0), beginPos);
	setBeginPosition(row(align, 1), beginPos);
	setClippedEndPosition(row(align, 0), toSourcePosition(row(align, 0), endPos));
	setClippedEndPosition(row(align, 1), toSourcePosition(row(align, 1), endPos));
}

template<typename TPos, typename TCoord>
void
_possibleEnds(String<Triple<TPos, TPos, TCoord> > & possEnds,
			  String<Pair<TPos, TCoord> > const & bestEnds) {
SEQAN_CHECKPOINT
	typedef Pair<TPos, TCoord>	TBestEnd;
	typedef Triple<TPos, TPos, TCoord> TPossibleEnd;
	typedef typename Iterator<String<TBestEnd> >::Type	TBestEndIterator;

	TBestEndIterator lengthsIt = begin(bestEnds);
	TBestEndIterator lengthsEnd = end(bestEnds);
	TPos len = 0;

	/*for (unsigned i = 0; i < length(bestEnds); ++i) {
		std::cout << i << ": " << bestEnds[i].i1;
		std::cout << " (" << bestEnds[i].i2.i1 << "," << bestEnds[i].i2.i2 << ")" << std::endl;
	}
	std::cout << std::endl;*/

	if (lengthsIt < lengthsEnd) {
		bool match = false;
		if ((*lengthsIt).i1 == 0) match = true;

		TBestEndIterator prevLength = begin(bestEnds);
		++lengthsIt;
		++len;

		for (; lengthsIt < lengthsEnd; ++lengthsIt, ++prevLength, ++len) {
			if ((*lengthsIt).i1 == (*prevLength).i1 && (*lengthsIt).i1 <= length(bestEnds)) match = true;
			else {
				if (match) {
					appendValue(possEnds, TPossibleEnd(len-1, (*prevLength).i1, (*prevLength).i2));
				}
				match = false;
			}
		}

		if (match) {
			appendValue(possEnds, TPossibleEnd(len-1, (*prevLength).i1, (*prevLength).i2));
		}
	}

	/*for (unsigned i = 0; i < length(possEnds); ++i) {
		std::cout << possEnds[i].i1 << ": " << possEnds[i].i2;
		std::cout << " (" << possEnds[i].i3.i1 << "," << possEnds[i].i3.i2 << ")" << std::endl;
	}*/
}

template <typename TTrace, typename TPos, typename TCoord, typename TStringSet, typename TScore, typename TDiagonal>
inline void
_align_banded_nw_best_ends(TTrace& trace,
						   String<Pair<TPos, TCoord> > & bestEnds,
						   TStringSet const& str,
						   TScore const & sc,
						   TDiagonal const diagL,
						   TDiagonal const diagU)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TTrace>::Type TSize;

	SEQAN_ASSERT_GEQ(diagU, diagL);

	// Initialization
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1) + 1;
	TSize len2 = length(str2) + 1;
	TSize diagonalWidth = (TSize) (diagU - diagL + 1);
	TSize hi_diag = diagonalWidth;
	TSize lo_diag = 0;
	if (diagL > 0) lo_diag = 0;
	else lo_diag = (diagU < 0) ? hi_diag : (TSize) (1 - diagL); 
	TSize lo_row = (diagU <= 0) ? -diagU : 0;
	TSize hi_row = len2;
	if (len1 - diagL < hi_row) hi_row = len1 - diagL;
	TSize height = hi_row - lo_row;

	typedef String<TScoreValue> TRow;
	TRow mat, len;
	resize(mat, diagonalWidth);
	resize(len, diagonalWidth);
	resize(trace, height * diagonalWidth);
	fill(bestEnds, _max(len1+len2-2, diagonalWidth), Pair<TPos, TCoord>(_max(len1+len2-1, diagonalWidth+1), TCoord()));
	
	//// Debug stuff
	//String<TScoreValue> originalMat;
	//resize(originalMat, len1 * len2);
	//TSize count = 0;
	//for(TSize i=0; i<len2; ++i) {
	//	for(TSize j=0; j<len1; ++j) {
	//		value(originalMat, i * len1 + j) = count;
	//		std::cerr << count << ',';
	//		++count;
	//	}
	//	std::cerr << std::endl;
	//}
	//std::cerr << std::endl;

	// Classical DP with affine gap costs
	typedef typename Iterator<TRow, Standard>::Type TRowIter;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	TSize actualCol = 0;
	TSize actualRow = 0;
	TScoreValue verti_val = 0;
	TScoreValue hori_val = 0;
	TScoreValue hori_len = len1+len2+1;

	for(TSize row = 0; row < height; ++row) {
		actualRow = row + lo_row;
		if (lo_diag > 0) --lo_diag;
		if ((TDiagonal)actualRow >= (TDiagonal)len1 - diagU) --hi_diag;
		TTraceIter traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
		TRowIter matIt = begin(mat, Standard()) + lo_diag;
		TRowIter lenIt = begin(len, Standard()) + lo_diag;
		hori_val = InfimumValue<TScoreValue>::VALUE;
		hori_len = len1+len2+1;
		for(TSize col = lo_diag; col<hi_diag; ++col, ++matIt, ++traceIt, ++lenIt) {
			actualCol = col + diagL + actualRow;
			if (actualCol >= len1) break;
			//std::cerr << row << ',' << col << ':' << value(originalMat, actualRow * len1 + actualCol) << std::endl;

			if ((actualRow != 0) && (actualCol != 0)) {
				// Get the new maximum for mat
				*matIt += score(const_cast<TScore&>(sc), ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
				*traceIt = Diagonal;
				++(*lenIt);
				if ((verti_val = (col < diagonalWidth - 1) ? *(matIt+1) + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : InfimumValue<TScoreValue>::VALUE) > *matIt) {
					*matIt = verti_val;
					*traceIt = Vertical;
					*lenIt = *(lenIt+1) + 1;
				}						
				if ((hori_val = (col > 0) ? hori_val + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2) : InfimumValue<TScoreValue>::VALUE) > *matIt) {
					*matIt = hori_val;
					*traceIt = Horizontal;
					*lenIt = hori_len + 1;
				}
				hori_val = *matIt;
				hori_len = *lenIt;
			} else {			
				// Usual initialization for first row and column
				if (actualRow == 0) {
					*matIt = actualCol * scoreGapExtendHorizontal(sc, ((int) actualCol - 1), -1, str1, str2);
					*lenIt = actualCol;
				}
				else {
					*matIt = actualRow * scoreGapExtendVertical(sc, -1, ((int) actualRow - 1), str1, str2);
					*lenIt = actualRow;
					hori_val = *matIt;
					hori_len = actualRow;
				}
			}
			TPos errors = (*matIt - (*lenIt * scoreMatch(const_cast<TScore&>(sc)))) /
						(scoreGap(const_cast<TScore&>(sc)) - scoreMatch(const_cast<TScore&>(sc)));
			if (*lenIt >= 0 && errors < value(bestEnds, *lenIt).i1)
				value(bestEnds, *lenIt) = Pair<TPos, TCoord>(errors, TCoord(row, col));
			//std::cerr << row << ',' << col << ':' << *matIt << std::endl;
		}
	}
}

template<typename TInfixA, typename TInfixB, typename TSeed>
void
_reverseLeftExtension(TInfixA const & a,
					  TInfixB const & b,
					  TSeed & seed,
					  TSeed & seedOld) {
SEQAN_CHECKPOINT
	TInfixB infixA = infix(host(a), getBeginDim0(seed), getBeginDim0(seedOld));
	TInfixB infixB = infix(host(b), getBeginDim1(seed), getBeginDim1(seedOld));
	reverseInPlace(infixA);
	reverseInPlace(infixB);
}

template<typename TMatrix, typename TPossEnd, typename TInfixA, typename TInfixB, typename TSeed, typename TScore>
void
_fillMatrixBestEndsLeft(TMatrix & matrixLeft,
							String<TPossEnd> & possibleEndsLeft,
							TInfixA const & a,
							TInfixB const & b,
							TSeed & seed,
							TSeed & seedOld,
							TScore const & scoreMatrix) {
SEQAN_CHECKPOINT
	typedef StringSet<TInfixB>		TInfixSet;
	typedef typename TPossEnd::T1	TPos;
	typedef typename TPossEnd::T3	TCoord;
	typedef Pair<TPos, TCoord>		TBestEnd; // number of errors and coordinate in alignment matrix for a given trace length

	TInfixB infixA = infix(host(a), getBeginDim0(seed), getBeginDim0(seedOld));
	TInfixB infixB = infix(host(b), getBeginDim1(seed), getBeginDim1(seedOld));

	reverseInPlace(infixA);
	reverseInPlace(infixB);

	TInfixSet str;
	appendValue(str, infixA);
	appendValue(str, infixB);

	String<TBestEnd> bestEnds;
	_align_banded_nw_best_ends(matrixLeft, bestEnds, str, scoreMatrix, 
							   getUpperDiagonal(seedOld) - getUpperDiagonal(seed),
							   getUpperDiagonal(seedOld) - getLowerDiagonal(seed));
	_possibleEnds(possibleEndsLeft, bestEnds);
}

template<typename TMatrix, typename TPossEnd, typename TInfixA, typename TInfixB, typename TSeed, typename TScore>
void
_fillMatrixBestEndsRight(TMatrix & matrixRight,
							String<TPossEnd> & possibleEndsRight,
							TInfixA const & a,
							TInfixB const & b,
							TSeed & seed,
							TSeed & seedOld,
							TScore const & scoreMatrix) {
SEQAN_CHECKPOINT
	typedef StringSet<TInfixB>		TInfixSet;
	typedef typename TPossEnd::T1	TPos;
	typedef typename TPossEnd::T3	TCoord;
	typedef Pair<TPos, TCoord>		TBestEnd; // number of errors and coordinate in alignment matrix for a given trace length

	TInfixSet str;
	appendValue(str, infix(host(a), getEndDim0(seedOld), getEndDim0(seed)));
	appendValue(str, infix(host(b), getEndDim1(seedOld), getEndDim1(seed)));

	String<TBestEnd> bestEnds;
	_align_banded_nw_best_ends(matrixRight, bestEnds, str, scoreMatrix, 
							   getLowerDiagonal(seedOld) - getUpperDiagonal(seed),
							   getLowerDiagonal(seedOld) - getLowerDiagonal(seed)); 
	_possibleEnds(possibleEndsRight, bestEnds);
}

template<typename TMatrix, typename TCoord, typename TInfixA, typename TInfixB, typename TSeed, typename TPos, typename TAlign>
void
_tracebackLeft(TMatrix const & matrixLeft,
			   TCoord const & coordinate,
			   TInfixA const & a,
			   TInfixB const & b,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const endLeftA,
			   TPos const endLeftB,
			   TAlign & align) {
SEQAN_CHECKPOINT
	typedef StringSet<TInfixB>							TInfixSet;
	typedef typename Size<TInfixSet>::Type				TSize;
	typedef typename Iterator<String<TraceBack> >::Type	TIterator;

	TInfixSet str;
	TInfixB infixA = infix(host(a), getBeginDim0(seed), getBeginDim0(seedOld));
	TInfixB infixB = infix(host(b), getBeginDim1(seed), getBeginDim1(seedOld));
	appendValue(str, infixA);
	appendValue(str, infixB);

	bool overallMaxValue[2]; // only needed for standard traceback function
	overallMaxValue[0] = 1; overallMaxValue[1] = 0;
	TPos overallMaxIndex[4];
	overallMaxIndex[0] = coordinate.i1;
	overallMaxIndex[1] = coordinate.i2;

	_Align_Traceback<TPos> traceBack;
	_align_banded_nw_trace(traceBack, str, matrixLeft, overallMaxValue, overallMaxIndex,
				   getUpperDiagonal(seedOld) - getUpperDiagonal(seed), getUpperDiagonal(seedOld) - getLowerDiagonal(seed));
	
	reverseInPlace(traceBack.sizes);
	reverseInPlace(traceBack.tvs);

	TIterator tvsIt = end(traceBack.tvs) - 1;
	TIterator tvsBegin = begin(traceBack.tvs);
	TSize newLen = length(traceBack.tvs);
	while (tvsIt >= tvsBegin && (*tvsIt) != (TraceBack)0) {
		--newLen;
		--tvsIt;
	}
	resize(traceBack.tvs, newLen);
	resize(traceBack.sizes, newLen);

	Align<TInfixB> infixAlign;
	resize(rows(infixAlign), 2);
	assignSource(row(infixAlign, 0), infix(str[0], length(str[0]) - endLeftA, length(str[0])));
	assignSource(row(infixAlign, 1), infix(str[1], length(str[1]) - endLeftB, length(str[1])));

	_pump_trace_2_Align(infixAlign, traceBack);
	integrateAlign(align, infixAlign);
}


template<typename TMatrix, typename TCoord, typename TInfixA, typename TInfixB, typename TSeed, typename TPos, typename TAlign>
void
_tracebackRight(TMatrix const & matrixRight,
			   TCoord const & coordinate,
			   TInfixA const & a,
			   TInfixB const & b,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const endRightA,
			   TPos const endRightB,
			   TAlign & align) {
SEQAN_CHECKPOINT
	typedef StringSet<TInfixB>							TInfixSet;
	typedef typename Size<TInfixSet>::Type				TSize;
	typedef typename Iterator<String<TraceBack> >::Type	TIterator;
		
	TInfixSet str;
	appendValue(str, infix(host(a), getEndDim0(seedOld), getEndDim0(seed)));
	appendValue(str, infix(host(b), getEndDim1(seedOld), getEndDim1(seed)));

	bool overallMaxValue[2]; // only needed for standard traceback function
	overallMaxValue[0] = 1; overallMaxValue[1] = 0;
	TPos overallMaxIndex[4];
	overallMaxIndex[0] = coordinate.i1;
	overallMaxIndex[1] = coordinate.i2;

	_Align_Traceback<TPos> traceBack;
	_align_banded_nw_trace(traceBack, str, matrixRight, overallMaxValue, overallMaxIndex,
				   getLowerDiagonal(seedOld) - getUpperDiagonal(seed), getLowerDiagonal(seedOld) - getLowerDiagonal(seed));

	TSize skipLen = 0;
	TIterator tvsIt = begin(traceBack.tvs);
	TIterator tvsEnd = end(traceBack.tvs);
	while (tvsIt < tvsEnd && (*tvsIt) != (TraceBack)0) {
		++skipLen;
		++tvsIt;
	}
	replace(traceBack.tvs, 0, skipLen, String<TraceBack>());
	replace(traceBack.sizes, 0, skipLen, String<TPos>());

	Align<TInfixB> infixAlign;
	resize(rows(infixAlign), 2);
	assignSource(row(infixAlign, 0), infix(str[0], 0, endRightA));
	assignSource(row(infixAlign, 1), infix(str[1], 0, endRightB));

	_pump_trace_2_Align(infixAlign, traceBack);
	integrateAlign(align, infixAlign);
}

template<typename TInfixA, typename TInfixB, typename TSeed, typename TPos, typename TDir, typename TScore,
		 typename TSize, typename TEps, typename TAlign>
bool
_bestExtension(TInfixA const & a,
			   TInfixB const & b,
			   TSeed & seed,
			   TSeed & seedOld,
			   TPos const alignLen,
			   TPos const alignErr,
			   TScore const & scoreMatrix,
			   TDir const direction,
			   TSize const minLength,
			   TEps const eps,
			   TAlign & align) {
SEQAN_CHECKPOINT
	typedef String<TraceBack>				TAlignmentMatrix;
	typedef Pair<TPos>						TCoord; // coordinate in alignment matrix
	typedef Triple<TPos, TPos, TCoord>		TEndInfo;  // trace length, number of errors, coordinate in alignment matrix
	typedef typename Iterator<String<TEndInfo> >::Type	TEndIterator;
	typedef typename Value<TScore>::Type	TScoreValue;
	typedef StringSet<TInfixB>	TInfixSet;

	// variables for banded alignment, possible ends of match
	TAlignmentMatrix matrixRight, matrixLeft;
	String<TEndInfo> possibleEndsLeft, possibleEndsRight;

	// fill banded matrix and gaps string for ...
	if (direction == EXTEND_BOTH || direction == EXTEND_LEFT) { // ... extension to the left
		_fillMatrixBestEndsLeft(matrixLeft, possibleEndsLeft, a, b, seed, seedOld, scoreMatrix);
		// Caution: left extension infix is now reversed in host(a and b) !!!
	} else appendValue(possibleEndsLeft, TEndInfo(0, 0, TCoord(0, 0)));
	if (direction == EXTEND_BOTH || direction == EXTEND_RIGHT) { // ... extension to the right
		_fillMatrixBestEndsRight(matrixRight, possibleEndsRight, a, b, seed, seedOld, scoreMatrix);
	} else appendValue(possibleEndsRight, TEndInfo(0, 0, TCoord(0, 0)));

	// longest eps match on poss ends string
	Pair<TEndIterator> endPair = longestEpsMatch(possibleEndsLeft, possibleEndsRight, alignLen, alignErr, minLength, eps);

	if (endPair == Pair<TEndIterator>(0, 0)) { // no eps-match found
		if (direction != 1) 
			_reverseLeftExtension(a, b, seed, seedOld); // back to original orientation
		return false;
	}

	// end positions of maximal eps-match
	TPos endLeftA = 0, endLeftB = 0;
	TPos endRightA = 0, endRightB = 0;
	if((*endPair.i1).i1 != 0) { // extension to the left
		endLeftB = (*endPair.i1).i3.i1;
		if (getUpperDiagonal(seedOld) - getLowerDiagonal(seed) <= 0)
			endLeftB -= (TPos)(getUpperDiagonal(seedOld) - getLowerDiagonal(seed));
		endLeftA = (TPos)((*endPair.i1).i3.i2 + endLeftB + getUpperDiagonal(seedOld) - getUpperDiagonal(seed));
	}
	if((*endPair.i2).i1 != 0) { // ... extension to the right
		endRightB = (*endPair.i2).i3.i1;
		if (getLowerDiagonal(seedOld) - getLowerDiagonal(seed) <= 0)
			endRightB -= (TPos)(getLowerDiagonal(seedOld) - getLowerDiagonal(seed));
		endRightA = (TPos)((*endPair.i2).i3.i2 + endRightB + getLowerDiagonal(seedOld) - getUpperDiagonal(seed));
	}

	// set begin and end positions of align
	setClippedBeginPosition(row(align, 0), getBeginDim0(seedOld) - endLeftA);
	setClippedBeginPosition(row(align, 1), getBeginDim1(seedOld) - endLeftB);
	setBeginPosition(row(align, 0), 0);
	setBeginPosition(row(align, 1), 0);
	setClippedEndPosition(row(align, 0), getEndDim0(seedOld) + endRightA);
	setClippedEndPosition(row(align, 1), getEndDim1(seedOld) + endRightB);

	// traceback through matrix from begin/end pos on ...
	if((*endPair.i1).i1 != 0) { // ... extension to the left
		_tracebackLeft(matrixLeft, (*endPair.i1).i3, a, b, seed, seedOld, endLeftA, endLeftB, align);
	}
	if((*endPair.i2).i1 != 0) { // ... extension to the right
		_tracebackRight(matrixRight, (*endPair.i2).i3, a, b, seed, seedOld, endRightA, endRightB, align);
	}

	if (direction == EXTEND_BOTH || direction == EXTEND_LEFT) 
		_reverseLeftExtension(a, b, seed, seedOld); // back to original orientation

	return true;
}

template<typename TScoreValue, typename TScore, typename TInfixA, typename TInfixB, typename TSize, typename TEps, typename TAlign>
bool
_extendAndExtract(Align<TInfixB> const & localAlign,
				  TScoreValue scoreDropOff,
				  TScore const & scoreMatrix,
				  TInfixA const & a,
				  TInfixB const & b,
				  ExtensionDirection direction,
				  TSize minLength,
				  TEps eps,
				  TAlign & align) {
SEQAN_CHECKPOINT
    typedef typename Position<TInfixB>::Type TPos;
    typedef Seed<Simple> TSeed;

	integrateAlign(align, localAlign);

	// begin and end position of local alignment (seed)
	TPos seedBeginA = clippedBeginPosition(row(localAlign, 0)) + beginPosition(a);
	TPos seedBeginB = clippedBeginPosition(row(localAlign, 1)) + beginPosition(b);
	TPos seedEndA = clippedEndPosition(row(localAlign, 0)) + beginPosition(a);
	TPos seedEndB = clippedEndPosition(row(localAlign, 1)) + beginPosition(b);

	if (direction == EXTEND_NONE) {
		// set begin and end positions of align
		setClippedBeginPosition(row(align, 0), seedBeginA);
		setClippedBeginPosition(row(align, 1), seedBeginB);
		setBeginPosition(row(align, 0), 0);
		setBeginPosition(row(align, 1), 0);
		setClippedEndPosition(row(align, 0), seedEndA);
		setClippedEndPosition(row(align, 1), seedEndB);

		if ((TSize)length(row(align, 0)) < minLength)
			return false;

		/*TODO longest eps match?*/
	} else {
		// gapped X-drop extension of seed alignments
		TSeed seed(seedBeginA, seedBeginB, seedEndA, seedEndB);
		TSeed seedOld(seed);
		extendSeed(seed, host(a), host(b), direction, scoreMatrix, scoreDropOff, GappedXDrop());

		if (static_cast<__int64>(getSeedSize(seed)) < minLength - (int)floor(minLength*eps))
			return false;

		TPos alignLen = length(row(localAlign, 0)); // TODO not length(row0)!
		TPos alignErr = 0;
		for (TPos i = 0; i < alignLen; ++i) {
			if (!isMatch(localAlign, i)) ++alignErr;
		}

		if (!_bestExtension(a, b, seed, seedOld, alignLen, alignErr, scoreMatrix, direction, minLength, eps, align))
			return false;
	}
	return true;
}

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
	setClippedBeginPosition(row(align, 0), beginPosition(a) + clippedBeginPosition(row(bandedAlign, 0)));
	setClippedBeginPosition(row(align, 1), beginPosition(b) + clippedBeginPosition(row(bandedAlign, 1)));
	setBeginPosition(row(align, 0), 0);
	setBeginPosition(row(align, 1), 0);
	setClippedEndPosition(row(align, 0), beginPosition(a) + clippedEndPosition(row(bandedAlign, 0)));
	setClippedEndPosition(row(align, 1), beginPosition(b) + clippedEndPosition(row(bandedAlign, 1)));

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
	if (!_extendAndExtract(bandedAlign, scoreDropOff, scoreMatrix, a, b, EXTEND_BOTH, minLength, eps, align))
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

		//std::cout << beginPosition(infix(finder)) << ",";
		//std::cout << endPosition(infix(finder)) << "  ";
		//std::cout << beginPosition(infix(pattern, value(host(needle(pattern)), pattern.curSeqNo))) << ",";
		//std::cout << endPosition(infix(pattern, value(host(needle(pattern)), pattern.curSeqNo))) << std::endl;

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
		if (length(*it) > 0) {
			maskOverlaps(*it, minLength);		// remove overlaps and duplicates
			compactMatches(*it, numMatches);	// keep only the <numMatches> longest matches
		}
	}
    
    return numSwiftHits;
}
