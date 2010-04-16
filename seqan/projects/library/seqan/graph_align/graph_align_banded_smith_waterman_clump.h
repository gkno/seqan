 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_BANDED_SMITH_WATERMAN_CLUMP_H
#define SEQAN_HEADER_GRAPH_ALIGN_BANDED_SMITH_WATERMAN_CLUMP_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

//template<typename TString>
//inline void
//_initAlign(Align<TString>& align,
//           StringSet<TString, Dependent<> > const& str) {
//    resize(rows(align), 2);
//    assignSource(row(align, 0), str[0]);
//    assignSource(row(align, 1), str[1]);
//}

template<typename TString>
inline void
_initAlign(Graph<Alignment<StringSet<TString, Dependent<> > > >& align,
           StringSet<TString, Dependent<> > const& str) {
	assignStringSet(align, str);
}

//template<typename TString, typename TPos>
//inline void
//_finishAlign(Align<TString>& align,
//             TPos const begin1,
//             TPos const end1,
//             TPos const begin2,
//             TPos const end2) {
//    setSourceBeginPosition(row(align, 0), begin1);
//	setSourceBeginPosition(row(align, 1), begin2);
//	setBeginPosition(row(align, 0), 0);
//	setBeginPosition(row(align, 1), 0);
//	setSourceEndPosition(row(align, 0), end1);
//	setSourceEndPosition(row(align, 1), end2);
//}

template<typename TString, typename TPos>
inline void
_finishAlign(Graph<Alignment<StringSet<TString, Dependent<> > > >& g,
             TPos const begin1,
             TPos const end1,
             TPos const begin2,
             TPos const end2) {
    // Nothing to be done?
    //stringSet(g)[0] = infix(stringSet(g)[0], begin1, end1);
    //stringSet(g)[1] = infix(stringSet(g)[1], begin2, end2);
}

////////////////////////////////////////////////////////////////////////////////
//template <typename TString, typename TId, typename TPos, typename TTraceValue>
//inline void
//_align_trace_print(Align<TString>& align,
//				   StringSet<TString, Dependent<> > const&,
//				   TId,
//				   TPos const pos1,
//				   TId,
//				   TPos const pos2,
//				   TPos const segLen,
//				   TTraceValue const tv)
//{
//	SEQAN_CHECKPOINT
//
//	// TraceBack values
//	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
//
//	if (segLen == 0 || tv == Diagonal) return;
//
//    if (tv == Horizontal) {
//        insertGaps(row(align, 0), toViewPosition(row(align, 0), pos1), segLen);
//    } else if (tv == Vertical) {
//        insertGaps(row(align, 1), toViewPosition(row(align, 1), pos2), segLen);
//    }
//}

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TAlign, typename TTrace, typename TVal, typename TIndexPair, typename TDiagonal, typename TForbidden>
inline void
_align_banded_sw_trace(TStringSet const& str,
                       TAlign& align,
                       TTrace const& trace,
                       TVal const initialDir,
                       TIndexPair const& indexPair,
                       TDiagonal diagL,
                       TDiagonal diagU,
                       TForbidden& forbidden) {
    SEQAN_CHECKPOINT
    typedef typename Value<TStringSet>::Type TString;
    typedef typename Id<TStringSet>::Type TId;
    typedef typename Size<TTrace>::Type TSize;
    typedef typename Value<TTrace>::Type TTraceValue;

    // Traceback values
	TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

    // Initialization
    _initAlign(align, str);
    TString const& str1 = str[0];
    TString const& str2 = str[1];
    TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
    TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
    TSize len1 = length(str1);
    TSize len2 = length(str2);
    TSize lo_row = (diagU <= 0) ? -diagU : 0;
    TSize diagonalWidth = (TSize) (diagU - diagL + 1);

    // Start the trace from the cell with the max value
    TSize row = indexPair[0];
    TSize col = indexPair[1];
    TSize endRow = indexPair[0] + lo_row;
    TSize endCol = indexPair[1] + diagL + endRow;
	TSize actualRow = row + lo_row;
    TSize actualCol = col + diagL + actualRow;
    if ((actualCol == 0) || (actualRow == 0)) return;
	if (actualCol < len1) _align_trace_print(align, str, id1, actualCol, id2, actualRow, len1 - actualCol, Horizontal);
	if (actualRow < len2) _align_trace_print(align, str, id1, actualCol, id2, actualRow, len2 - actualRow, Vertical);
	_setForbiddenCell(forbidden, row+1, col+1, diagonalWidth);
	
    TTraceValue traceValue = initialDir;
    TTraceValue nextTraceValue;
    if (traceValue == Diagonal) nextTraceValue = trace[(--row) * diagonalWidth + col];
    else if (traceValue == Vertical) nextTraceValue = trace[(--row) * diagonalWidth + (++col)];
    else if (traceValue == Horizontal) nextTraceValue = trace[row * diagonalWidth + (--col)];
    actualRow = row + lo_row;
    actualCol = col + diagL + actualRow;
	if (nextTraceValue == Diagonal) _setForbiddenCell(forbidden, row+1, col+1, diagonalWidth);
    TSize segLen = 1;
	
    while (nextTraceValue != Stop) {
        if (traceValue == nextTraceValue) {
            ++segLen;
        } else {
            _align_trace_print(align, str, id1, actualCol, id2, actualRow, segLen, traceValue);
            segLen = 1;
        }
        traceValue = nextTraceValue;
        if (traceValue == Diagonal) nextTraceValue = trace[(--row) * diagonalWidth + col];
        else if (traceValue == Vertical) nextTraceValue = trace[(--row) * diagonalWidth + (++col)];
        else if (traceValue == Horizontal) nextTraceValue = trace[row * diagonalWidth + (--col)];
        actualRow = row + lo_row;
        actualCol = col + diagL + actualRow;
		if (nextTraceValue == Diagonal) _setForbiddenCell(forbidden, row+1, col+1, diagonalWidth);
    }
    if (segLen) _align_trace_print(align, str, id1, actualCol, id2, actualRow, segLen, traceValue);
    
	// Handle the remaining sequence
	if (actualCol != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) actualCol, Horizontal);
	if (actualRow != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) actualRow, Vertical);

    _finishAlign(align, actualCol, endCol, actualRow, endRow);
}

////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScore, typename TIndexPair, typename TDiagonal, typename TForbidden>
inline typename Value<TScore>::Type
_align_banded_sw(TTrace& trace,
                 TStringSet const& str,
                 TScore const& sc,
                 typename Value<TTrace>::Type& initialDir,
                 TIndexPair& indexPair,
                 TDiagonal diagL,
                 TDiagonal diagU,
                 TForbidden forbidden) {
    SEQAN_CHECKPOINT
    typedef typename Value<TTrace>::Type TTraceValue;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Value<TStringSet>::Type TString;
    typedef typename Size<TTrace>::Type TSize;

    // TraceBack values for Smith Waterman
    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2; TTraceValue Stop = 3;

    // Initialization
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

    typedef String<TScoreValue> TColumn;
    TColumn mat;
    fill(mat, diagonalWidth, 0);
    fill(trace, height * diagonalWidth, Stop);
    
    // Record the max score
    TScoreValue score_max = 0;
    indexPair[0] = 0; indexPair[1] = 0;
    initialDir = Stop;

    // Classical DP
    //TScoreValue max_val = 0;
    TSize actualCol, actualRow;
    TScoreValue verti_val, hori_val;

    // Initialize first row
    typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
    typedef typename Iterator<TColumn, Standard>::Type TMatIter;
    TTraceIter traceIt;
    TMatIter matIt;
    
    if (lo_diag > 0) --lo_diag;
    if (lo_row >= len1 - diagU) --hi_diag;

    for (TSize row = 1; row < height; ++row) {
        actualRow = row + lo_row;
        if (lo_diag > 0) --lo_diag;
        if (actualRow >= len1 - diagU) --hi_diag;
        traceIt = begin(trace, Standard()) + row * diagonalWidth + lo_diag;
        matIt = begin(mat, Standard()) + lo_diag;
        TScoreValue diagValMat = *matIt;
        hori_val = 0;
        
        for (TSize col = lo_diag; col < hi_diag; ++col, ++matIt, ++traceIt) {
            actualCol = col + diagL + actualRow;
            if (actualCol != 0) {
                if (_isClumping(forbidden, col+1, row+1, diagonalWidth)) {
                    //*traceIt = Stop; // TODO: Do we need this line?
                    *matIt = 0;
                } else {
                    // Get the new maximum for mat
                    *matIt += score(const_cast<TScore&>(sc), ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
                    *traceIt = Diagonal;
                }

                // Get the new maximum for vertical
                if (col < diagonalWidth - 1) {
                    verti_val = *(matIt+1) + scoreGapExtendVertical(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
                } else {
                    verti_val = -2;
                }
                if (verti_val > *matIt) {
                    *matIt = verti_val;
                    *traceIt = Vertical;
                }

                // Get the new maximum for horizontal
                if (col > 0) {
                    hori_val = hori_val + scoreGapExtendHorizontal(sc, ((int) actualCol - 1), ((int) actualRow - 1), str1, str2);
                } else {
                    hori_val = -2;
                }
                if (hori_val > *matIt) {
                    *matIt = hori_val;
                    *traceIt = Horizontal;
                }

                // Check if new maximum is greater than 0
                if (0 > *matIt) {
                    *matIt = 0;
                    *traceIt = Stop;
                    
                }

                // Record the new best score
                if (*matIt > score_max) {
                    indexPair[0] = row; indexPair[1] = col;
                    score_max = *matIt;
                    initialDir = *traceIt;
                }
            } else {
                // Init first col
                *matIt = 0;
                *traceIt = Stop;
            }
            hori_val = *matIt;
        }
    }
 //   // Debug code
 //   std::cerr << std::endl;
	//for(TSize i= 0; i<height;++i) {
	//	for(TSize j= 0; j<diagonalWidth;++j) {
	//		std::cerr << (TSize) getValue(trace, i*diagonalWidth + j) << ',';
	//	}
	//	std::cerr << std::endl;
	//}
	//std::cerr << "Max score: " << indexPair[0] << ',' << indexPair[1] << ':' << score_max << " (" << (TSize) initialDir << ")" << std::endl;
    return score_max;
}

////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TForbidden, typename TScore, typename TDiagonal>
inline typename Value<TScore>::Type
_localAlignment(TStringSet& str,
                TAlign& align,	
                TForbidden& forbidden,
                TScore const& sc,
                TDiagonal diag1,
                TDiagonal diag2,
                BandedSmithWatermanClump) {
    SEQAN_CHECKPOINT
    typedef typename Value<TScore>::Type TScoreValue;
    typedef typename Size<TStringSet>::Type TSize;
    typedef unsigned char TDir;

    // Maximum value
    TScoreValue maxScore;
    TSize indexPair[2];

    // Create the trace
    TDir initialDir;
    String<TDir> trace;
    maxScore = _align_banded_sw(trace, str, sc, initialDir, indexPair, diag1, diag2, forbidden);

    // Follow the trace and create the alignment
    _align_banded_sw_trace(str, align, trace, initialDir, indexPair, diag1, diag2, forbidden);

    return maxScore;
}

template<typename TString, typename TAlignments, typename TScores, typename TScoreValue, typename TSpec2, typename TDiagonal>
inline void
_localAlignment(StringSet<TString, Dependent<> > const& str,
				TAlignments& alignments,
				TScores& scores,
				Score<TScoreValue, TSpec2> const& sc,
				TScoreValue minScore,
				TDiagonal diag1,
				TDiagonal diag2,
				BandedSmithWatermanClump) {
	SEQAN_CHECKPOINT
	typedef typename Value<TAlignments>::Type TAlign;
	typedef typename Size<TString>::Type TSize;
  
	// For clumpping remember the used positions
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
    TSize diagonalWidth = (TSize) (diag2 - diag1 + 1);
    
    TSize lo_row = (diag2 < 0) ? -diag2 : 0;
    TSize hi_row = (len1 - diag1 < len2) ? len1 - diag1 : len2;
    TSize height = hi_row - lo_row;

	String<bool> forbidden;
	fill(forbidden, (height+1) * diagonalWidth, false);//(len1+1) * (len2+1), false);

	// Stop looking for local alignments, if their score is too low
	while (true) {
		// Create the local alignment
        TAlign align;
		TScoreValue local_score = _localAlignment(str, align, forbidden, sc, diag1, diag2, BandedSmithWatermanClump());

        //// Debug code
        //for (int i = 0; i < height; ++i) {
        //    for (int j = 0; j < diagonalWidth; ++j) {
        //        std::cerr << forbidden[i * diagonalWidth + j] << ", ";
        //    }
        //    std::cerr << std::endl;
        //}
        std::cerr << align << std::endl;

		if (local_score >= minScore) {
            appendValue(alignments, align);
            appendValue(scores, local_score);
		} else break;
	}
}

template<typename TString, typename TAlignments, typename TScores, typename TScoreValue, typename TSpec2, typename TDiagonal, typename TTag>
inline void
multiLocalAlignment(StringSet<TString, Dependent<> > const& str,
                    TAlignments& alignments,
                    TScores& scores,
                    Score<TScoreValue, TSpec2> const& sc,
                    TScoreValue minScore,
                    TDiagonal diag1,
                    TDiagonal diag2,
                    TTag) {
	// Make multiple local alignment and save them in alignments container
	_localAlignment(str, alignments, scores, sc, minScore, diag1, diag2, TTag());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...