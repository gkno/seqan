#ifndef SEQAN_HEADER_GRAPH_ALIGN_HIRSCHBERG_H
#define SEQAN_HEADER_GRAPH_ALIGN_HIRSCHBERG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Hirschberg Alignment - TODO!!!
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace>
void
_align_hirschberg_trace(TAlign& align,		 
						TStringSet const& str,
						TTrace& trace)
{
	SEQAN_CHECKPOINT
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Iterator<TTrace>::Type TTraceIter;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Id<TStringSet>::Type TId;

	TTraceIter traceIter = begin(trace, Standard() );
	TTraceIter traceIterEnd = end(trace, Standard() );

	/*	
	// Debug
	unsigned int count = 0;
	for(;traceIter != traceIterEnd;goNext(traceIter)) {
	  std::cout << count << ',' << getValue(traceIter) << std::endl;
	  ++count;
	}
	*/

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};
	
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);	 
	TSize len1 = length(str[0]);
	TSize currentPointer = getValue(trace, len1);
	TSize movePointer = getValue(trace, len1 - 1);
	do {
		if ((currentPointer - movePointer) == 1) {
			TSize segLen = 0;
			while ((currentPointer - movePointer) == 1) {
				++segLen;
				--len1;
				--currentPointer;
				if (len1>0) movePointer = getValue(trace, len1 - 1);
				else break;
			}
			_align_trace_print(align, str, id1, len1, id2, currentPointer, segLen, (Byte) Diagonal);
			//std::cout << "Diagonal " << segLen << std::endl;
		} else if ((currentPointer - movePointer) == 0) {
			TSize segLen = 0;
			while ((currentPointer - movePointer) == 0) {
				++segLen;
				--len1;
				if (len1>0) movePointer = getValue(trace, len1 - 1);
				else break;
			}
			_align_trace_print(align, str, id1, len1, id2, (TSize) 0, segLen, (Byte) Horizontal);
			//std::cout << "Horizontal " << segLen << std::endl;
		} else {
			_align_trace_print(align, str, id1, (TSize) 0, id2, (TSize) movePointer + 1, (TSize) currentPointer - (movePointer + 1), (Byte) Vertical);
			//std::cout << "Vertical " << currentPointer - (movePointer + 1)  << std::endl;
			currentPointer = movePointer + 1;
		}
	} while (len1 != 0);
	if (getValue(trace, 0) != 0) {
		_align_trace_print(align, str, id1, (TSize) 0, id2, (TSize) 0, (TSize) getValue(trace, 0), (Byte) Vertical);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScoreValue>
TScoreValue
_align_hirschberg(TTrace& trace,
				  TStringSet const& str,
				  Score<TScoreValue, Simple> const& sc,
				  Hirschberg) 
{
	SEQAN_CHECKPOINT
	typedef String<TScoreValue> TColumn;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	typedef typename Value<TStringSet const>::Type TString;
	typedef typename Infix<TString const>::Type TInfix;
	typedef typename Size<TString>::Type TSize;
	
	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue gapOpen = scoreGapOpen(sc);
	TScoreValue upperLeft = 0;

	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
	resize(trace, len1+1);
	typedef std::pair<TSize, TSize> TPoint;
	typedef std::list<TPoint> TMidPointQueue;
	TMidPointQueue midpoints;
	TSize x1 = 0;
	TSize x2 = 0;
	TSize y1 = 0;
	TSize y2 = 0;
	bool firstRun = true;
	TTraceValue tvMat, tvHorizontal, tvVertical;
	do {
		// Initialization
		upperLeft = 0;
		tvMat = (Byte) Diagonal;
	  
		// Iterator until the number of midpoints equals the sequence length
		firstRun = true;
		typename TMidPointQueue::iterator it = midpoints.begin();
		do {
			// Step1: Get the alignment window
			if (firstRun) {
				x1 = 0;
				y1 = 0;
				firstRun = false;
			} else {
				x1 = x2;
				y1 = y2;
				++it;
			}

			if (it == midpoints.end()) {
				x2 = len1;
				y2 = len2;
			} else {
				TPoint p = *it;
				x2 = p.first;
				y2 = p.second;
			}

			// Step2: Do the alignment
			TSize middle = x1 + ((x2-x1)/2);
			TInfix infix1 = infix(str[0], x1, x2);
			TInfix infix2 = infix(str[1], y1, y2);
			TSize inf_middle = (TSize) ((x2-x1)/2);
		

			// The DP Matrix for diagonal walks
			TColumn mat;
			// The DP Matrix for gaps from the left
			TColumn horizontal;
			// The DP Matrix for gaps from the top
			TScoreValue vert;

			// The pointers
			unsigned int verticalPointer;
			String<unsigned int> horizontalPointer;
			String<unsigned int> diagonalPointer;
			unsigned int diagonalPointerOld;
			unsigned int tmpPointer;

			// Initialization
			TSize inf_len1 = x2-x1;
			TSize inf_len2 = y2-y1;
			TScoreValue tmp = 0;
			TScoreValue inf = infimumValue<TScoreValue>() / 2;
			resize(mat, (inf_len2+1));   // One column for the diagonal matrix
			vert = inf;
			fill(horizontal, (inf_len2+1), inf);  
			resize(horizontalPointer, (inf_len2+1));  
			resize(diagonalPointer, (inf_len2+1));  
				
			// Classical DP
			assignValue(mat, 0, upperLeft);
			verticalPointer = 0;
			assignValue(horizontalPointer, 0, 0);
			assignValue(diagonalPointer, 0, 0);
			for(TSize row = 1; row <= inf_len2; ++row) {
				if (tvMat == (Byte) Vertical) {
					assignValue(mat, row, upperLeft + gap + (row - 1) * gap);
					assignValue(horizontal, row, getValue(mat, row));
				} else {
					assignValue(mat, row, upperLeft + gapOpen + (row - 1) * gap);
					assignValue(horizontal, row, getValue(mat, row) + gapOpen - gap);
				}
				assignValue(horizontalPointer, row, row);
				assignValue(diagonalPointer, row, row);
			}
			for(TSize col = 1; col <= inf_len1; ++col) {
				diagonalPointerOld = getValue(diagonalPointer,0);
				assignValue(diagonalPointer, 0, getValue(horizontalPointer,0));
				verticalPointer = getValue(horizontalPointer,0);
				TScoreValue diagValMat = getValue(mat, 0);
				if (tvMat == (Byte) Horizontal) {
					assignValue(mat, 0, upperLeft + gap + (col - 1) * gap);
					vert = getValue(mat, 0);
				} else {			    
					assignValue(mat, 0, upperLeft + gapOpen + (col - 1) * gap);
					vert = getValue(mat, 0) + gapOpen - gap;
				}
				// If there are no rows we step horizontally !!!
				tvMat = (Byte) Horizontal;
				for(TSize row = 1; row <= inf_len2; ++row) {
					// Get the new maximum for vertical
					if ((tmp = getValue(mat, row - 1) + gapOpen) > vert + gap) {
						vert = tmp;
						tvVertical = (Byte) Diagonal;
					} else {
						vert = vert + gap;
						tvVertical = (Byte) Vertical;
					}
					
					// Get the new maximum for left
					if ((tmp = getValue(mat, row) + gapOpen) > getValue(horizontal, row) + gap) {
						assignValue(horizontal, row, tmp);
						tvHorizontal = (Byte) Diagonal;
					} else {
						assignValue(horizontal, row, getValue(horizontal, row) + gap);
						tvHorizontal = (Byte) Horizontal;
					}
					
					// Get the new maximum for mat
					TScoreValue sc_ = score(const_cast<Score<TScoreValue, Simple>&>(sc), infix1[col-1], infix2[row-1]);
					tmp = diagValMat + sc_;
					tvMat = (Byte) Diagonal;
					if (vert > tmp) {
						tmp = vert;
						tvMat = (Byte) Vertical;
					}
					if (getValue(horizontal, row) > tmp) {
						tmp = getValue(horizontal,row);
						tvMat = (Byte) Horizontal;
					}
					
					// Assign the new diagonal values
					diagValMat = getValue(mat, row);
					assignValue(mat, row, tmp);
					
					// Take care of the pointers
					if (col > inf_middle) {
						if (tvVertical == (Byte) Diagonal) assignValue(verticalPointer, row, getValue(diagonalPointer, row - 1));
						if (tvHorizontal == (Byte) Diagonal) assignValue(horizontalPointer, row, getValue(diagonalPointer, row));
						
						tmpPointer = getValue(diagonalPointer, row);
						if (tvMat == (Byte) Diagonal) assignValue(diagonalPointer, row, diagonalPointerOld);
						else if (tvMat == (Byte) Horizontal) assignValue(diagonalPointer, row, getValue(horizontalPointer, row));
						else if (tvMat == (Byte) Vertical) assignValue(diagonalPointer, row, verticalPointer);
						diagonalPointerOld = tmpPointer;
					}
					//std::cout << getValue(diagonalPointer,row) << ',' << getValue(horizontalPointer, row) << ',' << verticalPointer << ';';
				}
				//std::cout << std::endl;
			}
			tmp = getValue(mat, inf_len2);
			tmpPointer = getValue(diagonalPointer, inf_len2);
			if (getValue(horizontal, inf_len2) ==  tmp) {
				tmpPointer = getValue(horizontalPointer, inf_len2);
			}
			else if (vert == tmp) {
				tmpPointer = verticalPointer;
			}
			TSize cut = y1 + tmpPointer;
			upperLeft = tmp;

			/*
			// Debug code
			std::cout << infix1 << std::endl;
			std::cout << infix2 << std::endl;
			std::cout << "Trace-Point" << middle << ',' << cut << std::endl;
			std::cout << "Max-Score" << tmp << std::endl;
			std::cout << "--" << std::endl;
			*/

			// Step3: Insert the new midpoint
			if ((x2 - x1) > 1) {
				it = midpoints.insert(it, std::make_pair(middle,cut));
				++it;
			}
		} while (it != midpoints.end());
	} while ((midpoints.size() < len1 - 1));
	TTraceIter traceIter = begin(trace, Standard() );
	typename TMidPointQueue::iterator it = midpoints.begin();
	if (it->second > 0) assignValue(traceIter, it->second - 1); // Start point of the trace
	else assignValue(traceIter, 0); 
	goNext(traceIter);
	for(;it!=midpoints.end();++it) {
		assignValue(traceIter, (it->second));
		goNext(traceIter);
	}
	assignValue(traceIter, len2); // End point of the trace
	
	return upperLeft;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 Score<TScoreValue, Simple> const& sc,
				 Hirschberg)
{
	SEQAN_CHECKPOINT
	String<unsigned int, External<> > trace;
	TScoreValue maxScore = _align_hirschberg(trace, str, sc, Hirschberg());

	// Follow the trace and create the graph
	_align_hirschberg_trace(align, str, trace);

	return maxScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TStringSet const& str,
				 Score<TScoreValue, Simple> const& sc,
				 Hirschberg)
{
	SEQAN_CHECKPOINT
	String<unsigned int> trace;
	return _align_hirschberg(trace, str, sc, Hirschberg());	
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
