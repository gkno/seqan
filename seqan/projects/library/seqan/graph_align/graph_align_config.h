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

 ============================================================================
  $Id: graph_align_config.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_CONFIG_H
#define SEQAN_HEADER_GRAPH_ALIGN_CONFIG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//	Graph - AlignConfig
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Class.AlignConfig:
..cat:Alignments
..summary:The AlignConfig class encapsulates how DP is carried out. 
It indicates at what ends gaps are free, the so-called free ends-space alignments.
..signature:AlignConfig<bool TTop, bool TLeft, bool TRight, bool TBottom, TSpec>
..param.TTop:If true then 0's in top row.
...default:$false$
..param.TLeft:If true then 0's in the left row.
...default:$false$
..param.TRight:If true then maximum is also searched in the last column.
...default:$false$
..param.TBottom:If true then maximum is also searched in the last row.
...default:$false$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
*/
template<bool TTop = false, bool TLeft = false, bool TRight = false, bool TBottom = false, typename TSpec = Default>
class AlignConfig;

// 1 config
template<typename TSpec>
class AlignConfig<false, false, false, false, TSpec> 
{
};

// 2 config
template<typename TSpec>
class AlignConfig<false, false, false, true, TSpec> 
{
};

// 3 config
template<typename TSpec>
class AlignConfig<false, false, true, false, TSpec> 
{
};

// 4 config
template<typename TSpec>
class AlignConfig<false, false, true, true, TSpec> 
{
};

// 5 config
template<typename TSpec>
class AlignConfig<false, true, false, false, TSpec> 
{
};

// 6 config
template<typename TSpec>
class AlignConfig<false, true, false, true, TSpec> 
{
};

// 7 config
template<typename TSpec>
class AlignConfig<false, true, true, false, TSpec> 
{
};

// 8 config
template<typename TSpec>
class AlignConfig<false, true, true, true, TSpec> 
{
};

// 9 config
template<typename TSpec>
class AlignConfig<true, false, false, false, TSpec> 
{
};

// 10 config
template<typename TSpec>
class AlignConfig<true, false, false, true, TSpec> 
{
};

// 11 config
template<typename TSpec>
class AlignConfig<true, false, true, false, TSpec> 
{
};

// 12 config
template<typename TSpec>
class AlignConfig<true, false, true, true, TSpec> 
{
};

// 13 config
template<typename TSpec>
class AlignConfig<true, true, false, false, TSpec> 
{
};

// 14 config
template<typename TSpec>
class AlignConfig<true, true, false, true, TSpec> 
{
};

// 15 config
template<typename TSpec>
class AlignConfig<true, true, true, false, TSpec> 
{
};

// 16 config
template<typename TSpec>
class AlignConfig<true, true, true, true, TSpec> 
{
};


//////////////////////////////////////////////////////////////////////////////
//	FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TRight, bool TBottom, typename TSpec, typename TElement, typename TCost>
inline void
_initFirstColumn(AlignConfig<TTop, false, TRight, TBottom, TSpec> const,
				 TElement& element,
				 TCost const cost)
{
	SEQAN_CHECKPOINT
	element = cost;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TRight, bool TBottom, typename TSpec, typename TElement, typename TCost>
inline void
_initFirstColumn(AlignConfig<TTop, true, TRight, TBottom, TSpec> const,
				 TElement& element,
				 TCost const)
{
	SEQAN_CHECKPOINT
	element = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TLeft, bool TRight, bool TBottom, typename TSpec, typename TElement, typename TCost>
inline void
_initFirstRow(AlignConfig<false, TLeft, TRight, TBottom, TSpec> const,
			  TElement& element,
			  TCost const cost)
{
	SEQAN_CHECKPOINT
	element = cost;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TLeft, bool TRight, bool TBottom, typename TSpec, typename TElement, typename TCost>
inline void
_initFirstRow(AlignConfig<true, TLeft, TRight, TBottom, TSpec> const,
			  TElement& element,
			  TCost const)
{
	SEQAN_CHECKPOINT
	element = 0;
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastRow(AlignConfig<TTop, TLeft, TRight, false, TSpec> const,
		 TValue1&,
		 TIndex1&,
		 TValue2 const,
		 TIndex2 const)
{
	SEQAN_CHECKPOINT
	// Nop
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastRow(AlignConfig<TTop, TLeft, TRight, true, TSpec> const,
		 TValue1& maxValue,
		 TIndex1& maxIndex,
		 TValue2 const val,
		 TIndex2 const index)
{
	SEQAN_CHECKPOINT
	if (val > maxValue[0]) {
		maxValue[0] = val;
		maxIndex[0] = index;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TColumn>
inline void
_lastColumn(AlignConfig<TTop, TLeft, false, TBottom, TSpec> const,
			TValue1& maxValue,
			TIndex1&,
			TColumn const& column)
{
	SEQAN_CHECKPOINT
	maxValue[1] = column[length(column) - 1];
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TColumn>
inline void
_lastColumn(AlignConfig<TTop, TLeft, true, TBottom, TSpec> const,
			TValue1& maxValue,
			TIndex1& maxIndex,
			TColumn const& column)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TColumn>::Type TSize;
	typedef typename Iterator<TColumn, Standard>::Type TColIter;
	TSize limit = length(column) - 1;
	maxValue[1] = column[limit];
	TColIter itCol = begin(column, Standard());
	TColIter itColEnd = end(column, Standard());
	for(TSize i = 0;itCol != itColEnd; ++i, ++itCol) {
		if (*itCol > maxValue[1]) {
			maxValue[1] = *itCol;
			maxIndex[1] = i;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TScoreValue, bool TTop, bool TLeft, typename TSpec, typename TValue, typename TIndex, typename TSize>
inline TScoreValue
_maxOfAlignment(AlignConfig<TTop, TLeft, false, false, TSpec> const,
				TValue& maxValue,
				TIndex&,
				TSize const,
				TSize const)
{
	SEQAN_CHECKPOINT
	return maxValue[1];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TScoreValue, bool TTop, bool TLeft, typename TSpec, typename TValue, typename TIndex, typename TSize>
inline TScoreValue
_maxOfAlignment(AlignConfig<TTop, TLeft, true, false, TSpec> const,
				TValue& maxValue,
				TIndex& maxIndex,
				TSize const len1,
				TSize const)
{
	SEQAN_CHECKPOINT
	maxIndex[0] = len1;
	return maxValue[1];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TScoreValue, bool TTop, bool TLeft, typename TSpec, typename TValue, typename TIndex, typename TSize>
inline TScoreValue
_maxOfAlignment(AlignConfig<TTop, TLeft, false, true, TSpec> const,
				TValue& maxValue,
				TIndex& maxIndex,
				TSize const,
				TSize const len2)
{
	SEQAN_CHECKPOINT
	maxIndex[1] = len2;
	return maxValue[0];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TScoreValue, bool TTop, bool TLeft, typename TSpec, typename TValue, typename TIndex, typename TSize>
inline TScoreValue
_maxOfAlignment(AlignConfig<TTop, TLeft, true, true, TSpec> const,
				TValue& maxValue,
				TIndex& maxIndex,
				TSize const len1,
				TSize const len2)
{
	SEQAN_CHECKPOINT
	// Find the maximum
	if (maxValue[1] > maxValue[0]) maxIndex[0] = len1;
	else maxIndex[1] = len2;
	return (maxValue[0] > maxValue[1]) ? maxValue[0] : maxValue[1];
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastColumn(AlignConfig<TTop, TLeft, false, TBottom, TSpec> const,
			TValue1& maxValue,
			TIndex1& maxIndex,
			TValue2 const val,	
			TIndex2 const row,
			TIndex2 const col)
{
	SEQAN_CHECKPOINT
	maxValue[1] = val; maxIndex[2] = row; maxIndex[3] = col;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastColumn(AlignConfig<TTop, TLeft, true, TBottom, TSpec> const,
			TValue1& maxValue,
			TIndex1& maxIndex,
			TValue2 const val,
			TIndex2 const row,
			TIndex2 const col)
{
	SEQAN_CHECKPOINT
	if (val > maxValue[1]) {maxValue[1] = val; maxIndex[2] = row; maxIndex[3] = col; }
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastRow(AlignConfig<TTop, TLeft, TRight, false, TSpec> const,		
		 TValue1& maxValue,
		 TIndex1& maxIndex,
		 TValue2 const val,
		 TIndex2 const row,
		 TIndex2 const col)
{
	SEQAN_CHECKPOINT
	maxValue[0] = val; maxIndex[0] = row; maxIndex[1] = col;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_lastRow(AlignConfig<TTop, TLeft, TRight, true, TSpec> const,
		 TValue1& maxValue,
		 TIndex1& maxIndex,
		 TValue2 const val,
		 TIndex2 const row,
		 TIndex2 const col)
{
	SEQAN_CHECKPOINT
	if (val > maxValue[0]) {maxValue[0] = val; maxIndex[0] = row; maxIndex[1] = col; }
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
