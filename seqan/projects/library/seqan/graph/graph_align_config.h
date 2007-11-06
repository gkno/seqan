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
  $Id$
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
..include:graph.h
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

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue1, typename TIndex1, typename TValue2, typename TIndex2>
inline void
_processLastRow(AlignConfig<TTop, TLeft, TRight, false, TSpec> const,
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
_processLastRow(AlignConfig<TTop, TLeft, TRight, true, TSpec> const,
				TValue1& maxValue,
				TIndex1& maxIndex,
				TValue2 const val,
				TIndex2 const index)
{
	SEQAN_CHECKPOINT
	if (val > maxValue.first) {
		maxValue.first = val;
		maxIndex.first = index;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TColumn>
inline void
_processLastColumn(AlignConfig<TTop, TLeft, false, TBottom, TSpec> const,
				   TValue1& maxValue,
				   TIndex1& maxIndex,
				   TColumn const& column)
{
	SEQAN_CHECKPOINT
	maxIndex.second = length(column) - 1;
	maxValue.second = getValue(column, maxIndex.second);
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TBottom, typename TSpec, typename TValue1, typename TIndex1, typename TColumn>
inline void
_processLastColumn(AlignConfig<TTop, TLeft, true, TBottom, TSpec> const,
				   TValue1& maxValue,
				   TIndex1& maxIndex,
				   TColumn const& column)
{
	SEQAN_CHECKPOINT
	maxIndex.second = length(column) - 1;
	maxValue.second = getValue(column, maxIndex.second);
	unsigned int limit = maxIndex.second;
	for(unsigned int i = 1; i<limit; ++ i) {
		if (getValue(column, i) > maxValue.second) {
			maxValue.second = getValue(column, i);
			maxIndex.second = i;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, typename TSpec, typename TValue>
inline typename TValue::first_type
_retrieveMaxOfAlignment(AlignConfig<TTop, TLeft, TRight, false, TSpec> const,
						TValue& maxValue)
{
	SEQAN_CHECKPOINT
	return maxValue.second;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, typename TSpec, typename TValue>
inline typename TValue::first_type
_retrieveMaxOfAlignment(AlignConfig<TTop, TLeft, false, true, TSpec> const,
						TValue& maxValue)
{
	SEQAN_CHECKPOINT
	return maxValue.first;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, typename TSpec, typename TValue>
inline typename TValue::first_type
_retrieveMaxOfAlignment(AlignConfig<TTop, TLeft, true, true, TSpec> const,
						TValue& maxValue)
{
	SEQAN_CHECKPOINT
	if (maxValue.first > maxValue.second) return maxValue.first;
	else return maxValue.second;
}

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
