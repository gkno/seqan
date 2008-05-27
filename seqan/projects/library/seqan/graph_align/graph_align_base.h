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
  $Id: graph_align_base.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_BASE_H
#define SEQAN_HEADER_GRAPH_ALIGN_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{




//////////////////////////////////////////////////////////////////////////////
// Alignment: Simple Traceback Alphabet
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_Byte_2_TraceBack
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_TraceBack<T>::VALUE[256] = 
{
	0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.TraceBack:
..cat:Alphabets
..summary: Trace back values.
..general:Class.SimpleType
..signature:TraceBack
..remarks:
...text:The @Metafunction.ValueSize@ of $TraceBack$ is 3. 
The values are defined in the following way: 0=Diagonal Move, 1=Horizontal Move, 2=Vertical Move
..see:Metafunction.ValueSize
*/
struct _TraceBack {};
typedef SimpleType<unsigned char, _TraceBack> TraceBack;

template <> struct ValueSize< TraceBack > { enum { VALUE = 3 }; };
template <> struct BitsPerValue< TraceBack > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<TraceBack, Byte> { typedef TraceBack Type; };
inline void assign(TraceBack & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_TraceBack<>::VALUE[c_source];
}



//////////////////////////////////////////////////////////////////////////////
// Alignment: Extended Traceback Alphabet (Gotoh)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template <typename T = void>
struct _Translate_Table_Byte_2_TraceBackGotoh
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_TraceBackGotoh<T>::VALUE[256] = 
{
	0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,   11,   12,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.TraceBackGotoh:
..cat:Alphabets
..summary: Trace back values for gotoh.
..general:Class.SimpleType
..signature:TraceBackGotoh
..remarks:
...text:The @Metafunction.ValueSize@ of $TraceBackGotoh$ is 12. 
The values are defined in the following way: Move in Diagonal Matrix, Move in Horizontal Matrix, Move in Vertical Matrix
The values are: 
0=Diag, Diag, Diag; 1=Diag, Diag, Vert; 2=Diag, Hori, Diag; 3=Diag, Hori, Vert; 
4=Hori, Diag, Diag; 5=Hori, Diag, Vert; 6=Hori, Hori, Diag; 7=Hori, Hori, Vert; 
8=Vert, Diag, Diag; 9=Vert, Diag, Vert; 10=Vert, Hori, Diag; 11=Vert, Hori, Vert; 
12 = Stop (For SmithWaterman Traceback)
..see:Metafunction.ValueSize
*/
struct _TraceBackGotoh {};
typedef SimpleType<unsigned char, _TraceBackGotoh> TraceBackGotoh;

template <> struct ValueSize< TraceBackGotoh > { enum { VALUE = 13 }; };
template <> struct BitsPerValue< TraceBackGotoh > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<TraceBackGotoh, Byte> { typedef TraceBackGotoh Type; };
inline void assign(TraceBackGotoh & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_TraceBackGotoh<>::VALUE[c_source];
}





//////////////////////////////////////////////////////////////////////////////
// Alignment: Trace-back, internal functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringSet, typename TId, typename TPos, typename TTraceValue>
inline void
_align_trace_print(TFile& file,
				   TStringSet const& str,
				   TId const id1,
				   TPos const pos1,
				   TId const id2,
				   TPos const pos2,
				   TPos const segLen,
				   TTraceValue const tv)
{
	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	if (segLen == 0) return;

	if (tv == (Byte) Horizontal) {
		for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
			_streamPut(file, '(');
			_streamPut(file, (str[0])[i]);
			_streamPut(file, ',');
			_streamPut(file, gapValue<char>());
			_streamPut(file, ')');
			_streamPut(file, '\n');
		}
	}
	else if (tv == (Byte) Vertical) {
		for (int i = pos2 + segLen - 1; i>= (int) pos2;--i) {
			_streamPut(file, '(');
			_streamPut(file, gapValue<char>());
			_streamPut(file, ',');
			_streamPut(file, (str[1])[i]);
			_streamPut(file, ')');
			_streamPut(file, '\n');
		}
	}
	else if (tv == (Byte) Diagonal) {
		int j = pos2 + segLen - 1;
		for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
			_streamPut(file, '(');
			_streamPut(file, (str[0])[i]);
			_streamPut(file, ',');
			_streamPut(file, (str[1])[j]);
			_streamPut(file, ')');
			_streamPut(file, '\n');
			--j;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TStringSet2, typename TId, typename TPos, typename TTraceValue>
inline void
_align_trace_print(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				   TStringSet2 const&,
				   TId const id1,
				   TPos const pos1,
				   TId const id2,
				   TPos const pos2,
				   TPos const segLen,
				   TTraceValue const tv)
{
	SEQAN_CHECKPOINT

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	if (segLen == 0) return;

	if (tv == (Byte) Horizontal) addVertex(g, id1, pos1, segLen);
	else if (tv == (Byte) Vertical) addVertex(g, id2, pos2, segLen);
	else if (tv == (Byte) Diagonal) addEdge(g, addVertex(g, id1, pos1, segLen), addVertex(g, id2, pos2, segLen));
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFragment, typename TStringSet, typename TId, typename TPos, typename TTraceValue>
inline void
_align_trace_print(String<TFragment>& matches,
				   TStringSet const&,
				   TId const id1,
				   TPos const pos1,
				   TId const id2,
				   TPos const pos2,
				   TPos const seqLen,
				   TTraceValue const tv)
{
	SEQAN_CHECKPOINT

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	if (seqLen == 0) return;

	if (tv == (Byte) Horizontal) {
		// Nop, no match
	} else if (tv == (Byte) Vertical) {
		// Nop, no match
	} else if (tv == (Byte) Diagonal) {
		appendValue(matches, TFragment(id1, pos1, id2, pos2, seqLen));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TVertexDescriptor, typename TSpec, typename TStringSet, typename TId, typename TPos, typename TTraceValue>
inline void
_align_trace_print(String<String<TVertexDescriptor, TSpec> >& nodeString,
				   TStringSet const& str,
				   TId const,
				   TPos const pos1,
				   TId const,
				   TPos const pos2,
				   TPos const segLen,
				   TTraceValue const tv)
{
	typedef String<TVertexDescriptor, TSpec> TVertexDescriptorString;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Iterator<TVertexDescriptorString>::Type TStringIter;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	if (segLen == 0) return;
	// Number of vertex descriptors in the first string at any position (e.g., group of 5 sequences = group of 5 vertex descriptors)
	TSize len1 = length(getValue(getValue(str,0), 0));
	// Number of vertex descriptors in the second string at any position (e.g., group of 5 sequences = group of 5 vertex descriptors)
	TSize len2 = length(getValue(getValue(str,1), 0));

	// Resize the node string
	TSize index = length(nodeString);
	resize(nodeString, index + segLen);

	if (tv == (Byte) Horizontal) {
		for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
			fill(value(nodeString, index), len1 + len2, nilVertex);
			TStringIter it = begin(value(nodeString, index));
			for(TPos all = 0;all<len1;++all) {
				*it = getValue(getValue(getValue(str,0),i), all);
				goNext(it);
			}
			++index;
		}
	}
	else if (tv == (Byte) Vertical) {
		for (int i = pos2 + segLen - 1; i>= (int) pos2;--i) {
			fill(value(nodeString, index), len1 + len2, nilVertex);
			TStringIter it = begin(value(nodeString, index));
			it+=len1;
			for(TPos all = 0;all<len2;++all) {
				*it = getValue(getValue(getValue(str,1),i), all);
				goNext(it);
			}
			++index;
		}
	}
	else if (tv == (Byte) Diagonal) {
		int j = pos2 + segLen - 1;
		for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
			resize(value(nodeString, index), len1 + len2);
			TStringIter it = begin(value(nodeString, index));
			for(TPos all = 0;all<len1;++all) {
				*it = getValue(getValue(getValue(str,0),i), all);
				goNext(it);
			}
			for(TPos all = 0;all<len2;++all) {
				*it = getValue(getValue(getValue(str,1),j), all);
				goNext(it);
			}
			++index;
			--j;
		}
	}
}


template <typename TSourceValue, typename TType, typename TTargetValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2>
inline void
convertScoringMatrix(Score<TSourceValue, TType> const & in,
					 Score<TTargetValue, ScoreMatrix<TSequenceValue, TSpec> > & out,
					 TVal1 gapExtend,
					 TVal2 gapOpen)
{
	typedef typename Size<Score<TSourceValue, TType> >::Type TSize;
	out.data_gap_extend = (TTargetValue) gapExtend;
	out.data_gap_open = (TTargetValue) gapOpen;
	TSize alphSize = ValueSize<TSequenceValue>::VALUE;
	for(TSize row = 0; row<alphSize; ++row) {
		for(TSize col = 0; col<alphSize; ++col) {
			setScore(out,  TSequenceValue(row), TSequenceValue(col), score(in, TSequenceValue(row), TSequenceValue(col)));
		}
	}

}



template <typename TSourceValue, typename TType, typename TTargetValue, typename TSequenceValue, typename TSpec>
inline void
convertScoringMatrix(Score<TSourceValue, TType> const & in,
					 Score<TTargetValue, ScoreMatrix<TSequenceValue, TSpec> > & out)
{
	convertScoringMatrix(in, out, scoreGapExtend(in), scoreGapOpen(in));
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
