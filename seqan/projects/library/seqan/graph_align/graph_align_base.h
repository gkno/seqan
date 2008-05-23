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
//Alignment Tags: see basic/basic_tag.h

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




//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TComponentMap, typename TOrderMap, typename TComponentLength> 
inline bool
convertAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				 TComponentMap& component,
				 TOrderMap& order,
				 TComponentLength& compLength)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TComponentMap>::Type TComponent;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Check for empty graph
	if (empty(g)) return false;

	// Connected Components
	TSize numComponents = connected_components(g, component);

	// Make a directed graph to represent the ordering of the components
	// Note: Multiple vertices might have the same component
	Graph<Directed<void, WithoutEdgeId> > componentGraph;
	reserve(_getVertexString(componentGraph), numComponents);
	for(TSize i = 0; i<numComponents;++i) addVertex(componentGraph);
	
	TSize nseq = length(value(g.data_sequence));
	String<std::set<TComponent> > componentsPerSeq;
	typedef String<String<TComponent> > TOrderedComponents;
	TOrderedComponents orderedComponentsPerSeq;
	resize(componentsPerSeq, nseq);
	resize(orderedComponentsPerSeq, nseq);
	typename TPosToVertexMap::const_iterator it1 = g.data_pvMap.begin();
	typename TPosToVertexMap::const_iterator it1End = g.data_pvMap.end();
	for(;it1!=it1End;++it1) {
		// If sections are not assigned to a vertex -> no alignment
		if (it1->second == nilVertex) return false;
		
		// Remember the sequence that component belongs to 
		TSize currentSeq = idToPosition(value(g.data_sequence), it1->first.first);

		// Append component
		TComponent c = getProperty(component, it1->second);
		if ((value(componentsPerSeq,currentSeq)).empty()) {
			String<TComponent> tmp;
			value(orderedComponentsPerSeq, currentSeq) = tmp;
		}
		appendValue(value(orderedComponentsPerSeq, currentSeq), c);
		// If two components appear twice in the same sequence -> no alignment
		if (!((value(componentsPerSeq,currentSeq)).insert(c)).second) return false;	
	}
	clear(componentsPerSeq);

	// Draw edges for the components within a sequence
	typedef typename Iterator<TOrderedComponents>::Type TIterTOrderedComponents;
	TIterTOrderedComponents itBegin = begin(orderedComponentsPerSeq);
	TIterTOrderedComponents itEnd = end(orderedComponentsPerSeq);
	for(;itBegin != itEnd; ++itBegin) {
		TSize n = length(*itBegin);
		for(TSize i = 0; i<n-1; ++i) {
			addEdge(componentGraph, value((*itBegin), i), value((*itBegin), i+1));
		}
	}
	
	// Make a topological sort of the component graph
	topological_sort(componentGraph, order);

	//// Debug code
	//std::cout << "Topological sort: " << std::endl;
	//for(TSize i = 0; i<length(order);++i) {
	//	std::cout << order[i] << ',';
	//}
	//std::cout << std::endl;

	// Walk through all sequences and check the component order
	unsigned int compIndex = 0;
	unsigned int compIndexLen = length(order);
	typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();
	TIdType currentSeq = it->first.first;
	for(; it != g.data_pvMap.end(); ++it) {
		if (it->first.first != currentSeq) {
			compIndex = 0;
			currentSeq = it->first.first;
		}
		unsigned int c = getProperty(component, it->second);
		compLength.insert(std::make_pair(c, fragmentLength(g, it->second)));
		while ((compIndex < compIndexLen) && (order[compIndex] != c)) ++compIndex;
		// Crossing components -> no alignment
		if (compIndex >= compIndexLen) return false;
		// Next component
		++compIndex;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.convertAlignment:
..cat:Graph.Alignment Graph
..summary:Converts an alignment graph into an alignment matrix.
..signature:convertAlignment(g, matrix)
..param.g:In-parameter: An alignment graph.
...type:Spec.Alignment Graph
..param.matrix:Out-parameter: A string that represents an alignment matrix.
..returns: A bool that is true iff the alignment graph is a valid alignment
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix> 
inline bool
convertAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				 TMatrix& mat)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TIdType;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;
	typedef std::map<unsigned int, unsigned int> TComponentLength;
	
	// Strongly Connected Components, topological sort, and length of each component
	String<unsigned int> component;
	String<unsigned int> order;
	TComponentLength compLength;

	if (!convertAlignment(g, component, order, compLength)) return false;

	// Create the matrix
	TSize len = 0;
	TSize nseq = length(stringSet(g));
	for(TComponentLength::iterator cIt=compLength.begin(); cIt != compLength.end(); ++cIt) len+=cIt->second;
	char gapChar = gapValue<char>();
	fill(mat, len * nseq, gapChar);

	// Fill the matrix
	TSize row = 0;
	TSize col = 0;
	typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();
	unsigned int compIndex = 0;
	unsigned int compIndexLen = length(order);
	TIdType currentSeq = it->first.first;
	for(; it != g.data_pvMap.end(); ++it) {
		if (it->first.first != currentSeq) {
			SEQAN_TASSERT(col <= len);
			//std::cout << std::endl;
			++row;col=0;
			compIndex = 0;
			currentSeq = it->first.first;
		}
		unsigned int c = getProperty(component, it->second);
		while ((compIndex < compIndexLen) && (order[compIndex] != c)) {
			for(TSize i=0;i<compLength[order[compIndex]];++i) {
				//std::cout << gapValue<char>();
				assignValue(mat, row*len + col, gapValue<char>() );
				++col;
			}
			++compIndex;
		}
		String<TValue> str = label(g,it->second);
		//std::cout << str;
		for(TSize i=0;i<length(str);++i) {
			assignValue(mat, row*len + col, (TValue) getValue(str, i));
			++col;
		}
		++compIndex;
	}
	SEQAN_TASSERT(row + 1 == nseq);
	//std::cout << std::endl;

	return true;
}


template<typename TStringSet, typename TCargo, typename TSpec>
inline void
rebuildGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Size<TGraph>::Type TSize;

	// Initialization
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	TSize nseq = length(stringSet(g));

	// Collect all character pairs
	typedef std::pair<TSize, TSize> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);	
	TEdgeIterator itE(g);
	for(;!atEnd(itE);++itE) {
		TVertexDescriptor sV = sourceVertex(itE);
		TVertexDescriptor tV = targetVertex(itE);
		TSize seq1 = idToPosition(stringSet(g), sequenceId(g, sV));
		TSize seq2 = idToPosition(stringSet(g), sequenceId(g, tV));
		TSize index = 0;
		TSize pos1 = 0;
		TSize pos2 = 0;
		if (seq1 < seq2) {
			index = seq1 * nseq + seq2;
			pos1 = fragmentBegin(g, sV);
			pos2 = fragmentBegin(g, tV);
		} else {
			index = seq2 * nseq + seq1;
			pos1 = fragmentBegin(g, tV);
			pos2 = fragmentBegin(g, sV);
		}
		for(TSize i = 0; i<fragmentLength(g, sV); ++i) {
			resPair[index].insert(std::make_pair(pos1 + i, pos2 + i));
		}
	}

	// Rebuild the graph with maximal segments
	for(TSize i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		TSize startMatch1 = pos->first;
		TSize startMatch2 = pos->second;
		TSize len = 1;
		++pos;
		while(pos != posEnd) {
			if ((startMatch1 + len == pos->first) && (startMatch2 + len == pos->second)) ++len;
			else {
				appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
				startMatch1 = pos->first;
				startMatch2 = pos->second;
				len = 1;
			}
			++pos;
		}
		appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
	}
	clearVertices(g);
	matchRefinement(matches,stringSet(g),g);
}





}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
