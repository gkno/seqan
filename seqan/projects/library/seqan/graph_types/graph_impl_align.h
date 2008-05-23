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
  $Id: graph_impl_align.h 2103 2008-05-23 07:57:13Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_ALIGN_H
#define SEQAN_HEADER_GRAPH_IMPL_ALIGN_H

#include <map>

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TId = unsigned int, typename TSize = unsigned int>
class FragmentInfo {
public:
	TId data_seq_id;
	TSize data_begin;
	TSize data_length;

	FragmentInfo() :
		data_seq_id(0),
		data_begin(0),
		data_length(0)
	{
	}

	FragmentInfo(TId id, TSize beg, TSize len) :
		data_seq_id(id),
		data_begin(beg),
		data_length(len)
	{
	}

};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Alignment Graph:
..cat:Graph
..general:Class.Graph
..summary:An alignment graph.
..description:
...image:alignmentGraph|An alignment graph with 3 sequences.
..signature:Graph<Alignment<TStringSet, TCargo, TSpec> > 
..param.TStringSet:The type of the string set containing the sequence information.
...default:$Class.StringSet$
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use @Metafunction.Cargo@ to get the cargo type of an undirected graph.
...default:$void$
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
...remarks:Use WithoutEdgeId here to omit edge ids.
Note: If edges do not store ids external property maps do not work.
...default:$Default$, see @Tag.Default@.
*/
template<typename TString, typename TSpecial, typename TCargo, typename TSpec>
class Graph<Alignment<StringSet<TString, Dependent<TSpecial> >, TCargo, TSpec> > 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef typename Size<Graph>::Type TSize;
		typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;

		// Alignment graph
		Graph<Undirected<TCargo, TSpec> > data_align;

		// Sequences
		Holder<StringSet<TString, Dependent<TSpecial> > > data_sequence;
		
		// Alignment specific members
		String<FragmentInfo<TIdType, TSize> > data_fragment;

		// STL Map to retrieve a vertex given SeqId, Position
		TPosToVertexMap data_pvMap;


		Graph() {
		}


		template <typename TDefault>
		Graph(StringSet<TString, Dependent<TDefault> > const& sSet) {
			SEQAN_CHECKPOINT
			data_sequence = sSet;

			// Cover all sequences with nil vertices
			TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
			for(TSize k=0; k<length(sSet);++k) {
				data_pvMap.insert(std::make_pair(std::make_pair(positionToId(sSet,k), length(sSet[k])), nilVertex));
			}
		}

		template <typename TDefault>
		Graph(StringSet<TString, Owner<TDefault> > const& sSet) {
			SEQAN_CHECKPOINT
			StringSet<TString, Dependent<> > depStr(sSet);
			data_sequence = depStr;

			// Cover all sequences with nil vertices
			TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
			for(TSize k=0; k<length(sSet);++k) {
				data_pvMap.insert(std::make_pair(std::make_pair(positionToId(const_cast<StringSet<TString, Owner<TDefault> >&>(sSet),k), length(sSet[k])), nilVertex));
			}
		}


		~Graph() {
			SEQAN_CHECKPOINT
			clear(*this);
		}

		Graph(Graph const & _other) 
		{
			SEQAN_CHECKPOINT
			_copyGraph(_other, *this);	
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			clear(*this);
			_copyGraph(_other, *this);
			return *this;
		}
};

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type*>&
_getVertexString(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	return const_cast<String<TEdgeStump*>&>(g.data_align.data_vertex);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename VertexIdHandler<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
_getVertexIdManager(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexIdHandler<TGraph>::Type TVertexIdManager;
	return const_cast<TVertexIdManager&>(g.data_align.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename EdgeIdHandler<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
_getEdgeIdManager(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) {
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename EdgeIdHandler<TGraph>::Type TEdgeIdManager;
	return const_cast<TEdgeIdManager&>(g.data_align.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& source,
		   Graph<Alignment<TStringSet, TCargo, TSpec> >& dest,
		   bool) 
{
	SEQAN_CHECKPOINT
	clear(dest);
	dest.data_align = source.data_align;
	dest.data_sequence = source.data_sequence;
	dest.data_fragment = source.data_fragment;
	dest.data_pvMap = source.data_pvMap;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& source,
		   Graph<Alignment<TStringSet, TCargo, TSpec> >& dest) 
{
	_copyGraph(source, dest, false); // Never transpose, underlying graph is undirected
}


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
transpose(Graph<Alignment<TStringSet, TCargo, TSpec> > const& source,
		  Graph<Alignment<TStringSet, TCargo, TSpec> >& dest)
{
	SEQAN_CHECKPOINT
	// Alignment graph, no transpose just copy
	_copyGraph(source, dest, false);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
transpose(Graph<Alignment<TStringSet, TCargo, TSpec> >&)
{
	SEQAN_CHECKPOINT
	// Nothing to do in an alignment graph
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
numEdges(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return numEdges(g.data_align);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
numVertices(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return numVertices(g.data_align);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline bool 
empty(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return empty(g.data_align);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Alignment<TStringSet, TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearEdges(g.data_align);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Alignment<TStringSet, TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	clear(g.data_fragment);
	g.data_pvMap.clear();
	clearVertices(g.data_align);


	// Don't forget to cover the sequences with nil vertices again
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	if(!empty(value(g.data_sequence))) {
		for(TSize k=0; k<length(stringSet(g));++k) {
			g.data_pvMap.insert(std::make_pair(std::make_pair(positionToId(stringSet(g),k), length(stringSet(g)[k])), nilVertex));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
clear(Graph<Alignment<TStringSet, TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	// Only clear also removes the sequences
	clear(value(g.data_sequence));
	clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
outDegree(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g, 
		  TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return outDegree(g.data_align, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
inDegree(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g, 
		 TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return inDegree(g.data_align, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
degree(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
	   TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return degree(g.data_align, vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TId, typename TPos, typename TLength> 
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
addVertex(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
		  TId id,
		  TPos begin,
		  TLength len)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef FragmentInfo<TIdType, TSize> TFragmentInfo;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;

	//for(TPosToVertexMap::const_iterator p = g.data_pvMap.begin(); p != g.data_pvMap.end(); ++p) {
	//	std::cout << p->first.first << ',' << p->first.second << ':' << p->second << std::endl;
	//}

	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Store the new fragment
	typename TPosToVertexMap::iterator interval = g.data_pvMap.lower_bound(std::make_pair(id, begin + len));
	// Segment does not belong to Sequence anymore
	SEQAN_TASSERT(interval != g.data_pvMap.end());
	// Segment end must be assigned to nil so far
	SEQAN_TASSERT(interval->second == nilVertex);
	// Segment must belong to the whole old interval
	SEQAN_TASSERT(*interval == *g.data_pvMap.upper_bound(std::make_pair(id, begin)));

	// Insert new vertex
	TVertexDescriptor vd = addVertex(g.data_align);
	if (length(g.data_fragment) <= vd) resize(g.data_fragment, vd + 1, Generous());
	assignProperty(g.data_fragment, vd, TFragmentInfo(id, begin, len));

	// Update position to vertex map
	// Does the end of the new fragment coincides with the end of the interval?
	if ( (TSize) begin + len == (TSize) interval->first.second) {
		// Does the beginning of the new fragment coincides with the beginning of the interval?
		if ((begin == 0) ||
			(g.data_pvMap.find(std::make_pair(id, begin)) != g.data_pvMap.end())) {
			// Replace interval
			interval->second = vd;
		} else {
			// Split interval once
			g.data_pvMap.insert(std::make_pair(std::make_pair(interval->first.first,begin), interval->second));
			g.data_pvMap.erase(interval);
			g.data_pvMap.insert(std::make_pair(std::make_pair(id,begin+len), vd));
		}
	} else {
		// Does the beginning of the new fragment coincides with the beginning of the interval?
		if ((begin == 0) ||
			(g.data_pvMap.find(std::make_pair(id, begin)) != g.data_pvMap.end())) {
			// Split interval once
			// Just insert here because we store interval ends
			g.data_pvMap.insert(std::make_pair(std::make_pair(id,begin+len), vd));
		} else {
			// Split interval twice
			TIdType tmp = interval->first.second;
			g.data_pvMap.insert(std::make_pair(std::make_pair(interval->first.first,begin), interval->second));
			g.data_pvMap.erase(interval);
			g.data_pvMap.insert(std::make_pair(std::make_pair(id,begin+len), vd));
			g.data_pvMap.insert(std::make_pair(std::make_pair(id,tmp), nilVertex));
		}
	}
	return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVD>
inline void 
removeVertex(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
			 TVD const v) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;

	// Clear the interval
	typename TPosToVertexMap::iterator interval = g.data_pvMap.lower_bound(std::make_pair(sequenceId(g,v), fragmentBegin(g,v) + fragmentLength(g,v)));
	SEQAN_TASSERT(interval != g.data_pvMap.end());
	interval->second = getNil<TVertexDescriptor>();

	// Remove the vertex
	removeVertex(g.data_align,v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
addEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	return addEdge(g.data_align, source, target);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TCargo2> 
inline typename EdgeDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
addEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TCargo2 const cargo) 
{
	SEQAN_CHECKPOINT
	return addEdge(g.data_align, source, target, cargo);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
		   TVertexDescriptor const source, 
		   TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	removeEdge(g.data_align, source, target);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void 
removeEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	removeEdge(g.data_align, sourceVertex(g.data_align,edge), targetVertex(g.data_align,edge));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeOutEdges(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
			   TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	removeOutEdges(g.data_align, v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline void 
removeInEdges(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
			  TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	removeInEdges(g.data_align,v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TEdgeDescriptor> 
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
targetVertex(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return targetVertex(g.data_align, edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
sourceVertex(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return sourceVertex(g.data_align, edge);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g, 
				   TMatrix& mat) 
{
	getAdjacencyMatrix(g.data_align, mat);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename EdgeDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
findEdge(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		 TVertexDescriptor const v,
		 TVertexDescriptor const w)
{
	SEQAN_CHECKPOINT
	return findEdge(g.data_align, v, w);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
	  TIDString const &,
	  Raw)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef FragmentInfo<TIdType, TSize> TSegment;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const, Rooted>::Type TIterConst;


	
	String<char> align;
	if (!convertAlignment(g, align)) {
		_streamWrite(target,"Adjacency list:\n");
		for(TIterConst it = begin(g.data_align.data_vertex);!atEnd(it);goNext(it)) {
			TVertexDescriptor sourceV = position(it);
			_streamPutInt(target, sourceV);
			TSegment seg = getProperty(g.data_fragment, sourceV);
			_streamWrite(target," (SeqId:");
			_streamPutInt(target, seg.data_seq_id);
			_streamWrite(target," ,Begin:");
			_streamPutInt(target, seg.data_begin);
			_streamWrite(target," ,Length:");
			_streamPutInt(target, seg.data_length);
			_streamWrite(target,") -> ");
			TEdgeStump* current = getValue(it);
			while(current!=0) {
				TVertexDescriptor adjV = getTarget(current);
				if (adjV != sourceV) {
					_streamPutInt(target, adjV);
					_streamPut(target, ',');
					current=getNextS(current);
				} else {
					adjV = getSource(current);
					_streamPutInt(target, adjV);
					_streamPut(target, ',');
					current=getNextT(current);
				}
			}
			_streamPut(target, '\n');
		}
		_streamWrite(target,"Edge list:\n");
		for(TIterConst it = begin(g.data_align.data_vertex);!atEnd(it);goNext(it)) {
			TVertexDescriptor sourceV = position(it);
			TEdgeStump* current = getValue(it);
			while(current!=0) {
				TVertexDescriptor targetV = getTarget(current);
				if (sourceV != targetV) {
					_streamWrite(target,"Source: ");
					_streamPutInt(target, sourceV);		
					_streamPut(target, ',');
					_streamWrite(target,"Target: ");
					_streamPutInt(target, targetV);
					_streamPut(target, ' ');
					_streamWrite(target,"(Id: ");
					_streamPutInt(target, _getId(current));
					_streamPut(target, ',');
					_streamWrite(target," Cargo-Type: ");
					_streamWrite(target, typeid(getCargo(current)).name());
					_streamPut(target, ')');
					_streamPut(target, '\n');
					current=getNextS(current);
				} else {
					current=getNextT(current);
				}
			}
		}
	} else {
		TSize nseq = length(stringSet(g));
		TSize colLen = length(align) / nseq;
		
		TSize baseCount=0;
		TSize leftSpace=6;
		TSize xPos = 0;
		_streamWrite(target,"Alignment matrix:\n");
		while (xPos < colLen) {
			TSize windowSize = 50;
			if ((xPos + windowSize)>colLen) windowSize = colLen - xPos;
			
			// Print header line
			TSize offset=0;
			// Larger numbers need to be further left
			if (baseCount != 0) offset = (unsigned int) floor(log((double)baseCount) / log((double)10));
			for(TSize j = 0;j<leftSpace-offset;++j) {
				_streamPut(target, ' ');
			}
			_streamPutInt(target, baseCount);
			baseCount+=windowSize;
			_streamPut(target, ' ');
			for(TSize col = 1;col<=windowSize;++col) {
				if ((col % 10)==0) _streamPut(target, ':');
				else if ((col % 5)==0) _streamPut(target, '.');
				else _streamPut(target, ' ');
			}
			_streamPut(target, ' ');
			_streamPut(target, '\n');

			// Print sequences
			for(TSize row=0;row<2*nseq-1;++row) {
				for(TSize col = 0;col<leftSpace+2;++col) _streamPut(target, ' ');
				if ((row % 2)==0) {
					for(TSize col = xPos;col<xPos+windowSize;++col) {
						_streamPut(target, getValue(align, (row/2)*colLen+col));
					}
				} else {
					for(TSize col = xPos;col<xPos+windowSize;++col) {
						if ((getValue(align,((row-1)/2)*colLen + col) != gapValue<char>()) &&
							(getValue(align,((row+1)/2)*colLen + col) != gapValue<char>()) &&
							(getValue(align,((row-1)/2)*colLen + col) == getValue(align,((row+1)/2)*colLen + col))) {
								_streamPut(target, '|');
						} else {
							_streamPut(target, ' ');
						}
					}
				}
				_streamPut(target, '\n');
			}
			_streamPut(target, '\n');
			xPos+=windowSize;
		}
		_streamPut(target, '\n');
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec, typename TNames>
inline void
write(TFile & file,
	  Graph<TSpec> const& g,
	  TNames const& names,
	  FastaFormat) 
{
	typedef Graph<TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;

	String<char> align;
	if (convertAlignment(g, align)) {	
		TSize nseq = length(stringSet(g));
		TSize colLen = length(align) / nseq;
		
		for(TSize i = 0; i<nseq; ++i) {
			_streamPut(file, '>');
			_streamWrite(file,names[i]);
			_streamPut(file, '\n');
			TSize col = 0;
			while(col < colLen) {
				TSize max = 0;
				if ((colLen - col) < 60) max = colLen - col;
				else max = 60;
				for(TSize finger = col; finger<col+max; ++finger) {
					_streamPut(file, getValue(align, i*colLen + finger));
				}
				col += max;
				_streamPut(file, '\n');
			}
		}
	}
}



//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec, typename TNames>
inline void
write(TFile & file,
	  Graph<TSpec> const& g,
	  TNames const& names,
	  MsfFormat) 
{
	typedef Graph<TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;

	String<char> align;
	if (convertAlignment(g, align)) {	
		TSize nseq = length(stringSet(g));
		TSize colLen = length(align) / nseq;
		
		_streamWrite(file,"PileUp\n");
		_streamPut(file, '\n');
		_streamWrite(file," MSF: ");
		_streamPutInt(file, colLen);
		_streamWrite(file," Type: P");
		_streamWrite(file," Check: 0 ..");
		_streamPut(file, '\n');
		_streamPut(file, '\n');
		TSize offset = 0;
		for(TSize i = 0; i<nseq; ++i) {
			_streamWrite(file," Name: ");
			_streamWrite(file,names[i]);
			_streamWrite(file," oo  Len:  ");
			TSize len = length(names[i]);
			if (len > offset) offset = len;
			_streamPutInt(file, colLen);
			_streamWrite(file," Check: 0");
			_streamWrite(file," Weight: 1.00");
			_streamPut(file, '\n');
		}
		offset += 5;
		_streamPut(file, '\n');
		_streamWrite(file,"//\n");
		_streamPut(file, '\n');
		_streamPut(file, '\n');
		TSize col = 0;
		while(col < colLen) {
			TSize max = 0;
			for(TSize i = 0; i<nseq; ++i) {
				if ((colLen - col) < 50) max = colLen - col;
				else max = 50;
				_streamWrite(file,names[i]);
				for(TSize j = 0; j<offset - length(names[i]); ++j) {
					_streamPut(file, ' ');
				}
				for(TSize finger = col; finger<col+max; ++finger) {
					if ((finger - col) % 10 == 0) _streamPut(file, ' ');
					if (getValue(align, i*colLen + finger) == '-') _streamPut(file, '.');
					else _streamPut(file, getValue(align, i*colLen + finger));
				}
				_streamPut(file, '\n');
			}
			col += max;
			_streamPut(file, '\n');
			_streamPut(file, '\n');
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TStringSet, typename TSpec, typename TEdge>
inline void
__writeCargo(TFile & file,
			 Graph<Alignment<TStringSet, void, TSpec> > const&,
			 TEdge const&)
{
	_streamPutInt(file, 0);
}

//////////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TEdge>
inline void
__writeCargo(TFile & file,
			 Graph<Alignment<TStringSet, TCargo, TSpec> > const&,
			 TEdge const& edge)
{
	_streamPutInt(file, getCargo(edge));
}

//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
inline void
write(TFile & file,
	  Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
	  TNames const& names,
	  CgVizFormat) 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef FragmentInfo<TId, TSize> TSegment;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const, Rooted>::Type TIterConst;
	
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);

	_streamWrite(file,"{DATA Data\n");
	_streamWrite(file,"\t[__GLOBAL__] tracks=");
	_streamPutInt(file, nseq);
	_streamPut(file, '\n');
	for(TSize i = 0; i<nseq; ++i) {
		_streamWrite(file,"\tfasta_id=\"");
		_streamWrite(file,names[i]);
		_streamWrite(file,"\" sequence=\"");
		_streamWrite(file,str[i]);
		_streamWrite(file,"\" track=");
		_streamPutInt(file, i);
		_streamWrite(file," type=\"");
		_streamWrite(file,"DNA");
		_streamWrite(file,"\": ");
		_streamPutInt(file, 0);
		_streamPut(file, ' ');
		_streamPutInt(file, length(str[i])- 1);
		_streamPut(file, '\n');
	}
	_streamWrite(file,"}\n");
	for(TSize i = 0; i<nseq; ++i) {
		_streamWrite(file,"{DATA ");
		_streamPutInt(file, i);
		//_streamWrite(file,names[i]);
		_streamWrite(file,"-seqlen\n");
		_streamWrite(file,"\t[__GLOBAL__]\n");
		_streamWrite(file,"\tlength=");
		_streamPutInt(file, length(str[i]));
		_streamWrite(file,":\t");
		_streamPutInt(file, 0);
		_streamPut(file, ' ');
		_streamPutInt(file, length(str[i])- 1);
		_streamPut(file, '\n');
		_streamWrite(file,"}\n");
	}
	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			_streamWrite(file,"{DATA ");
			_streamPutInt(file, i);
			//_streamWrite(file,names[i]);
			_streamWrite(file,"-vs-");
			_streamPutInt(file, j);
			//_streamWrite(file,names[j]);
			_streamPut(file, '\n');
			_streamWrite(file,"\t[__GLOBAL__]\n");
			for(TIterConst it = begin(g.data_align.data_vertex);!atEnd(it);goNext(it)) {
				TVertexDescriptor sourceV = position(it);
				TId id1 = sequenceId(g, sourceV);
				if ((positionToId(str, id1) != i) &&
					(positionToId(str, id1) != j)) continue;
				TEdgeStump* current = getValue(it);
				while(current!=0) {
					TVertexDescriptor targetV = getTarget(current);
					TId id2 = sequenceId(g, targetV);
					if (sourceV != targetV) {
						if ((positionToId(str, id2) != i) &&
							(positionToId(str, id2) != j)) {
								current=getNextS(current);
								continue;
						}
						_streamWrite(file,"\t");
						_streamWrite(file,"source=");
						_streamPutInt(file, sourceV);		
						_streamPut(file, ' ');
						_streamWrite(file,"target=");
						_streamPutInt(file, targetV);
						_streamPut(file, ' ');
						_streamWrite(file,"edgeId=");
						_streamPutInt(file, _getId(current));
						_streamPut(file, ' ');
						_streamWrite(file,"cargo=");
						__writeCargo(file,g,current);
						_streamPut(file, ' ');
						_streamWrite(file,"label=");
						_streamWrite(file,label(g,sourceV));
						_streamPut(file, ' ');
						_streamWrite(file,"labelOpp=");
						_streamWrite(file,label(g,targetV));
						_streamPut(file, ':');
						_streamPut(file, '\t');
						_streamPutInt(file, fragmentBegin(g, sourceV));
						_streamPut(file, ' ');
						_streamPutInt(file, fragmentBegin(g, targetV));
						_streamPut(file, ' ');
						_streamPutInt(file, fragmentBegin(g, sourceV) + fragmentLength(g, sourceV));
						_streamPut(file, ' ');
						_streamPutInt(file, fragmentBegin(g, targetV) + fragmentLength(g, targetV));
						_streamPut(file, '\n');
						current=getNextS(current);
					} else {
						current=getNextT(current);
					}
				}
			}
			_streamWrite(file,"}\n");	
		}
	}
}

//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TAlignmentMatrix, typename TOldBegEndPos, typename TReadBegEndPos, typename TGappedConsensus>
inline void
write(TFile & file,
	  Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
	  TAlignmentMatrix const& mat,
	  TOldBegEndPos const& oldBegEndPos,
	  TReadBegEndPos const& readBegEndPos,
	  TGappedConsensus const& gappedConsensus,
	  FastaReadFormat) 
{
	typedef typename Size<TAlignmentMatrix>::Type TSize;
	typedef typename Value<TAlignmentMatrix>::Type TValue;

	// Initialization
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	TSize len = length(gappedConsensus);
	TSize maxCoverage = length(mat) / len;
	TValue gapChar = gapValue<TValue>();
	
	// Print the alignment matrix
	TSize winSize = 60;
	int offset = 2;
	TSize column = 0;
	while (column<len) {
		TSize window_end = column + winSize;
		if (window_end >= len) window_end = len;
		// Position
		for(int i = 0; i<offset - 2; ++i) _streamPut(file,' ');
		_streamWrite(file,"Pos: ");
		_streamPutInt(file,column);
		_streamPut(file,'\n');
		// Ruler
		for(int i = 0; i<offset + 3; ++i) _streamPut(file,' ');
		for(TSize local_col = 1; local_col<window_end - column + 1; ++local_col) {
			if ((local_col % 10)==0) _streamPut(file, ':');
			else if ((local_col % 5)==0) _streamPut(file, '.');
			else _streamPut(file, ' ');
		}
		_streamPut(file,'\n');
		// Reads
		for(TSize row = 0; row<maxCoverage; ++row) {
			TSize tmp = row;
			int off = 0;
			while (tmp / 10 != 0) {
				tmp /= 10;
				++off;
			}
			for(int i = 0; i<offset - off; ++i) _streamPut(file,' ');
			_streamPutInt(file, row);
			_streamPut(file,':');
			_streamPut(file,' ');
			for(TSize local_col = column; local_col<window_end; ++local_col) {
				_streamPut(file, getValue(mat, row*len+local_col));
			}
			_streamPut(file,'\n');
		}
		_streamPut(file,'\n');
		// Consensus
		for(int i = 0; i<offset; ++i) _streamPut(file,' ');
		_streamWrite(file,"C: ");
		for(unsigned int local_col = column; local_col<window_end; ++local_col) _streamPut(file, gappedConsensus[local_col]);
		_streamPut(file,'\n');
		_streamPut(file,'\n');
		column+=winSize;
	}
	_streamPut(file,'\n');
	_streamPut(file,'\n');

	// Print all reads
	for(TSize i = 0; i<nseq; ++i) {
		_streamWrite(file,"typ:R");
		_streamPutInt(file, i);
		_streamPut(file,'\n');
		_streamWrite(file,"seq:");
		if ((oldBegEndPos[i]).i1 > (oldBegEndPos[i]).i2) reverseComplementInPlace(str[i]);
		_streamWrite(file,str[i]);
		_streamPut(file,'\n');
		_streamWrite(file,"Pos:");
		if ((oldBegEndPos[i]).i1 > (oldBegEndPos[i]).i2) {
			_streamPutInt(file, (readBegEndPos[i]).i2);
			_streamPut(file,',');
			_streamPutInt(file, (readBegEndPos[i]).i1);
		} else {
			_streamPutInt(file, (readBegEndPos[i]).i1);
			_streamPut(file,',');
			_streamPutInt(file, (readBegEndPos[i]).i2);
		}
		_streamPut(file,'\n');
		
		std::stringstream gapCoords;
		TSize letterCount = 0;
		TSize gapCount = 0;
		for(TSize column = (readBegEndPos[i]).i1; column<(readBegEndPos[i]).i2; ++column) {
			if (value(mat, (readBegEndPos[i]).i3 * len + column) == gapChar) {
				++gapCount;
				gapCoords << letterCount << ' ';
			} else ++letterCount;
		}
		_streamWrite(file,"dln:");
		_streamPutInt(file, gapCount);
		_streamPut(file,'\n');
		_streamWrite(file,"del:");
		_streamWrite(file, gapCoords.str().c_str());
		_streamPut(file,'\n');
		_streamPut(file,'\n');
	}
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignStringSet:
..cat:Graph.Alignment Graph
..summary:Assigns a new string set to an alignment graph.
..signature:assignStringSet(g, str)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..param.str:A string set.
..see:Function.getStringSet
..see:Function.stringSet
*/
template<typename TString, typename TDefault, typename TCargo, typename TSpec, typename TDefault2>
inline void
assignStringSet(Graph<Alignment<StringSet<TString, Dependent<TDefault> >, TCargo, TSpec> >& g,
				StringSet<TString, Dependent<TDefault2> >& sStr)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<StringSet<TString, Dependent<TDefault> >, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	clear(g);
	g.data_sequence = (StringSet<TString, Dependent<TDefault> >) sStr;
	for(unsigned k=0; k<length(sStr);++k) {
		g.data_pvMap.insert(std::make_pair(std::make_pair(positionToId(sStr,k), length(sStr[k])), nilVertex));
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TDefault, typename TCargo, typename TSpec, typename TDefault2>
inline void
assignStringSet(Graph<Alignment<StringSet<TString, Dependent<TDefault> >, TCargo, TSpec> >& g,
				StringSet<TString, Owner<TDefault2> >& sStr)
{
	SEQAN_CHECKPOINT
	StringSet<TString, Dependent<> > depStr(sStr);
	assignStringSet(g, depStr);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.getStringSet:
..cat:Graph.Alignment Graph
..summary:Gets the string set of an alignment graph.
..signature:getStringSet(g)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..returns:A string set.
..see:Function.assignStringSet
..see:Function.stringSet
*/
template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Host<Graph<Alignment<TStringSet, TCargo, TSpec> > const>::Type&
getStringSet(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return value(g.data_sequence);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.stringSet:
..cat:Graph.Alignment Graph
..summary:Gets the string set of an alignment graph.
..signature:stringSet(g)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..returns:A reference to a string set.
..see:Function.assignStringSet
..see:Function.getStringSet
*/
template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Host<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
stringSet(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return const_cast<TStringSet&>(value(g.data_sequence));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.label:
..cat:Graph.Alignment Graph
..summary:Gets the label that is associated with this vertex descriptor.
..signature:label(g, v)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:The label.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Infix<typename Value<TStringSet>::Type>::Type
label(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
	  TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef FragmentInfo<TIdType, TSize> TSegment;
	TSegment seg = getProperty(g.data_fragment, v);
	//std::cout << seg.data_seq_id << ",";
	//std::cout << seg.data_begin << ",";
	//std::cout << seg.data_length << ",";
	//std::cout << getValueById(value(g.data_sequence), seg.data_seq_id) << std::endl;
	//std::cout << infix(getValueById(value(g.data_sequence), seg.data_seq_id), seg.data_begin, seg.data_begin + seg.data_length) << std::endl;
	return infix(getValueById(value(g.data_sequence), seg.data_seq_id), seg.data_begin, seg.data_begin + seg.data_length);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.sequenceId:
..cat:Graph.Alignment Graph
..summary:Gets the sequence id that is associated with this vertex descriptor.
..signature:sequenceId(g, v)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:The sequence id.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
sequenceId(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_seq_id);
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.fragmentBegin:
..cat:Graph.Alignment Graph
..summary:Gets the begin position for this vertex descriptor in the sequence.
..signature:fragmentBegin(g, v)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:The begin position.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
fragmentBegin(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			  TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_begin);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.fragmentLength:
..cat:Graph.Alignment Graph
..summary:Gets the length of the label of a given vertex descriptor in the sequence.
..signature:fragmentLength(g, v)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:The length of the fragment represented by this vertex descriptor.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
fragmentLength(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			   TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_length);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.findVertex:
..cat:Graph.Alignment Graph
..summary:Finds a vertex given a sequence id and a position.
..signature:fragmentLength(g, v)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:The length of the fragment represented by this vertex descriptor.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId, typename TPos> 
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
findVertex(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
		   TSeqId id,
		   TPos pos)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	if (pos >= (TPos) length(getValueById(stringSet(g),id))) return getNil<TVertexDescriptor>();
	else return g.data_pvMap.upper_bound(std::make_pair(id, pos))->second;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId, typename TPosition, typename TSeqId2, typename TPosition2> 
inline void
getProjectedPosition(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					 TSeqId const id1,
					 TPosition const pos1,
					 TSeqId2& id2,
					 TPosition2& pos2)
{
	SEQAN_CHECKPOINT
	SEQAN_TASSERT(length(stringSet(g)) == 2);

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	TStringSet& str = stringSet(g);
	TVertexDescriptor sV = findVertex(g, id1, pos1);

	// Case 1: No projection possible
	if (sV == getNil<TVertexDescriptor>()) {
		if ( (TSeqId) positionToId(str, 0) == id1) id2 = (TSeqId2) positionToId(str,1);
		else id2 = (TSeqId2) positionToId(str,0);
		pos2 = 0;
		return;
	}

	// Case 2: Projection is possible
	TEdgeStump* current = getValue(g.data_align.data_vertex, sV);
	if(current != (TEdgeStump*) 0) {
		TVertexDescriptor tV = target(current);
		if (tV == sV) tV = source(current);
		pos2 = (TPosition2) (fragmentBegin(g,tV) + (pos1 - fragmentBegin(g, sV)));
		id2 = (TSeqId2) sequenceId(g, tV);
		return;
	} else {
		// If no out-going edge, get the preceding or following vertex
		if (fragmentBegin(g, sV) == 0) {
			getProjectedPosition(g, id1, fragmentBegin(g,sV) + fragmentLength(g, sV), id2, pos2);
			return;
		} else {
			getProjectedPosition(g, id1, fragmentBegin(g,sV) - 1, id2, pos2);
			return;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getFirstCoveredPosition:
..cat:Graph.Alignment Graph
..summary:Finds the first position in a sequence that is not assigned to a nil vertex.
..signature:getFirstCoveredPosition(g, id)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..param.id:A sequence id.
..returns:A sequence position
..see:Function.getLastCoveredPosition
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId> 
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
getFirstCoveredPosition(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						TSeqId const id)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;


	typename TPosToVertexMap::const_iterator it = g.data_pvMap.upper_bound(std::make_pair(id, 0));

	// Case 1: id is not covered
	if (it == g.data_pvMap.end()) return length(getValueById(stringSet(g), id));

	// Case 2: We found a nil vertex, go one forward
	if (it->second == getNil<TVertexDescriptor>()) {
		++it;
		if (it == g.data_pvMap.end()) return length(getValueById(stringSet(g), id));
	}

	// Now we have the right vertex, return the beginning if the sequence id still fits
	if (it->first.first != id) return length(getValueById(stringSet(g), id));
	else return fragmentBegin(g, it->second);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getLastCoveredPosition:
..cat:Graph.Alignment Graph
..summary:Finds the last position in a sequence that is not assigned to a nil vertex.
..signature:getLastCoveredPosition(g, id)
..param.g:An alignment graph.
...type:Spec.Alignment Graph
..param.id:A sequence id.
..returns:A sequence position
..see:Function.getFirstCoveredPosition
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId> 
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
getLastCoveredPosition(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TSeqId id)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;

	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	typename TPosToVertexMap::const_iterator it = g.data_pvMap.lower_bound(std::make_pair(id, length(getValueById(stringSet(g), id))));

	// Case 1: No last position (all nil)
	if ((it == g.data_pvMap.begin()) && (it->second == nilVertex)) return 0;

	// Case 2: Found a nil position but there is a vertex before
	if (it->second == nilVertex) {
		--it;
	}

	// If the sequence id still matches return the position behind the last position belonging to this vertex
	if (it->first.first != id) return 0;
	return fragmentBegin(g, it->second) + fragmentLength(g, it->second);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
