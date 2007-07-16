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
..include:graph.h
*/
template<typename TStringSet, typename TCargo, typename TSpec>
class Graph<Alignment<TStringSet, TCargo, TSpec> > 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef typename Size<Graph>::Type TSize;
		typedef typename Size<TStringSet>::Type TStringSetSize;
		typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;

		// Alignment graph
		Graph<Undirected<TCargo, TSpec> > data_align;

		// Sequences
		Holder<TStringSet> data_sequence;
		
		// Alignment specific members
		String<FragmentInfo<TIdType, TSize> > data_fragment;

		// STL Map to retrieve a vertex given SeqId, Position
		TPosToVertexMap data_pvMap;


//____________________________________________________________________________

		Graph() {
		}


		Graph(TStringSet sSet) {
			SEQAN_CHECKPOINT
			data_sequence = sSet;
			TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
			for(TStringSetSize k=0; k<length(sSet);++k) {
				data_pvMap.insert(std::make_pair(std::make_pair(positionToId(sSet,k), length(sSet[k])), nilVertex));
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

	clearVertices(g.data_align);
	clear(g.data_fragment);
	g.data_pvMap.clear();
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
void 
removeEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
		   TVertexDescriptor const source, 
		   TVertexDescriptor const target) 
{
	SEQAN_CHECKPOINT
	removeEdge(g.data_align, source, target);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
void 
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
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;


	
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
			for(unsigned int j = 0;j<leftSpace-offset;++j) {
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TStringSet2>
inline void
assignStringSet(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				TStringSet2& sStr)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	clear(g);
	g.data_sequence = (TStringSet) sStr;
	for(unsigned k=0; k<length(sStr);++k) {
		g.data_pvMap.insert(std::make_pair(std::make_pair(positionToId(sStr,k), length(sStr[k])), nilVertex));
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Host<Graph<Alignment<TStringSet, TCargo, TSpec> > const>::Type&
getStringSet(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return value(g.data_sequence);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline typename Host<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
stringSet(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return const_cast<TStringSet&>(value(g.data_sequence));
}

//////////////////////////////////////////////////////////////////////////////

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
	//std::cout << infix(getString(value(g.data_sequence), seg.data_seq_id), seg.data_begin, seg.data_begin + seg.data_length) << std::endl;
	return infix(getValueById(value(g.data_sequence), seg.data_seq_id), seg.data_begin, seg.data_begin + seg.data_length);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
sequenceId(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_seq_id);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
fragmentBegin(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			  TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_begin);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
fragmentLength(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			   TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_fragment[v].data_length);
}

//////////////////////////////////////////////////////////////////////////////

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
	if (sV == getNil<TVertexDescriptor>()) {
		if ( (TSeqId) positionToId(str, 0) == id1) id2 = (TSeqId2) positionToId(str,1);
		else id2 = (TSeqId2) positionToId(str,0);
		pos2 = 0;
		return;
	}

	TEdgeStump* current = getValue(g.data_align.data_vertex, sV);
	if(current != (TEdgeStump*) 0) {
		TVertexDescriptor tV = target(current);
		if (tV == sV) tV = source(current);
		pos2 = (TPosition2) (fragmentBegin(g,tV) + (pos1 - fragmentBegin(g, sV)));
		id2 = (TSeqId2) sequenceId(g, tV);
		return;
	} else {
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId> 
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
getFirstCoveredPosition(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
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
	typename TPosToVertexMap::const_iterator it = g.data_pvMap.upper_bound(std::make_pair(id, 0));

	if (it == g.data_pvMap.end()) return length(getValueById(stringSet(g), id));
	if (it->second == nilVertex) {
		++it;
		if (it == g.data_pvMap.end()) return length(getValueById(stringSet(g), id));
	}
	if (it->first.first != id) return length(getValueById(stringSet(g), id));
	else return fragmentBegin(g, it->second);
}

//////////////////////////////////////////////////////////////////////////////

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

	if ((it == g.data_pvMap.begin()) && (it->second == nilVertex)) return 0;
	if (it->second == nilVertex) {
		--it;
	}
	if (it->first.first != id) return 0;
	return fragmentBegin(g, it->second) + fragmentLength(g, it->second);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix> 
inline bool
convertAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				 TMatrix& mat)
{
	SEQAN_CHECKPOINT
	if (empty(g)) return false;

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef FragmentInfo<TIdType, TSize> TFragmentInfo;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;
	typedef std::map<unsigned int, unsigned int> TComponentLength;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Strongly Connected Components
	String<unsigned int> component;
	strongly_connected_components(g, component);

	// Make a directed graph to represent the ordering of the components
	Graph<Directed<> > componentGraph;
	for(TSize i = 0; i<length(component);++i) addVertex(componentGraph);

	// Walk through all sequences and add edges
	typename TPosToVertexMap::const_iterator it1 = g.data_pvMap.begin();
	while (it1!=g.data_pvMap.end()) {
		// If sections are not assigned to a vertex -> no alignment
		if (it1->second == nilVertex) return false;
		typename TPosToVertexMap::const_iterator it2 = it1;
		++it2;
		TIdType currentSeq = it1->first.first;
		while ((it2!=g.data_pvMap.end()) && 
			(it2->first.first == currentSeq)) {
			if (it2->second == nilVertex) return false;
			unsigned int c1 = getProperty(component, it1->second);
			unsigned int c2 = getProperty(component, it2->second);
			// If two components appear twice in the same sequence -> no alignment
			if (c1 == c2) return false;
			else {
				if (findEdge(componentGraph, c1, c2) == 0) addEdge(componentGraph, c1, c2);
			}
			++it2;
		}
		++it1;
	}

	// Make a topological sort of the component graph
	String<unsigned int> order;
	topological_sort(componentGraph, order);

	// Walk through all sequences and check the component order
	// Also store the length of each component
	TComponentLength compLength;
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

	
	// Create the matrix
	TSize len = 0;
	TSize nseq = length(stringSet(g));
	for(TComponentLength::iterator cIt=compLength.begin(); cIt != compLength.end(); ++cIt) len+=cIt->second;
	char gapChar = gapValue<char>();
	fill(mat, len * nseq, gapChar);

	// Fill the matrix
	TSize row = 0;
	TSize col = 0;
	it = g.data_pvMap.begin();
	compIndex = 0;
	currentSeq = it->first.first;
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


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TLibraries> 
inline void
combineGraphs(Graph<Alignment<TStringSet, TCargo, TSpec> >& outGraph,
			  TLibraries& libs)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

	// Clear out library
	clearVertices(outGraph);
	TSize numLibs = length(libs);

	// All the matches with score values
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TCargo, Block<> > score_values;

	// Max score and index start position for every library
	String<TCargo> max_scores;
	String<TCargo> index_start;
	resize(max_scores, numLibs);
	resize(index_start, numLibs);

	// Get all matches
	TSize count = 0;
	for(TSize i = 0; i<numLibs; ++i) {
		assignValue(index_start, i, count);
		TCargo maxCargoLib = 0;
		TGraph* lib = getValue(libs, i);
		TEdgeIterator it(*lib);
		for(;!atEnd(it);++it) {
			if (getCargo(*it) > maxCargoLib) maxCargoLib = getCargo(*it);
			TVertexDescriptor sV = sourceVertex(it);
			TVertexDescriptor tV = targetVertex(it);
			push_back(matches, TFragment( (unsigned int) sequenceId(*lib, sV), (unsigned int) fragmentBegin(*lib,sV), (unsigned int) sequenceId(*lib, tV),  (unsigned int)  fragmentBegin(*lib,tV),  (unsigned int)  fragmentLength(*lib,tV)));
			push_back(score_values, cargo(*it));
			++count;
		}
		assignValue(max_scores, i, maxCargoLib);
	}

	// Match refinement
	TStringSet& str = stringSet(outGraph);	
	matchRefinement(matches,str,outGraph);

	// Adapt edge weights
	count = 0;
	TSize currentLib = 0;
	TFragmentStringIter endIt = end(matches);
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
		if ((currentLib < numLibs - 1) && (count >= (TSize) getValue(index_start, currentLib+1))) ++currentLib;
		double scaling = (double) 100 / (double) getValue(max_scores, currentLib);
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			SEQAN_TASSERT(pos2 < pos2 + fragmentLength(*it, id2))
			TVertexDescriptor p1 = findVertex(outGraph, id1, pos1);
			TVertexDescriptor p2 = findVertex(outGraph, id2, pos2);
			TEdgeDescriptor e = findEdge(outGraph, p1, p2);
			TSize fragLen = fragmentLength(outGraph, p1); 
			double newVal = (double) fragLen / (double) (end1 - pos1);
			newVal *= scaling;
			newVal *= (double) getValue(score_values, position(it));
			if (e != 0) cargo(e) += (TCargo) newVal;
			else addEdge(outGraph, p1, p2, (TCargo) newVal);
			SEQAN_TASSERT(fragLen == fragmentLength(outGraph, p2))
			pos1 += fragLen;
			pos2 += fragLen;
		}
		++count;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrixAlign, typename TMatrix> 
inline void
_mutualInformationContent(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						  TMatrixAlign const& align,
						  TMatrix& mat)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Value<TMatrixAlign>::Type TChar;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;
				
	TSize nseq = length(stringSet(g));
	TSize colLen = length(align) / nseq;
	TSize value_size = ValueSize<TAlphabet>::VALUE;
	TChar gapChar = gapValue<TChar>();
	

	resize(mat, colLen * colLen);
	for(TSize i = 0; i < colLen; ++i) {
		for(TSize j = i + 1; j < colLen; ++j) {
			String<TValue> f_i;
			fill(f_i, value_size, 0);
			String<TValue> f_j;
			fill(f_j, value_size, 0);
			String<TValue> f_ij;
			fill(f_ij, value_size * value_size, 0);
			TSize count = 0;
			for(TSize seq = 0; seq < nseq; ++seq) {
				if ((getValue(align, seq * colLen + i) != gapChar) && (getValue(align, seq * colLen + j) != gapChar)) {
					TAlphabet c1 = (TAlphabet) getValue(align, seq * colLen + i);
					TAlphabet c2 = (TAlphabet) getValue(align, seq * colLen + j);
					value(f_i, (Byte) c1) += 1;
					value(f_j, (Byte) c2) += 1;
					value(f_ij, ((Byte) c1) * value_size + (Byte) c2) += 1;
					++count;
				}
			}
			// Calculate h_ij
			TValue h_ij = 0;
			if (count != 0) {
				for(TSize x = 0; x < value_size;++x) {
					value(f_i, x) /= count;
					//std::cout << getValue(f_i, x) << std::endl;
					for(TSize y = 0; y < value_size;++y) {
						if (x == 0) {
							value(f_j, y) /= count;
							//std::cout << getValue(f_j, y) << std::endl;
						}
						value(f_ij, x * value_size + y) /= count;
						//std::cout << getValue(f_ij, x * value_size + y) << std::endl;
						if (getValue(f_ij, x * value_size + y) != 0) {
							TValue ratio = getValue(f_ij, x * value_size + y) / (getValue(f_i, x) * getValue(f_j, y));
							h_ij += getValue(f_ij, x * value_size + y) * (log((TValue) ratio) / log((TValue) 2));						}
					}
				}
			}
			assignValue(mat, i * colLen + j, h_ij);
			assignValue(mat, j * colLen + i, h_ij);
			//std::cout << i << ',' << j << ':' << h_ij << std::endl;
		}
	}

	//// Debug code
	//for(TSize i = 0; i < colLen; ++i) {
	//	for(TSize j = i+1; j < colLen; ++j) {
	//		std::cout << i << ',' << j << ':' << getValue(mat, i * colLen + j) << std::endl;
	//	}
	//}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix> 
inline void
mutualInformationContent(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						 TMatrix& mat)
{
	SEQAN_CHECKPOINT	
	String<char> align;
	if (convertAlignment(g, align)) {
		_mutualInformationContent(g, align, mat);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
