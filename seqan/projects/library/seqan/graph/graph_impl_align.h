#ifndef SEQAN_HEADER_GRAPH_IMPL_ALIGN_H
#define SEQAN_HEADER_GRAPH_IMPL_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TId = unsigned int, typename TSize = unsigned int>
class SegmentInfo {
public:
	TId data_seq_id;
	TSize data_begin;
	TSize data_length;

	SegmentInfo() :
		data_seq_id(0),
		data_begin(0),
		data_length(0)
	{
	}

	SegmentInfo(TId id, TSize beg, TSize len) :
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
		typedef typename Size<Graph>::Type TSize;

		// Alignment graph
		Graph<Undirected<TCargo, TSpec> > data_align;

		// Sequences
		Holder<TStringSet> data_sequence;
		
		// Alignment specific members
		String<SegmentInfo<TIdType, TSize> > data_segment;

//____________________________________________________________________________

		Graph(TStringSet sSet) {
			SEQAN_CHECKPOINT
			data_sequence = sSet;
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

private:
		Graph() {
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
		   bool transpose) 
{
	SEQAN_CHECKPOINT
	clear(dest);
	dest.data_align = source.data_align;
	dest.data_sequence = source.data_sequence;
	dest.data_segment = source.data_segment;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Alignment<TStringSet, TCargo, TSpec> > const& source,
		   Graph<Alignment<TStringSet, TCargo, TSpec> >& dest) 
{
	_copyGraph(source, dest, false); // Never transpose
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
transpose(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
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
	clearVertices(g.data_align);
	clear(g.data_segment);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline void 
clear(Graph<Alignment<TStringSet, TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
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
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef SegmentInfo<TIdType, TSize> TSegmentInfo;

	TVertexDescriptor vd = addVertex(g.data_align);
	if (length(g.data_segment) <= vd) resize(g.data_segment, vd + 1, Generous());
	assignProperty(g.data_segment, vd, TSegmentInfo(id, begin, len));
	return vd;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	removeVertex(g.data_align,v);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
addEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
		TVertexDescriptor source, 
		TVertexDescriptor target) 
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
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef SegmentInfo<TIdType, TSize> TSegment;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdgeStump;
	typedef typename Iterator<String<TEdgeStump*> const>::Type TIterConst;
	_streamWrite(target,"Adjacency list:\n");
	for(TIterConst it = begin(g.data_align.data_vertex);!atEnd(it);goNext(it)) {
		TVertexDescriptor sourceV = position(it);
		_streamPutInt(target, sourceV);
		TSegment seg = getProperty(g.data_segment, sourceV);
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
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TStringSet2>
inline void
assignStringSet(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				TStringSet2& sStr)
{
	SEQAN_CHECKPOINT
	clear(g);
	g.data_sequence = (TStringSet) sStr;
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
inline typename Value<TStringSet>::Type
label(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
	  TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef SegmentInfo<TIdType, TSize> TSegment;
	TSegment seg = getProperty(g.data_segment, v);
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
	return const_cast<typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_segment[v].data_seq_id);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
segmentBegin(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			 TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Position<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_segment[v].data_begin);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&
segmentLength(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			  TVertexDescriptor const v)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Size<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type&>(g.data_segment[v].data_length);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
