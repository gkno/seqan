#ifndef SEQAN_HEADER_GRAPH_IMPL_ALIGN_H
#define SEQAN_HEADER_GRAPH_IMPL_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TSize>
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
			SEQAN_CHECKPOINT
		}
};



//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Undirected<TCargo, TSpec> > >::Type*> const&
_getVertexString(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g) {
	return g.data_align.data_vertex;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline String<typename EdgeType<Graph<Undirected<TCargo, TSpec> > >::Type*>&
_getVertexString(Graph<Alignment<TStringSet, TCargo, TSpec> >& g) {
	return g.data_align.data_vertex;
}

/////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type, Default> const &
_getVertexIdManager(Graph<Alignment<TStringSet,TCargo, TSpec> > const& g) {
	return g.data_align.data_id_managerV;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec>
inline IdManager<typename Id<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type, Default> &
_getVertexIdManager(Graph<Alignment<TStringSet, TCargo, TSpec> >& g) {
	return g.data_align.data_id_managerV;
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TVertexDescriptor> 
inline typename EdgeDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
addEdge(Graph<Alignment<TStringSet, TCargo, TSpec> >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TCargo const cargo) 
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


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
