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
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef typename Size<Graph>::Type TSize;
		typedef typename Size<TStringSet>::Type TStringSetSize;
		typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;

		// Alignment graph
		Graph<Undirected<TCargo, TSpec> > data_align;

		// Sequences
		Holder<TStringSet> data_sequence;
		
		// Alignment specific members
		String<SegmentInfo<TIdType, TSize> > data_segment;

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
		   bool transpose) 
{
	SEQAN_CHECKPOINT
	clear(dest);
	dest.data_align = source.data_align;
	dest.data_sequence = source.data_sequence;
	dest.data_segment = source.data_segment;
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
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	clearVertices(g.data_align);
	clear(g.data_segment);
	g.data_pvMap.clear();
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	if(!empty(value(g.data_sequence))) {
		for(TSize k=0; k<length(stringSet(g));++k) {
			g.data_pvMap.insert(std::make_pair(std::make_pair(positionToId(stringSet(g),k), length(stringSet(g)[0])), nilVertex));
		}
	}
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
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Size<TGraph>::Type TSize;
	typedef SegmentInfo<TIdType, TSize> TSegmentInfo;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;

	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Store the new segment
	typename TPosToVertexMap::iterator interval = g.data_pvMap.lower_bound(std::make_pair(id, begin + len));
	// Segment does not belong to Sequence anymore
	SEQAN_TASSERT(interval != g.data_pvMap.end());
	// Segment end must be assigned to nil so far
	SEQAN_TASSERT(interval->second == nilVertex);
	// Segment must belong to the whole old interval
	SEQAN_TASSERT(interval == g.data_pvMap.upper_bound(std::make_pair(id, begin)));

	// Insert new vertex
	TVertexDescriptor vd = addVertex(g.data_align);
	if (length(g.data_segment) <= vd) resize(g.data_segment, vd + 1, Generous());
	assignProperty(g.data_segment, vd, TSegmentInfo(id, begin, len));

	// Update position to vertex map
	// Does the end of the new segment coincides with the end of the interval?
	if ( (TSize) begin + len == (TSize) interval->first.second) {
		// Does the beginning of the new segment coincides with the beginning of the interval?
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
		// Does the beginning of the new segment coincides with the beginning of the interval?
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
	typename TPosToVertexMap::iterator interval = g.data_pvMap.lower_bound(std::make_pair(sequenceId(g,v), segmentBegin(g,v) + segmentLength(g,v)));
	SEQAN_TASSERT(interval != g.data_pvMap.end());
	interval->second = getNil<TVertexDescriptor>();

	// Remove the vertex
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

	Matrix<char> align;
	if (convertAlignment(g, align)) {
		_streamWrite(target,"Alignment matrix:\n");
		TSize colLen = length(align, 0);
		TSize nseq = length(align, 1);
		for (TSize row=0;row<nseq;++row) {
			for(TSize col=0;col<colLen;++col) {
				_streamPut(target, getValue(align, row*colLen+col));
			}
			_streamPut(target, '\n');
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

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TSeqId, typename TPos> 
inline typename VertexDescriptor<Graph<Alignment<TStringSet, TCargo, TSpec> > >::Type 
findVertex(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
		   TSeqId id,
		   TPos pos)
{
	SEQAN_CHECKPOINT
	return g.data_pvMap.upper_bound(std::make_pair(id, pos))->second;
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
	typedef SegmentInfo<TIdType, TSize> TSegmentInfo;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;
	typedef std::map<unsigned int, unsigned int> TComponentLength;

	// Strongly Connected Components
	String<unsigned int> component;
	strongly_connected_components(g, component);

	typedef std::list<unsigned int> TComps;
	typedef std::set<unsigned int> TUniqueComps;
	TComps comps;
	TUniqueComps bag;
	TComponentLength compLength;

	// Walk through the first sequence
	typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();
	TIdType firstSeq = it->first.first;
	for(; it != g.data_pvMap.end(); ++it) {
		if (it->first.first != firstSeq) break;
		unsigned int c = getProperty(component, it->second);
		if (!bag.insert(c).second) return false;
		compLength.insert(std::make_pair(c, segmentLength(g, it->second)));
		comps.push_back(c);
	}

	// Walk through all other sequences
	TComps::iterator pos = comps.begin();
	TIdType currentSeq = it->first.first;
	for(; it != g.data_pvMap.end(); ++it) {
		if (it->first.first != currentSeq) {
			pos = comps.begin();
			currentSeq = it->first.first;
		}
		unsigned int c = getProperty(component, it->second);
		// Have we seen this component before?
		if (bag.find(c) != bag.end()) {
			while ((pos != comps.end()) && (*pos != c)) ++ pos;
			// Crossing components
			if (pos == comps.end()) return false;
		} else {
			bag.insert(c);
			compLength.insert(std::make_pair(c, segmentLength(g, it->second)));
			pos = comps.insert(pos, c);
		}
		++pos;
	}

	// Create the matrix
	TSize len = 0;
	TSize nseq = length(stringSet(g));
	for(TComponentLength::iterator cIt=compLength.begin(); cIt != compLength.end(); ++cIt) len+=cIt->second;
	setDimension(mat, 2);setLength(mat, 0, len);setLength(mat, 1, nseq);
	resize(mat);

	// Fill the matrix
	TSize row = 0;
	TSize col = 0;
	it = g.data_pvMap.begin();
	pos = comps.begin();
	currentSeq = it->first.first;
	for(; it != g.data_pvMap.end(); ++it) {
		if (it->first.first != currentSeq) {
			SEQAN_TASSERT(col == len);
			//std::cout << std::endl;
			++row;col=0;
			pos = comps.begin();
			currentSeq = it->first.first;
		}
		unsigned int c = getProperty(component, it->second);
		while ((pos != comps.end()) && (*pos != c)) {
			for(TSize i=0;i<compLength[*pos];++i) {
				//std::cout << '-';
				assignValue(mat, row*len + col, '-');
				++col;
			}
			++pos;
		}
		String<TValue> str = label(g,it->second);
		//std::cout << str;
		for(TSize i=0;i<length(str);++i) {
			assignValue(mat, row*len + col, (TValue) getValue(str, i));
			++col;
		}

		++pos;
	}
	SEQAN_TASSERT(row + 1 == nseq);
	//std::cout << std::endl;

	return true;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
