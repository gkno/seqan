#ifndef SEQAN_HEADER_GRAPH_IMPL_AUTOMATON_H
#define SEQAN_HEADER_GRAPH_IMPL_AUTOMATON_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TEdge, typename TAlphabet>
class AutomatonEdgeArray {
public:
	TEdge data_edge[ValueSize<TAlphabet>::VALUE];

	AutomatonEdgeArray() 
	{
		typedef typename VertexDescriptor<TEdge>::Type TVertexDescriptor;
		typedef typename Size<TAlphabet>::Type TSize;
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
		for(TSize i=0;i<table_length;++i) {
			assignTarget(&data_edge[i], nilVal);
		}
	}
};

////////////////
// Automaton
////////////////
/**
.Spec.Automaton:
..cat:Graph
..general:Class.Graph
..summary:A rooted graph with edge labels.
..signature:Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>
..param.TAlphabet:The alphabet type that is used for the edge labels.
...metafunction:Metafunction.Alphabet
...remarks:Use Alphabet<Automaton<TAlphabet, TCargo, TEdgeSpec> > to get the alphabet type of an automaton.
...default:$char$
..param.TCargo:The cargo type that can be attached to edges.
...metafunction:Metafunction.Cargo
...remarks:Use Alphabet<Automaton<TAlphabet, TCargo, TEdgeSpec> >@ to get the cargo type of an automaton.
...default:$void$
..param.TEdgeSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..include:graph.h
*/
template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
class Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdge;
		typedef typename VertexDescriptor<Graph>::Type TVertexDescriptor;
		typedef typename IdHandler<TEdge, TIdType>::Type TEdgeIdManager;

		String<AutomatonEdgeArray<TEdge, TAlphabet> > data_vertex;		// List of tables
		IdManager<TIdType> data_id_managerV;
		TEdgeIdManager data_id_managerE;
		TVertexDescriptor root;
	

//____________________________________________________________________________


		Graph() : root(0) {
			SEQAN_CHECKPOINT
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
			_copyGraph(_other, *this);
			return *this;
		}
};

//////////////////////////////////////////////////////////////////////////////
// Automaton specific graph functions
//////////////////////////////////////////////////////////////////////////////
template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& dest,
		   bool transpose)
{
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Cargo<TEdge>::Type TCargoType;
	typedef typename Size<TAlphabet>::Type TSize;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	
	clear(dest);
	resize(dest.data_vertex, length(source.data_vertex));
	dest.root = source.root;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	for(TIterConst it = begin(source.data_vertex);!atEnd(it);goNext(it)) {
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TVertexDescriptor sourceVertex = position(it);
		for(TSize i=0;i<table_length;++i) {
			TEdgeDescriptor const edSource = (TEdgeDescriptor) &source.data_vertex[sourceVertex].data_edge[i];
			TVertexDescriptor targetVertex = getTarget(edSource);
			if (targetVertex == nilVal) continue;
			TEdgeDescriptor edTarget;
			if (!transpose) {
				edTarget = &dest.data_vertex[sourceVertex].data_edge[i];
				assignTarget(edTarget, targetVertex);
			} else {
				edTarget = &dest.data_vertex[targetVertex].data_edge[i];
				assignTarget(edTarget, sourceVertex);
			}
			_assignId(edTarget, _getId(edSource));
			assignCargo(edTarget, getCargo(edSource));
		}
	}
	dest.data_id_managerV = source.data_id_managerV;
	dest.data_id_managerE = source.data_id_managerE;
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
_copyGraph(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& source,
		   Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& dest)
{
	_copyGraph(source,dest,false);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void 
transpose(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& source,
		  Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& dest)
{
	SEQAN_CHECKPOINT
	_copyGraph(source, dest, true);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void 
transpose(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> dest;
	_copyGraph(g, dest, true);
	g = dest;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
numEdges(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g)
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerE);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
numVertices(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline bool 
empty(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g) 
{
	SEQAN_CHECKPOINT
	return (!(idCount(g.data_id_managerV)));
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearEdges(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	clear(g.data_vertex);
	releaseAll(g.data_id_managerE);
	resize(g.data_vertex, getIdUpperBound(g.data_id_managerV));
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void
clearVertices(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	clearEdges(g);
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline void 
clear(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
outDegree(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
		  TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	for(TSize i=0;i<table_length;++i) {
		if ( (TVertexDescriptor) getTarget(&g.data_vertex[vertex].data_edge[i])!=nilVal) ++count;
	}
	return count;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
inDegree(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
		 TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		for(TSize i=0;i<table_length;++i) {
			TEdgeDescriptor ed = (TEdgeDescriptor) &(getValue(it)).data_edge[i];
			if ( (TVertexDescriptor) getTarget(ed)==nilVal) continue;
			if ( (TVertexDescriptor) getTarget(ed)==vertex) ++count;			
		}
	}
	return count;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
degree(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g, 
	   TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return (inDegree(g,vertex)+outDegree(g,vertex));
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
addVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph >::Type TEdge;

	TVertexDescriptor vd = obtainId(g.data_id_managerV);
	if (vd == length(g.data_vertex)) {
		appendValue(g.data_vertex, AutomatonEdgeArray<TEdge, TAlphabet>()); 
	} else {
		value(g.data_vertex, vd) =  AutomatonEdgeArray<TEdge, TAlphabet>();
	}
	return vd;
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	removeOutEdges(g,v); // Remove all outgoing edges
	removeInEdges(g,v); // Remove all incoming edges
	releaseId(g.data_id_managerV, v); // Release id
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TLabel const label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(label);
	TEdgeDescriptor e = &g.data_vertex[source].data_edge[(TSize) letter];
	TId id = obtainId(g.data_id_managerE);
	_assignId(e, id);
	assignTarget(e, target);
	return e;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TLabel, typename TEdgeCargo>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g, 
				  TVertexDescriptor const source, 
				  TVertexDescriptor const target,
				  TLabel const label,
				  TEdgeCargo const cargo)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TAlphabet letter(label);
	TEdgeDescriptor e = addEdge(g,source,target,label);
	assignCargo(e,cargo);
	TId id = obtainId(g.data_id_managerE);
	_assignId(e, id);
	return e;
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline void
removeEdge(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, _getId(edge)) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	releaseId(g.data_id_managerE, _getId(edge));
	assignTarget(edge, _get_nil<TVertexDescriptor>());
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline void
removeEdge(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
	       TVertexDescriptor const source,
	       TVertexDescriptor const target,
		   TLabel const label)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename Size<TGraph>::Type TSize;

	TAlphabet letter(label);
	removeEdge(g, &g.data_vertex[source].data_edge[(TSize) letter]);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void
removeOutEdges(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
			   TVertexDescriptor const vertex)
{
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	for(TSize i=0;i<table_length;++i) {
		TEdgeDescriptor ed = &g.data_vertex[vertex].data_edge[i];
		if ( (TVertexDescriptor) getTarget(ed) == nilVal) continue;
		assignTarget(ed, nilVal);
		releaseId(g.data_id_managerE, _getId(ed));
	}
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void
removeInEdges(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
			  TVertexDescriptor const vertex)
{
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > >::Type TIter;
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		for(TSize i=0;i<table_length;++i) {
			TEdgeDescriptor ed = &(value(it)).data_edge[i];
			if ( (TVertexDescriptor) getTarget(ed) == nilVal) continue;
			if ( (TVertexDescriptor) getTarget(ed) != vertex) continue;
			assignTarget(ed, nilVal);
			releaseId(g.data_id_managerE, _getId(ed));	
		}
	}
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
targetVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return (getTarget(edge));
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
sourceVertex(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
		
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TGraph>::Type TSize;

	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > >::Type TIter;
	for(TIter it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		for(TSize i=0;i<table_length;++i) {
			TEdgeDescriptor ed = &(value(it)).data_edge[i];
			if (getTarget(ed) == nilVal) continue;
			if (ed==edge) return position(it);
		}
	}
	// We should never reach this point
	SEQAN_ASSERT(false)
	return 0;
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Size<TMatrix>::Type TMatrixSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TMatrixSize len = getIdUpperBound(g.data_id_managerV);
	setDimension(mat, 2);
	setLength(mat, 0, len);
	setLength(mat, 1, len);
	resize(mat);
	for (TMatrixSize i=0;i<len*len;++i) value(mat,i) = 0;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		for(TSize i=0;i<table_length;++i) {
			if (((getValue(it)).data_edge[i].data_target!=nilVal))
			{
				TVertexDescriptor const source = position(it);
				TVertexDescriptor const target = (getValue(it)).data_edge[i].data_target;
				assignValue(mat, source*len+target, getValue(mat, source*len+target)+1);
			}
		}
	}
}


template<typename TFile, typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();

	_streamWrite(target,"Automaton - State: (Input / NextState)\n");
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (!idInUse(g.data_id_managerV, position(it))) continue;
		TVertexDescriptor sourceVertex = position(it);
		_streamPutInt(target, sourceVertex);
		_streamWrite(target,": ");
		for(TSize i=0;i<table_length;++i) {
			_streamPut(target, ' ');
			_streamPut(target, '(');
			_streamPut(target, TAlphabet(i));
			_streamPut(target, ' ');
			_streamPut(target, '/');
			_streamPut(target, ' ');
			if (g.data_vertex[sourceVertex].data_edge[i].data_target ==  nilVal) _streamWrite(target,"nil");
			else _streamPutInt(target, g.data_vertex[sourceVertex].data_edge[i].data_target);
			_streamPut(target, ')');
			_streamPut(target, ' ');
		}
		_streamPut(target, '\n');
	}
}


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline void
assignRoot(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g,
		TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	g.root = vertex;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
getRoot(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g)
{
	SEQAN_CHECKPOINT
	return g.root;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type&
root(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>& g)
{
	SEQAN_CHECKPOINT
	return g.root;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor>
inline bool
isRoot(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
	   TVertexDescriptor v)
{
	SEQAN_CHECKPOINT
	return ( (TVertexDescriptor) g.root == v);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
getSuccessor(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
			 TVertexDescriptor vertex,
			 TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	typedef Graph<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> > TGraph;
	typedef typename Size<TAlphabet>::Type TSize;
	TAlphabet letter(c);
	return g.data_vertex[vertex].data_edge[(TSize) letter].data_target;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
getPredecessor(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
					 TVertexDescriptor vertex,
					 TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	TAlphabet letter(c);
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	for(TIterConst it = begin(g.data_vertex);!atEnd(it);goNext(it)) {
		if (( (TVertexDescriptor) (getValue(it)).data_edge[(TSize) letter].data_target!=nilVal) &&
			( (TVertexDescriptor) (getValue(it)).data_edge[(TSize) letter].data_target==vertex))
		{
			return position(it);
	    }
	}
	// We should never reach this point
	SEQAN_ASSERT(false)
	return 0;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TIterator>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
			TVertexDescriptor const vertex,
			TIterator beginIt,
			TIterator endIt)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	TVertexDescriptor succ = vertex;
	while (beginIt!=endIt) {
		TVertexDescriptor tmp = getSuccessor(g,succ,*beginIt);
		if (tmp == nilVal) break;
		succ = tmp;
		++beginIt;
	}
	return succ;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
			TVertexDescriptor const vertex,
			TCharacters const& chars)
{
	SEQAN_CHECKPOINT
	return parseString(g,vertex,begin(chars),end(chars));
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
			TVertexDescriptor const vertex,
			TCharacters const* chars)
{
	SEQAN_CHECKPOINT
	return parseString(g,vertex,chars,chars+length(chars));
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
