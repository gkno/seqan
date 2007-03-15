#ifndef SEQAN_HEADER_GRAPH_IMPL_AUTOMATON_H
#define SEQAN_HEADER_GRAPH_IMPL_AUTOMATON_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Automaton
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template<typename TEdge, typename TAlphabet>
class AutomatonEdgeArray {
public:
	TEdge data_edge[ValueSize<TAlphabet>::VALUE];

	AutomatonEdgeArray() 
	{
		typedef typename VertexDescriptor<TEdge>::Type TVertexDescriptor;
		typedef typename Size<TAlphabet>::Type TSize;
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
		for(TSize i=0;i<table_length;++i) {
			assignTarget(&data_edge[i], nilVal);
		}
	}
};

/**
.Spec.Automaton:
..cat:Graph
..general:Class.Graph
..summary:An automaton that stores the edges in a transition table.
..signature:Graph<Automaton<TAlphabet, TCargo, TSpec> > 
..param.TAlphabet:The alphabet type that is used for the transition labels.
...metafunction:Metafunction.Alphabet
...remarks:Use @Metafunction.Alphabet@ to get the type of the labels in an automaton.
...default:$char$
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
template<typename TAlphabet, typename TCargo, typename TSpec>
class Graph<Automaton<TAlphabet, TCargo, TSpec> > 
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
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& source,
		   Graph<Automaton<TAlphabet, TCargo, TSpec> >& dest,
		   bool transpose)
{
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Cargo<TEdge>::Type TCargoType;
	typedef typename Size<TAlphabet>::Type TSize;
	typedef typename Iterator<String<AutomatonEdgeArray<TEdge, TAlphabet> > const>::Type TIterConst;
	
	clear(dest);
	resize(dest.data_vertex, length(source.data_vertex));
	dest.root = source.root;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
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


template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
_copyGraph(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& source,
		   Graph<Automaton<TAlphabet, TCargo, TSpec> >& dest)
{
	_copyGraph(source,dest,false);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void 
transpose(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& source,
		  Graph<Automaton<TAlphabet, TCargo, TSpec> >& dest)
{
	SEQAN_CHECKPOINT
	_copyGraph(source, dest, true);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void 
transpose(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	Graph<Automaton<TAlphabet, TCargo, TSpec> > dest;
	_copyGraph(g, dest, true);
	g = dest;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
numEdges(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerE);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
numVertices(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return idCount(g.data_id_managerV);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline bool 
empty(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g) 
{
	SEQAN_CHECKPOINT
	return (!(idCount(g.data_id_managerV)));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
clearEdges(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	clear(g.data_vertex);
	releaseAll(g.data_id_managerE);
	resize(g.data_vertex, getIdUpperBound(g.data_id_managerV));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void
clearVertices(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	clearEdges(g);
	releaseAll(g.data_id_managerV);
	clear(g.data_vertex);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline void 
clear(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	clearVertices(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
outDegree(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
		  TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	for(TSize i=0;i<table_length;++i) {
		if ( (TVertexDescriptor) getTarget(&g.data_vertex[vertex].data_edge[i])!=nilVal) ++count;
	}
	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
inDegree(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
		 TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TSize count=0;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
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

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline typename Size<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
degree(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g, 
	   TVertexDescriptor const vertex) 
{
	SEQAN_CHECKPOINT
	return (inDegree(g,vertex)+outDegree(g,vertex));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
addVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
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

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
removeVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g, 
			 TVertexDescriptor const v) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)

	removeOutEdges(g,v); // Remove all outgoing edges
	removeInEdges(g,v); // Remove all incoming edges
	releaseId(g.data_id_managerV, v); // Release id
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g, 
		TVertexDescriptor const source, 
		TVertexDescriptor const target,
		TLabel const label) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TEdgeDescriptor e = findEdge(g, source, label);
	_assignId(e, obtainId(g.data_id_managerE));
	assignTarget(e, target);
	return e;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel, typename TEdgeCargo>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
addEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g, 
				  TVertexDescriptor const source, 
				  TVertexDescriptor const target,
				  TLabel const label,
				  TEdgeCargo const cargo)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
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

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline void
removeEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
	       TVertexDescriptor const source,
	       TVertexDescriptor const target,
		   TLabel const label)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, source) == true)
	SEQAN_ASSERT(idInUse(g.data_id_managerV, target) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;

	TAlphabet letter(label);
	removeEdge(g, &g.data_vertex[source].data_edge[(TSize) letter]);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline void
removeEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		   TEdgeDescriptor const edge)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, _getId(edge)) == true)
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	releaseId(g.data_id_managerE, _getId(edge));
	assignTarget(edge, getNil<TVertexDescriptor>());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeOutEdges(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
			   TVertexDescriptor const vertex)
{
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	for(TSize i=0;i<table_length;++i) {
		TEdgeDescriptor ed = &g.data_vertex[vertex].data_edge[i];
		if ( (TVertexDescriptor) getTarget(ed) == nilVal) continue;
		assignTarget(ed, nilVal);
		releaseId(g.data_id_managerE, _getId(ed));
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
removeInEdges(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
			  TVertexDescriptor const vertex)
{
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
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

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
targetVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
	return (getTarget(edge));
}

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeDescriptor>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
sourceVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
			 TEdgeDescriptor const edge) 
{
	SEQAN_CHECKPOINT
		
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TGraph>::Type TSize;

	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
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

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TCargo, typename TSpec, typename TMatrix>
inline void
getAdjacencyMatrix(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
				   TMatrix& mat) 
{
	SEQAN_CHECKPOINT

	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Size<TMatrix>::Type TMatrixSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
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

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TAlphabet, typename TCargo, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
	  TIDString const &,
	  Raw)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeType<TGraph>::Type TEdge;
	typedef typename Size<TAlphabet>::Type TSize;
	TSize table_length = ValueSize<TAlphabet>::VALUE;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

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

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignRoot:
..cat:Spec.Automaton
..summary:Assigns a new root vertex to the automaton.
..signature:assignRoot(g, vertex)
..param.g:A automaton.
...type:Spec.Automaton
..param.vertex:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:void
..see:Function.getRoot
..see:Function.root
..see:Function.isRoot
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void
assignRoot(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		TVertexDescriptor const vertex)
{
	SEQAN_CHECKPOINT
	g.root = vertex;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.root:
..cat:Spec.Automaton
..cat:Spec.Tree
..summary:Gets a reference to the root of the tree.
..signature:root(g)
..param.g:A tree or an automaton.
...type:Spec.Automaton
...type:Spec.Tree
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.getRoot
..see:Function.assignRoot
..see:Function.isRoot
*/

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type&
root(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	return g.root;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getRoot:
..cat:Spec.Automaton
..cat:Spec.Tree
..summary:Get method for the root of a tree or an automaton.
..signature:getRoot(g)
..param.g:An automaton or a tree.
...type:Spec.Automaton
...type:Spec.Tree
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.assignRoot
..see:Function.root
..see:Function.isRoot
*/

template<typename TAlphabet, typename TCargo, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
getRoot(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g)
{
	SEQAN_CHECKPOINT
	return g.root;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.isRoot:
..cat:Spec.Automaton
..cat:Spec.Tree
..summary:Tests whether a given vertex is the root or not.
..signature:isRoot(g, v)
..param.g:An automaton or a tree.
...type:Spec.Automaton
...type:Spec.Tree
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:True if vertex is the root.
..see:Function.getRoot
..see:Function.assignRoot
..see:Function.root
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor>
inline bool
isRoot(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
	   TVertexDescriptor v)
{
	SEQAN_CHECKPOINT
	return ( (TVertexDescriptor) g.root == v);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.findEdge:
..cat:Spec.Automaton
..cat:Spec.Tree
..summary:Finds an edge. 
In an automaton an edge is uniquely defined by a vertex and a label.
In a tree two adjacent vertices uniquely define an edge.
..signature:findEdge(g, v, c)
..signature:findEdge(g, v, w)
..param.g:An automaton or a tree.
...type:Spec.Automaton
...type:Spec.Tree
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.c:An edge label.
...type:Metafunction.Alphabet
..param.w:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..returns:An edge descriptor
...type:Metafunction.EdgeDescriptor
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
findEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		 TVertexDescriptor const v,
		 TLabel const c)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;

	TAlphabet label(c);
	return &g.data_vertex[v].data_edge[(TSize) label];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TLabel>
inline typename EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
findEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
		 TVertexDescriptor const v,
		 TLabel const c)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, v) == true)
	
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	TGraph* graph = const_cast<TGraph*>(&g);
	return findEdge(*graph, v, c);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.getSuccessor:
..cat:Spec.Automaton
..cat:Spec.Word Graph
..summary:Gets the successor for a given vertex and an edge label.
For an automaton a single character is required whereas for a word graph getSuccessor takes a string.
..signature:getSuccessor(g, v, c)
..signature:getSuccessor(g, v, str)
..param.g:An automaton or a word graph.
...type:Spec.Automaton
...type:Spec.Word Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.c:A label.
...type:Metafunction.Alphabet
..param.str:A string.
...type:Class.String
..returns:A vertex descriptor or nil if successor is not defined.
...type:Metafunction.VertexDescriptor
..see:Function.parseString
..see:Function.getNil
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TChar>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
getSuccessor(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			 TVertexDescriptor vertex,
			 TChar const c) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	return getTarget(findEdge(g, vertex, c));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.parseString:
..cat:Spec.Automaton
..cat:Spec.Word Graph
..summary:Parses a string one character at a time and moves accordingly in the automaton.
..signature:parseString(g, v, beginIt, endIt)
..signature:parseString(g, v, str)
..param.g:An automaton or a word graph.
...type:Spec.Automaton
...type:Spec.Word Graph
..param.v:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..param.beginIt:An iterator to the beginning of the string.
...type:Metafunction.Iterator.
..param.endIt:An iterator to the end of the string.
...type:Metafunction.Iterator.
..param.str:A string.
...type:Class.String
..returns:A vertex descriptor.
...type:Metafunction.VertexDescriptor
..see:Function.getSuccessor
*/

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TIterator>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			TVertexDescriptor const vertex,
			TIterator beginIt,
			TIterator endIt)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(g.data_id_managerV, vertex) == true)
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor succ = vertex;
	while (beginIt!=endIt) {
		TVertexDescriptor tmp = getSuccessor(g,succ,*beginIt);
		if (tmp == nilVal) break;
		succ = tmp;
		++beginIt;
	}
	return succ;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo,  typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			TVertexDescriptor const vertex,
			TCharacters const& chars)
{
	SEQAN_CHECKPOINT
	return parseString(g,vertex,begin(chars),end(chars));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TCharacters>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type 
parseString(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
			TVertexDescriptor const vertex,
			TCharacters const* chars)
{
	SEQAN_CHECKPOINT
	return parseString(g,vertex,chars,chars+length(chars));
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
