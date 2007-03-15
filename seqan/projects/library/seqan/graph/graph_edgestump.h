#ifndef SEQAN_HEADER_GRAPH_EDGESTUMP_H
#define SEQAN_HEADER_GRAPH_EDGESTUMP_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStump
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Class.EdgeStump:
..cat:Graph
..summary:The EdgeStump class encapsulates a single edge. 
It represents either a list node in the adjacency list of a graph or an array field if edges are stored in an array.
..signature:EdgeStump<TCargo, bool TList, bool TSource, bool TId, TSpec>
..param.TCargo:The cargo type of an edge.
...metafunction:Metafunction.Cargo
...remarks:The cargo can be used to store arbitrary information with an edge.
...default:$void$
..param.TList:Boolean value that indicates whether it is a list node or not.
...remarks:If it is a list node it has one or two next pointers.
...default:$true$
..param.TSource:Boolean value that indicates whether the source is stored in the EdgeStump or not.
...remarks:If this value is true and it is a list node an additional source next pointer is present.
...default:$false$
..param.TId:Boolean value that indicates whether an id is stored in the EdgeStump or not.
Note: Without edge ids external property maps do not work for edges!
...default:$true$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:The default EdgeStump in all graph types does not consider a cargo. 
However, in default usage every graph does store an edge id. 
Edge ids are used to append additional properties to edges with the help of external property maps.
..include:graph.h
*/
template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, true, false, false, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
		TCargo data_cargo;
		EdgeStump* data_nextT;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, true, false, true, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TId data_id;
		TCargo data_cargo;
		EdgeStump* data_nextT;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, true, true, false, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
		TVertexDescriptor data_source;
		TCargo data_cargo;
		EdgeStump* data_nextT;
		EdgeStump* data_nextS;
};

//////////////////////////////////////////////////////////////////////////////


template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, true, true, true, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TVertexDescriptor data_source;
		TId data_id;
		TCargo data_cargo;
		EdgeStump* data_nextT;
		EdgeStump* data_nextS;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, false, false, false, TSpec>
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
		TCargo data_cargo;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, false, false, true, TSpec>
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TId data_id;
		TCargo data_cargo;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, false, true, false, TSpec>
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
		TVertexDescriptor data_source;
		TCargo data_cargo;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
class EdgeStump<TCargo, false, true, true, TSpec>
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TVertexDescriptor data_source;
		TId data_id;
		TCargo data_cargo;
};

//////////////////////////////////////////////////////////////////////////////
//	Graph - Cargoless EdgeStump
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, true, false, false, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
		EdgeStump* data_nextT;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, true, false, true, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TId data_id;
		EdgeStump* data_nextT;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, true, true, false, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
		TVertexDescriptor data_source;
		EdgeStump* data_nextT;
		EdgeStump* data_nextS;
};

//////////////////////////////////////////////////////////////////////////////


template<typename TSpec>
class EdgeStump<void, true, true, true, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TVertexDescriptor data_source;
		TId data_id;
		EdgeStump* data_nextT;
		EdgeStump* data_nextS;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, false, false, false, TSpec>
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, false, false, true, TSpec>
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TId data_id;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, false, true, false, TSpec>
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
		TVertexDescriptor data_source;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class EdgeStump<void, false, true, true, TSpec>
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TVertexDescriptor data_source;
		TId data_id;
};

//////////////////////////////////////////////////////////////////////////////
// EdgeStump - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cargo.param.T.type:Class.EdgeStump

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> > {
	typedef TCargo Type;
};

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> const> {
	typedef TCargo const Type;
};

template<bool TList, bool TSource, bool TId, typename TSpec>
struct Cargo<EdgeStump<void, TList, TSource, TId, TSpec> > {
	typedef void* Type;
};

template<bool TList, bool TSource, bool TId, typename TSpec>
struct Cargo<EdgeStump<void, TList, TSource, TId, TSpec> const> {
	typedef void* Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.EdgeStump

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct Spec<EdgeStump<TCargo, TList, TSource, TId, TSpec> > 
{
	typedef TSpec Type;
};

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
struct Spec<EdgeStump<TCargo, TList, TSource, TId, TSpec> const> 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IdHandler.param.T.type:Class.EdgeStump


template<typename TCargo, bool TList, bool TSource, typename TSpec, typename TIdType>
struct IdHandler<EdgeStump<TCargo, TList, TSource, false, TSpec>, TIdType> {
	// Dummy IdManager
	typedef IdManager<void> Type;
};


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.getCargo:
..cat:Graph
..summary:Get method for the edge cargo.
..signature:getCargo(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Returns the cargo.
..remarks:If cargo is not present the return value is (void*) 0.
..see:Function.cargo
..see:Funktion.assignCargo
*/


template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> const>::Type&
getCargo(EdgeStump<TCargo, TList, TSource, TId, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type&
getCargo(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<void, TList, TSource, TId, TSpec> const>::Type
getCargo(EdgeStump<void, TList, TSource, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<void, TList, TSource, TId, TSpec> >::Type
getCargo(EdgeStump<void, TList, TSource, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.cargo:
..cat:Graph
..summary:Access to the cargo.
..signature:cargo(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Returns a reference to the cargo.
..remarks:If cargo is not present the return value is (void*) 0.
..see:Function.getCargo
..see:Funktion.assignCargo
*/


template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> const>::Type&
cargo(EdgeStump<TCargo, TList, TSource, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type& 
cargo(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}


//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<void, TList, TSource, TId, TSpec> >::Type
cargo(EdgeStump<void, TList, TSource, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec>
inline typename Cargo<EdgeStump<void, TList, TSource, TId, TSpec> const>::Type
cargo(EdgeStump<void, TList, TSource, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignCargo:
..cat:Graph
..summary:Assigns a new cargo to the edge.
..signature:assignCargo(es, cargo)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..param.cargo:New cargo object.
...remarks:Type of the new cargo object must match Cargo<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type.
..returns:void
..remarks:In cargoless EdgeStumps this operation is a NOP.
..see:Function.cargo
..see:Funktion.getCargo
*/

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TList, bool TSource, bool TId, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStump<void, TList, TSource, TId, TSpec>* es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignTarget:
..cat:Graph
..summary:Assigns a target vertex to an edge.
..signature:assignTarget(es, t)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..param.t:Target vertex.
..returns:void
..see:Function.target
..see:Funktion.getTarget
*/


template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es, 
			 TVertexDescriptor const t) 
{
	SEQAN_CHECKPOINT
	es->data_target = t;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.target:
..cat:Graph
..summary:Accesses the target of an EdgeStump.
..signature:target(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Reference to the target vertex.
..see:Function.assignTarget
..see:Funktion.getTarget
*/

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type&
target(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type
target(EdgeStump<TCargo, TList, TSource, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getTarget:
..cat:Graph
..summary:Get method for the target.
..signature:getTarget(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Target vertex.
..see:Function.assignTarget
..see:Funktion.target
*/

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> const>::Type
getTarget(EdgeStump<TCargo, TList, TSource, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, TSource, TId, TSpec> >::Type
getTarget(EdgeStump<TCargo, TList, TSource, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignSource:
..cat:Graph
..summary:Assigns a source vertex to an edge.
..remarks:A source vertex is not required in an edge stump.
However, EdgeStumps can be configured to contain a source vertex, e.g., in undirected graphs.
..signature:assignSource(es, s)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..param.t:Source vertex.
..returns:void
..see:Function.source
..see:Funktion.getSource
*/

template<typename TCargo, bool TList, bool TId, typename TSpec, typename TVertexDescriptor>
inline void 
assignSource(EdgeStump<TCargo, TList, true, TId, TSpec>* es, 
			 TVertexDescriptor const s) 
{
	SEQAN_CHECKPOINT
	es->data_source = s;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec, typename TVertexDescriptor>
inline void 
assignSource(EdgeStump<TCargo, TList, false, TId, TSpec>* es, 
			 TVertexDescriptor const s) 
{
	SEQAN_CHECKPOINT
	// NOP
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.source:
..cat:Graph
..summary:Accesses the source of an EdgeStump.
..remarks:A source vertex is not required in an edge stump.
However, EdgeStumps can be configured to contain a source vertex, e.g., in undirected graphs.
..signature:source(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Reference to the source vertex.
..see:Function.assignSource
..see:Funktion.getSource
*/

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, true, TId, TSpec> >::Type&
source(EdgeStump<TCargo, TList, true, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, true, TId, TSpec> >::Type
source(EdgeStump<TCargo, TList, true, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, false, TId, TSpec> >::Type
source(EdgeStump<TCargo, TList, false, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No source available
	return 0;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, false, TId, TSpec> >::Type
source(EdgeStump<TCargo, TList, false, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No source available
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getSource:
..cat:Graph
..summary:Get method for the source.
..remarks:A source vertex is not required in an edge stump.
However, EdgeStumps can be configured to contain a source vertex, e.g., in undirected graphs.
..signature:getSource(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Source vertex.
..see:Function.assignSource
..see:Funktion.source
*/

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, true, TId, TSpec> const>::Type
getSource(EdgeStump<TCargo, TList, true, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, true, TId, TSpec> >::Type
getSource(EdgeStump<TCargo, TList, true, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, false, TId, TSpec> const>::Type
getSource(EdgeStump<TCargo, TList, false, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// Nop
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TId, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TList, false, TId, TSpec> >::Type
getSource(EdgeStump<TCargo, TList, false, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// Nop
	return 0;
}




//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignNextT:
..cat:Graph
..summary:Assigns another EdgeStump to the next target pointer.
..signature:assignNextT(es, es2)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..param.es2:Pointer to the following EdgeStump.
...type:Class.EdgeStump
..returns:void
..see:Function.nextT
..see:Funktion.getNextT
*/

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline void 
assignNextT(EdgeStump<TCargo, true, TSource, TId, TSpec>* es, 
			EdgeStump<TCargo, true, TSource, TId, TSpec>* es2) 
{
	SEQAN_CHECKPOINT
	es->data_nextT = es2;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.nextT:
..cat:Graph
..summary:Accesses the next target pointer.
..signature:nextT(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Reference to the next target pointer.
..see:Function.assignNextT
..see:Funktion.getNextT
*/

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, TSource, TId, TSpec>* &
nextT(EdgeStump<TCargo, true, TSource, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_nextT;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, TSource, TId, TSpec>* &
nextT(EdgeStump<TCargo, true, TSource, TId, TSpec> const* es) 
{
	return es->data_nextT;
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.getNextT:
..cat:Graph
..summary:Get method for the next target pointer.
..signature:getNextT(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Pointer to the next edge stump in target list.
..see:Function.assignNextT
..see:Funktion.nextT
*/

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, TSource, TId, TSpec>*
getNextT(EdgeStump<TCargo, true, TSource, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_nextT;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TSource, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, TSource, TId, TSpec>*
getNextT(EdgeStump<TCargo, true, TSource, TId, TSpec> const* es) 
{
	return es->data_nextT;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignNextS:
..cat:Graph
..summary:Assigns another EdgeStump to the next source pointer.
..signature:assignNextS(es, es2)
..remarks:EdgeStumps can be configured to have no source. Then there is no next source pointer.
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..param.es2:Pointer to the following EdgeStump.
...type:Class.EdgeStump
..returns:void
..see:Function.nextS
..see:Funktion.getNextS
*/

template<typename TCargo, bool TId, typename TSpec>
inline void 
assignNextS(EdgeStump<TCargo, true, true, TId, TSpec>* es, 
			EdgeStump<TCargo, true, true, TId, TSpec>* es2) 
{
	SEQAN_CHECKPOINT
	es->data_nextS = es2;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline void 
assignNextS(EdgeStump<TCargo, true, false, TId, TSpec>* es, 
			EdgeStump<TCargo, true, false, TId, TSpec>* es2) 
{
	SEQAN_CHECKPOINT
	// Nop
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.nextS:
..cat:Graph
..summary:Accesses the next source pointer.
..signature:nextS(es)
..remarks:EdgeStumps can be configured to have no source. Then there is no next source pointer.
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Reference to the next source pointer.
..see:Function.assignNextS
..see:Funktion.getNextS
*/

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, true, TId, TSpec>* &
nextS(EdgeStump<TCargo, true, true, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_nextS;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, true, TId, TSpec>* &
nextS(EdgeStump<TCargo, true, true, TId, TSpec> const* es) 
{
	return es->data_nextS;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, false, TId, TSpec>*
nextS(EdgeStump<TCargo, true, false, TId, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// Nop
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, false, TId, TSpec>*
nextS(EdgeStump<TCargo, true, false, TId, TSpec> const* es) 
{
	// Nop
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getNextS:
..cat:Graph
..summary:Get method for the next source pointer.
..remarks:EdgeStumps can be configured to have no source. Then there is no next source pointer.
..signature:getNextS(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Pointer to the next edge stump in source list.
..see:Function.assignNextS
..see:Funktion.nextS
*/

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, true, TId, TSpec>*
getNextS(EdgeStump<TCargo, true, true, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_nextS;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TId, typename TSpec>
inline EdgeStump<TCargo, true, false, TId, TSpec>*
getNextS(EdgeStump<TCargo, true, false, TId, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No source pointer
	return 0;
}




//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec, typename TId2>
void 
_assignId(EdgeStump<TCargo, TList, TSource, true, TSpec>* es, 
		  TId2 const id) 
{
	SEQAN_CHECKPOINT
	es->data_id = id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec, typename TId2>
void 
_assignId(EdgeStump<TCargo, TList, TSource, false, TSpec>* es, 
		  TId2 const id) 
{
	SEQAN_CHECKPOINT
	// No id -> does nothing
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
inline typename Id<EdgeStump<TCargo, TList, TSource, true, TSpec> const>::Type
_getId(EdgeStump<TCargo, TList, TSource, true, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
inline typename Id<EdgeStump<TCargo, TList, TSource, true, TSpec> >::Type
_getId(EdgeStump<TCargo, TList, TSource, true, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
inline typename Id<EdgeStump<TCargo, TList, TSource, false, TSpec> >::Type 
_getId(EdgeStump<TCargo, TList, TSource, false, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, bool TList, bool TSource, typename TSpec>
inline typename Id<EdgeStump<TCargo, TList, TSource, false, TSpec> >::Type 
_getId(EdgeStump<TCargo, TList, TSource, false, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
