#ifndef SEQAN_HEADER_GRAPH_EDGESTUMP_H
#define SEQAN_HEADER_GRAPH_EDGESTUMP_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStump
//////////////////////////////////////////////////////////////////////////////

/**
.Class.EdgeStump:
..cat:Graph
..summary:The EdgeStump class encapsulates a single edge. 
It represents a list node in the adjacency list of a directed graph.
..signature:EdgeStump<TCargo,TSpec>
..param.TCargo:The cargo type of an edge.
...metafunction:Metafunction.Cargo
...remarks:The cargo can be used to store arbitrary information with an edge.
...default:$void$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:The default EdgeStump does not consider a cargo. 
However, it does store an edge id.
Edge ids are used to append additional properties to edges with the help of external property maps.
..include:graph.h
*/
template<typename TCargo, typename TSpec>
class EdgeStump
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TId data_id;
		TCargo data_cargo;
		EdgeStump* data_next;
};

/**
.Spec.Cargoless EdgeStump:
..cat:Graph
..general:Class.EdgeStump
..summary:EdgeStump that does not reserve space for a cargo.
..signature:EdgeStump<void, TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:This is the default EdgeStump.
Edge ids can be used to append additional properties to edges with the help of external property maps.
..include:graph.h
*/
template<typename TSpec>
class EdgeStump<void, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TId data_id;
		EdgeStump* data_next;
};

/**
.Spec.Id-free EdgeStump:
..cat:Graph
..general:Class.EdgeStump
..summary:EdgeStump that does not store an id.
..signature:EdgeStump<TCargo, WithoutEdgeId>
..param.TCargo:The cargo type of an edge.
...metafunction:Metafunction.Cargo
...remarks:The cargo can be used to store arbitrary information with an edge.
..remarks:Caution: If edge ids are omitted external property maps do not work.
Necessary edge information can only be stored as an edge cargo.
..include:graph.h
*/
template<typename TCargo>
class EdgeStump<TCargo, WithoutEdgeId> 
{
	public:
		typedef typename VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef typename Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		TCargo data_cargo;
		EdgeStump* data_next;
};

/**
.Spec.Minimal EdgeStump:
..cat:Graph
..general:Class.EdgeStump
..summary:EdgeStump without a cargo and without an id.
..signature:EdgeStump<void, WithoutEdgeId>
..remarks:Edges solely connect vertices. No additional edge information can be stored.
..include:graph.h
*/
template<>
class EdgeStump<void, WithoutEdgeId> 
{
	public:
		typedef VertexDescriptor<EdgeStump>::Type TVertexDescriptor;
		typedef Id<EdgeStump>::Type TId;
		TVertexDescriptor data_target;
		EdgeStump* data_next;
};

//////////////////////////////////////////////////////////////////////////////
// EdgeStump - Metafunctions
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cargo.param.T.type:Class.EdgeStump

template<typename TCargo, typename TSpec>
struct Cargo<EdgeStump<TCargo, TSpec> > {
	typedef TCargo Type;
};

template<typename TCargo, typename TSpec>
struct Cargo<EdgeStump<TCargo, TSpec> const> {
	typedef TCargo const Type;
};

template<typename TSpec>
struct Cargo<EdgeStump<void, TSpec> > {
	typedef void* Type;
};

template<typename TSpec>
struct Cargo<EdgeStump<void, TSpec> const> {
	typedef void* Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.EdgeStump

template<typename TCargo, typename TSpec>
struct Spec<EdgeStump<TCargo, TSpec> > 
{
	typedef TSpec Type;
};

template<typename TCargo, typename TSpec>
struct Spec<EdgeStump<TCargo, TSpec> const> 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IdHandler.param.T.type:Class.EdgeStump

template<typename TCargo, typename TIdType>
struct IdHandler<EdgeStump<TCargo, WithoutEdgeId>, TIdType> {
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

///.Function.getCargo.param.es.type:Class.EdgeStump

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TSpec> const>::Type&
getCargo(EdgeStump<TCargo, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TSpec> >::Type&
getCargo(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStump<void, TSpec> const>::Type
getCargo(EdgeStump<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStump<void, TSpec> >::Type
getCargo(EdgeStump<void, TSpec>* es) 
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

///.Function.cargo.param.es.type:Class.EdgeStump

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TSpec> const>::Type&
cargo(EdgeStump<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TSpec> >::Type& 
cargo(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStump<void, TSpec> >::Type
cargo(EdgeStump<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStump<void, TSpec> const>::Type
cargo(EdgeStump<void, TSpec> const* es) 
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
...remarks:Type of the new cargo object must match Cargo<EdgeStump<TCargo,TSpec> >::Type.
..returns:void
..remarks:In Cargoless EdgeStumps this operation is a NOP.
..see:Function.cargo
..see:Funktion.getCargo
*/

///.Function.assignCargo.param.es.type:Class.EdgeStump

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStump<TCargo, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStump<void, TSpec>* es, 
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

///.Function.assignTarget.param.es.type:Class.EdgeStump

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeStump<TCargo, TSpec>* es, 
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

///.Function.target.param.es.type:Class.EdgeStump

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TSpec> >::Type&
target(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TSpec> >::Type
target(EdgeStump<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getTarget:
..cat:Graph
..summary:Get method for the target.
..signature:target(es)
..param.es:Pointer to the EdgeStump.
...type:Class.EdgeStump
..returns:Target vertex.
..see:Function.assignTarget
..see:Funktion.target
*/

///.Function.getTarget.param.es.type:Class.EdgeStump

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TSpec> const>::Type
getTarget(EdgeStump<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TSpec> >::Type
getTarget(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TId>
void 
_assignId(EdgeStump<TCargo, TSpec>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	es->data_id = id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TId>
void 
_assignId(EdgeStump<TCargo, WithoutEdgeId>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	// No id -> does nothing
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStump<TCargo, TSpec> const>::Type
_getId(EdgeStump<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStump<TCargo, TSpec> >::Type
_getId(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
inline typename Id<EdgeStump<TCargo, WithoutEdgeId> >::Type 
_getId(EdgeStump<TCargo, WithoutEdgeId> const* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
inline typename Id<EdgeStump<TCargo, WithoutEdgeId> >::Type 
_getId(EdgeStump<TCargo, WithoutEdgeId>* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
