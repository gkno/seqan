#ifndef SEQAN_HEADER_GRAPH_EDGESTUMPU_H
#define SEQAN_HEADER_GRAPH_EDGESTUMPU_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStumpU
//////////////////////////////////////////////////////////////////////////////

/**
.Class.EdgeStumpU:
..cat:Graph
..summary:The EdgeStumpU class encapsulates a single undirected edge. 
It represents a list node in the adjacency list of a undirected graph.
..signature:EdgeStumpU<TCargo,TSpec>
..param.TCargo:The cargo type of an edge.
...metafunction:Metafunction.Cargo
...remarks:The cargo can be used to store arbitrary information with an edge.
...default:$void$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:The default EdgeStumpU does not consider a cargo. 
However, it does store an edge id.
Edge ids are used to append additional properties to edges with the help of external property maps.
..include:graph.h
*/
template<typename TCargo, typename TSpec>
class EdgeStumpU
{
	public:
		typedef typename VertexDescriptor<EdgeStumpU>::Type TVertexDescriptor;
		typedef typename Id<EdgeStumpU>::Type TId;
		TVertexDescriptor data_source;		// Smaller vertex id is always source
		TVertexDescriptor data_target;		// Larger vertex id is always target
		TId data_id;
		TCargo data_cargo;
		EdgeStumpU* data_next_source;
		EdgeStumpU* data_next_target;
};

/**
.Spec.Cargoless EdgeStumpU:
..cat:Graph
..general:Class.EdgeStumpU
..summary:EdgeStumpU that does not reserve space for a cargo.
..signature:EdgeStumpU<void, TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:This is the default EdgeStumpU.
Edge ids can be used to append additional properties to edges with the help of external property maps.
..include:graph.h
*/
template<typename TSpec>
class EdgeStumpU<void, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStumpU>::Type TVertexDescriptor;
		typedef typename Id<EdgeStumpU>::Type TId;
		TVertexDescriptor data_source;
		TVertexDescriptor data_target;
		TId data_id;
		EdgeStumpU* data_next_source;
		EdgeStumpU* data_next_target;
};

/**
.Spec.Id-free EdgeStumpU:
..cat:Graph
..general:Class.EdgeStumpU
..summary:EdgeStumpU that does not store an id.
..signature:EdgeStumpU<TCargo, WithoutEdgeId>
..param.TCargo:The cargo type of an edge.
...metafunction:Metafunction.Cargo
...remarks:The cargo can be used to store arbitrary information with an edge.
..remarks:Caution: If edge ids are omitted external property maps do not work.
Necessary edge information can only be stored as an edge cargo.
..include:graph.h
*/
template<typename TCargo>
class EdgeStumpU<TCargo, WithoutEdgeId> 
{
	public:
		typedef typename VertexDescriptor<EdgeStumpU>::Type TVertexDescriptor;
		typedef typename Id<EdgeStumpU>::Type TId;
		TVertexDescriptor data_source;
		TVertexDescriptor data_target;
		TCargo data_cargo;
		EdgeStumpU* data_next_source;
		EdgeStumpU* data_next_target;
};

/**
.Spec.Minimal EdgeStumpU:
..cat:Graph
..general:Class.EdgeStumpU
..summary:EdgeStumpU without a cargo and without an id.
..signature:EdgeStumpU<void, WithoutEdgeId>
..remarks:Edges solely connect vertices. No additional edge information can be stored.
..include:graph.h
*/
template<>
class EdgeStumpU<void, WithoutEdgeId> 
{
	public:
		typedef VertexDescriptor<EdgeStumpU>::Type TVertexDescriptor;
		typedef Id<EdgeStumpU>::Type TId;
		TVertexDescriptor data_source;
		TVertexDescriptor data_target;
		EdgeStumpU* data_next_source;
		EdgeStumpU* data_next_target;
};

//////////////////////////////////////////////////////////////////////////////
// EdgeStumpU - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cargo.param.T.type:Class.EdgeStumpU
template<typename TCargo, typename TSpec>
struct Cargo<EdgeStumpU<TCargo, TSpec> > {
	typedef TCargo Type;
};

template<typename TCargo, typename TSpec>
struct Cargo<EdgeStumpU<TCargo, TSpec> const> {
	typedef TCargo const Type;
};

template<typename TSpec>
struct Cargo<EdgeStumpU<void, TSpec> > {
	typedef void* Type;
};

template<typename TSpec>
struct Cargo<EdgeStumpU<void, TSpec> const> {
	typedef void* Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec>
struct Spec<EdgeStumpU<TCargo, TSpec> > 
{
	typedef TSpec Type;
};

template<typename TCargo, typename TSpec>
struct Spec<EdgeStumpU<TCargo, TSpec> const> 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IdHandler.param.T.type:Class.EdgeStumpU

template<typename TCargo, typename TIdType>
struct IdHandler<EdgeStumpU<TCargo, WithoutEdgeId>, TIdType> {
	// Dummy IdManager
	typedef IdManager<void> Type;
};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

///.Function.getCargo.param.es.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpU<TCargo, TSpec> const>::Type&
getCargo(EdgeStumpU<TCargo, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpU<TCargo, TSpec> >::Type&
getCargo(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpU<void, TSpec> const>::Type
getCargo(EdgeStumpU<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpU<void, TSpec> >::Type
getCargo(EdgeStumpU<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.cargo.param.es.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpU<TCargo, TSpec> const>::Type&
cargo(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpU<TCargo, TSpec> >::Type& 
cargo(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpU<void, TSpec> >::Type
cargo(EdgeStumpU<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpU<void, TSpec> const>::Type
cargo(EdgeStumpU<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.assignCargo.param.es.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpU<TCargo, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpU<void, TSpec>* es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

//////////////////////////////////////////////////////////////////////////////

///.Function.assignTarget.param.es.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeStumpU<TCargo, TSpec>* es, 
			 TVertexDescriptor const t) 
{
	SEQAN_CHECKPOINT
	es->data_target = t;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.target.param.es.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type&
target(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type
target(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.getTarget.param.es.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> const>::Type
getTarget(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type
getTarget(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignSource:
..cat:Graph
..summary:Assigns a source vertex to an edge.
..signature:assignSource(es, s)
..param.es:Pointer to the EdgeStumpU.
...type:Class.EdgeStumpU
..param.s:Source vertex.
..returns:void
..see:Function.source
..see:Funktion.getSource
*/

///.Function.assignSource.param.es.type:Class.EdgeStump

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignSource(EdgeStumpU<TCargo, TSpec>* es, 
			 TVertexDescriptor const s) 
{
	SEQAN_CHECKPOINT
	es->data_source = s;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.source:
..cat:Graph
..summary:Accesses the source of an EdgeStumpU.
..signature:source(es)
..param.es:Pointer to the EdgeStumpU.
...type:Class.EdgeStumpU
..returns:Reference to the source vertex.
..see:Function.assignSource
..see:Funktion.getSource
*/

///.Function.source.param.es.type:Class.EdgeStumpU


template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type&
source(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type
source(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getSource:
..cat:Graph
..summary:Get method for the source vertex.
..signature:getSource(es)
..param.es:Pointer to the EdgeStumpU.
...type:Class.EdgeStumpU
..returns:Source vertex.
..see:Function.assignSource
..see:Funktion.source
*/

///.Function.getSource.param.es.type:Class.EdgeStumpU

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> const>::Type
getSource(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type
getSource(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TId>
void 
_assignId(EdgeStumpU<TCargo, TSpec>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	es->data_id = id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TId>
void 
_assignId(EdgeStumpU<TCargo, WithoutEdgeId>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	// No id -> does nothing
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpU<TCargo, TSpec> const>::Type
_getId(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpU<TCargo, TSpec> >::Type
_getId(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
inline typename Id<EdgeStumpU<TCargo, WithoutEdgeId> >::Type 
_getId(EdgeStumpU<TCargo, WithoutEdgeId> const* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
inline typename Id<EdgeStumpU<TCargo, WithoutEdgeId> >::Type 
_getId(EdgeStumpU<TCargo, WithoutEdgeId>* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
