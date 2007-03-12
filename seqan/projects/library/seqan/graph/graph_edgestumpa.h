#ifndef SEQAN_HEADER_GRAPH_EDGESTUMPA_H
#define SEQAN_HEADER_GRAPH_EDGESTUMPA_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStumpA
//////////////////////////////////////////////////////////////////////////////


/**
.Class.EdgeStumpA:
..cat:Graph
..summary:The EdgeStumpA class encapsulates a single automaton edge. 
It represents an array element in the transition table of an automaton.
..signature:EdgeStumpA<TCargo,TSpec>
..param.TCargo:The cargo type of an edge.
...metafunction:Metafunction.Cargo
...remarks:The cargo can be used to store arbitrary information with an edge.
...default:$void$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:The default EdgeStumpA does not consider a cargo. 
However, it does store an edge id.
Edge ids are used to append additional properties to edges with the help of external property maps.
..include:graph.h
*/
template<typename TCargo, typename TSpec>
class EdgeStumpA
{
	public:
		typedef typename VertexDescriptor<EdgeStumpA>::Type TVertexDescriptor;
		typedef typename Id<EdgeStumpA>::Type TId;
		TVertexDescriptor data_target;
		TId data_id;
		TCargo data_cargo;
};

/**
.Spec.Cargoless EdgeStumpA:
..cat:Graph
..general:Class.EdgeStumpA
..summary:EdgeStumpA that does not reserve space for a cargo.
..signature:EdgeStumpA<void, TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:This is the default EdgeStumpA.
Edge ids can be used to append additional properties to edges with the help of external property maps.
..include:graph.h
*/
template<typename TSpec>
class EdgeStumpA<void, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStumpA>::Type TVertexDescriptor;
		typedef typename Id<EdgeStumpA>::Type TId;
		TVertexDescriptor data_target;
		TId data_id;
};

/**
.Spec.Id-free EdgeStumpA:
..cat:Graph
..general:Class.EdgeStumpA
..summary:EdgeStumpA that does not store an id.
..signature:EdgeStumpA<TCargo, WithoutEdgeId>
..param.TCargo:The cargo type of an edge.
...metafunction:Metafunction.Cargo
...remarks:The cargo can be used to store arbitrary information with an edge.
..remarks:Caution: If edge ids are omitted external property maps do not work.
Necessary edge information can only be stored as an edge cargo.
..include:graph.h
*/
template<typename TCargo>
class EdgeStumpA<TCargo, WithoutEdgeId> 
{
	public:
		typedef typename VertexDescriptor<EdgeStumpA>::Type TVertexDescriptor;
		typedef typename Id<EdgeStumpA>::Type TId;
		TVertexDescriptor data_target;
		TCargo data_cargo;
};

/**
.Spec.Minimal EdgeStumpA:
..cat:Graph
..general:Class.EdgeStumpA
..summary:EdgeStumpA without a cargo and without an id.
..signature:EdgeStumpA<void, WithoutEdgeId>
..remarks:Edges solely connect vertices. No additional edge information can be stored.
..include:graph.h
*/
template<>
class EdgeStumpA<void, WithoutEdgeId> 
{
	public:
		typedef VertexDescriptor<EdgeStumpA>::Type TVertexDescriptor;
		typedef Id<EdgeStumpA>::Type TId;
		TVertexDescriptor data_target;
};

//////////////////////////////////////////////////////////////////////////////
// EdgeStumpA - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cargo.param.T.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec>
struct Cargo<EdgeStumpA<TCargo, TSpec> > {
	typedef TCargo Type;
};

template<typename TCargo, typename TSpec>
struct Cargo<EdgeStumpA<TCargo, TSpec> const> {
	typedef TCargo const Type;
};

template<typename TSpec>
struct Cargo<EdgeStumpA<void, TSpec> > {
	typedef void* Type;
};

template<typename TSpec>
struct Cargo<EdgeStumpA<void, TSpec> const> {
	typedef void* Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec>
struct Spec<EdgeStumpA<TCargo, TSpec> > 
{
	typedef TSpec Type;
};

template<typename TCargo, typename TSpec>
struct Spec<EdgeStumpA<TCargo, TSpec> const> 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.IdHandler.param.T.type:Class.EdgeStumpA

template<typename TCargo, typename TIdType>
struct IdHandler<EdgeStumpA<TCargo, WithoutEdgeId>, TIdType> {
	// Dummy IdManager
	typedef IdManager<void> Type;
};


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Function.getCargo.param.es.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpA<TCargo, TSpec> const>::Type&
getCargo(EdgeStumpA<TCargo, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpA<TCargo, TSpec> >::Type&
getCargo(EdgeStumpA<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpA<void, TSpec> const>::Type
getCargo(EdgeStumpA<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpA<void, TSpec> >::Type
getCargo(EdgeStumpA<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}


//////////////////////////////////////////////////////////////////////////////

///.Function.cargo.param.es.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpA<TCargo, TSpec> >::Type& 
cargo(EdgeStumpA<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpA<TCargo, TSpec> const>::Type& 
cargo(EdgeStumpA<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpA<void, TSpec> >::Type
cargo(EdgeStumpA<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpA<void, TSpec> const>::Type
cargo(EdgeStumpA<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.assignCargo.param.es.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpA<TCargo, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpA<void, TSpec>* es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

//////////////////////////////////////////////////////////////////////////////

///.Function.assignTarget.param.es.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeStumpA<TCargo, TSpec>* es,
			 TVertexDescriptor const t) 
{
	SEQAN_CHECKPOINT
	es->data_target = t;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeStumpA<void, TSpec>* es, 
	     TVertexDescriptor t) 
{
	SEQAN_CHECKPOINT
	es->data_target = t;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.target.param.es.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpA<TCargo, TSpec> >::Type&
target(EdgeStumpA<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpA<TCargo, TSpec> >::Type
target(EdgeStumpA<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}


//////////////////////////////////////////////////////////////////////////////

///.Function.getTarget.param.es.type:Class.EdgeStumpA

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpA<TCargo, TSpec> const>::Type 
getTarget(EdgeStumpA<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpA<TCargo, TSpec> >::Type 
getTarget(EdgeStumpA<TCargo, TSpec>* es) 
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
_assignId(EdgeStumpA<TCargo, TSpec>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	es->data_id = id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TId>
void 
_assignId(EdgeStumpA<TCargo, WithoutEdgeId>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	// No id -> does nothing
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpA<TCargo, TSpec> const>::Type
_getId(EdgeStumpA<TCargo, TSpec> const* es) 
{
	return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpA<TCargo, TSpec> >::Type
_getId(EdgeStumpA<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
inline typename Id<EdgeStumpA<TCargo, WithoutEdgeId> >::Type 
_getId(EdgeStumpA<TCargo, WithoutEdgeId> const* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo>
inline typename Id<EdgeStumpA<TCargo, WithoutEdgeId> >::Type 
_getId(EdgeStumpA<TCargo, WithoutEdgeId>* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
