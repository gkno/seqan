#ifndef SEQAN_HEADER_GRAPH_EDGESTUMP_H
#define SEQAN_HEADER_GRAPH_EDGESTUMP_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStump
//////////////////////////////////////////////////////////////////////////////
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
// EdgeStump Specific Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename T>
struct Cargo;

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

/////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TSpec> const>::Type&
getCargo(EdgeStump<TCargo, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TSpec> >::Type&
getCargo(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TSpec>
inline typename Cargo<EdgeStump<void, TSpec> const>::Type
getCargo(EdgeStump<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeStump<void, TSpec> >::Type
getCargo(EdgeStump<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TSpec> const>::Type&
cargo(EdgeStump<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStump<TCargo, TSpec> >::Type& 
cargo(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}


template<typename TSpec>
inline typename Cargo<EdgeStump<void, TSpec> >::Type
cargo(EdgeStump<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeStump<void, TSpec> const>::Type
cargo(EdgeStump<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStump<TCargo, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStump<void, TSpec>* es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeStump<TCargo, TSpec>* es, 
			 TVertexDescriptor const t) 
{
	SEQAN_CHECKPOINT
	es->data_target = t;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TSpec> >::Type&
target(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TSpec> >::Type
target(EdgeStump<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TSpec> const>::Type
getTarget(EdgeStump<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStump<TCargo, TSpec> >::Type
getTarget(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec, typename TId>
void 
_assignId(EdgeStump<TCargo, TSpec>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	es->data_id = id;
}

template<typename TCargo, typename TId>
void 
_assignId(EdgeStump<TCargo, WithoutEdgeId>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	// No id -> does nothing
}

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStump<TCargo, TSpec> const>::Type
_getId(EdgeStump<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStump<TCargo, TSpec> >::Type
_getId(EdgeStump<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

template<typename TCargo>
inline typename Id<EdgeStump<TCargo, WithoutEdgeId> >::Type 
_getId(EdgeStump<TCargo, WithoutEdgeId> const* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}

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
