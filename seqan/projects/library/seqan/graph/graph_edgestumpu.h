#ifndef SEQAN_HEADER_GRAPH_EDGESTUMPU_H
#define SEQAN_HEADER_GRAPH_EDGESTUMPU_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStumpU
//////////////////////////////////////////////////////////////////////////////
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
// EdgeStumpU Specific Metafunctions
//////////////////////////////////////////////////////////////////////////////
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

/////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpU<TCargo, TSpec> const>::Type&
getCargo(EdgeStumpU<TCargo, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpU<TCargo, TSpec> >::Type&
getCargo(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TSpec>
inline typename Cargo<EdgeStumpU<void, TSpec> const>::Type
getCargo(EdgeStumpU<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeStumpU<void, TSpec> >::Type
getCargo(EdgeStumpU<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpU<TCargo, TSpec> const>::Type&
cargo(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpU<TCargo, TSpec> >::Type& 
cargo(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}


template<typename TSpec>
inline typename Cargo<EdgeStumpU<void, TSpec> >::Type
cargo(EdgeStumpU<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeStumpU<void, TSpec> const>::Type
cargo(EdgeStumpU<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpU<TCargo, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpU<void, TSpec>* es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeStumpU<TCargo, TSpec>* es, 
			 TVertexDescriptor const t) 
{
	SEQAN_CHECKPOINT
	es->data_target = t;
}

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignSource(EdgeStumpU<TCargo, TSpec>* es, 
			 TVertexDescriptor const s) 
{
	SEQAN_CHECKPOINT
	es->data_source = s;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type&
target(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type
target(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type&
source(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type
source(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> const>::Type
getTarget(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type
getTarget(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> const>::Type
getSource(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpU<TCargo, TSpec> >::Type
getSource(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_source;
}

template<typename TCargo, typename TSpec, typename TId>
void 
_assignId(EdgeStumpU<TCargo, TSpec>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	es->data_id = id;
}

template<typename TCargo, typename TId>
void 
_assignId(EdgeStumpU<TCargo, WithoutEdgeId>* es, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	// No id -> does nothing
}

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpU<TCargo, TSpec> const>::Type
_getId(EdgeStumpU<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpU<TCargo, TSpec> >::Type
_getId(EdgeStumpU<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_id;
}

template<typename TCargo>
inline typename Id<EdgeStumpU<TCargo, WithoutEdgeId> >::Type 
_getId(EdgeStumpU<TCargo, WithoutEdgeId> const* es) 
{
	SEQAN_CHECKPOINT
	// No real id
	return 0;
}

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
