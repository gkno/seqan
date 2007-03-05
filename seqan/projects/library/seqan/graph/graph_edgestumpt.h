#ifndef SEQAN_HEADER_GRAPH_EDGESTUMPT_H
#define SEQAN_HEADER_GRAPH_EDGESTUMPT_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStumpT
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TSpec>
class EdgeStumpT
{
	public:
		typedef typename VertexDescriptor<EdgeStumpT>::Type TVertexDescriptor;
		typedef typename Id<EdgeStumpT>::Type TId;
		TVertexDescriptor data_parent;		// Parent vertex descriptor
		TVertexDescriptor data_child;		// Child vertex descriptor
		TCargo data_cargo;
		EdgeStumpT* data_next;
};


template<typename TSpec>
class EdgeStumpT<void, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeStumpT>::Type TVertexDescriptor;
		typedef typename Id<EdgeStumpT>::Type TId;
		TVertexDescriptor data_parent;		// Parent vertex descriptor
		TVertexDescriptor data_child;		// Child vertex descriptor
		EdgeStumpT* data_next;
};


//////////////////////////////////////////////////////////////////////////////
// EdgeStumpt Specific Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TSpec>
struct Cargo<EdgeStumpT<TCargo, TSpec> > {
	typedef TCargo Type;
};

template<typename TCargo, typename TSpec>
struct Cargo<EdgeStumpT<TCargo, TSpec> const> {
	typedef TCargo const Type;
};

template<typename TSpec>
struct Cargo<EdgeStumpT<void, TSpec> > {
	typedef void* Type;
};

template<typename TSpec>
struct Cargo<EdgeStumpT<void, TSpec> const> {
	typedef void* Type;
};

/////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpT<TCargo, TSpec> const>::Type&
getCargo(EdgeStumpT<TCargo, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpT<TCargo, TSpec> >::Type&
getCargo(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TSpec>
inline typename Cargo<EdgeStumpT<void, TSpec> const>::Type
getCargo(EdgeStumpT<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeStumpT<void, TSpec> >::Type
getCargo(EdgeStumpT<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpT<TCargo, TSpec> const>::Type&
cargo(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpT<TCargo, TSpec> >::Type& 
cargo(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}


template<typename TSpec>
inline typename Cargo<EdgeStumpT<void, TSpec> >::Type
cargo(EdgeStumpT<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeStumpT<void, TSpec> const>::Type
cargo(EdgeStumpT<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpT<TCargo, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpT<void, TSpec>* es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignParent(EdgeStumpT<TCargo, TSpec>* es, 
			 TVertexDescriptor const p) 
{
	SEQAN_CHECKPOINT
	es->data_parent = p;
}

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignChild(EdgeStumpT<TCargo, TSpec>* es,
			TVertexDescriptor const ch) 
{
	SEQAN_CHECKPOINT
	es->data_child = ch;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type&
parent(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_parent;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type
parent(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_parent;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type&
chid(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_child;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type
child(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_child;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> const>::Type
getParent(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_parent;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type
getParent(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_parent;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> const>::Type
getChild(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_child;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type
getChild(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_child;
}

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpT<TCargo, TSpec> const>::Type
_getId(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// Edge id = child id for tree edges
	return es->data_child;
}

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpT<TCargo, TSpec> >::Type
_getId(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// Edge id = child id for tree edges
	return es->data_child;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
