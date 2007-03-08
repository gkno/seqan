#ifndef SEQAN_HEADER_GRAPH_EDGESTUMPT_H
#define SEQAN_HEADER_GRAPH_EDGESTUMPT_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStumpT
//////////////////////////////////////////////////////////////////////////////

// Default Edge Stump Tree: No cargo but edge id
template<typename TCargo = void, typename TSpec = Default>
class EdgeStumpT;


/**
.Class.EdgeStumpT:
..cat:Graph
..summary:The EdgeStumpT class encapsulates a single tree edge. 
It represents a list node in the adjacency list of a tree.
..signature:EdgeStumpT<TCargo,TSpec>
..param.TCargo:The cargo type of a tree edge.
...metafunction:Metafunction.Cargo
...remarks:The cargo can be used to store arbitrary information with an edge.
...default:$void$
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:The default EdgeStumpT does not consider a cargo. 
However, every tree edge has an id. 
This id can be used to append additional properties to edges using external property maps.
..include:graph.h
*/
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

/**
.Spec.Cargoless EdgeStumpT:
..cat:Graph
..general:Class.EdgeStumpT
..summary:EdgeStumpT that does not reserve space for a cargo.
..signature:EdgeStumpT<void, TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Default$, see @Tag.Default@.
..remarks:This is the default EdgeStumpT.
Edge ids can be used to append additional properties to edges with the help of external property maps.
..include:graph.h
*/
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
// EdgeStumpT - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Cargo.param.T.type:Class.EdgeStumpT

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

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Spec.param.T.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec>
struct Spec<EdgeStumpT<TCargo, TSpec> > 
{
	typedef TSpec Type;
};

template<typename TCargo, typename TSpec>
struct Spec<EdgeStumpT<TCargo, TSpec> const> 
{
	typedef TSpec Type;
};

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

///.Function.getCargo.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpT<TCargo, TSpec> const>::Type&
getCargo(EdgeStumpT<TCargo, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpT<TCargo, TSpec> >::Type&
getCargo(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpT<void, TSpec> const>::Type
getCargo(EdgeStumpT<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpT<void, TSpec> >::Type
getCargo(EdgeStumpT<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.cargo.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpT<TCargo, TSpec> const>::Type&
cargo(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeStumpT<TCargo, TSpec> >::Type& 
cargo(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TSpec>
inline typename Cargo<EdgeStumpT<void, TSpec> >::Type
cargo(EdgeStumpT<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Cargo<EdgeStumpT<void, TSpec> const>::Type
cargo(EdgeStumpT<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.assignCargo.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpT<TCargo, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeStumpT<void, TSpec>* es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignParent:
..cat:Graph
..summary:Assigns a parent vertex to an edge.
..signature:assignParent(es, p)
..param.es:Pointer to the EdgeStumpT.
...type:Class.EdgeStumpT
..param.p:Parent vertex.
..returns:void
..see:Function.parent
..see:Funktion.getParent
*/

///.Function.assignParent.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignParent(EdgeStumpT<TCargo, TSpec>* es, 
			 TVertexDescriptor const p) 
{
	SEQAN_CHECKPOINT
	es->data_parent = p;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.parent:
..cat:Graph
..summary:Accesses the parent vertex of an edge.
..signature:parent(es)
..param.es:Pointer to the EdgeStumpT.
...type:Class.EdgeStumpT
..returns:Reference to the parent vertex.
..see:Function.assignParent
..see:Funktion.getParent
*/

///.Function.parent.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type&
parent(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_parent;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type
parent(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_parent;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getParent:
..cat:Graph
..summary:Get method for the parent vertex of an edge.
..signature:getParent(es)
..param.es:Pointer to the EdgeStumpT.
...type:Class.EdgeStumpT
..returns:Parent vertex.
..see:Function.parent
..see:Funktion.assignParent
*/

///.Function.parent.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> const>::Type
getParent(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_parent;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type
getParent(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_parent;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignChild:
..cat:Graph
..summary:Assigns a child vertex to an edge.
..signature:assignChild(es, ch)
..param.es:Pointer to the EdgeStumpT.
...type:Class.EdgeStumpT
..param.ch:Child vertex.
..returns:void
..see:Function.child
..see:Funktion.getChild
*/

///.Function.assignChild.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignChild(EdgeStumpT<TCargo, TSpec>* es,
			TVertexDescriptor const ch) 
{
	SEQAN_CHECKPOINT
	es->data_child = ch;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.child:
..cat:Graph
..summary:Accesses the child vertex of an edge.
..signature:child(es)
..param.es:Pointer to the EdgeStumpT.
...type:Class.EdgeStumpT
..returns:Reference to the child vertex.
..see:Function.assignChild
..see:Funktion.getChild
*/

///.Function.child.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type&
child(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_child;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type
child(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_child;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getChild:
..cat:Graph
..summary:Get method for the child vertex of an edge.
..signature:getChild(es)
..param.es:Pointer to the EdgeStumpT.
...type:Class.EdgeStumpT
..returns:Child vertex.
..see:Function.assignChild
..see:Funktion.getChild
*/

///.Function.getChild.param.es.type:Class.EdgeStumpT

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> const>::Type
getChild(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_child;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeStumpT<TCargo, TSpec> >::Type
getChild(EdgeStumpT<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_child;
}

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStumpT<TCargo, TSpec> const>::Type
_getId(EdgeStumpT<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// Edge id = child id for tree edges
	return es->data_child;
}

//////////////////////////////////////////////////////////////////////////////

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
