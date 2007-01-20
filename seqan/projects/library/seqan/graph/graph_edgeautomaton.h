#ifndef SEQAN_HEADER_GRAPH_EDGEAUTOMATON_H
#define SEQAN_HEADER_GRAPH_EDGEAUTOMATON_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeStump
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TSpec>
class EdgeAutomaton
{
	public:
		typedef typename VertexDescriptor<EdgeAutomaton>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
		TCargo data_cargo;
};


template<typename TSpec>
class EdgeAutomaton<void, TSpec> 
{
	public:
		typedef typename VertexDescriptor<EdgeAutomaton>::Type TVertexDescriptor;
		TVertexDescriptor data_target;
};

//////////////////////////////////////////////////////////////////////////////
// EdgeStump Specific Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename T>
struct Cargo;

template<typename TCargo, typename TSpec>
struct Cargo<EdgeAutomaton<TCargo, TSpec> > {
	typedef TCargo Type;
};

template<typename TCargo, typename TSpec>
struct Cargo<EdgeAutomaton<TCargo, TSpec> const> {
	typedef TCargo const Type;
};

template<typename TSpec>
struct Cargo<EdgeAutomaton<void, TSpec> > {
	typedef void* Type;
};

template<typename TSpec>
struct Cargo<EdgeAutomaton<void, TSpec> const> {
	typedef void* Type;
};

/////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeAutomaton<TCargo, TSpec> const>::Type&
getCargo(EdgeAutomaton<TCargo, TSpec> const* es)
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeAutomaton<TCargo, TSpec> >::Type&
getCargo(EdgeAutomaton<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TSpec>
inline typename Cargo<EdgeAutomaton<void, TSpec> const>::Type
getCargo(EdgeAutomaton<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeAutomaton<void, TSpec> >::Type
getCargo(EdgeAutomaton<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeAutomaton<TCargo, TSpec> >::Type& 
cargo(EdgeAutomaton<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeAutomaton<TCargo, TSpec> const>::Type& 
cargo(EdgeAutomaton<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_cargo;
}


template<typename TSpec>
inline typename Cargo<EdgeAutomaton<void, TSpec> >::Type
cargo(EdgeAutomaton<void, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeAutomaton<void, TSpec> const>::Type
cargo(EdgeAutomaton<void, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeAutomaton<TCargo, TSpec>* es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es->data_cargo = t;
}

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeAutomaton<void, TSpec>* es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeAutomaton<TCargo, TSpec>* es, 
			 TVertexDescriptor const t) 
{
	SEQAN_CHECKPOINT
	es->data_target = t;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeAutomaton<TCargo, TSpec> >::Type&
target(EdgeAutomaton<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeAutomaton<TCargo, TSpec> >::Type
target(EdgeAutomaton<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeAutomaton<TCargo, TSpec> const>::Type 
getTarget(EdgeAutomaton<TCargo, TSpec> const* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeAutomaton<TCargo, TSpec> >::Type 
getTarget(EdgeAutomaton<TCargo, TSpec>* es) 
{
	SEQAN_CHECKPOINT
	return es->data_target;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
