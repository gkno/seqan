#ifndef SEQAN_HEADER_GRAPH_EDGEAUTOMATON_H
#define SEQAN_HEADER_GRAPH_EDGEAUTOMATON_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
//	Graph - EdgeAutomaton
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
// EdgeAutomaton Specific Metafunctions
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
getCargo(EdgeAutomaton<TCargo, TSpec> const& es)
{
	SEQAN_CHECKPOINT
	return es.data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeAutomaton<TCargo, TSpec> >::Type&
getCargo(EdgeAutomaton<TCargo, TSpec>& es) 
{
	SEQAN_CHECKPOINT
	return es.data_cargo;
}

template<typename TSpec>
inline typename Cargo<EdgeAutomaton<void, TSpec> const>::Type
getCargo(EdgeAutomaton<void, TSpec> const& es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeAutomaton<void, TSpec> >::Type
getCargo(EdgeAutomaton<void, TSpec>& es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeAutomaton<TCargo, TSpec> >::Type& 
cargo(EdgeAutomaton<TCargo, TSpec>& es) 
{
	SEQAN_CHECKPOINT
	return es.data_cargo;
}

template<typename TCargo, typename TSpec>
inline typename Cargo<EdgeAutomaton<TCargo, TSpec> const>::Type& 
cargo(EdgeAutomaton<TCargo, TSpec> const& es) 
{
	SEQAN_CHECKPOINT
	return es.data_cargo;
}


template<typename TSpec>
inline typename Cargo<EdgeAutomaton<void, TSpec> >::Type
cargo(EdgeAutomaton<void, TSpec>& es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TSpec>
inline typename Cargo<EdgeAutomaton<void, TSpec> const>::Type
cargo(EdgeAutomaton<void, TSpec> const& es) 
{
	SEQAN_CHECKPOINT
	// No real cargo
	return 0;
}

template<typename TCargo, typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeAutomaton<TCargo, TSpec>& es,
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	es.data_cargo = t;
}

template<typename TSpec, typename TCargo2>
inline void 
assignCargo(EdgeAutomaton<void, TSpec>& es, 
			TCargo2 const& t) 
{
	SEQAN_CHECKPOINT
	// No real cargo
}

template<typename TCargo, typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeAutomaton<TCargo, TSpec>& es, 
	     TVertexDescriptor const t) 
{
	SEQAN_CHECKPOINT
	es.data_target = t;
}

template<typename TSpec, typename TVertexDescriptor>
inline void 
assignTarget(EdgeAutomaton<void, TSpec>& es, 
	     TVertexDescriptor t) 
{
	SEQAN_CHECKPOINT
	es.data_target = t;
}


template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeAutomaton<TCargo, TSpec> >::Type&
target(EdgeAutomaton<TCargo, TSpec>& es) 
{
	SEQAN_CHECKPOINT
	return es.data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeAutomaton<TCargo, TSpec> >::Type
target(EdgeAutomaton<TCargo, TSpec> const& es) 
{
	SEQAN_CHECKPOINT
	return es.data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeAutomaton<TCargo, TSpec> const>::Type 
getTarget(EdgeAutomaton<TCargo, TSpec> const& es) 
{
	SEQAN_CHECKPOINT
	return es.data_target;
}

template<typename TCargo, typename TSpec>
inline typename VertexDescriptor<EdgeAutomaton<TCargo, TSpec> >::Type 
getTarget(EdgeAutomaton<TCargo, TSpec>& es) 
{
	SEQAN_CHECKPOINT
	return es.data_target;
}

template<typename TCargo, typename TSpec, typename TId>
void 
_assignId(EdgeAutomaton<TCargo, TSpec>& es, 
		  TId const id) 
{
	// Should be never called!!!
	SEQAN_ASSERT(false)
	return 0;
}

template<typename TCargo, typename TSpec>
inline typename Id<EdgeStump<TCargo, TSpec> const>::Type
_getId(EdgeAutomaton<TCargo, TSpec> const& es) 
{
	// Should be never called!!!
	SEQAN_ASSERT(false)
	return 0;
}

//_getId must be called on EdgeDescriptor
template<typename TVertexDescriptor, typename TAlphabet>
inline typename Id<TVertexDescriptor>::Type
_getId(Pair<TVertexDescriptor, TAlphabet> const& ed)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alph_size = ValueSize<TAlphabet>::VALUE;
	return (ed.i1 * alph_size + (TSize) ed.i2);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
