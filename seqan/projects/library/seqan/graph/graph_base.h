#ifndef SEQAN_HEADER_GRAPH_BASE_H
#define SEQAN_HEADER_GRAPH_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// General Metafunctions
//////////////////////////////////////////////////////////////////////////////

// IdType	
template<typename T>
struct Id {
	typedef unsigned int Type;
};

template<typename T>
struct Id<T const> {
	typedef unsigned int Type;
};

// VertexDescriptor = Id (for all graphs)
template<typename TGraph>
struct VertexDescriptor {
	typedef typename Id<TGraph>::Type Type;
};

// EdgeDescriptor: Varies for different graph types
template<typename T>
struct EdgeDescriptor;


//////////////////////////////////////////////////////////////////////////////
// Default IdManager
//////////////////////////////////////////////////////////////////////////////
template<typename TIdType = unsigned int, typename TSpec = Default>
class IdManager;

// Tag for no edge id
struct WithoutEdgeId_;
typedef Tag<WithoutEdgeId_> const WithoutEdgeId;


//////////////////////////////////////////////////////////////////////////////
// Default Edge Stump: No cargo but edge id
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo = void, typename TSpec = Default>
class EdgeStump;


//////////////////////////////////////////////////////////////////////////////
// Default Automaton Edge: No cargo and no edge id
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo = void, typename TSpec = Default>
class EdgeAutomaton;


//////////////////////////////////////////////////////////////////////////////
// Default EdgeList: Uses a list of edge stumps
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo = void, typename TSpec = Default>
struct EdgeList;


//////////////////////////////////////////////////////////////////////////////
// Default Automaton
//////////////////////////////////////////////////////////////////////////////
template<typename TAlphabet = char, typename TCargo = void, typename TSpec = Default>
struct Automaton;


//////////////////////////////////////////////////////////////////////////////
// Default Graph
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges = EdgeList<>, typename TSpec = Default>
class Graph;


//////////////////////////////////////////////////////////////////////////////
// Graph - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec>
struct Size<Graph<TEdges, TSpec> > {
	typedef size_t Type;
};

template<typename TEdges, typename TSpec>
struct Size<Graph<TEdges, TSpec> const> : Size<Graph<TEdges, TSpec> >
{
};


template<typename T>
struct Alphabet;

template<typename TAlphabet, typename TCargo, typename TEdgeSpec>
struct Alphabet<Automaton<TAlphabet, TCargo, TEdgeSpec> > {
	typedef TAlphabet Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec>
struct Alphabet<Automaton<TAlphabet, TCargo, TEdgeSpec> const> {
	typedef TAlphabet Type;
};

template<typename T>
struct EdgeType;

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeStump<TCargo, TEdgeSpec> Type;
};

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeStump<TCargo, TEdgeSpec> const Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> > {
	typedef EdgeAutomaton<TCargo, TEdgeSpec> Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeType<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const> {
	typedef EdgeAutomaton<TCargo, TEdgeSpec> const Type;
};

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> > 
{
	typedef typename EdgeType<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> >::Type* Type;
};

template<typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const>
{
	typedef typename EdgeType<Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const>::Type* Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> > 
{
	typedef Pair<typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type, TAlphabet> Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec>
struct EdgeDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const>
{
	typedef Pair<typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> >::Type, TAlphabet> Type;
};

template<typename T, typename TIdType>
struct IdHandler {
	typedef IdManager<TIdType> Type;
};

template<typename TCargo, typename TIdType>
struct IdHandler<EdgeStump<TCargo, WithoutEdgeId>, TIdType> {
	// Dummy IdManager
	typedef IdManager<void> Type;
};


//////////////////////////////////////////////////////////////////////////////
// Utility functions
//////////////////////////////////////////////////////////////////////////////
template <typename T>
inline T const
_get_nil(T *)
{
	return ~0;
}

template <typename T>
inline T const
_get_nil()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return _get_nil(_tag);
}

template <typename T>
inline T const
_get_infinity()
{
SEQAN_CHECKPOINT
	T * _tag = 0;
	return supremumValueImpl(_tag);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
