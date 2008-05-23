 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_iterator_vertex.h 1757 2008-02-27 16:26:20Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ITERATOR_VERTEX_H
#define SEQAN_HEADER_GRAPH_ITERATOR_VERTEX_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph VertexIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Vertex Iterator:
..cat:Graph
..summary:Vertex iterator for @Class.Graph@.
..signature:Iterator<TGraph, VertexIterator>
..param.TGraph:A graph.
...type:Class.Graph
..general:Class.Iter
..see:Spec.Edge Iterator
..see:Spec.Out-Edge Iterator
..see:Spec.Adjacency Iterator
..see:Spec.Bfs Iterator
..see:Spec.Dfs Preorder Iterator
*/
template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > > 
{
public:
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph const* data_host;
	TVertexDescriptor data_pos;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph) : 
		data_host(&_graph), 
		data_pos(getIdLowerBound(_getVertexIdManager(*data_host))) 
	{
		SEQAN_CHECKPOINT
	}

	Iter(Iter const& _iter) : 
		data_host(_iter.data_host), 
		data_pos(_iter.data_pos) 
	{
		SEQAN_CHECKPOINT
	}

	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host = _other.data_host;
		data_pos = _other.data_pos;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Graph InternalVertexIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.Graph

template<typename TGraph>
struct Iterator<TGraph, VertexIterator>
{	
	typedef Iter<TGraph, GraphIterator<InternalVertexIterator<VertexIterator> > > Type;
};

template<typename TGraph>
struct Iterator<TGraph const, VertexIterator>
{	
	typedef Iter<TGraph const, GraphIterator<InternalVertexIterator<VertexIterator> > > Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >
{
	typedef typename VertexDescriptor<TGraph>::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >
{
	typedef typename VertexDescriptor<TGraph const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type& Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Spec<Iter<TGraph const, GraphIterator<InternalVertexIterator<TIteratorSpec> > > >
{
	typedef TIteratorSpec Type;
};


//////////////////////////////////////////////////////////////////////////////
// Graph InternalVertexIterator - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/*DISABLED: see sequence_interface.h and basic_iterator.h
.Function.getValue:
..cat:Graph
..summary:The vertex or edge the iterator points to.
..signature:getValue(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..returns:A vertex descriptor or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..see:Function.value
*/
/**
.Function.getValue:
..cat:Graph
..param.object
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
*/

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	return it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Function.value:
..cat:Graph
..summary:The vertex the iterator points to.
..signature:value(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..returns:A vertex descriptor or edge descriptor.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..see:Function.getValue
*/
/**
.Function.value:
..cat:Graph
..param.object
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
*/

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	return it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.hostGraph:
..cat:Graph
..summary:The graph this iterator is working on.
..signature:hostGraph(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..returns:A pointer to the host graph.
*/

template<typename TGraph, typename TSpec>
inline typename Host<Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
} 

//////////////////////////////////////////////////////////////////////////////

/*DISABLED: see basic_iterator.h
.Function.atBegin:
..cat:Graph
..summary:Determines whether the iterator is at the beginning or not.
..signature:atBegin(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..returns:True if the iterator is at the beginning, false otherwise
..see:Function.goBegin
*/
/**
.Function.atBegin:
..cat:Graph
..param.iterator
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..see:Function.goBegin
*/

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	return (getValue(it) == getIdLowerBound(_getVertexIdManager(*it.data_host)));	
}

//////////////////////////////////////////////////////////////////////////////

/*DISABLED: see basic_iterator.h
.Function.goBegin:
..cat:Graph
..summary:Resets the iterator to the beginning.
..signature:goBegin(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..returns:void
..see:Function.atBegin
*/
/**
.Function.goBegin:
..cat:Graph
..param.iterator
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
*/

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = getIdLowerBound(_getVertexIdManager(*it.data_host));
}

//////////////////////////////////////////////////////////////////////////////

/*DISABLED: see basic_iterator.h
.Function.atEnd:
..cat:Graph
..summary:Determines whether the iterator is at the end or not.
..signature:atEnd(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..returns:True if the iterator is at the end, false otherwise
..see:Function.goEnd
*/
/**
.Function.atEnd:
..cat:Graph
..param.iterator
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..see:Function.goEnd
*/

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	return (getValue(it) >= getIdUpperBound(_getVertexIdManager(*it.data_host)));	
}

//////////////////////////////////////////////////////////////////////////////

/*DISABLED: see basic_iterator.h
.Function.goEnd:
..cat:Graph
..summary:Resets the iterator to the end.
..signature:goEnd(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..returns:void
..see:Function.atEnd
*/
/**
.Function.goEnd:
..cat:Graph
..param.iterator
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
*/
template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = getIdUpperBound(_getVertexIdManager(*it.data_host));
}

//////////////////////////////////////////////////////////////////////////////

/*DISABLED: see basic_iterator.h
.Function.goNext:
..cat:Graph
..summary:Moves the iterator to the next vertex or next edge.
..signature:goNext(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
..returns:void
..remarks:This method does nothing if the iterator is already at the end.
..see:Function.goPrevious
*/
/**
.Function.goNext:
..cat:Graph
..param.iterator
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
...type:Spec.Bfs Iterator
...type:Spec.Dfs Preorder Iterator
*/

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) ++it.data_pos;
	while ((!atEnd(it)) && (!idInUse(_getVertexIdManager(*it.data_host), it.data_pos))) ++it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

/*DISABLED: see basic_iterator.h
.Function.goPrevious:
..cat:Graph
..summary:Moves the iterator to the preceding vertex or the preceding edge.
..signature:goPrevious(it)
..param.it:A vertex or edge iterator.
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
..returns:void
..remarks:This method does nothing if the iterator is already at the beginning.
..see:Function.goNext
*/
/**
.Function.goPrevious:
..cat:Graph
..param.iterator
...type:Spec.Vertex Iterator
...type:Spec.Out-Edge Iterator
...type:Spec.Edge Iterator
...type:Spec.Adjacency Iterator
*/

template<typename TGraph, typename TSpec>
inline void
goPrevious(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atBegin(it)) --it.data_pos;
	while ((!atBegin(it)) && (!idInUse(_getVertexIdManager(*it.data_host), it.data_pos))) --it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >&
operator --(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >
operator --(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > > ret = it;
	goPrevious(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos==it2.data_pos);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<InternalVertexIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos!=it2.data_pos);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
