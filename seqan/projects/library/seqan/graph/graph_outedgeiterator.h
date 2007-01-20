#ifndef SEQAN_HEADER_GRAPH_OUTEDGEITERATOR_H
#define SEQAN_HEADER_GRAPH_OUTEDGEITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph OutEdgeIterator
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > > 
{
public:
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph const* data_host;
	TVertexDescriptor data_source;
	TEdgeDescriptor data_edge;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph, TVertexDescriptor const v) : 
		data_host(&_graph),
		data_source(v),
		data_edge(getValue(_graph.data_vertex,v))
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(Iter const& _iter) : 
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_edge(_iter.data_edge)
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
		data_source = _other.data_source;
		data_edge = _other.data_edge;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph OutEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TIteratorSpec>
struct Iterator<TGraph, OutEdgeIterator<TIteratorSpec> >
{	
	typedef Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Iterator<TGraph const, OutEdgeIterator<TIteratorSpec> >
{	
	typedef Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename EdgeDescriptor<TGraph >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename EdgeDescriptor<TGraph const>::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};




//////////////////////////////////////////////////////////////////////////////
// Graph OutIterator - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

template<typename TGraph, typename TSpec>
inline typename HostGraph<Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
}


template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == getValue(it.data_host->data_vertex, it.data_source));
}

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = getValue(it.data_host->data_vertex,it.data_source);
}

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == 0);
}

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = 0;
}

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) it.data_edge = it.data_edge->data_next;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

template<typename TGraph, typename TSpec>
inline void
goPrevious(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef typename EdgeType<TGraph>::Type TEdge;
	TEdge* current = getValue(it.data_host->data_vertex, it.data_source);
	if (current == it.data_edge) return;
	while (current->data_next != it.data_edge) current = current->data_next;
	it.data_edge = current;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >&
operator --(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >
operator --(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > > ret = it;
	goPrevious(it);
	return ret;
}

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_edge==it2.data_edge) && 
			(it1.data_source==it2.data_source));
}

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_edge!=it2.data_edge) || 
			(it1.data_source!=it2.data_source));
}


// Functions for speed-up
template<typename TGraph, typename TSpec>
inline typename VertexDescriptor<TGraph>::Type 
sourceVertex(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_source;
}

template<typename TGraph, typename TSpec>
inline typename VertexDescriptor<TGraph>::Type 
targetVertex(Iter<TGraph, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(*it.data_host, it.data_edge);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
