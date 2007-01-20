#ifndef SEQAN_HEADER_GRAPH_ADJACENCYITERATOR_H
#define SEQAN_HEADER_GRAPH_ADJACENCYITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph AdjacencyIterator
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > > 
{
public:
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, OutEdgeIterator<> >::Type TOutEdgeIterator;
	TOutEdgeIterator data_edge_it;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph, TVertexDescriptor const v) : 
		data_edge_it(_graph, v)
	{
		SEQAN_CHECKPOINT
	}
	
	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter(Iter const& _iter) : data_edge_it(_iter.data_edge_it)
	{
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_edge_it = _other.data_edge_it;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph AdjacencyIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TIteratorSpec>
struct Iterator<TGraph, AdjacencyIterator<TIteratorSpec> >
{	
	typedef Iter<TGraph, GraphIterator<AdjacencyIterator<TIteratorSpec> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Iterator<TGraph const, AdjacencyIterator<TIteratorSpec> >
{	
	typedef Iter<TGraph const, GraphIterator<AdjacencyIterator<TIteratorSpec> > > Type;
};


template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<AdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<AdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<AdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<AdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph const, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<AdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<AdjacencyIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph const, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////
// Graph AdjacencyIterator - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return getValue(it);
}

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

template<typename TGraph, typename TSpec>
inline typename HostGraph<Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return hostGraph(it.data_edge_it);
} 

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return atBegin(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goBegin(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (atEnd(it.data_edge_it));
}

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goEnd(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goNext(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

template<typename TGraph, typename TSpec>
inline void
goPrevious(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goPrevious(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >&
operator --(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >
operator --(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > > ret = it;
	goPrevious(it);
	return ret;
}

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_edge_it==it2.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<AdjacencyIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_edge_it!=it2.data_edge_it);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
