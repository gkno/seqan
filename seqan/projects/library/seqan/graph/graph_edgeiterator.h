#ifndef SEQAN_HEADER_GRAPH_EDGEITERATOR_H
#define SEQAN_HEADER_GRAPH_EDGEITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph EdgeIterator
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > > 
{
public:
	typedef typename Iterator<TGraph, VertexIterator<> >::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator<> >::Type TOutEdgeIterator;
	TVertexIterator data_vertex_it;
	TOutEdgeIterator data_edge_it;

	Iter()	
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph) : 
		data_vertex_it(_graph),
		data_edge_it(_graph, getIdLowerBound(_graph.data_id_managerV))  
	{
		SEQAN_CHECKPOINT
		while((!atEnd(data_vertex_it)) && 
				(atEnd(data_edge_it))) 
		{
				++data_vertex_it;
				typedef typename Iterator<TGraph, OutEdgeIterator<> >::Type TOutEdgeIterator;
				data_edge_it = TOutEdgeIterator(hostGraph(*this), getValue(data_vertex_it));			
		}
	}

	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter(Iter const& _iter) : 
		data_vertex_it(_iter.data_vertex_it),
		data_edge_it(_iter.data_edge_it)
	{
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_vertex_it = _other.data_vertex_it;
		data_edge_it = _other.data_edge_it;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph EdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TIteratorSpec>
struct Iterator<TGraph, EdgeIterator<TIteratorSpec> >
{	
	typedef Iter<TGraph, GraphIterator<EdgeIterator<TIteratorSpec> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Iterator<TGraph const, EdgeIterator<TIteratorSpec> >
{	
	typedef Iter<TGraph const, GraphIterator<EdgeIterator<TIteratorSpec> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<EdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<EdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<EdgeIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<EdgeIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<EdgeIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<EdgeIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////
// Graph EdgeIterator - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	return getValue(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	return value(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline typename Reference<Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

template<typename TGraph, typename TSpec>
inline typename HostGraph<Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return hostGraph(it.data_vertex_it);
} 

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (atBegin(it.data_vertex_it) &&
			atBegin(it.data_edge_it));
}

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef typename Iterator<TGraph, OutEdgeIterator<> >::Type TOutEdgeIterator;
	goBegin(it.data_vertex_it);
	it.data_edge_it = TOutEdgeIterator(hostGraph(it), getIdLowerBound(hostGraph(it).data_id_managerV));	
	while((!atEnd(it.data_vertex_it)) && 
			(atEnd(it.data_edge_it))) 
	{
		++it.data_vertex_it;
		it.data_edge_it = TOutEdgeIterator(hostGraph(it), getValue(it.data_vertex_it));			
	}
}

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return atEnd(it.data_vertex_it);
}

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goEnd(it.data_edge_it);
	goEnd(it.data_vertex_it);
}

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) {
		++it.data_edge_it;
		if(!atEnd(it.data_edge_it)) return;
		else {
			while((!atEnd(it.data_vertex_it)) && 
					(atEnd(it.data_edge_it))) {
				++it.data_vertex_it;
				if ((atEnd(it.data_vertex_it)) ||
					(!idInUse(hostGraph(it).data_id_managerV, getValue(it.data_vertex_it)))) continue;
				typedef typename Iterator<TGraph, OutEdgeIterator<> >::Type TOutEdgeIterator;
				it.data_edge_it = TOutEdgeIterator(hostGraph(it), getValue(it.data_vertex_it));			
			}
		}
	}
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

template<typename TGraph, typename TSpec>
inline void
goPrevious(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atBegin(it)) {
		if ((atEnd(it)) ||
			(atBegin(it.data_edge_it))) 
		{
			--it.data_vertex_it;
			typedef typename Iterator<TGraph, OutEdgeIterator<> >::Type TOutEdgeIterator;
			it.data_edge_it = TOutEdgeIterator(hostGraph(it), getValue(it.data_vertex_it));
			while((!atBegin(it.data_vertex_it)) && 
					(atEnd(it.data_edge_it))) {
				--it.data_vertex_it;
				typedef typename Iterator<TGraph, OutEdgeIterator<> >::Type TOutEdgeIterator;
				it.data_edge_it = TOutEdgeIterator(hostGraph(it), getValue(it.data_vertex_it));			
			}
			goEnd(it.data_edge_it);
		}
		--it.data_edge_it;
	}
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >&
operator --(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >
operator --(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > > ret = it;
	goPrevious(it);
	return ret;
}

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_vertex_it==it2.data_vertex_it) && 
			(it1.data_edge_it==it2.data_edge_it));
}

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_vertex_it!=it2.data_vertex_it) || 
			(it1.data_edge_it!=it2.data_edge_it));
}

// Functions for speed-up
template<typename TGraph, typename TSpec>
inline typename VertexDescriptor<TGraph>::Type 
sourceVertex(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return sourceVertex(it.data_edge_it);
}

template<typename TGraph, typename TSpec>
inline typename VertexDescriptor<TGraph>::Type 
targetVertex(Iter<TGraph, GraphIterator<EdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(it.data_edge_it);
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
