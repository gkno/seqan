#ifndef SEQAN_HEADER_GRAPH_BFSITERATOR_H
#define SEQAN_HEADER_GRAPH_BFSITERATOR_H

#include <deque>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph BfsIterator
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
class Iter<TGraph, GraphIterator<BfsIterator<TSpec> > > 
{
public:
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph const* data_host;
	TVertexDescriptor data_source;
	String<bool> data_tokenMap;
	std::deque<TVertexDescriptor> data_queue;

	void _init() {
		initVertexMap(*data_host,data_tokenMap);
		typedef typename Iterator<String<bool> >::Type TIter;
		TIter it = begin(data_tokenMap);
		for(;!atEnd(it);goNext(it)) {
			assignValue(it,false);
		}
		assignProperty(data_tokenMap, data_source, true);
		data_queue.clear();
		data_queue.push_back(data_source);
	}

	Iter()
	{
		SEQAN_CHECKPOINT
	}
	
	Iter(TGraph const& _graph, TVertexDescriptor v) : 
		data_host(&_graph),
		data_source(v)
	{
		SEQAN_CHECKPOINT
		_init();
	}
	
	
	~Iter() {
		SEQAN_CHECKPOINT
	}

	Iter(Iter const& _iter) :
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_tokenMap(_iter.data_tokenMap),
		data_queue(_iter.data_queue)
	{
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) {
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host=_other.data_host;
		data_source=_other.data_source;
		data_tokenMap=_other.data_tokenMap;
		data_queue=_other.data_queue;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph BfsIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TIteratorSpec>
struct Iterator<TGraph, BfsIterator<TIteratorSpec> >
{	
	typedef Iter<TGraph, GraphIterator<BfsIterator<TIteratorSpec> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Iterator<TGraph const, BfsIterator<TIteratorSpec> >
{	
	typedef Iter<TGraph const, GraphIterator<BfsIterator<TIteratorSpec> > > Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph, GraphIterator<BfsIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Value<Iter<TGraph const, GraphIterator<BfsIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<TGraph const, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph, GraphIterator<BfsIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct Reference<Iter<TGraph const, GraphIterator<BfsIterator<TIteratorSpec> > > >
{
	typedef typename Reference<Iter<TGraph const, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph, GraphIterator<BfsIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TGraph, typename TIteratorSpec>
struct GetValue<Iter<TGraph const, GraphIterator<BfsIterator<TIteratorSpec> > > >
{
	typedef typename GetValue<Iter<TGraph const, GraphIterator<VertexIterator<TIteratorSpec> > > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////
// Graph BfsIterator - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<BfsIterator<TSpec> > > >::Type
getValue(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_queue.front();
}

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<BfsIterator<TSpec> > > >::Type
value(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return getValue(it);
}

template<typename TGraph, typename TSpec>
inline typename GetValue<Iter<TGraph, GraphIterator<BfsIterator<TSpec> > > >::Type
operator * (Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

template<typename TGraph, typename TSpec>
inline typename HostGraph<Iter<TGraph, GraphIterator<BfsIterator<TSpec> > > >::Type const&
hostGraph(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
}

template<typename TGraph, typename TSpec>
inline bool
atBegin(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (it.data_queue.empty()) return false;
	else return (it.data_queue.front() == it.data_source);
}

template<typename TGraph, typename TSpec>
inline void
goBegin(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it._init();
}

template<typename TGraph, typename TSpec>
inline bool
atEnd(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_queue.empty());
}

template<typename TGraph, typename TSpec>
inline void
goEnd(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_queue.clear();
}

template<typename TGraph, typename TSpec>
inline void
goNext(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (it.data_queue.empty()) return;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor u = it.data_queue.front();
	it.data_queue.pop_front();
	typedef typename Iterator<TGraph, AdjacencyIterator<> >::Type TAdjacencyIterator;
	TAdjacencyIterator itad(*it.data_host,u);
	for(;!atEnd(itad);goNext(itad)) {
		TVertexDescriptor v = getValue(itad);
		if (getProperty(it.data_tokenMap, v) == false) {
			assignProperty(it.data_tokenMap, v, true);
			it.data_queue.push_back(v);
		}
	}
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >&
operator ++(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}

template<typename TGraph, typename TSpec>
inline Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >
operator ++(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TGraph, GraphIterator<BfsIterator<TSpec> > > ret = it;
	goNext(it);
	return ret;
}

template<typename TGraph, typename TSpec>
inline bool
operator ==(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_source==it2.data_source) &&
			(it1.data_tokenMap==it2.data_tokenMap) &&
			(it1.data_queue==it2.data_queue));
}

template<typename TGraph, typename TSpec>
inline bool
operator !=(Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it1,
			Iter<TGraph, GraphIterator<BfsIterator<TSpec> > >& it2)
{
	SEQAN_CHECKPOINT
	return ((it1.data_source!=it2.data_source) ||
			(it1.data_tokenMap!=it2.data_tokenMap) ||
			(it1.data_queue!=it2.data_queue));
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
