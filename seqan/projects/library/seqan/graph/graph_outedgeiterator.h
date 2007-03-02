#ifndef SEQAN_HEADER_GRAPH_OUTEDGEITERATOR_H
#define SEQAN_HEADER_GRAPH_OUTEDGEITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// OutEdgeIterator for EdgeList
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
class Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > > 
{
public:
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> TGraph;
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
// OutEdgeIterator for EdgeListU
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
class Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > > 
{
public:
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> TGraph;
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
// Graph OutEdgeIterator for Automaton
//////////////////////////////////////////////////////////////////////////////
template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIteratorSpec>
class Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > 
{
public:
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TAlphabet>::Type TSize;
	TGraph const* data_host;
	TVertexDescriptor data_source;
	TSize data_pos;
	TSize data_begin;
	TSize data_end;


	Iter()	
	{
		SEQAN_CHECKPOINT
	}

	Iter(TGraph const& _graph, TVertexDescriptor const v) : 
		data_host(&_graph),
		data_source(v)
	{
		SEQAN_CHECKPOINT
		TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
		TSize table_length = ValueSize<TAlphabet>::VALUE;
		TSize pos = 0;
		while (	(pos < table_length) &&
				(_graph.data_vertex[v].data_edge[pos].data_target == nilVal))
		{
				++pos;
		}
		data_pos = pos;
		data_begin = pos;
		data_end = table_length;
	}

	Iter(Iter const& _iter) : 
		data_host(_iter.data_host),
		data_source(_iter.data_source),
		data_pos(_iter.data_pos),
		data_begin(_iter.data_begin),
		data_end(_iter.data_end)
	{
		SEQAN_CHECKPOINT
	}

	~Iter() 
	{
		SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & _other) 
	{
		SEQAN_CHECKPOINT
		if (this == &_other) return *this;
		data_host = _other.data_host;
		data_source = _other.data_source;
		data_pos = _other.data_pos;
		data_begin = _other.data_begin;
		data_end = _other.data_end;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Graph OutEdgeIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Iterator<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, OutEdgeIterator<TIteratorSpec> >
{	
	typedef Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Iterator<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> const, OutEdgeIterator<TIteratorSpec> >
{	
	typedef Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Iterator<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, OutEdgeIterator<TIteratorSpec> >
{	
	typedef Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Iterator<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> const, OutEdgeIterator<TIteratorSpec> >
{	
	typedef Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIteratorSpec>
struct Iterator<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>, OutEdgeIterator<TIteratorSpec> >
{	
	typedef Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TIteratorSpec>
struct Iterator<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const, OutEdgeIterator<TIteratorSpec> >
{	
	typedef Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > Type;
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

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Reference<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Reference<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Reference<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Reference<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Reference<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct Reference<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct GetValue<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct GetValue<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct GetValue<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct GetValue<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type& Type;
};


template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct GetValue<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TIteratorSpec>
struct GetValue<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >
{
	typedef typename Value<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec> const, GraphIterator<OutEdgeIterator<TIteratorSpec> > > >::Type Type;
};



//////////////////////////////////////////////////////////////////////////////
// Graph OutIterator - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename GetValue<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
getValue(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	return TEdgeDescriptor(it.data_source,TAlphabet(it.data_pos));
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_edge;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename Reference<Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > > >::Type
value(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	return TEdgeDescriptor(it.data_source,TAlphabet(it.data_pos));
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

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == getValue(it.data_host->data_vertex, it.data_source));
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == getValue(it.data_host->data_vertex, it.data_source));
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
atBegin(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_begin == it.data_pos);
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = getValue(it.data_host->data_vertex,it.data_source);
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = getValue(it.data_host->data_vertex,it.data_source);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goBegin(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = it.data_begin;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == 0);
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_edge == 0);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
atEnd(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return (it.data_end == it.data_pos);
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = 0;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_edge = 0;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goEnd(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = it.data_end;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) it.data_edge = it.data_edge->data_next;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) {
		if (it.data_source == getSource(it.data_edge)) it.data_edge = it.data_edge->data_next_source;
		else it.data_edge = it.data_edge->data_next_target;
	}
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goNext(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	if (it.data_pos < it.data_end) ++it.data_pos;
	while (	(it.data_pos < it.data_end) &&
			(it.data_host->data_vertex[it.data_source].data_edge[it.data_pos].data_target == nilVal))
	{
				++it.data_pos;
	}
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

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef typename EdgeType<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> >::Type TEdge;
	TEdge* current = getValue(it.data_host->data_vertex, it.data_source);
	if (current == it.data_edge) return;
	while ((current != 0) && (current->data_next != it.data_edge)) current = current->data_next;
	it.data_edge = current;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> TGraph;
	typedef typename EdgeType<TGraph>::Type TEdge;
	TEdge* current = getValue(it.data_host->data_vertex, it.data_source);
	if (current == it.data_edge) return;
	while (current != 0) {
		if (it.data_source == getSource(current)) {
			if (it.data_edge == current->data_next_source) break;
			else current = current->data_next_source;
		} else {
			if (it.data_edge == current->data_next_target) break;
			else current = current->data_next_target;
		}
	}
	it.data_edge = current;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline void
goPrevious(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVal = _get_nil<TVertexDescriptor>();
	if (it.data_pos != it.data_begin) --it.data_pos;
	while (	(it.data_pos != it.data_begin) &&
			(it.data_host->data_vertex[it.data_source].data_edge[it.data_pos].data_target == nilVal))
	{
				--it.data_pos;
	}
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

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_edge==it2.data_edge) && 
			(it1.data_source==it2.data_source));
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_edge==it2.data_edge) && 
			(it1.data_source==it2.data_source));
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
operator ==(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_pos==it2.data_pos) && 
			(it1.data_source==it2.data_source));
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_edge!=it2.data_edge) || 
			(it1.data_source!=it2.data_source));
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_edge!=it2.data_edge) || 
			(it1.data_source!=it2.data_source));
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline bool
operator !=(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it1,
			Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it2)
{
SEQAN_CHECKPOINT
	return ((it1.data_pos!=it2.data_pos) || 
			(it1.data_source!=it2.data_source));
}


// Functions for speed-up
template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> >::Type 
sourceVertex(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> , GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_source;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> >::Type 
sourceVertex(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> , GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_source;
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec> >::Type
sourceVertex(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_source;
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> >::Type 
targetVertex(Iter<Graph<EdgeList<TCargo, TEdgeSpec>, TGraphSpec> , GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return targetVertex(*it.data_host, it.data_edge);
}

template<typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> >::Type 
targetVertex(Iter<Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> , GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	typedef Graph<EdgeListU<TCargo, TEdgeSpec>, TGraphSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor target = targetVertex(*it.data_host, it.data_edge);
	if (target != it.data_source) return target;
	else return sourceVertex(*it.data_host, it.data_edge);
}

template<typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TGraphSpec, typename TSpec>
inline typename VertexDescriptor<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec> >::Type
targetVertex(Iter<Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TGraphSpec>, GraphIterator<OutEdgeIterator<TSpec> > >& it)
{
	SEQAN_CHECKPOINT
	return it.data_host->data_vertex[it.data_source].data_edge[it.data_pos].data_target;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
