#ifndef SEQAN_HEADER_BLAST_HIT_ITERATOR_H
#define SEQAN_HEADER_BLAST_HIT_ITERATOR_H

//SEQAN_NO_DDDOC: do not generate documentation for this file

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Hit Iterator
//////////////////////////////////////////////////////////////////////////////




/**
.Spec.HitIterator:
..cat:Blast
..summary:Hit iterator for @Class.BlastReport@.
..signature:Iterator<TBlastReport, HitIterator>
..param.TBlastReport:A Blast report.
...type:Class.Blast
..general:Class.Iter
*/
template<typename TBlastReport>
class Iter<TBlastReport, SimpleBlastIterator<HitIterator> > 
{
public:

	TBlastReport * data_host;
	unsigned int data_pos;

	Iter()	
	{
	SEQAN_CHECKPOINT
	}
	
	Iter(TBlastReport & blast) : 
		data_host(&blast), 
		data_pos(0) 
	{
	SEQAN_CHECKPOINT
	}

	Iter(Iter const& it) : 
		data_host(it.data_host), 
		data_pos(it.data_pos) 
	{
	SEQAN_CHECKPOINT
	}

	~Iter() 
	{
	SEQAN_CHECKPOINT
	}

	Iter const&	operator = (Iter const & other) 
	{
	SEQAN_CHECKPOINT
		if (this == &other) return *this;
		data_host = other.data_host;
		data_pos = other.data_pos;
		return *this;
	}
//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// SimpleBlastIterator<TSpec> - Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.Blast
template<typename TBlastObject, typename TIteratorSpec>
struct Host<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >
{	
	typedef TBlastObject Type;
};



//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
struct Reference<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >
{
	typedef typename Value<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >::Type& Type;
};

template<typename TBlastObject, typename TIteratorSpec>
struct Reference<Iter<TBlastObject const, SimpleBlastIterator<TIteratorSpec> > >
{
	typedef typename Value<Iter<TBlastObject const, SimpleBlastIterator<TIteratorSpec> > >::Type& Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
struct GetValue<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >
{
	typedef typename Value<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >::Type Type;
};

template<typename TBlastObject, typename TIteratorSpec>
struct GetValue<Iter<TBlastObject const, SimpleBlastIterator<TIteratorSpec> > >
{
	typedef typename Value<Iter<TBlastObject const, SimpleBlastIterator<TIteratorSpec> > >::Type Type;
};



//////////////////////////////////////////////////////////////////////////////
// SimpleBlastIterator<HitIterator> - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.Blast

template<typename TBlastHsp, typename TInfoSpec>
struct Iterator<BlastReport<TBlastHsp,StoreReport<TInfoSpec> >, HitIterator>
{	
	typedef Iter<BlastReport<TBlastHsp,StoreReport<TInfoSpec> >, SimpleBlastIterator<HitIterator> > Type;
};

template<typename TBlastHsp, typename TInfoSpec>
struct Iterator<BlastReport<TBlastHsp,StoreReport<TInfoSpec> > const, HitIterator>
{	
	typedef Iter<BlastReport<TBlastHsp,StoreReport<TInfoSpec> > const, SimpleBlastIterator<HitIterator> > Type;
};




//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TSpec>
struct Value<Iter<BlastReport<TBlastHsp,TSpec>, SimpleBlastIterator<HitIterator> > >
{
	typedef typename Hit<BlastReport<TBlastHsp,TSpec> >::Type Type;
};

template<typename TBlastHsp, typename TSpec>
struct Value<Iter<BlastReport<TBlastHsp,TSpec> const, SimpleBlastIterator<HitIterator> > >
{
	typedef typename Hit<BlastReport<TBlastHsp,TSpec> const>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Blast SimpleBlastIterator<TSpec> - FUNCTIONS
// TSpecs: HspIterator and HitIterator
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline typename Reference<Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > >::Type
operator * (Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
	SEQAN_CHECKPOINT
	return value(it);
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.atBegin:
..cat:Blast
..summary:Determines whether the iterator is at the beginning or not.
..signature:atBegin(it)
..param.it:A hit or hsp iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:True if the iterator is at the beginning, false otherwise
..see:Function.goBegin
*/
template<typename TBlastObject, typename TIteratorSpec>
inline bool
atBegin(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == 0);	
}

//////////////////////////////////////////////////////////////////////////////



/**
.Function.goBegin:
..cat:Blast
..param.iterator
...type:Spec.HitIterator
...type:Spec.HspIterator
*/
template<typename TBlastObject, typename TIteratorSpec>
inline void
goBegin(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = 0;
}



//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..cat:Blast
..summary:Moves the iterator to the next hit or hsp.
..signature:goNext(it)
..param.it:A hit or hsp iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:void
..remarks:This method does nothing if the iterator is already at the end.
..see:Function.goPrevious
*/
template<typename TBlastObject, typename TIteratorSpec>
inline void
goNext(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(it)) ++it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >&
operator ++(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
SEQAN_CHECKPOINT
	goNext(it);
	return it;
}



//////////////////////////////////////////////////////////////////////////////

/**
.Function.goPrevious:
..cat:Blast
..summary:Moves the iterator to the preceding hit or hsp.
..signature:goPrevious(it)
..param.it:A hit or hsp iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:void
..remarks:This method does nothing if the iterator is already at the beginning.
..see:Function.goNext
*/

template<typename TBlastObject, typename TIteratorSpec>
inline void
goPrevious(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
	SEQAN_CHECKPOINT
	if (!atBegin(it)) --it.data_pos;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >&
operator --(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it)
{
SEQAN_CHECKPOINT
	goPrevious(it);
	return it;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >
operator --(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it, int)
{
	SEQAN_CHECKPOINT
	Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> > ret = it;
	goPrevious(it);
	return ret;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline bool
operator ==(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it1,
			Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos==it2.data_pos);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastObject, typename TIteratorSpec>
inline bool
operator !=(Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it1,
			Iter<TBlastObject, SimpleBlastIterator<TIteratorSpec> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos!=it2.data_pos);
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Blast SimpleBlastIterator<HitIterator> - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

/**
.Function.getValue:
..cat:Blast
..summary:The hit the iterator points to.
..signature:getValue(it)
..param.it:A Blast search result iterator.
...type:Spec.Hit Iterator
...type:Spec.Hsp Iterator
..returns:A BlastHit or BlastHsp.
...type:Metafunction.Hit
...type:Metafunction.Hsp
..see:Function.value
*/
template<typename TBlastReport>
inline typename GetValue<Iter<TBlastReport, SimpleBlastIterator<HitIterator> > >::Type
getValue(Iter<TBlastReport, SimpleBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	return it.data_host->hits[it.data_pos];
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.value:
..cat:Blast
..summary:The hit the iterator points to.
..signature:value(it)
..param.it:A Blast search result iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:A BlastHit or BlastHsp.
...type:Metafunction.Hit
...type:Metafunction.Hsp
..see:Function.getValue
*/
template<typename TBlastReport>
inline typename Reference<Iter<TBlastReport, SimpleBlastIterator<HitIterator> > >::Type 
value(Iter<TBlastReport, SimpleBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	return it.data_host->hits[it.data_pos];
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.hostReport:
..cat:Blast
..summary:The BlastReport this iterator is working on.
..signature:hostReport(it)
..param.it:A Blast hit iterator.
...type:Spec.HitIterator
..returns:A pointer to the host BlastReport.
*/
template<typename TBlastReport>
inline typename Host<Iter<TBlastReport, SimpleBlastIterator<HitIterator> > >::Type/* const*/ &
hostReport(Iter<TBlastReport, SimpleBlastIterator<HitIterator> > & it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
} 

//////////////////////////////////////////////////////////////////////////////


/**
.Function.atEnd:
..cat:Blast
..summary:Determines whether the iterator is at the end or not.
..signature:atEnd(it)
..param.it:A hit or hsp iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:True if the iterator is at the end, false otherwise
..see:Function.goEnd
*/
template<typename TBlastReport>
inline bool
atEnd(Iter<TBlastReport, SimpleBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == length(it.data_host->hits));	
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.goEnd:
..cat:Blast
..summary:Resets the iterator to the end.
..signature:goEnd(it)
..param.it:A hit or hsp iterator.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:void
..see:Function.atEnd
*/
template<typename TBlastReport>
inline void
goEnd(Iter<TBlastReport, SimpleBlastIterator<HitIterator> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = length(it.data_host->hits);
}


//////////////////////////////////////////////////////////////////////////////




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
