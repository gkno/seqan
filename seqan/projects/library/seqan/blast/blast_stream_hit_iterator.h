#ifndef SEQAN_HEADER_BLAST_STREAM_HIT_ITERATOR_H
#define SEQAN_HEADER_BLAST_STREAM_HIT_ITERATOR_H

//SEQAN_NO_DDDOC: do not generate documentation for this file

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Stream Hit Iterator
//////////////////////////////////////////////////////////////////////////////




/**
.Spec.HitIterator:
..cat:Blast
..summary:Hit iterator for @Class.BlastReport@.
..signature:Iterator<TBlastReport, HitIterator>
..param.TBlastReport:A Blast report with @Spec.StreamReport@.
...type:Class.Blast
..general:Class.Iter
*/
template<typename TBlastHsp, typename TFile>
class Iter<BlastReport<TBlastHsp, StreamReport<TFile> >, StreamBlastIterator<HitIterator> > 
{
public:
	typedef BlastReport<TBlastHsp, StreamReport<TFile> > TBlastReport;
	typedef typename Hit<TBlastReport>::Type TBlastHit;
	typedef typename Position<TFile>::Type TPosition;

	TBlastHit data_hit;
	TBlastReport* data_host;
	TPosition data_pos, data_next_pos;
	bool data_at_end;

	Iter()	
	{
	SEQAN_CHECKPOINT
	}
	
	Iter(TBlastReport & blast)  
	{
	SEQAN_CHECKPOINT
		data_host = &blast; 
		data_pos = blast.first_hit_pos;
		data_next_pos = data_pos;
		if(blast.hits_found)
			data_at_end = false;
		else
			data_at_end = true;
		data_hit.data_host = &blast;

	}

	Iter(Iter const& other) : 
		data_host(other.data_host), 
		data_pos(other.data_pos), 
		data_next_pos(other.data_next_pos), 
		data_at_end(other.data_at_end),
		data_hit(other.data_hit)
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
		data_next_pos = other.data_next_pos;
		data_hit = other.data_hit;
		return *this;
	}
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Blast StreamHitIterator - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Host.param.T.type:Class.Blast
template<typename TBlastReport>
struct Host<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >
{	
	typedef TBlastReport Type;
};


///.Metafunction.Iterator.param.T.type:Class.Blast

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastReport<TBlastHsp,StreamReport<TFile> >, HitIterator>
{	
	typedef Iter<BlastReport<TBlastHsp,StreamReport<TFile> >, StreamBlastIterator<HitIterator> > Type;
};

template<typename TBlastHsp, typename TFile>
struct Iterator<BlastReport<TBlastHsp,StreamReport<TFile> > const, HitIterator>
{	
	typedef Iter<BlastReport<TBlastHsp,StreamReport<TFile> > const, StreamBlastIterator<HitIterator> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TSpec>
struct Value<Iter<BlastReport<TBlastHsp,TSpec>, StreamBlastIterator<HitIterator> > >
{
	typedef typename Hit<BlastReport<TBlastHsp,TSpec> >::Type Type;
};

template<typename TBlastHsp, typename TSpec>
struct Value<Iter<BlastReport<TBlastHsp,TSpec> const, StreamBlastIterator<HitIterator> > >
{
	typedef typename Hit<BlastReport<TBlastHsp,TSpec> const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastReport>
struct Reference<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >
{
	typedef typename Value<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type& Type;
};

template<typename TBlastReport>
struct Reference<Iter<TBlastReport const, StreamBlastIterator<HitIterator> > >
{
	typedef typename Value<Iter<TBlastReport const, StreamBlastIterator<HitIterator> > >::Type& Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastReport>
struct GetValue<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >
{
	typedef typename Value<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type Type;
};

template<typename TBlastReport>
struct GetValue<Iter<TBlastReport const, StreamBlastIterator<HitIterator> > >
{
	typedef typename Value<Iter<TBlastReport const, StreamBlastIterator<HitIterator> > >::Type Type;
};

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Blast StreamBlastIterator<HitIterator> - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


/**
.Function.atBegin:
..cat:Blast
..summary:Determines whether the iterator is at the beginning or not.
..signature:atBegin(file,it)
..param.it:A hit or hsp iterator on a Blast report with @Spec.StreamReport@.
..param.it:A stream.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:True if the iterator is at the beginning, false otherwise
..see:Function.goBegin
*/
template<typename TBlastReport, typename TFile>
inline bool
atBegin(TFile &,
		Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	return (it.data_pos == it.data_host->first_hit_pos);	
}

//////////////////////////////////////////////////////////////////////////////



/**
.Function.goBegin:
..cat:Blast
..summary:Resets the iterator to the beginning.
..signature:goBegin(file,it)
..param.it:A hit or hsp iterator on a Blast report with @Spec.StreamReport@.
..param.it:A stream.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:void
..see:Function.atBegin
*/
template<typename TBlastReport, typename TFile>
inline void
goBegin(TFile &,
		Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
	SEQAN_CHECKPOINT
	it.data_pos = it.data_host->first_hit_pos;
	it.data_next_pos = it.data_pos;
	if(it.data_host->hits_found)
		it.data_at_end = false;
	else
		it.data_at_end = true;
}



//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..cat:Blast
..summary:Moves the iterator to the next hit or hsp.
..signature:goNext(file,it)
..param.it:A hit or hsp iterator on a Blast report with @Spec.StreamReport@.
..param.it:A stream.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:void
..remarks:This method does nothing if the iterator is already at the end.
..see:Function.goPrevious
*/
template<typename TBlastReport, typename TFile>
inline void
goNext(TFile & file, 
	   Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
	SEQAN_CHECKPOINT
	if (!atEnd(file,it)) 
	{
		if(it.data_pos == it.data_next_pos)
			getNextHitFilePos(file,it);
		if(it.data_pos == it.data_next_pos)
			it.data_at_end = true;
		else
            it.data_pos = it.data_next_pos;
	}
}

//////////////////////////////////////////////////////////////////////////////


//template<typename TBlastReport>
//inline void
//goPrevious(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
//{
//	SEQAN_CHECKPOINT
//	if (!atBegin(it)) --it.data_pos;
//}
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
//
//template<typename TBlastReport>
//inline Iter<TBlastReport, StreamBlastIterator<HitIterator> >
//operator --(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it, int)
//{
//	SEQAN_CHECKPOINT
//	Iter<TBlastReport, StreamBlastIterator<HitIterator> > ret = it;
//	goPrevious(it);
//	return ret;
//}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastReport>
inline bool
operator ==(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it1,
			Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos==it2.data_pos && it1.data_host==it2.data_host);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TBlastReport>
inline bool
operator !=(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it1,
			Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it2)
{
SEQAN_CHECKPOINT
	return (it1.data_pos!=it2.data_pos && it1.data_host!=it2.data_host);
}

//////////////////////////////////////////////////////////////////////////////


/**
.Function.getValue:
..cat:Blast
..summary:The hit the iterator points to.
..signature:getValue(file,it)
..param.it:A hit or hsp iterator on a Blast report with @Spec.StreamReport@.
..param.it:A stream.
...type:Spec.Hit Iterator
...type:Spec.Hsp Iterator
..returns:A BlastHit or BlastHsp.
...type:Metafunction.Hit
...type:Metafunction.Hsp
..see:Function.value
*/
template<typename TBlastReport, typename TFile>
inline typename GetValue<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type
getValue(TFile & file,
		 Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	if(it.data_pos != (it.data_hit).begin_pos)
	{
		_streamSeekG(file,it.data_pos);
		it.data_host->act_c = '>';
		_parseBlastHit(file,it.data_host->act_c,it.data_hit);
		//if(_parseBlastHit(file,it.data_host->act_c,it.data_hit) == it.data_pos)
		//	it.data_at_end = true;
	}
	return it.data_hit;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.value:
..cat:Blast
..summary:The hit the iterator points to.
..signature:value(file,it)
..param.it:A hit or hsp iterator on a Blast report with @Spec.StreamReport@.
..param.it:A stream.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:A BlastHit or BlastHsp.
...type:Metafunction.Hit
...type:Metafunction.Hsp
..see:Function.getValue
*/
template<typename TBlastReport, typename TFile>
inline typename Reference<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type
value(TFile & file,
	  Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
	if(it.data_pos != (it.data_hit).begin_pos)
		it.data_hit = getValue(file,it);
	return it.data_hit;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TBlastReport>
inline typename Host<Iter<TBlastReport, StreamBlastIterator<HitIterator> > >::Type &
hostReport(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
	SEQAN_CHECKPOINT
	return *it.data_host;
} 

//////////////////////////////////////////////////////////////////////////////


/**
.Function.atEnd:
..cat:Blast
..summary:Determines whether the iterator is at the end or not.
..signature:atEnd(file,it)
..param.it:A hit or hsp iterator on a Blast report with @Spec.StreamReport@.
..param.it:A stream.
...type:Spec.HitIterator
...type:Spec.HspIterator
..returns:True if the iterator is at the end, false otherwise
..see:Function.goEnd
*/

template<typename TBlastReport, typename TFile>
inline bool
atEnd(TFile &,
	  Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
{
SEQAN_CHECKPOINT
//	return (it.data_last_pos != it.data_pos);	
	return it.data_at_end;	
}

//////////////////////////////////////////////////////////////////////////////

//template<typename TBlastReport>
//inline void
//goEnd(Iter<TBlastReport, StreamBlastIterator<HitIterator> >& it)
//{
//	SEQAN_CHECKPOINT
//	it.data_pos = doof;
//}


//////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
inline void
getNextHitFilePos(TFile & file,
				  Iter<BlastReport<TBlastHsp, StreamReport<TFile> >, StreamBlastIterator<HitIterator> >& it)
{
	typedef typename Position<TFile>::Type TPosition;

	_streamSeekG(file,it.data_pos);
	char c = '>';
	it.data_host->act_c = c;

	_parse_skipWhitespace(file,c);
	_parse_skipLine(file,c);

//	String<char> delim = ">";
	if(_parse_untilBeginLine(file,c,'>'))
	{
		TPosition event_pos = _streamTellG(file);
		if(c=='>' && ( !it.data_host->next_report || (it.data_host->next_report && event_pos < it.data_host->next_report_pos)))
			it.data_next_pos = event_pos;
		else
	        _streamSeekG(file,it.data_host->next_report_pos);

	}//end hit
	else
        _streamSeekG(file,it.data_pos);

}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
