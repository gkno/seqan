#ifndef SEQAN_HEADER_STREAM_ALGORITHMS_H
#define SEQAN_HEADER_STREAM_ALGORITHMS_H


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////


/**
.Internal._streamPutInt:
..summary:Converts an integer to a character and writes it to stream.
..cat:Streams
..signature:_streamPutInt(stream, integer)
..param.target:An output stream.
...type:Adaption."std::iostream"
..param.number:A number that is written to $stream$.
*/
template <typename TStream>
inline void
_streamPutInt(TStream & target,
			  int number)
{
SEQAN_CHECKPOINT
	char str[BitsPerValue<int>::VALUE];
	sprintf(str, "%d", number);
	_streamWrite(target, str);
}


//////////////////////////////////////////////////////////////////////////////



/**
.Internal._streamWrite:
..summary:Writes a sequence to stream.
..cat:Streams
..signature:_streamWrite(stream, sequence)
..param.stream:An input stream.
..param.sequence:A sequence that is written to $stream$.
*/

template <typename TTarget, typename TSource>
inline void
_streamWrite(TTarget & target,
			 TSource const & source)
{
SEQAN_CHECKPOINT
	typename Iterator<TSource const, Standard>::Type it = begin(source, Standard());
	typename Iterator<TSource const, Standard>::Type it_end = end(source, Standard());

	for (; it < it_end; ++it)
	{
		typename GetValue<TSource const>::Type val_ = getValue(it);
		_streamPut(target, val_);
	}
}

//____________________________________________________________________________

template <typename TTarget, typename TSourceValue>
inline void
_streamWrite(TTarget & target,
			 TSourceValue const * source)
{
SEQAN_CHECKPOINT

	for (; !atEnd(source); ++source)
	{
		_streamPut(target, *source);
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Internal._streamWriteRange:
..summary:Writes a range to stream.
..cat:Streams
..signature:_streamWriteRange(stream, begin_iterator, end_iterator)
..param.stream:An input stream.
..param.sequence:A sequence that is written to $stream$.
*/

template <typename TTarget, typename TIterator>
inline void
_streamWriteRange(TTarget & target,
				  TIterator begin_,
				  TIterator end_)
{
SEQAN_CHECKPOINT
	for (; begin_ != end_; ++begin_)
	{
		_streamPut(target, *begin_);
	}
}


	

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
