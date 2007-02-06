#ifndef SEQAN_HEADER_ALIGN_BASE_H
#define SEQAN_HEADER_ALIGN_BASE_H

#include <cmath>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment Specific Metafunctions
//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Cols:
..summary:Type of column container of an alignment.
..signature:Cols<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The type of the container that allows access to the columns of $T$.
*/
template <typename T> struct Cols;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Col:
..summary:Type of a column in an alignment.
..signature:Col<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The column type of $T$.
..remarks:The returned type is equivalent to $Value<Cols<T>::Type>::Type$.
..see:Metafunction.Cols
..see:Metafunction.Value
*/
template <typename T>
struct Col:
	Value<typename Cols<T>::Type>
{
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Rows:
..summary:Type of row container of an alignment.
..signature:Rows<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The type of the container that allows access to the rows of $T$.
..see:Metafunction.Cols
*/
template <typename T> struct Rows;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Row:
..summary:Type of a row in an alignment.
..signature:Row<T>::Type
..param.T:An alignment type.
...type:Class.Align
..returns.param.Type:The row type of $T$.
..remarks:The returned type is equivalent to $Value<Rows<T>::Type>::Type$.
..see:Metafunction.Rows
..see:Metafunction.Value
*/
template <typename T>
struct Row:
	Value<typename Rows<T>::Type>
{
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Align
//////////////////////////////////////////////////////////////////////////////
//Default implementation: array of Gaps objects

/**
.Class.Align:
..cat:Alignments
..summary:An alignment of sequences.
..signature:Align<TSource, TSpec>
..param.TSource:Type of the ungapped sequences.
...metafunction:Metafunction.Source
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
..remarks:The default implementation of $Align$ stores the alignment in a set of @Class.Gaps.Gaps<TSource.TSpec>@ objects.
 Hence, the default implementation is row-based, so it will be faster to access the alignment row-wise than column-wise.
*/

template <typename TSource, typename TSpec>
class Align
{
public:
	typedef Gaps<TSource, TSpec> TGaps;
	typedef String<TGaps> TRows;
	typedef typename Size<TRows>::Type TRowsSize;

	TRows data_rows;

//____________________________________________________________________________
public:
	Align()
	{
SEQAN_CHECKPOINT
	}
	Align(Align const & _other):
		data_rows(_other.data_rows)
	{
SEQAN_CHECKPOINT
	}
	~Align()
	{
SEQAN_CHECKPOINT
	}

	Align const &
	operator = (Align const & _other)
	{
SEQAN_CHECKPOINT
		data_rows = _other.data_rows;
		return *this;
	}

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
//ALIGN INTERFACE
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct Value<Align<TSource, TSpec> >:
	Value<Gaps<TSource, TSpec> >
{
};
template <typename TSource, typename TSpec>
struct Value<Align<TSource, TSpec> const>:
	Value<Gaps<TSource, TSpec> const>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.GetValue.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct GetValue<Align<TSource, TSpec> >:
	GetValue<Gaps<TSource, TSpec> >
{
};
template <typename TSource, typename TSpec>
struct GetValue<Align<TSource, TSpec> const>:
	GetValue<Gaps<TSource, TSpec> const>
{
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Reference.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct Reference<Align<TSource, TSpec> >:
	Reference<Gaps<TSource, TSpec> >
{
};
template <typename TSource, typename TSpec>
struct Reference<Align<TSource, TSpec> const>:
	Reference<Gaps<TSource, TSpec> const>
{
};

//////////////////////////////////////////////////////////////////////////////

// struct Cols<Align>: see below (in align_cols_base)

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Rows.param.T.type:Class.Align

template <typename TSource, typename TSpec>
struct Rows<Align<TSource, TSpec> >
{
	typedef String<Gaps<TSource, TSpec> > Type;
};
template <typename TSource, typename TSpec>
struct Rows<Align<TSource, TSpec> const>
{
	typedef String<Gaps<TSource, TSpec> > const Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.rows:
..cat:Alignments
..summary:The container of rows in an alignment. 
..signature:Rows rows(align)
..param.align:An alignment.
...type:Class.Align
..returns:The container of rows in $align$. 
...metafunction:Metafunction.Rows
..see:Function.cols
..see:Metafunction.Rows
*/
template <typename TSource, typename TSpec>
inline typename Rows< Align<TSource, TSpec> >::Type &
rows(Align<TSource, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.data_rows;
}
template <typename TSource, typename TSpec>
inline typename Rows< Align<TSource, TSpec> const >::Type &
rows(Align<TSource, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return me.data_rows;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.row:
..cat:Alignments
..summary:A row in an alignment. 
..signature:Row & row(align, position)
..param.align:An alignment.
...type:Class.Align
..param.position:A position in the @Function.rows@ container of $align$.
..returns:The row in @Function.rows@ container of $align$ at the given $position$. 
...metafunction:Metafunction.Row
..remarks:This function is equivalent to $value(rows(align), position)$.
..see:Function.rows
..see:Function.col
..see:Metafunction.Row
*/
template <typename TSource, typename TSpec, typename TPosition>
inline typename Row< Align<TSource, TSpec> >::Type &
row(Align<TSource, TSpec> & me, 
	TPosition _pos)
{
SEQAN_CHECKPOINT
	return value(rows(me), _pos);
}
template <typename TSource, typename TSpec, typename TPosition>
inline typename Row< Align<TSource, TSpec> const>::Type &
row(Align<TSource, TSpec> const & me, 
	TPosition _pos)
{
SEQAN_CHECKPOINT
	return value(rows(me), _pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.cols:
..cat:Alignments
..summary:The container of columns in an alignment. 
..signature:Cols cols(align)
..param.align:An alignment.
...type:Class.Align
..returns:The container of columns in $align$. 
...metafunction:Metafunction.Cols
..see:Metafunction.Cols
*/
template <typename TSource, typename TSpec>
inline typename Cols< Align<TSource, TSpec> >::Type
cols(Align<TSource, TSpec> & me)
{
SEQAN_CHECKPOINT
	return typename Cols< Align<TSource, TSpec> >::Type(me);
}
template <typename TSource, typename TSpec>
inline typename Cols< Align<TSource, TSpec> const>::Type
cols(Align<TSource, TSpec> const & me)
{
SEQAN_CHECKPOINT
	return typename Cols< Align<TSource, TSpec> const>::Type(me);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.col:
..cat:Alignments
..summary:A column in an alignment. 
..signature:Col & col(align, position)
..param.align:An alignment.
...type:Class.Align
..param.position:A position in the @Function.cols@ container of $align$.
..returns:The column in @Function.cols@ container of $align$ at the given $position$. 
...metafunction:Metafunction.Col
..remarks:This function is equivalent to $value(cols(align), position)$.
..see:Function.cols
..see:Metafunction.Col
*/
template <typename TSource, typename TSpec, typename TPosition>
inline typename Col< Align<TSource, TSpec> >::Type
col(Align<TSource, TSpec> & me,
	TPosition _pos)
{
SEQAN_CHECKPOINT
	return value(cols(me), _pos);
}
template <typename TSource, typename TSpec, typename TPosition>
inline typename Col< Align<TSource, TSpec> const>::Type
col(Align<TSource, TSpec> const & me,
	TPosition _pos)
{
SEQAN_CHECKPOINT
	return value(cols(me), _pos);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.detach.param.object.type:Class.Align

template <typename TSource, typename TSpec>
inline void
detach(Align<TSource, TSpec> & me)
{
SEQAN_CHECKPOINT
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Rows<TAlign>::Type TRows;
	typedef typename Iterator<TRows>::Type TRowsIterator;

	TRowsIterator it = begin(rows(me));
	TRowsIterator it_end = end(rows(me));

	while (it != it_end)
	{
		detach(*it);
		++it;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSource, typename TSpec, typename TIDString>
inline void
write(TFile & target,
	  Align<TSource, TSpec> const & source,
	  TIDString const &,
	  Raw)
{
SEQAN_CHECKPOINT
	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;

	TRowsPosition row_count = length(rows(source));
	TPosition begin_ = beginPosition(cols(source));
	TPosition end_ = endPosition(cols(source));
	
	unsigned int baseCount=0;
	unsigned int leftSpace=6;
	while(begin_ < end_) {
		unsigned int windowSize_ = 50;
		if ((begin_ + windowSize_)>end_) windowSize_=end_ - begin_;

		// Print header line
		unsigned int offset=0;
		if (baseCount != 0) offset = (unsigned int) floor(log((double)baseCount) / log((double)10));
		for(unsigned int j = 0;j<leftSpace-offset;++j) {
			_streamPut(target, ' ');
		}
		_streamPutInt(target, baseCount);
		baseCount+=windowSize_;
		_streamPut(target, ' ');
		for(TPosition i = 1;i<=windowSize_;++i) {
			if ((i % 10)==0) _streamPut(target, ':');
			else if ((i % 5)==0) _streamPut(target, '.');
			else _streamPut(target, ' ');
		}
		_streamPut(target, ' ');
		_streamPut(target, '\n');
		
		// Print sequences
		for(TRowsPosition i=0;i<2*row_count-1;++i) {
			for(unsigned int j = 0;j<leftSpace+2;++j) _streamPut(target, ' ');
			if ((i % 2)==0) {
				TRow& row_ = row(source, i/2);
				typedef typename Iterator<typename Row<TAlign>::Type const>::Type TIter;
				TIter begin1_ = iter(row_, begin_);
				TIter end1_ = iter(row_, begin_ + windowSize_);
				for (; begin1_ != end1_; ++begin1_) {
					if (isGap(begin1_)) _streamPut(target, gapValue<char>());
					else _streamPut(target, *begin1_);
				}
				//_streamWriteRange(target, iter(row_, begin_), iter(row_, begin_ + windowSize_));
			} else {
				for(unsigned int j = 0;j<windowSize_;++j) {
					if ((!isGap(row(source, (i-1)/2), begin_+j)) &&
						(!isGap(row(source, (i+1)/2), begin_+j)) &&
						(row(source, (i-1)/2)[begin_+j]==row(source, (i+1)/2)[begin_+j]))
					{
						_streamPut(target, '|');
					} else {
						_streamPut(target, ' ');
					}
				} 
			}
			_streamPut(target, '\n');
		}
		_streamPut(target, '\n');
		begin_+=50;
	}
	_streamPut(target, '\n');
}

//////////////////////////////////////////////////////////////////////////////
// stream operators
//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TSource, typename TSpec>
inline TStream &
operator << (TStream & target, 
			 Align<TSource, TSpec> const & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

//////////////////////////////////////////////////////////////////////////////
/*
template <typename TStream, typename TSource, typename TSpec>
inline TStream &
operator >> (TStream & source, 
			 Align<TSource, TSpec> & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
*/
//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
