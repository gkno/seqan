#ifndef SEQAN_HEADER_SHAPE_QGRAM3_H
#define SEQAN_HEADER_SHAPE_QGRAM3_H

namespace SEQAN_NAMESPACE_MAIN
{

struct GappedShape;

/**
Spec.GappedShape
..cat:Index
..general:Class.Shape
..summary:For gapped q-grams.
*/
template<typename TValue>
class Shape<TValue,GappedShape>
{

public:

	unsigned span;
	unsigned num_gaps;
	unsigned shape_len;
	String<unsigned> shape;
	
	Shape()
	{
SEQAN_CHECKPOINT
	}

	Shape(unsigned num_gaps, unsigned span):
	span(span),
	num_gaps(num_gaps),
	shape_len(span - num_gaps)
	{
SEQAN_CHECKPOINT
		resize(shape, span - num_gaps);
		for(unsigned i = 0; i < span - num_gaps; ++i)
			shape[i] = 1;
	}

	Shape(Shape const & rhs):
	span(rhs.span),
	num_gaps(rhs.num_gaps),
	shape_len(rhs.shape_len),
	shape(rhs.shape)
	{
SEQAN_CHECKPOINT
	}	


	unsigned & operator[] (unsigned offset)
	{
	SEQAN_CHECKPOINT
		return shape[offset];
	}



	~Shape()
	{
		SEQAN_CHECKPOINT
	}
	
};






template <typename TValue>
void
_setShape(Shape<TValue,GappedShape> & me, unsigned span, 
		  unsigned num_gaps, unsigned shape_len)
{
SEQAN_CHECKPOINT
	me.num_gaps = num_gaps;
	me.span = span;
	me.shape_len = shape_len;
	resize(me.shape, shape_len);
	for (unsigned i = 0; i < shape_len; ++i)
		me.shape[i] = 1;
	me.shape[0] = 0;

}



template <typename TValue>
inline unsigned &
shapeCountBlanks(Shape<TValue, GappedShape> & me)
{
SEQAN_CHECKPOINT
	return me.num_gaps;
}



template <typename TValue, typename TSequenceIterator>
typename Size<Shape<TValue, GappedShape> >::Type
inline hash(Shape<TValue, GappedShape> & shape, TSequenceIterator qgram_it)	
{
SEQAN_CHECKPOINT

	typedef typename Size<Shape<TValue,GappedShape> >::Type TSize;
	typedef typename Value<Shape<TValue,GappedShape> >::Type TValue;
	TValue pos = 0;			
	TSize * pshape = begin(shape.shape);
	TSize * pshape_end = end(shape.shape);//
			
	while(pshape < pshape_end)
	{
		qgram_it += *pshape;
		pos = pos * ValueSize<Shape<TValue, SimpleShape> >::VALUE + (int)*qgram_it;
		++pshape;
	}

	return pos;

}


template <typename TSequenceIterator, typename TValue>
typename Size<Shape<TValue, GappedShape> >::Type
inline hashNext(
	Shape<TValue, GappedShape> & shape,
	TSequenceIterator qgram_1, 
	TSequenceIterator qgram_2, 
	typename Size<Shape<TValue, GappedShape> >::Type x)
{
SEQAN_CHECKPOINT
	return hash(shape, qgram_2);
}



/**.Function.stringToShape:
..cat:Index
..summary:Takes a shape given as a string of 'x' (relevant position) and '_' 
(irrelevant position) and converts it into a Shape object.
..signature:stringToShape(object,shape_string)
..param.object:Shape object that is manipulated.
...type:Class.Shape
..param.shape_string:A string of 'x' and '_' representing relevant and irrelevant positions (blanks) respectively.
...type:Class.String
*/
template <typename TValue>
void 
stringToShape(Shape<TValue,GappedShape> & shape,String<char> const & shape_string)
{
SEQAN_CHECKPOINT
	unsigned count_gaps = 0;
	unsigned shape_len;
	typename Iterator<String<char> const>::Type string_it, string_it_last, string_end;
	string_it = begin(shape_string);
	string_it_last = string_it;
	++string_it;
	string_end = end(shape_string);
	while(string_it < string_end) 
	{
		if (*string_it == '_')
			++count_gaps;
		++string_it_last;
		++string_it;
	}
	shape_len = length(shape_string) - count_gaps;
	_setShape(shape,length(shape_string),count_gaps,shape_len);
	string_it = begin(shape_string);
	unsigned j = 0;
	while(string_it < string_end && j < shape_len) 
	{
		if(*string_it != '_')
			++j;
		else
			++shape[j];
		++string_it;
	}
}


}	// namespace seqan

#endif
