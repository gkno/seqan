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
template<typename TString>
class Shape<TString,GappedShape>
{

public:

	int span;
	int num_gaps;
	int shape_len;
	String<int> shape;
	
	Shape()
	{
SEQAN_CHECKPOINT
	}

	Shape(int num_gaps, int span):
	span(span),
	shape_len(span - num_gaps),
	num_gaps(num_gaps)
	{
SEQAN_CHECKPOINT
		resize(shape, span - num_gaps);
		for(int i = 0; i < span - num_gaps; ++i)
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


	int & operator[] (int offset)
	{
	SEQAN_CHECKPOINT
		return shape[offset];
	}



	~Shape()
	{
		SEQAN_CHECKPOINT
	}
	
};






template <typename TString>
void
_setShape(Shape<TString,GappedShape> & me, int span, 
		  int num_gaps, int shape_len)
{
SEQAN_CHECKPOINT
	me.num_gaps = num_gaps;
	me.span = span;
	me.shape_len = shape_len;
	resize(me.shape, shape_len);
	for (int i = 0; i < shape_len; ++i)
		me.shape[i] = 1;
	me.shape[0] = 0;

}



template <typename TString>
inline int &
shapeCountBlanks(Shape<TString, GappedShape>  & me)
{
SEQAN_CHECKPOINT
	return me.num_gaps;
}



template <typename TString, typename TSequenceIterator>
typename Size<Shape<TString, GappedShape> >::Type
hash(Shape<TString, GappedShape> & shape, TSequenceIterator qgram_it)	
{
SEQAN_CHECKPOINT

	typedef typename Size<Shape<TString,GappedShape> >::Type TSize;
	typedef typename Value<Shape<TString,GappedShape> >::Type TValue;
	TValue pos = 0;			
	TSize * pshape = begin(shape.shape);
	TSize * pshape_end = end(shape.shape);//
			
	while(pshape < pshape_end)
	{
		qgram_it += *pshape;
		pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		++pshape;
	}

	return pos;

}


template <typename TSequenceIterator, typename TString>
typename Size<Shape<TString, GappedShape> >::Type
hashNext(Shape<TString, GappedShape> & shape, 
		 TSequenceIterator & qgram_1, TSequenceIterator & qgram_2, 
		 typename Size<Shape<TString, GappedShape> >::Type & x)
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
template <typename TString>
void 
stringToShape(Shape<TString,GappedShape> & shape,String<char> & shape_string)
{
SEQAN_CHECKPOINT
	int count_gaps = 0;
	int shape_len;
	typename Iterator<String<char> >::Type string_it, string_it_last, string_end;
	string_it = begin(shape_string);
	string_it_last = string_it;
	++string_it;
	string_end = end(shape_string);
	while(string_it < string_end) 
	{
		if(*string_it == '_')
				++count_gaps;
		++string_it_last;
		++string_it;
	}
	shape_len = length(shape_string) - count_gaps;
	_setShape(shape,length(shape_string),count_gaps,shape_len);
	string_it = begin(shape_string);
	int j = 0;
	while(string_it < string_end && j < shape_len) 
	{
		if(*string_it != '_')
		{
			++j;
		}
		else
		{
			++shape[j];
		}
		++string_it;
	}

}


}	// namespace seqan

#endif
