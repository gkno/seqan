#ifndef SEQAN_HEADER_SHAPE_QGRAM1_H
#define SEQAN_HEADER_SHAPE_QGRAM1_H

namespace SEQAN_NAMESPACE_MAIN
{

struct GappedShape1;

template<typename TString>
class Shape<TString, GappedShape1>
{

public:

	//q = qgram inkl. lücken
	int q;
	int num_gaps;
	int /*const*/ term;
	String<int> shape;
	
	Shape()
	{
	}

	// Konstruktor
	Shape(int num_gaps, int q, int alp_size):
	q(q),
	num_gaps(num_gaps),
	//q = qgram ohne lücken
	term(intPow(alp_size, q-1-num_gaps))
	{
		resize(shape, num_gaps + 1);
		for (int i = 0; i < num_gaps + 1; ++i)
		shape[i] = -1;

	}

	// Kopierkonstruktor
	Shape(Shape const & rhs):
	q(rhs.q),
	num_gaps(rhs.num_gaps),
	term(rhs.term),
	shape(rhs.shape)
	{
	}	


	~Shape(){}
	
	int & operator[] (int offset){ return shape[offset];}
};


template <typename TString>
void
setShape(Shape<TString,GappedShape1> & me,int num_gaps, int q, int alp_size)
{
	me.num_gaps = num_gaps;
	me.q = q;
	//q = qgram ohne lücken
	me.term = intPow(alp_size, q-1-num_gaps);
	resize(me.shape, num_gaps + 1);
	for (int i = 0; i < num_gaps + 1; ++i)
		me.shape[i] = -1;


}


template <typename TString>
inline int &
numGaps(Shape<TString, GappedShape1>  & me)
{
	return me.num_gaps;
}




template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, GappedShape1> & shape, TSequenceIterator qgram_it)	
{

	typedef typename Size<Shape<TString,GappedShape3> >::Type TValue;
	typedef typename Position<TString>::Type TPosition;
	TPosition pos = 0;			// resultierende Position des q-Grams
	TValue * pshape = begin(shape.shape);
	TValue i = 0;
	TValue span = seedSize(shape);

	while(i < span)
	{
		if(i == *pshape)
			++pshape;
		else
			pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		++i;
		++qgram_it;
	}

	return pos;
}



template <typename TSequenceIterator, typename TString>
int
hash_next(Shape<TString, GappedShape1> & shape, TSequenceIterator qgram_1_it, TSequenceIterator qgram_2_it, int & x)
{
	return hash(shape, qgram_2_it);
}



//returns a Shape object that represents the string given, e.g. "xxx__x__x_xxx"
//'_' means gap
template <typename TString>
void
stringToShape(String<char> & shape_string, Shape<TString,GappedShape1> & shape)
{
	int count_gaps = 0;
	typename Iterator<String<char> >::Type string_it, string_it_last, string_end;
	string_it = begin(shape_string);
	string_end = end(shape_string);
	while(string_it < string_end) 
	{
		if(*string_it == '_')
			++count_gaps;
		++string_it;
	}
	//q = qgram ohne lücken
	setShape(shape,count_gaps,length(shape_string),ValueSize<Shape<TString,GappedShape1> >::VALUE);
	string_it = begin(shape_string);
	int j = 0;
	for(int i = 0; i < length(shape_string); ++i) 
	{
		if(shape_string[i] == '_')
		{
			shape[j] = i;
			++j;
		}
	}
	
}



}	// namespace seqan

#endif
