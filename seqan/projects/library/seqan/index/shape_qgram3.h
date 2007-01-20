#ifndef SEQAN_HEADER_SHAPE_QGRAM3_H
#define SEQAN_HEADER_SHAPE_QGRAM3_H

namespace SEQAN_NAMESPACE_MAIN
{

struct GappedShape3;


// Objekt Shape, klassiefiziert die Seedart. Speichert die Anzahl und die genauen Stellen innerhalb eines q-Grams,
// die bei der Seed-Suche irrelevant sind.
// Weiterhin speichert es die Alpabhet-größe, die Größe des q-Grams q und berechnet damit einen für die 
// Berechnung des hash-Wertes des entsprechenden q-Grams benötigten term.
// Ist die Anzhal irrelevanter Stellen = 0, handelt es sich um exakte Seeds, sonst um exakte gegappte seeds.
template<typename TString>
class Shape<TString,GappedShape3>
{

public:

	int q;
	int num_gaps;
	int /*const*/ term;
	//shape_len ist die länge des shape vektors
	int shape_len;
	//[num_springen num_springen...]
	String<int> shape;
	
	Shape()
	{
	}
	// Konstruktor
	Shape(int num_gaps, int q, int alp_size, int shape_len):
	q(q),
	shape_len(shape_len),
	num_gaps(num_gaps),
	term(intPow(alp_size, q-1-num_gaps))
	{
		resize(shape, shape_len);
		for(int i = 0; i < shape_len; ++i)
			shape[i] = 1;
	}

	// Kopierkonstruktor
	Shape(Shape const & rhs):
	q(rhs.q),
	num_gaps(rhs.num_gaps),
	shape_len(rhs.shape_len),
	term(rhs.term),
	shape(rhs.shape)
	{
	}	


	int & operator[] (int offset){ return shape[offset];}



	~Shape(){}
	
};

template <typename TString>
void
setShape(Shape<TString,GappedShape3> & me,int num_gaps, int q, int alp_size, int shape_len)
{
	me.num_gaps = num_gaps;
	me.q = q;
	me.shape_len = shape_len;
	me.term = intPow(alp_size, q-1-num_gaps);
	resize(me.shape, shape_len);
	for (int i = 0; i < shape_len; ++i)
		me.shape[i] = 1;

}



template <typename TString>
inline int &
numGaps(Shape<TString, GappedShape3>  & me)
{
	return me.num_gaps;
}


template <typename TString>
inline int &
shapeLen(Shape<TString, GappedShape3>  & me)
{
	return me.shape_len;
}

template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, GappedShape3> & shape, TSequenceIterator qgram_it)	
{

	typedef typename Size<Shape<TString,GappedShape3> >::Type TValue;
	typedef typename Position<TString>::Type TPosition;
	TPosition pos = 0;			// resultierende Position des q-Grams
	TValue * pshape = begin(shape.shape);
	TValue * pshape_end = end(shape.shape);//
		
	while(pshape < pshape_end)
	{
		pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		qgram_it += *pshape;
		++pshape;
	}

	return pos;

}


template <typename TSequenceIterator, typename TString>
int
hash_next(Shape<TString, GappedShape3> & shape, TSequenceIterator & qgram_1, TSequenceIterator & qgram_2, int & x)
{
	return hash(shape, qgram_2);
}



template <typename TString>
void 
stringToShape(String<char> & shape_string,Shape<TString,GappedShape3> & shape)
{
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
	setShape(shape,count_gaps,length(shape_string),ValueSize<Shape<TString,GappedShape3> >::VALUE,shape_len);
	string_it_last = begin(shape_string);
	string_it = string_it_last +1;
	int j = 0;
	while(string_it < string_end) 
	{
		if(*string_it != '_')
		{
			++j;
		}
		else
		{
			++shape[j];
		}
		++string_it_last;
		++string_it;
	}
	shape[shape_len-1] = 0;

}


}	// namespace seqan

#endif
