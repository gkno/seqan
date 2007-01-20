#ifndef SEQAN_HEADER_SHAPE_BASE_H
#define SEQAN_HEADER_SHAPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

struct SimpleShape;
struct GappedShape1;
struct GappedShape2;
struct GappedShape3;

// Objekt Shape, klassiefiziert die Seedart. Speichert die Anzahl und die genauen Stellen innerhalb eines q-Grams,
// die bei der Seed-Suche irrelevant sind.
// Weiterhin speichert es die Alpabhet-größe, die Größe des q-Grams q und berechnet damit einen für die 
// Berechnung des hash-Wertes des entsprechenden q-Grams benötigten term.
// Ist die Anzhal irrelevanter Stellen = 0, handelt es sich um exakte Seeds, sonst um exakte gegappte seeds.


template <typename TString = Dna, typename TSpec = GappedShape3>
class Shape;

template<typename TString>
class Shape<TString, SimpleShape>
{

public:

	int q;
	int const term;
	

	// Konstruktor
	Shape( int q, int alp_size):
	q(q),
	term(intPow(alp_size, q-1))
	{
	}

	// Kopierkonstruktor
	Shape(Shape const & rhs):
	q(rhs.q),
	term(rhs.term),
	{
	}	


	~Shape(){}
	
};





template <typename TString, typename TSpec>
inline int &
powTerm(Shape<TString, TSpec>  & me)
{
	return me.term;
}

template <typename TString, typename TSpec>
inline int &
seedSize(Shape<TString, TSpec> & me)
{
	return me.q;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hashfunktion, weist einem q-Gram eine Position zu (alphabetisch)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename TString, typename TSequenceIterator>
int 
hash(Shape<TString, SimpleShape> & shape, TSequenceIterator qgram_it)	
{

//	typename Size<Shape<TString,SimpleShape> >::Type const alp_size = ValueSize<Shape<TString, SimpleShape> >::VALUE;

	int pos = 0;			// resultierende Position des q-Grams
	int i = 0;
	int span = seedSize(shape);
	while(i < span)
	{
		pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*qgram_it;
		++qgram_it;
		++i;
	}


	return pos;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hash_next, weist einem q-Gram eine Position zu (alphabetisch) (bei exakten Seeds abgeleitete Berechnung,
// bei exakte gegappte Seeds Aufruf von hash)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename TSequenceIterator, typename TString>
int
hash_next(Shape<TString, SimpleShape> & shape, TSequenceIterator qgram_1_it, TSequenceIterator qgram_2_it, int & x)
{

//	typename Size<Shape<TString,SimpleShape> >::Type const alp_size = ValueSize<Shape<TString, SimpleShape> >::VALUE;
	int term = powTerm(shape);
	int span = seedSize(shape);
	int y;
    y = (x - (int)*qgram_1_it * term) * ValueSize<Shape<TString, SimpleShape> >::VALUE + (int)*(qgram_2_it+span);
	
	return y;
}




//////////// METAFUNCTIONS /////////////////


///.Metafunction.Size.param.T.type:Class.Shape
template <typename TString, typename TSpec>
struct Size<Shape<TString,TSpec> >
{
	typedef int Type;
};

///.Metafunction.Value.param.T.type:Class.Shape
template <typename TString, typename TSpec>
struct Value<Shape<TString,TSpec> >
{
	typedef int Type;
};

template <typename TString, typename TSpec>
struct ValueSize< Shape<TString, TSpec> > {
	enum { VALUE = ValueSize<TString>::VALUE }; 
};


}	// namespace seqan

#endif
