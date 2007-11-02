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
// Weiterhin speichert es die Alpabhet-groesse, die Groesse des q-Grams q und berechnet damit einen fuer die 
// Berechnung des hash-Wertes des entsprechenden q-Grams benoetigten term.
// Ist die Anzhal irrelevanter Stellen = 0, handelt es sich um exakte Seeds, sonst um exakte gegappte seeds.


template <typename TString = Dna, typename TSpec = GappedShape3>
class Shape;

template<typename TString>
class Shape<TString, SimpleShape>
{

public:

	int q;
	int term;
	

	// Konstruktor
	Shape( int q, int alp_size):
	q(q),
	term(_intPow(alp_size, q-1))
	{
	}

	// Kopierkonstruktor
	Shape(Shape const & rhs):
	q(rhs.q),
	term(rhs.term)
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
struct ValueSize< Shape<TString, TSpec> > 
{
	enum { VALUE = ValueSize<TString>::VALUE }; 
};

///hinzugefuegt
template<typename TShape>
struct ShapeType;

template <typename TString, typename TSpec>
struct ShapeType<Shape<TString,TSpec> >
{
	typedef TString Type;
};
///

}	// namespace seqan

#endif
