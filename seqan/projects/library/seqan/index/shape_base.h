#ifndef SEQAN_HEADER_SHAPE_BASE_H
#define SEQAN_HEADER_SHAPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct SimpleShape;
	struct GappedShape;
/*
	typedef Tag<_SimpleShape>	SimpleShape;
	typedef Tag<_GappedShape>	GappedShape;
*/

/**
.Class.Shape:
..cat:Index
..summary:Stores a shape for an ungapped or gapped q-gram.
..signature:Shape<TString, TSpec>
..remarks:The @Metafunction.ValueSize@ of Shape is the ValueSize of TString which is the alphabet size.
..param.TString:The type of string the q-gram shape is applied to (e.g. $Dna$).
..param.TSpec:The specializing type.
...default:$GappedShape$, for gapped q-grams.
*/
template <typename TString = Dna, typename TSpec = GappedShape>
class Shape;




inline int intPow(int a, int b)
{
SEQAN_CHECKPOINT
	int ret = 1;
	while (b != 0)
	{
		if (b & 1) ret *= a;
		a *= a;
		b >>= 1;
	}	
	return ret;
}



/**
Spec.SimpleShape
..cat:Index
..general:Class.Shape
..summary:For ungapped q-grams.
*/
template<typename TString>
class Shape<TString, SimpleShape>
{

public:

	int span;
	int /*const*/ term;
	
/**
.Shape#Shape:
..class:Class.Shape
..summary:Constructor
..signature:Shape<TString, TSpec> ()
..signature:Shape<TString, TSpec> (shape)
..signature:Shape<TString, SimpleShape> (span)
..signature:Shape<TString, GappedShape> (num_blanks, span)
..param.shape:Other Shape object. (copy constructor)
..param.span:Total length of the q-gram (including blanks).
..param.num_blanks:Number of blanks (irrelevant positions).
*/
	Shape()
	{
SEQAN_CHECKPOINT
	}

	Shape( int span):
	span(span),
	term(intPow( ValueSize<TString>::VALUE, span-1))
	{
SEQAN_CHECKPOINT
	}

	// Kopierkonstruktor
	Shape(Shape const & rhs):
	span(rhs.span),
	term(rhs.term)
	{
SEQAN_CHECKPOINT
	}	


	~Shape()
	{
SEQAN_CHECKPOINT
	}
	
};



//------------------- Functions -----------------------------
/**.Function.shapeSpan:
..cat:Index
..summary:Span of a q-gram.
..signature:shapeSpan(object)
..param.object.type:Class.Shape
..returns:Span of object.
*/
template <typename TString, typename TSpec>
inline int &
shapeSpan(Shape<TString, TSpec> & me)
{
SEQAN_CHECKPOINT
	return me.span;
}

/**.Function.shapeCountBlanks:
..cat:Index
..summary:Number of blanks (irrelevant positions) in a shape.
..signature:shapeCountBlanks(object)
..param.object.type:Class.Shape
..returns:Number of blanks in object.
*/
template <typename TString>
inline int 
shapeCountBlanks(Shape<TString, SimpleShape>  & me)
{
SEQAN_CHECKPOINT
	return 0;
}

/**.Function.hash:
..cat:Index
..summary:Computes a hash value for a q-gram (alphabetical position value).
..signature:hash(shape,q_gram)
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.q_gram:Sequence iterator pointing to the first position of the q-gram.
..returns:Hash value of the q_gram.
*/
template <typename TString, typename TSequenceIterator>
typename Size<Shape<TString, SimpleShape> >::Type
hash(Shape<TString, SimpleShape> & shape, TSequenceIterator qgram_it)	
{
SEQAN_CHECKPOINT
	typename Size<Shape<TString, SimpleShape> >::Type pos = 0;		
	int i = 0;
	int span = shapeSpan(shape);
	while(i < span)
	{
		pos = pos * ValueSize<Shape<TString, SimpleShape> >::VALUE + (unsigned)*qgram_it;
		++qgram_it;
		++i;
	}


	return pos;
}

/**
.Function.hashNext:
..cat:Index
..summary:Computes the hash value for the next q-gram given the q-gram before (that starts one 
position before the next q-gram) and its hash value.
..signature:hashNext(shape,q_gram1,qgram2,hash1)
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.q_gram1:Sequence iterator pointing to first position of the q-gram that the 
hash value is supposed to be determined for.
..param.q_gram1:Sequence iterator pointing to first position of the last q-gram.
..param.hash1:Hash value of the last q-gram.
..returns:Hash value of q_gram1.
*/
template <typename TSequenceIterator, typename TString>
typename Size<Shape<TString, SimpleShape> >::Type
inline hashNext(
	Shape<TString, SimpleShape> & shape, 
	TSequenceIterator qgram_1_it, 
	TSequenceIterator qgram_2_it, 
	typename Size<Shape<TString, SimpleShape> >::Type x)
{
SEQAN_CHECKPOINT
	int span = shapeSpan(shape);
    return 
		(x - (int)*qgram_1_it * shape.term) * ValueSize<Shape<TString, SimpleShape> >::VALUE
		+ (int)*(qgram_2_it + span - 1);
}



template <typename TString>
void
stringToShape(Shape<TString,SimpleShape> & shape, String<char> const & shape_string)
{
	shape.span = length(shape_string);
	shape.term = intPow(ValueSize<TString>::VALUE, length(shape_string)-1);
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

///.Metafunction.ValueSize.param.T.type:Class.Shape
template <typename TString, typename TSpec>
struct ValueSize< Shape<TString, TSpec> > {
	enum { VALUE = ValueSize<TString>::VALUE }; 
};


}	// namespace seqan

#endif
