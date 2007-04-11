#ifndef SEQAN_HEADER_SHAPE_BASE_H
#define SEQAN_HEADER_SHAPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct SimpleShape;
	struct GappedShape;

	template <unsigned q>
	struct FixedShape;

	template <unsigned q>
	struct FixedGappedShape;

/*
	typedef Tag<_SimpleShape>	SimpleShape;
	typedef Tag<_GappedShape>	GappedShape;
*/

/**
.Class.Shape:
..cat:Index
..summary:Stores a shape for an ungapped or gapped q-gram.
..signature:Shape<TValue, TSpec>
..remarks:The @Metafunction.ValueSize@ of Shape is the ValueSize of TValue which is the alphabet size.
..param.TValue:The type of string the q-gram shape is applied to (e.g. $Dna$).
..param.TSpec:The specializing type.
...default:$GappedShape$, for gapped q-grams.
*/
template <typename TValue = Dna, typename TSpec = GappedShape>
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
template<typename TValue>
class Shape<TValue, SimpleShape>
{

public:

	unsigned span;
	int /*const*/ term;
	
/**
.Shape#Shape:
..class:Class.Shape
..summary:Constructor
..signature:Shape<TValue, TSpec> ()
..signature:Shape<TValue, TSpec> (shape)
..signature:Shape<TValue, SimpleShape> (span)
..signature:Shape<TValue, GappedShape> (num_blanks, span)
..param.shape:Other Shape object. (copy constructor)
..param.span:Total length of the q-gram (including blanks).
..param.num_blanks:Number of blanks (irrelevant positions).
*/
	Shape()
	{
SEQAN_CHECKPOINT
	}

	Shape( unsigned span):
	span(span),
	term(intPow( ValueSize<TValue>::VALUE, span-1))
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
template <typename TValue, typename TSpec>
inline unsigned &
shapeSpan(Shape<TValue, TSpec> & me)
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
template <typename TValue, typename TSpec>
inline unsigned 
shapeCountBlanks(Shape<TValue, TSpec> &)
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
template <typename TValue, typename TIter>
typename Size<Shape<TValue, SimpleShape> >::Type
hash(Shape<TValue, SimpleShape> & shape, TIter qgram_it)	
{
SEQAN_CHECKPOINT
	typename Size<Shape<TValue, SimpleShape> >::Type pos = 0;		
	unsigned i = 0;
	unsigned span = shapeSpan(shape);
	while(i < span)	{
		pos = pos * ValueSize<Shape<TValue, SimpleShape> >::VALUE + (unsigned)*qgram_it;
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
template <typename TIter, typename TValue>
typename Size<Shape<TValue, SimpleShape> >::Type
inline hashNext(
	Shape<TValue, SimpleShape> & shape, 
	TIter it1, 
	TIter it2, 
	typename Size<Shape<TValue, SimpleShape> >::Type x)
{
SEQAN_CHECKPOINT
	unsigned span = shapeSpan(shape);
    return 
		(x - (int)*it1 * shape.term) * ValueSize<Shape<TValue, SimpleShape> >::VALUE
		+ (int)*(it2 + span - 1);
}



template <typename TValue>
void
stringToShape(Shape<TValue,SimpleShape> & shape, String<char> const & shape_string)
{
	shape.span = length(shape_string);
	shape.term = intPow(ValueSize<TValue>::VALUE, length(shape_string)-1);
}





	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with fixed length q
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue, unsigned q>
	class Shape<TValue, FixedShape<q> >
	{
	public:
		unsigned long term;
		
		Shape() {}
		Shape(Shape const &other):
			term(other.term) {}
	};

	template <typename TValue, unsigned q>
	inline unsigned
	shapeSpan(Shape<TValue, FixedShape<q> > & me)
	{
	SEQAN_CHECKPOINT
		return q;
	}

	template <typename TValue, unsigned q, typename TIter>
	inline typename Value< Shape<TValue, FixedShape<q> > >::Type
	hash(Shape<TValue, FixedShape<q> > &shape, TIter it)	
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, FixedShape<q> > >::Type THValue;
		THValue hValue = 0;
		for(unsigned i = 0; i < q; ++i, ++it)
			hValue = hValue * ValueSize<Shape<TValue, FixedShape<q> > >::VALUE + (THValue)*it;
		return hValue;
	}

	template <typename TValue, unsigned q, typename TIter>
	inline typename Value< Shape<TValue, FixedShape<q> > >::Type
	hashNext(Shape<TValue, FixedShape<q> > &shape, 	TIter it1, TIter it2, 
		typename Value< Shape<TValue, FixedShape<q> > >::Type x)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, FixedShape<q> > >::Type THValue;
		return 
			(x - (THValue)*it1 * shape.term) * ValueSize<Shape<TValue, FixedShape<q> > >::VALUE
			+ (THValue)*(it2 + q - 1);
	}




//////////// METAFUNCTIONS /////////////////


///.Metafunction.Value.param.T.type:Class.Shape
template <typename TValue, typename TSpec>
struct Value<Shape<TValue,TSpec> >
{
	typedef unsigned Type;
};

///.Metafunction.ValueSize.param.T.type:Class.Shape
template <typename TValue, typename TSpec>
struct ValueSize< Shape<TValue, TSpec> >:
	public ValueSize<TValue> {};


}	// namespace seqan

#endif
