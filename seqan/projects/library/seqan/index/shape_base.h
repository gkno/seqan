#ifndef SEQAN_HEADER_SHAPE_BASE_H
#define SEQAN_HEADER_SHAPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct SimpleShape;
	struct GappedShape;

	template <unsigned q>
	struct FixedGappedShape;


	inline int _intPow(int a, int b)
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
.Class.Shape:
..cat:Index
..summary:Stores a shape for an ungapped or gapped q-gram.
..signature:Shape<TValue, TSpec>
..remarks:The @Metafunction.ValueSize@ of Shape is the ValueSize of TValue which is the alphabet size.
..param.TValue:The type of string the q-gram shape is applied to (e.g. $Dna$).
..param.TSpec:The specializing type.
...default:$GappedShape$, for gapped q-grams.
*/
	template <typename TValue = Dna, typename TSpec = SimpleShape>
	class Shape;

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Shape
	template <typename TValue, typename TSpec>
	struct Value<Shape<TValue,TSpec> >
	{
		typedef unsigned Type;
	};

///.Metafunction.Size.param.T.type:Class.Shape
	template <typename TValue, typename TSpec>
	struct Size<Shape<TValue,TSpec> >
	{
		typedef unsigned Type;
	};

///.Metafunction.ValueSize.param.T.type:Class.Shape
	template <typename TValue, typename TSpec>
	struct ValueSize< Shape<TValue, TSpec> >:
		public ValueSize<TValue> {};


//////////////////////////////////////////////////////////////////////////////

/**
Spec.SimpleShape
..cat:Index
..general:Class.Shape
..summary:For ungapped q-grams.
*/

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with variable length
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue>
	class Shape<TValue, SimpleShape>
	{
	public:
//____________________________________________________________________________

		unsigned					span;
		typename Value<Shape>::Type	hValue;
		typename Value<Shape>::Type	leftFactor;
		TValue						leftChar;
//____________________________________________________________________________
		
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
		Shape() {}

		Shape(unsigned _span)
		{
			resize(*this, _span);
		}

		Shape(Shape const &other):
			span(other.span),
			hValue(other.hValue),
			leftFactor(other.leftFactor),
			leftChar(other.leftChar) {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with fixed length q
	//////////////////////////////////////////////////////////////////////////////

	template <unsigned q>
	struct FixedShape;

	template <typename TValue, unsigned q>
	class Shape<TValue, FixedShape<q> >
	{
	public:
//____________________________________________________________________________

		enum					  { span = q };
		typename Value<Shape>::Type	hValue;
		enum					  { leftFactor = Power<ValueSize<TValue>::VALUE, q - 1>::VALUE };
		TValue						leftChar;
//____________________________________________________________________________
		
		Shape() {}
		Shape(Shape const &other):
			hValue(other.hValue),
			leftChar(other.leftChar) {}
	};



//////////////////////////////////////////////////////////////////////////////

	template <typename TValue, typename TSpec>
	inline typename Value< Shape<TValue, TSpec> >::Type
	value(Shape<TValue, TSpec> &me)
	{
		return me.hValue;
	}

//____________________________________________________________________________

/**.Function.shapeSpan:
..cat:Index
..summary:Span of a q-gram.
..signature:shapeSpan(object)
..param.object.type:Class.Shape
..returns:Span of object.
*/
	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, TSpec> >::Type
	length(Shape<TValue, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		return me.span;
	}

//____________________________________________________________________________

/**.Function.shapeCountBlanks:
..cat:Index
..summary:Number of blanks (irrelevant positions) in a shape.
..signature:shapeCountBlanks(object)
..param.object.type:Class.Shape
..returns:Number of blanks in object.
*/
	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, TSpec> >::Type
	shapeCountBlanks(Shape<TValue, TSpec> &)
	{
	SEQAN_CHECKPOINT
		return 0;
	}

//____________________________________________________________________________

	template <typename TValue, typename TSize>
	inline typename Size< Shape<TValue, SimpleShape> >::Type
	resize(Shape<TValue, SimpleShape> & me, TSize new_length)
	{
	SEQAN_CHECKPOINT
		me.leftFactor = _intPow(ValueSize<TValue>::VALUE, new_length - 1);
		return me.span = new_length;
	}

	template <typename TValue, typename TStringValue, typename TStringSpec>
	inline void
	stringToShape(Shape<TValue, SimpleShape> & me, String<TStringValue, TStringSpec> const & shape_string)
	{
	SEQAN_CHECKPOINT
		resize(me, length(shape_string));
	}

//____________________________________________________________________________

/**.Function.hash:
..cat:Index
..summary:Computes a hash value for a q-gram (alphabetical position value).
..signature:hash(shape,q_gram)
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.q_gram:Sequence iterator pointing to the first position of the q-gram.
..returns:Hash value of the q_gram.
*/

	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash(Shape<TValue, TSpec> &me, TIter it)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;
		typedef typename Size< Shape<TValue, TSpec> >::Type		TSize;

		me.hValue = (THValue)(me.leftChar = *it);
		for(TSize i = 1; i < me.span; ++i) {
			++it;
			me.hValue = me.hValue * ValueSize<Shape<TValue, TSpec> >::VALUE + (THValue)*it;
		}
		return me.hValue;
	}

//____________________________________________________________________________

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

	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hashNext(Shape<TValue, TSpec> &me, TIter &it)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type THValue;
		me.hValue = 
			(me.hValue - (THValue)me.leftChar * me.leftFactor) 
				* ValueSize<Shape<TValue, TSpec> >::VALUE
			+ (THValue)*(it + me.span - 1);
		me.leftChar = *it;
		return me.hValue;
	}

}	// namespace seqan

#endif
