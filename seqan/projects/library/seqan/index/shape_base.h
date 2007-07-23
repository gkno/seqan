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


	template <typename TValue, typename TExponent>
	inline TValue _intPow(TValue a, TExponent b)
	{
	SEQAN_CHECKPOINT
		TValue ret = 1;
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

///.Metafunction.Length.param.T.type:Class.Shape
    template <typename TValue, unsigned q>
	struct Length< Shape<TValue, FixedShape<q> > >
	{
		enum { VALUE = q };
	};
    template <typename TValue, unsigned q>
	struct Length< Shape<TValue, FixedGappedShape<q> > >
	{
		enum { VALUE = q };
	};

///.Metafunction.ValueSize.param.T.type:Class.Shape
	template <typename TValue, typename TSpec>
	struct ValueSize< Shape<TValue, TSpec> > {
		enum { VALUE = Power<
						ValueSize<TValue>::VALUE, 
						Length< Shape<TValue, TSpec> >::VALUE >::VALUE };
	};


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
		typename Value<Shape>::Type	XValue;
		typename Value<Shape>::Type	leftFactor;
		typename Value<Shape>::Type	leftFactor2;
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
			XValue(other.XValue),
			leftFactor(other.leftFactor),
			leftFactor2(other.leftFactor2),
			leftChar(other.leftChar) {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// ungapped shape with fixed length q
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue, unsigned q>
	class Shape<TValue, FixedShape<q> >
	{
	public:
//____________________________________________________________________________

		enum { span = q };
		enum { leftFactor = Power<ValueSize<TValue>::VALUE, q - 1>::VALUE };
		enum { leftFactor2 = (Power<ValueSize<TValue>::VALUE, q>::VALUE - 1) / (ValueSize<TValue>::VALUE - 1) };
		// Sigma^(q-1) + Sigma^(q-2) + ... + Sigma + 1

		typename Value<Shape>::Type	hValue;		// current hash value
		typename Value<Shape>::Type	XValue;		// Sum_{i=0..q-1} (x_i + 1)
		TValue						leftChar;	// left-most character
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
	length(Shape<TValue, TSpec> const &me)
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
	shapeCountBlanks(Shape<TValue, TSpec> const &)
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
		me.leftFactor = _intPow((unsigned)ValueSize<TValue>::VALUE, new_length - 1);
		me.leftFactor2 = (_intPow((unsigned)ValueSize<TValue>::VALUE, new_length) - 1) / (ValueSize<TValue>::VALUE - 1);
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
..summary:Computes a hash value for a q-gram.
..signature:hash(shape, it)
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the q-gram.
..returns:Hash value of the q-gram.
*/

	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash(Shape<TValue, TSpec> &me, TIter it)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;
		typedef typename Size< Shape<TValue, TSpec> >::Type		TSize;

		me.hValue = _ord(me.leftChar = *it);
		for(TSize i = 1; i < me.span; ++i) {
			++it;
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + _ord(*it);
		}
		return me.hValue;
	}

//____________________________________________________________________________

/**
.Function.hashNext:
..cat:Index
..summary:Computes the hash value for the next q-gram.
..signature:hashNext(shape, it)
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the next q-gram.
..returns:Hash value of the q-gram.
..remarks:@Function.hash@ has to be called before with $shape$ on the left adjacent q-gram.
*/

	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hashNext(Shape<TValue, TSpec> &me, TIter &it)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;
		me.hValue = 
			(me.hValue - _ord(me.leftChar) * me.leftFactor) * ValueSize<TValue>::VALUE
			+ _ord(*(it + me.span - 1));
		me.leftChar = *it;
		return me.hValue;
	}

//____________________________________________________________________________

/**.Function.hash2:
..cat:Index
..summary:Computes a hash value for a q-gram, even if it is shorter than q.
..signature:hash2(shape, it, charsLeft)
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the q-gram.
..param.shape:Number of characters left to read before the next border. The characters after a border are considered as 0's.
..returns:Hash value of the q-gram.
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = me.XValue = _ord(me.leftChar = *it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				// update sum of x_i
				me.XValue += _ord(*it);
				// shift hash
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
			}
		} else
			me.hValue = me.XValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + me.XValue;
		return me.hValue += iEnd;
	}

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2Upper(Shape<TValue, TSpec> &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		THValue hValue, XValue;
		TSize i = 0;
		if (iEnd > 0) {
			hValue = XValue = _ord(*it);
			for(i = 1; i < iEnd; ++i) {
				++it;
				// update sum of x_i
				XValue += _ord(*it);
				// shift hash
				hValue = hValue * ValueSize<TValue>::VALUE + XValue;
			}
		} else
			hValue = XValue = 0;

		if (charsLeft <= me.span) {
			++XValue;
			++hValue;
		}

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			hValue = hValue * ValueSize<TValue>::VALUE + XValue;
		return hValue += iEnd;
	}

//____________________________________________________________________________

/**
.Function.hash2Next:
..cat:Index
..summary:Computes the hash value for the next q-gram, even if it is shorter than q.
..signature:hash2Next(shape, it)
..param.shape:Shape to be used for hashing the q-gram.
...type:Class.Shape
..param.it:Sequence iterator pointing to the first character of the next q-gram.
..returns:Hash value of the q-gram.
..remarks:@Function.hash@ has to be called before with $shape$ on the left adjacent q-gram.
*/

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, TSpec> >::Type
	hash2Next(Shape<TValue, TSpec> &me, TIter &it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		// remove first, shift left, and add next character
		typedef typename Value< Shape<TValue, TSpec> >::Type	THValue;

		if (charsLeft >= me.span) {
			// update sum of x_i
			me.XValue = me.XValue + _ord(*(it + me.span - 1)) - _ord(me.leftChar);
			// shift hash
			me.hValue = (me.hValue - _ord(me.leftChar) * me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
						- me.span * (ValueSize<TValue>::VALUE - 1);
		} else {
			// update sum of x_i
			me.XValue -= _ord(me.leftChar);
			// shift hash
			me.hValue = (me.hValue - _ord(me.leftChar) * me.leftFactor2) * ValueSize<TValue>::VALUE + me.XValue
				        - charsLeft * (ValueSize<TValue>::VALUE - 1) - ValueSize<TValue>::VALUE;
		}

		me.leftChar = *it;
		return me.hValue;
	}

}	// namespace seqan

#endif
