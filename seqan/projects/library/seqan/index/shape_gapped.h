 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_SHAPE_GAPPED_H
#define SEQAN_HEADER_SHAPE_GAPPED_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// HardwiredShape allows compiler-time defined gapped shape

	// Pxx = spaces between '1's
	template <
		int P00 = 0, int P01 = 0, int P02 = 0, int P03 = 0, int P04 = 0, 
		int P05 = 0, int P06 = 0, int P07 = 0, int P08 = 0, int P09 = 0,
		int P10 = 0, int P11 = 0, int P12 = 0, int P13 = 0, int P14 = 0,
		int P15 = 0, int P16 = 0, int P17 = 0, int P18 = 0, int P19 = 0	
	>
	struct HardwiredShape {
		static const int DIFFS[];
	};

	template <
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	const int HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	>::DIFFS[] = {
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	};


//////////////////////////////////////////////////////////////////////////////
// Length meta-function for fixed gapped shapes

	template <>
	struct Length< HardwiredShape<
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0> >
	{
		enum { VALUE = 1 };
	};

	template <
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	struct Length< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19> >
	{
		enum { VALUE = Length< HardwiredShape<
			P01,P02,P03,P04,P05,
			P06,P07,P08,P09,P10,
			P11,P12,P13,P14,P15,
			P16,P17,P18,P19, 0 > >::VALUE + P00 };
	};

	template <
	    typename TValue,
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	struct Length< Shape<TValue, FixedGappedShape< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > > >:
	Length< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > {};


//////////////////////////////////////////////////////////////////////////////
// Weight meta-function for fixed gapped shapes

	template <>
	struct Weight< HardwiredShape<
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0> >
	{
		enum { VALUE = 1 };
	};

	template <
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	struct Weight< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19> >
	{
		enum { VALUE = Weight< HardwiredShape<
			P01,P02,P03,P04,P05,
			P06,P07,P08,P09,P10,
			P11,P12,P13,P14,P15,
			P16,P17,P18,P19, 0 > >::VALUE + 1 };
	};

	template <
	    typename TValue,
		int P00, int P01, int P02, int P03, int P04, 
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19	
	>
	struct Weight< Shape<TValue, FixedGappedShape< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > > >:
	Weight< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > {};


//////////////////////////////////////////////////////////////////////////////

/**
Spec.GappedShape
..cat:Index
..general:Class.Shape
..summary:For gapped q-grams.
*/

	//////////////////////////////////////////////////////////////////////////////
	// variable gapped shape
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue>
	class Shape<TValue, GappedShape>
	{
	public:
//____________________________________________________________________________

		unsigned span;
		unsigned weight;
		String<int> diffs;
	
		typename Value<Shape>::Type	hValue;		// current hash value

		Shape()	{}

		// c'tor for ungapped shapes
		Shape(unsigned _span):
			span(_span),
			weight(_span)
		{
		SEQAN_CHECKPOINT
			resize(diffs, _span);
			for(unsigned i = 0; i < _span; ++i)
				diffs[i] = 1;
		}

		Shape(Shape const &other):
			span(other.span),
			weight(other.weight),
			diffs(other.diffs),
			hValue(other.hValue) {}	
	};

	//////////////////////////////////////////////////////////////////////////////
	// fixed gapped shape
	//////////////////////////////////////////////////////////////////////////////

	template <typename TValue, typename TSpec>
	class Shape<TValue, FixedGappedShape<TSpec> >
	{
	public:
//____________________________________________________________________________

		typedef FixedGappedShape<TSpec>	TShapeSpec;

		enum { span = Length<Shape>::VALUE };
		enum { weight = Weight<Shape>::VALUE };
		const int *diffs;
	
		typename Value<Shape>::Type	hValue;		// current hash value

		Shape():
			diffs(TSpec::DIFFS) {}

		Shape(Shape const &other):
			diffs(other.diffs),	
			hValue(other.hValue) {}
	};

//////////////////////////////////////////////////////////////////////////////

	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	shapeCountBlanks(Shape<TValue, FixedGappedShape<TSpec> > const & me)
	{
	SEQAN_CHECKPOINT
		return me.span - me.weight;
	}

	template <typename TValue, typename TSpec>
	inline typename Size< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	shapeWeight(Shape<TValue, FixedGappedShape<TSpec> > const & me)
	{
	SEQAN_CHECKPOINT
		return me.weight;
	}
//____________________________________________________________________________

	template <typename TValue, typename TIter>
	inline typename Value< Shape<TValue, GappedShape> >::Type
	hash(Shape<TValue, GappedShape> &me, TIter it)	
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, GappedShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, GappedShape> >::Type	TSize;

		me.hValue = _ord(*it);
		for(TSize i = 1; i < me.span; ++i) {
			goFurther(it, me.diffs[i]);
			me.hValue = me.hValue * ValueSize<TValue>::VALUE + _ord(*it);
		}
		return me.hValue;
	}

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	hash(Shape<TValue, FixedGappedShape<TSpec> > &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type	THValue;

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = _ord(*it);
			for(TSize i = 1; i < iEnd; ++i) {
				goFurther(it, me.diffs[i]);
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + _ord(*it);
			}
		} else
			return me.hValue = 0;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

	template <typename TValue, typename TSpec, typename TIter, typename TSize>
	inline typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	hashUpper(Shape<TValue, FixedGappedShape<TSpec> > &me, TIter it, TSize charsLeft)
	{
	SEQAN_CHECKPOINT
		typedef typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type	THValue;

		TSize iEnd = me.span;
		if (iEnd > charsLeft) iEnd = charsLeft;

		TSize i = 0;
		if (iEnd > 0) {
			me.hValue = _ord(*it);
			for(TSize i = 1; i < iEnd; ++i) {
				goFurther(it, me.diffs[i]);
				me.hValue = me.hValue * ValueSize<TValue>::VALUE + _ord(*it);
			}
			++me.hValue;
		} else
			return me.hValue = 1;

		// fill shape with zeros
		for(; i < (TSize)me.span; ++i)
			me.hValue *= ValueSize<TValue>::VALUE;
		return me.hValue;
	}

//____________________________________________________________________________

	template <typename TValue, typename TSpec, typename TIter>
	inline void
	_hashHardwiredShape(Shape<TValue, TSpec> &, TIter &, HardwiredShape<
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0 > const)
	{
	}

	template <
		         int P01, int P02, int P03, int P04,
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19,
		typename TValue, typename TSpec, typename TIter
	>
	inline void
	_hashHardwiredShape(Shape<TValue, TSpec> &me, TIter &it, HardwiredShape<
		 1 ,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 > const)
	{
		typedef typename Value<Shape<TValue, TSpec> >::Type	THValue;

		++it;
		me.hValue = me.hValue * ValueSize<TValue>::VALUE + _ord(*it);

		_hashHardwiredShape(me, it, HardwiredShape<
			P01,P02,P03,P04,P05,
			P06,P07,P08,P09,P10,
			P11,P12,P13,P14,P15,
			P16,P17,P18,P19, 0 >());
	}

	template <
		int P00, int P01, int P02, int P03, int P04,
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19,
		typename TValue, typename TSpec, typename TIter
	>
	inline void
	_hashHardwiredShape(Shape<TValue, TSpec> &me, TIter &it, HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 > const)
	{
		typedef typename Value<Shape<TValue, TSpec> >::Type	THValue;

		it += P00;
		me.hValue = me.hValue * ValueSize<TValue>::VALUE + _ord(*it);

		_hashHardwiredShape(me, it, HardwiredShape<
			P01,P02,P03,P04,P05,
			P06,P07,P08,P09,P10,
			P11,P12,P13,P14,P15,
			P16,P17,P18,P19, 0 >());
	}

	template <
		int P00, int P01, int P02, int P03, int P04,
		int P05, int P06, int P07, int P08, int P09,
		int P10, int P11, int P12, int P13, int P14,
		int P15, int P16, int P17, int P18, int P19,
		typename TValue, typename TIter
	>
	typename Value< Shape<TValue, FixedGappedShape< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > > >::Type
	hash(Shape<TValue, FixedGappedShape< HardwiredShape<
		P00,P01,P02,P03,P04,
		P05,P06,P07,P08,P09,
		P10,P11,P12,P13,P14,
		P15,P16,P17,P18,P19 
	> > > &me, TIter it)
	{
	SEQAN_CHECKPOINT
		typedef HardwiredShape<
			P00,P01,P02,P03,P04,
			P05,P06,P07,P08,P09,
			P10,P11,P12,P13,P14,
			P15,P16,P17,P18,P19 >								TSpec;
		typedef FixedGappedShape<TSpec>							TShape;
		typedef typename Value< Shape<TValue, TShape> >::Type	THValue;
		typedef typename Size< Shape<TValue, TShape> >::Type	TSize;

		me.hValue = _ord(*it);
		_hashHardwiredShape(me, it, TSpec());
		return me.hValue;
	}

//____________________________________________________________________________

	template <typename TValue, typename TSpec, typename TIter>
	inline typename Value< Shape<TValue, FixedGappedShape<TSpec> > >::Type
	hashNext(Shape<TValue, FixedGappedShape<TSpec> > &me, TIter it)	
	{
	SEQAN_CHECKPOINT
		return hash(me, it);
	}


//____________________________________________________________________________

/**.Function.stringToShape:
..cat:Index
..summary:Takes a shape given as a string of '1' (relevant position) and '0' 
(irrelevant position) and converts it into a Shape object.
..signature:stringToShape(shape, shapeString)
..param.shape:Shape object that is manipulated.
...type:Class.Shape
..param.shapeString:A string of '1' and '0' representing relevant and irrelevant positions (blanks) respectively. This string must begin with a '1'.
...type:Class.String
*/

	template <typename TValue, typename TSpec, typename TShapeString>
	inline void
	stringToShape(
		Shape<TValue, FixedGappedShape<TSpec> > &me, 
		TShapeString const &shapeString)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TShapeString const>::Type		TIter;
		typedef typename Iterator<String<int> >::Type			TShapeIter;

		me.span = length(shapeString);

		unsigned oneCount = 0;
		TIter it = begin(shapeString, Standard());
		TIter itEnd = end(shapeString, Standard());
		for(; it != itEnd; ++it)
			if (*it == '1')
				++oneCount;

		me.weight = oneCount;
		resize(me.diffs, oneCount);

		unsigned diff = 0;
		it = begin(shapeString, Standard());
		TShapeIter itS = begin(me.diffs, Standard());
		for(; it != itEnd; ++it) {
			if (*it == '1') {
				*itS = diff;
				++itS;
				diff = 0;
			}
			++diff;
		}
	}

}	// namespace seqan

#endif
