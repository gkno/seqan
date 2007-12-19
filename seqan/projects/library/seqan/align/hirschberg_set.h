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

#ifndef SEQAN_HEADER_HIRSCHBERG_SET_H
#define SEQAN_HEADER_HIRSCHBERG_SET_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Internal._HirschbergSet:
..cat:Classes
..summary:Represents a part of the Dynamic-Programming Matrix, defined by the start and end position in the two sequences
..signature:_HirschbergSet()
..signature:_HirschbergSet(begin1,end1,begin2,end2)
..remarks:This class is used to define the part of the sequences, that need to be inspected during in the alignment. It is
used of @Tag.MyersHirschberg@ and @Tag.Hirschberg@ 
..memfunc:Internal._begin1
..memfunc:Internal._begin2
..memfunc:Internal._setBegin1
..memfunc:Internal._setBegin2
..memfunc:Internal._end1
..memfunc:Internal._end2
..memfunc:Internal._setEnd1
..memfunc:Internal._setEnd2
..memfunc:Internal._score
..memfunc:Internal._setScore
*/
class _HirschbergSet
{

private:
	int x1,x2,y1,y2;
	int score;

public:
	_HirschbergSet()
		: x1(0),x2(0),y1(0),y2(0)
	{
SEQAN_CHECKPOINT
	}

	_HirschbergSet(int a1,int a2,int b1,int b2,int sc)
		: x1(a1),x2(a2),y1(b1),y2(b2),score(sc)
	{
SEQAN_CHECKPOINT
		SEQAN_ASSERT(a1 <= b1);
		SEQAN_ASSERT(a2 <= b2);
	}

	_HirschbergSet & 
	operator = (_HirschbergSet const & other_)
	{
SEQAN_CHECKPOINT
		_setBegin1(*this,_begin1(other_));
		_setEnd1(*this,_end1(other_));
		_setBegin2(*this,_begin2(other_));
		_setEnd2(*this,_end2(other_));
		_setScore(*this,_score(other_));
		return *this;
	}

	// //////////////////////////////////////////////////////////////////////////////////////////////
	// Accessor Methods
	// //////////////////////////////////////////////////////////////////////////////////////////////

	// Sequence 1
/**
.Internal.begin1:
..cat:Functions
..summary:Returns the begin position of the _HirschberSet in the first sequence
..signature:begin1(_HirschbergSet)
..param._HirschbergSet: Reference to the _HirschbergSet
..param._HirschbergSet.type:Internal._HirschbergSet
*/
/**
.Internal.setBegin1:
..cat:Functions
..summary:Sets the begin position of the _HirschberSet in the first sequence
..signature:setBegin1(_HirschbergSet,new_begin)
..param._HirschbergSet: Reference to the _HirschbergSet
..param.new_begin: The new begin position in the first sequence
*/

///.Internal.setBegin1.param._HirschbergSet.type:Internal._HirschbergSet

/**
.Internal.end1:
..cat:Functions
..summary:Returns the end position of the _HirschberSet in the first sequence
..signature:end1(_HirschbergSet)
..param._HirschbergSet: Reference to the _HirschbergSet
..param._HirschbergSet.type:Internal._HirschbergSet
*/
/**
.Internal.setEnd1:
..cat:Functions
..summary:Sets the end position of the _HirschberSet in the first sequence
..signature:setEnd1(_HirschbergSet,new_end)
..param._HirschbergSet: Reference to the _HirschbergSet
..param._HirschbergSet.type:Internal._HirschbergSet
..param.new_begin: The new end position in the first sequence
*/
	friend inline int &
	_begin1(_HirschbergSet & me)
	{
SEQAN_CHECKPOINT
		return me.x1;
	}


	friend inline int const & 
	_begin1(_HirschbergSet const & me)
	{
SEQAN_CHECKPOINT
		return me.x1;
	}

	friend inline void
	_setBegin1(_HirschbergSet & me, int const & new_begin)
	{
SEQAN_CHECKPOINT
		me.x1 = new_begin;
	}

	friend inline int &
	_end1(_HirschbergSet & me)
	{
SEQAN_CHECKPOINT
		return me.x2;
	}

	friend inline int const & 
	_end1(_HirschbergSet const & me)
	{
SEQAN_CHECKPOINT
		return me.x2;
	}

	friend inline void
	_setEnd1(_HirschbergSet & me, int const & new_end)
	{
SEQAN_CHECKPOINT
		me.x2 = new_end;
	}
	// Sequence 2
/**
.Internal.begin2:
..cat:Functions
..summary:Returns the begin position of the _HirschberSet in the second sequence
..signature:begin2(_HirschbergSet)
..param._HirschbergSet: Reference to the _HirschbergSet
..param._HirschbergSet.type:Internal._HirschbergSet
*/
/**
.Internal.setBegin2:
..cat:Functions
..summary:Sets the begin position of the _HirschberSet in the second sequence
..signature:setBegin2(_HirschbergSet,new_begin)
..param._HirschbergSet: Reference to the _HirschbergSet
..param._HirschbergSet.type:Internal._HirschbergSet
..param.new_begin: The new begin position in the first sequence
*/
/**
.Internal.end2:
..cat:Functions
..summary:Returns the end position of the _HirschberSet in the second sequence
..signature:end2(_HirschbergSet)
..param._HirschbergSet: Reference to the _HirschbergSet
..param._HirschbergSet.type:Internal._HirschbergSet
*/
/**
.Internal.setEnd2:
..cat:Functions
..summary:Sets the end position of the _HirschberSet in the second sequence
..signature:setEnd2(_HirschbergSet,new_end)
..param._HirschbergSet: Reference to the _HirschbergSet
..param._HirschbergSet.type:Internal._HirschbergSet
..param.new_begin: The new end position in the second sequence
*/
	friend inline int &
	_begin2(_HirschbergSet & me)
	{
SEQAN_CHECKPOINT
		return me.y1;
	}

	friend inline int const &
	_begin2(_HirschbergSet const & me)
	{
SEQAN_CHECKPOINT
		return me.y1;
	}

	friend inline void
	_setBegin2(_HirschbergSet & me, int const & new_begin)
	{
SEQAN_CHECKPOINT
		me.y1 = new_begin;
	}

	friend inline int &
	_end2(_HirschbergSet & me)
	{
SEQAN_CHECKPOINT
		return me.y2;
	}

	friend inline int const &
	_end2(_HirschbergSet const & me)
	{
SEQAN_CHECKPOINT
		return me.y2;
	}

	friend inline void
	_setEnd2(_HirschbergSet & me, int const & new_end)
	{
SEQAN_CHECKPOINT
		me.y2 = new_end;
	}

	// //////////////////////////////////////////////////////////////////////////////////////////////
	// Score Methods
	// //////////////////////////////////////////////////////////////////////////////////////////////

	friend inline int &
	_score(_HirschbergSet & me)
	{
SEQAN_CHECKPOINT
		return me.score;
	}

	friend inline int const &
	_score(_HirschbergSet const & me)
	{
SEQAN_CHECKPOINT
		return me.score;
	}

	friend inline void
	_setScore(_HirschbergSet & me,int new_score)
	{
SEQAN_CHECKPOINT
		me.score = new_score;
	}

	// //////////////////////////////////////////////////////////////////////////////////////////////
	//  Debug Methods
	//		functions are only used for debugging or verbose output, therefore they
	//      are only active in SEQAN_DEBUG
	// //////////////////////////////////////////////////////////////////////////////////////////////
#ifdef SEQAN_DEBUG
	
	friend inline void
	print(_HirschbergSet const & me)
	{
		std::cout << me.x1 << " " << me.x2 << "\t" << me.y1 << " " << me.y2 << std::endl;
	}
#endif
};


inline bool
operator == (_HirschbergSet const & lhs, 
		_HirschbergSet const & rhs)
{
	return ((_begin1(lhs) == _begin1(rhs)) && (_end1(lhs) == _end1(rhs)) && (_begin2(lhs) == _begin2(rhs)) && (_end2(lhs) == _end2(rhs)));
}

}

#endif // #ifndef SEQAN_HEADER_HIRSCHBERG_SET_H
