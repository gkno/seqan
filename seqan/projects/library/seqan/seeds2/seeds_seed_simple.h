/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
 ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ============================================================================
  Specialization "Simple" for class Seed.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
#define SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// Already defined in module score.
// struct _Simple;
// typedef Tag<_Simple> Simple;

/**
.Spec.SimpleSeed:
..summary:Describes a seed with start and end position and diagonal upper and lower bounds.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, SimpleSeed>
..param.TPosition:The type number that schuld be used. Must have negative numbers (e.g. int/long).
.Memfunc.SimpleSeed#Seed:
..class:Spec.SimpleSeed
..summary:Constructor
..signature: Seed<TPosition, SimpleSeed> ()
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, length)
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, qEndPos, dEndPos)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.qEndPos: End in query sequence.
..param.dEndPos: End in database sequence.
..param.length: Length of the seed.
..include:seqan/seeds.h
*/
template <typename TPosition>
class Seed<TPosition, Simple>
{
public:
    TPosition _leftDim0;
    TPosition _leftDim1;
    TPosition _rightDim0;
    TPosition _rightDim1;
    TPosition _leftDiagonal;
    TPosition _rightDiagonal;

    Seed()
            : _leftDim0(0), _leftDim1(0), _rightDim0(0), _rightDim1(0),
              _leftDiagonal(0), _rightDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition leftDim0, TPosition leftDim1, TPosition seedLength)
            : _leftDim0(leftDim0),
              _leftDim1(leftDim1),
              _rightDim0(leftDim0 + seedLength - 1),
              _rightDim1(leftDim1 + seedLength - 1),
              _leftDiagonal(leftDim1 - leftDim0),
              _rightDiagonal(leftDim1 - leftDim0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition leftDim0, TPosition leftDim1, TPosition rightDim0,
         TPosition rightDim1)
            : _leftDim0(leftDim0),
              _leftDim1(leftDim1),
              _rightDim0(rightDim0),
              _rightDim1(rightDim1),
              _leftDiagonal(_max(leftDim1 - leftDim0, rightDim1 - rightDim0)),
              _rightDiagonal(_min(leftDim1 - leftDim0, rightDim1 - rightDim0))
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template<typename TPosition>
inline TPosition 
startDiagonal(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._leftDim1 - seed._leftDim0;
}

template<typename TPosition, typename TSpecSeed>
inline TPosition 
endDiagonal(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._rightDim1 - seed._rightDim0;
}

template< typename TPosition, typename TDimension>
inline TPosition 
leftPosition(Seed<TPosition, Simple> const & seed, TDimension dim)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(dim, static_cast<dim>(0));
    SEQAN_ASSERT_LEQ(dim, static_cast<dim>(1));
	return dim ? seed._leftDim1 : seed._leftDim0;
}

template< typename TPosition, typename TDimension, typename TPosition2> 
inline TPosition 
setLeftPosition(Seed<TPosition, Simple> & seed, 
				TSize dim,
				TPosition2 newLeftPosition)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(dim, static_cast<dim>(0));
    SEQAN_ASSERT_LEQ(dim, static_cast<dim>(1));
    if (dim)
        seed._leftDim1 = newLeftPosition;
	else
        seed._leftDim0 = newLeftPosition;
}

template< typename TPosition, typename TDimension>
inline TPosition 
rightPosition(Seed<TPosition, Simple> const & seed, TDimension dim)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(dim, static_cast<dim>(0));
    SEQAN_ASSERT_LEQ(dim, static_cast<dim>(1));
	return dim ? seed._rightDim1 : seed._rightDim0;
}

template< typename TPosition, typename TDimension, typename TPosition2> 
inline TPosition 
setRightPosition(Seed<TPosition, Simple> & seed, 
                 TSize dim,
                 TPosition2 newRightPosition)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(dim, static_cast<dim>(0));
    SEQAN_ASSERT_LEQ(dim, static_cast<dim>(1));
    if (dim)
        seed._rightDim1 = newRightPosition;
	else
        seed._rightDim0 = newRightPosition;
}

// TODO(holtgrew): Same as for Chained, maybe inherit from same base class Seed2D?
template <typename TPosition>
inline unsigned
dimension(Seed<TPosition, Simple> & seed)
{
    SEQAN_CHECKPOINT;
    return 2u;
}

template <typename TPosition>
inline TPosition 
leftDim0(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._leftDim0;
}

template <typename TPosition>
inline TPosition 
rightDim0(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;;
	return seed._rightDim0;
}

template <typename TPosition>
inline TPosition 
leftDim1(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._leftDim1;
}

template <typename TPosition>
inline TPosition 
rightDim1(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._rightDim1;
}

template<typename TPosition>
inline TPosition 
leftDiagonal(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._leftDiagonal;
}

template<typename TPosition>
inline TPosition 
rightDiagonal(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._rightDiagonal;
}

// TODO(holtgrew): Maybe rename to seedLength or so, name-clashes with containers' length().
template<typename TPosition>
inline TPosition 
length(Seed<TPosition, Simple> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._rightDim0 - seed._leftDim0 + 1;
}

template<typename TPosition>
inline void 
setLeftDim0(Seed<TPosition, Simple> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
	seed._leftDim0 = newLeftPosition;
}

template<typename TPosition>
inline void 
setRightDim0(Seed<TPosition, Simple> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
	seed._rightDim0 = newRightPosition;
}

template<typename TPosition>
inline void 
setLeftDim1(Seed<TPosition, Simple> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
	seed._leftDim1 = newLeftPosition;
}

template<typename TPosition>
inline void 
setRightDim1(Seed<TPosition, Simple> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
	seed._rightDim1 = newRightPosition;
}

template<typename TPosition>
inline void 
setLeftDiagonal(Seed<TPosition, Simple> & seed,
				TPosition newDiag)
{
	SEQAN_CHECKPOINT;
	seed._leftDiagonal = newDiag;
}

template<typename TPosition>
inline void 
setRightDiagonal(Seed<TPosition, Simple> & seed, 
				 TPosition newDiag)
{
	SEQAN_CHECKPOINT;
	seed._rightDiagonal = newDiag;
}

#endif  // SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
