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
  Specialization "Chained" for class Seed.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
#define SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// TODO(holtgrew): Document seed tags.
struct _Chained;
typedef Tag<_Chained> Chained;

/**
.Spec.ChainedSeed
..summary:Describes a seed with start and end position2 and diagonal upper and lower bounds. Additionaly diagonal segments
between start and end position2 are stored.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, ChainedSeed>
..param.TPosition:The type of number that schuld be used. Must have negative numbers (e.g. int/long).
.Memfunc.ChainedSeed#Seed:
..class:Spec.ChainedSeed
..summary:Constructor
..signature: Seed<TPosition, ChainedSeed> ()
..signature: Seed<TPosition, ChainedSeed> (qStartPos, dStartPos, length)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.length: Length of the seed.
..include:seqan/seeds.h
*/
template <typename TPosition>
class Seed<TPosition, Chained>
{
public:
    typedef SeedSegment<TPosition, TPosition> TSeedSegment;

    ::std::list<_TSeedSegment> _seedSegments;
    TPosition _leftDiagonal;
    TPosition _rightDiagonal;

    Seed() : _leftDiagonal(0), _rightDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition leftDim0, TPosition leftDim1, TPosition seedLength)
    {
        SEQAN_CHECKPONT;
        appendValue(_seedSegments, TSeedSegment(leftDim0, leftDim1, seedLength));
        _leftDiagonal = leftDim1 - leftDim0;
        _rightDiagonal = _leftDiagonal;
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// TODO(holtgrew): Document this.
template <typename TPosition>
struct SeedSegment<Seed<TPosition, Chained> >
{
    typedef SeedSegment<TPosition, TPosition> Type;
};

template <typename TPosition>
struct SeedSegment<Seed<TPosition, Chained> const>
{
    typedef SeedSegment<TPosition, TPosition> const Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template<typename TPosition>
inline TPosition 
startDiagonal(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;
    return front(seed._seedSegments).leftDim1 - front(seed._seedSegments).leftDim0;
}

template<typename TPosition, typename TSpecSeed>
inline TPosition 
endDiagonal(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;
    return back(seed._seedSegments).leftDim1 - back(seed._seedSegments).leftDim0;
}

template< typename TPosition, typename TDimension>
inline TPosition 
leftPosition(Seed<TPosition, Chained> const & seed, TDimension dim)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(dim, static_cast<dim>(0));
    SEQAN_ASSERT_LEQ(dim, static_cast<dim>(1));
    if (dim)
        return leftDim1(seed);
    else
        return leftDim0(seed);
}

template< typename TPosition, typename TDimension, typename TPosition2> 
inline TPosition 
setLeftPosition(Seed<TPosition, Chained> & seed, 
				TSize dim,
				TPosition2 newLeftPosition)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(dim, static_cast<dim>(0));
    SEQAN_ASSERT_LEQ(dim, static_cast<dim>(1));
    if (dim)
        return setLeftDim1(seed, newLeftPosition);
    else
        return setLeftDim0(seed, newLeftPosition);
}

template< typename TPosition, typename TDimension>
inline TPosition 
rightPosition(Seed<TPosition, Chained> const & seed, TDimension dim)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(dim, static_cast<dim>(0));
    SEQAN_ASSERT_LEQ(dim, static_cast<dim>(1));
    if (dim)
        return rightDim1(seed);
    else
        return rightDim0(seed);
}

template< typename TPosition, typename TDimension, typename TPosition2> 
inline TPosition 
setRightPosition(Seed<TPosition, Chained> & seed, 
                 TSize dim,
                 TPosition2 newRightPosition)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(dim, static_cast<dim>(0));
    SEQAN_ASSERT_LEQ(dim, static_cast<dim>(1));
    if (dim)
        return setRightDim1(seed, newRightPosition);
    else
        return setRightDim0(seed, newRightPosition);
}

// TODO(holtgrew): Same as for Simple, maybe inherit from same base class Seed2D?
template <typename TPosition>
inline unsigned
dimension(Seed<TPosition, Chained> & me)
{
    SEQAN_CHECKPOINT;
    return 2u;
}

template <typename TPosition>
inline TPosition 
leftDim0(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;
    return front(seed._seedSegments).leftDim0;
}

template <typename TPosition>
inline TPosition 
rightDim0(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;;
    return back(seed._seedSegments).leftDim1 + back(seed._seedSegments).length - 1;
}

template <typename TPosition>
inline TPosition 
leftDim1(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;
    return front(seed._seedSegments).leftDim1;
}

template <typename TPosition>
inline TPosition 
rightDim1(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;
    return back(seed._seedSegments).leftDim0 + back(seed._seedSegments).length - front(seed._seedSegments).leftDim0;
}

// TODO(holtgrew): Same as Simple Seed.
template<typename TPosition>
inline TPosition 
leftDiagonal(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._leftDiagonal;
}

// TODO(holtgrew): Same as Simple Seed.
template<typename TPosition>
inline TPosition 
rightDiagonal(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._rightDiagonal;
}

// TODO(holtgrew): Maybe rename to seedLength or so, name-clashes with containers' length().
template<typename TPosition>
inline TPosition 
length(Seed<TPosition, Chained> const & seed)
{
	SEQAN_CHECKPOINT;
    return back(seed._seedSegments).leftDim0 + back(seed._seedSegmetns).length - front(seed._seedSegments).leftDim0;
}

template<typename TPosition>
inline void 
setLeftDim0(Seed<TPosition, Chained> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
    TPosition lengthDiff = front(seed._seedSegments).leftDim0 - newLeftPosition;
    front(seed._seedSegments).leftDim0 = newLeftPosition;
    front(seed._seedSegments).leftDim1 -= newLeftPosition;
    front(seed._seedSegments).length += newLeftPosition;
}

template<typename TPosition>
inline void 
setRightDim0(Seed<TPosition, Chained> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
    back(seed._seedSegments).length = newRightPosition - back(seed._seedSegments).leftDim0 + 1;
}

template<typename TPosition>
inline void 
setLeftDim1(Seed<TPosition, Chained> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
    TPosition lengthDiff = front(seed._seedSegments).leftDim1 - newLeftPosition;
    front(seed._seedSegments).leftDim0 -= newLeftPosition;
    front(seed._seedSegments).leftDim1 = newLeftPosition;
    front(seed._seedSegments).length += newLeftPosition;
}

template<typename TPosition>
inline void 
setRightDim1(Seed<TPosition, Chained> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
    back(seed._seedSegments).length = newRightPosition - back(seed._seedSegments).leftDim1 + 1;
}

// TODO(holtgrew): Same as Simple Seed.
template<typename TPosition>
inline void 
setLeftDiagonal(Seed<TPosition, Chained> & seed,
				TPosition newDiag)
{
	SEQAN_CHECKPOINT;
    seed._leftDiagonal = newDiag;
}

// TODO(holtgrew): Same as Simple Seed.
template<typename TPosition>
inline void 
setRightDiagonal(Seed<TPosition, Chained> & seed, 
				 TPosition newDiag)
{
	SEQAN_CHECKPOINT;
    seed._rightDiagonal = newDiag;
}

// TODO(holtgrew): Rename to appendDiagonal?
/**
.Function.appendDiag
..summary: Adds diagonal to the seed.
..cat:Seed Handling
..signature:appendDiag(seed, diagonal)
..param.seed: The seed to which the diagonal should be added.
...type:Spec.ChainedSeed
..param.diag: The diagonal to add.
...type:Class.SeedSegment.
...remarks: A diagonal consists of three values: 1: start in 1. sequence, 2: start in 2. sequence, 3: length of match
*/
template <typename TPosition>
inline void
appendDiag(Seed<TPosition, Chained> & seed,
           typename SeedSegment<Seed<TPosition, Chained> >::Type const & diagonal)
{
    SEQAN_CHECKPOINT;
    appendValue(seed._seedSegments, diagonal);
}

#endif  // SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
