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

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// TODO(holtgrew): Document seed tags.
struct _Chained;
typedef Tag<_Chained> ChainedSeed;  // TODO(holtgrew): Chained already taken as template in file. Maybe prefer non-parameterized types for simpler names.

/**
.Spec.ChainedSeedSeed
..summary:Describes a seed with start and end position2 and diagonal upper and lower bounds. Additionaly diagonal segments
between start and end position2 are stored.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, ChainedSeedSeed>
..param.TPosition:The type of number that schuld be used. Must have negative numbers (e.g. int/long).
.Memfunc.ChainedSeedSeed#Seed:
..class:Spec.ChainedSeedSeed
..summary:Constructor
..signature: Seed<TPosition, ChainedSeedSeed> ()
..signature: Seed<TPosition, ChainedSeedSeed> (qStartPos, dStartPos, length)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.length: Length of the seed.
..include:seqan/seeds.h
*/
template <typename TPosition>
class Seed<TPosition, ChainedSeed>
{
public:
    typedef SeedDiagonal<TPosition, TPosition> _TSeedDiagonal;

    ::std::list<_TSeedDiagonal> _seedDiagonals;
    TPosition _leftDiagonal;
    TPosition _rightDiagonal;

    Seed() : _leftDiagonal(0), _rightDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition leftDim0, TPosition leftDim1, TPosition seedLength)
            : _leftDiagonal(leftDim1 - leftDim0),
              _rightDiagonal(leftDim1 - leftDim0)
    {
        SEQAN_CHECKPOINT;
        appendValue(_seedDiagonals, _TSeedDiagonal(leftDim0, leftDim1, seedLength));
    }

	Seed(TPosition leftDim0, TPosition leftDim1, TPosition rightDim0, TPosition rightDim1)
    {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_FAIL("Implement me!");
        (void)leftDim0;
        (void)leftDim1;
        (void)rightDim0;
        (void)rightDim1;
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

/**
.Metafunction.Value.param.T:Class.Seed
 */
template <typename TPosition>
struct Value<Seed<TPosition, ChainedSeed> >
{
    typedef SeedDiagonal<TPosition, TPosition> Type;
};

template <typename TPosition>
struct Value<Seed<TPosition, ChainedSeed> const>
{
    typedef SeedDiagonal<TPosition, TPosition> const Type;
};

// ===========================================================================
// Functions
// ===========================================================================

// From base class Seed<TPosition, TSpec>

template<typename TPosition>
inline TPosition 
startDiagonal(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    // TODO(holtgrew): front() is not specialized for std::list, apparently!
    return seed._seedDiagonals.front().leftDim1 - seed._seedDiagonals.front().leftDim0;
}

template<typename TPosition>
inline TPosition 
endDiagonal(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    // TODO(holtgrew): back() is not specialized for std::list, apparently!
//     return back(seed._seedDiagonals).leftDim1 + back(seed._seedDiagonals).length - 1;
    return seed._seedDiagonals.back().leftDim1 + seed._seedDiagonals.back().length - 1;
}

template< typename TPosition, typename TDimension>
inline TPosition 
leftPosition(Seed<TPosition, ChainedSeed> const & seed, TDimension dim)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return 0;
}

template< typename TPosition, typename TDimension>
inline TPosition 
rightPosition(Seed<TPosition, ChainedSeed> const & seed, TDimension dim)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return 0;
}

template <typename TPosition>
inline TPosition 
leftDim0(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return 0;
}

template <typename TPosition>
inline TPosition 
rightDim0(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return 0;
}

template <typename TPosition>
inline TPosition 
leftDim1(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    // TODO(holtgrew): front() is not specialized for std::list, apparently!
    return seed._seedDiagonals.front().leftDim1;
}

template <typename TPosition>
inline TPosition 
rightDim1(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    // TODO(holtgrew): back() is not specialized for std::list, apparently!
//     return back(seed._seedDiagonals).leftDim0 + back(seed._seedDiagonals).length - front(seed._seedDiagonals).leftDim0;
    return seed._seedDiagonals.back().leftDim0 + seed._seedDiagonals.back().length - seed._seedDiagonals.front().leftDim0;
}

// TODO(holtgrew): Same as Simple Seed.
template<typename TPosition>
inline TPosition 
leftDiagonal(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._leftDiagonal;
}

// TODO(holtgrew): Same as Simple Seed.
template<typename TPosition>
inline TPosition 
rightDiagonal(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._rightDiagonal;
}

// TODO(holtgrew): Maybe rename to seedLength or so, name-clashes with containers' length().
template<typename TPosition>
inline TPosition 
length(Seed<TPosition, ChainedSeed> const & seed)
{
	SEQAN_CHECKPOINT;
    return back(seed._seedDiagonals).leftDim0 + back(seed._seedSegmetns).length - front(seed._seedDiagonals).leftDim0;
}

template<typename TPosition>
inline void 
setLeftDim0(Seed<TPosition, ChainedSeed> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
    TPosition lengthDiff = front(seed._seedDiagonals).leftDim0 - newLeftPosition;
    front(seed._seedDiagonals).leftDim0 = newLeftPosition;
    front(seed._seedDiagonals).leftDim1 -= newLeftPosition;
    front(seed._seedDiagonals).length += newLeftPosition;
}

template<typename TPosition>
inline void 
setRightDim0(Seed<TPosition, ChainedSeed> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
    back(seed._seedDiagonals).length = newRightPosition - back(seed._seedDiagonals).leftDim0 + 1;
}

template<typename TPosition>
inline void 
setLeftDim1(Seed<TPosition, ChainedSeed> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
    TPosition lengthDiff = front(seed._seedDiagonals).leftDim1 - newLeftPosition;
    front(seed._seedDiagonals).leftDim0 -= newLeftPosition;
    front(seed._seedDiagonals).leftDim1 = newLeftPosition;
    front(seed._seedDiagonals).length += newLeftPosition;
}

template<typename TPosition>
inline void 
setRightDim1(Seed<TPosition, ChainedSeed> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
    back(seed._seedDiagonals).length = newRightPosition - back(seed._seedDiagonals).leftDim1 + 1;
}

// TODO(holtgrew): Same as Simple Seed.
template<typename TPosition>
inline void 
setLeftDiagonal(Seed<TPosition, ChainedSeed> & seed,
				TPosition newDiag)
{
	SEQAN_CHECKPOINT;
    seed._leftDiagonal = newDiag;
}

// TODO(holtgrew): Same as Simple Seed.
template<typename TPosition>
inline void 
setRightDiagonal(Seed<TPosition, ChainedSeed> & seed, 
				 TPosition newDiag)
{
	SEQAN_CHECKPOINT;
    seed._rightDiagonal = newDiag;
}

// TODO(holtgrew): Same as for Chained, maybe inherit from same base class Seed2D?
template <typename TPosition>
inline unsigned
dimension(Seed<TPosition, ChainedSeed> & seed)
{
    SEQAN_CHECKPOINT;
    return 2u;
}


// From specialization Seed<TPosition, ChainedSeed>

// TODO(holtgrew): Rename to appendDiagonal?
/**
.Function.appendDiag
..summary: Adds diagonal to the seed.
..cat:Seed Handling
..signature:appendDiag(seed, diagonal)
..param.seed: The seed to which the diagonal should be added.
...type:Spec.ChainedSeedSeed
..param.diag: The diagonal to add.
...type:Class.SeedDiagonal.
...remarks: A diagonal consists of three values: 1: start in 1. sequence, 2: start in 2. sequence, 3: length of match
*/
template <typename TPosition>
inline void
appendDiag(Seed<TPosition, ChainedSeed> & seed,
           typename Value<Seed<TPosition, ChainedSeed> >::Type const & diagonal)
{
    SEQAN_CHECKPOINT;
    appendValue(seed._seedDiagonals, diagonal);
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
