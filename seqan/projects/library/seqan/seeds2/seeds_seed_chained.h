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

struct _Chained;
typedef Tag<_Chained> ChainedSeed;  // TODO(holtgrew): Chained already taken as template in file. Maybe prefer non-parameterized types for simpler names.

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
template <typename TConfig>
class Seed<ChainedSeed, TConfig>
        : TConfig::TScoreMixin
{
public:
    typedef typename TConfig::TPosition TPosition;
    typedef typename TConfig::TSize TSize;
    typedef typename TConfig::TDiagonal TDiagonal;

    typedef SeedDiagonal<TPosition, TSize> TSeedDiagonal;

    ::std::list<TSeedDiagonal> _seedDiagonals;
    TDiagonal _lowerDiagonal;
    TDiagonal _upperDiagonal;

    Seed() : _lowerDiagonal(0), _upperDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition leftDim0, TPosition leftDim1, TPosition seedLength)
            : _lowerDiagonal(leftDim1 - leftDim0),
              _upperDiagonal(leftDim1 - leftDim0)
    {
        SEQAN_CHECKPOINT;
        appendValue(_seedDiagonals, TSeedDiagonal(leftDim0, leftDim1, seedLength));
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

/**
.Metafunction.Value.param.T:Spec.ChainedSeed
 */
template <typename TConfig>
struct Value<Seed<ChainedSeed, TConfig> >
{
    typedef Seed<ChainedSeed, TConfig> _TSeed;
    typedef typename Position<_TSeed>::Type _TPosition;
    typedef typename Size<_TSeed>::Type _TSize;

    typedef SeedDiagonal<_TPosition, _TSize> Type;
};

template <typename TConfig>
struct Value<Seed<ChainedSeed, TConfig> const>
{
    typedef Seed<ChainedSeed, TConfig> _TSeed;
    typedef typename Position<_TSeed>::Type _TPosition;
    typedef typename Size<_TSeed>::Type _TSize;

    typedef SeedDiagonal<_TPosition, _TSize> const Type;
};

/**
.Metafunction.Iterator.param.T:Spec.ChainedSeed
 */
template <typename TConfig>
struct Iterator<Seed<ChainedSeed, TConfig>, Standard>
{
    typedef Seed<ChainedSeed, TConfig> _TSeed;
    typedef typename Value<_TSeed>::Type _TSeedDiagonal;
    typedef typename ::std::list<_TSeedDiagonal>::iterator Type;
};

template <typename TConfig>
struct Iterator<Seed<ChainedSeed, TConfig> const, Standard>
{
    typedef Seed<ChainedSeed, TConfig> _TSeed;
    typedef typename Value<_TSeed>::Type _TSeedDiagonal;
    typedef typename ::std::list<_TSeedDiagonal>::const_iterator Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getLeftDim0(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return front(seed._seedDiagonals).leftDim0;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getRightDim0(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return back(seed._seedDiagonals).leftDim0 + back(seed._seedDiagonals).length - 1;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getLeftDim1(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return front(seed._seedDiagonals).leftDim1;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getRightDim1(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return back(seed._seedDiagonals).leftDim1 + back(seed._seedDiagonals).length - 1;
}

/**
.Function.appendDiagonal
..summary: Adds diagonal to the Chained Seed.
..cat:Seed Handling
..signature:appendDiag(seed, diagonal)
..param.seed: The seed to which the diagonal should be added.
...type:Spec.ChainedSeed
..param.diag: The diagonal to add.
...type:Class.SeedDiagonal.
...remarks: A diagonal consists of three values: 1: start in 1. sequence, 2: start in 2. sequence, 3: length of match
*/
template <typename TConfig>
inline void
appendDiagonal(Seed<ChainedSeed, TConfig> & seed,
               typename Value<Seed<ChainedSeed, TConfig> >::Type const & diagonal)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Check for consistency!
    appendValue(seed._seedDiagonals, diagonal);
}

/**
.Function.begin.param.object.type:Spec.ChainedSeed
*/
template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> >::Type
begin(Seed<ChainedSeed, TConfig> & seed, Standard const &)
{
    SEQAN_CHECKPOINT;
    return seed._seedDiagonals.begin();
}

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> const>::Type
begin(Seed<ChainedSeed, TConfig> const & seed, Standard const &)
{
    SEQAN_CHECKPOINT;
    return seed._seedDiagonals.begin();
}

/**
.Function.end.param.object.type:Spec.ChainedSeed
*/
template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> >::Type
end(Seed<ChainedSeed, TConfig> & seed, Standard const &)
{
    SEQAN_CHECKPOINT;
    return seed._seedDiagonals.end();
}

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> const>::Type
end(Seed<ChainedSeed, TConfig> const & seed, Standard const &)
{
    SEQAN_CHECKPOINT;
    return seed._seedDiagonals.end();
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
