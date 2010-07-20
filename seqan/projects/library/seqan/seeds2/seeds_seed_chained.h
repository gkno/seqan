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
        : public TConfig::TScoreMixin
{
public:
    typedef typename TConfig::TPosition TPosition;
    typedef typename TConfig::TSize TSize;
    typedef typename TConfig::TDiagonal TDiagonal;

    typedef typename TConfig::TScoreMixin TScoreMixin;

    typedef SeedDiagonal<TPosition, TSize> TSeedDiagonal;

    ::std::list<TSeedDiagonal> _seedDiagonals;
    TDiagonal _lowerDiagonal;
    TDiagonal _upperDiagonal;

    Seed() : TScoreMixin(), _lowerDiagonal(0), _upperDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition beginDim0, TPosition beginDim1, TPosition seedLength)
            : TScoreMixin(),
              _lowerDiagonal(beginDim1 - beginDim0),
              _upperDiagonal(beginDim1 - beginDim0)
    {
        SEQAN_CHECKPOINT;
        appendValue(_seedDiagonals, TSeedDiagonal(beginDim0, beginDim1, seedLength));
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
.Metafunction.Size.param.T:Spec.ChainedSeed
 */
template <typename TConfig>
struct Size<Seed<ChainedSeed, TConfig> >
{
    typedef typename TConfig::TSize Type;
};

template <typename TConfig>
struct Size<Seed<ChainedSeed, TConfig> const>
        : Size<Seed<ChainedSeed, TConfig> > {};

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

/**
.Metafunction.Reference.param.T:Spec.ChainedSeed
 */
template <typename TConfig>
struct Reference<Seed<ChainedSeed, TConfig> >
{
    typedef Seed<ChainedSeed, TConfig> _TSeed;
    typedef typename Value<_TSeed>::Type _TSeedDiagonal;
    typedef _TSeedDiagonal & Type;
};

template <typename TConfig>
struct Reference<Seed<ChainedSeed, TConfig> const>
{
    typedef Seed<ChainedSeed, TConfig> _TSeed;
    typedef typename Value<_TSeed>::Type _TSeedDiagonal;
    typedef _TSeedDiagonal & Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TStream, typename TConfig>
inline TStream &
operator<<(TStream & stream, Seed<ChainedSeed, TConfig> const & seed)
{
    typedef Seed<ChainedSeed, TConfig> const TSeed;
    typedef typename Iterator<TSeed>::Type TIterator;

    stream << "Seed<ChainedSeed, TConfig>([";
    for (TIterator it = begin(seed); it != end(seed); ++it) {
        if (it != begin(seed))
            stream << ", ";
        stream << *it;
    }
    stream << "])";
    return stream;
}

template <typename TConfig>
inline bool
operator==(Seed<ChainedSeed, TConfig> const & a, Seed<ChainedSeed, TConfig> const & b)
{
    SEQAN_CHECKPOINT;
    return a._seedDiagonals == b._seedDiagonals &&
            a._upperDiagonal == b._upperDiagonal &&
            a._lowerDiagonal == b._lowerDiagonal;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getBeginDim0(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return front(seed._seedDiagonals).beginDim0;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getEndDim0(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return back(seed._seedDiagonals).beginDim0 + back(seed._seedDiagonals).length;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getBeginDim1(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return front(seed._seedDiagonals).beginDim1;
}

template <typename TConfig>
inline typename Position<Seed<ChainedSeed, TConfig> >::Type
getEndDim1(Seed<ChainedSeed, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return back(seed._seedDiagonals).beginDim1 + back(seed._seedDiagonals).length;
}

/**
.Function.length.param.object.type:Spec.ChainedSeed
*/
template <typename TConfig>
inline typename Size<Seed<ChainedSeed, TConfig> >::Type
length(Seed<ChainedSeed, TConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    return length(seed._seedDiagonals);
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
.Function.truncateDiagonals
..summary:Removes diagonals from the given first one to the end of the seed's diagonals.
..cat:Seed Handling
..signature:truncateDiagonals(seed, first)
..param.seed: The seed to which the diagonal should be added.
...type:Spec.ChainedSeed
..param.first: Iterator the first diagonal to remove.
*/
template <typename TConfig>
inline void
truncateDiagonals(Seed<ChainedSeed, TConfig> & seed,
                  typename Iterator<Seed<ChainedSeed, TConfig> >::Type const & first)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Add erase() to std::list adaptors?
    seed._seedDiagonals.erase(first, seed._seedDiagonals.end());
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
.Function.front.param.object.type:Spec.ChainedSeed
*/
template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> >::Type
front(Seed<ChainedSeed, TConfig> & seed)
{
    SEQAN_CHECKPOINT;
    return front(seed._seedDiagonals);
}

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> const>::Type
front(Seed<ChainedSeed, TConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    return front(seed._seedDiagonals);
}

/**
.Function.back.param.object.type:Spec.ChainedSeed
*/
template <typename TConfig>
inline typename Reference<Seed<ChainedSeed, TConfig> >::Type
back(Seed<ChainedSeed, TConfig> & seed)
{
    SEQAN_CHECKPOINT;
    return back(seed._seedDiagonals);
}

template <typename TConfig>
inline typename Iterator<Seed<ChainedSeed, TConfig> const>::Type
back(Seed<ChainedSeed, TConfig> const & seed)
{
    SEQAN_CHECKPOINT;
    return back(seed._seedDiagonals);
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

// Basic Functions

template <typename TConfig>
void
move(Seed<ChainedSeed, TConfig> & target, Seed<ChainedSeed, TConfig> & source)
{
    SEQAN_CHECKPOINT;
    std::swap(target._seedDiagonals, source._seedDiagonals);
    target._lowerDiagonal = source._lowerDiagonal;
    target._upperDiagonal = source._upperDiagonal;
}

// Debug Output

template <typename TStream, typename TConfig>
inline void
_write(TStream & stream, Seed<ChainedSeed, TConfig> const & seed, _Tikz const &)
{
    // Overall seed.
    stream << "\\draw[seed] (" << getBeginDim1(seed) << ", -" << getBeginDim0(seed) << ") -- (" << (getEndDim1(seed) - 1) << ", -" << (getEndDim0(seed) - 1) << ");" << std::endl;
    // Diagonals.
    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Iterator<TSeed const, Standard>::Type TIterator;
    for (TIterator it = begin(seed); it != end(seed); ++it)
        stream << "\\draw[seed diagonal] (" << it->beginDim1 << ", -" << it->beginDim0 << ") -- (" << (it->beginDim1 + it->length - 1) << ", -" << (it->beginDim0 + it->length - 1) << ");" << std::endl;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_CHAINED_H_
