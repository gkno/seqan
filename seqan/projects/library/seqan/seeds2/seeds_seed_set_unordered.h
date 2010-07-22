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
  The Unordered specialization of the class SeedSet.  Seeds are stored
  in a string and not kept in a particular order.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_SET_UNORDERED_H_
#define SEQAN_SEEDS_SEEDS_SEED_SET_UNORDERED_H_

#include <cmath>

namespace seqan {

template <typename T>
T _abs(T const & x)
{
    if (x < static_cast<T>(0))
        return -x;
    else
        return x;
}

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

struct _Unordered;
typedef Tag<_Unordered> Unordered;

// TODO(holtgrew): Maybe allow iterating over seeds that have reached a certain quality (length/score).

template <typename TSeedSpec, typename TSeedSetConfig>
class SeedSet<TSeedSpec, Unordered, TSeedSetConfig>
        : public TSeedSetConfig::TQualityThresholdMixin
{
public:
    typedef typename TSeedSetConfig::TQualityThresholdMixin _TQualityThresholdMixin;
    typedef typename TSeedSetConfig::TSeedConfig TSeedConfig;
    typedef Seed<TSeedSpec, TSeedConfig> TSeed;

    typedef Allocator<SinglePool<sizeof(TSeed) > > TSeedAllocator;

    typedef String<TSeed *> TAllSeeds;
    typedef std::set<TSeed *> THighQualitySeeds;

    TSeedAllocator _seedAllocator;
    TAllSeeds _allSeeds;
    // TODO(holtgrew): High quality seeds are only necessary if qualities are considered.
    THighQualitySeeds _highQualitySeeds;

    SeedSet()
            : _TQualityThresholdMixin()
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TSeedSpec, typename TSeedSetConfig>
struct Position<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >
{
    typedef SeedSet<TSeedSpec, Unordered, TSeedSetConfig> _TSeedSet;
    typedef String<typename _TSeedSet::TSeed> _TSeedString;
    typedef typename Position<_TSeedString>::Type Type;
};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Position<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>
        : Position<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> > {};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Size<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >
{
    typedef SeedSet<TSeedSpec, Unordered, TSeedSetConfig> _TSeedSet;
    typedef String<typename _TSeedSet::TSeed> _TSeedString;
    typedef typename Size<_TSeedString>::Type Type;
};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Size<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>
        : Size<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> > {};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >
{
    typedef typename TSeedSetConfig::TSeedConfig _TSeedConfig;
    typedef Seed<TSeedSpec, _TSeedConfig> _TSeed;
    typedef _TSeed Type;
};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>
{
    typedef typename TSeedSetConfig::TSeedConfig _TSeedConfig;
    typedef Seed<TSeedSpec, _TSeedConfig> _TSeed;
    typedef _TSeed const Type;
};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Reference<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >
{
    typedef typename TSeedSetConfig::TSeedConfig _TSeedConfig;
    typedef Seed<TSeedSpec, _TSeedConfig> _TSeed;
    typedef _TSeed & Type;
};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Reference<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>
{
    typedef typename TSeedSetConfig::TSeedConfig _TSeedConfig;
    typedef Seed<TSeedSpec, _TSeedConfig> _TSeed;
    typedef _TSeed const & Type;
};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig>, Standard>
{
    typedef typename TSeedSetConfig::TSeedConfig _TSeedConfig;
    typedef Seed<TSeedSpec, _TSeedConfig> _TSeed;
    typedef SeedSet<TSeedSpec, Unordered, TSeedSetConfig> _TSeedSet;
//     typedef String<_TSeed> _TSeedString;
//     typedef typename Iterator<_TSeedString, Standard>::Type _TIterator;
//     typedef _TIterator Type;
    typedef Iter<_TSeedSet, Indirect<typename std::set<_TSeed *>::iterator> > Type;
};

template <typename TSeedSpec, typename TSeedSetConfig>
struct Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const, Standard>
{
    typedef typename TSeedSetConfig::TSeedConfig _TSeedConfig;
    typedef Seed<TSeedSpec, _TSeedConfig> _TSeed;
    typedef SeedSet<TSeedSpec, Unordered, TSeedSetConfig> _TSeedSet;
//     typedef String<_TSeed const> _TSeedString;
//     typedef typename Iterator<_TSeedString const, Standard>::Type _TIterator;
//     typedef _TIterator Type;
    typedef Iter<_TSeedSet const, Indirect<typename std::set<_TSeed *>::const_iterator> > Type;
};

// ===========================================================================
// Functions
// ===========================================================================

// Standard Container Functions

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Size<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type
length(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Do not use dot function.
    return seedSet._highQualitySeeds.size();
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Size<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>::Type
length(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const & seedSet)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Do not use dot function.
    return seedSet._highQualitySeeds.size();
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type
begin(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type TIterator;
    // TODO(holtgrew): Do not use dot-method.
    return TIterator(seedSet._highQualitySeeds.begin());
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>::Type
begin(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const & seedSet, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>::Type TIterator;
    // TODO(holtgrew): Do not use dot-method.
    return TIterator(seedSet._highQualitySeeds.begin());
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type
end(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type TIterator;
    // TODO(holtgrew): Do not use dot-method.
    return TIterator(seedSet._highQualitySeeds.end());
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>::Type
end(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const & seedSet, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>::Type TIterator;
    // TODO(holtgrew): Do not use dot-method.
    return TIterator(seedSet._highQualitySeeds.end());
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type
front(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Do not use dot-method.
    return **seedSet._highQualitySeeds.begin();
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>::Type
front(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const & seedSet)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Do not use dot-method.
    return **seedSet._highQualitySeeds.begin();
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type
back(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Do not use dot-method.
    return **seedSet._highQualitySeeds.rbegin();
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const>::Type
back(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> const & seedSet)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Do not use dot-method.
    return **seedSet._highQualitySeeds.rbegin();
}

// SeedSet Functions

// TODO(holtgrew): Add bulk-addSeeds functions.

// Returns true iff b can be merged into a where a is the one to the
// upper left, b the one to the lower right.
template <typename TSeedSpec, typename TSeedConfig, typename TThreshold>
inline bool
_seedsMergeable(Seed<TSeedSpec, TSeedConfig> const & a,
                Seed<TSeedSpec, TSeedConfig> const & b,
                TThreshold const & threshold,
                Merge const &)
{
    // TODO(holtgrew): TThreshold could be Position<TSeed>::Type.
    SEQAN_CHECKPOINT;
    // If the two seeds do not overlap, they cannot be merged.
    if (getBeginDim0(b) >= getEndDim0(a) && getBeginDim1(b) >= getEndDim1(a))
        return false;
    // If the distance between the diagonals exceeds the threshold
    // then the seeds cannot be merged.
    typedef typename _MakeUnsigned<TThreshold>::Type TUnsignedThreshold;
    if (static_cast<TUnsignedThreshold>(_abs(getEndDiagonal(a) - getStartDiagonal(b))) >= static_cast<TUnsignedThreshold>(threshold))
        return false;
    // Otherwise, the seeds can be merged.
    return true;
}

template <typename TSeedConfig>
inline void
_mergeSeeds(Seed<Simple, TSeedConfig> & seed,
            Seed<Simple, TSeedConfig> const & other,
            Merge const &)
{
    _updateMergeScore(seed, other);

    // Merging simple seeds simply works by updating their coordinates
    // to the minimal/maximal values.
    setBeginDim0(seed, _min(getBeginDim0(seed), getBeginDim0(other)));
    setBeginDim1(seed, _min(getBeginDim1(seed), getBeginDim1(other)));
    setEndDim0(seed, _max(getEndDim0(seed), getEndDim0(other)));
    setEndDim1(seed, _max(getEndDim1(seed), getEndDim1(other)));
    setLowerDiagonal(seed, _min(getLowerDiagonal(seed), getLowerDiagonal(other)));
    setUpperDiagonal(seed, _max(getUpperDiagonal(seed), getUpperDiagonal(other)));
}

template <typename TSeedConfig>
inline void
_mergeSeeds(Seed<ChainedSeed, TSeedConfig> & seed,
            Seed<ChainedSeed, TSeedConfig> const & other,
            Merge const &)
{
    // For chained seeds, we first remove all diagonals from seed
    // until the last diagonal of seed starts truly before other.
    // Then, we possibly shorten the last diagonal.  Finally, we copy
    // over all diagonals from other.

    // Remove diagonals.
    typedef Seed<ChainedSeed, TSeedConfig> TSeed;
    typedef typename Iterator<TSeed, Standard>::Type TIterator;
    TIterator it;
    TIterator lastKept = begin(seed);
    for (it = begin(seed); it != end(seed); ++it) {
        if (it->beginDim0 < getBeginDim0(other) &&
            it->beginDim1 < getBeginDim1(other))
            break;
        lastKept = it;
    }
    truncateDiagonals(seed, it);

    // Shorten last diagonal if necessary.
    if (it->beginDim0 + it->length > getBeginDim0(other) &&
        it->beginDim1 + it->length > getBeginDim1(other)) {
        it->length = _min(getBeginDim0(other) - it->beginDim0, getBeginDim1(other) - it->beginDim1);
    } else if (it->beginDim0 + it->length > getBeginDim0(other)) {
        it->length = getBeginDim0(other) - it->beginDim0;
    } else if (it->beginDim1 + it->length > getBeginDim1(other)) {
        it->length = getBeginDim1(other) - it->beginDim1;
    }

    // Copy over other diagonals.
    typedef typename Iterator<TSeed const, Standard>::Type TConstIterator;
    for (TConstIterator it = begin(other, Standard()); it != end(other, Standard()); ++it)
        appendDiagonal(seed, *it);
}

template <typename TSeedSpec, typename TSeedSetConfig, typename TThreshold>
inline bool
addSeed(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet,
        typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type const & seed,
        TThreshold const & maxDiagDist,
        Merge const &)
{
    SEQAN_CHECKPOINT;

    typedef SeedSet<TSeedSpec, Unordered, TSeedSetConfig> TSeedSet;
    typedef typename TSeedSet::TAllSeeds TAllSeeds;
    typedef typename Iterator<TAllSeeds, Standard>::Type TSeedPtrIterator;
    typedef typename Value<TSeedSet>::Type TSeed;

    // Iterate over all seeds and search for the first one in this
    // arbitrary order that is mergeable with parameter seed within a
    // maximal diagonal distance maxDiagDist.  We allow either seed to
    // be the left one.
    //
    // TODO(holtgrew): Search for *closest* overlapping one instead!
    for (TSeedPtrIterator it = begin(seedSet._allSeeds); it != end(seedSet._allSeeds); ++it) {
        bool mergedThisSeed = false;
        if (_seedsMergeable(*value(it), seed, maxDiagDist, Merge())) {
            // Merge seed into *it.
            _mergeSeeds(*value(it), seed, Merge());
            mergedThisSeed = true;
        } else if (_seedsMergeable(seed, *value(it), maxDiagDist, Merge())) {
            // Merge *it into seed.
            TSeed tmp;
            move(tmp, *value(it));
            assign(*value(it), seed);
            _mergeSeeds(*value(it), tmp, Merge());
            mergedThisSeed = true;
        }
        if (mergedThisSeed) {
            // Possibly add resulting seed to the subset of seeds with
            // high enough score.
			typedef typename TSeedSetConfig::TQualityThreshold TQualityThreshold;
            if (_qualityReached(*value(it), seedSet, TQualityThreshold()))
                // TODO(holtgrew): Do not use dot-methods.
                seedSet._highQualitySeeds.insert(value(it));
            return true;
        }
    }
    return false;
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline void
addSeed(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet,
        typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type const & seed,
        Chaos const &)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Write me!");
}


template <typename TSeedSpec, typename TSeedSetConfig>
inline void
addSeed(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet,
        typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type const & seed,
        SimpleChain const &)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Write me!");
}

template <typename TSeedSpec, typename TSeedSetConfig>
inline void
addSeed(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet,
        typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type const & seed,
        Single const &)
{
    SEQAN_CHECKPOINT;
    typedef SeedSet<TSeedSpec, Unordered, TSeedSetConfig> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    // Allocate space for new seed in allocator and assign seed from
    // parameter to this space.
    //
    // TODO(holtgrew): Move would be faster if it is a chained seed with many diagonals.
    TSeed * tmp;
    allocate(seedSet._seedAllocator, tmp, 1);
    *tmp = seed;

    appendValue(seedSet._allSeeds, tmp);
	typedef typename TSeedSetConfig::TQualityThreshold TQualityThreshold;
    if (_qualityReached(seed, seedSet, TQualityThreshold()))
        // TODO(holtgrew): Do not use dot-methods.
        seedSet._highQualitySeeds.insert(tmp);
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SET_UNORDERED_H_

