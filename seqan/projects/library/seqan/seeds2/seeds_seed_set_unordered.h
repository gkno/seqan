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

template <typename TSeedSpec, typename TSeedSetConfig, typename TDistanceThreshold, typename TBandwidth, typename TCombination>
bool
_findSeedForCombination(
        typename Iterator<typename SeedSet<TSeedSpec, Unordered, TSeedSetConfig>::TAllSeeds, Standard>::Type & mergePartner,
        bool & seedIsOnTheLeft,
        SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet,
        typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type const & seed,
        TDistanceThreshold const & maxDistance,
        TBandwidth const & bandwidth,
        TCombination const & tag)
{
    SEQAN_CHECKPOINT;

    typedef typename Iterator<typename SeedSet<TSeedSpec, Unordered, TSeedSetConfig>::TAllSeeds, Standard>::Type TSeedPtrIterator;

    // Iterate over all seeds and search for the first one in this
    // arbitrary order that is combineable with parameter seed within
    // a maximal diagonal distance maxDistance.  We allow either seed
    // to be the left one.
    //
    // TODO(holtgrew): Search for *closest* overlapping one instead!
    for (TSeedPtrIterator it = begin(seedSet._allSeeds); it != end(seedSet._allSeeds); ++it) {
        if (_seedsCombineable(*value(it), seed, maxDistance, bandwidth, tag)) {
//			std::cout << "Combineable: " << (*value(it)) << " and " << seed << std::endl;
            // seed is to be merged into *it.
            mergePartner = it;
            seedIsOnTheLeft = false;
            return true;
        } else if (_seedsCombineable(seed, *value(it), maxDistance, bandwidth, tag)) {
//			std::cout << "Combineable: " << seed << " and " << (*value(it)) << std::endl;
            // *it is to be merged into seed.
            mergePartner = it;
            seedIsOnTheLeft = true;
            return true;
        }
    }

    // Found no seed to combine with.
    return false;
}

// TODO(holtgrew): Score not needed for Merge!

template <typename TSeedSpec, typename TSeedSetConfig, typename TDistanceThreshold, typename TBandwidth, typename TScoreValue, typename TSequence0, typename TSequence1, typename TCombination>
inline bool
addSeed(SeedSet<TSeedSpec, Unordered, TSeedSetConfig> & seedSet,
        typename Value<SeedSet<TSeedSpec, Unordered, TSeedSetConfig> >::Type const & seed,
        TDistanceThreshold const & maxDiagDist,
        TBandwidth const & bandwidth,
        Score<TScoreValue, Simple> const & scoringScheme,
        TSequence0 const & sequence0,
        TSequence1 const & sequence1,
        TCombination const & tag)
{
    SEQAN_CHECKPOINT;

    typedef SeedSet<TSeedSpec, Unordered, TSeedSetConfig> TSeedSet;
    typedef typename TSeedSet::TAllSeeds TAllSeeds;
    typedef typename Iterator<TAllSeeds, Standard>::Type TSeedPtrIterator;
    typedef typename Value<TSeedSet>::Type TSeed;

    // Try to find a seed for recombination.
    TSeedPtrIterator it = 0;
    bool seedIsOnTheLeft = false;
    bool foundSeed = _findSeedForCombination(it, seedIsOnTheLeft, seedSet, seed, maxDiagDist, bandwidth, tag);

    // If we could find a seed: Combine them.
    if (foundSeed) {
        // Swap seed and *value(it) if seed is on the left and
        // *value(it) is on the right.  Then, merge both.
        if (!seedIsOnTheLeft) {
            _combineSeeds(*value(it), seed, scoringScheme, sequence0, sequence1, tag);
        } else {
            TSeed tmp;
            move(tmp, *value(it));
            assign(*value(it), seed);
            _combineSeeds(*value(it), tmp, scoringScheme, sequence0, sequence1, tag);
        }
        // If the new seed has a high enough quality, add it to the
        // set of high-scoring seeds.
        typedef typename TSeedSetConfig::TQualityThreshold TQualityThreshold;
        if (_qualityReached(*value(it), seedSet, TQualityThreshold())) {
            // TODO(holtgrew): Do not use dot-methods.
            seedSet._highQualitySeeds.insert(value(it));
        }
        return true;
    }
    return false;
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

    // Allocate space for new seed in allocator and copy construct the
    // seed there.
    //
    // TODO(holtgrew): Move would be faster if it is a chained seed with many diagonals.
    TSeed * tmp;
    allocate(seedSet._seedAllocator, tmp, 1);
    new(tmp) TSeed(seed);

    appendValue(seedSet._allSeeds, tmp);
	typedef typename TSeedSetConfig::TQualityThreshold TQualityThreshold;
    if (_qualityReached(seed, seedSet, TQualityThreshold()))
        // TODO(holtgrew): Do not use dot-methods.
        seedSet._highQualitySeeds.insert(tmp);
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SET_UNORDERED_H_

