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

// TODO(holtgrew): These functions should use "combine" instead of "merge".

// TODO(holtgrew): Add bulk-addSeeds functions.

// Returns true iff b can be merged into a where a is the one to the
// upper left, b the one to the lower right.
template <typename TSeedSpec, typename TSeedConfig, typename TThreshold>
inline bool
_seedsCombineable(Seed<TSeedSpec, TSeedConfig> const & a,
                  Seed<TSeedSpec, TSeedConfig> const & b,
                  TThreshold const & maxDiagonalDistance,
                  Nothing const & /*maxBandwidth*/,
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
    if (static_cast<TUnsignedThreshold>(_abs(getEndDiagonal(a) - getStartDiagonal(b))) >= static_cast<TUnsignedThreshold>(maxDiagonalDistance))
        return false;
    // Otherwise, the seeds can be merged.
    return true;
}


// Returns true iff b can be simple-chained to a where a is the one to
// the upper left, b the one to the lower right.
template <typename TSeedSpec, typename TSeedConfig, typename TThreshold>
inline bool
_seedsCombineable(Seed<TSeedSpec, TSeedConfig> const & a,
                  Seed<TSeedSpec, TSeedConfig> const & b,
                  TThreshold const & maxGapSize,
                  Nothing const & /*maxBandwidth*/,
                  SimpleChain const &)
{
    // TODO(holtgrew): We should be able to configure whether we want to have Manhattan, euclidean, minimal edit distance, for seeds.
    // TODO(holtgrew): TThreshold could be Position<TSeed>::Type.
    SEQAN_CHECKPOINT;

    // b has to be right of a for the two seeds to be chainable.
    if (getBeginDim0(b) < getEndDim0(a) || getBeginDim1(b) < getEndDim1(a))
        return false;

    // Distance is maximal distance, this corresponds to going the
    // distacen in the smaller distance with matches/mismatches and
    // the rest with indels.
    TThreshold distance = _max(getBeginDim0(b) - getEndDim0(a), getBeginDim1(b) - getEndDim1(a));
    // Compare distance with threshold.
    return distance <= maxGapSize;
}


// Returns true iff b can be Chaos chained to a where a is the one to
// the upper left, b the one to the lower right.
//
// TODO(holtgrew): Replace bandwidth with diagonalDistance.
template <typename TSeedSpec, typename TSeedConfig, typename TDistanceThreshold, typename TBandwidthThreshold>
inline bool
_seedsCombineable(Seed<TSeedSpec, TSeedConfig> const & a,
                  Seed<TSeedSpec, TSeedConfig> const & b,
                  TDistanceThreshold const & maxGapSize,
                  TBandwidthThreshold const & bandwidth,
                  Chaos const &)
{
    SEQAN_CHECKPOINT;

    // b has to be right of a for the two seeds to be chainable.
    if (getBeginDim0(b) < getEndDim0(a) || getBeginDim1(b) < getEndDim1(a))
        return false;

    // TODO(holtgrew): Ignore gray areas from CHAOS paper "so sequence order does not matter".

    // The diagonal distance has to be smaller than the bandwidth.
    // TODO(holtgrew): s/getStartDiagonal/getBeginDiagonal/
    TBandwidthThreshold diagonalDistance = _abs(getEndDiagonal(b) - getStartDiagonal(a));
    if (diagonalDistance > bandwidth)
        return false;

    // Distance is maximal distance, this corresponds to going the
    // distacen in the smaller distance with matches/mismatches and
    // the rest with indels.
    TDistanceThreshold distance = _max(getBeginDim0(b) - getEndDim0(a), getBeginDim1(b) - getEndDim1(a));
    // Compare distance with threshold.
    return distance <= maxGapSize;
}


// Updating the coordinates of seeds is the same for merging and
// simple chaining.  Only the score computation differs.
template <typename TSeedConfig>
inline void
_updateSeedsCoordinatesMergeOrSimpleChain(
        Seed<Simple, TSeedConfig> & seed,
        Seed<Simple, TSeedConfig> const & other)
{
    SEQAN_CHECKPOINT;

    setBeginDim0(seed, _min(getBeginDim0(seed), getBeginDim0(other)));
    setBeginDim1(seed, _min(getBeginDim1(seed), getBeginDim1(other)));
    setEndDim0(seed, _max(getEndDim0(seed), getEndDim0(other)));
    setEndDim1(seed, _max(getEndDim1(seed), getEndDim1(other)));
    setLowerDiagonal(seed, _min(getLowerDiagonal(seed), getLowerDiagonal(other)));
    setUpperDiagonal(seed, _max(getUpperDiagonal(seed), getUpperDiagonal(other)));
}


template <typename TSeedConfig, typename TScoreValue>
inline void
_combineSeeds(Seed<Simple, TSeedConfig> & seed,
              Seed<Simple, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & /*scoringScheme*/,
              Nothing const & /*sequence0*/,
              Nothing const & /*sequence1*/,
              Merge const &)
{
    SEQAN_CHECKPOINT;

    _updateSeedsScoreMerge(seed, other);
    _updateSeedsCoordinatesMergeOrSimpleChain(seed, other);
}


template <typename TSeedConfig, typename TScoreValue>
inline void
_combineSeeds(Seed<Simple, TSeedConfig> & seed,
              Seed<Simple, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & scoringScheme,
              Nothing const & /*sequence0*/,
              Nothing const & /*sequence1*/,
              SimpleChain const &)
{
    SEQAN_CHECKPOINT;

    typedef Seed<Simple, TSeedConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;

    _updateSeedsScoreSimpleChain(seed, other, scoringScheme);
    _updateSeedsCoordinatesMergeOrSimpleChain(seed, other);
}


template <typename TSeedConfig, typename TScoreValue, typename TSequence0, typename TSequence1>
inline void
_combineSeeds(Seed<Simple, TSeedConfig> & seed,
              Seed<Simple, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & scoringScheme,
              TSequence0 const & sequence0,
              TSequence1 const & sequence1,
              Chaos const &)
{
    SEQAN_CHECKPOINT;

    typedef Seed<Simple, TSeedConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;

    // TODO(holtgrew): Assert seed left of other.

    // Compute gaps in both dimensions, the remaining gap is the
    // vertical/horizontal distance we will not fill with CHAOS
    // chaining.
    //
    // TODO(holtgrew): We need + 1 here, do we need it anywhere else?
    TPosition gapDim0 = getBeginDim0(other) - getEndDim0(seed) + 1;
    TPosition gapDim1 = getBeginDim1(other) - getEndDim1(seed) + 1;
    TPosition minGap = _min(gapDim0, gapDim1);
    TPosition maxGap = _max(gapDim0, gapDim1);
    TPosition remainingGap = maxGap - minGap;

    // Compute new score using the CHAOS method.
    //
    // First, compute the score when force-aligning from seed.
    TPosition posLeft0 = getEndDim0(seed);
    TPosition posLeft1 = getEndDim1(seed);
    TScoreValue tmpScore = 0;
    for (TPosition i = 0; i < minGap; ++i)
        tmpScore += score(scoringScheme, sequence0[posLeft0 + i], sequence1[posLeft1 + i]);

    SEQAN_ASSERT_GT(getBeginDim0(other), static_cast<TPosition>(0));
    SEQAN_ASSERT_GT(getBeginDim1(other), static_cast<TPosition>(0));
    TPosition posRight0 = getBeginDim0(other);
    TPosition posRight1 = getBeginDim1(other);

    // Now, try to put the gap at each position and get the position
    // with the highest score.  If there are two such positions, the
    // first one found is returned which is the one that is furthest
    // away from seed.
    TPosition bestGapPos = 0;  // delta to lowermost position
    TScoreValue bestScore = tmpScore;
    for (TPosition i = 1; i < minGap; ++i) {
        tmpScore -= score(scoringScheme, sequence0[posLeft0 + minGap - 1 - i], sequence1[posLeft1 + minGap - 1 - i]);
        tmpScore += score(scoringScheme, sequence0[posRight0 - i], sequence1[posRight1 - i]);
        if (tmpScore > bestScore) {
            // Found a better score.
            bestScore = tmpScore;
            bestGapPos = i;
        }
    }

    // Now, the best gap is when extending the lower right seed
    // (other) by bestGapPos to the upper right.  However, this is
    // ignored for simple seeds: We simply update the score and are
    // done.
    _updateSeedsScoreChaos(seed, other, bestScore + remainingGap * scoreGap(scoringScheme));

    // For simple seeds, the coordinate computation is the same as for
    // merge/simple chain.
    //
    // TODO(holtgrew): Adjust the name of updateSeedsCoordinatesMergeOrSimpleChain to reflect this.
    _updateSeedsCoordinatesMergeOrSimpleChain(seed, other);
}


template <typename TSeedConfig, typename TScoreValue>
inline void
_combineSeeds(Seed<ChainedSeed, TSeedConfig> & seed,
              Seed<ChainedSeed, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & /*scoringScheme*/,
              Nothing const & /*sequence0*/,
              Nothing const & /*sequence1*/,
              Merge const &)
{
    SEQAN_CHECKPOINT;
    // For chained seeds, we first remove all diagonals from seed
    // until the last diagonal of seed starts truly before other.
    // Then, we possibly shorten the last diagonal.  Finally, we copy
    // over all diagonals from other.

    _updateSeedsScoreMerge(seed, other);

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


template <typename TSeedConfig, typename TScoreValue>
inline void
_combineSeeds(Seed<ChainedSeed, TSeedConfig> & seed,
              Seed<ChainedSeed, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & scoringScheme,
              Nothing const & /*sequence0*/,
              Nothing const & /*sequence1*/,
              SimpleChain const &)
{
    SEQAN_CHECKPOINT;
    // Simply copy over the diagonals of the seed (other) into the
    // left one (seed) after updating the score.

    _updateSeedsScoreSimpleChain(seed, other, scoringScheme);

    // Copy over other diagonals.
    typedef Seed<ChainedSeed, TSeedConfig> TSeed;
    typedef typename Iterator<TSeed const, Standard>::Type TConstIterator;
    for (TConstIterator it = begin(other, Standard()); it != end(other, Standard()); ++it)
        appendDiagonal(seed, *it);
}


template <typename TSeedConfig, typename TScoreValue, typename TSequence0, typename TSequence1>
inline void
_combineSeeds(Seed<ChainedSeed, TSeedConfig> & seed,
              Seed<ChainedSeed, TSeedConfig> const & other,
              Score<TScoreValue, Simple> const & scoringScheme,
              TSequence0 const & /*sequence0*/,
              TSequence1 const & /*sequence1*/,
              Chaos const &)
{
    SEQAN_CHECKPOINT;
    // Simply copy over the diagonals of the seed (other) into the
    // left one (seed) after updating the score.

    SEQAN_ASSERT_FAIL("Write me!");
}


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
            // seed is to be merged into *it.
            mergePartner = it;
            seedIsOnTheLeft = false;
            return true;
        } else if (_seedsCombineable(seed, *value(it), maxDistance, bandwidth, tag)) {
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
        if (_qualityReached(*value(it), seedSet, TQualityThreshold()))
            // TODO(holtgrew): Do not use dot-methods.
            seedSet._highQualitySeeds.insert(value(it));
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

