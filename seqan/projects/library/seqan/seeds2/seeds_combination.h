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
  Algorithms for combining (i.e. merging and chaining) seeds.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_COMBINATION_H_
#define SEQAN_SEEDS_SEEDS_COMBINATION_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Tag.Local Chaining
..cat:Seed Handling
..summary:The local chaining algorithms to use when adding a seed to a @Class.SeedSet@.
..see:Class.SeedSet
..see:Function.addSeed
..tag.Merge:Merge with existing seed.
..tag.Chaos:CHAOS chaining.
..tag.SimpleChain:Simple chaining.
..tag.Single:Add single seed without merging and chaining.
..include:seqan/seeds.h
*/
struct _Merge;
typedef Tag<_Merge> Merge;

struct _Chaos;
typedef Tag<_Chaos> Chaos;

struct _SimpleChain;
typedef Tag<_SimpleChain> SimpleChain;

struct _Single;
typedef Tag<_Single> Single;

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

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

    // b has to be right of a for the two seeds to be mergeable.
    if (getBeginDim0(b) < getBeginDim0(a) || getBeginDim1(b) < getBeginDim1(a))
        return false;
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
    TPosition gapDim0 = getBeginDim0(other) - getEndDim0(seed);
    TPosition gapDim1 = getBeginDim1(other) - getEndDim1(seed);
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

    // std::cout << "Merging chained seeds " << seed << " and " << other << std::endl;
    SEQAN_ASSERT_LEQ_MSG(getBeginDim0(seed), getBeginDim0(other), "Monotony in both dimensions required for merging.");
    SEQAN_ASSERT_LEQ_MSG(getBeginDim1(seed), getBeginDim1(other), "Monotony in both dimensions required for merging.");
    
    _updateSeedsScoreMerge(seed, other);

    // Remove diagonals.
    typedef Seed<ChainedSeed, TSeedConfig> TSeed;
    typedef typename Iterator<TSeed, Standard>::Type TIterator;
    TIterator it;
    // TODO(holtgrew): Could use back() instead of lastKept.
    TIterator lastKept = begin(seed);
    for (it = begin(seed); it != end(seed); ++it) {
        if (it->beginDim0 >= getBeginDim0(other) && it->beginDim1 >= getBeginDim1(other))
            break;
        lastKept = it;
    }
    if (it != end(seed))
        truncateDiagonals(seed, it);
    // std::cout << "Seed after truncating diagonals: " << seed << std::endl;

    // Shorten last diagonal if necessary.
    if (lastKept->beginDim0 + lastKept->length > getBeginDim0(other) && lastKept->beginDim1 + lastKept->length > getBeginDim1(other)) {
        lastKept->length = _min(getBeginDim0(other) - lastKept->beginDim0, getBeginDim1(other) - lastKept->beginDim1);
    } else if (lastKept->beginDim0 + lastKept->length > getBeginDim0(other)) {
        lastKept->length = getBeginDim0(other) - lastKept->beginDim0;
    } else if (lastKept->beginDim1 + lastKept->length > getBeginDim1(other)) {
        lastKept->length = getBeginDim1(other) - lastKept->beginDim1;
    }

    // Copy over other diagonals.
    typedef typename Iterator<TSeed const, Standard>::Type TConstIterator;
    for (TConstIterator it = begin(other, Standard()); it != end(other, Standard()); ++it)
        appendDiagonal(seed, *it);

    // std::cout << "Chained seed after merging: " << seed << std::endl;
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

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_SEEDS_COMBINATION_UNORDERED_H_
