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
  The class SeedSet.  ScoringScheme and related tags are defined in
  seeds_scoring_scheme.h
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_
#define SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

struct _MinSeedSize;
typedef Tag<_MinSeedSize> MinSeedSize;

struct _MinScore;
typedef Tag<_MinScore> MinScore;

// TODO(holtgrew): Put configs and mixins into their own header.
// TODO(holtgrew): Add metafunctions for mixin classes?

struct DefaultSeedSetConfig
{
    typedef DefaultSeedConfig TSeedConfig;
    typedef Nothing TQualityThreshold;
    typedef Nothing TQualityThresholdMixin;
};

// TODO(holtgrew): Rename _qualityReached to _qualityAboveThreshold
// Quality is always reached for seed sets without a quality threshold.
template <typename TSeed, typename TSeedSet>
bool _qualityReached(TSeed const & /*seed*/, TSeedSet const & /*seedSet*/, Nothing const &)
{
    SEQAN_CHECKPOINT;
    return true;
}

template <typename TSize>
struct _MinSeedSizeMixin
{
    TSize _minSeedSizeThreshold;

    _MinSeedSizeMixin() : _minSeedSizeThreshold(InfimumValue<TSize>::VALUE) {}
};

template <typename TSeedSet, typename TSize>
void setMinSeedSizeThreshold(TSeedSet & seedSet, TSize const & size)
{
    SEQAN_CHECKPOINT;
    seedSet._minSeedSizeThreshold = size;
}

template <typename TSize>
TSize getMinSeedSizeThreshold(_MinSeedSizeMixin<TSize> const & mixin)
{
    SEQAN_CHECKPOINT;
    return mixin._minSeedSizeThreshold;
}

template <typename TSeed, typename TSeedSet>
bool _qualityReached(TSeed const & seed, TSeedSet const & seedSet, MinSeedSize const &)
{
    SEQAN_CHECKPOINT;
    return getSeedSize(seed) >= getMinSeedSizeThreshold(seedSet);
}

struct DefaultSeedSetConfigLength
{
    typedef DefaultSeedConfig TSeedConfig;
    typedef MinSeedSize TQualityThreshold;
    typedef _MinSeedSizeMixin<TSeedConfig::TSize> TQualityThresholdMixin;
};

template <typename TScore>
struct _MinScoreMixin
{
    TScore _minScoreThreshold;
    _MinScoreMixin() : _minScoreThreshold(InfimumValue<TScore>::VALUE) {}
};

template <typename TSeedSet, typename TScore>
void setMinScoreThreshold(TSeedSet & seedSet, TScore const & score)
{
    SEQAN_CHECKPOINT;
    seedSet._minScoreThreshold = score;
}

template <typename TScore>
TScore getMinScoreThreshold(_MinScoreMixin<TScore> const & mixin)
{
    SEQAN_CHECKPOINT;
    return mixin._minScoreThreshold;
}


template <typename TSeed, typename TSeedSet>
bool _qualityReached(TSeed const & seed, TSeedSet const & seedSet, MinScore const &)
{
    SEQAN_CHECKPOINT;
    return getScore(seed) >= getMinScoreThreshold(seedSet);
}


struct DefaultSeedSetConfigScore
{
    typedef DefaultSeedConfigScore TSeedConfig;
    typedef MinScore TQualityThreshold;
    typedef _MinScoreMixin<TSeedConfig::TScoreValue> TQualityThresholdMixin;
};


/**
.Class.SeedSet:
..summary:Handles a set of seeds with local chaining on adding seeds.
..cat:Seed Handling
..signature:SeedSet<TSeedSpec, TScored, TSpec[, TSeedConfig]>
..param.TSeedSpec:Specialization of the seed to use.
..param.TScored:Either UnScored or a seed set scoring scheme specification.
..param.TSpec:Specialization of the seed set.
..param.TSeedConfig:Configuration for the seeds.  Sensible defaults are chosen based on the other template parameters.
..include:seqan/seeds.h
*/
template <typename TSeedSpec, typename TSpec, typename TSeedSetConfig = DefaultSeedSetConfig>
class SeedSet;

// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.Position.param.T:Class.SeedSet
///.Metafunction.Size.param.T:Class.SeedSet
///.Metafunction.Value.param.T:Class.SeedSet
///.Metafunction.GetValue.param.T:Class.SeedSet
///.Metafunction.Reference.param.T:Class.SeedSet
///.Metafunction.Iterator.param.T:Class.SeedSet
/**
.Metafunction.Seed:
..summary:Return the type of the underlying seed.
..signature:Seed<TSeedSet>::Value
..param.TSeedSet:The SeedSet to retrieve the seed type from.
..include:seqan/seeds.h
*/
/**
.Metafunction.ScoringScheme:
..summary:Return the type of the seed set's scorings cheme.
..signature:ScoringScheme<TSeedSet>::Value
..param.TSeedSet:The SeedSet to retrieve the scoring scheme type from.
..include:seqan/seeds.h
*/

// ===========================================================================
// Functions
// ===========================================================================

// Basic Container Functions

/**
.Function.begin.param.object.type:Class.SeedSet
.Function.end.param.object.type:Class.SeedSet
.Function.length.param.object.type:Class.SeedSet
.Function.front.param.object.type:Class.SeedSet
.Function.back.param.container.type:Class.SeedSet
..include:seqan/seeds2.h
 */
// TODO(holtgrew): dddoc {begin,end,length,front,back}All(T)

// SeedSet Functions

/**
.Function.addSeed:
..summary:Adds a seed to an existing set.
..cat:Seed Handling
..signature:addSeed(set, qPos, dPos, length, tag)
..signature:addSeed(set, qPos, dPos, qrPos, drPos, tag)
..signature:addSeed(set, seed, tag)
..param.set:The set to which the new seed should be added.
...type:Class.SeedSet
..param.qPos: Start position in sequence1.
..param.dPos: Start position in sequence2.
..param.length: Length of the seed.
..param.tag: The algorithm that should be used to add the new seed.
...type:Tag.Seed Adding
...remark: Note that not every algorithm can be used with each type of @Class.Seed@.
..param.qrPos: End position in sequence1.
..param.drPos: End Position in sequence2.
..param.seed: The new Seed.
...type:Class.Seed
...remarks: The seed is copied and then added.
..returns:Boolean if succesfully added.
...remarks:Always true for Tag Single.
..include:seqan/seeds2.h
*/

/**
.Function.addSeeds:
..summary:Adds several seeds to an existing set. If a merging or chaining algorithm is used seeds are added if the merging or chaining fails.
..cat:Seed Handling
..signature:addSeed(set, container, tag)
..signature:addSeed(set, begin, end, tag)
..param.set:The set to which the new seed sould be added.
...type:Class.SeedSet
..param.container: Content is copied to set.
...type:Concept.Container
..param.begin: Iterator pointing to the first value to add.
..param.end: Iterator pointing just behind the last value to add.
..param.tag: The algorithm that should be used to add the new @Class.Seed@.
...type:Tag.Local Chaining
...remarks: Note that not every algorithm can be used with each specialization of @Class.Seed@.
..include:seqan/seeds2.h
*/

// Debugging / TikZ Output

template <typename TStream, typename TQuerySequence, typename TDatabaseSequence, typename TSeedSetSpec, typename TSeedSpec, typename TSeedConfig>
inline void
_write(TStream & stream,
       TQuerySequence & sequence0,
       TDatabaseSequence & sequence1,
       SeedSet<TSeedSpec, TSeedSetSpec, TSeedConfig> const & seedSet,
       _Tikz const &)
{
    typedef SeedSet<TSeedSpec, TSeedSetSpec, TSeedConfig> TSeedSet;

    stream << "\\begin{tikzpicture}[" << std::endl
           << "    seed/.style={very thick}," << std::endl
           << "    seed diagonal/.style={red,<->}" << std::endl
           << "    ]" << std::endl;

    // Draw sequences.
    stream << "  \\draw";
    // Draw query / sequence 0;
    for (unsigned i = 0; i < length(sequence0); ++i)
        stream << std::endl << "    (0, -" << i << ") node {" << sequence0[i] << "}";
    stream << std::endl;
    // Draw database / sequence 1.
    for (unsigned i = 0; i < length(sequence1); ++i)
        stream << std::endl << "    (" << i << ", 0) node {" << sequence1[i] << "}";
    stream << ";" << std::endl;

    // Draw seeds.
    typedef typename Iterator<TSeedSet const, Standard>::Type TIterator;
    for (TIterator it = begin(seedSet); it != end(seedSet); ++it)
        _write(stream, value(it), _Tikz());
    stream << "\\end{tikzpicture}" << std::endl;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SET_BASE_H_
