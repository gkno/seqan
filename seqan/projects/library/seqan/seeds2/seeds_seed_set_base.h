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


// struct _Scored;
// typedef Tag<_Scored> Scored;

// struct _UnScored;
// typedef Tag<_UnScored> UnScored;


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
template <typename TSeedSpec, typename TSpec, typename TSeedConfig = DefaultSeedConfig>
class SeedSet;

// ===========================================================================
// Metafunctions
// ===========================================================================

/**
.Metafunction.Position.param.T:Class.SeedSet
*/
/**
.Metafunction.Size.param.T:Class.SeedSet
*/
/**
.Metafunction.Value.param.T:Class.SeedSet
*/
/**
.Metafunction.GetValue.param.T:Class.SeedSet
*/
/**
.Metafunction.Reference.param.T:Class.SeedSet
*/
/**
.Metafunction.Iterator.param.T:Class.SeedSet
*/
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
.Function.back.param.object.type:Class.SeedSet
 */

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
