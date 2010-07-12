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
  Tags for the scoring scheme selection for SeedSet specializations.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SCORING_SCHEME_H_
#define SEQAN_SEEDS_SEEDS_SCORING_SCHEME_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Tag.Gap Costs
..cat:Seed Handling
..summary:The distance function to use for gap costs.
..see:Class.SeedSet
..see:Tag.Scoring Scheme
..tag.NoGapCosts:Gaps do not inflict any cost.
..tag.ManhattanDistance:The manhattan distance is used.
..tag.QueryDistance:Only gaps in the query sequence are evaluated.
..tag.DatabaseDistance:Only gaps in the database sequence are evaluated.
..include:seqan/seeds.h
*/

struct _NoGapCosts;
typedef Tag<_NoGapCosts> NoGapCosts;

struct _ManhattanDistance;
typedef Tag<_ManhattanDistance> ManhattanDistance;

struct _QueryDistance;
typedef Tag<_QueryDistance> QueryDistance;

struct _DatabaseDistance;
typedef Tag<_DatabaseDistance> DatabaseDistance;


/**
.Tag.Quality Factor
..cat:Seed Handling
..summary:Select quality factor evaluation.
..seed:Class.SeedSet
..see:Tag.Scoring Scheme
..tag.SeedScore:Use the seed score for the quality evaluation.
..tag.SeedLength:Use the seed length for the quality evaluation.
.include:seqan/seeds.h
 */
struct _SeedScore;
typedef Tag<_SeedScore> SeedScore;

struct _SeedLength;
typedef Tag<_SeedLength> SeedLength;


/**
.Tag.Scoring Scheme
..summary:Information about a seed set's scoring scheme.
..cat:Seed Handling
..signature:ScoringScheme<TGapCosts, TQualityFactor, TScore>
..param.TQualityFactor:Selects the seed evaluation function, for example for thresholds.
..param.TGapCosts:The cost function for gaps. NoGapCosts, Manhattan, QueryDistance, DatabaseDistance.
..param.TScore:The @Class.Score@ to use for mismatch scoring. void for no scores
..include:seqan/seeds.h
*/
template <typename TQualityFactor, typename TGapCosts, typename TScore>
struct ScoringScheme {};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SCORING_SCHEME_H_
