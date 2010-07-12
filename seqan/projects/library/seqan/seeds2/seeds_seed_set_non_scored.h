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
  Non-Scored SeedSet specialization.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_SET_NON_SCORED_H_
#define SEQAN_SEEDS_SEEDS_SEED_SET_NON_SCORED_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Spec.Scored SeedSet
..summary:SeedSet that uses non-scored Seeds.
..cat:Seed Handling
..general:Class.SeedSet
..signature:SeedSet<TPosition, TSeedSpec, TScoringScheme>
..param.TPosition: Type that saves the positions and upper/lower bounds.
...remarks: Positive and negative values are needed.
..param.TSeedSpec:The @Class.Seed@ specialization.
..param.TScoringScheme:The scoring sheme to use.
..include:seqan/seeds.h
*/

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SET_NON_SCORED_H_
