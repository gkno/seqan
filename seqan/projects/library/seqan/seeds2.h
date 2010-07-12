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
  Module for two-dimensional seeding and chaining.
 ==========================================================================*/

#ifndef SEQAN_HEADER_SEEDS_H
#define SEQAN_HEADER_SEEDS_H

// ===========================================================================
// Preliminaries
// ===========================================================================
#include <list>

#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/score.h>
#include <seqan/align.h>
#include <seqan/map.h>

// ===========================================================================
// Seeds Module
// ===========================================================================

// Basic definitions
// #include <seqan/seeds2/seeds_base.h>

// Class Seed and specializations
#include <seqan/seeds2/seeds_seed_base.h>
#include <seqan/seeds2/seeds_seed_simple.h>
#include <seqan/seeds2/seeds_seed_diagonal.h>
#include <seqan/seeds2/seeds_seed_chained.h>

// Algorithms on Seed objects.
// #include <seqan/seeds2/seeds_extension.h>

// Class SeedSet and specializations
// #include <seqan/seeds2/seeds_scoring_scheme.h>
// #include <seqan/seeds2/seeds_seed_set_base.h>
// #include <seqan/seeds2/seeds_seed_set_scored.h>
// #include <seqan/seeds2/seeds_seed_set_non_scored.h>

// Algorithms on SeedSet objects
// #include <seqan/seeds2/seeds_add_seed.h>

// Global chaining algorithms
// #include <seqan/seeds2/seeds_global_chaining_base.h>
// #include <seqan/seeds2/seeds_global_chaining_gusfield.h>
// #include <seqan/seeds2/seeds_global_chaining_lagan.h>

#endif  // #ifndef SEQAN_HEADER_SEEDS_H
