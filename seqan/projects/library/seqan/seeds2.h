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

#include <algorithm>
#include <list>
#include <new>

#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/score.h>
#include <seqan/align.h>
#include <seqan/map.h>

// #ifdef SEQAN_SWITCH_USE_FORWARDS
// #include <seqan/seeds2/seeds2_generated_forwards.h>
// #endif

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

// Seed extension algorithms.
#include <seqan/seeds2/seeds_extension.h>

// Algorithms for chaining and merging seeds.
#include <seqan/seeds2/seeds_combination.h>

// Class SeedSet, specializations, iterators.
#include <seqan/seeds2/basic_iter_indirect.h>
#include <seqan/seeds2/seeds_seed_set_base.h>
#include <seqan/seeds2/seeds_seed_set_unordered.h>

// Banded and unbanded DP algorithms.
#include <seqan/seeds2/align_dynprog_linear.h>
#include <seqan/seeds2/align_dynprog_banded_linear.h>
#include <seqan/seeds2/align_dynprog_affine.h>
#include <seqan/seeds2/align_dynprog_banded_affine.h>

// Banded chain alignment.
#include <seqan/seeds2/align_chain_banded.h>

// Global chaining algorithms
// #include <seqan/seeds2/seeds_global_chaining_base.h>
// #include <seqan/seeds2/seeds_global_chaining_gusfield.h>
// #include <seqan/seeds2/seeds_global_chaining_lagan.h>

#endif  // #ifndef SEQAN_HEADER_SEEDS_H
