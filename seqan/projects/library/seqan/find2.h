/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007-2010

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
  Umbrella header for the find module.
 ==========================================================================*/

#ifndef SEQAN_SEQAN_FIND_H_
#define SEQAN_SEQAN_FIND_H_

// Prerequisites.
#include <cmath>

#include <deque>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
//#include <seqan/modifier.h>
//#include <seqan/score.h>
//#include <seqan/graph_types.h>
//#include <seqan/graph_algorithms.h>
//#include <seqan/map.h>
//#include <seqan/find.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/find2/find2_generated_forwards.h>
#endif

#include <seqan/find2/find_base.h>
#include <seqan/find2/find_finder_default.h>

// Exact pattern matching.
#include <seqan/find2/find_exact_simple.h>

// Complex pattern matching.
#include <seqan/find2/find_pattern_wild_shiftand.h>

// Multiple exact pattern search.
#include <seqan/find2/find_multiple_exact_shiftand.h>

// Approximate pattern matching with Hamming distance.
#include <seqan/find2/find_hamming_simple.h>

// Approximate matching with linear/affine gap costs, edit distance etc.
#include <seqan/find2/find_approx_dpsearch.h>

#endif  // SEQAN_SEQAN_FIND_H_
