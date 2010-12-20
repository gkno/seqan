 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

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
  Umbrella header for the random module.
 ==========================================================================*/

#ifndef SEQAN_RANDOM_H_
#define SEQAN_RANDOM_H_

//____________________________________________________________________________
// Prerequisites

#include <cmath>
#include <cstdlib>

#include <seqan/basic.h>
#include <seqan/misc/misc_random.h>

//____________________________________________________________________________
// Module Headers

// Basic Definitions
#include <seqan/random/random_base.h>

// Random Number Generation
#include <seqan/random/random_mt19937.h>

// RNG With Special Distributions.
#include <seqan/random/random_uniform.h>
#include <seqan/random/random_normal.h>
#include <seqan/random/random_lognormal.h>  // uses normal.h, uniform.h
#include <seqan/random/random_geometric.h>
#include <seqan/random/random_rng_functor.h>

// Functions with randomness.
#include <seqan/random/random_shuffle.h>

//____________________________________________________________________________

#endif  // SEQAN_RANDOM_H_
