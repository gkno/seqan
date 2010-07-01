/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  Basic definitions for the module random.
  ===========================================================================
*/

#ifndef SEQAN_RANDOM_RANDOM_BASE_H_
#define SEQAN_RANDOM_RANDOM_BASE_H_

namespace seqan {

/**
.Class.RNG:
..summary:Random Number Generator
..signature:RNG<TSpec>
..param.TSpec:Random Number Generator specialization.
*/
template <typename TSpec>
class RNG;

/**
.Class.PDF:
..summary:ProbabilityDensityFunction
..signature:PDF<TSpec>
..param.TSpec:Specialization.
*/
template <typename TSpec>
class PDF;

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_BASE_H_
