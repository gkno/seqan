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
  Author: Carsten Kemena <carsten.kemena@crg.es>
 ============================================================================
  Code for performing banded alignment around a seed.  The code is based on
  Carsten Kemena's original implementation and was adjusted to the new seeds
  and seed set interface.
  
  This file contains a facade function.  The real implementation is
  split into the headers align_seed_banded_{linear,affine}.h which
  have the implementations for the linear and affine gap costs.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_SEED_BANDED_H_
#define SEQAN_SEEDS_ALIGN_SEED_BANDED_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

/**
.Function.bandedAlignment:
..summary:Calculates a banded alignment around a Seed.
..cat:Seed Handling
..signature:bandedAlignmnet(align, seed, k, score)
..param.align:An alignment dataStructure containing two sequences (sections). The sequence (parts) must be equal to the part which is covered by the seed.
...type:Class.Align
..param.seed: The seed for wich the banded alignmnent shall be constructed.
...type:Class.Seed
..param.k: A value describing the additional extension of the diagonal band.
..param.score: The score matrix used.
...type:Spec.Simple Score
..returns:The score of the optimal banded alignment given in align.
..remarks:Use the function @Function.globalAlignment@ with the tag $BandedNeedlemanWunsch$ or $BandedGotoh$ for more general banded alignment without a seed.
..see:Function.globalAlignment
..see:Tag.Global Alignment Algorithms
..include:seqan/seeds.h
*/
template <typename TSource, typename TAlignSpec, typename TSeedSpec, typename TSeedConfig, typename TBandwidth, typename TScoreValue>
TScoreValue
bandedAlignment(
        Align<TSource, TAlignSpec> & alignment,
        Seed<TSeedSpec, TSeedConfig> const & seed,
        TBandwidth k,
        Score<TScoreValue, Simple> const & scoringScheme)
{
    SEQAN_CHECKPOINT;
    
    // TODO(holtgrew): For consistency with global alignment, there should probably be a version that accepts a set of two seeds? Could be more robust than having to give exactly the right infixes in the alignment.

    // Choose banded alignment around the given seed, depending on whether the
    // used scoring scheme is linear.
	if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
		return bandedAlignment(alignment, seed, k, scoringScheme, NeedlemanWunsch());
	else
		return bandedAlignment(alignment, seed, k, scoringScheme, Gotoh());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_SEED_BANDED_H_

