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
  This header contains the public interface to the banded chain
  alignment.  The implementations for linear and affine gap costs are
  in align_chain_banded_linear.h and align_chain_banded_affine.h.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_H_
#define SEQAN_SEEDS_ALIGN_CHAIN_BANDED_H_

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
.Function.bandedChainAlignment:
..summary:Calculates a banded alignment around a chain of seeds. 
..cat:Seed Handling
..signature:bandedChainAlignment(seedChain, k, alignment, scoreMatrix);
..param.seedChain:A chain of seeds, must be ascendingly sorted in both dimensions.
..param.k:Half of the width of the band.
..param.alignment:The alignment where the result is stored.
...type:Class.Align
..param.scoreMatrix: The score matrix.
...type:Spec.Simple Score
...remarks: Depending on the score matrix the Needleman-Wunsch or the Gotoh algorithm is used. For a description of the algorithm see the masters thesis of C. Kemena, Section 5.3.3 LAGAN Alignment.
..returns: The score of the alignment.
*/
// TODO(holtgrew): wholeAlignment is the result and should be the first parameter.
// TODO(holtgrew): Adjust the documentation to the parameter names.
template<typename TContainer, typename TValue, typename TScoreValue, typename TAlign, bool START1_FREE, bool START0_FREE, bool END1_FREE, bool END0_FREE>
TScoreValue
bandedChainAlignment(TContainer const & seedChain, 
					 TValue k,
					 TAlign & alignment, 
					 Score<TScoreValue, Simple> const & scoringScheme,
                     AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const & alignConfig)
{
    SEQAN_CHECKPOINT;

	if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
		return _bandedChainAlignment(alignment, seedChain, k, scoringScheme, alignConfig, NeedlemanWunsch());
	else
		return _bandedChainAlignment(alignment, seedChain, k, scoringScheme, alignConfig, Gotoh());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_H_
