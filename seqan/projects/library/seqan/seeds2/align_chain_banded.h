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
..param.seedChain:A chain of seeds.
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
template<typename TContainer, typename TValue, typename TScoreValue, typename TAlign>
TScoreValue
bandedChainAlignment(TContainer const & seedChain, 
					 TValue k,
					 TAlign & alignment, 
					 Score<TScoreValue, Simple> const & scoringScheme)
{
    SEQAN_CHECKPOINT;

    // The assumption is that seedChain is sorted by dimension 0,
    // either ascendingly or descendingly.  In the first case, we have
    // to reverse the order since the banded chain alignment assumes
    // reverse order.  Because of the assumption about being sorted,
    // we only test the first two entries.
    if (length(seedChain) >= 2u && getBeginDim0(seedChain[0]) < getBeginDim0(seedChain[1])) {
        // TODO(holtgrew): If we have to copy, why not use a reverse string modifier so it is reversed on the fly?
        TContainer copyChain(seedChain);
        std::reverse(begin(copyChain), end(copyChain));

        // TODO(holtgrew): For consistency, we should have bandedChainAlignment(..., {Gotoh, NeedlemanWunsch}) functions.

        if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
            return _bandedChainAlignment_NW(alignment, copyChain, k, scoringScheme);
        else
            return _bandedChainAlignment_Gotoh(copyChain, k, alignment, scoringScheme);
    }
    // Chain of seeds is properly sorted, simply kick of the banded
    // chain alignment.
	if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
		return _bandedChainAlignment_NW(alignment, seedChain, k, scoringScheme);
	else
		return _bandedChainAlignment_Gotoh(seedChain, k, alignment, scoringScheme);
}
    
}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_H_
