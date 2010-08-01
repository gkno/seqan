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
  Author: Birte Kehr <bkehr@fu-berlin.de>
  Author: Carsten Kemena <carsten.kemena@crg.es>
 ============================================================================
  Banded chain alignment around a chain of seeds; Case with linear gap costs.
  Based on the original code by Carsten Kemena, adapted to the new seeds
  interface.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_
#define SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_

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

template <typename TContainer, typename TBandwidth, typename TScoreValue, typename TAlign, bool START1_FREE, bool START0_FREE, bool END1_FREE, bool END0_FREE, typename TGlobalAlignmentTag>
TScoreValue
_bandedChainAlignment(
        TAlign & alignment,
        TContainer const & seedChain,
        TBandwidth k,
        Score<TScoreValue, Simple> const & scoringScheme,
        AlignConfig<START1_FREE, START0_FREE, END1_FREE, END0_FREE> const & alignConfig,
        TGlobalAlignmentTag const &)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(seedChain), 0u);

    //
    // Function-Wide Typedefs
    //
    typedef typename Source<TAlign>::Type TSequence;
    typedef typename Infix<TSequence>::Type TSegment;
    typedef typename Value<TContainer>::Type TSeed;
    typedef typename Iterator<TContainer, Standard>::Type TIterator;
    typedef Score<TScoreValue, Simple> TScoringScheme;
    typedef _AlignmentChain<TSegment, TScoringScheme, TGlobalAlignmentTag> TAlignmentChain;

    //
    // Initialization
    //
	TSegment seq1 = sourceSegment(row(alignment, 0));
	TSegment seq2 = sourceSegment(row(alignment, 1));
    TAlignmentChain alignmentChain(k, scoringScheme, seq1, seq2);

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _bandedChainAlignment" << std::endl;
    std::cout << "| seq1:" << seq1 << std::endl;
    std::cout << "| seq2:" << seq2 << std::endl;
    std::cout << "`--" << std::endl;

    //
    // Compute Alignment using Alignment Chain
    //
    // Compute alignment for leading rectangle and the first seed.
    _alignLeadingRectangle(alignmentChain, front(seedChain), alignConfig);
    _alignSeed(alignmentChain, front(seedChain));

    // For all seeds from the second one from the left to the
    // rightmost one: Align rectangle left of it and then the seed.
    for (TIterator it = begin(seedChain) + 1; it != end(seedChain); ++it) {
        _alignRectangle(alignmentChain, value(it - 1), value(it));
        _alignSeed(alignmentChain, value(it));
    }

    // Compute alignment for the trailing rectangle.
    _alignTrailingRectangle(alignmentChain, back(seedChain));

    // Glue all alignments together.
    return _glueAlignmentChain(alignment, alignmentChain, seedChain, alignConfig);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_
