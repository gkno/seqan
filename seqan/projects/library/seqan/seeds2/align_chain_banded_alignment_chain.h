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
  This header defines the data structure alignment pool which supports
  the banded chain alignment algorithm.  It allows to store multiple
  alignments for each seed in the chain and each gap as well as the leading
  and trailing gap.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_ALIGNMENT_CHAIN_H_
#define SEQAN_SEEDS_ALIGN_CHAIN_BANDED_ALIGNMENT_CHAIN_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// Note that begin and end gaps are always free for alignment chains.
// For seeds, it does not make sense not to make them free and since
// we do not force the leading and trailing characters to align in the
// gapped X-drop extension, it makes even less sense.
template <typename TSegment, typename TScoringScheme, typename TAlignmentTag>
class _AlignmentChain
{
public:
    typedef typename Value<TScoringScheme>::Type TScoreValue;

    // The scoring scheme used for the alignment.
    Holder<TScoringScheme> scoringScheme_;

    // The first sequence / dimension 0 / query sequence, vertical down
    Holder<TSegment> sequence0_;

    // The second sequence / dimension 1 / database sequence, horizontal ltr
    Holder<TSegment> sequence1_;

    // The alignment matrices for the different part.  The first
    // element is the alignment matrix for the rectangle to the lower
    // right, followed by the one for the first seed in the chain
    // (lower-rightmost one), the next rectangle and so on.
    String<Matrix<TScoreValue, 2> > alignmentMatrices_;

    // TODO(holtgrew): Default constructor + setSequence{0,1} missing for now.

    _AlignmentChain(TScoringScheme scoringScheme, TSegment sequence0, TSegment sequence1)
            : scoringScheme_(scoringScheme), sequence0_(sequence0), sequence1_(sequence1)
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TSegment, typename TScoringScheme, typename TAlignmentTag>
struct ScoringScheme<_AlignmentChain<TSegment, TScoringScheme, TAlignmentTag> >
{
    typedef TScoringScheme Type;
};

template <typename TSegment, typename TScoringScheme, typename TAlignmentTag>
struct ScoringScheme<_AlignmentChain<TSegment, TScoringScheme, TAlignmentTag> const>
        : ScoringScheme<_AlignmentChain<TSegment, TScoringScheme, TAlignmentTag> > {};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignTrailingRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;

    typedef Seed<TSeedSpec, TSeedConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef typename Suffix<TSequence>::Type TSuffix;
    typedef Matrix<TScoreValue, 2> TMatrix;

    TSuffix segment0 = suffix(value(alignmentChain.sequence0_), getEndDim0(seed));
    TSuffix segment1 = suffix(value(alignmentChain.sequence1_), getEndDim1(seed));

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignTrailingRectangle" << std::endl;
    std::cout << "| segment0: '" << segment0 << "'" << std::endl;
    std::cout << "| segment1: '" << segment1 << "'" << std::endl;
    std::cout << "`--" << std::endl;

    // Append a matrix to the chain's matrices.  _needleman_wunsch()
    // below will resize the matrix appropriately.
    TMatrix tmpMatrix;
    appendValue(alignmentChain.alignmentMatrices_, tmpMatrix);
    _needleman_wunsch(back(alignmentChain.alignmentMatrices_), segment0, segment1, value(alignmentChain.scoringScheme_));
}


template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignLeadingRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;
}


template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & rightSeed,
        Seed<TSeedSpec, TSeedConfig> const & leftSeed)
{
    SEQAN_CHECKPOINT;
}

template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignSeed(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;

    // Copy values from rectangle alignment matrix to the lower right.

    // Perform banded alignment around the seed, using the values previously copied.
}

template <typename TAlignment, typename TSequence, typename TScoringScheme, typename TAlignmentTag>
typename Value<typename ScoringScheme<Alignment<TSequence, TScoringScheme, TAlignmentTag> >::Type >::Type
_glueAlignmentChain(
        TAlignment & alignment,
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> const & alignmentChain)
{
    SEQAN_CHECKPOINT;
    return 0;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_ALIGN_CHAIN_BANDED_ALIGNMENT_CHAIN_H_
