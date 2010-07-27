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
    // TODO(holtgrew): Underscores should be in front of variable names, right?
    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef typename Position<TSegment>::Type TPosition;

    // The scoring scheme used for the alignment.
    TScoringScheme scoringScheme_;

    // The first sequence / dimension 0 / query sequence, vertical down
    Holder<TSegment> sequence0_;

    // The second sequence / dimension 1 / database sequence, horizontal ltr
    Holder<TSegment> sequence1_;

    // The bandwidth for the banded alignment.
    TPosition bandwidth_;

    // The alignment matrices for the different part.  The first
    // element is the alignment matrix for the rectangle to the lower
    // right, followed by the one for the first seed in the chain
    // (lower-rightmost one), the next rectangle and so on.
    String<Matrix<TScoreValue, 2> > alignmentMatrices_;

    // TODO(holtgrew): Default constructor + setSequence{0,1} missing for now.

    _AlignmentChain(TPosition bandwidth, TScoringScheme scoringScheme, TSegment sequence0, TSegment sequence1)
            : bandwidth_(bandwidth), scoringScheme_(scoringScheme), sequence0_(sequence0), sequence1_(sequence1)
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

    // Get dimensions of overlaps from the end of the seed in both dimensions.
    TPosition delta0 = getUpperDiagonal(seed) - getStartDiagonal(seed) + alignmentChain.bandwidth_;
    TPosition delta1 = getLowerDiagonal(seed) - getEndDiagonal(seed) + alignmentChain.bandwidth_;
    // Get suffixes to compute the alignment matrix for.
    TSuffix segment0 = suffix(value(alignmentChain.sequence0_), getEndDim0(seed) - delta0);
    TSuffix segment1 = suffix(value(alignmentChain.sequence1_), getEndDim1(seed) - delta1);

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignTrailingRectangle" << std::endl;
    std::cout << "| segment0: '" << segment0 << "'" << std::endl;
    std::cout << "| segment1: '" << segment1 << "'" << std::endl;
    std::cout << "`--" << std::endl;

    // Append a matrix to the chain's matrices.  _needleman_wunsch()
    // below will resize the matrix appropriately.
    TMatrix tmpMatrix;
    appendValue(alignmentChain.alignmentMatrices_, tmpMatrix);
    // TODO(holtgrew): Should be end gap free!
    _needleman_wunsch(back(alignmentChain.alignmentMatrices_), segment0, segment1, alignmentChain.scoringScheme_);

    TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    std::cout << ",-- Last rectangle, NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    for (unsigned i = 0; i < length(matrix, 0); ++i) {
        std::cout << "|\t";
        for (unsigned j = 0; j < length(matrix, 1); ++j) {
            std::cout << value(matrix, i, j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "`--" << std::endl;
}


template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignLeadingRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;

    typedef typename Prefix<TSequence>::Type TPrefix;

    TPrefix segment0 = prefix(value(alignmentChain.sequence0_), getBeginDim0(seed));
    TPrefix segment1 = prefix(value(alignmentChain.sequence1_), getBeginDim1(seed));

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignLeadingRectangle" << std::endl;
    std::cout << "| segment0: '" << segment0 << "'" << std::endl;
    std::cout << "| segment1: '" << segment1 << "'" << std::endl;
    std::cout << "`--" << std::endl;
}


template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & rightSeed,
        Seed<TSeedSpec, TSeedConfig> const & leftSeed)
{
    SEQAN_CHECKPOINT;

    typedef typename Infix<TSequence>::Type TInfix;

    TInfix segment0 = infix(value(alignmentChain.sequence0_), getEndDim0(leftSeed), getBeginDim1(rightSeed));
    TInfix segment1 = infix(value(alignmentChain.sequence1_), getEndDim1(leftSeed), getBeginDim1(rightSeed));
    
    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignRectangle" << std::endl;
    std::cout << "| segment0: '" << segment0 << "'" << std::endl;
    std::cout << "| segment1: '" << segment1 << "'" << std::endl;
    std::cout << "`--" << std::endl;
}

template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignSeed(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef typename Position<TSequence>::Type TPosition;
    typedef typename Infix<TSequence>::Type TInfix;
    typedef Matrix<TScoreValue, 2> TMatrix;
    // TODO(holtgrew): Const matrices broken.
    typedef typename Iterator<TMatrix /*const*/, Standard>::Type TConstMatrixIterator;

    TInfix segment0 = infix(value(alignmentChain.sequence0_), getBeginDim0(seed), getEndDim0(seed));
    TInfix segment1 = infix(value(alignmentChain.sequence1_), getBeginDim1(seed), getEndDim1(seed));
    
    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignSeed" << std::endl;
    std::cout << "| segment0: '" << segment0 << "'" << std::endl;
    std::cout << "| segment1: '" << segment1 << "'" << std::endl;
    std::cout << "`--" << std::endl;

    // Copy values from rectangle alignment matrix to the lower right.
    //
    // Compute overlapping sizes.
    TPosition delta0 = getUpperDiagonal(seed) - getStartDiagonal(seed) + alignmentChain.bandwidth_;
    TPosition delta1 = getLowerDiagonal(seed) - getEndDiagonal(seed) + alignmentChain.bandwidth_;
    // Allocate memory for copying.
    String<TScoreValue> scoreValues;
    reserve(scoreValues, delta0 + delta1 + 1);
    // In the banded seed alignment, the values are first copied in
    // from top to bottom along the vertical edge of the empty
    // triangle.  Then, the bottom border is filled from right to
    // left.  We copy out the values from the matrix in this order.
    {
        std::cout << "delta0 = " << delta0 << ", delta1 = " << delta1 << std::endl;
        std::cout << "copied out scores: ";
        SEQAN_ASSERT_GT(length(alignmentChain.alignmentMatrices_), 0u);
        TMatrix /*const*/ & matrix = back(alignmentChain.alignmentMatrices_);
        TConstMatrixIterator it = begin(matrix);
        setPosition(it, delta1);
        for (TPosition i = 0; i < delta0; ++i) {
            std::cout << *it << " ";
            appendValue(scoreValues, *it);
            goNext(it, 0);
        }
        std::cout << *it << " ";
        appendValue(scoreValues, *it);
        for (TPosition i = 0; i < delta1; ++i) {
            goPrevious(it, 1);
            std::cout << *it << " ";
            appendValue(scoreValues, *it);
        }
        std::cout << std::endl;
    }

    // Perform banded alignment around the seed, using the values previously copied.
    {
        TMatrix matrix;
        appendValue(alignmentChain.alignmentMatrices_, matrix);
        _bandedAlignment_NW_align(back(alignmentChain.alignmentMatrices_),
                                  seed,
                                  alignmentChain.bandwidth_,
                                  segment0,
                                  segment1,
                                  alignmentChain.scoringScheme_,
                                  scoreValues);
    }
    TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    std::cout << ",-- Banded Seed NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    for (unsigned i = 0; i < length(matrix, 0); ++i) {
        std::cout << "|\t";
        for (unsigned j = 0; j < length(matrix, 1); ++j) {
            std::cout << value(matrix, i, j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "`--" << std::endl;
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
