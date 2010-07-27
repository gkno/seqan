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
  This header defines the data structure AlignmentChain which supports
  the banded chain alignment algorithm.  It allows to store multiple
  alignments for each seed in the chain and each gap as well as the leading
  and trailing gap.

  The data structure supports the consecutive alignment of seeds and
  the rectangles between them.  The functions encapsulate the
  computation of dimensions for the rectangles and seeds.  They also
  encapsulate the dimension computation for parts to be copied.

  The alignment and especially the copying of the data itself is
  delegated to the alignment algorithms in
  align_dynprog[_banded]_{linear,affine}.h.

  The functions should not have to be adapted to accomodate less
  copying since they are not concerned with the copying themselves.
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
    // element is the alignment matrix for the rectangle to the upper
    // left right, followed by the one for the first seed in the chain
    // (left-uppermost one), the next rectangle and so on.
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
_alignLeadingRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Prefix<TSequence>::Type TPrefix;
    typedef typename Size<TSequence>::Type TSize;

    // Compute the overlaps of the sequence prefixes with the seed
    // banded alignment matrix.
    //
    // TODO(holtgrew): Are these computations correct?
    TSize rightOverlap0 = getLowerDiagonal(seed) - getStartDiagonal(seed) + alignmentChain.bandwidth_;
    TSize rightOverlap1 = getUpperDiagonal(seed) - getStartDiagonal(seed) + alignmentChain.bandwidth_;
    // The rectangle must overlap with the band around the seed (given
    // the bandwidth is > 0).
    SEQAN_ASSERT_GT(rightOverlap0, 0u);
    SEQAN_ASSERT_GT(rightOverlap1, 0u);
    // Limit overlap so it does not go over the alignment matrix
    // (given the bandwidth is > 0).
    rightOverlap0 = _min(rightOverlap0, length(value(alignmentChain.sequence0_)) - getBeginDim0(seed));
    rightOverlap1 = _min(rightOverlap1, length(value(alignmentChain.sequence1_)) - getBeginDim1(seed));
    // Get prefixes of the sequences to be aligned for the first
    // rectangle.  They use the overlap computed above.
    TPrefix prefix0 = prefix(value(alignmentChain.sequence0_), getBeginDim0(seed) + rightOverlap0);
    TPrefix prefix1 = prefix(value(alignmentChain.sequence1_), getBeginDim1(seed) + rightOverlap1);

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignLeadingRectangle" << std::endl;
    std::cout << "| prefix0: '" << prefix0 << "' (overlap = " << rightOverlap0 << ")" << std::endl;
    std::cout << "| prefix1: '" << prefix1 << "' (overlap = " << rightOverlap1 << ")" << std::endl;
    std::cout << "`--" << std::endl;

    // Append a new alignment matrix for the DP programming.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Compute DP matrix with free begin gaps in both sequences but no
    // free end gaps.
    _align_dynProg(back(alignmentChain.alignmentMatrices_), prefix0, prefix1, alignmentChain.scoringScheme_, AlignConfig<true, true, false, false>(), TAlignmentTag());

    // TODO(holtgrew): Temporary debug code.
    TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    std::cout << ",-- First rectangle, NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
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
_alignTrailingRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Suffix<TSequence>::Type TSuffix;
    typedef typename Size<TSequence>::Type TSize;

    // Compute the overlaps of the sequence suffixes with the seed
    // banded alignment matrix.
    //
    // TODO(holtgrew): Are these computations correct?
    TSize leftOverlap0 = getUpperDiagonal(seed) - getEndDiagonal(seed) + alignmentChain.bandwidth_;
    TSize leftOverlap1 = getEndDiagonal(seed) - getLowerDiagonal(seed) + alignmentChain.bandwidth_;
    // The rectangle must overlap with the band around the seed (given
    // the bandwidth is > 0).
    SEQAN_ASSERT_GT(leftOverlap0, 0u);
    SEQAN_ASSERT_GT(leftOverlap1, 0u);
    // Limit overlaps so they do not go beyond the alignment matrix.
    leftOverlap0 = _min(leftOverlap0, getEndDim0(seed));
    leftOverlap1 = _min(leftOverlap1, getEndDim1(seed));
    // Get suffixes of the sequences to be aligned for the last
    // rectangle.  They use the overlap computed above.
    TSuffix suffix0 = suffix(value(alignmentChain.sequence0_), getEndDim0(seed) - leftOverlap0);
    TSuffix suffix1 = suffix(value(alignmentChain.sequence1_), getEndDim1(seed) - leftOverlap1);

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignTrailingRectangle" << std::endl;
    std::cout << "| suffix0: '" << suffix0 << "' (overlap = " << leftOverlap0 << ")" << std::endl;
    std::cout << "| suffix1: '" << suffix1 << "' (overlap = " << leftOverlap1 << ")" << std::endl;
    std::cout << "`--" << std::endl;

    // Get a description of the lower horizontal stripe of the
    // previous banded DP matrix to copy into the rectangular DP
    // programming matrix.
    _DPMatrixRectangle<TScoreValue> bandedDPMatrixRectangle(back(alignmentChain.alignmentMatrices_), leftOverlap0, leftOverlap1);

    // Append a new alignment matrix for the DP programming.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Compute DP matrix copying in the values from the lower stripe
    // of a banded DP alignment matrix.
    _align_dynProg(back(alignmentChain.alignmentMatrices_), suffix0, suffix1, alignmentChain.scoringScheme_, bandedDPMatrixRectangle, TAlignmentTag());

    // TODO(holtgrew): Temporary debug code.
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
_alignRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & leftSeed,
        Seed<TSeedSpec, TSeedConfig> const & rightSeed)
{
    SEQAN_CHECKPOINT;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename Size<TSequence>::Type TSize;

    // Compute overlap with left seed in both dimensions.
    //
    // TODO(holtgrew): Are these computations correct?
    TSize leftOverlap0 = getUpperDiagonal(leftSeed) - getEndDiagonal(leftSeed) + alignmentChain.bandwidth_;
    TSize leftOverlap1 = getEndDiagonal(leftSeed) - getLowerDiagonal(leftSeed) + alignmentChain.bandwidth_;
    // The rectangle must overlap with the band around the seed.
    SEQAN_ASSERT_GEQ(leftOverlap0, 0u);
    SEQAN_ASSERT_GEQ(leftOverlap1, 0u);
    // TODO(holtgrew): Are these computations correct?
    TSize rightOverlap0 = getLowerDiagonal(rightSeed) - getStartDiagonal(rightSeed) + alignmentChain.bandwidth_;
    TSize rightOverlap1 = getUpperDiagonal(rightSeed) - getStartDiagonal(rightSeed) + alignmentChain.bandwidth_;
    // The rectangle must overlap with the band around the seed.
    SEQAN_ASSERT_GEQ(rightOverlap0, 0u);
    SEQAN_ASSERT_GEQ(rightOverlap1, 0u);

    // Limit overlaps so they do not go beyond the alignment matrix.
    leftOverlap0 = _min(leftOverlap0, getEndDim0(leftSeed));
    leftOverlap1 = _min(leftOverlap1, getEndDim1(leftSeed));
    rightOverlap0 = _min(rightOverlap0, length(value(alignmentChain.sequence0_)) - getBeginDim0(rightSeed));
    rightOverlap1 = _min(rightOverlap1, length(value(alignmentChain.sequence1_)) - getBeginDim1(rightSeed));

    // Get infixes of the sequences to be aligned for the rectangle in
    // the middle of the chain.  They use the overlaps computed above.
    TInfix infix0 = infix(value(alignmentChain.sequence0_), getEndDim0(leftSeed) - leftOverlap0, getBeginDim0(rightSeed) + rightOverlap0);
    TInfix infix1 = infix(value(alignmentChain.sequence1_), getEndDim1(leftSeed) - leftOverlap1, getBeginDim1(rightSeed) + rightOverlap1);

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignRectangle" << std::endl;
    std::cout << "| infix0: '" << infix0 << "' overlaps = (" << leftOverlap0 << ", " << rightOverlap0 << ")"<< std::endl;
    std::cout << "| infix1: '" << infix1 << "' overlaps = (" << leftOverlap1 << ", " << rightOverlap1 << ")" << std::endl;
    std::cout << "`--" << std::endl;

    // Get a description of the lower horizontal stripe of the
    // previous banded DP matrix to copy into the rectangular DP
    // programming matrix.  It is the same as the overlap with this DP
    // matrix.
    _DPMatrixRectangle<TScoreValue> bandedDPMatrixRectangle(back(alignmentChain.alignmentMatrices_), leftOverlap0, leftOverlap1);

    // Append a new alignment matrix for the DP programming.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Compute DP matrix copying in the values from the lower stripe
    // of a banded DP alignment matrix.
    _align_dynProg(back(alignmentChain.alignmentMatrices_), infix0, infix1, alignmentChain.scoringScheme_, bandedDPMatrixRectangle, TAlignmentTag());

    // TODO(holtgrew): Temporary debug code.
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
_alignSeed(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & seed)
{
    SEQAN_CHECKPOINT;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename Size<TSequence>::Type TSize;

    // Get the infixes corresponding to the seed.
    TInfix infix0 = infix(value(alignmentChain.sequence0_), getBeginDim0(seed), getEndDim0(seed));
    TInfix infix1 = infix(value(alignmentChain.sequence1_), getBeginDim1(seed), getEndDim1(seed));

    // Compute the overlapping part of the seed's banded alignment
    // matrix and the rectangle left of it.
    //
    // TODO(holtgrew): Are these computations correct?
    TSize leftOverlap0 = getStartDiagonal(seed) - getLowerDiagonal(seed) + alignmentChain.bandwidth_;
    TSize leftOverlap1 = getUpperDiagonal(seed) - getStartDiagonal(seed) + alignmentChain.bandwidth_;
    // The overlaps should be greater than 0 (given the bandwidth is > 0).
    SEQAN_ASSERT_GT(leftOverlap0, 0u);
    SEQAN_ASSERT_GT(leftOverlap1, 0u);
    // Get the DP matrix rectangle description corresponding to this overlap.
    _DPMatrixRectangle<TScoreValue> bandedDPMatrixRectangle(back(alignmentChain.alignmentMatrices_), leftOverlap0, leftOverlap1);

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignSeed" << std::endl;
    std::cout << "| infix0: '" << infix0 << "'" << std::endl;
    std::cout << "| infix1: '" << infix1 << "'" << std::endl;
    std::cout << "`--" << std::endl;

    // _align_banded_dynProg // XXX TODO(holtgrew): Actually call banded dynamic programming.

    // TODO(holtgrew): Temporary debug code.
    TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    std::cout << ",-- Seed Alignment, Banded NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
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
    SEQAN_ASSERT_FAIL("Write _glueAlignmentChain()!");
    return 0;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_ALIGN_CHAIN_BANDED_ALIGNMENT_CHAIN_H_
