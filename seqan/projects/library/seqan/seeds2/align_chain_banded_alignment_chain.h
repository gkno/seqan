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

    _AlignmentChain(TPosition bandwidth, TScoringScheme const & scoringScheme, TSegment /*const*/ & sequence0, TSegment /*const*/ & sequence1)
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

// Given a seed and an alignment chain, the overlap of the seed with
// the rectangle to the lower right is written to overlap0 and
// overlap1.
template <typename TSize, typename TSequence, typename TScoringScheme, typename TSeedSpec, typename TSeedConfig, typename TAlignmentTag>
inline void
_computeLowerRightOverlap(TSize & overlap0,
                          TSize & overlap1,
                          Seed<TSeedSpec, TSeedConfig> const & seed,
                          _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> const & alignmentChain)
{
    SEQAN_CHECKPOINT;
    overlap0 = getUpperDiagonal(seed) - getEndDiagonal(seed) + alignmentChain.bandwidth_;
    overlap1 = getEndDiagonal(seed) - getLowerDiagonal(seed) + alignmentChain.bandwidth_;
    // Don't overlap more than the seed...
    overlap0 = _min(overlap0, getEndDim0(seed) - getBeginDim0(seed));
    overlap1 = _min(overlap1, getEndDim1(seed) - getBeginDim1(seed));
}


// Given a seed and an alignment chain, the overlap of the seed with
// the rectangle to the upper left is written to overlap0 and
// overlap1.
template <typename TSize, typename TSequence, typename TScoringScheme, typename TSeedSpec, typename TSeedConfig, typename TAlignmentTag>
inline void
_computeUpperLeftOverlap(TSize & overlap0,
                         TSize & overlap1,
                         Seed<TSeedSpec, TSeedConfig> const & seed,
                         _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> const & alignmentChain)
{
    SEQAN_CHECKPOINT;
    overlap0 = getStartDiagonal(seed) - getLowerDiagonal(seed) + alignmentChain.bandwidth_;
    overlap1 = getUpperDiagonal(seed) - getStartDiagonal(seed) + alignmentChain.bandwidth_;
    // Don't overlap more than the seed...
    overlap0 = _min(overlap0, getEndDim0(seed) - getBeginDim0(seed));
    overlap1 = _min(overlap1, getEndDim1(seed) - getBeginDim1(seed));
}


// Performs the alignment matrix filling step for the leading
// rectangle in the alignment chain.
template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignLeadingRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & rightSeed)
{
    SEQAN_CHECKPOINT;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename Size<TSequence>::Type TSize;

    // Compute overlap of the rectangle with the seed to the lower right.
    TSize rightOverlap0, rightOverlap1;
    _computeUpperLeftOverlap(rightOverlap0, rightOverlap1, rightSeed, alignmentChain);

    // Get infixes of the sequences that correspond to this rectangle.
    TInfix prefix0 = prefix(value(alignmentChain.sequence0_), getBeginDim0(rightSeed) + rightOverlap0);
    TInfix prefix1 = prefix(value(alignmentChain.sequence1_), getBeginDim1(rightSeed) + rightOverlap1);

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignLeadingRectangle" << std::endl;
    std::cout << "| prefix0: '" << prefix0 << "'" << std::endl;
    std::cout << "| prefix1: '" << prefix1 << "'" << std::endl;
    std::cout << "`--" << std::endl;

    // Append a new alignment matrix to the chain.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Resize the alignment matrix to the appropriate size.
    _align_resizeMatrix(back(alignmentChain.alignmentMatrices_), prefix0, prefix1, TAlignmentTag());
    // Initialize the matrix gutter.
    _align_initGutter(back(alignmentChain.alignmentMatrices_), alignmentChain.scoringScheme_, AlignConfig<true, true, false, false>(), TAlignmentTag());
    // Fill the Matrix using standard dynamic programming.
    _align_fillMatrix(back(alignmentChain.alignmentMatrices_), prefix0, prefix1, alignmentChain.scoringScheme_, TAlignmentTag());

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


/*
template <typename TSequence, typename TScoringScheme, typename TAlignmentTag, typename TSeedSpec, typename TSeedConfig>
void
_alignTrailingRectangle(
        _AlignmentChain<TSequence, TScoringScheme, TAlignmentTag> & alignmentChain,
        Seed<TSeedSpec, TSeedConfig> const & leftSeed)
{
    SEQAN_CHECKPOINT;

    typedef typename Value<TScoringScheme>::Type TScoreValue;
    typedef Matrix<TScoreValue, 2> TMatrix;
    typedef typename Suffix<TSequence>::Type TSuffix;
    typedef typename Size<TSequence>::Type TSize;

    // Compute overlap of the rectangle with the seed to the upper left.
    TSize leftOverlap0, leftOverlap1;
    _computeLowerRightOverlap(leftOverlap0, leftOverlap1, leftSeed, alignmentChain);

    // Get infixes of the sequences that correspond to this rectangle.
    TInfix suffix0 = suffix(value(alignmentChain.sequence0_), getEndDim0(leftSeed));
    TInfix suffix1 = suffix(value(alignmentChain.sequence1_), getEndDim1(leftSeed));

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignTrailingRectangle" << std::endl;
    std::cout << "| suffix0: '" << suffix0 << "'" << std::endl;
    std::cout << "| suffix1: '" << suffix1 << "'" << std::endl;
    std::cout << "`--" << std::endl;

    // Append a new alignment matrix to the chain.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Resize the alignment matrix to the appropriate size.
    _align_resizeMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, TAlignmentTag());
    // Copy over the data from the banded seed alignment matrix.
    _align_initGutterFromBanded(back(alignmentChain.alignmentMatrices_), value(end(alignmentChain.alignmentMatrices_) - 2), leftOverlap0, leftOverlap1, TAlignmentTag());
    // Initialize the rest of the gutter.
    _align_initGutter(back(alignmentChain.alignmentMatrices_), AlignConfig<true, true, false, false>(), leftOverlap0, leftOverlap1, TAlignmentTag());
    // Fill the Matrix using standard dynamic programming.
    _align_fillMatrix(back(alignmentChain.alignmentMatrices_), suffix0, suffix1, alignmentChain.scoringScheme_, TAlignmentTag());

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

    // Compute overlap of the rectangle with the seed to the upper left and to the lower right.
    TSize leftOverlap0, leftOverlap1, rightOverlap0, rightOverlap1;
    _computeLowerRightOverlap(leftOverlap0, leftOverlap1, leftSeed, alignmentChain);
    _computeUpperLeftOverlap(rightOverlap0, rightOverlap1, rightSeed, alignmentChain);

    // Get infixes of the sequences that correspond to this rectangle.
    TInfix infix0 = infix(value(alignmentChain.sequence0_), getEndDim0(leftSeed) - leftOverlap0, getBeginDim0(rightSeed) + rightOverlap0);
    TInfix infix1 = infix(value(alignmentChain.sequence1_), getEndDim1(leftSeed) - leftOverlap1, getBeginDim1(rightSeed) + rightOverlap1);

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignRectangle" << std::endl;
    std::cout << "| infix0: '" << infix0 << "'" << std::endl;
    std::cout << "| infix1: '" << infix1 << "'" << std::endl;
    std::cout << "`--" << std::endl;

    // Append a new alignment matrix to the chain.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Resize the alignment matrix to the appropriate size.
    _align_resizeMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, TAlignmentTag());
    // Copy over the data from the banded seed alignment matrix.
    _align_initGutterFromBanded(back(alignmentChain.alignmentMatrices_), value(end(alignmentChain.alignmentMatrices_) - 2), leftOverlap0, leftOverlap1, TAlignmentTag());
    // Initialize the rest of the gutter.
    _align_initGutter(back(alignmentChain.alignmentMatrices_), AlignConfig<false, false, false, false>(), leftOverlap0, leftOverlap1, TAlignmentTag());
    // Fill the Matrix using standard dynamic programming.
    _align_fillMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, alignmentChain.scoringScheme_, TAlignmentTag());

    // TODO(holtgrew): Temporary debug code.
    TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    std::cout << ",-- Middle rectangle, NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
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
    typedef Seed<TSeedSpec, TSeedConfig> TSeed;
    typedef typename Diagonal<TSeed>::Type TDiagonal;

    // Compute overlap with the rectangle to the upper left.
    TSize leftOverlap0, leftOverlap1, rightOverlap0, rightOverlap1;
    _computeUpperLeftOverlap(leftOverlap0, leftOverlap1, seed, alignmentChain);

    // Get the infixes corresponding to the seed.
    TInfix infix0 = infix(value(alignmentChain.sequence0_), getBeginDim0(seed), getEndDim0(seed));
    TInfix infix1 = infix(value(alignmentChain.sequence1_), getBeginDim1(seed), getEndDim1(seed));

    // TODO(holtgrew): Temporary debug code.
    std::cout << ",-- _alignSeed" << std::endl;
    std::cout << "| infix0: '" << infix0 << "'" << std::endl;
    std::cout << "| infix1: '" << infix1 << "'" << std::endl;
    std::cout << "`--" << std::endl;

    // Append a new alignment matrix to the chain.
    appendValue(alignmentChain.alignmentMatrices_, TMatrix());
    // Resize the alignment matrix to the appropriate size.
    _alignBanded_resizeMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, TAlignmentTag());
    // Copy over the data from the unbanded banded rectangle alignment matrix.
    _alignBanded_initGutterFromUnbanded(back(alignmentChain.alignmentMatrices_), value(end(alignmentChain.alignmentMatrices_) - 2), leftOverlap0, leftOverlap1, TAlignmentTag());
    // Initialize the rest of the gutter.
    _alignBanded_initGutter(back(alignmentChain.alignmentMatrices_), AlignConfig<false, false, false, false>(), leftOverlap0, leftOverlap1, TAlignmentTag());
    // Fill the Matrix using banded dynamic programming.
    _align_fillMatrix(back(alignmentChain.alignmentMatrices_), infix0, infix1, alignmentChain.scoringScheme_, TAlignmentTag());

    // TODO(holtgrew): Temporary debug code.
    TMatrix & matrix = back(alignmentChain.alignmentMatrices_);
    std::cout << ",-- Seed Alignment, Banded NW Matrix " << length(matrix, 0) << " x " << length(matrix, 1) << std::endl;
    for (unsigned i = 0; i < length(matrix, 0); ++i) {
        std::cout << "|\t";
        for (unsigned j = 0; j < i; ++j)
            std::cout << "\t";
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
*/
}  // namespace seqan

#endif  // SEQAN_SEEDS_ALIGN_CHAIN_BANDED_ALIGNMENT_CHAIN_H_
