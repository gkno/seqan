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

// TODO(holtgrew): Remove includes here, for debug only.
#include <iostream>
#include <fstream>

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

template <typename TStream, typename TAlignment, typename TAlignmentChain, typename TSeedChain>
void
write(TStream & stream, TAlignment const & alignment, TAlignmentChain const & alignmentChain, TSeedChain const & seedChain, _Tikz const &)
{
    typedef typename Row<TAlignment const>::Type TRow;
    typedef typename Position<TAlignment>::Type TPosition;
    typedef typename Iterator<TRow, Standard>::Type TRowIterator;
    typedef typename Iterator<TSeedChain, Standard>::Type TSeedIterator;
  
    stream << "\\documentclass{article}" << std::endl
           << "\\usepackage{tikz,nicefrac,amsmath,pifont}" << std::endl
           << "\\usetikzlibrary{arrows,snakes,backgrounds,patterns,matrix,shapes,fit,calc,shadows,plotmarks}" << std::endl
           << "\\begin{document}" << std::endl
           // << "\\begin{tikzpicture}[scale=.5,core/.style={fill=blue!20, fill opacity=.8},band/.style={fill=blue!10, fill opacity=.8}]" << std::endl;
           << "\\begin{tikzpicture}[scale=.5,core/.style={},band/.style={},rectangle/.style={fill=green!10, fill opacity=.1}]" << std::endl;

    // Grid.
    stream << "\\draw[help lines] (0, -" << length(row(alignment, 0)) << ") grid (" << length(row(alignment, 1)) << ", 0);" << std::endl;
    // Top row of characters.
    stream << "\\draw";
    int i = 0;
    for (TRowIterator it = begin(row(alignment, 0)); it != end(row(alignment, 0)); ++it, ++i) {
        stream << " (-1," << -0.5-i << ") node {" << *it << "}";
    }
    stream << ";" << std::endl;
    // Left row of characters.
    stream << "\\draw";
    i = 0;
    for (TRowIterator it = begin(row(alignment, 1)); it != end(row(alignment, 1)); ++it, ++i) {
        stream << " (" << 0.5+i << ", 1) node {" << *it << "}";
    }
    stream << ";" << std::endl;
    // Draw seeds.
    for (TSeedIterator it = begin(seedChain); it != end(seedChain); ++it) {
        int height = getEndDim0(*it) - getBeginDim0(*it);
        int width = getEndDim1(*it) - getBeginDim1(*it);
        int ext = _min(height, width);
        int s0 = getBeginDim0(*it);
        int s1 = getBeginDim1(*it);
        int e0 = getEndDim0(*it);
        int e1 = getEndDim1(*it);
        int ld = getLowerDiagonal(*it);
        int ud = getUpperDiagonal(*it);
        int sd = getStartDiagonal(*it);
        int ed = getEndDiagonal(*it);
        int b = alignmentChain.bandwidth_;
        // Draw seed "core."
        stream << "\\draw[core] (" << s1 << ", -" << s0 << ") -- (" << s1 + ext << ", -" << s0 + ext << ") -- (" << e1 << ", -" << e0 << ") -- (" << e1 - ext << ", -" << e0 - ext << ") -- cycle;" << std::endl;
        // Draw left band.
        stream << "\\draw[band] (" << s1 - ((ld - sd) + b) << ", -" << s0 << ") -- (" << s1 - ((ld - sd) + b) + ext << ", -" << s0 + ext << ") -- (" << s1 + ext << ", -" << s0 + ext << ") -- (" << s1 << ", -" << s0 << ") -- cycle;" << std::endl;
        // Draw right band.
        stream << "\\draw[band] (" << e1 << ", -" << e0 << ") -- (" << e1 + ((ud - ed) + b) << ", -" << e0 << ") -- (" << e1 + ((ud - ed) + b) - ext << ", -" << e0 - ext << ") -- (" << e1 - ext << ", -" << e0 - ext << ") -- cycle;" << std::endl;
    }

    // Draw rectangles.
    // First rectangle.
    {
        int s0 = getBeginDim0(front(seedChain));
        int s1 = getBeginDim1(front(seedChain));
        int o0, o1;
        _computeLowerRightOverlap(o0, o1, front(seedChain), alignmentChain);
        stream << "\\draw[rectangle] (0, 0) rectangle (" << s1 + o1 << ", -" << s0 + o0 << ");" << std::endl;
    }
    // Middle rectangles.
    if (length(seedChain) > 1u) {
        for (TSeedIterator it = begin(seedChain) + 1; it != end(seedChain); ++it) {
            int e0 = getEndDim0(value(it - 1));
            int e1 = getEndDim1(value(it - 1));
            int s0 = getBeginDim0(value(it));
            int s1 = getBeginDim1(value(it));
            int ou0, ou1, ol0, ol1;
            _computeLowerRightOverlap(ou0, ou1, value(it-1), alignmentChain);
            _computeUpperLeftOverlap(ol0, ol1, value(it), alignmentChain);
            printf("... e0=%d e1=%d s0=%d s1=%d ou0=%d ou1=%d ol0=%d ol1=%d\n", e0, e1, s0, s1, ou0, ou1, ol0, ol1);
            stream << "\\draw[rectangle] (" << e1 - ou1 << ", -" << e0 - ou0 << ") rectangle (" << s1 + ol1 << ", -" << s0 + ol0 << ");" << std::endl;
        }
    }
    // Last rectangle.
    {
        int e0 = getEndDim0(back(seedChain));
        int e1 = getEndDim1(back(seedChain));
        int o0, o1;
        _computeLowerRightOverlap(o0, o1, back(seedChain), alignmentChain);
        stream << "\\draw[rectangle] (" << e1 - o1 << ", -" << e0 - o0 << ") rectangle (" << length(row(alignment, 1)) << ", -" << length(row(alignment, 0)) << ");" << std::endl;
    }

    stream << "\\end{tikzpicture}" << std::endl
           << "\\end{document}" << std::endl;
}

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

    // // TODO(holtgrew): Temporary debug code.
    // std::cout << ",-- _bandedChainAlignment" << std::endl;
    // std::cout << "| seq1:" << seq1 << std::endl;
    // std::cout << "| seq2:" << seq2 << std::endl;
    // std::cout << "`--" << std::endl;

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

    // // Write alignment to file.
    // {
    //     std::ofstream file;
    //     file.open("/tmp/example.tex", std::ios::out);
    //     write(file, alignment, alignmentChain, seedChain, _Tikz());
    // }

    // Glue all alignments together.
    return _glueAlignmentChain(alignment, alignmentChain, seedChain, alignConfig);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_ALIGN_CHAIN_BANDED_LINEAR_H_
