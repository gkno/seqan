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

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

}  // namespace seqan

#endif  // SEQAN_SEEDS_ALIGN_CHAIN_BANDED_ALIGNMENT_CHAIN_H_
