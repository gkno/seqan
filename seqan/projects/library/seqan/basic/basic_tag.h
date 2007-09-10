 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_TAG_H
#define SEQAN_HEADER_BASIC_TAG_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Tag.DotDrawing
..summary:Switch to trigger drawing in dot format.
..value.DotDrawing:Graphs in dot format.
*/

struct DotDrawing_;
typedef Tag<DotDrawing_> const DotDrawing;


/**
.Tag.HammingDistance
..summary:Switch to trigger Hamming distance, which is a measure of character substitutions.
*/

/**
.Tag.LevenshteinDistance
..summary:Switch to trigger Levenshtein distance, which is a measure of edit operations (character substitutions, deletions or insertions).
*/

struct _HammingDistance;
struct _LevenshteinDistance;

typedef Tag<_HammingDistance>		HammingDistance;
typedef Tag<_LevenshteinDistance>	LevenshteinDistance;
typedef Tag<_LevenshteinDistance>	EditDistance;

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
