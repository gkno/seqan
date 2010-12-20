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
..include:seqan/basic.h
*/

struct DotDrawing_;
typedef Tag<DotDrawing_> const DotDrawing;


/**
.Tag.HammingDistance
..summary:Switch to trigger Hamming distance, which is a measure of character substitutions.
..include:seqan/basic.h
*/

/**
.Tag.LevenshteinDistance
..summary:Switch to trigger Levenshtein distance, which is a measure of edit operations (character substitutions, deletions or insertions).
..remarks:$EditDistance$ is a synonym for $LevenshteinDistance$.
..see:Spec.EditDistance
..include:seqan/basic.h
*/

struct _HammingDistance;
struct _LevenshteinDistance;

typedef Tag<_HammingDistance>		HammingDistance;
typedef Tag<_LevenshteinDistance>	LevenshteinDistance;
typedef Tag<_LevenshteinDistance>	EditDistance; 


//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Alignment: Tags
//////////////////////////////////////////////////////////////////////////////
//Sollte eigentlich nach align/, aber da jetzt ja so viele
//alignment algorithmen in graph/ gelandet sind...

/**
.Tag.Global Alignment Algorithms:
..summary:Global alignment algorithm used by globalAlignment.
..see:Function.globalAlignment
..see:Tag.Local Alignment Algorithms
..include:seqan/basic.h
*/

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.NeedlemanWunsch:
	Dynamic programming algorithm for alignments by Needleman and Wunsch.
..include:seqan/basic.h
*/

struct _NeedlemanWunsch;
typedef Tag<_NeedlemanWunsch> NeedlemanWunsch;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.BandedNeedlemanWunsch:
	The Needleman-Wunsch alignment algorithm in a banded version.
..include:seqan/basic.h
*/
struct _BandedNeedlemanWunsch;
typedef Tag<_BandedNeedlemanWunsch> BandedNeedlemanWunsch;


//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.Gotoh:
	Gotoh's affine gap cost alignment algorithm.
..include:seqan/basic.h
*/
struct _Gotoh;
typedef Tag<_Gotoh> Gotoh;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.BandedGotoh:
	Gotoh's affine gap cost alignment algorithm in a banded version.
..include:seqan/basic.h
*/
struct _BandedGotoh;
typedef Tag<_BandedGotoh> BandedGotoh;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.MyersBitVector:
	Myers' bit vector alignment algorithm for edit distance.
	Note that this algorithm does not returns the alignment itself, but only computes the score.
..include:seqan/basic.h
*/
struct MyersBitVector_;
typedef Tag<MyersBitVector_> const MyersBitVector;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.MyersHirschberg:
	Myers' bit vector algorithm for edit distance combined with Hirschberg's linear space alignment algorithm.
..include:seqan/basic.h
*/
struct MyersHirschberg_;
typedef Tag<MyersHirschberg_> const MyersHirschberg;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.Hirschberg:
	Hirschberg's linear space global alignment algorithm.
..include:seqan/basic.h
*/
struct Hirschberg_;
typedef Tag<Hirschberg_> const Hirschberg;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Global Alignment Algorithms.value.Lcs:
	Longest common subsequence algorithm.
..include:seqan/basic.h
*/
struct Lcs_;
typedef Tag<Lcs_> const Lcs;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms:
..summary:Local alignment algorithm used by localAlignment.
..see:Function.localAlignment
..include:seqan/basic.h
*/

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.SmithWaterman:
	Triggers a Smith Waterman local alignment algorithm.
..include:seqan/basic.h
*/
struct _SmithWaterman;
typedef Tag<_SmithWaterman> const SmithWaterman;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.BandedSmithWaterman:
	Triggers a banded version of the Smith Waterman local alignment algorithm.
..include:seqan/basic.h
*/
struct _BandedSmithWaterman;
typedef Tag<_BandedSmithWaterman> const BandedSmithWaterman;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.WatermanEggert:
	Local alignment algorithm by Waterman and Eggert with "declumping" (i.e. only non-overlapping local alignments are computed).
.Tag.Local Alignment Algorithms.value.SmithWatermanClump:
	Same as $WatermanEggert$.
..include:seqan/basic.h
*/
struct _SmithWatermanClump;
typedef Tag<_SmithWatermanClump> const SmithWatermanClump;
typedef Tag<_SmithWatermanClump> const WatermanEggert;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Local Alignment Algorithms.value.BandedWatermanEggert:
	Triggers a banded version of the local alignment algorithm by Waterman and Eggert with "declumping".
.Tag.Local Alignment Algorithms.value.BandedSmithWatermanClump:
	Same as $BandedWatermanEggert$.
..include:seqan/basic.h
*/
struct _BandedWatermanEggert;
typedef Tag<_BandedWatermanEggert> const BandedSmithWatermanClump;
typedef Tag<_BandedWatermanEggert> const BandedWatermanEggert;

//////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Tag.RNA Folding Algorithms.value.Nussinov:
	Nussinov style RNA folding algorithm
..include:seqan/basic.h
*/
struct Nussinov_;
typedef Tag<Nussinov_> const Nussinov;

//////////////////////////////////////////////////////////////////////////////

struct TagBlat_;
typedef Tag<TagBlat_> const Blat;


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
