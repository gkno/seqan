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

  Based on the code by Carsten Kemena <carsten.kemena@crg.es>, debugged
  by Birte Kehr <birte.kehr@fu-berlin.de>.
  ============================================================================
  Generic parts of the seed extension algorithms.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_EXTENSION_H_
#define SEQAN_SEEDS_SEEDS_EXTENSION_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Tag.Seed Extension
..cat:Seed Handling
..summary:The algorithms used to extend a seed.
..see:Function.extendSeed
..see:Function.extendSeeds
..see:Function.extendSeedScore
..see:Function.extendSeedsScore
..tag.MatchExtend:Extends a seed until a mismatch occurs.
..tag.UngappedXDrop:Ungapped extension of a seed until score drops below a Value.
..tag.GappedXDrop:Gapped extension of a seed until score drops below a Value. Only @Spec.SimpleSeed@s.
..include:seqan/seeds.h
*/
struct _MatchExtend;
typedef Tag<_MatchExtend> const MatchExtend;

struct _UnGappedXDrop;
typedef Tag<_UnGappedXDrop> const UnGappedXDrop;

struct _GappedXDrop;
typedef Tag<_GappedXDrop> const GappedXDrop;

enum ExtensionDirection
{
    EXTEND_LEFT,
    EXTEND_RIGHT,
    EXTEND_BOTH
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

/**
.Function.extendSeed
..summary:Extends a seed.
..cat:Seed Handling
..signature:extendSeed(seed, query, database, direction, MatchExtend)
..signature:extendSeed(seed, query, database, direction, scoreDropOff, scoreMatrix, {UngappedXDrop, GappedXDrop})
..param.seed: The seed to extend.
...type:Class.Seed
..param.query: The query sequence.
...type:Class.String
..param.query: The database sequence.
...type:Class.String
..param.direction: Defines the direction in which the seed should be extended. 0 = left, 1 = right, 2 = both
..param.scoreDropOff: The score drop after which the extension should stop. The extension stops if this value is exceeded.
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.scoreMatrix: The scoring scheme.
...type:Spec.Simple Score
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.tag: The algorithm to use.
...type:Tag.Seed Extension.MatchExtend
...type:Tag.Seed Extension.UngappedXDrop
...type:Tag.Seed Extension.GappedXDrop
*/
template <typename TConfig, typename TQuery, typename TDatabase>
inline void 
extendSeed(Seed<Simple, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
		   MatchExtend const &)
{
    // For Simple Seeds, we can simply update the begin and end values
    // in each dimension.
	SEQAN_CHECKPOINT;

    typedef Seed<Simple, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
		TPosition posDim0 = getBeginDim0(seed) ;
		TPosition posDim1 = getBeginDim1(seed);
		while (posDim0 >= 1 && posDim1 >= 1 && query[posDim0 - 1] == database[posDim1 - 1]) {
			--posDim0;
			--posDim1;
		}
		setBeginDim0(seed, posDim0);
		setBeginDim1(seed, posDim1);
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
		TPosition posDim0 = getEndDim0(seed) ;
		TPosition posDim1 = getEndDim1(seed);
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && query[posDim0] == database[posDim1]) {
			++posDim0;
			++posDim1;
		}
		setEndDim0(seed, posDim0);
		setEndDim1(seed, posDim1);
	}
}


template <typename TConfig, typename TQuery, typename TDatabase>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
		   MatchExtend const &)
{
    // For Chained Seeds, we extend the first and the last Seed
    // Diagonal.
	SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(seed), 0u);

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Value<TSeed>::Type TSeedDiagonal;
    typedef typename Position<TSeedDiagonal>::Type TPosition;
    typedef typename Size<TSeedDiagonal>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        TSeedDiagonal & diag = front(seed);
		TPosition posDim0 = diag.beginDim0;
		TPosition posDim1 = diag.beginDim1;
        TSize diagonalLength = diag.length;
		while (posDim0 >= 1 && posDim1 >= 1 && query[posDim0 - 1] == database[posDim1 - 1]) {
			--posDim0;
			--posDim1;
            ++diagonalLength;
		}
        diag.beginDim0 = posDim0;
        diag.beginDim1 = posDim1;
        diag.length = diagonalLength;
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
        TSeedDiagonal & diag = back(seed);
		TPosition posDim0 = diag.beginDim0 + diag.length;
		TPosition posDim1 = diag.beginDim1 + diag.length;
        TSize diagonalLength = diag.length;
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && query[posDim0] == database[posDim1]) {
			++posDim0;
			++posDim1;
            ++diagonalLength;
		}
        diag.length = diagonalLength;
	}
}

/**
.Function.extendSeeds
..summary: Extension of seeds.
..cat:Seed Handling
..signature:extendSeeds(container, query, database, direction, MatchExtend)
..signature:extendSeeds(begin, end, query, database, direction, MatchExtend)
..signature:extendSeeds(container, query, database, direction, scoreDropOff, scoreMatrix, {UngappedXDrop, GappedXDrop})
..signature:extendSeeds(begin, end, query, database, direction, scoreDropOff, scoreMatrix, {UngappedXDrop, GappedXDrop})
..param.container: The container with the @Class.Seed@ objects to extend.
...type:Concept.Container
..param.begin: Iterator pointing to the first value to add.
..param.end: Iterator pointing just behind the last value to add.
..param.query: The query sequence.
...type:Class.String
..param.query: The database sequence.
...type:Class.String
..param.direction: Defines the direction in which the seed should be extended. 0 = left, 1 = right, 2 = both
..param.scoreDropOff: The score drop after which the extension should stop. The extension stops if this value is exceeded.
..param.scoreMatrix: The scoring sheme.
...type:Spec.Simple Score
..param.tag: The algorithm to use.
...type:Tag.Seed Extension.MatchExtend
...type:Tag.Seed Extension.UngappedXDrop
...type:Tag.Seed Extension.GappedXDrop
*/

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_EXTENSION_H_
