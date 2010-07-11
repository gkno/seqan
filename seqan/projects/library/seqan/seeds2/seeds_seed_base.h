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
  The class Seed.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_BASE_H_
#define SEQAN_SEEDS_SEEDS_SEED_BASE_H_

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Tag.Seed Specs
..summary:Specialization tags for @Class.Seed@.
..cat:Seed Handling
..tag.Simple:Simple seed that only stores coordinates of the begin end end point and the diagonals.
..tag.Chained:Seed that stores the dot-plot diagonals (i.e. matches) that the seed consists of.
..see:Spec.SimpleSeed
..see:Spec.ChainedSeed
..include:seqan/seeds.h
*/

/**
.Class.Seed:
..summary:Describe a seed.
..cat:Seed Handling
..signature:Seed<TPosition, TSpec>
..param.TPosition:The type for storing position.
..param.TSpec:The seed specialization type.
..include:seqan/seeds.h
 */
// TODO(holtgrew): Maybe add TSize?
template <typename TPosition, typename TSpec>
class Seed;

// ===========================================================================
// Metafunctions
// ===========================================================================

/**
.Metafunction.Position.param.T:type:Class.Seed
 */
template <typename TPosition, typename TSpec>
struct Position<Seed<TPosition, TSpec> >
{
    typedef TPosition Type;
};

template <typename TPosition, typename TSpec>
struct Position<Seed<TPosition, TSpec> const>
        : Position<Seed<TPosition, TSpec> > {};

// ===========================================================================
// Functions
// ===========================================================================

// TODO(holtgrew): Probably, all the getter should begin with "get"

/**
.Function.startDiagonal:
..summary: Returns the diagonal of the start point.
..cat:Seed Handling
..signature:startDiagonal(seed)
..param.seed:The seed whose start diagonal should be returned.
...type:Class.Seed
..returns:The diagonal of the start point.
..include:seqan/seeds.h
*/

/**
.Function.endDiagonal:
..summary: Returns the diagonal of the end point.
..cat:Seed Handling
..signature:endDiagonal(seed)
..param.seed:The seed whose end diagonal should be returned.
...type:Class.Seed
..returns:The diagonal of the end point.
..include:seqan/seeds.h
*/

/**
.Function.leftPosition:
..summary:The begin position of segment in a seed.
..cat:Seed Handling
..signature:leftPosition(seed, dim)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..returns:Begin position of the $dim$-th segment in $seed$.
..include:seqan/seeds.h
*/

/**
.Function.setLeftPosition:
..summary:Sets begin position of segment in a seed.
..cat:Seed Handling
..signature:setLeftPosition(seed, dim, new_pos)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..param.new_pos:The new begin position of the $dim$-th segment in $seed$.
..see:Function.leftPosition
..include:seqan/seeds.h
*/

/**
.Function.rightPosition:
..summary:The end position of segment in a seed.
..cat:Seed Handling
..signature:rightPosition(seed, dim)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..returns:End position of the $dim$-th segment in $seed$.
..see:Function.leftPosition
..include:seqan/seeds.h
*/

/**
.Function.setRightPosition:
..summary:Sets end position of segment in a seed.
..cat:Seed Handling
..signature:setRightPosition(seed, dim, new_pos)
..param.seed:A seed.
...type:Class.Seed
..param.dim:Number of segment.
...remarks:$dim <$ @Function.dimension.dimension(seed)@
..param.new_pos:The new end position of the $dim$-th segment in $seed$.
..see:Function.rightPosition
..see:Function.setLeftPosition
..include:seqan/seeds.h
*/

/**
.Function.dimension:
..summary:Dimension of a seed.
..cat:Seed Handling
..signature:dimension(seed)
..param.seed:A seed.
...type:Class.Seed
..returns:The number of reference sequences for $seed$.
..include:seqan/seeds.h
*/

/**
.Function.leftDim0:
..summary: Returns the first position of the seed in the query.
..cat:Seed Handling
..signature:leftDim0(seed)
..param.seed:The seed whose query position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
..include:seqan/seeds.h
*/

/**
.Function.rightDim0:
..summary: Returns the last position of the seed in the query.
..cat:Seed Handling
..signature:rightDim0(seed)
..param.seed:The seed whose last in the query position should be returned.
...type:Class.Seed
..returns: End of the seed.
..include:seqan/seeds.h
*/

/**
.Function.leftDim1:
..summary: Returns the first position of the seed in the database.
..cat:Seed Handling
..signature:leftDim1(seed)
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
..include:seqan/seeds.h
*/

/**
.Function.rightDim1:
..summary: Returns the last position of the seed in the database.
..cat:Seed Handling
..signature:rightDim1(seed)
..param.seed:The seed whose last in the database position should be returned.
...type:Class.Seed
..returns: End of the seed.
..include:seqan/seeds.h
*/

/**
.Function.leftDiagonal:
..summary: Returns the most left diagonal of the seed (maximum diagonal value).
..cat:Seed Handling
..signature:leftDiagonal(seed)
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns:The most left diagonal.
..include:seqan/seeds.h
*/

/**
.Function.rightDiagonal:
..summary: Returns the most right diagonal of the seed (minimum diagonal value).
..cat:Seed Handling
..signature:rightDiagonal(seed)
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns:The most right diagonal.
..include:seqan/seeds.h
*/

/**
.Function.length.param.object.type:Class.Seed
*/

/**
.Function.setLeftDim0:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setLeftDim0(seed, start)
..param.seed:The seed whose start position should be updated.
...type:Class.Seed
..param.start:The query position where the seed should start.
..include:seqan/seeds.h
*/

/**
.Function.setRightDim0:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setRightDim0(seed, end)
..param.seed:The seed whose end position should be updated.
...type:Class.Seed
..param.end:The query position where the seed should end.
..include:seqan/seeds.h
*/

/**
.Function.setLeftDim1:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setLeftDim1(seed, start)
..param.seed:The seed whose start position should be updated.
...type:Class.Seed
..param.start:The database position where the seed should start.
..include:seqan/seeds.h
*/

/**
.Function.setRightDim1:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setRightDim1(seed, end)
..param.seed:The seed whose end position should be updated.
...type:Class.Seed
..param.end:The database position where the seed should end.
..include:seqan/seeds.h
*/

/**
.Function.setLeftDiagonal:
..summary: Sets a new value for the most left diagonal.
..cat:Seed Handling
..signature:setLeftDiagonal(seed, diag)
..param.seed:The seed whose left diagonal value should be updated.
...type:Class.Seed
..param.diag:The new value for the most left diagonal.
..include:seqan/seeds.h
*/

/**
.Function.setRightDiagonal:
..summary: Sets a new value for the most right diagonal.
..cat:Seed Handling
..signature:setRightDiagonal(seed, diag)
..param.seed:The seed whose right diagonal value should be updated.
...type:Class.Seed
..param.diag:The new value for the most right diagonal.
..include:seqan/seeds.h
*/

/**
.Function.extendSeed
..summary:Extends a seed.
..cat:Seed Handling
..signature:extendSeed(seed, query, database, direction, tag)
..signature:extendSeed(seed, scoreDropOff, scoreMatrix, query, database, direction, tag)
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

#endif  // SEQAN_SEEDS_SEEDS_SEED_BASE_H_
