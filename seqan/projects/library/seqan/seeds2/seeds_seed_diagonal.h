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
  SeedDiagonal objects store information about seed parts/segments of a
  ChainedSeed.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_
#define SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/*
.Class.SeedDiagonal:
..summary:Store the information about a seed segment.
..signature:SeedDiagonal<TPosition, TSize>
..cat:Seed Handling
..param.TPosition:The type to use for positions.
..param.TSize:The type to use for the seed length.
*/
template <typename TPosition, typename TSize>
class SeedDiagonal
{
public:
/**
.Memvar.SeedDiagonal#leftDim0:The position in the query (dimension 0).
*/
    TPosition leftDim0;
/**
 .Memvar.SeedDiagonal#leftDim0:The position in the database sequence (dimension 1).
*/
    TPosition leftDim1;
/**
.Memvar.SeedDiagonal#length:The length of the diagonal.
*/
    TSize length;
    
    SeedDiagonal()
            : leftDim0(0), leftDim1(0), length(0)
    { SEQAN_CHECKPOINT; }

    SeedDiagonal(TPosition _leftDim0, TPosition _leftDim1, TSize _length)
            : leftDim0(_leftDim0), leftDim1(_leftDim1), length(_length)
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================
  
/**
.Metafunction.Position.param.T:SeedDiagonal
*/
template <typename TPosition, typename TSize>
struct Position<SeedDiagonal<TPosition, TSize> >
{
    typedef TPosition Type;
};

template <typename TPosition, typename TSize>
struct Position<SeedDiagonal<TPosition, TSize> const>
  : Position<SeedDiagonal<TPosition, TSize> > {};

/**
.Metafunction.Size.param.T:SeedDiagonal
*/
template <typename TPosition, typename TSize>
struct Size<SeedDiagonal<TPosition, TSize> >
{
  typedef TSize Type;
};

template <typename TPosition, typename TSize>
struct Size<SeedDiagonal<TPosition, TSize> const>
  : Size<SeedDiagonal<TPosition, TSize> > {};
  
// ===========================================================================
// Functions
// ===========================================================================

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_DIAGONAL_H_
