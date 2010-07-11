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
  Specialization "Chained" for class Seed.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_SEGMENT_H_
#define SEQAN_SEEDS_SEEDS_SEED_SEGMENT_H_

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// TODO(holtgrew): Maybe rename to SeedDiagonal?
// TODO(holtgrew): Dddoc this!

/*
Store the information about a seed segment.
 */
template <typename TPosition, typename TSize>
class SeedSegment
{
public:
    TPosition leftDim0;
    TPosition leftDim1;
    TSize length;

    SeedSegment()
            : leftDim0(0), leftDim1(0), length(0)
    {}

    SeedSegment(TPosition _leftDim0, TPosition _leftDim1, TSize _length)
            : leftDim0(_leftDim0), leftDim1(_leftDim1), length(_length)
    {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TPosition, typename TSize>
inline TSize
length(SeedSegment<TPosition, TSize> const & seedSegment)
{
    SEQAN_CHECKPOINT;
    return seedSegment.length;
}

template <typename TPosition, typename TSize>
inline void
setLength(SeedSegment<TPosition, TSize> & seedSegment, TSize newLength)
{
    SEQAN_CHECKPOINT;
    seedSegment.length = newLength;
}

template <typename TPosition, typename TSize>
inline TSize
leftDim0(SeedSegment<TPosition, TSize> const & seedSegment)
{
    SEQAN_CHECKPOINT;
    return seedSegment.leftDim0;
}

template <typename TPosition, typename TSize>
inline void
setLeftDim0(SeedSegment<TPosition, TSize> & seedSegment, TPosition newDim)
{
    SEQAN_CHECKPOINT;
    seedSegment.leftDim0 = newDim;
}

template <typename TPosition, typename TSize>
inline TSize
leftDim1(SeedSegment<TPosition, TSize> const & seedSegment)
{
    SEQAN_CHECKPOINT;
    return seedSegment.leftDim1;
}

template <typename TPosition, typename TSize>
inline void
setLeftDim1(SeedSegment<TPosition, TSize> & seedSegment, TPosition newDim)
{
    SEQAN_CHECKPOINT;
    seedSegment.leftDim1 = newDim;
}

#endif  // SEQAN_SEEDS_SEEDS_SEED_SEGMENT_H_
