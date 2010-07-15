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
  Specialization "Simple" for class Seed.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
#define SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// Already defined in module score.
// struct _Simple {};
// typedef Tag<_Simple> Simple;  // Already defined in module base.
struct Simple {}; // .. but type is not complete yet

/**
.Spec.SimpleSeed:
..summary:Describes a seed with start and end position and diagonal upper and lower bounds.
..cat:Seed Handling
..general:Class.Seed
..signature:Seed<TPosition, SimpleSeed>
..param.TPosition:The type number that schuld be used. Must have negative numbers (e.g. int/long).
.Memfunc.SimpleSeed#Seed:
..class:Spec.SimpleSeed
..summary:Constructor
..signature: Seed<TPosition, SimpleSeed> ()
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, length)
..signature: Seed<TPosition, SimpleSeed> (qStartPos, dStartPos, qEndPos, dEndPos)
..param.qStartPos: Start in query sequence.
..param.dStartPos: Start in database sequence.
..param.qEndPos: End in query sequence.
..param.dEndPos: End in database sequence.
..param.length: Length of the seed.
..include:seqan/seeds.h
*/
template <typename TConfiguration>
class Seed<Simple, TConfiguration>
        : TConfiguration::TScoreMixin
{
public:
    typedef typename TConfiguration::TPosition TPosition;
    typedef typename TConfiguration::TSize TSize;
    typedef typename TConfiguration::TDiagonal TDiagonal;

    TPosition _leftDim0;
    TPosition _leftDim1;
    TPosition _rightDim0;
    TPosition _rightDim1;
    TDiagonal _lowerDiagonal;
    TDiagonal _upperDiagonal;

    Seed()
            : _leftDim0(0), _leftDim1(0), _rightDim0(0), _rightDim1(0),
              _lowerDiagonal(0), _upperDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition leftDim0, TPosition leftDim1, TPosition seedLength)
            : _leftDim0(leftDim0),
              _leftDim1(leftDim1),
              _rightDim0(leftDim0 + seedLength - 1),
              _rightDim1(leftDim1 + seedLength - 1),
              _lowerDiagonal(leftDim1 - leftDim0),
              _upperDiagonal(leftDim1 - leftDim0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition leftDim0, TPosition leftDim1, TPosition rightDim0,
         TPosition rightDim1)
            : _leftDim0(leftDim0),
              _leftDim1(leftDim1),
              _rightDim0(rightDim0),
              _rightDim1(rightDim1),
              _lowerDiagonal(_min(leftDim1 - leftDim0, rightDim1 - rightDim0)),
              _upperDiagonal(_max(leftDim1 - leftDim0, rightDim1 - rightDim0))
    { SEQAN_CHECKPOINT; }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getLeftDim0(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._leftDim0;
}

/**
.Function.setLeftDim0:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setLeftDim0(seed, start)
..param.seed:The seed whose start position should be updated.
...type:Spec.SimpleSeed
..param.start:The query position where the seed should start.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setLeftDim0(Seed<Simple, TConfig> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
	seed._leftDim0 = newLeftPosition;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getRightDim0(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;;
	return seed._rightDim0;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getLeftDim1(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._leftDim1;
}

/**
.Function.setLeftDim1:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setLeftDim1(seed, start)
..param.seed:The seed whose start position should be updated.
...type:Spec.SimpleType
..param.start:The database position where the seed should start.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setLeftDim1(Seed<Simple, TConfig> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
	seed._leftDim1 = newLeftPosition;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getRightDim1(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._rightDim1;
}

/**
.Function.setRightDim0:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setRightDim0(seed, end)
..param.seed:The seed whose end position should be updated.
...type:Spec.SimpleSeed
..param.end:The query position where the seed should end.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setRightDim0(Seed<Simple, TConfig> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
	seed._rightDim0 = newRightPosition;
}

/**
.Function.setRightDim1:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setRightDim1(seed, end)
..param.seed:The seed whose end position should be updated.
...type:Spec.Simple Seed
..param.end:The database position where the seed should end.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setRightDim1(Seed<Simple, TConfig> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
	seed._rightDim1 = newRightPosition;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
