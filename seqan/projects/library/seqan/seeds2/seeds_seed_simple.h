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

    TPosition _beginDim0;
    TPosition _beginDim1;
    TPosition _endDim0;
    TPosition _endDim1;
    TDiagonal _lowerDiagonal;
    TDiagonal _upperDiagonal;

    Seed()
            : _beginDim0(0), _beginDim1(0), _endDim0(0), _endDim1(0),
              _lowerDiagonal(0), _upperDiagonal(0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition beginDim0, TPosition beginDim1, TPosition seedLength)
            : _beginDim0(beginDim0),
              _beginDim1(beginDim1),
              _endDim0(beginDim0 + seedLength),
              _endDim1(beginDim1 + seedLength),
              _lowerDiagonal(beginDim1 - beginDim0),
              _upperDiagonal(beginDim1 - beginDim0)
    { SEQAN_CHECKPOINT; }

    Seed(TPosition beginDim0, TPosition beginDim1, TPosition endDim0,
         TPosition endDim1)
            : _beginDim0(beginDim0),
              _beginDim1(beginDim1),
              _endDim0(endDim0),
              _endDim1(endDim1),
              _lowerDiagonal(_min(beginDim1 - beginDim0, endDim1 - endDim0)),
              _upperDiagonal(_max(beginDim1 - beginDim0, endDim1 - endDim0))
    { SEQAN_CHECKPOINT; }

    template <typename TSeed2>
    Seed(TSeed2 const & other)
            : _beginDim0(getBeginDim0(other)),
              _beginDim1(getBeginDim1(other)),
              _endDim0(getEndDim0(other)),
              _endDim1(getEndDim1(other)),
              _lowerDiagonal(getLowerDiagonal(other)),
              _upperDiagonal(getUpperDiagonal(other))
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
getBeginDim0(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._beginDim0;
}

/**
.Function.setBeginDim0:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setBeginDim0(seed, start)
..param.seed:The seed whose start position should be updated.
...type:Spec.SimpleSeed
..param.start:The query position where the seed should start.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setBeginDim0(Seed<Simple, TConfig> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
	seed._beginDim0 = newLeftPosition;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getEndDim0(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;;
	return seed._endDim0;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getBeginDim1(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._beginDim1;
}

/**
.Function.setBeginDim1:
..summary: Updates the start point of the seed.
..cat:Seed Handling
..signature:setBeginDim1(seed, start)
..param.seed:The seed whose start position should be updated.
...type:Spec.SimpleType
..param.start:The database position where the seed should start.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setBeginDim1(Seed<Simple, TConfig> & seed, 
            TPosition newLeftPosition)
{
	SEQAN_CHECKPOINT;
	seed._beginDim1 = newLeftPosition;
}

template <typename TConfig>
inline typename Position<Seed<Simple, TConfig> >::Type
getEndDim1(Seed<Simple, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
	return seed._endDim1;
}

/**
.Function.setEndDim0:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setEndDim0(seed, end)
..param.seed:The seed whose end position should be updated.
...type:Spec.SimpleSeed
..param.end:The query position where the seed should end.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setEndDim0(Seed<Simple, TConfig> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
	seed._endDim0 = newRightPosition;
}

/**
.Function.setEndDim1:
..summary: Updates the end point of the seed.
..cat:Seed Handling
..signature:setEndDim1(seed, end)
..param.seed:The seed whose end position should be updated.
...type:Spec.Simple Seed
..param.end:The database position where the seed should end.
..include:seqan/seeds.h
*/
template <typename TConfig, typename TPosition>
inline void 
setEndDim1(Seed<Simple, TConfig> & seed, 
             TPosition newRightPosition)
{
	SEQAN_CHECKPOINT;
	seed._endDim1 = newRightPosition;
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SIMPLE_H_
