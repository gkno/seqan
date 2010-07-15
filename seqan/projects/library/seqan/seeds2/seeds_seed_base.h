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

namespace seqan {

// TODO(holtgrew): Add possibility to generate TikZ code.

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

// Mixin member for mixing in scores into seeds.
template <typename TScore>
struct _ScoreMixin
{
    TScore score;
};

// Default configuration for seeds without score.
struct DefaultSeedConfig
{
    typedef size_t TPosition;
    typedef size_t TSize;
    typedef _MakeSigned<size_t>::Type TDiagonal;
    typedef False THasScore;
    typedef Nothing TScoreValue;
    typedef Nothing TScoreMixin;
};

// Default configuration for seeds with score.
struct DefaultSeedConfigScore
{
    typedef size_t TPosition;
    typedef size_t TSize;
    typedef _MakeSigned<size_t>::Type TDiagonal;
    typedef True THasScore;
    typedef int TScoreValue;
    typedef _ScoreMixin<int> TScoreMixin;
};

/**
.Class.Seed:
..summary:Describe a seed.
..cat:Seed Handling
..signature:Seed<TSpec, TConfig>
..param.TSpec:The seed specialization type.
..param.TConfig:The configuration object to use for this seed.
..include:seqan/seeds.h
 */
template <typename TSpec, typename TConfig = DefaultSeedConfig>
class Seed;

// ===========================================================================
// Metafunctions
// ===========================================================================

/**
.Metafunction.Position.param.T:type:Class.Seed
 */
template <typename TSpec, typename TConfig>
struct Position<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TPosition Type;
};

template <typename TSpec, typename TConfig>
struct Position<Seed<TSpec, TConfig> const>
        : Position<Seed<TSpec, TConfig> > {};

/**
.Metafunction.Size.param.T:type:Class.Seed
 */
template <typename TSpec, typename TConfig>
struct Size<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TSize Type;
};

template <typename TSpec, typename TConfig>
struct Size<Seed<TSpec, TConfig> const>
        : Size<Seed<TSpec, TConfig> > {};

/**
.Metafunction.Diagonal:
..cat:Seed Handling
..summary:Returns type of the value for the diagonal of a seed.
..signature:Diagonal<T>::Type
..param.T:Type of the seed to retrieve the diagonal for.
...type:Class.Seed
 */
template <typename T>
struct Diagonal;

template <typename TSpec, typename TConfig>
struct Diagonal<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TDiagonal Type;
};

template <typename TSpec, typename TConfig>
struct Diagonal<Seed<TSpec, TConfig> const>
        : Diagonal<Seed<TSpec, TConfig> > {};

/**
.Metafunction.HasScore:
..cat:Seed Handling
..summary:Returns True if the seed stores a score, False otherwise.
..signature:HasScore<T>::Type
..param.T:Type of the seed to retrieve whether it has a score for.
...type:Class.Seed
 */
template <typename T>
struct HasScore;

template <typename TSpec, typename TConfig>
struct HasScore<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::THasScore Type;
};

template <typename TSpec, typename TConfig>
struct HasScore<Seed<TSpec, TConfig> const>
        : HasScore<Seed<TSpec, TConfig> > {};

/**
.Metafunction.SeedScore:
..cat:Seed Handling
..summary:Returns type of the value for the score of a seed.
..signature:SeedScore<T>::Type
..param.T:Type of the seed to retrieve the score for.
...type:Class.Seed
 */
template <typename T>
struct SeedScore;

template <typename TSpec, typename TConfig>
struct SeedScore<Seed<TSpec, TConfig> >
{
    typedef typename TConfig::TScoreValue Type;
};

template <typename TSpec, typename TConfig>
struct SeedScore<Seed<TSpec, TConfig> const>
        : SeedScore<Seed<TSpec, TConfig> > {};

// ===========================================================================
// Functions
// ===========================================================================

/**
.Function.getLeftDim0:
..summary: Returns the first position of the seed in the query.
..cat:Seed Handling
..signature:leftDim0(seed)
..param.seed:The seed whose query position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
..include:seqan/seeds.h
*/

/**
.Function.getRightDim0:
..summary: Returns the last position of the seed in the query.
..cat:Seed Handling
..signature:rightDim0(seed)
..param.seed:The seed whose last in the query position should be returned.
...type:Class.Seed
..returns: End of the seed.
..include:seqan/seeds.h
*/

/**
.Function.getLeftDim1:
..summary: Returns the first position of the seed in the database.
..cat:Seed Handling
..signature:leftDim1(seed)
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns: Begin of the seed.
..include:seqan/seeds.h
*/

/**
.Function.getRightDim1:
..summary: Returns the last position of the seed in the database.
..cat:Seed Handling
..signature:rightDim1(seed)
..param.seed:The seed whose last in the database position should be returned.
...type:Class.Seed
..returns: End of the seed.
..include:seqan/seeds.h
*/

/**
.Function.getLowerDiagonal:
..summary: Returns the most left diagonal of the seed (maximum diagonal value).
..cat:Seed Handling
..signature:lowerDiagonal(seed)
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns:The most left diagonal.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
getLowerDiagonal(Seed<TSpec, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._lowerDiagonal;
}

/**
.Function.setLowerDiagonal:
..summary: Sets a new value for the most left diagonal.
..cat:Seed Handling
..signature:setLeftDiagonal(seed, diag)
..param.seed:The seed whose left diagonal value should be updated.
...type:Class.Seed
..param.diag:The new value for the most left diagonal.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig, typename TPosition>
inline void 
setLowerDiagonal(Seed<TSpec, TConfig> & seed, 
				 TPosition newDiag)
{
	SEQAN_CHECKPOINT;
	seed._lowerDiagonal = newDiag;
}

/**
.Function.getUpperDiagonal:
..summary: Returns the most right diagonal of the seed (minimum diagonal value).
..cat:Seed Handling
..signature:upperDiagonal(seed)
..param.seed:The seed whose database position should be returned.
...type:Class.Seed
..returns:The most right diagonal.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
getUpperDiagonal(Seed<TSpec, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
    return seed._upperDiagonal;
}

/**
.Function.setUpperDiagonal:
..summary: Sets a new value for the most right diagonal.
..cat:Seed Handling
..signature:setRightDiagonal(seed, diag)
..param.seed:The seed whose right diagonal value should be updated.
...type:Class.Seed
..param.diag:The new value for the most right diagonal.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig, typename TPosition>
inline void 
setUpperDiagonal(Seed<TSpec, TConfig> & seed, 
				 TPosition newDiag)
{
	SEQAN_CHECKPOINT;
	seed._upperDiagonal = newDiag;
}

// Computed values, based on properties returned by getters.

/**
.Function.getStartDiagonal:
..summary: Returns the diagonal of the start point.
..cat:Seed Handling
..signature:startDiagonal(seed)
..param.seed:The seed whose start diagonal should be returned.
...type:Class.Seed
..returns:The diagonal of the start point.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
getStartDiagonal(Seed<TSpec, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
    return getLeftDim1(seed) - getLeftDim0(seed);
}

/**
.Function.getEndDiagonal:
..summary: Returns the diagonal of the end point.
..cat:Seed Handling
..signature:endDiagonal(seed)
..param.seed:The seed whose end diagonal should be returned.
...type:Class.Seed
..returns:The diagonal of the end point.
..include:seqan/seeds.h
*/
template <typename TSpec, typename TConfig>
inline typename Diagonal<Seed<TSpec, TConfig> >::Type
getEndDiagonal(Seed<TSpec, TConfig> const & seed)
{
	SEQAN_CHECKPOINT;
    return getRightDim1(seed) - getRightDim0(seed);
}

// TODO(holtgrew): COULD introduce {get,set}{Left,Right}(dim, value)

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_BASE_H_
