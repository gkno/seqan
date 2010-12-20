// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================


//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file


#ifndef SEQAN_SEEDHANDLINGTAGS_H
#define SEQAN_SEEDHANDLINGTAGS_H

namespace SEQAN_NAMESPACE_MAIN
{
/**
.Tag.Seed Adding
..summary:The algorithm used to add a seed to a SeedSet.
..see:Function.addSeed
..see:Function.addSeeds
..cat:Seed Handling
..tag.Single:
	Simple adding of a seed to the set. No chaining or merging.
..tag.Blat:
	Chaining of seeds. Gap is filled with smaller matching segments.
..tag.Chaos:
	Chaining of seeds. One gap is introduced at the best position.
..tag.SimpleChain:
	Chaining of seeds.
..include:seqan/seeds.h
*/
struct _addSeeding_Single;
typedef Tag<_addSeeding_Single> const Single;

//also defined in blast_base.h (103)
//moved to basic_tag
/*
struct _Chain_Blat;
typedef Tag<_Chain_Blat> const Blat;
*/

struct _Chain_Chaos;
typedef Tag<_Chain_Chaos> const Chaos;

struct _Chain_Simple;
typedef Tag<_Chain_Simple> const SimpleChain;


template <typename T>
struct BLOCK_SIZE
{
};

/**
.Tag.SeedSet
..cat:Seed Handling
..summary: Tags for the behaviour of a SeedSet
..tag.DefaultScore: 
	Enables scoring of seeds in a SeedSet.
..tag.DefaultNoScore:
	Disables scoring of seeds in a SeedSet.
..include:seqan/seeds.h
*/



struct _Gap_Cost_Manhatten;
typedef Tag<_Gap_Cost_Manhatten> const Manhattan;

struct _Gap_Cost_QueryDistance;
typedef Tag<_Gap_Cost_QueryDistance> const QueryDistance;

struct _Gap_Cost_DatabaseDistance;
typedef Tag<_Gap_Cost_DatabaseDistance> const DatabaseDistance;

struct _Gap_Cost_NoGapCost;
typedef Tag<_Gap_Cost_NoGapCost> const NoGapCost;


struct _Good_Seed;
typedef Tag<_Good_Seed> const SeedScore;

struct _Good_Seed2;
typedef Tag<_Good_Seed2> const SeedLength;


template<typename T1, typename T2, typename T3>
struct Scoring_Scheme;

typedef Tag<Scoring_Scheme<SeedScore, Manhattan, int> > const DefaultScore;


typedef Tag<Scoring_Scheme<SeedLength, Manhattan, void> > const DefaultNoScore;


template<typename T>
struct GapCosts
{
	typedef T Type;
};

template <typename TGapCosts, typename TQualityFactor, typename TScore>
struct GapCosts<const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, TScore> > >
{
	typedef TGapCosts Type;
};

template<typename T>
struct QualityFactor
{
	typedef T Type;
};

template <typename TGapCosts, typename TQualityFactor, typename TScore>
struct QualityFactor<const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, TScore> > >
{
	typedef  TQualityFactor Type;
};

template<typename T>
struct ScoreType
{
	typedef T Type;
};

template <typename TGapCosts, typename TQualityFactor, typename TScore>
struct ScoreType<const Tag<Scoring_Scheme<TQualityFactor, TGapCosts, TScore> > >
{
	typedef  TScore Type;
};

/* see find_pattern_base.h
template <typename T>
struct ScoringScheme
{
	typedef T Type;
};
*/

} //namespace Seqan

#endif //#ifndef SEQAN_HEADER_
