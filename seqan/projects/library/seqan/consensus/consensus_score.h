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
==========================================================================*/

#ifndef SEQAN_HEADER_SEQAN_CONSENSUS_SCORE_H
#define SEQAN_HEADER_SEQAN_CONSENSUS_SCORE_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Consensus score tags
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

struct ConsensusScore_;
typedef Tag<ConsensusScore_> const ConsensusScore;

//////////////////////////////////////////////////////////////////////////////

struct FractionalScore_;
typedef Tag<FractionalScore_> const FractionalScore;

//////////////////////////////////////////////////////////////////////////////

template<typename TScore1, typename TScore2>
struct WeightedConsensusScore;




//////////////////////////////////////////////////////////////////////////////
// Scoring classes
//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
// ConsensusScore
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, ConsensusScore>
{
public:
	TValue data_penalty; // Not a consensus element
	TValue data_match;	 // Is a consensus element

	String<bool> consensus_set;		// Is the alphabet character part of the consensus set in the given column
	int column;


public:
	Score(): data_penalty(-1), data_match(0), column(-1) {}

};



template <typename TValue, typename TPos1, typename TSeq1>
inline void
_update(Score<TValue, ConsensusScore>& me,
		TPos1 const pos1,
		TSeq1 const& seq1)
{
	typedef typename Value<TSeq1>::Type TProfileType;
	typedef typename Size<TProfileType>::Type TSize;
	TProfileType const& myC = value(seq1, pos1);
	TSize maxCount = 0;
	TSize alphSize = ValueSize<TProfileType>::VALUE;
	resize(me.consensus_set, alphSize);
	for(TSize i = 0; i<alphSize; ++i) {
		if (myC.count[i] > maxCount) maxCount = myC.count[i];
	}
	for(TSize i = 0; i<alphSize; ++i) {
		if (myC.count[i] == maxCount) value(me.consensus_set, i) = true;
		else value(me.consensus_set, i) = false;
	}
	me.column = pos1;
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, ConsensusScore> const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &seq1,
	TSeq2 const &)
{
	if ((int) pos2 < (int) 0) return me.data_penalty;
	if ((int) pos1 != me.column) _update(const_cast<Score<TValue, ConsensusScore>&>(me), pos1, seq1);
	return (getValue(me.consensus_set, ValueSize<typename Value<TSeq1>::Type>::VALUE - 1)) ? me.data_match : me.data_penalty;
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, ConsensusScore> const & me,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return me.data_penalty;
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, ConsensusScore> const & me,
	  TPos1 pos1,
	  TPos2 pos2,
	  TSeq1 const &seq1,
	  TSeq2 const &seq2)
{
	if ((int) pos1 != me.column) _update(const_cast<Score<TValue, ConsensusScore>&>(me), pos1, seq1);
	return (getValue(me.consensus_set, seq2[pos2].count[0])) ? me.data_match : me.data_penalty;
}









//////////////////////////////////////////////////////////////////////////////
// FractionalScore
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, FractionalScore>
{
public:
	TValue data_penalty; 
	
	int sum;		// Total number of profile characters in the given column
	int column;


public:
	Score(): data_penalty(-1), column(-1) {}

};


template <typename TValue, typename TPos1, typename TSeq1>
inline void
_update(Score<TValue, FractionalScore>& me,
		TPos1 const pos1,
		TSeq1 const& seq1)
{
	typedef typename Size<TSeq1>::Type TSize;
	me.sum = 0;
	for(TSize i = 0; i < (TSize) ValueSize<typename Value<TSeq1>::Type>::VALUE; ++i) me.sum += (value(seq1, pos1)).count[i];
	me.column = pos1;
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, FractionalScore> const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &seq1,
	TSeq2 const &)
{
	if ((int) pos2 < 0) return me.data_penalty;
	if ((int) pos1 != me.column) _update(const_cast<Score<TValue, FractionalScore>&>(me), pos1, seq1);
	return me.data_penalty * ((me.sum - (TValue) (value(seq1, pos1)).count[ValueSize<typename Value<TSeq1>::Type>::VALUE - 1]) / me.sum);
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, FractionalScore> const & me,
	TPos1,
	TPos2,
	TSeq1 const &,
	TSeq2 const &)
{
	return me.data_penalty;
}


template <typename TValue, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, FractionalScore> const & me,
	  TPos1 pos1,
	  TPos2 pos2,
	  TSeq1 const &seq1,
	  TSeq2 const &seq2)
{
	if ((int) pos1 != me.column) _update(const_cast<Score<TValue, FractionalScore>&>(me), pos1, seq1);
	return me.data_penalty * ((me.sum - (TValue) seq1[pos1].count[seq2[pos2].count[0]]) / me.sum);
}








//////////////////////////////////////////////////////////////////////////////
// WeightedConsensusScore
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TScore1, typename TScore2>
class Score<TValue, WeightedConsensusScore<TScore1, TScore2> >
{
public:
	TScore1 sc1;
	TScore2 sc2;

public:
	Score() {}

};


template <typename TValue, typename TScore1, typename TScore2, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
	Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const &seq1,
	TSeq2 const &seq2)
{
	return (0.5 * (TValue) scoreGapExtendHorizontal(me.sc1, pos1, pos2, seq1, seq2) + 0.5 * (TValue) scoreGapExtendHorizontal(me.sc2, pos1, pos2, seq1, seq2));
}


template <typename TValue, typename TScore1, typename TScore2, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
	Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
	TPos1 pos1,
	TPos2 pos2,
	TSeq1 const & seq1,
	TSeq2 const & seq2)
{
	return (0.5 * (TValue) scoreGapExtendVertical(me.sc1, pos1, pos2, seq1, seq2) + 0.5 * (TValue) scoreGapExtendVertical(me.sc2, pos1, pos2, seq1, seq2));
}


template <typename TValue, typename TScore1, typename TScore2, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, WeightedConsensusScore<TScore1, TScore2> > const & me,
	  TPos1 pos1,
	  TPos2 pos2,
	  TSeq1 const &seq1,
	  TSeq2 const &seq2)
{
	return (0.5 * (TValue) score(me.sc1, pos1, pos2, seq1, seq2) + 0.5 * (TValue) score(me.sc2, pos1, pos2, seq1, seq2));
}

}

#endif

