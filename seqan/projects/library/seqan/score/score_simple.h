#ifndef SEQAN_HEADER_SCORE_SIMPLE_H
#define SEQAN_HEADER_SCORE_SIMPLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Simple
..cat:Scoring
..general:Class.Score
..summary:Simple scoring scheme that has scores for matches, mismatches, opening gaps and extending gaps.
*/
template <typename TValue>
class Score<TValue, Simple>
{
public:
	TValue data_match;
	TValue data_mismatch;
	TValue data_gap_extend;
	TValue data_gap_open;

public:
	Score(TValue _match = 0, TValue _mismatch = -1, TValue _gap_extend = -1, TValue _gap_open = 0):
		data_match(_match),
		data_mismatch(_mismatch),
		data_gap_extend(_gap_extend),
		data_gap_open(_gap_open)
	{
	}

/**.Memfunc.Score#Score:
..class:Class.Score
..summary:Constructor
..signature:Score<TValue, Simple> ()
..signature:Score<TValue, Simple> (score)
..signature:Score<TValue, Simple> (match, mismatch,gap_extend,gap_open)
..param.score:Other Score object. (copy constructor)
..param.match:TValue object.
..param.mismatch:TValue object.
..param.gap_extend:TValue object.
..param.gap_open:TValue object.
...remarks:gap_open is the extra score for opening gaps, so whenever a new gap is openend, the score for the first gap will be gap_open + gap_extend.
*/	Score(Score const & other):
		data_match(other.data_match),
		data_mismatch(other.data_mismatch),
		data_gap_extend(other.data_gap_extend),
		data_gap_open(other.data_gap_open)
	{
	}
	~Score()
	{
	}

	Score & operator = (Score const & other)
	{
		data_match = other.data_match;
		data_mismatch = other.data_mismatch;
		data_gap_extend = other.data_gap_extend;
		data_gap_open = other.data_gap_open;
		return *this;
	}

//____________________________________________________________________________
};
//////////////////////////////////////////////////////////////////////////////

//Shortcut:

typedef Score< > SimpleScore;

//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
/**.Function.scoreMatch:
..class:Class.Score
..cat:Alignments
..summary:Match score.
..signature:scoreMatch(object)
..param.object.type:Spec.Simple
..returns:Match score.
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
*/
template <typename TValue>
inline TValue &
scoreMatch(Score<TValue, Simple> & me)
{
	return me.data_match;
}
template <typename TValue>
inline TValue const &
scoreMatch(Score<TValue, Simple> const & me)
{
	return me.data_match;
}

/**.Function.scoreMismatch:
..class:Class.Score
..cat:Alignments
..summary:Mismatch score.
..signature:scoreMismatch(object)
..param.object.type:Spec.Simple
..returns:Mismatch score.
..see:Function.scoreMatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
*/
template <typename TValue>
inline TValue &
scoreMismatch(Score<TValue, Simple> & me)
{
	return me.data_mismatch;
}
template <typename TValue>
inline TValue const &
scoreMismatch(Score<TValue, Simple> const & me)
{
	return me.data_mismatch;
}
/**.Function.scoreGapExtend:
..class:Class.Score
..cat:Alignments
..summary:Score for extending gaps.
..signature:scoreGapExtend(object)
..param.object.type:Spec.Simple
..returns:Score for extending gaps.
..see:Function.scoreMismatch
..see:Function.scoreMatch
..see:Function.scoreGapOpen
*/
template <typename TValue>
inline TValue &
scoreGapExtend(Score<TValue, Simple> & me)
{
	return me.data_gap_extend;
}
template <typename TValue>
inline TValue const &
scoreGapExtend(Score<TValue, Simple> const & me)
{
	return me.data_gap_extend;
}
/**.Function.scoreGapOpen:
..class:Class.Score
..cat:Alignments
..summary:Score for opening a gap.
..signature:scoreGapOpen(object)
..param.object.type:Spec.Simple
..returns:Score for opening a gap.
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreMatch
*/
template <typename TValue>
inline TValue &
scoreGapOpen(Score<TValue, Simple> & me)
{
	return me.data_gap_open;
}
template <typename TValue>
inline TValue const &
scoreGapOpen(Score<TValue, Simple> const & me)
{
	return me.data_gap_open;
}

//////////////////////////////////////////////////////////////////////////////

// compute score of aligning two characters

template <typename TValue, typename T>
inline TValue
score(Score<TValue, Simple> & me,
	  T const & left,
	  T const & right)
{
	if (left == right) return scoreMatch(me);
	else return scoreMismatch(me);
}


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
