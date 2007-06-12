#ifndef SEQAN_HEADER_SCORE_MANHATTAN_H
#define SEQAN_HEADER_SCORE_MANHATTAN_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Score spec for manhattan distance

template <typename TValue>
class Score<TValue, Manhattan>
{
private:
	TValue data_match;
	TValue data_mismatch;
	TValue data_gap;

public:
	Score( TValue _match = 0, TValue _misalign = 1 ):
		data_match( _match ),
		data_mismatch( _misalign ),
		data_gap( _misalign )
	{
	}

	Score( TValue score ):
		data_match( 0 ),
		data_mismatch( score ),
		data_gap( score )
	{
	}

	Score(Score const & other):
		data_match(other.data_match),
		data_mismatch(other.data_mismatch),
		data_gap(other.data_gap)
	{
	}

	~Score()
	{
	}

	Score & operator = (Score const & other)
	{
		data_match = other.data_match;
		data_mismatch = other.data_mismatch;
		data_gap = other.data_gap;
		return *this;
	}

//____________________________________________________________________________

	friend inline TValue 
	scoreMatch(Score & me)
	{
		return me.data_match;
	}
	friend inline TValue const 
	scoreMatch(Score const & me)
	{
		return me.data_match;
	}

	friend inline TValue 
	scoreMismatch(Score & me)
	{
		return me.data_mismatch;
	}
	friend inline TValue const 
	scoreMismatch(Score const & me)
	{
		return me.data_mismatch;
	}

	friend inline TValue 
	scoreGapExtend(Score & me)
	{
		return me.data_gap;
	}
	friend inline TValue const &
	scoreGapExtend(Score const & me)
	{
		return me.data_gap;
	}

	template <typename T>
	friend inline T &
	scoreGapOpen(Score & me)
	{
		return me.data_gap;
	}
	friend inline TValue const &
	scoreGapOpen(Score const & me)
	{
		return me.data_gap;
	}

//____________________________________________________________________________
};
//////////////////////////////////////////////////////////////////////////////

//Shortcut:


//template< typename TValue >
//typedef typename Score< TValue, Manhattan > ManhattanScore;

//////////////////////////////////////////////////////////////////////////////

// compute score of aligning two characters

template <typename TValue, typename T>
inline TValue
score(Score<TValue, Manhattan> & me,
	  T const & left,
	  T const & right)
{
	if (left == right) return scoreMatch(me);
	else return scoreMismatch(me);
}


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
