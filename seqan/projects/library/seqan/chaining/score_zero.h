#ifndef SEQAN_HEADER_SCORE_ZERO_H
#define SEQAN_HEADER_SCORE_ZERO_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, Zero>
{
private:

public:
	Score()
	{
	}

	Score(Score const & other)
	{
	}

	~Score()
	{
	}

	Score & operator = (Score const & other)
	{
		if( this == & other )
			return *this;
		return *this;
	}

//____________________________________________________________________________

	friend inline TValue 
	scoreMatch(Score & me)
	{
		return 0;
	}
	friend inline TValue const 
	scoreMatch(Score const & me)
	{
		return 0;
	}

	friend inline TValue 
	scoreMismatch(Score & me)
	{
		return 0;
	}
	friend inline TValue const 
	scoreMismatch(Score const & me)
	{
		return 0;
	}

	friend inline TValue 
	scoreGapExtend(Score & me)
	{
		return 0;
	}
	friend inline TValue const 
	scoreGapExtend(Score const & me)
	{
		return 0;
	}

	template <typename T>
	friend inline T 
	scoreGapOpen(Score<T, Zero> & me)
	{
		return 0;
	}
	friend inline TValue const 
	scoreGapOpen(Score<TValue, Zero> const & me)
	{
		return 0;
	}

//____________________________________________________________________________
};
//////////////////////////////////////////////////////////////////////////////

//Shortcut:

//template< typename TValue >
//typedef typename Score<TValue, Zero> ZeroScore;

//////////////////////////////////////////////////////////////////////////////

// compute score of aligning two characters

template <typename TValue, typename T>
inline TValue
score(Score<TValue, Zero> const & me,
	  T const & left,
	  T const & right)
{
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
//compute score for chaining two fragments

template <typename TValue, typename TFragment>
inline TValue
scoreChainGap(Score<TValue, Zero> const &,
			  TFragment &,
			  TFragment &)
{
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
