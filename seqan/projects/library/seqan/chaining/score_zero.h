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

 ============================================================================
  $Id$
 ==========================================================================*/

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

	Score(Score const & )
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
	scoreMatch(Score &)
	{
		return 0;
	}
	friend inline TValue const 
	scoreMatch(Score const &)
	{
		return 0;
	}

	friend inline TValue 
	scoreMismatch(Score & me)
	{
		return 0;
	}
	friend inline TValue const 
	scoreMismatch(Score const &)
	{
		return 0;
	}

	friend inline TValue 
	scoreGap(Score &)
	{
		return 0;
	}
	friend inline TValue const 
	scoreGap(Score const &)
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
