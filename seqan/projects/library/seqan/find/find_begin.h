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
  $Id: $
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_BEGIN_H
#define SEQAN_HEADER_FIND_BEGIN_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Forward declarations of Pattern Specs that could be used prefix search

//see find_score.h
template <typename TScore, typename TSpec, typename TFindBeginPatternSpec>
struct DPSearch;

//see finder_myers_ukkonen.h
//template <typename TSpec, typename TFindBeginPatternSpec>
//struct Myers;


//____________________________________________________________________________

template <typename TScore = EditDistanceScore>
struct DefaultFindBeginPatternSpec
{
	typedef DPSearch<TScore, FindPrefix, void> Type;
};
//template <>
//struct DefaultFindBeginPatternSpec <EditDistanceScore>
//{
//	typedef Myers<FindPrefix, void> Type;
//};

//____________________________________________________________________________
//must be implemented for all patterns that do approximate infix finding

template <typename TPattern>
struct FindBeginPatternSpec 
{
	typedef void Type; //void means: no find begin (see _FindBegin)
};

//____________________________________________________________________________
// Metafunction FindBeginPattern: pattern used to find the begin of an approximate match

template <typename TPattern>
struct FindBeginPattern
{
	typedef typename Needle<TPattern>::Type TNeedle;
	typedef ModifiedString<TNeedle, ModReverse> TReverseNeedle;
	typedef typename FindBeginPatternSpec<TPattern>::Type TFindBeginPatternSpec;
	typedef Pattern<TReverseNeedle, TFindBeginPatternSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////
// Base class for approximate finders that manages findBegin stuff

template <typename TPattern, typename TFindBeginPatternSpec = typename FindBeginPatternSpec<TPattern>::Type >
struct _FindBegin
{
private:
	typedef typename FindBeginPattern<TPattern>::Type TFindBeginPattern;

public:
	TFindBeginPattern data_findBeginPattern;
};

template <typename TPattern>
struct _FindBegin <TPattern, void>
{
//need no findBegin if FindBeginPatternSpec is void
};



//////////////////////////////////////////////////////////////////////////////

template <typename TFindBeginPatternSpec>
struct _FindBegin_Impl
{
	template <typename TPattern>
	static inline void
	_findBeginInit(TPattern & pattern)
	{
//		setNeedle(pattern.data_findBeginPattern, reverseString(needle(pattern)));
		setHost(needle(pattern.data_findBeginPattern), needle(pattern));
		setScoringScheme(pattern.data_findBeginPattern, scoringScheme(pattern));
	}
//____________________________________________________________________________

	template <typename TFinder, typename TPattern, typename TLimit>
	static inline bool
	findBegin(TFinder & finder, TPattern & pattern, TLimit limit)
	{
		typedef typename FindBeginPattern<TPattern>::Type TFindBeginPattern;
		typedef typename Haystack<TFinder>::Type THaystack;
		typedef ModifiedString<THaystack, ModReverse> TReverseHaystack;
		typedef Finder<TReverseHaystack> TBeginFinder;
		typedef typename Position<THaystack>::Type TPosition;

		TFindBeginPattern & find_begin_pattern = pattern.data_findBeginPattern;
		setScoreLimit(find_begin_pattern, limit);

		//build begin_finder
		TBeginFinder begin_finder;
		THaystack & hayst = haystack(finder);
		setContainer(host(hostIterator(begin_finder)), hayst);
		TPosition begin_finder_beginPosition = position(finder);
		TPosition begin_finder_position;

		if (!finder._beginFind_called)
		{//start finding
			finder._beginFind_called = true;
			_setFinderLength(finder, 0);
			clear(begin_finder);
			begin_finder_position = begin_finder_beginPosition;
		}
		else
		{//resume finding
			_finderSetNonEmpty(begin_finder);
			SEQAN_ASSERT(length(finder) > 0)
			begin_finder_position = endPosition(finder) - length(finder);
		}
		setPosition(host(hostIterator(begin_finder)), begin_finder_position);
		_setFinderEnd(begin_finder);
		_setFinderLength(begin_finder, length(finder));

		bool begin_found = find(begin_finder, find_begin_pattern);
		if (begin_found)
		{//new begin found: report in finder
			_setFinderLength(finder, begin_finder_beginPosition - position(host(hostIterator(begin_finder)))+1);
		}
		return begin_found;
	}
	template <typename TFinder, typename TPattern>
	static inline bool
	findBegin(TFinder & finder, TPattern & pattern)
	{
		return findBegin(finder, pattern, scoreLimit(pattern));
	}

//____________________________________________________________________________

	template <typename TPattern>
	static inline typename Value<typename ScoringScheme<TPattern>::Type>::Type
	getBeginScore(TPattern & pattern)
	{
		return getScore(pattern.data_findBeginPattern);
	}
};

//default implementation of findBegin emulates the behaviour 
//of findBegin for approximate string matching algorithms
//for exact algorithms:
//the first call returns true ("begin position found")
//the next calls return false ("no further begin position found")
//no searching needed, since the exact algorithms set the 
//match length manually via "_setFinderLength" (it's the needle length)

template <>
struct _FindBegin_Impl<void>
{
	template <typename TPattern>
	static inline void
	_findBeginInit(TPattern &)
	{
	}

	template <typename TFinder, typename TPattern, typename TLimit>
	static inline bool
	findBegin(TFinder & finder, TPattern &, TLimit)
	{
		if (!finder._beginFind_called)
		{
			finder._beginFind_called = true;
			return true;
		}
		return false;
	}
	template <typename TFinder, typename TPattern>
	static inline bool
	findBegin(TFinder & finder, TPattern & pattern)
	{
		return findBegin(finder, pattern, 0);
	}

	template <typename TPattern>
	static inline typename Value<typename ScoringScheme<TPattern>::Type>::Type
	getBeginScore(TPattern & pattern)
	{
		return getScore(pattern);
	}
};


//////////////////////////////////////////////////////////////////////////////

//initialize pattern of begin finding
template <typename TPattern>
inline void
_findBeginInit(TPattern & pattern)
{
	typedef typename FindBeginPatternSpec<TPattern>::Type TFindBeginPatternSpec;
	return _FindBegin_Impl<TFindBeginPatternSpec>::_findBeginInit(pattern);
}

//find begin main interface 
template <typename TFinder, typename TPattern>
inline bool
findBegin(TFinder & finder,
		  TPattern & pattern)
{
	typedef typename FindBeginPatternSpec<TPattern>::Type TFindBeginPatternSpec;
	return _FindBegin_Impl<TFindBeginPatternSpec>::findBegin(finder, pattern);
}
template <typename TFinder, typename TPattern, typename TLimit>
inline bool
findBegin(TFinder & finder,
		  TPattern & pattern,
		  TLimit limit)
{
	typedef typename FindBeginPatternSpec<TPattern>::Type TFindBeginPatternSpec;
	return _FindBegin_Impl<TFindBeginPatternSpec>::findBegin(finder, pattern, limit);
}

template <typename TPattern>
inline typename Value<typename ScoringScheme<TPattern>::Type>::Type
getBeginScore(TPattern & pattern)
{
	typedef typename FindBeginPatternSpec<TPattern>::Type TFindBeginPatternSpec;
	return _FindBegin_Impl<TFindBeginPatternSpec>::getBeginScore(pattern);
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
