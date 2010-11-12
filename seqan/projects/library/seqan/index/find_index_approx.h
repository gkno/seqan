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

#ifndef SEQAN_HEADER_FIND_INDEX_APPROX_H
#define SEQAN_HEADER_FIND_INDEX_APPROX_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TSAValue1, typename TSAValue2>
struct AddResultsFunctor
{
	String<Pair<TSAValue1, String<TSAValue2> > > results;
    
    template <typename TIter1, typename TIter2>
    void operator() (TIter1 &iter1, TIter2 &iter2)
    {
		unsigned ofs = length(results);
		resize(results, ofs + countOccurrences(iter1));
		for (unsigned i = 0; i < countOccurrences(iter1); ++i)
		{
			results[ofs + i].i1 = getOccurrences(iter1)[i];
			results[ofs + i].i2 = getOccurrences(iter2);
		}
    }
};

template <
	bool enumerateA,
	bool enumerateB,
	typename TOnFoundFunctor, 
	typename TTreeIteratorA, 
	typename TIterPosA, 
	typename TTreeIteratorB, 
	typename TIterPosB, 
	typename TErrors >
inline void 
_approximateTreeSearch(
	TOnFoundFunctor	&onFoundFunctor, 
	TTreeIteratorA	iterA, 
	TIterPosA		iterPosA, 
	TTreeIteratorB	iterB_, 
	TIterPosB		iterPosB, 
	TErrors			errorsLeft)
{
	if (enumerateA && !goDown(iterA)) 
	{
        onFoundFunctor(iterA, iterB_);
		return;
	}
	if (enumerateB && !goDown(iterB_)) return;
	
	do 
	{
		TTreeIteratorB iterB = iterB_;
		do 
		{
			TErrors e = errorsLeft;
			TIterPosA ipA = iterPosA;
			TIterPosB ipB = iterPosB;
			
			if (ipB == repLength(iterB)) continue;
			
			while (true)
			{
				if (ipA == repLength(iterA))
				{
					if (ipB == repLength(iterB))
						_approximateTreeSearch<true,true>(onFoundFunctor, iterA, ipA, iterB, ipB, e);
					else
						_approximateTreeSearch<true,false>(onFoundFunctor, iterA, ipA, iterB, ipB, e);
					break;
				} 
				else if (ipB == repLength(iterB))
				{
					_approximateTreeSearch<false,true>(onFoundFunctor, iterA, ipA, iterB, ipB, e);
					break;
				}

				if (representative(iterA)[ipA] != representative(iterB)[ipB])
					if (e-- == 0) break;
				
				++ipA;
				++ipB;
			}
		} while (enumerateB && goRight(iterB));
	} while (enumerateA && goRight(iterA));
}


template <typename TSAValue>
struct AddSingleResultsFunctor
{
	String<TSAValue> results;
    
    template <typename TPattern, typename TIter>
    void operator() (TPattern &pattern, TIter &iter)
    {
        append(results, getOccurrences(iter));
    }
};


template <
	typename TOnFoundFunctor, 
	typename TString, 
	typename TStringPos, 
	typename TTreeIterator, 
	typename TIterPos, 
	typename TErrors >
inline void 
_approximateStringSearch(
	TOnFoundFunctor	&onFoundFunctor, 
	TString const &string, 
	TStringPos stringPos, 
	TTreeIterator iter, 
	TIterPos iterPos, 
	TErrors errorsLeft)
{
	if (errorsLeft == 0)
	{
		if (goDown(iter, suffix(string, stringPos)))
			onFoundFunctor(string, iter);
		return;
	}
	
	if (!goDown(iter)) return;
	do 
	{
		TErrors e = errorsLeft;
		TStringPos sp = stringPos;
		TIterPos ip = iterPos;
		
		if (ip == repLength(iter)) continue;
		
		while (true)
		{
			if (representative(iter)[ip] != string[sp])
				if (e-- == 0) break;
			
			if (++sp == length(string))
			{
                onFoundFunctor(string, iter);
				break;
			}
			
			if (++ip == repLength(iter))
			{
				_approximateStringSearch(onFoundFunctor, string, sp, iter, ip, e);
				break;
			}
		}
	} while (goRight(iter));
}

template <
    typename TOnFoundFunctor, 
    typename TString, 
    typename TTreeIterator, 
    typename TErrors>
inline void 
approximateStringSearch(
    TOnFoundFunctor &onFoundFunctor, 
    TString const &string, 
    TTreeIterator &iter, 
    TErrors errorsLeft)
{
	if (length(string) <= errorsLeft)
	{
        onFoundFunctor(string, iter);
		return;
	}	
	_approximateStringSearch(onFoundFunctor, string, 0u, iter, repLength(iter), errorsLeft);
}

template <
    typename TOnFoundFunctor, 
    typename TTreeIteratorA, 
    typename TTreeIteratorB, 
    typename TErrors>
inline void 
approximateTreeSearch(
    TOnFoundFunctor &onFoundFunctor, 
    TTreeIteratorA const &iterA, 
    TTreeIteratorB const &iterB, 
    TErrors errorsLeft)
{
	_approximateTreeSearch<true,true>(
        onFoundFunctor, 
        iterA, 
        repLength(iterA), 
        iterB, 
        repLength(iterB), 
        errorsLeft);
}

}
#endif
