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

#ifndef SEQAN_HEADER_REPEAT_BASE_H
#define SEQAN_HEADER_REPEAT_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	template <typename TPos, typename TPeriod>
	struct Repeat {
		TPos		beginPosition;
		TPos		endPosition;
		TPeriod		period;
	};

	template <typename TPos, typename TPeriod>
	struct Value< Repeat<TPos, TPeriod> > {
		typedef TPos Type;
	};

	template <typename TPos, typename TPeriod>
	struct Size< Repeat<TPos, TPeriod> > {
		typedef TPeriod Type;
	};



	template <typename TSize>
	struct RepeatFinderParams {
		TSize minRepeatLen;
		TSize maxPeriod;
	};

	// custom TSpec for our customized wotd-Index
	struct TRepeatFinder;

	template <typename TText>
	struct Cargo<Index<TText, Index_Wotd<TRepeatFinder> > > 
	{
		typedef Index<TText, Index_Wotd<TRepeatFinder> >	TIndex;
		typedef typename Size<TIndex>::Type					TSize;
		typedef RepeatFinderParams<TSize>					Type;
	};


	// node predicate
	template <typename TText, typename TSpec>
	bool nodePredicate(Iter<Index<TText, Index_Wotd<TRepeatFinder> >, TSpec> &it) 
	{
		return countOccurrences(it) * repLength(it) >= cargo(container(it)).minRepeatLen;
	}

	// monotonic hull
	template <typename TText, typename TSpec>
	bool nodeHullPredicate(Iter<Index<TText, Index_Wotd<TRepeatFinder> >, TSpec> &it) 
	{
		return repLength(it) <= cargo(container(it)).maxPeriod;
	}

	template <typename TPos>
	struct _RepeatLess : public ::std::binary_function<TPos, TPos, bool>
	{
		// key less
		inline bool operator() (TPos const &a, TPos const &b) {
			return posLess(a, b);
		}
	};

	// main function
	template <typename TRepeatStore, typename TText, typename TRepeatSize, typename TPeriodSize>
	void findRepeats(TRepeatStore &repString, TText const &text, TRepeatSize minRepeatLen, TPeriodSize maxPeriod) 
	{
		typedef Index<TText, Index_Wotd<TRepeatFinder> >					TIndex;
		typedef typename Size<TIndex>::Type									TSize;
		typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type	TNodeIterator;
		typedef typename Fibre<TIndex, Fibre_SA>::Type const				TSA;
		typedef typename Infix<TSA>::Type									TOccString;
		typedef typename Iterator<TOccString>::Type							TOccIterator;

		typedef typename Value<TRepeatStore>::Type							TRepeat;
		typedef typename Value<TOccString>::Type							TOcc;

		typedef ::std::map<TOcc,TRepeat,_RepeatLess<TOcc> >					TRepeatList;

		TIndex		index(text);
		TRepeatList list;

		// set repeat finder parameters
		cargo(index).minRepeatLen = minRepeatLen;
		cargo(index).maxPeriod = maxPeriod;

		TNodeIterator nodeIt(index);
		TOccIterator itA, itB, itRepBegin, itEnd;
		TRepeat rep;
		while (!atEnd(nodeIt)) 
		{
			// get occurrences
			TOccString occ = getOccurrences(nodeIt);
			itA = begin(occ, Standard());
			itEnd = end(occ, Standard());
			itRepBegin = itB = itA;

			// get representative/repeat length
			TSize repLen = repLength(nodeIt);
			if (repLen > 0) 
			{
				TSize minRepeats = (minRepeatLen + repLen - 1) / repLen - 1;

				for (++itB; itB != itEnd; ++itB)
				{
					if (posSub(*itB, *itA) != repLen || getSeqNo(*itA) != getSeqNo(*itB))
					{
						if ((TSize)(itA - itRepBegin) >= minRepeats) 
						{
							// insert repeat
							rep.beginPosition = *itRepBegin;
							rep.endPosition = posAdd(*itA, repLen);
							rep.period = repLen;
//							::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
							list.insert(::std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
						}
						itRepBegin = itB;
					}
					itA = itB;
				}

				if ((TSize)(itA - itRepBegin) >= minRepeats) 
				{
					// insert repeat
					rep.beginPosition = *itRepBegin;
					rep.endPosition = posAdd(*itA, repLen);
					rep.period = repLen;
//					::std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<::std::endl;
					list.insert(::std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
				}
			}

			goNext(nodeIt);
		}

		// copy low-complex regions to result string
		resize(repString, list.size());
		typename TRepeatList::const_iterator lit = list.begin();
		typename TRepeatList::const_iterator litEnd = list.begin();
		for (TSize i = 0; lit != litEnd; ++lit, ++i)
			repString[i] = (*lit).second;
	}

}	// namespace seqan

#endif
