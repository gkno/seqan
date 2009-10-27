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

#ifndef SEQAN_HEADER_STORE_INTERVALTREE_H
#define SEQAN_HEADER_STORE_INTERVALTREE_H
//#define DEBUG_TREE

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Create IntervallTreeStores
//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig>
inline void
createIntervalTreeStore(FragmentStore<TSpec, TConfig> & me, const bool &unknownO)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAnnotationStore 	TAnnotationStore;
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos 		TContigPos;
	typedef typename Value<TAnnotationStore>::Type 				TAnnotationStoreElement;
	typedef typename TAnnotationStoreElement::TId 				TId;
	typedef typename Iterator<TAnnotationStore>::Type 			TAnnotationIterator;
	
	typedef	typename FragmentStore<TSpec, TConfig>::TIntervalTreeStore 	TIntervalTreeStore;
	typedef typename Value<TIntervalTreeStore>::Type 			TIntervalTree;
	typedef typename TIntervalTree::TInterval 				TInterval;
	typedef 	 String<TInterval>					TIntervals;
	typedef typename Iterator<String<TIntervals> >::Type			TCIter;
	
	static const TId INVALID_ID = TAnnotationStoreElement::INVALID_ID;
	
	// get intervals for each contig (R- and F-strand):
	if (!empty(me.annotationStore) && !unknownO)
	{
		resize(me.intervalTreeStore_F, length(me.contigStore) );
		resize(me.intervalTreeStore_R, length(me.contigStore) );

		String<TIntervals> contigIntervals_F;
		String<TIntervals> contigIntervals_R;
		resize(contigIntervals_F, length(me.contigStore));
		resize(contigIntervals_R, length(me.contigStore));
		
		TAnnotationIterator itAnno = begin(me.annotationStore);
		TAnnotationIterator itAnnoEnd = end(me.annotationStore);
		TId beginPos;
		TId endPos;
		TInterval interval;
		for ( ; itAnno != itAnnoEnd; goNext(itAnno))
		{
			if (getValue(itAnno).contigId != INVALID_ID)
			{
				beginPos = getValue(itAnno).beginPos;
				endPos = getValue(itAnno).endPos;
				if (beginPos != INVALID_ID && beginPos <= endPos)
				{
					interval.i1 = beginPos;
					interval.i2 = endPos;
					interval.cargo = position(itAnno, me.annotationStore);
					appendValue(value(contigIntervals_F,  getValue(itAnno).contigId), interval, Generous());		
				}
				else if (beginPos != INVALID_ID  && beginPos > endPos)
				{
					interval.i1 = endPos;					
					interval.i2 = beginPos;
					interval.cargo = position(itAnno, me.annotationStore);
					appendValue(value(contigIntervals_R, getValue(itAnno).contigId), interval, Generous() );
				}
			}
		}
	
		// build trees for each contig and each strand:
		TCIter itF = begin(contigIntervals_F);
		TCIter itFEnd = end(contigIntervals_F);
		TCIter itR = begin(contigIntervals_R);
		for ( ; itF != itFEnd; goNext(itF), goNext(itR))
		{
			TIntervalTree intervalTree_F(getValue(itF)/*,RandomCenter()*/); // TODO  
			TIntervalTree intervalTree_R(getValue(itR)/*,RandomCenter()*/);
				
			assignValue(me.intervalTreeStore_F, position(itF, contigIntervals_F), intervalTree_F);  
			assignValue(me.intervalTreeStore_R, position(itR, contigIntervals_R), intervalTree_R);
		}
	}
	
	// if read orientation is not known:
	// get intervals for each contig:
	if (!empty(me.annotationStore) && unknownO)
	{
		resize(me.intervalTreeStore_F, length(me.contigStore) );
		clear(me.intervalTreeStore_R);

		String<TIntervals> contigIntervals;
		resize(contigIntervals, length(me.contigStore));
		TAnnotationIterator itAnno = begin(me.annotationStore);
		TAnnotationIterator itAnnoEnd = end(me.annotationStore);
		TId beginPos;
		TId endPos;
		TInterval interval;
		for ( ; itAnno != itAnnoEnd; goNext(itAnno))
		{
			if (getValue(itAnno).contigId != INVALID_ID)
			{
				beginPos = getValue(itAnno).beginPos;
				endPos = getValue(itAnno).endPos;
				if (beginPos != INVALID_ID)
				{
					if (beginPos <= endPos)
					{
						interval.i1 = beginPos;
						interval.i2 = endPos;
					}
					else
					{
						interval.i1 = endPos;
						interval.i2 = beginPos;
					}
					interval.cargo = position(itAnno, me.annotationStore);
					appendValue(value(contigIntervals,  getValue(itAnno).contigId), interval, Generous());		
				}
			}
		}
		// build trees for each contig:
		TCIter itC = begin(contigIntervals);
		TCIter itCEnd = end(contigIntervals);
		for ( ; itC != itCEnd; goNext(itC))
		{
			TIntervalTree intervalTree(getValue(itC)/*,RandomCenter()*/); // TODO
				
			assignValue(me.intervalTreeStore_F, position(itC, contigIntervals), intervalTree); 
		}
	}
}





//////////////////////////////////////////////////////////////////////////////
// find intervals containing one interval
//////////////////////////////////////////////////////////////////////////////


template<typename TIntervalTree, typename TInterval, typename TCargo>
inline void
findIntervalsForInterval(String<TCargo> & result, TIntervalTree & intervalTree, TInterval & interval, const unsigned & offsetInterval)
{
	String<TCargo> result1;
	String<TCargo> result2;
	
	findIntervals(intervalTree.g, intervalTree.pm, interval.i1 + offsetInterval, result1);
	findIntervals(intervalTree.g, intervalTree.pm, interval.i2 - offsetInterval, result2);
	
	interSec(result, result1, result2); 
	
}



//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
