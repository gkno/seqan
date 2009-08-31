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

#include <iostream>

#ifndef SEQAN_HEADER_STORE_IO_SAM_OUT_H
#define SEQAN_HEADER_STORE_IO_SAM_OUT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// write
    
    template<typename TFile, typename TSpec, typename TConfig>
    inline void write(TFile & target,
                      FragmentStore<TSpec, TConfig> & store,
                      SAM)
    {
        // write header
        
        // write aligments
        _writeAlignments(target, store, SAM());
    }
    
//////////////////////////////////////////////////////////////////////////////
// _writeAlignments

    template<typename TFile, typename TSpec, typename TConfig>
	inline void _writeAlignments(TFile & target,
                                 FragmentStore<TSpec, TConfig> & store,
                                 SAM)
    {
		typedef FragmentStore<TSpec, TConfig>							TFragmentStore;

		typedef typename TFragmentStore::TReadStore						TReadStore;
		typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
		typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
		typedef typename TFragmentStore::TContigStore					TContigStore;

        typedef typename Value<TReadStore>::Type						TRead;
        typedef typename Value<TContigStore>::Type						TContig;
        typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;

        typedef typename Value<TReadSeqStore>::Type						TReadSeq;
		typedef typename TContig::TContigSeq							TContigSeq;
        typedef typename Iterator<TAlignedReadStore, Standard>::Type	TAlignIter;
        typedef typename Id<TAlignedRead>::Type							TId;

		typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >	TReadGaps;
		typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;

		String<int> mateIndex;	// store outer library size for each pair match (indexed by pairMatchId)
		calculateMateIndices(mateIndex, store);
		
        TAlignIter it = begin(store.alignedReadStore, Standard());
        TAlignIter itEnd = end(store.alignedReadStore, Standard());
		TAlignIter mit;
		CharString cigar;

        for(; it != itEnd; ++it)
		{
            TId alignedId = (*it).id;
			TId readId = (*it).readId;
			TId mateIdx = mateIndex[2*(*it).pairMatchId + getMateNo(store, (*it).readId)];

			__int64 mpos = 0;
			int isize = 0;
			unsigned short flag = 0;

			if ((*it).beginPos > (*it).endPos)
				flag |= 0x0010;
			
			if (mateIdx < length(store.alignedReadStore))
			{
				mit = begin(store.alignedReadStore, Standard()) + mateIdx;
				mpos = min((*mit).beginPos, (*mit).endPos);
				isize = (*mit).beginPos - (*it).beginPos;
				flag |= 0x0002;
				if ((*mit).beginPos > (*mit).endPos)
					flag |= 0x0020;				
			}
			unsigned char mateNo = getMateNo(store, readId);
			if (mateNo == 0) flag |= 0x0040;
			if (mateNo == 1) flag |= 0x0080;

			if (readId < length(store.readStore))
			{
				TRead &read = store.readStore[readId];
				if (read.matePairId != TRead::INVALID_ID)
					flag |= 0x0001;
			}

			// qname
			if (readId < length(store.readNameStore))
				_streamWrite(target, store.readNameStore[readId]);
            _streamPut(target, '\t');
            
            // flag
            _streamPutInt(target, flag);
            _streamPut(target, '\t');
            
			// rname
			if ((*it).contigId < length(store.contigNameStore))
				_streamWrite(target, store.contigNameStore[(*it).contigId]);
            _streamPut(target, '\t');
            
			// begin
            _streamPutInt(target, min((*it).beginPos, (*it).endPos));
            _streamPut(target, '\t');
            
			if (alignedId < length(store.alignQualityStore))
				_streamPutInt(target, store.alignQualityStore[alignedId].score);
            _streamPut(target, '\t');
            
            // cigar
			TContigGaps	contigGaps(store.contigStore[(*it).contigId].seq, store.contigStore[(*it).contigId].gaps);
			TReadGaps	readGaps(store.readSeqStore[readId], store.alignedReadStore[alignedId].gaps);
			setBeginPosition(contigGaps, (*it).beginPos);
			setEndPosition(contigGaps, (*it).endPos);
			getCigarString(cigar, contigGaps, readGaps);
			_streamWrite(target, cigar);
            _streamPut(target, '\t');
            
            // mrnm
			if ((*it).contigId == (*mit).contigId)
				_streamWrite(target, '=');
			else
				if ((*mit).contigId < length(store.contigNameStore))
					_streamWrite(target, store.contigNameStore[(*mit).contigId]);
            _streamPut(target, '\t');
            
            // mpos
			_streamPutInt(target, mpos);
            _streamPut(target, '\t');
            
            // isize
			_streamPutInt(target, isize);
            _streamPut(target, '\t');

            
			if (readId < length(store.readSeqStore))
				_streamWrite(target, store.readSeqStore[readId]);
            _streamPut(target, '\t');
            
            // qual
            _streamPut(target, '\t');
            
			if (alignedId < length(store.alignedReadTagStore))
				_streamWrite(target, store.alignedReadTagStore[alignedId]);
            
            _streamPut(target, '\n');
        }
        
    }
    
    
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
