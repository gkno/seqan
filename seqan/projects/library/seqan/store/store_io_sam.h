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

#ifndef SEQAN_HEADER_STORE_IO_SAM_H
#define SEQAN_HEADER_STORE_IO_SAM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.SAM:
	SAM mapping file.
*/
//struct TagSAM_;
//typedef Tag<TagSAM_> const SAM;

    
//////////////////////////////////////////////////////////////////////////////
// Cigar struct
//////////////////////////////////////////////////////////////////////////////
    
    template <typename _TOperation = char, typename _TCount = unsigned>
	struct CigarElement
	{
		typedef _TOperation TOperation;
		typedef _TCount		TCount;

		TOperation			operation;
		TCount				count;
		
		CigarElement(TOperation o, TCount c):
			operation(o),
			count(c) {}
	};
	
//////////////////////////////////////////////////////////////////////////////
// _getClippedLength
    
    template <typename TCigarString, typename TNum>
    inline void _getClippedLength(TCigarString const & cigar, TNum & sum)
    {
        typedef typename Iterator<TCigarString, Standard>::Type TCigarIter;
        
        TCigarIter it = begin(cigar, Standard());
        TCigarIter itEnd = end(cigar, Standard());
        
        sum = 0;        
        for (; it != itEnd; ++it)
            if (getValue(it).operation != 'S' && getValue(it).operation != 'H')
                sum += getValue(it).count;
    }

//////////////////////////////////////////////////////////////////////////////
// cigarToGapAnchorRead
    template<typename TCigarString, typename TGaps>
    inline void cigarToGapAnchorRead(TCigarString const & cigar, TGaps & gaps)
    {
		typename Iterator<TGaps>::Type it = begin(gaps);
		for (unsigned i = 0; i < length(cigar); ++i)
		{
			switch (cigar[i].operation)
			{
				case 'P':
				case 'D':
					insertGaps(it, cigar[i].count);
				case 'M':
				case 'I':
					it += cigar[i].count;
			}
		}
	}

    template<typename TCigarString, typename TGaps>
    inline void cigarToGapAnchorContig(TCigarString const & cigar, TGaps & gaps)
    {
		typename Iterator<TGaps>::Type it = begin(gaps);
		for (unsigned i = 0; i < length(cigar); ++i)
		{
			switch (cigar[i].operation)
			{
				case 'P':
				case 'I':
					insertGaps(it, cigar[i].count);
				case 'M':
				case 'D':
					it += cigar[i].count;
			}
		}
	}
/*
    template<typename TChar, typename TNum, typename TGapAnchor>
    inline void cigarToGapAnchorRead(Cigar<TChar, TNum> const & cigar, String<TGapAnchor> & gaps)
    {
        typedef typename Iterator< String<TNum> >::Type TNumIter;
        typedef typename Iterator< String<TChar> >::Type TCharIter;
        
        // Iterators on cigar
        TNumIter it_c = begin(cigar.operationCount);
        TCharIter it_t = begin(cigar.operationType);
        
        // delete information in gap anchor string
        gaps = String<TGapAnchor>();
        
        // boolean that keeps track on type of the last cigar operation
        // is true if match/mismatch or insertion
        // if false if deletion, padding or clipped
        bool was_mi = true;
        
        // positions in the un-gapped and gapped sequence
        int pos_seq = 0;
        int pos_gapped = 0;
        
        // if the CIGAR is empty
        if(length(cigar.operationCount) == 0){
            return;
        }
        
        // check if there is a (soft) clipped sequence in the begining
        if(value(it_t) == 'S' || value(it_t) == 's'){
            
            // skip first 'count' characters
            pos_seq += value(it_c);
            
            // insert in GapAnchor
            TGapAnchor anchor = TGapAnchor(pos_seq, pos_gapped);
            append(gaps, anchor);
        }
        
        while(it_c != end(cigar.operationCount)){
            
            // If the operation type is match/mismatch or a Deletion in the reference there is no gap in the read sequence.
            if(value(it_t) == 'M' || value(it_t) == 'I' || value(it_t) == 'm' || value(it_t) == 'i'){
                
                // If this is the first m/i operation after a d/p operation (first one after a gap)
                if(!was_mi){
                    
                    // insert in gap anchor
                    TGapAnchor anchor = TGapAnchor(pos_seq, pos_gapped);
                    append(gaps, anchor);
                    
                    // switch operation type to match/mismatch or insertion
                    was_mi = true;
                }
                
                // Increment positions in the gapped and un-gapped sequnece
                pos_gapped += value(it_c);
                pos_seq += value(it_c);
            }
            
            // If the operation type is insertion or skipped in the reference or padding there is a gap in the read sequence.
            if(value(it_t) == 'D' || value(it_t) == 'P' || value(it_t) == 'N' || value(it_t) == 'd' || value(it_t) == 'p' || value(it_t) == 'n'){

                // switch operation type to deletion or padding
                was_mi = false;
                
                // Increment only the position in the gapped sequence
                pos_gapped += value(it_c);
            }
            
            // Iterate.
            ++it_c;
            ++it_t;
        }
        
        // following (soft) klipped characters are encode by the exceeding length of the sequence
    }
*/
    //////////////////////////////////////////////////////////////////////////////
    // cigarToContigGaps
    
    template <typename TCigarString, typename TPos, typename TGapString>
    inline void 
	cigarToContigGaps (
		TCigarString const & cigar, 
		TPos beginPos, 
		bool reverse, 
		TGapString & gaps)
    {
        typedef typename Iterator<TCigarString, Standard>::Type TCigarIter;
        typedef typename Value<TGapString>::Type TPair;
        
        // Iterators on cigar
        TCigarIter it;
		TCigarIter itEnd;        
        int direction;
        
        // if the CIGAR is empty
        if (empty(cigar))
            return;
        
        if (reverse)
		{
            it = end(cigar, Standard()) - 1;
            itEnd = begin(cigar, Standard()) - 1;
            direction = -1;
        } else
		{
            it = begin(cigar, Standard());
            itEnd = end(cigar, Standard());
            direction = 1;
        }
                
        // positions and length of the current gap
        unsigned gapLength = 0;
		// the gap is inserted after pos later
		// here we insert gaps before pos
        TPos pos = beginPos - 1;
        
        // check if there is a (soft) clipped sequence in the begining
        if (getValue(it).operation == 'S')
		{ 
            // skip first 'count' characters
            pos += getValue(it).count;
        }
        
			std::cout<<std::endl;
        for (; it != itEnd; it += direction)
		{
			char op = getValue(it).operation;
			std::cout<< op<<':'<<getValue(it).count<<"  ";

            // If the operation type is deletion or padding there is a gap in the reference sequence.
            if (op == 'P' || op == 'I')
			{                
                // Increment only the gap length
                gapLength += getValue(it).count;
            } else
			{
                // If this is the first i/p/n operation after a m/d operation (first one after a gap)
                if (gapLength != 0)
				{
                    // insert in gap pair
					append(gaps, TPair(pos, gapLength), Generous());
					std::cout<<" p:"<<pos<<" l:"<<gapLength<<" r:"<<reverse;
                    
                    // set gap length back to 0
                    gapLength = 0;
                }
				
				// Increment position
				pos += getValue(it).count;
			}
        }
        
        // following (soft) klipped characters are encode by the exceeding length of the sequence
		if (gapLength != 0)
		{
			// insert in gap pair
			append(gaps, TPair(pos, gapLength), Generous());
		}
    }
    
//////////////////////////////////////////////////////////////////////////////
// appendAlignment
    
    template<typename TSpec, typename TConfig, typename TId, typename TPos, typename TGaps>
    inline typename Size<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type
	appendAlignment (
		FragmentStore<TSpec, TConfig> & fragStore, 
		TId readId, 
		TId contigId, 
		TPos beginPos, 
		TPos endPos, 
		TGaps const & gaps)
	{
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
        
        TId id = length(fragStore.alignedReadStore);
        TAlignedElement alignedElem = TAlignedElement(id, readId, contigId, beginPos, endPos, gaps);
        append(fragStore.alignedReadStore, alignedElem);
		
		return id;
    }
    
//////////////////////////////////////////////////////////////////////////////
// _appendRead
// 
// adds a new entry to the read store if neccessary. Otherwise it writes the 
// correct Id in the variable using to qname to identify it
// If needed a mate pair entry is created
    
    template<typename TSpec, typename TConfig, typename TId, typename TName, typename TString, typename TFlag, typename TContext>
    inline void 
	_appendRead (
		FragmentStore<TSpec, TConfig> & fragStore, 
		TId & readId, 
		TName & qname,
		TString & readSeq,
		TFlag & flag,
		TContext & context)
	{
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;

		// search for readId by name
        if (_getIdByName(fragStore.readNameStore, qname, readId, fragStore.readNameStoreCache))
		{
			if ((flag & 1) == 1)
			{
				// if the read is in the store and paired
				// check the mate pair store if it is the same mate of the pair
				// assuming that only one flag 0x040 or 0x0080 is 1
				int inPair = 1 - ((flag & 0x40) >> 6);	// bit 7 is set => inPair = 0
														// else inPair = 1 (even if bits 6 and 7 are not set)
				
				TId matePairId = fragStore.readStore[readId].matePairId;
				if (matePairId != TMatePairElement::INVALID_ID)
				{
					readId = fragStore.matePairStore[matePairId].readId[inPair];
					if (readId == TMatePairElement::INVALID_ID)
					{
						// create new entry in read and read name store
						// set sequence and mate pair ID in new read store element
						readId = appendRead(fragStore, readSeq, matePairId);
						// add the identifier to the read name store
						appendName(fragStore.readNameStore, qname, fragStore.readNameStoreCache);
						// set the ID in the mate pair store
						fragStore.matePairStore[matePairId].readId[inPair] = readId;
						return;
					}
				}
			} else 
				return;
        }

		// if the read name is not in the store
		// create new entry in read and read name store
		readId = length(fragStore.readStore);

		// if the read is paired
		if ((flag & 1) == 1)
		{
			TMatePairElement mateElem;
			// set the first or second read ID in the mate pair element
			TId matePairId = length(fragStore.matePairStore);
			mateElem.readId[(flag & 0x80) >> 7] = readId;
			// get a new mate pair ID and add the new mate pair element
			appendValue(fragStore.matePairStore, mateElem);
			// set the new mate pair ID in the read element
			appendRead(fragStore, readSeq, matePairId);
		} 
		// if read is not paired
		else
			appendRead(fragStore, readSeq);
		
		appendName(fragStore.readNameStore, qname, fragStore.readNameStoreCache);
    }
    
//////////////////////////////////////////////////////////////////////////////
// _appendContig
// 
// adds a new entry to the read store if neccessary. Otherwise it writes the 
// correct Id in the variable using to qname to identify it
// If needed a mate pair entry is created
    
    template<typename TSpec, typename TConfig, typename TId, typename TName>
    inline void 
	_appendContig (
		FragmentStore<TSpec, TConfig> & fragStore, 
		TId & contigId, 
		TName & rName)
    {
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
        
        if (!_getIdByName(fragStore.contigNameStore, rName, contigId, fragStore.contigNameStoreCache))
        {
			// if the contig is not in the store yet
            // set the ID on the last entry after appending
            contigId = length(fragStore.contigStore);
            // append contig store
            appendName(fragStore.contigNameStore, rName, fragStore.contigNameStoreCache);
            appendValue(fragStore.contigStore, TContigElement());
        }
    }

//////////////////////////////////////////////////////////////////////////////
// _generatePairMatchIds
//
// shift is the first ID for which mateInfo contains information about a aligned read
// the function generates new pair match IDs for these entries. It uses the ID of the first found mate for this.
    
	template <typename TPos, typename TId>
	struct _MatchMateInfo
	{
		TId		readId;
		TId		contigId;
		TId		pairMatchId;
		TPos	beginPos;
	};
    
    struct _MatchMateInfoLess
	{
        template <typename TMInfo>
        inline bool 
        operator() (TMInfo const &a, TMInfo const &b) const 
		{
            return (a.contigId < b.contigId) || (a.contigId == b.contigId && a.beginPos < b.beginPos);
        }
    };

    template<typename TSpec, typename TConfig, typename TMatchMateInfos>
    inline void 
	_generatePairMatchIds (
		FragmentStore<TSpec, TConfig> & fragStore,
		TMatchMateInfos & matchMateInfos)
    {
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        
        typedef typename TFragmentStore::TReadStore			TReadStore;
        typedef typename TFragmentStore::TAlignedReadStore	TAlignedReadStore;
        typedef typename TFragmentStore::TContigPos			TContigPos;

        typedef typename Value<TReadStore>::Type						TRead;
        typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
        typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIter;    
		typedef typename Iterator<TMatchMateInfos, Standard>::Type		TMIter;    
        typedef typename Id<TFragmentStore>::Type						TId;
                
        TIter it = begin(fragStore.alignedReadStore, Standard());
		TIter itEnd = end(fragStore.alignedReadStore, Standard());
		TMIter mit = begin(matchMateInfos, Standard());
		TMIter mitEnd = end(matchMateInfos, Standard());

		if (it == itEnd || mit == mitEnd) return;

        // sort the aligned read store by: begin position, contig name
        sortAlignedReads(fragStore.alignedReadStore, SortBeginPos());
        sortAlignedReads(fragStore.alignedReadStore, SortContigId());
		std::sort(mit, mitEnd, _MatchMateInfoLess());

		TMIter mitNext = mit;
		while (mitNext != mitEnd && (*mit).beginPos == (*mitNext).beginPos)
			++mitNext;
        
//		if (qname == "EAS54_61:4:143:69:578")
//		readId =0;

		TContigPos pos = _min((*it).beginPos, (*it).endPos);
		typename Size<TReadStore>::Type readStoreLength = length(fragStore.readStore);

		while (mit != mitEnd)
		{
			int cmp = 0;
			if ((*it).contigId < (*mit).contigId) cmp = -1;
			else if ((*it).contigId > (*mit).contigId) cmp = 1;
			else if (pos < (*mit).beginPos) cmp = -1;
			else if (pos > (*mit).beginPos) cmp = 1;

			if (cmp == 1)
			{
				mit = mitNext;
				while (mitNext != mitEnd && (*mit).beginPos == (*mitNext).beginPos)
					++mitNext;
				continue;
			}
			if (cmp == 0)
			{
				for (TMIter m = mit; m != mitNext; ++m)
					if ((*it).readId != (*m).readId && (*it).readId < readStoreLength && (*m).readId < readStoreLength)		// correct position found
						if (fragStore.readStore[(*it).readId].matePairId == fragStore.readStore[(*m).readId].matePairId)	// correct mate found
							(*it).pairMatchId = (*m).pairMatchId;															// link mate
			}
			if (++it == itEnd) break;
			pos = _min((*it).beginPos, (*it).endPos);
		}
    }    
/*
//////////////////////////////////////////////////////////////////////////////
// comparePosLengthPair
    
    struct comparePosLengthPair{
        template<typename TPos>
        inline bool 
        operator() (TPos const & first, TPos const & second) const {
            return (first.i1 < second.i1);
        }
    };
    
//////////////////////////////////////////////////////////////////////////////
// _writeContigGapInStore
    
    template<typename TFragmentStore, typename TContigsGapPairs>
    inline void
    _writeContigsGapsInStore(TFragmentStore & fragStore, TContigsGapPairs & contigsGaps)
    {
        typedef typename Value<TContigsGapPairs>::Type							TContigGapPairs;
		typedef typename Value<TContigGapPairs>::Type							TGapPair;
        typedef typename Iterator<TContigGapPairs, Standard>::Type				TGapsIter;
		
		typedef typename TFragmentStore::TContigStore							TContigStore;
        typedef typename Value<TContigStore>::Type								TContig;
        typedef typename TContig::TContigSeq									TContigSeq;
        typedef typename TContig::TGapAnchors									TGapAnchors;

		typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >	TContigGaps;
		typedef typename Iterator<TContigGaps>::Type							TContigGapsIter;
        
        // for all contigs
        for (unsigned i = 0; i < length(contigsGaps); ++i)
		{
			TContigGaps contigGaps(fragStore.contigStore[i].seq, fragStore.contigStore[i].gaps);
            TContigGapPairs	& minGapLengths = contigsGaps[i];
            
            // Sort the contig gaps according to their start position
            std::sort(begin(minGapLengths, Standard()), end(minGapLengths, Standard()), comparePosLengthPair());
            
            // create iterator over contigGaps
			TGapsIter gapItEnd = end(minGapLengths, Standard());
			for (TGapsIter gapIt = begin(minGapLengths, Standard()); gapIt != gapItEnd; ++gapIt)
			{
				// get corresponding gapAnchor and iterator on it
				TContigGapsIter gapsIter(contigGaps, positionSeqToGap(contigGaps, (*gapIt).i1));
				++gapsIter;
				unsigned gapsCount = countGaps(gapsIter);
				std::cout << "c:"<<i<<"p:"<<(*gapIt).i1<<"  soll:"<<(*gapIt).i2<<"  ist:"<<gapsCount<<std::endl;
				if (gapsCount < (*gapIt).i2)
					insertGaps(gapsIter, (*gapIt).i2 - gapsCount);				
			}
        }
    }
*/
//////////////////////////////////////////////////////////////////////////////
// parsing functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// _parse_readCigar
    
    template <typename TFile, typename TCigarString, typename TChar>
    inline void
    _parse_readCigar(TFile & file, TCigarString & cigar, TChar & c)
    {
		typedef typename Value<TCigarString>::Type	TCigarElement;
		typedef typename TCigarElement::TOperation	TOperation;
		typedef typename TCigarElement::TCount		TCount;

		clear(cigar);
		
        // if the CIGAR is not set and '*'
        if (c == '*')
		{
            c = _streamGet(file);
            return;
        }
        
        while (!_streamEOF(file)) 
		{
            TCount count = _parse_readNumber(file, c);
			if (c >= 'a' && c <= 'z')
				c = c + 'A' - 'a';
            append(cigar, TCigarElement(c, count));
            
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n') break;
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_readSamIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parse_readSamIdentifier(TFile & file, TString & str, TChar& c)
    {
        append(str, c);
        while (!_streamEOF(file)) 
		{
            c = _streamGet(file);
            if (c== ' ' || c== '\t' || c == '\n') break;
            append(str, c);
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_is_dna
    
    template<typename TChar>
    inline bool
    _parse_is_dna(TChar const & c)
    {
        return c == (TChar)(Dna5)c;
    }
    
//////////////////////////////////////////////////////////////////////////////
//_parse_readDnaSeq
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parse_readDnaSeq(TFile & file, TString & str, TChar & c)
    {
		TChar first = c;
		if (!_streamEOF(file)) 
			c = _streamGet(file);

        if (!_parse_is_dna(first))
			return;
        append(str, first, Generous());
        
        for (; !_streamEOF(file) && _parse_is_dna(c); c = _streamGet(file))
            append(str, c, Generous());
    }
        
//////////////////////////////////////////////////////////////////////////////
// _parse_is_PhredQual
    
    template <typename TChar>
    inline bool
    _parse_is_PhredQual(TChar c)
    {
        return c >= '!' && c <= '~';
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_readSeqQual
//
    
    template<typename TFile, typename TQualString, typename TChar>
    inline void
    _parse_readSeqQual(TFile & file, TQualString & str, TChar & c)
    {
        typedef typename Size<TQualString>::Type				TSize;
        typedef typename Iterator<TQualString, Standard>::Type	TIter;
        
        TIter itBegin = begin(str, Standard());
        TSize rest = length(str);
        
        int q = 0;
        for (TIter it = itBegin; rest != 0 && _parse_is_PhredQual(c); --rest, ++it)
        {
            q = c - 33;
			if (!_streamEOF(file)) 
				c = _streamGet(file);
			else
				if (rest > 1)
					rest = 1;
			
			if (q == '*' - 33 && !_parse_is_PhredQual(c) && it == itBegin)
				return;
			
			// As the DNA5Q data structure can only store quality values till 60, all
			// qualities over this threshold are changed to 60
            if (q > 60) q = 60;
            assignQualityValue(*it, q);
        }
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_readCharsTillEndOfLine
//
// Reads all symbols till the next '\n' and writes them in the CharString str
// the c is the first character after the '\n'.
    
    template<typename TFile, typename TChar>
    inline void
    _parse_readCharsTillEndOfLine(TFile & file, String<char> & str, TChar& c)
    {
        // read all chars till '\n'
        while (c != '\n')
        {
            append(str, c, Generous());
			if (_streamEOF(file)) return;
	        c = _streamGet(file);
        }
        
        // read the first char after the '\n'
		if (!_streamEOF(file))
	        c = _streamGet(file);
    }
	
	
//////////////////////////////////////////////////////////////////////////////
// read functions for SAM
//////////////////////////////////////////////////////////////////////////////
    
//////////////////////////////////////////////////////////////////////////////
// read
    
    template<typename TFile, typename TSpec, typename TConfig>
    inline void 
    read (
		TFile & file,
		FragmentStore<TSpec, TConfig> & fragStore,
		SAM)
    {
        typedef Value<FILE>::Type TValue;
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
		typedef typename TFragmentStore::TContigPos TContigPos;
        typedef typename Id<TFragmentStore>::Type TId;
        
        // data structure to temporarily store the gaps that need to be inserted in the contig sequences
        typedef _MatchMateInfo<TContigPos, TId> TMatchMateInfo;
        typedef String<TMatchMateInfo> TMatchMateInfos;
        typedef StringSet<String<typename TFragmentStore::TContigGapAnchor> > TContigAnchorGaps;

        // data structure to temporarily store information about match mates
        TMatchMateInfos matchMateInfos;
		TContigAnchorGaps contigAnchorGaps;
        
        if (_streamEOF(file)) return;

		// get first character from the stream
        char c = _streamGet(file);
        
        // Read in header section
        _readHeader(file, fragStore, c, SAM());
        
        // Read in alignments section
        _readAlignments(file, fragStore, contigAnchorGaps, matchMateInfos, c, SAM());
        
        // set the match mate IDs using the information stored in matchMateInfos
        _generatePairMatchIds(fragStore, matchMateInfos);
        
        // insert gaps in the contigs using the information stored in contigsGaps
//            _writeContigsGapsInStore(fragStore, contigsGaps);

		convertPairWiseToGlobalAlignment(fragStore, contigAnchorGaps);
    }
    
//////////////////////////////////////////////////////////////////////////////
// _readHeader

    template<typename TFile, typename TSpec, typename TConfig, typename TChar>
    inline void 
    _readHeader (
		TFile & file,
		FragmentStore<TSpec, TConfig> & fragStore,
		TChar & c,
		SAM)
    {
        // skip header for now
        while(c == '@'){
            _parse_skipLine(file, c);
        }
        
    }
    
//////////////////////////////////////////////////////////////////////////////
// _readAlignments
//
// reads in alignement sections from a SAM file
    
    template<typename TFile, typename TSpec, typename TConfig, typename TContigAnchorGaps, typename TMatchMateInfos, typename TChar>
    inline void 
    _readAlignments (
		TFile & file,
        FragmentStore<TSpec, TConfig> & fragStore,
        TContigAnchorGaps & contigAnchorGaps,	
        TMatchMateInfos & matchMateInfos,
        TChar & c,
        SAM)
    {
        // create dummy entries in SAM specific aligned read quality store and aligned read tag store
        // is needed so the ID in the aligned store can be use to access the other stores
        // even if there exists previous entries without
		typedef FragmentStore<TSpec, TConfig> TFragmentStore;
		typedef typename TFragmentStore::TAlignQualityStore TAlignQualityStore;
		typedef typename TFragmentStore::TNameStore TNameStore;
		typedef typename Value<TAlignQualityStore>::Type TAlignQuality;
		
        int diff = length(fragStore.alignedReadStore) - length(fragStore.alignQualityStore);
        for(int i = 0; i < diff; ++i)
		{
			TAlignQuality q;
			q.score = supremumValue(q.score);
            append(fragStore.alignQualityStore, q, Generous());
        }
        diff = length(fragStore.alignedReadStore) - length(fragStore.alignedReadTagStore);
        for(int i = 0; i < diff; ++i){
            appendValue(fragStore.alignedReadTagStore, "", Generous());
        }
        
        // read in alignments
        int k = 0;
		Nothing contextSAM;
		
		refresh(fragStore.contigNameStoreCache);
		refresh(fragStore.readNameStoreCache);

        while(!_streamEOF(file)){
//            std::cout << k << std::endl;
            _readOneAlignment(file, fragStore, contigAnchorGaps, matchMateInfos, c, SAM(), contextSAM);
            ++k;
//			if (k%1000==0) std::cout <<'.'<<std::flush;
        }
    }
    
    
//////////////////////////////////////////////////////////////////////////////
// _readOneAlignment
//
// reads in one alignement section from a SAM file
    
    template<typename TFile, typename TSpec, typename TConfig, typename TContigAnchorGaps, typename TMatchMateInfos, typename TChar, typename TContextSAM>
    inline void 
    _readOneAlignment (
		TFile & file,
		FragmentStore<TSpec, TConfig> & fragStore,
		TContigAnchorGaps & contigAnchorGaps,
		TMatchMateInfos & matchMateInfos,
		TChar & c,
		SAM,
		TContextSAM & contextSAM)
    {
        SEQAN_CHECKPOINT
        // Basic types
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Id<TFragmentStore>::Type TId;
        typedef typename Size<TFragmentStore>::Type TSize;
        
        // All fragment store element types
        typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
        typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
        typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
        typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
        typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignQualityElement;
        typedef typename TAlignedElement::TGapAnchors TReadGapAnchors;
		
        // Type for sequence in readstore
        typedef typename TFragmentStore::TReadSeq TReadSeq2;
        
        // Type for gap anchor
//        typedef typename TFragmentStore::TReadGapAnchor TReadGapAnchor2;
        typedef typename TFragmentStore::TContigPos TContigPos;
		typedef Gaps<TReadSeq2, AnchorGaps<TReadGapAnchors> >	TReadGaps;
		typedef Gaps<Nothing, AnchorGaps<typename Value<TContigAnchorGaps>::Type> >	TContigGapsPW;
        
        // Type to temporarily store information about match mates
        typedef typename Value<TMatchMateInfos>::Type TMatchMateInfo;
        
        // read fields of alignments line
        
        // read the query name
        String<char> qname;
        _parse_readSamIdentifier(file, qname, c);
        _parse_skipWhitespace(file, c);

		// read the flag
        int flag;
        flag = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
		bool reverse = (flag & (1 << 4)) == (1 << 4);

		// read reference name
        String<char> rname;
        _parse_readSamIdentifier(file, rname, c);
        _parse_skipWhitespace(file, c);

		// read begin position
        TContigPos beginPos;
        beginPos = _parse_readNumber(file, c);
        --beginPos; // SAM stores positions starting at 1 the fragment store starting at 0
        _parse_skipWhitespace(file, c);

        // read map quality
        TAlignQualityElement mapQ;
        mapQ.score = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);

		// read CIGAR
        String<CigarElement<> > cigar;
        _parse_readCigar(file, cigar, c);
        _parse_skipWhitespace(file, c);
        
        // calculate the end position
        TContigPos endPos;
        _getClippedLength(cigar, endPos);
        endPos = beginPos + endPos;
        // if the read is on the antisense strand switch begin and end position
        if (reverse)
		{
            TContigPos temp = beginPos;
            beginPos = endPos;
            endPos = temp;
        }

		// generate gap anchor string for the read
        TReadGapAnchors readGapAnchors;
        
        // read mate reference name
        String<char> mrnm;
        _parse_readSamIdentifier(file, mrnm, c);
        _parse_skipWhitespace(file, c);

		// read mate position
        TContigPos mPos;
        mPos = _parse_readNumber(file, c);
        --mPos; // SAM stores positions starting at 1 the fragment store starting at 0
        _parse_skipWhitespace(file, c);

		// read iSize
        _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);

		// read in sequence
        TReadSeq2 readSeq;
        _parse_readDnaSeq(file, readSeq, c);
		if (reverse)
			reverseComplementInPlace(readSeq);
        _parse_skipWhitespace(file, c);

		// and associated qualities
        _parse_readSeqQual(file, readSeq, c);

		// insert alignment gaps
		TReadGaps readGaps(readSeq, readGapAnchors);
        cigarToGapAnchorRead(cigar, readGaps);
        
        // read in SAM tags
        String<char> tags;
        _parse_skipSpace(file, c);
        _parse_readCharsTillEndOfLine(file, tags, c);
        
        
        // check if read sequence is already in the store.
        // if so get the ID, otherwise create new entries in the
        // read, read name and mate pair store
        
        TId readId = 0;
        _appendRead(fragStore, readId, qname, readSeq, flag, contextSAM);
        
        // check if the contig is already in the store
        // get its ID or create a new one otherwise
        TId contigId = 0;
        _appendContig(fragStore, contigId, rname);

		if (empty(cigar)) return;
		
        // insert gaps in the contigGaps
//        cigarToContigGaps(cigar, beginPos, reverse, value(contigsGaps, contigId));

        // create a new entry in the aligned read store
        TId pairMatchId = appendAlignment(fragStore, readId, contigId, beginPos, endPos, readGapAnchors);
		resize(contigAnchorGaps, length(fragStore.alignedReadStore), Generous());
		TContigGapsPW contigGaps(back(contigAnchorGaps));
        cigarToGapAnchorContig(cigar, contigGaps);
		
		// create entries in SAM specific stores
        appendValue(fragStore.alignQualityStore, mapQ, Generous());
        appendValue(fragStore.alignedReadTagStore, tags, Generous());
        
        // store additional data about match mate temporarily
        // used in the end of the read function to generate match mate IDs
		TId mcontigId = contigId;
        if (mrnm != "*")
		{
			if (mrnm != "=")
				_appendContig(fragStore, mcontigId, mrnm);

			if (flag & 0x40)	// store mate info only for the first read in the pair
			{
				TMatchMateInfo matchMateInfo = {readId, mcontigId, pairMatchId, mPos};
				append(matchMateInfos, matchMateInfo);
				back(fragStore.alignedReadStore).pairMatchId = pairMatchId;
			}
		}
    }
    
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
