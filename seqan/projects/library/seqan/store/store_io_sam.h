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
struct TagSAM_;
typedef Tag<TagSAM_> const SAM;

    
//////////////////////////////////////////////////////////////////////////////
// Cigar struct
//////////////////////////////////////////////////////////////////////////////
    
    template<typename TChar = char, typename TNum = int>
    struct Cigar
    {
        String<TChar> operationType;
        String<TNum> operationCount;
        
        Cigar() : operationType(0), operationCount(0) {}
    };

//////////////////////////////////////////////////////////////////////////////
// append
    
    template<typename TChar, typename TNum>
    inline void append(Cigar<TChar, TNum> & cigar, TChar const c, TNum const i)
    {
        int l = length(cigar.operationType);
        
        resize(cigar.operationType, l + 1);
        assignValue(cigar.operationType, l, c);
        
        resize(cigar.operationCount, l + 1);
        assignValue(cigar.operationCount, l - 1, i);
        
    }

//////////////////////////////////////////////////////////////////////////////
// _getClippedLength
    
    template<typename TChar, typename TNum>
    inline void _getClippedLength(Cigar<TChar, TNum> const & cigar, TNum & sum)
    {
        typedef typename Iterator< String<TNum> >::Type TNumIter;
        typedef typename Iterator< String<TChar> >::Type TCharIter;
        
        TNumIter it_c = begin(cigar.operationCount);
        TCharIter it_t = begin(cigar.operationType);
        
        sum = 0;
        
        while(it_c != end(cigar.operationCount)){
            if(value(it_t) != 'S' && value(it_t) != 's' && value(it_t) != 'H' && value(it_t) != 'h')
                sum += value(it_c);
            
            ++it_c;
            ++it_t;
        }
    }

//////////////////////////////////////////////////////////////////////////////
// cigarToGapAnchorRead
    
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
        
        
        // check if there is a (soft) clipped sequence in the begining
        if(value(it_t) == 'S' || value(it_t) == 's'){
            
            // skip first 'count' characters
            pos_seq += value(it_c);
            
            // insert in GapAnchor
            TGapAnchor anchor = TGapAnchor(pos_seq, pos_gapped);
            append(gaps, anchor);
        }
        
        while(it_c != end(cigar.operationCount)){
            
            // If the operation type is match/mismatch or insertion there is no gap in the read sequence.
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
            
            // If the operation type is deletion or padding there is a gap in the read sequence.
            // TODO: what to do with 'skipped in reference sequence' (N)?
            if(value(it_t) == 'D' || value(it_t) == 'P' || value(it_t) == 'd' || value(it_t) == 'p'){

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
    
//////////////////////////////////////////////////////////////////////////////
// some helping functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// _print_gapAnchor
    
    template<typename TChar, typename TGapAnchor>
    inline void _print_gapAnchor(String<TChar> seq, String<TGapAnchor> gaps)
    {
        typedef typename Iterator< String<TChar> >::Type TCharIter;
        typedef typename Iterator< String<TGapAnchor> >::Type TGapIter;
        
        // create iterators and set them on the begining of the strings
        TCharIter it_c = begin(seq);
        TGapIter  it_g = begin(gaps);
        
        //
        int gapped = value(it_g).gapPos, last_gapped = gapped;
        int ungapped = value(it_g).seqPos, last_ungapped = ungapped;
        
        for(int i = 0; i < ungapped; ++i)
            ++it_c;
        
        while(it_g != end(gaps)){
            gapped = value(it_g).gapPos;
            ungapped = value(it_g).seqPos;
            //std::cout << "(" << ungapped << ", " << gapped << ")";
            
            int count_chars = ungapped - last_ungapped;
            int count_gaps = (gapped - last_gapped) - count_chars;
            
            for(int i = 0; i < count_chars && it_c != end(seq); ++i, ++it_c){
                _streamPut(std::cout, value(it_c));
            }
            
            for(int i = 0; i < count_gaps; ++i){
                _streamPut(std::cout, '-');
            }
            
            // iterate
            last_gapped = gapped;
            last_ungapped = ungapped;
            ++it_g;
        }
        
        while(it_c != end(seq)){
            _streamPut(std::cout, value(it_c));
            ++it_c;
        }
    }

//////////////////////////////////////////////////////////////////////////////
// _getMate
//
// checks if a read name is already in read name store and if corresponding 
// begin position in the aligned read strore fits the mate position
// returns the ID or -1 otherwise
    
    template<typename TSpec, typename TConfig, typename TPos, typename TID>
    inline void 
    _getMate(FragmentStore<TSpec, TConfig> & fragStore, String<char> & readName, TPos & mPos, TID & mateID)
    {
        typedef Iterator<StringSet<CharString> >::Type TNameStoreIter;
        
        // set mate ID = -1. Will be replaced if mate is found
        mateID = -1;
        
        // Iterator over read names
        TNameStoreIter it_rm = begin(fragStore.readNameStore);
        
        while (it_rm != end(fragStore.readNameStore)){
            // if the read name was found
            if(readName == value(it_rm)){
                // check in the aligned read store if the begin position is correct
                int suggestedID = position(it_rm);
                
                if(value(fragStore.alignedReadStore, suggestedID).beginPos == mPos){
                    mateID = suggestedID;
                    return;
                }                
            }
            
            ++it_rm;
        }
        
    }
    
//////////////////////////////////////////////////////////////////////////////
// _getID
    
    template<typename TType, typename TID>
    inline bool 
    _getID(StringSet<TType> & store, TType & elem, TID & elem_id)
    {
        typedef Iterator<StringSet<CharString> >::Type TNameStoreIter;
        
        // set mate ID = -1. Will be replaced if mate is found
        elem_id = -1;
        
        // Iterator over read names
        TNameStoreIter iter = begin(store);
        
        while (iter != end(store)){
            // if the element was found
            if(elem == value(iter)){
                // set the ID
                elem_id = position(iter);
                // and break the loop
                return true;               
            }
            
            ++iter;
        }
        return false;
    }
    
//////////////////////////////////////////////////////////////////////////////
// parsing functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// _parse_readCigar
    
    template<typename TFile, typename TChar, typename TNum>
    inline void
    _parse_readCigar(TFile & file, Cigar<TChar, TNum> & cigar, TChar& c)
    {
        TChar type;
        TNum count;
        
        while (!_streamEOF(file)) {
            count = _parse_readNumber(file, c);
            type = c;
            append(cigar, type, count);
            
            c = _streamGet(file);
            if (c== ' ' || c== '\t' || c == '\n') break;
        }
    }

    
//////////////////////////////////////////////////////////////////////////////
// _parse_readSamIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parse_readSamIdentifier(TFile & file, TString & str, TChar& c)
    {
        append(str, c);
        while (!_streamEOF(file)) {
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
        return ((c == 'a') || (c == 'c') || (c == 'g') || (c == 't') || (c == 'n') || (c == 'A') || (c == 'C') || (c == 'G') || (c == 'T') || (c == 'N'));
    }
    
//////////////////////////////////////////////////////////////////////////////
//_parse_readDnaSeq
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parse_readDnaSeq(TFile & file, TString & str, TChar& c)
    {
        if (!_parse_is_dna(c)) return;
        append(str, c);
        
        while (!_streamEOF(file)) {
            c = _streamGet(file);
            if (!_parse_is_dna(c)) break;
            append(str, c);
        }
        
    }
        
//////////////////////////////////////////////////////////////////////////////
// _parse_is_PhredQual
    
    template<typename TChar>
    inline bool
    _parse_is_PhredQual(TChar const & c)
    {
        return ( ((unsigned) c > 32) && ((unsigned) c < 127) );
    }
    
//////////////////////////////////////////////////////////////////////////////
// _parse_readSeqQual
//
// As the DNA5Q data structure can only store quality values till 40, all
// qualities over this threshold are changed to 40
    
    template<typename TFile, typename TChar>
    inline void
    _parse_readSeqQual(TFile & file, String<Dna5Q> & str, TChar& c)
    {
        typedef typename Iterator< String<Dna5Q> >::Type TIter;
        
        TIter it = begin(str);
        int q = 0;
        
        do
        {
            if (!_parse_is_PhredQual(c)) break;
            
            q = c - 32;
            if (q > 39) q = 40;
            assignQualityValue(value(it), q);
            
            c = _streamGet(file);
        }
        while (!_streamEOF(file) && it != end(str));
    }
    
//////////////////////////////////////////////////////////////////////////////
// read functions for SAM
//////////////////////////////////////////////////////////////////////////////
    template<typename TFile, typename TSpec, typename TConfig, typename TChar>
    inline void 
    _readAlignments(TFile& file,
        FragmentStore<TSpec, TConfig>& fragStore,
        TChar & c,
        SAM)
    {
        while(!_streamEOF(file)){
            _readOneAlignment(file, fragStore, c, SAM());
        }
    }
    
    
//////////////////////////////////////////////////////////////////////////////    
    template<typename TFile, typename TSpec, typename TConfig, typename TChar>
    inline void 
    _readOneAlignment(TFile& file,
        FragmentStore<TSpec, TConfig>& fragStore,
        TChar & c,
        SAM)
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
        
        // Type for sequence in readstore
        typedef typename TFragmentStore::TReadSeq TReadSeq2;
        
        // Type for gap anchor
        typedef typename TFragmentStore::TReadGapAnchor TReadGapAnchor2;
        
        // TODO: typedef for TPos in aligned, replace int of begin and end position with that
        
        // read fields of alignments line
        
        // read teh query name
        String<char> qname;
        _parse_readSamIdentifier(file, qname, c);
        _parse_skipWhitespace(file, c);
        
        // read the flag
        int flag;
        flag = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
        
        // read reference name
        String<char> rname;
        _parse_readSamIdentifier(file, rname, c);
        _parse_skipWhitespace(file, c);
        
        // read begin position
        int beginPos;
        beginPos = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
        
        // read map quality
        int mapQ;
        mapQ = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
        
        // read CIGAR
        Cigar<> cigar = Cigar<>();
        _parse_readCigar(file, cigar, c);
        _parse_skipWhitespace(file, c);
        
        // generate gap anchor string for the read
        String<TReadGapAnchor2 > readGaps;
        cigarToGapAnchorRead(cigar, readGaps);
        
        // read mate reference name
        String<char> mrnm;
        _parse_readSamIdentifier(file, mrnm, c);
        _parse_skipWhitespace(file, c);
        
        // read mate position
        int mPos;
        mPos = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
        
        // read iSizs
        int iSize;
        iSize = _parse_readNumber(file, c);
        _parse_skipWhitespace(file, c);
        
        // read in sequence
        TReadSeq2 readSeq;
        _parse_readDnaSeq(file, readSeq, c);
        _parse_skipWhitespace(file, c);
        
        // and associated qualities
        _parse_readSeqQual(file, readSeq, c);
        
        // skip tags and go to the end of the line
        _parse_skipLine(file, c);
        
        
        std::cout << "qname: \t" << qname << std::endl;
        std::cout << "flag: \t" << flag << std::endl;
        std::cout << "rname: \t" << rname << std::endl;
        std::cout << "pos: \t" << beginPos << std::endl;
        std::cout << "mapQ: \t" << mapQ << std::endl;

        std::cout << "mrnm: \t" << mrnm << std::endl;
        std::cout << "flag: \t" << flag << std::endl;
        std::cout << "mPos: \t" << mPos << std::endl;
        std::cout << "iSize: \t" << iSize << std::endl;
        std::cout << "seq: \t" << readSeq << std::endl;
        
        if(qname == "EAS219_FC30151:5:29:817:854"){
            int wasd = 0;
        }
        
        // check if read sequence is already in the store.
        // if so get the ID, otherwise create new entries in the
        // read, read name and mate pair store
        
        TId read_ID = 0;
        if(_getID(fragStore.readNameStore, qname, read_ID)){
            // if the read is paired
            if((flag & 1) == 1){
                // check the mate pair store if it is the same mate of the pair
                // assuming that only one flag 0x040 or 0x0080 is 1
                int inPair = ((flag & (1 << 7)) >> 7);
                
                TId matePair_ID = value(fragStore.readStore, read_ID).matePairId;

                read_ID = value(fragStore.matePairStore, matePair_ID).readId[inPair];
                
                if(read_ID == TMatePairElement::INVALID_ID){
                    // create new entry in read and read name store
                    read_ID = length(fragStore.readStore);
                    // set sequence and mate pair ID in new read store element
                    TReadStoreElement readElem = TReadStoreElement();
                    readElem.seq = readSeq;
                    readElem.matePairId = matePair_ID;
                    append(fragStore.readStore, readElem);
                    
                    // add the identifyer to the read name store
                    appendValue(fragStore.readNameStore, qname);
                    
                    // set the ID in the mate pair store
                    value(fragStore.matePairStore, matePair_ID).readId[inPair] = read_ID;
                }
            }
        } else { // if the read name is not in the store
            // create new entry in read and read name store
            read_ID = length(fragStore.readStore);
            TReadStoreElement readElem = TReadStoreElement();
            readElem.seq = readSeq;
            
            // if the read is paired
            if((flag & 1) == 1){
                TMatePairElement mateElem = TMatePairElement();
                // set the first or second read ID in the mate pair element
                mateElem.readId[((flag & (1 << 7)) >> 7)] = read_ID;
                // get a new mate pair ID and add the new mate pair element
                TId matePair_ID = length(fragStore.matePairStore);
                append(fragStore.matePairStore, mateElem);
                // set the new mate pair ID in the read element
                readElem.matePairId = matePair_ID;
            }
            
            append(fragStore.readStore, readElem);
            appendValue(fragStore.readNameStore, CharString());
        }
        
        // check if the contig is already in the store
        // get its ID or create a new one otherwise
        // TODO: contig
        
        // create a new entry in the aligned read store
        
        TAlignedElement alignedElem = TAlignedElement();
        alignedElem.readId = read_ID;
        alignedElem.beginPos = beginPos;
        
        // calculate the end position
        int div;
        _getClippedLength(cigar, div);
        // if the read is on the antisense strand the end position lies before the begin position
        if((flag & (1 << 4)) == (1 << 4))
            div *= -1;
        
        alignedElem.endPos = beginPos + div;
        alignedElem.gaps = readGaps;

        std::cout << "======================================" << std::endl;
        
    }
    
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
