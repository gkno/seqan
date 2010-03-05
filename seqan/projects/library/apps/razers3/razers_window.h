/*==========================================================================
 RazerS - Fast Read Mapping with Controlled Loss Rate
 http://www.seqan.de/projects/razers.html
 
 ============================================================================
 Copyright (C) 2010
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your options) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#ifndef SEQAN_HEADER_RAZERS_WINDOW_H
#define SEQAN_HEADER_RAZERS_WINDOW_H

#include <iostream>

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // Find read matches in a single genome sequence
	// 
	// Creates finder on the contig given by the ID and a verifier.
	// Searches through the contig using the findWindowNext() function.
	// The results are dumped in the (aligned) store.
    template <
	typename TFragmentStore, 
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TPreprocessing,
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode >
    void _mapSingleReadsToContigWindow(
         TFragmentStore							& store,
         unsigned								  contigId,				// ... and its sequence number
         Pattern<TReadIndex, Swift<TSwiftSpec> >& swiftPattern,
         TPreprocessing							& preprocessing,
         TCounts								& cnts,
         char									  orientation,				// q-gram index of reads
         TRazerSOptions							& options,
         TRazerSMode                      const & mode)
    {
        // FILTRATION
        typedef typename TFragmentStore::TContigSeq				TContigSeq;
        typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
        typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;
        
        // VERIFICATION
        typedef MatchVerifier <
            TFragmentStore, 
            TRazerSOptions, 
            TRazerSMode,
            TPreprocessing, 
            TSwiftPattern,
            TCounts >											TVerifier;
        typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
        
        // HITS
        typedef typename TSwiftFinder::THitString               THitString;
		typedef typename Value<THitString>::Type                TSwiftHit;
        typedef typename Size<THitString>::Type               THitStringSize;
        
		// output what is done if verbous
        if (options._debugLevel >= 1)
        {
            ::std::cerr << ::std::endl << "Process genome seq #" << contigId;
            if (orientation == 'F') ::std::cerr << "[fwd]";
            else                    ::std::cerr << "[rev]";
        }
		// lock contig
        lockContig(store, contigId);
        TContigSeq &contigSeq = store.contigStore[contigId].seq;
        if (orientation == 'R')	reverseComplementInPlace(contigSeq);
        
		// Create finder and verifier
        TSwiftFinder	swiftFinder(contigSeq, options.repeatLength, 1);
        TVerifier		verifier(store, options, preprocessing, swiftPattern, cnts);
        
        // initialize verifier
        verifier.onReverseComplement = (orientation == 'R');
        verifier.genomeLength = length(contigSeq);
        verifier.m.contigId = contigId;
        
		// if the pattern can be initialized and there is a non-repeat region in the contig that fits a qgram.
        if(windowFindBegin(swiftFinder, swiftPattern, options.errorRate, options._debugLevel)){
            
			// while there is more contig sequence to search through
            while(windowFindNext(swiftFinder, swiftPattern, 1000, options._debugLevel)){
                // get the found hits from the finder
                THitString hits = getSwiftHits(swiftFinder);
                // verifiy them
                for(THitStringSize h = 0; h < length(hits); ++h){
                    verifier.m.readId = hits[h].ndlSeqNo;         //array oder jedesmal berechnen
                    matchVerify(verifier, getSwiftRange(hits[h], contigSeq), hits[h].ndlSeqNo, host(host(swiftPattern)), mode);
                    ++options.countFiltration;
                }
            }
            
            // clear finders
            windowFindEnd(swiftFinder, swiftPattern);
            
        }
        

        if (!unlockAndFreeContig(store, contigId))							// if the contig is still used
            if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation
    }
    
}

#endif