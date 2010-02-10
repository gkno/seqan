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

#ifndef SEQAN_HEADER_RAZERS_PARALLEL_2_H
#define SEQAN_HEADER_RAZERS_PARALLEL_2_H

//#include "tbb/pipeline.h"
//#include "tbb/spin_mutex.h"
//#include "tbb/task_scheduler_init.h"
#include <omp.h>

namespace SEQAN_NAMESPACE_MAIN
{
    
    //////////////////////////////////////////////////////////////////////////////
    // Find read matches in a single genome sequence
    template <
        typename TFragmentStore, 
        typename TReadIndex, 
        typename TSwiftSpec, 
        typename TPreprocessing,
        typename TCounts,
        typename TRazerSOptions,
        typename TRazerSMode >
    void _mapSingleReadsToContig(
             TFragmentStore                                     & store,
             unsigned                                             contigId,				// ... and its sequence number
             String<Pattern<TReadIndex, Swift<TSwiftSpec> > >   & swiftPatterns,
             TPreprocessing                                     & preprocessingBlocks,
             TCounts                                            & cnts,
             char                                                 orientation,			// q-gram index of reads
             TRazerSOptions                                     & options,
             TRazerSMode const                                  & mode)
    {
        // FILTRATION
        typedef typename TFragmentStore::TContigSeq				TContigSeq;
        typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
        typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;
        typedef typename Size<TSwiftFinder>::Type               TSwiftFinderSize;

        // HITS
        typedef typename TSwiftFinder::THitString               THitString;
		typedef typename Value<THitString>::Type                TSwiftHit;
        typedef typename Size<THitString>::Type               THitStringSize;
        
        // VERIFICATION
        typedef typename Value<TPreprocessing>::Type                                                                        TPreprocessingBlock;
        typedef MatchVerifier <TFragmentStore, TRazerSOptions, TRazerSMode, TPreprocessingBlock, TSwiftPattern, TCounts >   TVerifier;
        typedef typename Fibre<TReadIndex, Fibre_Text>::Type                                                                TReadSet;
        
        if (options._debugLevel >= 1)
        {
            ::std::cerr << ::std::endl << "Process genome seq #" << contigId;
            if (orientation == 'F') ::std::cerr << "[fwd]";
            else                    ::std::cerr << "[rev]";
        }
        lockContig(store, contigId);
        TContigSeq &contigSeq = store.contigStore[contigId].seq;
        if (orientation == 'R')	reverseComplementInPlace(contigSeq);
        
        // Finder and verifier strings of the same size as there are swift patterns
        // One Swift finder, Swift pattern and verifier work together
        String<TSwiftFinder> swiftFinders;
        resize(swiftFinders, length(swiftPatterns), Exact());
        TSwiftFinder swiftFinder(contigSeq, options.repeatLength, 1);
        
        String<TVerifier> verifier;
        resize(verifier, length(swiftPatterns), Exact());
        
        for(unsigned i = 0; i < length(swiftPatterns); ++i){
            // initialize finders
            swiftFinders[i] = swiftFinder;
            
            // initialize verifier
            TVerifier oneVerifier(store, options, preprocessingBlocks[i], swiftPatterns[i], cnts);
            oneVerifier.onReverseComplement = (orientation == 'R');
            oneVerifier.genomeLength = length(contigSeq);
            oneVerifier.m.contigId = contigId;
            
            verifier[i] = oneVerifier;
        }
        
        // set up finder
        for (unsigned i = 0; i < length(swiftFinders); ++i)
            windowFindBegin(swiftFinders[i], swiftPatterns[i], options.errorRate, options._debugLevel);
        
        // FIXME: razerSoption
        unsigned windowSize = 5000;
        unsigned offSet = 0;
        
        // go over contig sequence
        while(offSet < length(contigSeq))
        {
            // each filter runs over the window separately
            // FIXME: Parallelize
# pragma omp parallel for shared(i)
            for (TSwiftFinderSize i = 0; i < length(swiftFinders); ++i){
                
                // filter window and save hits
                windowFindNext(swiftFinders[i], swiftPatterns[i], windowSize, options._debugLevel);

                // verify hits
                THitString hits = getSwiftHits(swiftFinders[i]);
                for(THitStringSize h = 0; h < length(hits); ++h)
                    verifier[i].m.readId = hits[h].ndlSeqNo;
                    // FIXME: arguments are not correct
                    matchVerify(verifier[i], range(swiftFinders[i], contigSeq), verifier[i].m.readId, host(host(swiftPatterns[i])), mode);
                
            }
            
            // move window
            offSet += windowSize;
        }
        
        // clear finders
        for (unsigned i = 0; i < length(swiftFinders); ++i)
            windowFindEnd(swiftFinders[i], swiftPatterns[i]);
        
        
        /*while (find(swiftFinder, swiftPattern, options.errorRate, options._debugLevel)) 
        {
            verifier.m.readId = (*swiftFinder.curHit).ndlSeqNo;
            if (!options.spec.DONT_VERIFY)
                matchVerify(verifier, range(swiftFinder, contigSeq), verifier.m.readId, readSet, mode);
            ++options.countFiltration;
        }*/
        
        if (!unlockAndFreeContig(store, contigId))							// if the contig is still used
            if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation
    }
    
    
    //////////////////////////////////////////////////////////////////////////////
    // Find read matches in many genome sequences
    template <
        typename TFSSpec, 
        typename TFSConfig, 
        typename TCounts,
        typename TSpec, 
        typename TAlignMode,
        typename TGapMode,
        typename TScoreMode,
        typename TReadIndexString>
    int _mapSingleReadsParallel(
                        FragmentStore<TFSSpec, TFSConfig>					& store,
                        TCounts												& cnts,
                        RazerSOptions<TSpec>								& options,
                        RazerSMode<TAlignMode, TGapMode, TScoreMode>  const & mode,
                        TReadIndexString									& readIndices)
    {
        typedef FragmentStore<TFSSpec, TFSConfig>                   TFragmentStore;
        typedef typename IF<TYPECMP<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
        typedef typename Value<TReadIndexString>::Type              TReadIndex;
        
        // filter
        typedef Pattern<TReadIndex, Swift<TSwiftSpec> >             TSwiftPattern;
        
        // verifier
        typedef Pattern<TRead, MyersUkkonen>                        TMyersPattern;
        typedef typename Infix<String<TMyersPattern> >::Type        TVerifierBlock;
        typedef typename Position<TReadIndexString>::Type           TPos;

        // configure Swift patterns
        String<TSwiftPattern> swiftPatterns;
        resize(swiftPatterns, length(readIndices), Exact());
        for (TPos i = 0; i < length(readIndices); ++i) {
            TSwiftPattern swiftPattern(readIndices[i]);
            swiftPattern.params.minThreshold = options.threshold;
            swiftPattern.params.tabooLength = options.tabooLength;
            swiftPatterns[i] = swiftPattern;
        }
        
        // init edit distance verifiers
        unsigned readCount = length(store.readSeqStore);//countSequences(readIndices); // FIXME: function for String of StringSets has the return the total number of reads
        String<TMyersPattern> forwardPatterns;
        resize(forwardPatterns, length(store.readSeqStore), Exact()); // FIXME:
        
        options.compMask[4] = (options.matchN)? 15: 0;
        if (options.gapMode == RAZERS_GAPPED)
        {
            resize(forwardPatterns, readCount, Exact());
            for(unsigned i = 0; i < readCount; ++i)
            {
                setHost(forwardPatterns[i], store.readSeqStore[i]); // FIXME: store.readSeqStore insead of indexText(...)
                _patternMatchNOfPattern(forwardPatterns[i], options.matchN);
                _patternMatchNOfFinder(forwardPatterns[i], options.matchN);
            }
        }
        
        // divide the verification Patterns in the same blocks as the filters
        String<TVerifierBlock> forwardPatternsBlocks;
        resize(forwardPatternsBlocks, length(readIndices), Exact());
        TPos offSet = 0;
        for (TPos i = 0; i < length(readIndices); ++i){
            // block of the same size as the corresponding index (number of reads in this index)
            TVerifierBlock block = infix(forwardPatterns, offSet, length(host(readIndices[i])));
            offSet += length(host(readIndices[i]));
            
            forwardPatternsBlocks[i] = block;
        }
        
        if (options.maqMapping)
        {
            resize(cnts, 2);
            for (unsigned i = 0; i < length(cnts); ++i)
                fill(cnts[i], readCount, 31); //initialize with maxeditDist, 11:5 for count:dist
        }
        
        // clear stats
        options.countFiltration = 0;
        options.countVerification = 0;
        options.timeMapReads = 0;
        options.timeDumpResults = 0;
        SEQAN_PROTIMESTART(find_time);
        
        // iterate over genome sequences
        for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId)
        {
            // lock to prevent releasing and loading the same contig twice
            // (once per _mapSingleReadsToContig call)
            lockContig(store, contigId);
            if (options.forward)
                _mapSingleReadsToContig(store, contigId, swiftPatterns, forwardPatternsBlocks, cnts, 'F', options, mode);
            
            if (options.reverse)
                _mapSingleReadsToContig(store, contigId, swiftPatterns, forwardPatternsBlocks, cnts, 'R', options, mode);
            unlockAndFreeContig(store, contigId);
        }
        
        options.timeMapReads = SEQAN_PROTIMEDIFF(find_time);
        if (options._debugLevel >= 1)
            ::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
        if (options._debugLevel >= 2) {
            ::std::cerr << ::std::endl;
            ::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
            ::std::cerr << "Filtration counter:      " << options.countFiltration << ::std::endl;
            ::std::cerr << "Successful verfications: " << options.countVerification << ::std::endl;
        }
        return 0;
    }
    
    
    template <
        typename TFSSpec, 
        typename TFSConfig, 
        typename TCounts,
        typename TSpec, 
        typename TShape,
        typename TAlignMode,
        typename TGapMode,
        typename TScoreMode >
    int _mapSingleReadsParallel(
            FragmentStore<TFSSpec, TFSConfig>					& store,
            TCounts												& cnts,
            RazerSOptions<TSpec>								& options,
            TShape const										& shape,
            RazerSMode<TAlignMode, TGapMode, TScoreMode>  const & mode)
    {
        typedef FragmentStore<TFSSpec, TFSConfig>                       TFragmentStore;
        typedef typename TFragmentStore::TReadSeqStore                  TReadSeqStore;
        typedef typename Value<TReadSeqStore>::Type                     TRead;
        typedef StringSet<TRead>                                        TReadSet;
        typedef Index<TReadSet, Index_QGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
        typedef typename Size<TReadSeqStore>::Type                      TSize;
        
        // number of cores and how many blocks there should be
        // a block is a subset of reads that is filtered and verified in a row
        // depending on how long it takes to process the individual blocks a single
        // thread might work through more than otheres
        unsigned cores = 2; //omp_get_num_procs();
        unsigned noOfBlocks = cores * 5; // TODO: razerSoption
        
        // if there are not enough reads that the parallel version use the normal one
        if(length(store.readSeqStore) < 1000){ // TODO: razerSoption
            return _mapSingleReads(store, cnts, options, shape, mode);
        } 
        else {
            // number of reads per thread
            unsigned perBlock = length(store.readSeqStore) / noOfBlocks;
            unsigned perBlockRemainder = length(store.readSeqStore) / noOfBlocks;
            
            String<TIndex> indices;
            resize(indices, noOfBlocks, Exact());
            
            int readID = 0;
            // create swift indices that can work in parallel
            for (unsigned b = 0; b < noOfBlocks; ++b) {
                TReadSet readSet;
                resize(readSet, perBlock, Exact());
                
                // the last one gets some extra
                if((b == noOfBlocks - 1) && (perBlockRemainder != 0)){
                    resize(readSet, perBlock + perBlockRemainder, Exact());
                }
                
                // get a subset of reads from the store
                for (unsigned i = 0; i < length(readSet); ++i) {
                    assign(readSet[i], store.readSeqStore[readID]);
                    ++readID;
                }
                
                // configure q-gram index
                TIndex swiftIndex(readSet, shape);
                cargo(swiftIndex).abundanceCut = options.abundanceCut;
                cargo(swiftIndex)._debugLevel = options._debugLevel;
                
                indices[b] = swiftIndex;
            }
            
            return _mapSingleReadsParallel(store, cnts, options, mode, indices);
        }
        

        
    }
    
}

#endif
