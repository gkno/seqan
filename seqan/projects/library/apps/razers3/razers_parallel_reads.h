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

#include <iostream>
#include <omp.h>

//#define RAZERS_PARALLEL_READS_WHOLE_GENOME
#ifdef RAZERS_PARALLEL_READS_WHOLE_GENOME
#include "razers_window.h"
#endif

namespace SEQAN_NAMESPACE_MAIN
{
    
    // TODO: doc
    template <typename TSwiftPatterns>
    struct ParallelSwiftPatternHandler
    {
        SEQAN_CHECKPOINT
        
        TSwiftPatterns &swiftPatterns;
        ParallelSwiftPatternHandler(TSwiftPatterns &_swiftPatterns):
            swiftPatterns(_swiftPatterns) {}
    };
    
    // TODO: doc
    template<
        typename TReadSet,
        typename TShape>
    void intiIndex(Index<TReadSet, Index_QGram<TShape, OpenAddressing> > & index, TReadSet & _text, TShape const & _shape)
    {
        SEQAN_CHECKPOINT
        
        value(index.text) = _text;
        index.shape = _shape;
    }

    // TODO: doc
    template < typename TSwiftPatterns, typename TReadNo, typename TMaxErrors >
    inline void 
    setMaxErrors(ParallelSwiftPatternHandler<TSwiftPatterns> &swift, TReadNo readNo, TMaxErrors maxErrors)
    {
        SEQAN_CHECKPOINT
        
        int blockSize = length(host(host(swift.swiftPatterns[0])));
        int indexNo = readNo / blockSize;
        int localReadNo = readNo % blockSize;
        
        int minT = _qgramLemma(swift.swiftPatterns[indexNo], localReadNo, maxErrors);
        if (minT > 1)
        {
            //		::std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
            if (maxErrors < 0) minT = SupremumValue<int>::VALUE;
            setMinThreshold(swift.swiftPatterns[indexNo], localReadNo, (unsigned)minT);
        }
    }
	
	// iterates over the aligned read store and look for a certain begin position. should set break point in if clause
	template <typename TFragmentStore>
	inline void containsRead(TFragmentStore const & store)
	{
		typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
		typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
		typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
		
		for(unsigned i = 0; i < length(store.alignedReadStore); ++i){
			if(store.alignedReadStore[i].beginPos == 28927578){
				TAlignedReadStoreElem a = store.alignedReadStore[i];
				TAlignedQualStoreElem q = store.alignQualityStore[a.id];
				int k = 0; k = 1;
			}
		}
	}
	
	// checks if the IDs in the alignedReadStore are continuously increasing
	// prints the two values and the position in the store if not
	template <typename TFragmentStore>
	inline void consistencyTest(TFragmentStore const & store)
	{
		typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
		typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
		typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
		
		for(unsigned i = 1; i < length(store.alignedReadStore); ++i){
			if(store.alignedReadStore[i-1].id != (store.alignedReadStore[i].id -1)){
				std::cout << i <<  ": " << store.alignedReadStore[i-1].id << " , " << store.alignedReadStore[i].id << "(i-1, i)\n";
				int k = 0; k = 1;
			}
		}
	}
	
	
	// BLOCK_STORE
/**
.Function.appendBlockStores:
..cat:Razers
..summary:Appends the aligned read and quality stores from the block stores to the main store given as first argument
..signature:appendBlockStores(store, blockStores, swiftPatternHandler, cnts, options, mode)
..param.store:@Class.FragmentStore@
..param.blockStores:@Class.String@ of @Class.FragmentStore@
..param.swiftPatternHandler:
..param.cnts:Counts
..param.options:RazerSOptions
..param.mode:RazerSMode
*/
	template <
		typename TFragmentStore,
		typename TPatternHandler,
		typename TCounts,
		typename TRazerSOptions,
		typename TRazerSMode >
	inline void
	appendBlockStores(
		TFragmentStore			& store,
		String<TFragmentStore>	& blockStores,
		TPatternHandler			& swiftPatternHandler,
		TCounts					& cnts,
		TRazerSOptions			& options,
		TRazerSMode const		& mode
	){
		typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
		typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
		typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
		
		// first: update the IDs and calculate new size
		TAlignedReadStoreSize oldSize = length(store.alignedReadStore);
		// so the prefix increment can be used in the loops
		TAlignedReadStoreSize sizeSum = --oldSize;

		for(unsigned i = 0; i < options.numberOfBlocks; ++i)
			for(unsigned j = 0; j < length(blockStores[i].alignedReadStore); ++j)
				blockStores[i].alignedReadStore[j].id = ++sizeSum;
		++sizeSum;
		
		// second: resize first so copying happens at most once and not every for each block in the worst case
		resize(store.alignedReadStore, sizeSum, Generous());
		resize(store.alignQualityStore, sizeSum, Generous());			

		// third: append single block stores
		for(unsigned i = 0; i < options.numberOfBlocks; ++i){
			for(unsigned j = 0; j < length(blockStores[i].alignedReadStore); ++j){
				store.alignedReadStore[++oldSize] = blockStores[i].alignedReadStore[j];
				store.alignQualityStore[oldSize] = blockStores[i].alignQualityStore[j];
			}
		}
		
		// fourth: compact matches
		if (length(store.alignedReadStore) > options.compactThresh)
		{
			oldSize = length(store.alignedReadStore);
			
			if (TYPECMP<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
				maskDuplicates(store, mode);	// overlapping parallelograms cause duplicates
			
			compactMatches(store, cnts, options, mode, swiftPatternHandler, COMPACT);
			
			if (options._debugLevel >= 2)
				::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
		}
		
	}
    
/**
.Function._mapSingleReadsToContig:
..cat:Razers
..summary:Appends the aligned read and quality stores from the block stores to the main store given as first argument
..signature:_mapSingleReadsToContig(store, contigId, swiftPatterns, preprocessingBlocks, cnts, orientation, options, mode)
..param.store:@Class.FragmentStore@
..param.contigId: the ID of the contig within the fragment store to which the reads are mapped
..param.swiftPatternHandler: Handler for swift pattern
..param.preprocessingBlocks: String of infixes of a string with bit vector patterns for each read
..param.cnts:Counts for statistics
..param.options:RazerSOptions
..param.mode:RazerSMode
*/
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
             ParallelSwiftPatternHandler<String<Pattern<TReadIndex, Swift<TSwiftSpec> > > > & swiftPatternHandler,
             TPreprocessing                                     & preprocessingBlocks,
             TCounts                                            & cnts,
             char                                                 orientation,			// q-gram index of reads
             TRazerSOptions                                     & options,
             TRazerSMode const                                  & mode)
    {
        SEQAN_CHECKPOINT
        
        // FILTRATION
        typedef typename TFragmentStore::TContigSeq				TContigSeq;
        typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
        typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;
        typedef typename Size<TSwiftFinder>::Type               TSwiftFinderSize;

        // HITS
        typedef typename TSwiftFinder::THitString               THitString;
		typedef typename Value<THitString>::Type                TSwiftHit;
        typedef typename Size<THitString>::Type                 THitStringSize;
        
        // VERIFICATION
        typedef typename Value<TPreprocessing>::Type                                                                        TPreprocessingBlock;
        typedef String<Pattern<TReadIndex, Swift<TSwiftSpec> > >                                                            TSwiftPatterns;
        typedef ParallelSwiftPatternHandler<TSwiftPatterns>                                                                 TSwiftPatternHandler;
        typedef MatchVerifier <TFragmentStore, TRazerSOptions, TRazerSMode, TPreprocessingBlock, TSwiftPatternHandler, TCounts >   TVerifier;
        typedef typename Fibre<TReadIndex, Fibre_Text>::Type                                                                TReadSet;
		typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type												TAlignedReadStoreSize;
        
        if (options._debugLevel >= 1)
        {
            ::std::cerr << ::std::endl << "Process genome seq #" << contigId;
            if (orientation == 'F') ::std::cerr << "[fwd]" << std::endl;
            else                    ::std::cerr << "[rev]" << std::endl;
        }
        
        lockContig(store, contigId);
        TContigSeq &contigSeq = store.contigStore[contigId].seq;
        if (orientation == 'R')	reverseComplementInPlace(contigSeq);
        
        // Finder and verifier strings of the same size as there are swift patterns
        // One Swift finder, Swift pattern and verifier work together
        String<TSwiftFinder> swiftFinders;
        resize(swiftFinders, options.numberOfBlocks, Exact());
        TSwiftFinder swiftFinder(contigSeq, options.repeatLength, 1);
        		
        String<TVerifier> verifier;
        resize(verifier, options.numberOfBlocks, Exact());
		
		// BLOCK_STORE
		// a temporary store for every block. reduces the critical region in pushing verified hit to an atomic ID increment
		String<TFragmentStore> blockStores;
		resize(blockStores, options.numberOfBlocks, Exact());
        
        for(int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
            // initialize finders
            swiftFinders[blockId] = swiftFinder;
			
			// BLOCK_STORE
            // initialize verifier
			TVerifier oneVerifier(blockStores[blockId], options, preprocessingBlocks[blockId], swiftPatternHandler, cnts);
            oneVerifier.onReverseComplement = (orientation == 'R');
            oneVerifier.genomeLength = length(contigSeq);
            oneVerifier.m.contigId = contigId;
            
            verifier[blockId] = oneVerifier;
        }
        
        bool beginOk = true;
        
        // set up finder
        for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
            beginOk = beginOk & windowFindBegin(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options.errorRate);
			if(not beginOk) break;
        }
        
        // if started correctly
        if(beginOk){
            unsigned blockSize = length(host(host(swiftPatternHandler.swiftPatterns[0])));
#ifdef RAZERS_TIMER
            CharString contigAndDirection;
            append(contigAndDirection, store.contigNameStore[contigId]);
            append(contigAndDirection, "-");
            append(contigAndDirection, orientation);
            
            // for waiting times
            String<_proFloat> waitingTimes;
            resize(waitingTimes, options.numberOfCores, Exact());
			for(unsigned coreId = 0; coreId < options.numberOfCores; ++coreId)
				waitingTimes[coreId] = sysTime();
#endif
            // use only core many threads
            omp_set_num_threads(options.numberOfCores);
            
            // go over contig sequence
			#pragma omp parallel 
            while(true)
            {
                bool stop = false;
                
                // Parallelize loop body each thread gets one at a time
                // As the data structures are split up alread and each thread works on only one element of
                // the strings (finder, patterns) at a time the variables can be shared
                // TODO: maybe try schedule(guided)
                #pragma omp for schedule(dynamic, 1)
                for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){ //TSwiftFinderSize

#ifdef RAZERS_TIMER					
                    // start time for filtering
                    _proFloat startTime = sysTime();

                    Pair<int, int> posLength(0, 0);
                    // filter window and save hits
                    stop = !windowFindNext(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options.windowSize);

                    // get time for filtering
                    _proFloat filteringTime = sysTime() - startTime;
                    // start time for verification
                    startTime = sysTime();
#else
                    stop = !windowFindNext(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options.windowSize);
#endif
                    // verify hits
                    THitString hits = getSwiftHits(swiftFinders[blockId]);
                    for(THitStringSize h = 0; h < length(hits); ++h){
                        verifier[blockId].m.readId = (blockId * blockSize) + hits[h].ndlSeqNo;         //array oder jedesmal berechnen
                        matchVerify(verifier[blockId], getSwiftRange(hits[h], contigSeq), hits[h].ndlSeqNo, host(host(swiftPatternHandler.swiftPatterns[blockId])), mode);
                        ++options.countFiltration;
                    }
#ifdef RAZERS_TIMER
                    // get time for filtering
                    _proFloat verificationTime = sysTime() - startTime;
                    
                    #pragma omp critical
                    {
                        std::cout << "timing>\t";
                        std::cout << omp_get_thread_num() << "\t";
                        std::cout << blockId << "\t";
                        std::cout << contigAndDirection << "\t";
                        std::cout << posLength.i1 << "\t";
                        std::cout << posLength.i2 << "\t";
                        std::cout << filteringTime << "\t";
                        std::cout << length(hits) << "\t";
                        std::cout << verificationTime << "\n";
                    }
                    // set waiting time
                    waitingTimes[omp_get_thread_num()] = sysTime();
                    
#endif
                }
                
#ifdef RAZERS_TIMER
                _proFloat now = sysTime();
                #pragma omp critical
                {
                    for(unsigned k = 0; k < length(waitingTimes); ++k){
                        std::cout << "waiting>\t"  << k << "\t" << (now - waitingTimes[k]) << "\n";
					}
                }
#endif
                
                if(stop) break;
            }
            
            // clear finders
            for (unsigned int blockId = 0; blockId < options.numberOfBlocks; ++blockId)
                windowFindEnd(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId]);
            
        } // end: if started correctly
		
		// BLOCK_STORE
		// append the alignedReadStore and the alignQualityStore from the blocks to the main store
		appendBlockStores(store, blockStores, swiftPatternHandler, cnts, options, mode);
		
        if (!unlockAndFreeContig(store, contigId))							// if the contig is still used
            if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation
    }
    
    // TODO: doc
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
        SEQAN_CHECKPOINT
        
        typedef FragmentStore<TFSSpec, TFSConfig>                   TFragmentStore;
        typedef typename IF<TYPECMP<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
        typedef typename Value<TReadIndexString>::Type              TReadIndex;
        
        // filter
        typedef Pattern<TReadIndex, Swift<TSwiftSpec> >             TSwiftPattern;
		typedef String<TSwiftPattern>								TSwiftPatterns;
		typedef ParallelSwiftPatternHandler<TSwiftPatterns>         TSwiftPatternHandler;
        
        // verifier
        typedef Pattern<TRead, MyersUkkonen>                        TMyersPattern;
        typedef typename Infix<String<TMyersPattern> >::Type        TVerifierBlock;
        typedef typename Position<TReadIndexString>::Type           TPos;

        // configure Swift patterns
        TSwiftPatterns swiftPatterns;
        
        resize(swiftPatterns, options.numberOfBlocks, Exact());
        for (TPos blockId = 0; blockId < options.numberOfBlocks; ++blockId) {
            assign(swiftPatterns[blockId].data_host, readIndices[blockId]);
            swiftPatterns[blockId].params.minThreshold = options.threshold;
            swiftPatterns[blockId].params.tabooLength = options.tabooLength;
			swiftPatterns[blockId].params.printDots = (blockId == 0) and (options._debugLevel > 0);
        }
		// use pattern handler instead of pattern to have access to the global read ID
        TSwiftPatternHandler swiftPatternHandler(swiftPatterns);
                
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
        resize(forwardPatternsBlocks, options.numberOfBlocks, Exact());
        TPos offSet = 0;
        for (TPos blockId = 0; blockId < options.numberOfBlocks; ++blockId){
            // block of the same size as the corresponding index (number of reads in this index)
            unsigned blockSize = length(host(readIndices[blockId]));
            TVerifierBlock block = infix(forwardPatterns, offSet, offSet + blockSize);            
            forwardPatternsBlocks[blockId] = block;
            offSet += blockSize;
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
        

#ifdef RAZERS_TIMER
		// print header line for timer
		std::cout << "timing>\tthread\tblock\tcontigAndDircetion\tpos\tlength\tfilter.time\tverifications\tverification.time\n";
		std::cout << "waiting>\tthread\ttime\n";
#endif
        
#ifndef RAZERS_PARALLEL_READS_WHOLE_GENOME		
        // iterate over genome sequences
        for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId)
        {
            // lock to prevent releasing and loading the same contig twice
            // (once per _mapSingleReadsToContig call)
            lockContig(store, contigId);
            if (options.forward)
                _mapSingleReadsToContig(store, contigId, swiftPatternHandler, forwardPatternsBlocks, cnts, 'F', options, mode);
            
            if (options.reverse)
                _mapSingleReadsToContig(store, contigId, swiftPatternHandler, forwardPatternsBlocks, cnts, 'R', options, mode);
            unlockAndFreeContig(store, contigId);
        }
#else
		// TODO: something is wrong here + need to reverse the genome only once
		// a temporary store for every block. reduces the critical region in pushing verified hit to an atomic ID increment
		String<TFragmentStore> blockStores;
		resize(blockStores, options.numberOfBlocks, Exact());
		
		// reduce the size at which the alignment store is compacted
		unsigned mainCompactThresh = options.compactThresh;
		options.compactThresh = mainCompactThresh / options.numberOfBlocks;
		
		# pragma omp parallel for //TODO: correct pragma
		for(int blockID = 0; blockID < (int)options.numberOfBlocks; ++blockID){
			// only first one in verbous
			if(blockID != 0)
				options._debugLevel = 0;
			
			// iterate over genome sequences
			for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId)
			{				
				// lock to prevent releasing and loading the same contig twice
				// (once per _mapSingleReadsToContig call)
				lockContig(store, contigId);
				if (options.forward)
					_mapSingleReadsToContigWindow(store, blockStores[blockID], contigId, swiftPatterns[blockID], swiftPatternHandler, forwardPatternsBlocks[blockID], cnts, 'F', options, mode);
				
				if (options.reverse)
					_mapSingleReadsToContigWindow(store, blockStores[blockID], contigId, swiftPatterns[blockID], swiftPatternHandler, forwardPatternsBlocks[blockID], cnts, 'R', options, mode);
				unlockAndFreeContig(store, contigId);
			}
		}
		
		// reset compact threshold
		options.compactThresh = mainCompactThresh; 
		// combine blocks
		appendBlockStores(store, blockStores, swiftPatternHandler, cnts, options, mode);
#endif
        
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
    
    // TODO: doc
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
        SEQAN_CHECKPOINT
        
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
        unsigned cores = options.numberOfCores;
#ifdef RAZERS_PARALLEL_READS_WHOLE_GENOME
		// does not make sence to split up further in this case
		options.blocksPerCore = 1;
#endif
        options.numberOfBlocks = cores * options.blocksPerCore;
        
        if (options._debugLevel >= 1){
            ::std::cerr << ::std::endl << "Number of cores:                 \t" << cores << std::endl;
            ::std::cerr <<                "Number of blocks:                \t" << options.numberOfBlocks << " (" << cores << " x " << options.blocksPerCore << ")" << std::endl;
        }
        
        // compare with noOfBlocks, there needs to be at least one read per block
        if(length(store.readSeqStore) < options.numberOfBlocks)
            options.numberOfBlocks = length(store.readSeqStore);
        
        // if there are not enough reads that the parallel version makes sence use the normal one
        if(length(store.readSeqStore) < 10) // TODO: usefull number
            return _mapSingleReads(store, cnts, options, shape, mode);
        else {
            // number of reads per block
            unsigned perBlock = length(store.readSeqStore) / options.numberOfBlocks;
            unsigned perBlockRemainder = length(store.readSeqStore) % options.numberOfBlocks;
            
            String<TIndex> indices;
            resize(indices, options.numberOfBlocks, Exact());
            
            int readID = 0;
            // create swift indices that can work in parallel
            for (unsigned blockID = 0; blockID < options.numberOfBlocks; ++blockID) {
                TReadSet readSet;
                
                // the last one gets some extra
                if((blockID == options.numberOfBlocks - 1) && (perBlockRemainder != 0)){
                    resize(readSet, perBlock + perBlockRemainder, Exact());
                } else {
                    resize(readSet, perBlock, Exact());
                }
                
                // get a subset of reads from the store
                for (unsigned i = 0; i < length(readSet); ++i) {
                    assign(readSet[i], store.readSeqStore[readID]);
                    ++readID;
                }
                
                // configure q-gram index
                intiIndex(indices[blockID], readSet, shape);
#ifdef RAZERS_OPENADDRESSING
				indices[blockID].alpha = options.loadFactor;
#endif
                cargo(indices[blockID]).abundanceCut = options.abundanceCut;
                cargo(indices[blockID])._debugLevel = options._debugLevel;
            }
            
            return _mapSingleReadsParallel(store, cnts, options, mode, indices);
        }
                
    }
    
}

#endif