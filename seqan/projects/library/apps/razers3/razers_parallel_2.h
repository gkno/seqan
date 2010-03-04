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
    
    // TODO: doc
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
        unsigned noOfBlocks = length(swiftPatterns);
        
        String<TSwiftFinder> swiftFinders;
        resize(swiftFinders, noOfBlocks, Exact());
        TSwiftFinder swiftFinder(contigSeq, options.repeatLength, 1);
        
        TSwiftPatternHandler swiftPatternHandler(swiftPatterns);
        String<TVerifier> verifier;
        resize(verifier, noOfBlocks, Exact());
        
        for(unsigned i = 0; i < noOfBlocks; ++i){
            // initialize finders
            swiftFinders[i] = swiftFinder;
            
            // initialize verifier
            TVerifier oneVerifier(store, options, preprocessingBlocks[i], swiftPatternHandler, cnts);
            oneVerifier.onReverseComplement = (orientation == 'R');
            oneVerifier.genomeLength = length(contigSeq);
            oneVerifier.m.contigId = contigId;
            
            verifier[i] = oneVerifier;
        }
        
        bool beginOk = true;
        
        // set up finder
        for (unsigned i = 0; i < noOfBlocks; ++i){
            // FIXME: break after first one is false
            beginOk = beginOk & windowFindBegin(swiftFinders[i], swiftPatterns[i], options.errorRate, options._debugLevel);
        }
        
        // if started correctly
        if(beginOk){
            unsigned blockSize = length(host(host(swiftPatterns[0])));
            int numberOfBlocks = length(swiftFinders);
            // FIXME: int cores = omp_get_num_procs();
#ifdef RAZERS_TIMER
            CharString contigAndDirection;
            append(contigAndDirection, store.contigNameStore[contigId]);
            append(contigAndDirection, "-");
            append(contigAndDirection, orientation);
            
            // for waiting times
            String<_proFloat> waitingTimes;
            resize(waitingTimes, cores, Exact());
#endif
            // use only core many threads
            // FIXME: omp_set_num_threads(cores);
            
            // go over contig sequence
            while(true)
            {
                bool stop = false;
                
                // Parallelize loop body each thread gets one at a time
                // As the data structures are split up alread and each thread works on only one element of
                // the strings (finder, patterns) at a time the variables can be shared
                // TODO: maybe try schedule(guided)
                // FIXME: #pragma omp parallel for schedule(dynamic, 1) //private(options)
                for (int i = 0; i < numberOfBlocks; ++i){ //TSwiftFinderSize

#ifdef RAZERS_TIMER
                    // start time for filtering
                    _proFloat startTime = sysTime();

                    Pair<int, int> posLength(0, 0);
                    // filter window and save hits
                    stop = !windowFindNext(swiftFinders[i], swiftPatterns[i], options.windowSize, (options._debugLevel & (i == 0)), posLength);

                    // get time for filtering
                    _proFloat filteringTime = sysTime() - startTime;
                    // start time for verification
                    startTime = sysTime();
#else
                    stop = !windowFindNext(swiftFinders[i], swiftPatterns[i], options.windowSize, (options._debugLevel & (i == 0)));
#endif
                    // verify hits
                    THitString hits = getSwiftHits(swiftFinders[i]);
                    for(THitStringSize h = 0; h < length(hits); ++h){
                        verifier[i].m.readId = (i * blockSize) + hits[h].ndlSeqNo;         //array oder jedesmal berechnen
                        matchVerify(verifier[i], getSwiftRange(hits[h], contigSeq), hits[h].ndlSeqNo, host(host(swiftPatterns[i])), mode);
                        ++options.countFiltration;
                    }
#ifdef RAZERS_TIMER
                    // get time for filtering
                    _proFloat verificationTime = sysTime() - startTime;
                    
                    #pragma omp critical
                    {
                        std::cout << "timing>\t";
                        std::cout << omp_get_thread_num() << "\t";
                        std::cout << i << "\t";
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
                    std::cout << "waiting>";
                    for(unsigned k = 0; k < length(waitingTimes); ++k){
                        std::cout << "\t" << (now - waitingTimes[k]);
                    }
                    std::cout << "\n";
                }
#endif
                
                if(stop) break;
            }
            
            // clear finders
            for (unsigned i = 0; i < length(swiftFinders); ++i)
                windowFindEnd(swiftFinders[i], swiftPatterns[i]);
            
        } // end: if started correctly
                
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
        
        // verifier
        typedef Pattern<TRead, MyersUkkonen>                        TMyersPattern;
        typedef typename Infix<String<TMyersPattern> >::Type        TVerifierBlock;
        typedef typename Position<TReadIndexString>::Type           TPos;

        // configure Swift patterns
        typedef String<TSwiftPattern> TSwiftPatterns;
        TSwiftPatterns swiftPatterns;
        
        resize(swiftPatterns, length(readIndices), Exact());
        for (TPos i = 0; i < length(readIndices); ++i) {
            assign(swiftPatterns[i].data_host, readIndices[i]);
            swiftPatterns[i].params.minThreshold = options.threshold;
            swiftPatterns[i].params.tabooLength = options.tabooLength;
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
            unsigned blockSize = length(host(readIndices[i]));
            TVerifierBlock block = infix(forwardPatterns, offSet, offSet + blockSize);            
            forwardPatternsBlocks[i] = block;
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
#endif
        
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
        unsigned cores = 1;//omp_get_num_procs();
        unsigned noOfBlocks = cores * options.blocksPerCore;
        
        if (options._debugLevel >= 1){
            ::std::cerr << ::std::endl << "Number of cores:                 \t" << cores << std::endl;
            ::std::cerr <<                "Number of blocks:                \t" << noOfBlocks << std::endl;
        }
        
        // compare with noOfBlocks, there needs to be at least one read per block
        if(length(store.readSeqStore) < noOfBlocks)
            noOfBlocks = length(store.readSeqStore);
        
        // if there are not enough reads that the parallel version makes sence use the normal one
        if(length(store.readSeqStore) < 10) // TODO: usefull number
            return _mapSingleReads(store, cnts, options, shape, mode);
        else {
            // number of reads per block
            unsigned perBlock = length(store.readSeqStore) / noOfBlocks;
            unsigned perBlockRemainder = length(store.readSeqStore) % noOfBlocks;
            
            String<TIndex> indices;
            resize(indices, noOfBlocks, Exact());
            
            int readID = 0;
            // create swift indices that can work in parallel
            for (unsigned b = 0; b < noOfBlocks; ++b) {
                TReadSet readSet;
                
                // the last one gets some extra
                if((b == noOfBlocks - 1) && (perBlockRemainder != 0)){
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
                intiIndex(indices[b], readSet, shape);
                cargo(indices[b]).abundanceCut = options.abundanceCut;
                cargo(indices[b])._debugLevel = options._debugLevel;
            }
            
            return _mapSingleReadsParallel(store, cnts, options, mode, indices);
        }
                
    }
    
}

#endif