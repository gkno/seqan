 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

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

#ifndef SEQAN_HEADER_RAZERS_MATEPAIRS_PARALLEL_H
#define SEQAN_HEADER_RAZERS_MATEPAIRS_PARALLEL_H

#include <seqan/misc/misc_dequeue.h>


namespace SEQAN_NAMESPACE_MAIN
{

// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.


//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence

template <
	typename TFSSpec,
	typename TFSConfig,
	typename TThreadId,
	typename TReadIndex,
	typename TSwiftSpec,
	typename TVerifier,
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode>
void goOverContig(
		FragmentStore<TFSSpec, TFSConfig>					& mainStore,
		FragmentStore<TFSSpec, TFSConfig>					& threadStore,
		unsigned											  contigId,
		TThreadId											& threadID,
		ParallelSwiftPatternHandler<String<
			Pattern<TReadIndex, Swift<TSwiftSpec> > > >		& swiftPatternHandlerL,
		ParallelSwiftPatternHandler<String<
			Pattern<TReadIndex, Swift<TSwiftSpec> > > >		& swiftPatternHandlerR,
		TVerifier											& verifierL,
		TVerifier											& verifierR,
		TCounts												& cnts,
		char												  orientation,
		TRazerSOptions										& options,
		TRazerSMode											& mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>				TFragmentStore;
	typedef typename TFragmentStore::TMatePairStore			TMatePairStore;
	typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
	typedef typename Value<TMatePairStore>::Type			TMatePair;
	typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
	typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
	typedef typename Fibre<TReadIndex, FibreText>::Type	TReadSet;
	typedef typename Id<TAlignedRead>::Type					TId;

	typedef typename TFragmentStore::TContigSeq				TGenome;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Position<TGenome>::Type				TGPos;
	typedef typename MakeSigned_<TGPos>::Type				TSignedGPos;
	typedef typename Infix<TGenome>::Type					TGenomeInf;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >			TSwiftFinderR;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;
	
	// MATE-PAIR FILTRATION
	typedef Triple<__int64, TAlignedRead, TAlignQuality>	TDequeueValue;
	typedef Dequeue<TDequeueValue>							TDequeue;
	typedef typename TDequeue::TIter						TDequeueIterator;

	// VERIFICATION
	typedef String<TSwiftPattern>							TSwiftPatterns;

	const unsigned NOT_VERIFIED = 1u << (8*sizeof(unsigned)-1);
	
	TGenome &genome = mainStore.contigStore[contigId].seq;
	
	TSwiftPattern & swiftPatternL = swiftPatternHandlerL.swiftPatterns[threadID];
	TSwiftPattern & swiftPatternR = swiftPatternHandlerR.swiftPatterns[threadID];
	
	TReadSet	&readSetL = host(host(swiftPatternL));
	//TReadSet	&readSetR = host(host(swiftPatternR));
		
	if (empty(readSetL))
		return;
	
	// distance <= libLen + libErr + 2*(parWidth-readLen) - shapeLen
	// distance >= libLen - libErr - 2*parWidth + shapeLen
	TSize readLength = length(readSetL[0]);
	TSignedGPos maxDistance = options.libraryLength + options.libraryError - 2 * (int)readLength - (int)length(indexShape(host(swiftPatternL)));
	TSignedGPos minDistance = options.libraryLength - options.libraryError + (int)length(indexShape(host(swiftPatternL)));
	TGPos scanShift = (minDistance < 0)? 0: minDistance;

	// exit if contig is shorter than library size
	if (length(genome) <= scanShift)
		return;
	
	TGenomeInf genomeInf = infix(genome, scanShift, length(genome));
	TSwiftFinderL swiftFinderL(genome, options.repeatLength, 1);
	TSwiftFinderR swiftFinderR(genomeInf, options.repeatLength, 1);
	
	TDequeue fifo;						// stores left-mate potential matches
	String<__int64> lastPotMatchNo;		// last number of a left-mate potential
	__int64 lastNo = 0;					// last number over all left-mate pot. matches in the queue
	__int64 firstNo = 0;				// first number over all left-mate pot. match in the queue
	Pair<TGPos> gPair;

	resize(lastPotMatchNo, length(host(swiftPatternL)), (__int64)-1, Exact());

	TSize gLength = length(genome);

	TAlignedRead mR;
	TAlignQuality qR;
	TDequeueValue fL(-1, mR, qR);	// to supress uninitialized warnings
	
//	unsigned const preFetchMatches = 2048;

	// iterate all verification regions returned by SWIFT
	while (find(swiftFinderR, swiftPatternR, options.errorRate)) 
	{
		unsigned matePairId = swiftPatternR.curSeqNo;
		unsigned globalMatePairId = matePairId + options.blockSize * threadID; // local ID to global ID
		
		TGPos rEndPos = endPosition(swiftFinderR) + scanShift;
		TGPos doubleParWidth = 2 * (*swiftFinderR.curHit).bucketWidth;
		
		// remove out-of-window left mates from fifo
		while (!empty(fifo) && (TSignedGPos)front(fifo).i2.endPos + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
		{
			popFront(fifo);
			++firstNo;
		}
		// add within-window left mates to fifo
		while (empty(fifo) || (TSignedGPos)back(fifo).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
		{
			if (find(swiftFinderL, swiftPatternL, options.errorRate, false))
			{
				gPair = positionRange(swiftFinderL);
				if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
				{
					// link in
					unsigned globalCurSeqNo = swiftPatternL.curSeqNo + options.blockSize * threadID; // local ID to global ID
					fL.i1 = lastPotMatchNo[swiftPatternL.curSeqNo];
					lastPotMatchNo[swiftPatternL.curSeqNo] = lastNo++;

					fL.i2.readId = mainStore.matePairStore[globalCurSeqNo].readId[0] | NOT_VERIFIED;
					fL.i2.beginPos = gPair.i1;
					fL.i2.endPos = gPair.i2;
					
					pushBack(fifo, fL);
				}
			} else
				break;
		}
		
		int	bestLeftScore = MinValue<int>::VALUE;
		int bestLibSizeError = MaxValue<int>::VALUE;
		TDequeueIterator bestLeft = TDequeueIterator();
		
		TDequeueIterator it;
		unsigned leftReadId = mainStore.matePairStore[globalMatePairId].readId[0];
		__int64 lastPositive = (__int64)-1;
		for (__int64 i = lastPotMatchNo[matePairId]; firstNo <= i; i = (*it).i1)
		{
			it = &value(fifo, i - firstNo);
			{
				// verify left mate (equal seqNo), if not done already
				if ((*it).i2.readId & NOT_VERIFIED)
				{
					if((TSignedGPos)(*it).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
					{
						TId globalReadId = (*it).i2.readId & ~NOT_VERIFIED;
						if (matchVerify(verifierL, infix(genome, (TSignedGPos)(*it).i2.beginPos, (TSignedGPos)(*it).i2.endPos), 
								globalReadId, mainStore.readSeqStore, mode))
						{
							verifierL.m.readId = globalReadId;		// has been verified positively
							(*it).i2 = verifierL.m;
							(*it).i3 = verifierL.q;
							
							// short-cut negative matches
							if (lastPositive == (__int64)-1)
								lastPotMatchNo[matePairId] = i;
							else
								value(fifo, lastPositive - firstNo).i1 = i;
							lastPositive = i;
						} else
							(*it).i2.readId = ~NOT_VERIFIED;				// has been verified negatively
					} else
						lastPositive = i;
				} else
					lastPositive = i;
				
				if ((*it).i2.readId == leftReadId)
				{
					int score = (*it).i3.score;
					if (bestLeftScore <= score)
					{
						int libSizeError = options.libraryLength - (int)((__int64)mR.endPos - (__int64)(*it).i2.beginPos);
						if (libSizeError < 0) libSizeError = -libSizeError;
						if (bestLeftScore == score)
						{
							if (bestLibSizeError > libSizeError)
							{
								bestLibSizeError = libSizeError;
								bestLeft = it;
							}
						}
						else
						{
							bestLeftScore = score;
							bestLibSizeError = libSizeError;
							bestLeft = it;
							if (bestLeftScore == 0) break;	// TODO: replace if we have real qualities
						}
					}
				}
			}
		}
		
		// short-cut negative matches
		if (lastPositive == (__int64)-1)
			lastPotMatchNo[matePairId] = (__int64)-1;
		else
			value(fifo, lastPositive - firstNo).i1 = (__int64)-1;
		
		// verify right mate, if left mate matches
		if (bestLeftScore != MinValue<int>::VALUE)
		{
			if (matchVerify(verifierR, infix(swiftFinderR), 
					mainStore.matePairStore[globalMatePairId].readId[1], mainStore.readSeqStore, mode))
			{
				// distance between left mate beginning and right mate end
				__int64 dist = (__int64)verifierR.m.endPos - (__int64)(*bestLeft).i2.beginPos;
				if (dist <= options.libraryLength + options.libraryError &&
					options.libraryLength <= dist + options.libraryError)
				{
					mR = verifierR.m;
					qR = verifierR.q;

					fL.i2 = (*bestLeft).i2;
					fL.i3 = (*bestLeft).i3;

					// transform mate readNo to global readNo
					TMatePair &mp = mainStore.matePairStore[globalMatePairId];
					fL.i2.readId = mp.readId[0];
					mR.readId    = mp.readId[1];

					// transform coordinates to the forward strand
					if (orientation == 'F')
					{
						TSize temp = mR.beginPos;
						mR.beginPos = mR.endPos;
						mR.endPos = temp;
					} else 
					{
						fL.i2.beginPos = gLength - fL.i2.beginPos;
						fL.i2.endPos = gLength - fL.i2.endPos;
						TSize temp = mR.beginPos;
						mR.beginPos = gLength - mR.endPos;
						mR.endPos = gLength - temp;
						dist = -dist;
					}
					
					// set a unique pair id
					fL.i2.pairMatchId = mR.pairMatchId = options.nextPairMatchId;
					options.nextPairMatchId += options.numberOfCores;
					if (options.nextPairMatchId == TAlignedRead::INVALID_ID)
						options.nextPairMatchId = 0;
					
					// score the whole match pair
					fL.i3.pairScore = qR.pairScore = fL.i3.score + qR.score;
					
					if (!options.spec.DONT_DUMP_RESULTS)
					{
						fL.i2.id = length(threadStore.alignedReadStore);
						appendValue(threadStore.alignedReadStore, fL.i2, Generous());
						appendValue(threadStore.alignQualityStore, fL.i3, Generous());
						mR.id = length(threadStore.alignedReadStore);
						appendValue(threadStore.alignedReadStore, mR, Generous());
						appendValue(threadStore.alignQualityStore, qR, Generous());
						
						if (length(threadStore.alignedReadStore) > options.compactThresh)
						{
							typename Size<TAlignedReadStore>::Type oldSize = length(threadStore.alignedReadStore);
//									maskDuplicates(matches);	// overlapping parallelograms cause duplicates
							compactPairMatches(mainStore, threadStore, cnts, options, swiftPatternHandlerL, swiftPatternHandlerR);
							
							if (length(threadStore.alignedReadStore) * 4 > oldSize)			// the threshold should not be raised
								options.compactThresh += (options.compactThresh >> 1);	// if too many matches were removed
							
							if (options._debugLevel >= 2)
								::std::cerr << '(' << oldSize - length(threadStore.alignedReadStore) << " matches removed)";
						}
					}
					++options.countVerification;
				}
				++options.countFiltration;
			}
		} 	
	}
}

template <
	typename TFSSpec,
	typename TFSConfig,
	typename TReadIndex,
	typename TSwiftSpec,
	typename TPreprocessing,
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode>
void _mapMatePairReadsParallel(
	FragmentStore<TFSSpec, TFSConfig>					& store,
	String<FragmentStore<TFSSpec, TFSConfig> >			& threadStores,
	unsigned											  contigId,				// ... and its sequence number
	ParallelSwiftPatternHandler<String<
		Pattern<TReadIndex, Swift<TSwiftSpec> > > >		& swiftPatternHandlerL,
	ParallelSwiftPatternHandler<String<
		Pattern<TReadIndex, Swift<TSwiftSpec> > > >		& swiftPatternHandlerR,
	TPreprocessing										& preprocessing,
	TCounts												& cnts,
	char												  orientation,
	TRazerSOptions										& options,
	String<TRazerSOptions>								& threadOptions,
	TRazerSMode											& mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>				TFragmentStore;
	typedef typename TFragmentStore::TMatePairStore			TMatePairStore;
	typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
	typedef typename Value<TMatePairStore>::Type			TMatePair;
	typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
	typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
	typedef typename Fibre<TReadIndex, FibreText>::Type	TReadSet;
	typedef typename Id<TAlignedRead>::Type					TId;

	typedef typename TFragmentStore::TContigSeq				TGenome;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Position<TGenome>::Type				TGPos;
	typedef typename MakeSigned_<TGPos>::Type				TSignedGPos;
	typedef typename Infix<TGenome>::Type					TGenomeInf;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >			TSwiftFinderR;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;
	
	// MATE-PAIR FILTRATION
	typedef Triple<__int64, TAlignedRead, TAlignQuality>	TDequeueValue;
	typedef Dequeue<TDequeueValue>							TDequeue;
	typedef typename TDequeue::TIter						TDequeueIterator;

	// VERIFICATION
	typedef String<TSwiftPattern>							TSwiftPatterns;
	typedef ParallelSwiftPatternHandler<TSwiftPatterns>		TSwiftPatternHandler;
	typedef MatchVerifier <TFragmentStore, TRazerSOptions,
		TRazerSMode, TPreprocessing,
		TSwiftPatternHandler, TCounts >						TVerifier;
	
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}
	
	lockContig(store, contigId);
	TGenome &genome = store.contigStore[contigId].seq;
	if (orientation == 'R')	reverseComplement(genome);
		
	// Create a verifier for each thread. This way each thread gets its own store to dump the matches in.
	String<TVerifier> verifierL, verifierR;
	resize(verifierL, options.numberOfCores, Exact());
	resize(verifierR, options.numberOfCores, Exact());
	
	for(int threadId = 0; threadId < (int)options.numberOfCores; ++threadId){
		// initialize verifier
		TVerifier oneVerifierL(threadStores[threadId], threadOptions[threadId], preprocessing, swiftPatternHandlerL, cnts);
		oneVerifierL.oneMatchPerBucket = true;
		oneVerifierL.m.contigId = contigId;
		// assign it to the string
		verifierL[threadId] = oneVerifierL;
		
		// initialize verifier
		TVerifier oneVerifierR(threadStores[threadId], threadOptions[threadId], preprocessing, swiftPatternHandlerR, cnts);
		oneVerifierR.oneMatchPerBucket = true;
		oneVerifierR.m.contigId = contigId;
		// assign it to the string
		verifierR[threadId] = oneVerifierR;
	}
		
	#pragma omp parallel default(none) \
	shared(options, store, threadStores, contigId, swiftPatternHandlerL, swiftPatternHandlerR, verifierL, verifierR, cnts, orientation, threadOptions, mode)
	{
	#pragma omp for
		for (int threadID = 0; threadID < (int)options.numberOfCores; ++threadID) {
			goOverContig(store, threadStores[threadID], contigId, threadID, swiftPatternHandlerL, swiftPatternHandlerR, 
				verifierL[threadID], verifierR[threadID], cnts, orientation, threadOptions[threadID], mode);
		}
	}
		
	if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
		if (orientation == 'R')	reverseComplement(genome);	// we have to restore original orientation
	
}


template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TCounts,
	typename TSpec, 
	typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
	typename TReadIndexString>
int _mapMatePairReadsParallelCreatePatterns(
		FragmentStore<TFSSpec, TFSConfig>					& store,
		TCounts												& cnts,
		RazerSOptions<TSpec>								& options,
		RazerSMode<TAlignMode, TGapMode, TScoreMode>  const & mode,
		TReadIndexString									& readIndicesL, // TODO const? 
		TReadIndexString									& readIndicesR)
{
	SEQAN_CHECKPOINT

	typedef FragmentStore<TFSSpec, TFSConfig>					TFragmentStore;
	typedef typename If<
				IsSameType<TGapMode,RazerSGapped>::VALUE,
	 			SwiftSemiGlobal, 
				SwiftSemiGlobalHamming>::Type 					TSwiftSpec;
	typedef typename Value<TReadIndexString>::Type				TReadIndex;

	// filter
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >				TSwiftPattern;
	typedef String<TSwiftPattern>								TSwiftPatterns;
	typedef ParallelSwiftPatternHandler<TSwiftPatterns>			TSwiftPatternHandler;

	// verifier
	typedef Pattern<TRead, MyersUkkonen>						TMyersPattern;
	// typedef Pattern<TRead, Myers<FindInfix, False, void> >		TMyersPattern;
	typedef typename Position<TReadIndexString>::Type			TPos;
	
	//unsigned pairCount = length(store.matePairStore);
	
	// A temporary store for every thread. They are combined after all contigs are processed
	String<TFragmentStore> threadStores;
	resize(threadStores, options.numberOfCores, Exact());

	// configure Swift patterns
	String<TSwiftPattern> swiftPatternsL, swiftPatternsR;
	resize(swiftPatternsL, options.numberOfCores, Exact());
	resize(swiftPatternsR, options.numberOfCores, Exact());
	
	// configure Swift
	for (TPos threadID = 0; threadID < options.numberOfCores; ++threadID) {
		assign(swiftPatternsL[threadID].data_host, readIndicesL[threadID]);
		swiftPatternsL[threadID].params.minThreshold = options.threshold;
		swiftPatternsL[threadID].params.tabooLength = options.tabooLength;
		swiftPatternsL[threadID].params.printDots = 0; // only one should print the dots
		
		assign(swiftPatternsR[threadID].data_host, readIndicesR[threadID]);
		swiftPatternsR[threadID].params.minThreshold = options.threshold;
		swiftPatternsR[threadID].params.tabooLength = options.tabooLength;
		swiftPatternsR[threadID].params.printDots = (threadID == 0) and (options._debugLevel > 0);
	}
	
	// use pattern handler instead of pattern to have access to the global read ID
	TSwiftPatternHandler swiftPatternHandlerL(swiftPatternsL);
	TSwiftPatternHandler swiftPatternHandlerR(swiftPatternsR);
	
	// init edit distance verifiers
	String<TMyersPattern> forwardPatterns;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (options.gapMode == RAZERS_GAPPED)
	{
		unsigned readCount = length(store.readSeqStore);
		resize(forwardPatterns, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatterns[i], store.readSeqStore[i]);
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}
	
	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;
	SEQAN_PROTIMESTART(find_time);
	
	String<RazerSOptions<TSpec> > threadOptions; 
	resize(threadOptions, options.numberOfCores, options, Exact());
	for (TPos threadID = 0; threadID < options.numberOfCores; ++threadID)
		threadOptions[threadID].nextPairMatchId = threadID;
	
	
	for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId) {
		// lock to prevent releasing and loading the same contig twice
		// (once per _mapSingleReadsToContig call)
		lockContig(store, contigId);
		
		if (options.forward)
			_mapMatePairReadsParallel(store, threadStores, contigId, 
				swiftPatternHandlerL, swiftPatternHandlerR,
				forwardPatterns, cnts, 'F', options, threadOptions, mode);
		
		if (options.reverse)
			_mapMatePairReadsParallel(store, threadStores, contigId, 
				swiftPatternHandlerL, swiftPatternHandlerR,
				forwardPatterns, cnts, 'R', options, threadOptions, mode);
		
		unlockAndFreeContig(store, contigId);
	}

	// TODO: combine stores
	// TODO: get filter and verify count from thread options
	// TODO: remember matepair id

	// Get the matches from different thread stores and write them in the main store
	appendBlockStores(store, threadStores, options);
	for(int i = 0; i < (int) options.numberOfCores; ++i){
		options.countVerification += threadOptions[i].countVerification;
	}

	// output for verbose and very verbose options
	options.timeMapReads = SEQAN_PROTIMEDIFF(find_time);
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:  " << options.countFiltration << ::std::endl;
		::std::cerr << "Verfication counter: " << options.countVerification << ::std::endl;
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
int _mapMatePairReadsParallel(
	FragmentStore<TFSSpec, TFSConfig>					& store,
	TCounts												& cnts,
	RazerSOptions<TSpec>								& options,
	TShape const										& shape,
	RazerSMode<TAlignMode, TGapMode, TScoreMode>  const & mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	typedef typename Position<TReadSeqStore>::Type		TReadSeqStorePos;
	
	typedef typename Value<TReadSeqStore>::Type			TRead;
	typedef StringSet<TRead>							TReadSet;
#ifndef RAZERS_OPENADDRESSING
	typedef Index<TReadSet, IndexQGram<TShape> >	TIndex;			// q-gram index
#else
	typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >	TIndex;
#endif

	if (options._debugLevel >= 1){
		::std::cerr << ::std::endl << "Number of cores:                 \t" << options.numberOfCores << std::endl;
	}
	
	unsigned pairCount = length(store.matePairStore);
	
	// compare with noOfBlocks, there needs to be at least one read per block
	if(options.numberOfCores == 1 or pairCount < options.numberOfCores or pairCount < 100)
		return _mapMatePairReads(store, cnts, options, shape, mode);
	else {
		// number of reads per block
		options.blockSize = pairCount / options.numberOfCores;
		unsigned perBlockRemainder = pairCount % options.numberOfCores;
		
		// TODO: create two strings of indices
		String<TIndex> indicesL, indicesR;
		resize(indicesL, options.numberOfCores, Exact());
		resize(indicesR, options.numberOfCores, Exact());

		TReadSeqStorePos readID = 0;
		// create swift indices that can work in parallel
		for (unsigned threadID = 0; threadID < options.numberOfCores; ++threadID) {
			TReadSet readSetL, readSetR;

			// the last one gets some extra
			if((threadID == options.numberOfCores - 1) && (perBlockRemainder != 0)){
				resize(readSetL, options.blockSize + perBlockRemainder, Exact());
				resize(readSetR, options.blockSize + perBlockRemainder, Exact());
			}
			else{
				resize(readSetL, options.blockSize, Exact());
				resize(readSetR, options.blockSize, Exact());
			}
			
			// split mate-pairs over two indices
			for (unsigned i = 0; i < length(readSetL); ++i) {
				assign(readSetL[i], store.readSeqStore[store.matePairStore[readID].readId[0]]);
				assign(readSetR[i], store.readSeqStore[store.matePairStore[readID].readId[1]]);
				++readID;
			}
			reverseComplement(readSetR);
			
			// configure q-gram indices
			intiIndex(indicesL[threadID], readSetL, shape);
			intiIndex(indicesR[threadID], readSetR, shape);
#ifdef RAZERS_OPENADDRESSING
			indicesL[threadID].alpha = options.loadFactor;
			indicesR[threadID].alpha = options.loadFactor;
#endif
			reverse(indexShape(indicesR[threadID]));		// right mate qualities are reversed -> reverse right shape
			
			cargo(indicesL[threadID]).abundanceCut = options.abundanceCut;
			cargo(indicesL[threadID])._debugLevel = options._debugLevel;
			cargo(indicesR[threadID]).abundanceCut = options.abundanceCut;
			cargo(indicesR[threadID])._debugLevel = options._debugLevel;
		}
		return _mapMatePairReadsParallelCreatePatterns(store, cnts, options, mode, indicesL, indicesR);
	}
}


} // End namespace

#endif
