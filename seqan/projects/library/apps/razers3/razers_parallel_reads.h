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

#ifndef SEQAN_HEADER_RAZERS_PARALLEL_READS_H
#define SEQAN_HEADER_RAZERS_PARALLEL_READS_H

namespace SEQAN_NAMESPACE_MAIN
{
	
#ifndef _OPENMP
int omp_get_thread_num(){
	return 2;
}
#endif


// Pattern handler: is necessary to be able to calculate the overall read ID based on the local 
// read ID within one thread. Used in matchVerify.
template <typename TSwiftPatterns>
struct ParallelSwiftPatternHandler
{
	// SEQAN_CHECKPOINT

	TSwiftPatterns &swiftPatterns;
	ParallelSwiftPatternHandler(TSwiftPatterns &_swiftPatterns):
	swiftPatterns(_swiftPatterns) {}
};


// Initializes a q-gram index. Does the same as the constructor with the same parameters, 
// which cannot be called for a string of indices.
template<
	typename TReadSet,
	typename TShape>
void intiIndex(
		Index<TReadSet, Index_QGram<TShape, OpenAddressing> > & index,
		TReadSet & _text,
		TShape const & _shape)
{
	SEQAN_CHECKPOINT

	value(index.text) = _text;
	index.shape = _shape;
}


// Uses the read ID to find the correct SWIFT pattern in the string handled by the ParallelSwiftPatternHandler,
// and the correct local ID within this SWIFT pattern to update the max errors.
template < typename TSwiftPatterns, typename TReadNo, typename TMaxErrors >
inline void 
setMaxErrors(ParallelSwiftPatternHandler<TSwiftPatterns> &swift, TReadNo readNo, TMaxErrors maxErrors)
{
	SEQAN_CHECKPOINT

	int blockSize = length(host(host(swift.swiftPatterns[0])));
	int indexNo = readNo / blockSize;
	int localReadNo = readNo % blockSize;

	int minT = _qgramLemma(swift.swiftPatterns[indexNo], localReadNo, maxErrors);
	if (minT > 1){
		if (maxErrors < 0) minT = SupremumValue<int>::VALUE;
		setMinThreshold(swift.swiftPatterns[indexNo], localReadNo, (unsigned)minT);
	}
}


// Appends the content fo the aligned read / quality stores in the blockStores to the counterparts in the (main) store.
template <
	typename TFragmentStore,
	typename TRazerSOptions>
inline void
appendBlockStores(
		TFragmentStore			& store,
		String<TFragmentStore>	& blockStores,
		TRazerSOptions			& options)
{
	typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
	typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
	
	// first: update the IDs and calculate new size
	TAlignedReadStoreSize oldSize = length(store.alignedReadStore);
	// so the prefix increment can be used in the loops
	TAlignedReadStoreSize sizeSum = --oldSize;
	
	for(unsigned i = 0; i < options.numberOfCores; ++i)
		for(unsigned j = 0; j < length(blockStores[i].alignedReadStore); ++j)
			blockStores[i].alignedReadStore[j].id = ++sizeSum;
	++sizeSum;
	
	// second: resize first so copying happens at most once and not every for each block in the worst case
	resize(store.alignedReadStore, sizeSum, Generous());
	resize(store.alignQualityStore, sizeSum, Generous());			
	
	// third: append single block stores
	for(unsigned i = 0; i < options.numberOfCores; ++i){
		for(unsigned j = 0; j < length(blockStores[i].alignedReadStore); ++j){
			store.alignedReadStore[++oldSize] = blockStores[i].alignedReadStore[j];
			store.alignQualityStore[oldSize] = blockStores[i].alignQualityStore[j];
		}
		
		clear(blockStores[i].alignedReadStore);
		clear(blockStores[i].alignQualityStore);		
		
	}
	
}


inline void updateTav(String<int> & tav, int const myId){
	SEQAN_CHECKPOINT
	
	int k = 1;
	int n = (int)length(tav);
	// Search for the next block that is still working
	while(tav[(myId + k) % n] == -1)
		++k;
	// If there is none (all checked).
	if(k == n)
		return;
	
	tav[(myId + k) % n] += tav[myId];
	tav[myId] = -1;
	
	return;
}


template<typename TSwiftHit>
struct SwiftHitComparison:
	public ::std::binary_function<TSwiftHit, TSwiftHit, bool>
{
	inline bool
	operator() (TSwiftHit const & i1, TSwiftHit const & i2) {
		return (i1.ndlSeqNo) < (i2.ndlSeqNo);
	}
};


template<
	typename TPosition,
	typename TSpec,
	typename THstkPos,
	typename TAlignMode,
	typename TErrors>
void partitionHits(
		String<TPosition>									& positions,
		String<_SwiftHit<Tag<
			_SwiftSemiGlobal<TSpec> >, THstkPos> > 			& hits,
		int const											  noOfParts,
		RazerSMode<TAlignMode, RazerSGapped, TErrors> const)
{
	typedef String<_SwiftHit<Tag<_SwiftSemiGlobal<TSpec> >, THstkPos> >		THitString;
	typedef typename Value<THitString>::Type									TSwiftHit;
	typedef typename Iterator<THitString>::Type								THitStringIter;
	
	THitStringIter i1 = begin(hits);
	THitStringIter i2 = end(hits);
	std::sort(i1, i2, SwiftHitComparison<TSwiftHit>());

	resize(positions, noOfParts + 1, Exact());
	positions[0] = 0;
	positions[noOfParts] = length(hits);
	
	TPosition partSize = length(hits) / noOfParts;
	TPosition myPos;
	
	for (int i = 1; i < noOfParts+1; ++i){
		myPos = i * partSize;
		
		unsigned now = 0;
		unsigned last = hits[myPos - 1].ndlSeqNo;
		
		for(; myPos < length(hits); ++myPos){
			now = hits[myPos].ndlSeqNo;
			if(last != now){
 				break;
			}
			
			last = now;
		}
		positions[i] = myPos;
	}
}


template<
	typename TPosition,
	typename TSpec,
	typename THstkPos,
	typename TAlignMode,
	typename TErrors>
void partitionHits(
		String<TPosition>									& positions,
		String<_SwiftHit<Tag<
			_SwiftSemiGlobal<TSpec> >, THstkPos> > const	& hits,
		int const											  noOfParts,
		RazerSMode<TAlignMode, RazerSUngapped, TErrors> const)
{
	resize(positions, noOfParts + 1, Exact());
	
	TPosition partSize = length(hits) / noOfParts;
	
	positions[0] = 0;
	for (int i = 1; i < noOfParts; ++i){
		positions[i] = i * partSize;
	}
	positions[noOfParts] = length(hits);
}


template<
	typename TVerifier,
	typename THitString,
	typename TReadId,
	typename TContigSeq,
	typename TReadSet,
	typename TRazerSMode>
inline void verifyHits(
		TVerifier									& verifier,
		THitString									& hits,
		typename Position<THitString>::Type const	  startPos,
		typename Position<THitString>::Type const	  endPos,
		TReadId const								  offSet,
		TContigSeq									& contigSeq,
		TReadSet									& readSet,
		TRazerSMode									& mode)
{
	typedef typename Position<THitString>::Type		THitStringPos;
	
	for(THitStringPos h = startPos; h < endPos; ++h){
		TReadId absReadId = offSet + hits[h].ndlSeqNo;
		verifier.m.readId = absReadId;
		
		matchVerify(verifier, swiftInfix(hits[h], contigSeq), absReadId, readSet, mode);
	}
}


template <
	typename TContigSeq, 
	typename TReadIndex, 
	typename TSwiftSpec,
	typename TVerifier,
	typename TCounts,
	typename TRazerSOptions,
	typename TFragmentStore,
	typename TRazerSMode >
inline void goOverContig(
		ParallelSwiftPatternHandler<String<
			Pattern<TReadIndex, Swift<TSwiftSpec> > > >	& swiftPatternHandler,
		String<Finder<TContigSeq, Swift<TSwiftSpec> > >	& swiftFinders,
		String<TVerifier>								& verifier,
		TCounts											& ,//cnts,
		TRazerSOptions									& options,
		String<TFragmentStore>							& threadStores,
		TFragmentStore									& store,
		TRazerSMode										& mode)
   {
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >		TSwiftFinder;
	typedef typename TSwiftFinder::THitString			THitString;
	typedef typename Value<THitString>::Type			TSwiftHit;
	typedef typename Position<THitString>::Type			THitStringPos;
	typedef typename Iterator<THitString>::Type			THitStringIter;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	typedef typename TFragmentStore::TAlignedReadStore	TAlignedReadStore;
	typedef typename Size<TAlignedReadStore>::Type		TAlignedReadStoreSize;
	
	// Needed because the omp shared clause does not allow for "." in the variabel names.
	TReadSeqStore const & readSet = store.readSeqStore;
	TContigSeq const & contigSeq = host(swiftFinders[0]);
	
	// Number of Threads Allowed for Verification. Shared by all blocks
	String<int> tav;
	fill(tav, options.numberOfBlocks, 1, Exact());
	
	#pragma omp parallel num_threads((int)options.numberOfCores)
	{
	#pragma omp master
	{		
		for(int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId)
		{
			#pragma omp task default(none) \
				shared(tav, verifier, swiftFinders, swiftPatternHandler, contigSeq, readSet, threadStores, options, mode) \
				firstprivate(blockId)
			{
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				#ifdef RAZERS_PARALLEL_TIMER
				_proFloat myTime = sysTime();
				#endif
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				
				bool sequenceLeft = true;
				while(sequenceLeft){
					sequenceLeft = windowFindNext(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], 
						options.windowSize);
						
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					#ifdef RAZERS_PARALLEL_TIMER
					printf("filter: %f sec (pos: %lu, thread: %d)\n", sysTime() - myTime, swiftFinders[blockId].curPos, blockId);
					myTime = sysTime();
					#endif
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					
					THitString & hits = getSwiftHits(swiftFinders[blockId]);
					options.countFiltration += length(hits);
					
					if(length(hits) == 0)
						continue;
					#pragma omp flush(tav)
					
					// if(true){ 
					if(tav[blockId] == 1 or (int) length(hits) < tav[blockId] * 3){
						verifyHits(verifier[blockId], hits, 0, length(hits),
							(blockId * options.blockSize), host(swiftFinders[0]), readSet, mode);
					}
					else {
						String<THitStringPos> positions;
						partitionHits(positions, hits, tav[blockId], mode);
						// verify
						for(int relId = 0; relId < tav[blockId]; ++relId)
						{
							#pragma omp task default(none) \
								shared(verifier, hits, positions, contigSeq, readSet, threadStores, options, mode) \
								firstprivate(blockId, relId)
							{
								// calculate the absolute Id 
								int absVerifyId = (options.numberOfBlocks + blockId - relId) % options.numberOfBlocks;
								
								verifyHits(verifier[absVerifyId], hits, positions[relId], positions[relId + 1],
									(blockId * options.blockSize), contigSeq, readSet, mode);
							}
						}
						#pragma omp taskwait
						
					} // End else
					
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					#ifdef RAZERS_PARALLEL_TIMER
					#pragma omp critical(printf)
					printf("verify: %f sec (thread: %d)\n", sysTime() - myTime, blockId);
					#endif
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				} // End while
				
				#pragma omp critical(update_tav)
				updateTav(tav, blockId);
				
			} // End task
		} // End for
		
		#pragma omp taskwait
		
	} // End Master
	} // End Parallel
		
	// Clear finders
	for(int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
		windowFindEnd(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId]);
	}
	
}


template <
	typename TFragmentStore, 
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TPreprocessing,
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode >
   void _mapSingleReadsToContig(
	TFragmentStore										& store,
	String<TFragmentStore>								& threadStores,
	int													  contigId,
	ParallelSwiftPatternHandler<String<
		Pattern<TReadIndex, Swift<TSwiftSpec> > > >		& swiftPatternHandler,
	TPreprocessing										& preprocessing,
	TCounts												& cnts,
	char												  orientation,
	TRazerSOptions										& options,
	String<TRazerSOptions>								& threadOptions,
	TRazerSMode const									& mode)
   {
	SEQAN_CHECKPOINT

	// FILTRATION
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;
	typedef typename Size<TSwiftFinder>::Type               TSwiftFinderSize;

	// HITS
	typedef typename TSwiftFinder::THitString				THitString;
	typedef typename Value<THitString>::Type				TSwiftHit;
	typedef typename Size<THitString>::Type					THitStringSize;

	// VERIFICATION
	typedef String<Pattern<TReadIndex, Swift<TSwiftSpec> > >TSwiftPatterns;
	typedef ParallelSwiftPatternHandler<TSwiftPatterns>		TSwiftPatternHandler;
	typedef MatchVerifier <TFragmentStore, TRazerSOptions,
		TRazerSMode, TPreprocessing,
		TSwiftPatternHandler, TCounts >						TVerifier;
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
	typedef typename Size<TAlignedReadStore>::Type			TAlignedReadStoreSize;

	// for verbose options
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F') ::std::cerr << "[fwd]" << std::endl;
		else                    ::std::cerr << "[rev]" << std::endl;
	}
	
	// lock contig
	lockContig(store, contigId);
	TContigSeq &contigSeq = store.contigStore[contigId].seq;
	// if (orientation == 'R')	reverseComplementInPlace(contigSeq);
	if (orientation == 'R'){
		reverseComplementInPlace(contigSeq);
	}
	
	// Finder and verifier strings of the same size as there are swift patterns
	// One Swift finder, Swift pattern and verifier work together
	String<TSwiftFinder> swiftFinders;
	resize(swiftFinders, options.numberOfBlocks, Exact());
	// The finders are the same for each block as they only depend on the reference.
	// Separate ones are needed nevertheless as they form a pair with the patterns.
	TSwiftFinder swiftFinder(contigSeq, options.repeatLength, 1);
	for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId)
		swiftFinders[blockId] = swiftFinder;
	
	// Create a verifier for each thread. This way each thread gets its own store to dump the matches in.
	// As consequence the dumping does not need to be critical
	String<TVerifier> verifier;
	resize(verifier, options.numberOfBlocks, Exact());
	
	for(int threadId = 0; threadId < (int)options.numberOfBlocks; ++threadId){
		// initialize verifier
		TVerifier oneVerifier(threadStores[threadId], threadOptions[threadId], preprocessing, swiftPatternHandler, cnts);
		oneVerifier.onReverseComplement = (orientation == 'R');
		oneVerifier.genomeLength = length(contigSeq);
		oneVerifier.m.contigId = contigId;
		// assign it to the string
		verifier[threadId] = oneVerifier;
	}
	
	// Set up finder. beginOK is true after the loop if all finders are set up successfully.
	bool beginOk = true;
	
	for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
		beginOk = windowFindBegin(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], 
			options.errorRate);
	}	
	
	// Only if the finders are set up.
	if(beginOk){
		goOverContig(swiftPatternHandler, swiftFinders, verifier, cnts, options, threadStores, store, mode);
	}
	
	if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
	if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation
	
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
int _mapSingleReadsParallelCreatePatterns(
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
	// typedef Pattern<TRead, MyersUkkonen>                        TMyersPattern;
	typedef Pattern<TRead, Myers<FindInfix, False, void> >		TMyersPattern;
	typedef typename Infix<String<TMyersPattern> >::Type        TVerifierBlock;
	typedef typename Position<TReadIndexString>::Type           TPos;

	// A temporary store for every block. They are combined after all contigs are processed
	String<TFragmentStore> threadStores;
	resize(threadStores, options.numberOfBlocks, Exact());

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
	unsigned readCount = length(store.readSeqStore);
	String<TMyersPattern> forwardPatterns;
	resize(forwardPatterns, readCount, Exact());

	options.compMask[4] = (options.matchN)? 15: 0;
	if (options.gapMode == RAZERS_GAPPED){
		resize(forwardPatterns, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i){
			setHost(forwardPatterns[i], store.readSeqStore[i]);
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}
	
	if (options.maqMapping){
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
	
	String<RazerSOptions<TSpec> > threadOptions; 
	fill(threadOptions, options.numberOfBlocks, options, Exact());
	
	// iterate over genome sequences
	for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId){
		// lock to prevent releasing and loading the same contig twice
		// (once per _mapSingleReadsToContig call)
		lockContig(store, contigId);
		if (options.forward){
			_mapSingleReadsToContig(store, threadStores, contigId, swiftPatternHandler, forwardPatterns, cnts, 'F', options, threadOptions, mode);
		}
		if (options.reverse){
			_mapSingleReadsToContig(store, threadStores, contigId, swiftPatternHandler, forwardPatterns, cnts, 'R', options, threadOptions, mode);
		}
		unlockAndFreeContig(store, contigId);
	}

	// Get the matches from different thread stores and write them in the main store
	appendBlockStores(store, threadStores, options);
	for(int i = 0; i < (int) options.numberOfBlocks; ++i){
		options.countVerification += threadOptions[i].countVerification;
	}

	
	// output for verbose and very verbose options
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
		RazerSMode<TAlignMode, TGapMode, TScoreMode>  const	& mode)
{
	SEQAN_CHECKPOINT

	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type						TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, Index_QGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

	// number of cores and how many blocks there should be
	// a block is a subset of reads that is filtered and verified in a row
	// depending on how long it takes to process the individual blocks a single
	// thread might work through more than otheres
	unsigned cores = options.numberOfCores;

	options.numberOfBlocks = cores;

	if (options._debugLevel >= 1){
		::std::cerr << ::std::endl << "Number of cores:                 \t" << cores << std::endl;
	}

	// compare with noOfBlocks, there needs to be at least one read per block
	if(length(store.readSeqStore) < options.numberOfBlocks)
		options.numberOfBlocks = length(store.readSeqStore);

	// if there are not enough reads that the parallel version makes sence use the normal one
	if(length(store.readSeqStore) < 100 || options.numberOfBlocks == 1)
		return _mapSingleReads(store, cnts, options, shape, mode);
	else {
		// number of reads per block
		unsigned perBlock = length(store.readSeqStore) / options.numberOfBlocks;
		options.blockSize = perBlock;
		unsigned perBlockRemainder = length(store.readSeqStore) % options.numberOfBlocks;

		String<TIndex> indices;
		resize(indices, options.numberOfBlocks, Exact());

		int readID = 0;
		// create swift indices that can work in parallel
		for (unsigned blockID = 0; blockID < options.numberOfBlocks; ++blockID) {
			TReadSet readSet;

			// the last one gets some extra
			if((blockID == options.numberOfBlocks - 1) && (perBlockRemainder != 0))
				resize(readSet, perBlock + perBlockRemainder, Exact());
			else
				resize(readSet, perBlock, Exact());

			// get a subset of reads from the store
			for (unsigned i = 0; i < length(readSet); ++i) {
				assign(readSet[i], store.readSeqStore[readID]);
				++readID;
			}
		
			// configure q-gram indices
			intiIndex(indices[blockID], readSet, shape);
			#ifdef RAZERS_OPENADDRESSING
			indices[blockID].alpha = options.loadFactor;
			#endif
			cargo(indices[blockID]).abundanceCut = options.abundanceCut;
			cargo(indices[blockID])._debugLevel = options._debugLevel;
			// build index
			//indexRequire(indices[blockID], QGram_SADir());
		}
		return _mapSingleReadsParallelCreatePatterns(store, cnts, options, mode, indices);
	}	
}
	
} // End namespace

#endif
