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

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

// #if defined (RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS) || defined (RAZERS_PARALLEL_READS_INDEPENDENT)
// #include <parallel/algorithm>
// # endif

namespace SEQAN_NAMESPACE_MAIN
{
	
#ifndef _OPENMP
int omp_get_thread_num(){
	return 0;
}
#endif

/**
.Class.ParallelSwiftPatternHandler:
..summary:Holds a string of @Spec.Swift@ @Class.Pattern@s to allow the access of the overall read IDs.
..cat:Razers
..signature:ParallelSwiftPatternHandler<TSwiftPatterns>
..param.TSwiftPatterns:The patterns type. With the @Class.String@ around.
*/
template <typename TSwiftPatterns>
struct ParallelSwiftPatternHandler
{
	// SEQAN_CHECKPOINT

	TSwiftPatterns &swiftPatterns;
	ParallelSwiftPatternHandler(TSwiftPatterns &_swiftPatterns):
	swiftPatterns(_swiftPatterns) {}
};


/**
.Function.intiIndex:
..cat:Razers
..summary:Connects an empty @Spec.Index_QGram@ with a host text and a shape.
..signature:intiIndex(index, _text, _shape)
..param.index:@Spec.Index_QGram@
..param._text:In case of RazerS a @Class.StringSet@ of the reads
..param._shape:@Class.Shape@
*/
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
		//::std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
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
	
	
/**
.Function.appendBlockStores:
..cat:Razers
..summary:Appends the aligned read and quality stores from the block stores to the main store given as first argument
..signature:appendBlockStores(store, blockStores, swiftPatternHandler, cnts, options, mode)
..param.store:@Class.FragmentStore@
..param.blockStores:@Class.String@ of @Class.FragmentStore@
..param.swiftPatternHandler:@Class.ParallelSwiftPatternHandler@
..param.cnts:Counts
..param.options:RazerSOptions
..param.mode:RazerSMode
*/
template <
	typename TFragmentStore,
	// typename TPatternHandler,
	// typename TCounts,
	typename TRazerSOptions//,
	// typename TRazerSMode 
	>
inline void
appendBlockStores(
		TFragmentStore			& store,
		String<TFragmentStore>	& blockStores,
		// TPatternHandler			& swiftPatternHandler,
		// TCounts					& cnts,
		TRazerSOptions			& options//,
		// TRazerSMode const		& mode
		)
{
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
	// if (length(store.alignedReadStore) > options.compactThresh)
	// 	{
	// 		oldSize = length(store.alignedReadStore);
	// 		
	// 		if (TYPECMP<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
	// 			maskDuplicates(store, mode);	// overlapping parallelograms cause duplicates
	// 		
	// 		compactMatches(store, cnts, options, mode, swiftPatternHandler, COMPACT);
	// 		
	// 		if (options._debugLevel >= 2)
	// 			::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
	// 	}
	
}
	
	
template <
	typename TFragmentStore,
	typename TPatternHandler,
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode >
inline void
_appendBlockStore(
		TFragmentStore			& mainStore,
		TFragmentStore			& blockStore,
		TPatternHandler			& swiftPatternHandler,
		TCounts					& cnts,
		TRazerSOptions			& options,
		TRazerSMode const		& mode
		)
{
	typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
	// first: update the IDs
	TAlignedReadStoreSize id1 = length(mainStore.alignedReadStore) - 1;
	TAlignedReadStoreSize id2 = id1;
	
	for(unsigned i = 0; i < length(blockStore.alignedReadStore); ++i)
		blockStore.alignedReadStore[i].id = ++id1;
	
	// second: append block store and clear it
	++id1;
	resize(mainStore.alignedReadStore, id1, Generous());
	resize(mainStore.alignQualityStore, id1, Generous());

	for(unsigned i = 0; i < length(blockStore.alignedReadStore); ++i){
		mainStore.alignedReadStore[++id2] = blockStore.alignedReadStore[i];
		mainStore.alignQualityStore[id2] = blockStore.alignQualityStore[i];
	}

	clear(blockStore.alignedReadStore);
	clear(blockStore.alignQualityStore);
	
	// third: compact matches
	if (length(mainStore.alignedReadStore) > options.compactThresh)
	{
		TAlignedReadStoreSize oldSize = length(mainStore.alignedReadStore);
		
		if (TYPECMP<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
			maskDuplicates(mainStore, mode);	// overlapping parallelograms cause duplicates
		
		compactMatches(mainStore, cnts, options, mode, swiftPatternHandler, COMPACT);
		
		if (options._debugLevel >= 2)
			::std::cerr << '(' << oldSize - length(mainStore.alignedReadStore) << " matches removed)";
	}
}


/**
.Function.appendBlockStore:
..cat:Razers
..summary:Appends the aligned read and quality stores from the block stores to the main store given as first argument
..signature:appendBlockStores(store, blockStores, swiftPatternHandler, cnts, options, mode)
..param.store:@Class.FragmentStore@
..param.blockStores:@Class.String@ of @Class.FragmentStore@
..param.swiftPatternHandler:@Class.ParallelSwiftPatternHandler@
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
appendBlockStore(
		bool					& appendStoreLocked,
		TFragmentStore			& mainStore,
		TFragmentStore			& blockStore,
		TPatternHandler			& swiftPatternHandler,
		TCounts					& cnts,
		TRazerSOptions			& options,
		TRazerSMode const		& mode
		)
{
	// Request lock.
	bool iLocked = false;
	#pragma omp critical(appendStore)
	{
		if(not appendStoreLocked){
			appendStoreLocked = true;
			iLocked = true;
		}
	}
	
	// Free lock.
	if(iLocked){
		_appendBlockStore(mainStore, blockStore, swiftPatternHandler, cnts, options, mode);
		
		#pragma omp critical(appendStore)
		appendStoreLocked = false;
		iLocked = false;
	}
}


#ifdef RAZERS_PARALLEL_READS_INDEPENDENT


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


// End RAZERS_PARALLEL_READS_INDEPENDENT
#endif

#if defined (RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS) || defined (RAZERS_PARALLEL_READS_INDEPENDENT)


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
	int stage,
	typename TSize>
inline TSize toBucket(TSize id){
	SEQAN_CHECKPOINT
	
	// return (id >> (4*stage)) & 15;
	return (id >> (8*stage)) & 255;
}
	
	
template<
	int stage,
	typename THit>
inline void myRadixPass(
		String<THit> 			& sorted,
		String<THit> const 		& unsorted)
{
	SEQAN_CHECKPOINT
	
	typedef typename Position<String<THit> >::Type		TPos;
	typedef String<TPos>						TBucket;
	typedef String<TBucket>								TBuckets;
	typedef typename Position<TBuckets>::Type			TBucketPos;
	
	TPos threads = 16;
	TPos noOfBuckets = 256;
	
	TBucket emptyBucket;
	fill(emptyBucket, noOfBuckets, 0, Exact());
	TBuckets buckets;
	fill(buckets, threads, emptyBucket, Exact());
	
	TPos partSize = length(unsorted) / threads;
	for(TPos p = 0; p < threads; ++p){
		#pragma omp task shared(buckets, unsorted) if(length(unsorted) > 10)
		{
			TPos endPos = (p == threads-1) ? length(unsorted) : (p+1)*partSize;
			for(TPos i = p*partSize; i < endPos; ++i){
				++buckets[p][toBucket<stage>(unsorted[i].ndlSeqNo)];
			}
		}
		
	}
	#pragma omp taskwait
	
	// Exclusive partial sum
	TPos sum1 = -1, sum2 = -1;
	
	for(TBucketPos i = 0; i < noOfBuckets; ++i){
		for(TPos p = 0; p < threads; ++p){
			sum2 += buckets[p][i];
			buckets[p][i] = sum1;
			sum1 = sum2;
		}
	}
	
	for(TPos p = 0; p < threads; ++p){
		#pragma omp task shared(buckets, sorted, unsorted) if(length(unsorted) > 10)
		{
			TPos endPos = (p == threads-1) ? length(unsorted) : (p+1)*partSize;
			for(TPos i = p*partSize; i < endPos; ++i){
				sorted[++buckets[p][toBucket<stage>(unsorted[i].ndlSeqNo)]] = unsorted[i];
			}
		}
	}
	#pragma omp taskwait
}


template<typename THit>
inline void myRadixSort(String<THit> & hits)
{
	SEQAN_CHECKPOINT
	
	String<THit> sorted;
	resize(sorted, length(hits), Exact());
			
	myRadixPass<0>(sorted, hits);
	myRadixPass<1>(hits, sorted);
	myRadixPass<2>(sorted, hits);
	myRadixPass<3>(hits, sorted);
	// If more stages are added be aware that hits needs to be the first argument in the last call
	// or sorted needs to be copied in hits.
}


// End (RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS) || defined (RAZERS_PARALLEL_READS_INDEPENDENT)
#endif

#ifdef RAZERS_PARALLEL_READS_INDEPENDENT


template<
	typename TPosition,
	typename TSpec,
	typename THstkPos>
void partitionHits(
		String<TPosition>									& positions,
		String<_SwiftHit<Tag<
			_SwiftSemiGlobal<TSpec> >, THstkPos> > const	& hits,
		int const											  noOfParts)
{ 
	resize(positions, noOfParts + 1, Exact());
	positions[0] = 0;
	positions[noOfParts] = length(hits);
	
	TPosition partSize = length(hits) / noOfParts;
	TPosition myPos;
	
	for (int blockId = 1; blockId < noOfParts+1; ++blockId){
		myPos = blockId * partSize;
		
		unsigned now = 0;
		unsigned last = hits[myPos - 1].ndlSeqNo;;
		
		for(; myPos < length(hits); ++myPos){
			now = hits[myPos].ndlSeqNo;
			if(last != now)
				break;
			
			last = now;
		}
		positions[blockId] = myPos;
	}
}


template<
	typename TVerifier,
	typename THitString,
	typename TReadId,
	typename TContigSeq,
	typename TReadSet,
	typename TRazerSOptions,
	typename TRazerSMode>
inline void verifyHits(
		TVerifier									& verifier,
		THitString									& hits,
		typename Position<THitString>::Type const	  startPos,
		typename Position<THitString>::Type const	  endPos,
		TReadId const								  offSet,
		TContigSeq									& contigSeq,
		TReadSet									& readSet,
		TRazerSOptions								& options,
		TRazerSMode									& mode)
{
	typedef typename Position<THitString>::Type		THitStringPos;
	
	__int64 verified = 0;
	
	for(THitStringPos h = startPos; h < endPos; ++h){
		TReadId absReadId = offSet + hits[h].ndlSeqNo;
		verifier.m.readId = absReadId;
		
		if(matchVerify(verifier, swiftInfix(hits[h], contigSeq), absReadId, readSet, mode))
			++verified;
	}
	
	#pragma omp atomic
	options.countVerification += verified;
	
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
inline void goOverContigIndependent(
		ParallelSwiftPatternHandler<String<
			Pattern<TReadIndex, Swift<TSwiftSpec> > > >	& swiftPatternHandler,
		String<Finder<TContigSeq, Swift<TSwiftSpec> > >	& swiftFinders,
		String<TVerifier>								& verifier,
		TCounts											& cnts,
		TRazerSOptions									& options,
		String<TFragmentStore>							& ,//threadStores,
		TFragmentStore									& store,
		TRazerSMode										& mode)
   {
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >		TSwiftFinder;
	typedef typename TSwiftFinder::THitString			THitString;
	typedef typename Value<THitString>::Type			TSwiftHit;
	typedef typename Position<THitString>::Type			THitStringPos;
	typedef typename Iterator<THitString>::Type			THitStringIter;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	
	// Needed because the omp shared clause does not allow for "." in the variabel names.
	TReadSeqStore const & readSet = store.readSeqStore;
	TContigSeq const & contigSeq = host(swiftFinders[0]);
	
	// Number of Threads Allowed for Verification. Shared by all blocks
	// Keeps track of how many threads can be used for verification per block.
	String<int> tav;
	fill(tav, options.numberOfBlocks, 1, Exact());
	
	// Spin lock for appending new matches to the main store.
	bool appendStoreLocked = false;
	
	#pragma omp parallel num_threads((int)options.numberOfCores)
	{
	#pragma omp master
	{
		// For each block.
		for(int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId)
		{
			#pragma omp task default(none) \
				shared(appendStoreLocked, tav, store, verifier, swiftFinders, swiftPatternHandler, contigSeq, readSet, options, cnts, mode) \
				firstprivate(blockId)
			{
				
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
				#ifdef FLEX_TIMER
				printf("%d start task\n", blockId);
				_proFloat myTime = sysTime();
				#endif
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				
				// Iterate over the sequence.
				bool sequenceLeft = true;
				while(sequenceLeft){
					
					// Filter sequence
					sequenceLeft = windowFindNext(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], 
						options.windowSize);
					
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					#ifdef FLEX_TIMER
					#pragma omp critical(printf)
					printf("%d filtering took %f, finder at: %lu\n", blockId, sysTime() - myTime, swiftFinders[blockId].curPos);
					myTime = sysTime();
					#endif
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					
					THitString & hits = getSwiftHits(swiftFinders[blockId]);
					if(length(hits) == 0)
						continue;
					#pragma omp flush(tav)
					
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					#ifdef FLEX_TIMER
					#pragma omp critical(printf)
					printf("%d tav: %d\n", blockId, tav[blockId]);
					#endif
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					
					// Verify found hits. // TODO: find an upper bound for the length of hits to verify in parallel or not.
					// #ifdef RAZERS_WORKSREALING
					if(tav[blockId] == 1 or (int) length(hits) < tav[blockId]){
					// #else
					// if(true){
					// #endif
						verifyHits(verifier[blockId], hits, 0, length(hits),
							(blockId * options.blockSize), host(swiftFinders[0]), readSet, options, mode);
						
						appendBlockStore(appendStoreLocked, store, *verifier[blockId].store, swiftPatternHandler, cnts, options, mode);
					}
					else {
						// First, verify hits.
						int partSize = length(hits) / tav[blockId];
						for(int relId = 0; relId < (tav[blockId] - 1); ++relId)
						{
							#pragma omp task default(none) \
								shared(verifier, hits, contigSeq, readSet, options, mode) \
								firstprivate(partSize, blockId, relId)
							{
								// calculate the absolute Id 
								int absVerifyId = (options.numberOfBlocks + blockId - relId) % options.numberOfBlocks;
								
								verifyHits(verifier[absVerifyId], hits, (relId * partSize), ((relId+1) * partSize),
									(blockId * options.blockSize), contigSeq, readSet, options, mode);
							}
						}
						// No task for the last part. Saves scheduling time.
						int absVerifyId = (options.numberOfBlocks + blockId - (tav[blockId] - 1)) % options.numberOfBlocks;
						
						verifyHits(verifier[absVerifyId], hits, ((tav[blockId] - 1) * partSize), (tav[blockId] * partSize),
							(blockId * options.blockSize), contigSeq, readSet, options, mode);
						
						#pragma omp taskwait
						
						// Second, append hits to main store.
						for(int relId = 0; relId < tav[blockId]; ++relId)
						{
							// calculate the absolute Id 
							int absVerifyId = (options.numberOfBlocks + blockId - relId) % options.numberOfBlocks;
							appendBlockStore(appendStoreLocked, store, *verifier[absVerifyId].store, swiftPatternHandler, cnts, options, mode);
						}
					} // End else
					
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
					#ifdef FLEX_TIMER
					#pragma omp critical(printf)
					printf("%d verifying took %f\n", blockId, sysTime() - myTime);
					#endif
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

				} // End while
				
				#pragma omp critical(update_tav)
				updateTav(tav, blockId);
				
				// TODO: remove
				// #pragma omp critical(printf)
				// {
				// 	printf("%d done; tav: ", blockId);
				// 	for(int i = 0; i < (int)length(tav); ++i)
				// 		printf("%d, ", tav[i]);
				// 	printf("\n");
				// }
				
			} // End task
		} // End for
		
		#pragma omp taskwait
		
	} // End Master
	} // End Parallel
	
	// Clear thread stores and finders
	for(int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
		windowFindEnd(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId]);
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
inline void goOverContigIndependent2(
		ParallelSwiftPatternHandler<String<
			Pattern<TReadIndex, Swift<TSwiftSpec> > > >	& swiftPatternHandler,
		String<Finder<TContigSeq, Swift<TSwiftSpec> > >	& swiftFinders,
		String<TVerifier>								& verifier,
		TCounts											& ,//cnts,
		TRazerSOptions									& options,
		String<TFragmentStore>							& ,//threadStores,
		TFragmentStore									& store,
		TRazerSMode										& mode)
{
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >		TSwiftFinder;
	typedef typename TSwiftFinder::THitString			THitString;
	typedef typename Value<THitString>::Type			TSwiftHit;
	typedef typename Position<THitString>::Type			THitStringPos;
	typedef typename Iterator<THitString>::Type			THitStringIter;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	
	// TODO: (size of aligned read store)?
	
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
				shared(tav, verifier, swiftFinders, swiftPatternHandler, contigSeq, readSet, options, mode) \
				firstprivate(blockId)
			{
				// TODO: remove
				// printf("%d start task\n", blockId);
				// _proFloat myTime = sysTime();
				
				bool sequenceLeft = true;
				while(sequenceLeft){
					sequenceLeft = windowFindNext(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], 
						options.windowSize);
						
					// TODO: remove
					// #pragma omp critical(printf)
					// printf("%d filtering took %f, finder at: %lu\n", blockId, sysTime() - myTime, swiftFinders[blockId].curPos);
					// myTime = sysTime();
					
					THitString & hits = getSwiftHits(swiftFinders[blockId]);
					if(length(hits) == 0)
						continue;
					#pragma omp flush(tav)
					
					// TODO: remove
					// #pragma omp critical(printf)
					// printf("%d tav: %d\n", blockId, tav[blockId]);
					
					if(true){ //(tav[blockId] == 1 or (int) length(hits) < tav[blockId]){
						verifyHits(verifier[blockId], hits, 0, length(hits),
							(blockId * options.blockSize), host(swiftFinders[0]), readSet, options, mode);
					}
					else {
						// TODO: remove
						// myTime = sysTime();
						// sort
						myRadixSort(hits);
						
						// TODO: remove
						// #pragma omp critical(printf)
						// printf("%d tav: %d\n", blockId, tav[blockId]);
						// myTime = sysTime();
						
						// THitStringIter i1 = begin(hits);
						// THitStringIter i2 = end(hits);
						// __gnu_parallel::sort(i1, i2, SwiftHitComparison<TSwiftHit>());
						// partition
						String<THitStringPos> positions;
						partitionHits(positions, hits, tav[blockId]);
						// verify
						// TODO: no task for first ...
						for(int relId = 0; relId < tav[blockId]; ++relId)
						{
							#pragma omp task default(none) \
								shared(verifier, hits, positions, contigSeq, readSet, options, mode) \
								firstprivate(blockId, relId)
							{
								// calculate the absolute Id 
								int absVerifyId = (options.numberOfBlocks + blockId - relId) % options.numberOfBlocks;
								
								verifyHits(verifier[absVerifyId], hits, positions[relId], positions[relId + 1],
									(blockId * options.blockSize), contigSeq, readSet, options, mode);
							}
						}
						#pragma omp taskwait
					} // End else
					
					// TODO: remove
					// #pragma omp critical(printf)
					// printf("%d verifying took %f\n", blockId, sysTime() - myTime);
					// myTime = sysTime();
				} // End while
				
				#pragma omp critical(update_tav)
				updateTav(tav, blockId);
				
				// TODO: remove
				#pragma omp critical(printf)
				{
					printf("%d done; tav: ", blockId);
					for(int i = 0; i < (int)length(tav); ++i)
						printf("%d, ", tav[i]);
					printf("\n");
				}
				
			} // End task
		} // End for
		
		#pragma omp taskwait
		
	} // End Master
	} // End Parallel
	
	// Get the matches from different thread stores and write them in the main store
	// appendBlockStores(store, threadStores, swiftPatternHandler, cnts, options, mode);
	
	// Clear thread stores and finders
	for(int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
		// clear(threadStores[blockId].alignedReadStore);
		// clear(threadStores[blockId].alignQualityStore);
		
		windowFindEnd(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId]);
	}
	
}


// END RAZERS_PARALLEL_READS_INDEPENDENT
#endif

#ifdef RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS

	
/**
.Class.ConvertedSwiftHit:
..summary:Representation of a Swift hit containing the corresponding read ID and infix in the reference.
..cat:Razers
..signature:ConvertedSwiftHit<TText>
..param.TText:Text on which the @Spec.Infix@ is based on.
*/
	template <typename TText>
	struct ConvertedSwiftHit
	{
		typedef typename Infix<TText>::Type TInfix;

		// static const unsigned accuracy = 100;

		unsigned	readID;
		__int64		begin;
		__int64		end;
		
		
		ConvertedSwiftHit(){}
		
		template <typename TSpec, typename THstkPos, typename TText2>
		ConvertedSwiftHit(_SwiftHit<Tag<_SwiftSemiGlobal<TSpec> >, THstkPos> & other, TText2 & ref, unsigned const & offSet):
			readID(offSet + other.ndlSeqNo)
		{
			typedef typename Infix<TText2>::Type TInfix;
			TInfix i = swiftInfix(other, ref);
			
			begin = beginPosition(i);
			end = endPosition(i);
		}

	};
	
	// template <typename TText>
	// bool 
	// ConvertedSwiftHit::operator< (const ConvertedSwiftHit<TText> & i2){
	// 	return readID < i2.readID;
	// }
	
/**
.Class.ConvertedSwiftHitComparison:
..summary:Functor to compare two @Class.ConvertedSwiftHit@ by their read ID with a certain accuracy.
..cat:Razers
..signature:ConvertedSwiftHitComparison
..remarks:For the comparison the read IDs are divied by the accuracy value and the results are then compared.
*/
	template<typename TText>
	struct ConvertedSwiftHitComparison:
		public ::std::binary_function<ConvertedSwiftHit<TText>, ConvertedSwiftHit<TText>, bool>
	{
		
		const unsigned accuracy;
		
/**
.Memfunc.ConvertedSwiftHitComparison#ConvertedSwiftHitComparison:
..class:Class.ConvertedSwiftHitComparison
..summary:Constructor
..signature:ConvertedSwiftHitComparison(_accuracy)
..param._accuracy:unsigned.Value by which the read IDs are divided before they are compared
*/
		ConvertedSwiftHitComparison(unsigned _accuracy):
			accuracy(_accuracy)
		{}
		
		bool operator() (ConvertedSwiftHit<TText> const & i1, ConvertedSwiftHit<TText> const & i2) {
			return (i1.readID / accuracy) < (i2.readID / accuracy);
		}
	};

	template<typename TText>
	struct ConvertedSwiftHitComparison2 :
		public ::std::binary_function<ConvertedSwiftHit<TText>, ConvertedSwiftHit<TText>, bool>
	{

		inline bool
		operator() (ConvertedSwiftHit<TText> const & i1, ConvertedSwiftHit<TText> const & i2) {
			return (i1.readID) < (i2.readID);
		}
	};


/**
.Function._collectHits:
..cat:Razers
..summary:Searches for new hits with the finder and pattern.
Stops when the finder reaches its end or the threshold of total hits is surpassed.
..signature:_collectHits(hits, totalHits, threshold, finder, pattern, windowSize)
..param.hits:String of @Class.ConvertedSwiftHit@ to store the found hits
..param.totalHits:Variable shared by all threads to know when to stop.
..param.threshold:The function stops searching for new hits when totalHits surpasses this value.
..param.finder:Swift finder.
..param.pattern:Swift pattern.
..param.windowSize:Window size @Function.windowFindNext@ is called with.
*/
	template <
		typename TSwiftHitString,
		typename TSwiftFinder,
		typename TSwiftPattern,
		typename TRazerSOptions>
	bool collectHitsForBlock(
			TSwiftHitString				& hits,
			unsigned					& totalHits,
			unsigned const				  offSet,
			TSwiftFinder				& finder,
			TSwiftPattern				& pattern,
			TRazerSOptions const		& options
#ifdef FLEX_TIMER
			, int						& windowNumbers
#endif
			)
	{
		typedef typename Size<TSwiftHitString>::Type	TSize; 
		
		// Get new hits till all threads together have enough.
		bool sequenceLeft = true;
		while(totalHits < options.collect and sequenceLeft){
			// Search for new hits.
			sequenceLeft = windowFindNext(finder, pattern, options.windowSize);
			append(hits, getSwiftHits(finder));
			// Update the overall hit count so that the other threads know when to stop as well.
			#pragma omp atomic
			totalHits += length(getSwiftHits(finder));
			
#ifdef FLEX_TIMER
			++windowNumbers;
#endif
		}
		
		for(TSize i = 0; i < length(hits); ++i)
			hits[i].ndlSeqNo += offSet;
		
		return sequenceLeft;
	}
	
	// TODO:doku
	template <
		typename TContigSeq, 
		typename TReadIndex, 
		typename TSwiftSpec,
		typename THits,
		typename TBlockSeqLeft,
		typename TRazerSOptions>
	bool collectHits(
			ParallelSwiftPatternHandler<String<
				Pattern<TReadIndex, Swift<TSwiftSpec> > > >	& swiftPatternHandler,
			String<Finder<TContigSeq, Swift<TSwiftSpec> > >	& swiftFinders,
			THits											& hits,
			TBlockSeqLeft									& blockSeqLeft,
			TRazerSOptions									& options)
	{
		// Number of hits collected. Shared by all threads. Updated in an atomic expression.
		unsigned totalHits = 0;
		
		// Is set true if any of the blocks has sequence left. Value is returned by the function.
		bool sequenceLeft = false;
		
#ifdef FLEX_TIMER
		String<_proFloat> waitingTimes;
		resize(waitingTimes, options.numberOfBlocks, Exact());
		for(unsigned coreId = 0; coreId < options.numberOfBlocks; ++coreId)
			waitingTimes[coreId] = sysTime();
		
		String<int> windowNumbers;
		fill(windowNumbers, options.numberOfBlocks, 0, Exact());
#endif
		
		// Create tasks that each collect hits for one block
		for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
			// Only use block if there is sequnce left.
			// Otherwise findWindowNext does unpredicted things.
			if(blockSeqLeft[blockId]){
#ifndef FLEX_TIMER
				#pragma omp task default(none) \
						shared(blockSeqLeft, hits, totalHits, options, swiftFinders, swiftPatternHandler, sequenceLeft) \
						firstprivate(blockId)
				{
					blockSeqLeft[blockId] = collectHitsForBlock(hits[blockId], totalHits, blockId*options.blockSize,
						swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options);

					// If any of the block has sequnece left the value will be true in the end.
					#pragma omp atomic
					sequenceLeft |= blockSeqLeft[blockId];
				} // End task.
#else
				#pragma omp task default(none) \
						shared(blockSeqLeft, hits, totalHits, options, swiftFinders, swiftPatternHandler, windowNumbers, waitingTimes, sequenceLeft) \
						firstprivate(blockId)
				{
					blockSeqLeft[blockId] = collectHitsForBlock(hits[blockId], totalHits, blockId*options.blockSize,
						swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options, windowNumbers[blockId]);
					
					// If any of the block has sequnece left the value will be true in the end.
					#pragma omp atomic
					sequenceLeft |= blockSeqLeft[blockId];
					
					waitingTimes[blockId] = sysTime();
				} // End task.
#endif
			} // End if.
		} // End for loop.
		
		// Wait for all blocks to stop bofore continuing.
		#pragma omp taskwait
		
#ifdef FLEX_TIMER
		std::cout << std::endl;
		// for waiting times
		_proFloat now = sysTime();
		for(int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId)
			std::cout << "filter>\t" << blockId << "\t" << (now - waitingTimes[blockId]) 
					<< "\t" << length(hits[blockId]) << "\t" <<  windowNumbers[blockId] << std::endl;
#endif
		
		return sequenceLeft;
	}


// TODO: doc
	template<
		typename TPosition,
		typename TSwiftHit,
		typename TRazerSOptions>
	inline void partitionHits(
			String<TPosition>						& positions,
			TPosition const							& startPos,
			String<TSwiftHit> const					& hits,
			TRazerSOptions const					& options)
	{ 
		// typedef String<ConvertedSwiftHit<TText> >			TConvertedHitString;
		// typedef typename Size<TConvertedHitString>::Type	TConvertedHitStringSize;
		
		resize(positions, options.numberOfBlocks+1, Exact());
		positions[0] = startPos;
		positions[options.numberOfBlocks] = length(hits); // Last element in positions
		
		TPosition partSize = (length(hits) - startPos) / options.numberOfBlocks;
		TPosition myPos;
		
		for (int blockId = 1; blockId < (int)options.numberOfBlocks; ++blockId){
			myPos = startPos + blockId*partSize;
			
			unsigned now = 0;
			unsigned last = hits[myPos - 1].ndlSeqNo;
			
			for(; myPos < length(hits); ++myPos){
				now = hits[myPos].ndlSeqNo;
				if(last != now)
					break;
				
				last = now;
			}
			positions[blockId] = myPos;
			
		}
#ifdef FLEX_TIMER
		// for(int i = 1; i < (int)length(positions); ++i)
		// 	printf("%lu, ", positions[i] - positions[i-1]);
		// printf("\n");
#endif
	}


	// TODO: doc
	template<
		typename TPosition,
		typename TSwiftHit,
		typename TVerifier,
		typename TReads,
		typename TContigSeq,
		typename TRazerSMode>
	inline void verifyHits(
			String<TPosition> const					& positions,
			TPosition const							& shortestPos,
			StringSet<String<TSwiftHit > > 			& hits,
			TVerifier								& verifier,
			TReads 									& readSet,
			TContigSeq 								& contigSeq,
			TRazerSMode const						& mode)
	{
		// typedef String<ConvertedSwiftHit<TText> >	THitString;
		
#ifdef FLEX_TIMER
				String<_proFloat> waitingTimes;
				resize(waitingTimes, (int)length(positions)-1, Exact());
				for(int coreId = 0; coreId < (int)length(positions)-1; ++coreId)
					waitingTimes[coreId] = sysTime();
				_proFloat myTime = sysTime();
#endif
		
		
		for(int i = 1; i < (int)length(positions); ++i){
			#pragma omp task default(none) \
				shared(positions, shortestPos, verifier, contigSeq, hits, readSet, mode, waitingTimes) \
				firstprivate(i)
			{
				// TVerifier & myVerifier = verifier[i-1];
				// THitString & myHits = hits[i-1];
				
				// First verify hits in your own string
				TPosition myPos = 0;
				printf("block: %d, from:%lu, to: %lu, length: %lu\n", i, myPos, shortestPos, length(hits[i-1]));
				for(; myPos < shortestPos; ++myPos){
					verifier[i-1].m.readId = hits[i-1][myPos].ndlSeqNo;
					matchVerify(verifier[i-1], swiftInfix(hits[i-1][myPos], contigSeq), hits[i-1][myPos].ndlSeqNo, readSet, mode);	
				}
				
#ifdef FLEX_TIMER
				waitingTimes[i-1] = sysTime();
#endif
			}
		}
		#pragma omp taskwait
		
#ifdef FLEX_TIMER
		std::cout << std::endl;
		// for waiting times
		_proFloat now = sysTime();
		for(int blockId = 0; blockId < (int)length(positions)-1; ++blockId)
			printf("verify 1>%d\t%f\n", blockId, (now - waitingTimes[blockId]));
			
		printf("verification 2 took:%f\n", (sysTime() - myTime));
		myTime = sysTime();
#endif
		
		for(int i = 1; i < (int)length(positions); ++i){
			#pragma omp task default(none) \
				shared(positions, shortestPos, verifier, contigSeq, hits, readSet, mode, waitingTimes) \
				firstprivate(i)
			{
				// TVerifier & myVerifier = verifier[i-1];
				// THitString & sharedHits = hits[0];
				
				// Second verifY your Partition in the shared string
				TPosition myPos = positions[i-1];
				printf("block: %d, from:%lu, to: %lu, length: %lu\n", i, myPos, positions[i], length(hits[i-1]));
				for(; myPos < positions[i]; ++myPos){
					verifier[i-1].m.readId = hits[0][myPos].ndlSeqNo;
					matchVerify(verifier[i-1], swiftInfix(hits[0][myPos], contigSeq), hits[0][myPos].ndlSeqNo, readSet, mode);	
				}
				
#ifdef FLEX_TIMER
				waitingTimes[i-1] = sysTime();
#endif
			}
		}
		#pragma omp taskwait
	
#ifdef FLEX_TIMER
		std::cout << std::endl;
		// for waiting times
		now = sysTime();
		for(int blockId = 0; blockId < (int)length(positions)-1; ++blockId)
			printf("verify 2>%d\t%f\n", blockId, (now - waitingTimes[blockId]));
		printf("verification 2 took:%f\n", (sysTime() - myTime));
#endif
	}


// TODO: doc
	template<typename TPos>
	struct Task
	{
		int blockId;
		TPos begin;
		TPos end;
	};

/**
.Function.goOverContigFlex:
..cat:Razers
..summary:
..signature:goOverContigFlex(swiftPatternHandler, swiftFinders, verifier, cnts, options, threadStores, store, mode)
..param.swiftPatternHandler:@Class.ParallelSwiftPatternHandler@
..param.swiftFinders:String of Finders
..param.verifier:String of Verifier
..param.cnts:For the Statistics
..param.options::RazerSOptions
..param.threadStores:String of @Class.FragmentStore@. One for each thread. Only temporarily used.
..param.store:Main @Class.FragmentStore@. Contains the matches in the end
..param.mode:RazerSMode
*/
	template <
		typename TContigSeq, 
		typename TReadIndex, 
		typename TSwiftSpec,
		typename TVerifier,
		typename TCounts,
		typename TRazerSOptions,
		typename TFragmentStore,
		typename TRazerSMode >
	void goOverContigFlex(
			ParallelSwiftPatternHandler<String<
				Pattern<TReadIndex, Swift<TSwiftSpec> > > >	& swiftPatternHandler,
			String<Finder<TContigSeq, Swift<TSwiftSpec> > >	& swiftFinders,
			TVerifier										& verifier,
			TCounts											& cnts,
			TRazerSOptions									& options,
			String<TFragmentStore>							& threadStores,
			TFragmentStore									& store,
			TRazerSMode const								& mode)
    {
		
		typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
		typedef typename TSwiftFinder::THitString				THitString;
		typedef typename Value<THitString>::Type				TSwiftHit;
		typedef typename Position<THitString>::Type				THitStringSize;
		typedef typename Iterator<THitString>::Type				THitStringIter;
		typedef typename Host<TSwiftFinder>::Type				TText;
		typedef typename TFragmentStore::TReadSeqStore			TReadSeqStore;
		typedef typename Infix<TContigSeq>::Type				TContigInfix;
		
		// initialize before the while loop so that the memory does not need to be reallocated every time
		StringSet<THitString, Owner<> > hits;
		resize(hits, options.numberOfBlocks, Exact());
		
		// To keep track if any block has some sequence left to cover.
		bool sequenceLeft = true;
		// To keep track which block is done and which not.
		String<bool> blockSeqLeft;
		fill(blockSeqLeft, options.numberOfBlocks, true, Exact());
		
		// Needed because the omp shared clause does not allow for "." in the variabel names.
		TReadSeqStore & readSet = store.readSeqStore;
	
		// TODO: is that needed?
		omp_set_nested(true);
				
		// Collect hits, verify them and store them in the main store while there is sequence left to cover.
		while(sequenceLeft){
				
			#pragma omp parallel num_threads((int)options.numberOfCores)
			{
			#pragma omp master
			{
				
#ifdef FLEX_TIMER
				_proFloat myTime = sysTime();
#endif
				
				// First: collect hits
				sequenceLeft = collectHits(swiftPatternHandler, swiftFinders, hits, blockSeqLeft, options);
				
#ifdef FLEX_TIMER
				printf("\n");
				printf("filtering took:%f\n", (sysTime() - myTime));
				myTime = sysTime();
#endif
				// for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId)
				// 	printf("block: %d, length %lu\n", blockId, length(hits[blockId]));
				
				// Second: Combine hits in one string so they can be sorted
				THitStringSize shortestLength = length(hits[0]);
				for (int blockId = 1; blockId < (int)options.numberOfBlocks; ++blockId)
					if(length(hits[blockId]) < shortestLength)
						shortestLength = length(hits[blockId]);
				
				for (int blockId = 1; blockId < (int)options.numberOfBlocks; ++blockId)
					append(hits[0], hits[blockId]);
				
				// for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId)
				// 	printf("block: %d, length %lu\n", blockId, length(hits[blockId]));
				
#ifdef FLEX_TIMER
				printf("appending took:%f\n", (sysTime() - myTime));
				myTime = sysTime();
#endif
				
				// TODO: remove
				// for(TConvertedHitStringSize myPos = 0; myPos < length(allHits); ++myPos){
				// 	printf("%u: (%ld, %ld)\n", allHits[myPos].readID, allHits[myPos].begin, allHits[myPos].end);
				// }
				
				// Third: Sort hits.
				// Get iterators for sorting
				// THitStringIter i1 = begin(hits[0]) + shortestLength;
				// THitStringIter i2 = end(hits[0]);
				
				// __gnu_parallel::sort(i1, i2, SwiftHitComparison<TSwiftHit>());
				// std::sort(i1, i2, SwiftHitComparison<TSwiftHit>());
				myRadixSort(hits[0]);
				
#ifdef FLEX_TIMER
				printf("sorting took:%f\n", (sysTime() - myTime));
				myTime = sysTime();
#endif
				
				// TODO: remove
				// for(TConvertedHitStringSize myPos = 0; myPos < length(allHits); ++myPos){
				// 	printf("%u: (%ld, %ld)\n", allHits[myPos].readID, allHits[myPos].begin, allHits[myPos].end);
				// }
				
				// Fourth: Get boundaries at which the hit string is split.
				String<THitStringSize> positions;
				partitionHits(positions, shortestLength, hits[0], options);
				
#ifdef FLEX_TIMER
				printf("partitioning took:%f\n", (sysTime() - myTime));
#endif
				
				// Fifth: Verify hits.
				verifyHits(positions, shortestLength, hits, verifier, readSet, host(swiftFinders[0]), mode);
				
				
				// clear hit strings
				for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId)
					clear(hits[blockId]);
				
				// fourth: get the matches from different thread stores and write them in the main one
				appendBlockStores(store, threadStores, swiftPatternHandler, cnts, options, mode);
				
				// clear thread stores
				for(int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
					clear(threadStores[blockId].alignedReadStore);
					clear(threadStores[blockId].alignQualityStore);
				}

			} // End omp master
			} // End omp parallel
						
		} // End while
		
		// clear finders
		for (unsigned int blockId = 0; blockId < options.numberOfBlocks; ++blockId)
			windowFindEnd(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId]);
		
	}

// END RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
#endif 
#if defined (RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS) || defined (RAZERS_PARALLEL_READS_INDEPENDENT)


/**
.Function._mapSingleReadsToContigFlex:
..cat:Razers
..summary:Appends the aligned read and quality stores from the block stores to the main store given as first argument
..signature:_mapSingleReadsToContigFlex(store, contigId, swiftPatterns, preprocessing, cnts, orientation, options, mode)
..param.store:@Class.FragmentStore@
..param.contigId: the ID of the contig within the fragment store to which the reads are mapped
..param.swiftPatternHandler: Handler for swift pattern
..param.preprocessing: String with bit vector patterns for each read
..param.cnts:Counts for statistics
..param.options:RazerSOptions
..param.mode:RazerSMode
*/
	template <
		typename TFragmentStore, 
		typename TReadIndex, 
		typename TSwiftSpec, 
		typename TPreprocessing,
		typename TCounts,
		typename TRazerSOptions,
		typename TRazerSMode >
    void _mapSingleReadsToContigFlex(
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
			_proFloat complementTime = sysTime();
			reverseComplementInPlace(contigSeq);
			printf("reverseComplement took %f sec\n", sysTime()-complementTime);
		} //TODO: reverseComplement in parallel?
		
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
		
		// A temporary store for every block. They are emptied after each verification step
		// String<TFragmentStore> threadStores;
		// resize(threadStores, options.numberOfBlocks, Exact());
		
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
		
		// TODO: neccessary
		#pragma omp parallel num_threads((int)options.numberOfCores)
		{
			#pragma omp for reduction(&:beginOk)
			for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
				beginOk = windowFindBegin(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], 
					options.errorRate);
			}	
		}
		// Only if the finders are set up.
		if(beginOk){
			#ifdef RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
				goOverContigFlex(swiftPatternHandler, swiftFinders, verifier, cnts, options, threadStores, store, mode);
			#endif
			#ifdef RAZERS_PARALLEL_READS_INDEPENDENT
				goOverContigIndependent2(swiftPatternHandler, swiftFinders, verifier, cnts, options, threadStores, store, mode);
			#endif
		}
		
		if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
		if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation
		
	}

// END RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS or RAZERS_PARALLEL_READS_INDEPENDENT
#endif
#ifdef RAZERS_PARALLEL_READS_WAIT_AFTER_WINDOW 
	
	template <
	typename TContigSeq, 
	typename TReadIndex, 
	typename TSwiftSpec,
	typename TVerifier,
	typename TRazerSOptions,
#ifdef RAZERS_TIMER
	typename TFragmentStore,
#endif
	typename TRazerSMode >
    void _goOverContig(
			TContigSeq										& contigSeq,
			ParallelSwiftPatternHandler<String<
			Pattern<TReadIndex, Swift<TSwiftSpec> > > >		& swiftPatternHandler,
			String<Finder<TContigSeq, Swift<TSwiftSpec> > >	& swiftFinders,
			TVerifier										& verifier,
			TRazerSOptions									& options,
#ifdef RAZERS_TIMER
			TFragmentStore									& store,
			unsigned										  contigId,
			char											  orientation,
#endif
			TRazerSMode const                               & mode)
    {
		
		typedef Finder<TContigSeq, Swift<TSwiftSpec> >	TSwiftFinder;
		typedef typename TSwiftFinder::THitString		THitString;
		typedef typename Value<THitString>::Type		TSwiftHit;
		typedef typename Size<THitString>::Type			THitStringSize;
		
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
		// go over contig sequence
		#pragma omp parallel num_threads((int)options.numberOfCores)
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
				stop = !windowFindNext(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options.windowSize, posLength);
				
				// get time for filtering
				_proFloat filteringTime = sysTime() - startTime;
				// start time for verification
				startTime = sysTime();
				// number of matches
				int matches = 0;
				
				// verify hits
				THitString hits = getSwiftHits(swiftFinders[blockId]);
				for(THitStringSize h = 0; h < length(hits); ++h){
					verifier[blockId].m.readId = (blockId * blockSize) + hits[h].ndlSeqNo;         //array oder jedesmal berechnen
					matches += matchVerify(verifier[blockId], getSwiftRange(hits[h], contigSeq), hits[h].ndlSeqNo, host(host(swiftPatternHandler.swiftPatterns[blockId])), mode);
					++options.countFiltration;
				}

				// get time for filtering
				_proFloat verificationTime = sysTime() - startTime;
				
				#pragma omp critical
				{
					std::cout << "timing>\t";
					#ifdef _OPENMP
					std::cout << omp_get_thread_num() << "\t";
					#endif
					std::cout << blockId << "\t";
					std::cout << contigAndDirection << "\t";
					std::cout << posLength.i1 << "\t";
					std::cout << posLength.i2 << "\t";
					std::cout << filteringTime << "\t";
					std::cout << length(hits) << "\t";
					std::cout << verificationTime << "\t";
					std::cout << matches << "\n";
				}
				// set waiting time
				#ifdef _OPENMP
				waitingTimes[omp_get_thread_num()] = sysTime();
				#endif
// not RAZERS_TIMER
#else 
				stop = !windowFindNext(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options.windowSize);
				
				// verify hits
				THitString & hits = getSwiftHits(swiftFinders[blockId]);
				for(THitStringSize h = 0; h < length(hits); ++h){
					verifier[blockId].m.readId = (blockId * blockSize) + hits[h].ndlSeqNo;         //array oder jedesmal berechnen
					matchVerify(verifier[blockId], swiftInfix(hits[h], contigSeq), hits[h].ndlSeqNo, host(host(swiftPatternHandler.swiftPatterns[blockId])), mode);
					++options.countFiltration;
				}
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
		typedef typename Value<TPreprocessing>::Type                                TPreprocessingBlock;
		typedef String<Pattern<TReadIndex, Swift<TSwiftSpec> > >                    TSwiftPatterns;
		typedef ParallelSwiftPatternHandler<TSwiftPatterns>                         TSwiftPatternHandler;
		typedef MatchVerifier <TFragmentStore, TRazerSOptions,
		TRazerSMode, TPreprocessingBlock,
		TSwiftPatternHandler, TCounts >									TVerifier;
		typedef typename Fibre<TReadIndex, Fibre_Text>::Type                        TReadSet;
		typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type		TAlignedReadStoreSize;

		if (options._debugLevel >= 1)
		{
			::std::cerr << ::std::endl << "Process genome seq #" << contigId;
			if (orientation == 'F') ::std::cerr << "[fwd]" << std::endl;
			else                    ::std::cerr << "[rev]" << std::endl;
		}

		lockContig(store, contigId);
		TContigSeq &contigSeq = store.contigStore[contigId].seq;
		if (orientation == 'R'){
			_proFloat complementTime = sysTime();
			reverseComplementInPlace(contigSeq);
			printf("reverseComplement took %f sec\n", sysTime()-complementTime);
		}
			

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
#ifndef RAZERS_TIMER
			_goOverContig(contigSeq, swiftPatternHandler, swiftFinders, verifier, options, mode);
#else
			_goOverContig(contigSeq, swiftPatternHandler, swiftFinders, verifier, options, store, contigId, orientation, mode);		
#endif
			
			// append the alignedReadStore and the alignQualityStore from the blocks to the main store
			appendBlockStores(store, blockStores, swiftPatternHandler, cnts, options, mode);
        }
		
		if (!unlockAndFreeContig(store, contigId))							// if the contig is still used
			if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation
	}

// End: RAZERS_PARALLEL_READS_WAIT_AFTER_WINDOW
#endif
	
/**
.Function._mapSingleReadsParallelCreatePatterns:
..cat:Razers
..summary:Creates @Class.Pattern@s for the filtration (@Class.SwiftLocal@) and verification (@Class.Myers@)
..signature:_mapSingleReadsParallelCreatePatterns(store, cnts, options, mode, readIndices)
..param.store:@Class.FragmentStore@
..param.cnts:Counts for statistics
..param.options:RazerSOptions
..param.mode:RazerSMode
..param.readIndices:@Class.String@ over @Class.Index_QGram@
*/
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

		typedef FragmentStore<TFSSpec, TFSConfig>					TFragmentStore;
		typedef typename IF<TYPECMP<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
		typedef typename Value<TReadIndexString>::Type				TReadIndex;

		// filter
		typedef Pattern<TReadIndex, Swift<TSwiftSpec> >				TSwiftPattern;
		typedef String<TSwiftPattern>								TSwiftPatterns;
		typedef ParallelSwiftPatternHandler<TSwiftPatterns>			TSwiftPatternHandler;

		// verifier
		typedef Pattern<TRead, Myers<FindInfix, False, void> >		TMyersPattern;
		typedef typename Infix<String<TMyersPattern> >::Type		TVerifierBlock;
		typedef typename Position<TReadIndexString>::Type			TPos;

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
		
#ifdef RAZERS_PARALLEL_READS_WAIT_AFTER_WINDOW
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
#endif
#ifdef RAZERS_PARALLEL_READS_INDEPENDENT
		// A temporary store for every block. They are emptied after each verification step
		String<TFragmentStore> threadStores;
		resize(threadStores, options.numberOfBlocks, Exact());
		
		String<RazerSOptions<TSpec> > threadOptions;
		fill(threadOptions, options.numberOfBlocks, options, Exact());
#endif

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
		
#ifdef RAZERS_TIMER
		// print header line for timer
		std::cout << "timing>\tthread\tblock\tcontigAndDircetion\tpos\tlength\tfilter.time\tverifications\tverification.time\tmatches\n";
		std::cout << "waiting>\tthread\ttime\n";
#endif

		// iterate over genome sequences
		for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId){
			// lock to prevent releasing and loading the same contig twice
			// (once per _mapSingleReadsToContig call)
			lockContig(store, contigId);
			if (options.forward){
				#if defined (RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS) || defined (RAZERS_PARALLEL_READS_INDEPENDENT)
				_mapSingleReadsToContigFlex(store, threadStores, contigId, swiftPatternHandler, forwardPatterns, cnts, 'F', options, threadOptions, mode);
				#elif RAZERS_PARALLEL_READS_WAIT_AFTER_WINDOW
				_mapSingleReadsToContig(store, contigId, swiftPatternHandler, forwardPatternsBlocks, cnts, 'F', options, mode);
				#endif
			}
			if (options.reverse){
				#if defined (RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS) || defined (RAZERS_PARALLEL_READS_INDEPENDENT)
				_mapSingleReadsToContigFlex(store, threadStores, contigId, swiftPatternHandler, forwardPatterns, cnts, 'R', options, threadOptions, mode);
				#elif RAZERS_PARALLEL_READS_WAIT_AFTER_WINDOW
				_mapSingleReadsToContig(store, contigId, swiftPatternHandler, forwardPatternsBlocks, cnts, 'R', options, mode);
				#endif
			}
			unlockAndFreeContig(store, contigId);
		}
		
#ifdef RAZERS_PARALLEL_READS_INDEPENDENT
		appendBlockStores(store, threadStores, options);
		for(int i = 0; i < (int) options.numberOfBlocks; ++i){
			// options.countFiltration += threadOptions[i].countFiltration;
			options.countVerification += threadOptions[i].countVerification;
		}
#endif		

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
}

#endif