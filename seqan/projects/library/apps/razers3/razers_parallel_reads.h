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
#ifdef _OPENMP
#include <omp.h>
#endif

#define RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS

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
		SEQAN_CHECKPOINT

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
			
			if(store.alignedReadStore[i].readId == 36){
				#pragma omp critical(output1)
				std::cout << "inStore 37 > " << store.alignedReadStore[i].beginPos << ", " << store.alignedReadStore[i].endPos << std::endl;
			}
//			if(store.alignedReadStore[i].beginPos == 28927578){
//				TAlignedReadStoreElem a = store.alignedReadStore[i];
//				TAlignedQualStoreElem q = store.alignQualityStore[a.id];
//				int k = 0; k = 1;
//			}
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
		
		//delete:
//		for(unsigned i = 0; i < length(store.alignedReadStore); ++i){
//			if(store.alignedReadStore[i].readId == 36){
//				#pragma omp critical(output1)
//				std::cout << "inStore 37 > " << store.alignedReadStore[i].beginPos << ", " << store.alignedReadStore[i].endPos << " (not compacted)" << std::endl;
//			}
//		}
		
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

		//delete:
//		for(unsigned i = 0; i < length(store.alignedReadStore); ++i){
//			if(store.alignedReadStore[i].readId == 36){
//				#pragma omp critical(output1)
//				std::cout << "inStore 37 > " << store.alignedReadStore[i].beginPos << ", " << store.alignedReadStore[i].endPos << " (compacted)" << std::endl;
//			}
//		}
		
	}
	
#ifdef RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
	
/**
.Class._ConvertedSwiftHit:
..summary:Representation of a Swift hit containing the corresponding read ID and infix in the reference.
..cat:Razers
..signature:_ConvertedSwiftHit<TText>
..param.TText:Text on which the @Spec.Infix@ is based on.
*/
	template <typename TText>
	struct _ConvertedSwiftHit
	{
		typedef typename Infix<TText>::Type TInfix;

		static const unsigned accuracy = 100;

		unsigned	readID;
		TInfix		onContig;

		_ConvertedSwiftHit(){}
		
		// TODO: better template
		template <typename TSpec, typename THstkPos, typename TText2>
		_ConvertedSwiftHit(_SwiftHit<Tag<_SwiftSemiGlobal<TSpec> >, THstkPos> & other, TText2 & ref):
			readID(other.ndlSeqNo),
			onContig(getSwiftRange(other, ref))
		{}
		/*
		bool operator<(_ConvertedSwiftHit const & other){
			return (readID / accuracy) < (other.readID / accuracy);
		}
		*/
	};
	
	
	struct _ConvertedSwiftHitComparison {
		
		unsigned accuracy;
		
		_ConvertedSwiftHitComparison(unsigned _accuracy):
			accuracy(_accuracy)
		{}
		
		template<typename TText>
		bool operator() (_ConvertedSwiftHit<TText> i1, _ConvertedSwiftHit<TText> i2) {
			return (i1.readID / accuracy) < (i2.readID / accuracy);
		}
	};

/**
.Function._collectHits:
..cat:Razers
..summary:Searches for new hits with the finder and pattern.
Stops when the finder reaches its end or the threshold of total hits is surpassed.
..signature:_collectHits(hits, totalHits, threshold, finder, pattern, windowSize)
..param.hits:String of @Class._ConvertedSwiftHit@ to store the found hits
..param.totalHits:Variable shared by all threads to know when to stop.
..param.threshold:The function stops searching for new hits when totalHits surpasses this value.
..param.finder:Swift finder.
..param.pattern:Swift pattern.
..param.windowSize:Window size @Function.windowFindNext@ is called with.
*/
	template <
		typename TConvertedSwiftHitString,
		typename TSwiftFinder,
		typename TSwiftPattern>
	bool _collectHits(
			TConvertedSwiftHitString	& hits,
			unsigned					& totalHits,
			unsigned const				  threshold,
			TSwiftFinder				& finder,
			TSwiftPattern				& pattern,
			unsigned const				  windowSize)
	{
		typedef typename TSwiftFinder::THitString				TSwiftHitString;
		typedef typename Size<TConvertedSwiftHitString>::Type	TSize; 
		typedef typename Host<TSwiftFinder>::Type				TText;
		
		// Get new hits till all threads together have enough.
		bool sequenceLeft = true;
		while(totalHits < threshold and sequenceLeft){
			// Search for new hits.
			sequenceLeft = windowFindNext(finder, pattern, windowSize);
			TSwiftHitString oriHits = getSwiftHits(finder);
			// Resize the hit string so it can hold the new hits
			TSize oldSize = length(hits);
			resize(hits, (oldSize + length(oriHits)));
			// Convert hits and append them to the hit string for this block.
			// Conversion is neccessary so that the hits can be sorted by the readID.
			for(TSize hitID = 0; hitID < length(oriHits); ++hitID)
				hits[oldSize + hitID] = _ConvertedSwiftHit<TText>(oriHits[hitID], host(finder));
			
			// Update the overall hit count so that the other threads know when to stop as well.
			#pragma omp atomic
			totalHits += length(oriHits);
		}
		return sequenceLeft;
	}

/**
.Class._AssignmentDetail:
..summary:Representation of a verification task.
..cat:Razers
..signature:_AssignmentDetail<TPos>
..param.TPos: Type of the start and end positions
*/
	template <typename TPos>
	struct _TaskDetails{
		int blockId;
		TPos current;
		TPos end;
		bool running;
		bool split;
		int splitWith;
		bool sorted;
		
		_TaskDetails(){}
		
		_TaskDetails(int _blockId, TPos _end):
			blockId(_blockId),
			current(0),
			end(_end),
			running(true),
			split(false),
			splitWith(-1),
			sorted(false)
		{}
		
		_TaskDetails(int _blockId, TPos _start, TPos _end):
			blockId(_blockId),
			current(_start),
			end(_end),
			running(true),
			split(false),
			splitWith(-1),
			sorted(true)
		{}
		
	};

	// TODO: doku
	template <
		typename TText,
		typename TSize,
		typename TReadSet,
		typename TVerifier,
		typename TRazerSOptions,
		typename TRazerSMode>
	bool _verifyHits(
					 StringSet<String<_ConvertedSwiftHit<TText> > >	& hits,
					 int const										  myTaskId,
					 String<_TaskDetails<TSize> >					& tasks,
					 TReadSet										& readSet,
					 String<TVerifier> 								& verifier,
					 TRazerSOptions									& options,
					 TRazerSMode const								& mode,
					 int											  level) //TODO: remove
	{
		typedef String<_ConvertedSwiftHit<TText> >		THitString;
		typedef typename Position<THitString>::Type		THitStringPos;
		typedef typename Iterator<THitString>::Type		THitIter;
		typedef _TaskDetails<TSize>						TTask;
				
		TTask & task = tasks[myTaskId];
		TSize & current = task.current;
		THitString & myHits = hits[task.blockId];
		
		// go over all hits or till a different thread set a flag to stop and split
		for(; (current < task.end) and (not task.split); ++current){
			unsigned myId = (task.blockId * options.blockSize) + myHits[current].readID;
			int threadId = myTaskId; //omp_get_thread_num();
			verifier[threadId].m.readId = myId;
			
			matchVerify(verifier[threadId], myHits[current].onContig, myId, readSet, mode);
			
			#pragma omp atomic
			++options.countFiltration;
			
			#pragma omp flush(tasks)
		}
		
		// if split is true
		if(task.split){
			// set running false
			task.running = false;
			
			// TODO:replace with option
			unsigned accuracy = 5;
			
			// If not sorted yet, sort remaining hits.
			if(not task.sorted){
				// get iterators for sorting
				THitIter i1 = begin(myHits) + current;
				THitIter i2 = end(myHits);
				
				// TODO replace with parallel sort
				_ConvertedSwiftHitComparison comp(accuracy);
				sort(i1, i2, comp);
			}
			
			// search for first hit after the half for that the read Id changes significantly
			// split the hits at this position, so that the myers patterns shared by the 
			// verifier don't collide
			TSize left = task.end - task.current;
			if (left == 0) return true;
			TSize half = left / 2;
			unsigned now, last = myHits[task.current + half - 1].readID / accuracy;
			TSize i = task.current + half;
			for(; i < task.end; ++i){
				now = myHits[i].readID / accuracy;
				if(last != now)
					break;
				
				last = now;
			}
			
			// split and start new tasks
			int splitWith = task.splitWith;
			
			// First the other one so the task.end can be used without 
			// writing it in a temporary variable.

			_TaskDetails<TSize> me(task.blockId, task.current, i);
			_TaskDetails<TSize> you(task.blockId, i, task.end);
			
			#pragma omp flush(tasks)
			tasks[splitWith] = you;
			tasks[myTaskId]  = me;
			#pragma omp flush(tasks)
			
			#pragma omp parallel sections
			{
				#pragma omp section
				_verifyHits(hits, splitWith, tasks, readSet, verifier, options, mode, level+1);
				
				#pragma omp section
				_verifyHits(hits, myTaskId, tasks, readSet, verifier, options, mode, level+1);
			}
						
		}
		else // if ended by itself
		{
			
			// critical so two threads don't make the same one split
			#pragma omp critical(split_with)
			{
				// Inital value of maxHitsLeft is also the minimum number of hits needed
				// to trigger a the splitting.
				THitStringPos maxHitsLeft = 10;
				bool oneIsSplitting = false;
				
				int splitWith = -1;
				
				// check if another block has sufficient many hits left to steal work and non is in the process of splitting
				for(int taskId = 0; taskId < (int)options.numberOfBlocks; ++taskId){
					//
					if(not tasks[taskId].split){
						THitStringPos hitsLeft = tasks[taskId].end;
						if(hitsLeft > maxHitsLeft){
							maxHitsLeft = hitsLeft;
							splitWith = taskId;
						}
					}
					else
						oneIsSplitting = true;	
				}
				// set running false
				task.running = false;
			}
			
		}
		return true;
	}


	// TODO: doku
	template <
		typename TContigSeq, 
		typename TReadIndex, 
		typename TSwiftSpec,
		typename TVerifier,
		typename TCounts,
		typename TRazerSOptions,
		typename TFragmentStore,
		typename TRazerSMode >
	void _goOverContigFlex(
			//TContigSeq										& contigSeq,
			ParallelSwiftPatternHandler<String<
				Pattern<TReadIndex, Swift<TSwiftSpec> > > > & swiftPatternHandler,
			String<Finder<TContigSeq, Swift<TSwiftSpec> > >	& swiftFinders,
			TVerifier										& verifier,
			TCounts											& cnts,
			TRazerSOptions									& options,
			String<TFragmentStore>							& threadStores,
			TFragmentStore									& store,
			TRazerSMode const								& mode)
    {
		
		typedef Finder<TContigSeq, Swift<TSwiftSpec> >  TSwiftFinder;
		typedef typename TSwiftFinder::THitString		THitString;
		typedef typename Value<THitString>::Type		TSwiftHit;
		typedef typename Position<THitString>::Type		THitStringSize;
		typedef typename Host<TSwiftFinder>::Type		TText;
		typedef String<_ConvertedSwiftHit<TText> >		TConvertedHitString;
		typedef typename Size<TConvertedHitString>::Type	TConvertedHitStringSize;
		typedef typename TFragmentStore::TReadSeqStore	TReadSeqStore;
		
			// initialize before the while loop so that the memory does not need to be reallocated every time
			StringSet<TConvertedHitString, Owner<> > hits;
			resize(hits, options.numberOfBlocks, Exact());
			
			// Number of hits that are colleced before the verification is started
			unsigned threshold = 10000;
			// To keep track if any block has some sequence left to cover.
			bool sequenceLeft = true;
			// To keep track which block is done and which not.
			String<bool> blockSeqLeft;
			fill(blockSeqLeft, options.numberOfBlocks, true, Exact());
			
			// ??? is that needed?
			omp_set_nested(true);
			
			// Collect hits, verify them and store them in the main store while there is sequence left to cover.
			while(sequenceLeft)
			{
				// first: collect hits
				// Number of hits collected in this iteration. Shared by all threads
				// Updated in an atomic expression
				unsigned totalHits = 0;
				
				sequenceLeft = false;
				
				// Go over sequence with blocks in parallel
				#pragma omp parallel num_threads((int)options.numberOfCores)
				{
				#pragma omp for
				for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
					// Only use block if there is sequnce left.
					// Otherwise findwindow next does unpredicted things.
					if(blockSeqLeft[blockId]){
						blockSeqLeft[blockId] = _collectHits(hits[blockId], totalHits, threshold,
															 swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options.windowSize);
						// if any of the block has sequnece left the while loops again
						sequenceLeft |= blockSeqLeft[blockId];
					}
				}
				}
					
				/* Alternative to preceding parallel for
				 if used needs to be adjusted as the for loop is
				 #pragma omp parallel num_threads((int)options.numberOfCores)
				 {
				 int blockId = omp_get_thread_num();
				 bool left = _collectHits(hits[blockId], totalHits, threshold, swiftFinders[blockId], 
				 swiftPatternHandler.swiftPatterns[blockId], options.windowSize);
				 
				 #pragma omp atomic
				 sequenceLeft |= left;
				 }
				 */
					
				// third: verify the hits
				// create verification tasks
				String<_TaskDetails<TConvertedHitStringSize> > tasks;
				resize(tasks, options.numberOfBlocks, Exact());
				for (int taskId = 0; taskId < (int)options.numberOfBlocks; ++taskId)
					tasks[taskId] = _TaskDetails<TConvertedHitStringSize>(taskId, length(hits[taskId]));
									
				// start verification
				#pragma omp parallel num_threads((int)options.numberOfCores)
				{
				#pragma omp for
				for (int taskId = 0; taskId < (int)options.numberOfBlocks; ++taskId){
					_verifyHits(hits, taskId, tasks, store.readSeqStore, verifier, options, mode, 0);
				}					
				}
					
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
			}

		// clear finders
		for (unsigned int blockId = 0; blockId < options.numberOfBlocks; ++blockId)
			windowFindEnd(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId]);
	}

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
		int													  contigId,
		ParallelSwiftPatternHandler<String<
			Pattern<TReadIndex, Swift<TSwiftSpec> > > >		& swiftPatternHandler,
		TPreprocessing										& preprocessing,
		TCounts												& cnts,
		char												  orientation,
		TRazerSOptions										& options,
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
		if (orientation == 'R')	reverseComplementInPlace(contigSeq);

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
		resize(verifier, options.numberOfCores, Exact());
		
		// A temporary store for every block. They are emptied after each verification step
		String<TFragmentStore> threadStores;
		resize(threadStores, options.numberOfCores , Exact());
		for(int threadId = 0; threadId < (int)options.numberOfCores; ++threadId){
			// initialize verifier
			TVerifier oneVerifier(threadStores[threadId], options, preprocessing, swiftPatternHandler, cnts);
			oneVerifier.onReverseComplement = (orientation == 'R');
			oneVerifier.genomeLength = length(contigSeq);
			oneVerifier.m.contigId = contigId;
			// assign it to the string
			verifier[threadId] = oneVerifier;
		}

		// Set up finder. beginOK is true after the loop if all finders are set up successfully.
		bool beginOk = true;
		for (int blockId = 0; blockId < (int)options.numberOfBlocks; ++blockId){
			beginOk = beginOk & windowFindBegin(swiftFinders[blockId], swiftPatternHandler.swiftPatterns[blockId], options.errorRate);
			if(not beginOk) break;
		}

		// Only if the finders are set up.
		if(beginOk)
			_goOverContigFlex(swiftPatternHandler, swiftFinders, verifier,
								cnts, options, threadStores, store, mode);
		
		if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
		if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation

		//delete:
		//std::cout << "genomeLength: " << verifier[0].genomeLength << std::endl;
		
	}

#else // not RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
	
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
			Pattern<TReadIndex, Swift<TSwiftSpec> > > > & swiftPatternHandler,
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
				THitString hits = getSwiftHits(swiftFinders[blockId]);
				for(THitStringSize h = 0; h < length(hits); ++h){
					verifier[blockId].m.readId = (blockId * blockSize) + hits[h].ndlSeqNo;         //array oder jedesmal berechnen
					matchVerify(verifier[blockId], getSwiftRange(hits[h], contigSeq), hits[h].ndlSeqNo, host(host(swiftPatternHandler.swiftPatterns[blockId])), mode);
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
#endif // end else: not RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
	
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
		
#ifndef RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
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
				#ifdef RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
				_mapSingleReadsToContigFlex(store, contigId, swiftPatternHandler, forwardPatterns, cnts, 'F', options, mode);
				#else
				_mapSingleReadsToContig(store, contigId, swiftPatternHandler, forwardPatternsBlocks, cnts, 'F', options, mode);
				#endif
			}
			if (options.reverse){
				#ifdef RAZERS_PARALLEL_READS_FLEXIBLE_VERIFICATION_BLOCKS
				_mapSingleReadsToContigFlex(store, contigId, swiftPatternHandler, forwardPatterns, cnts, 'R', options, mode);
				#else
				_mapSingleReadsToContig(store, contigId, swiftPatternHandler, forwardPatternsBlocks, cnts, 'R', options, mode);
				#endif
			}
			unlockAndFreeContig(store, contigId);
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