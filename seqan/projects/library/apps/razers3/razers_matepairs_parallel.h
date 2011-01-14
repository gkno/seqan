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

#include "razers_parallel.h"
#include "razers_matepairs.h"

namespace SEQAN_NAMESPACE_MAIN
{

// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
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
	FragmentStore<TFSSpec, TFSConfig>		& store,
	unsigned								  contigId,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> >	& swiftPatternL,
	Pattern<TReadIndex, Swift<TSwiftSpec> >	& swiftPatternR,
#ifdef RAZERS_BANDED_MYERS
	TPreprocessing							&,
	TPreprocessing							&,
#else  // #ifdef RAZERS_BANDED_MYERS
	TPreprocessing							& preprocessingL,
	TPreprocessing							& preprocessingR,
#endif  // #ifdef RAZERS_BANDED_MYERS
	TCounts									& cnts,
	char									  orientation,			// q-gram index of reads
	TRazerSOptions							& options,
	TRazerSMode						  const & mode)
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
	typedef MatchVerifier <
		TFragmentStore, 
		TRazerSOptions, 
		TRazerSMode, 
		TSwiftPattern,
		TCounts,
        TPreprocessing>											TVerifier;

	const unsigned NOT_VERIFIED = 1u << (8*sizeof(unsigned)-1);

	// iterate all genomic sequences
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

	TReadSet	&readSetL = host(host(swiftPatternL));
	TReadSet	&readSetR = host(host(swiftPatternR));
	TVerifier	verifierL(store, options, swiftPatternL, cnts);
	TVerifier	verifierR(store, options, swiftPatternR, cnts);

#ifndef RAZERS_BANDED_MYERS
	verifierL.preprocessing = &preprocessingL;
	verifierR.preprocessing = &preprocessingR;
#endif  // #ifdef RAZERS_BANDED_MYERS

	verifierL.oneMatchPerBucket = true;
	verifierR.oneMatchPerBucket = true;
	verifierL.m.contigId = contigId;
	verifierR.m.contigId = contigId;

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
		TGPos rEndPos = endPosition(swiftFinderR) + scanShift;
		TGPos doubleParWidth = 2 * (*swiftFinderR.curHit).bucketWidth;
		
		// remove out-of-window left mates from fifo
		while (!empty(fifo) && (TSignedGPos)front(fifo).i2.endPos + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
		{
			popFront(fifo);
			++firstNo;
		}
/*		
		if (empty(fifo) || back(fifo).endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
			for (unsigned i = 0; i < preFetchMatches; ++i)
				if (find(swiftFinderL, swiftPatternL, options.errorRate, false))
					pushBack(fifo, mL);
				else
					break;
*/
		// add within-window left mates to fifo
		while (empty(fifo) || (TSignedGPos)back(fifo).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
		{
			if (find(swiftFinderL, swiftPatternL, options.errorRate))
			{
				gPair = positionRange(swiftFinderL);
				if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
				{
					// link in
					fL.i1 = lastPotMatchNo[swiftPatternL.curSeqNo];
					lastPotMatchNo[swiftPatternL.curSeqNo] = lastNo++;
					
					fL.i2.readId = store.matePairStore[swiftPatternL.curSeqNo].readId[0] | NOT_VERIFIED;
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
		unsigned leftReadId = store.matePairStore[matePairId].readId[0];
		__int64 lastPositive = (__int64)-1;
		for (__int64 i = lastPotMatchNo[matePairId]; firstNo <= i; i = (*it).i1)
		{
			it = &value(fifo, i - firstNo);

			// search left mate
//			if (((*it).i2.readId & ~NOT_VERIFIED) == leftReadId)
			{
				// verify left mate (equal seqNo), if not done already
				if ((*it).i2.readId & NOT_VERIFIED)
				{
//					if (matchVerify(
//							(*it).i2, (*it).i3, infix(genome, (TSignedGPos)(*it).i2.beginPos, (TSignedGPos)(*it).i2.endPos), 
//							matePairId, readSetL, forwardPatternsL, 
//							options, TSwiftSpec()))
					if ((TSignedGPos)(*it).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
					{
						if (matchVerify(verifierL, infix(genome, (TSignedGPos)(*it).i2.beginPos, (TSignedGPos)(*it).i2.endPos), 
								matePairId, readSetL, mode))
						{
							verifierL.m.readId = (*it).i2.readId & ~NOT_VERIFIED;		// has been verified positively
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
/*
				if ((*it).i2.readId == leftReadId)
				{
					bestLeft = it;
					bestLeftScore = (*it).i3.score;
					break;
				}
*/
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
//			if (matchVerify(
//					mR, qR, infix(swiftFinderR),
//					matePairId, readSetR, forwardPatternsR,
//					options, TSwiftSpec()))
			if (matchVerify(verifierR, infix(swiftFinderR), 
					matePairId, readSetR, mode))
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
					TMatePair &mp = store.matePairStore[matePairId];
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
					if (++options.nextPairMatchId == TAlignedRead::INVALID_ID)
						options.nextPairMatchId = 0;

					// score the whole match pair
					fL.i3.pairScore = qR.pairScore = fL.i3.score + qR.score;

					// both mates match with correct library size
/*								std::cout << "found " << matePairId << " on " << orientation << contigId;
					std::cout << " dist:" << dist;
					if (orientation=='F')
						std::cout << " \t_" << fL.i2.beginPos+1 << "_" << mR.endPos;
					else
						std::cout << " \t_" << mR.beginPos+1 << "_" << mL.endPos;
//							std::cout << " L_" << (*bestLeft).beginPos << "_" << (*bestLeft).endPos << "_" << (*bestLeft).editDist;
//							std::cout << " R_" << mR.beginPos << "_" << mR.endPos << "_" << mR.editDist;
					std::cout << std::endl;
*/
					if (!options.spec.DONT_DUMP_RESULTS)
					{
						fL.i2.id = length(store.alignedReadStore);
						appendValue(store.alignedReadStore, fL.i2, Generous());
						appendValue(store.alignQualityStore, fL.i3, Generous());
						mR.id = length(store.alignedReadStore);
						appendValue(store.alignedReadStore, mR, Generous());
						appendValue(store.alignQualityStore, qR, Generous());

						if (length(store.alignedReadStore) > options.compactThresh)
						{
							typename Size<TAlignedReadStore>::Type oldSize = length(store.alignedReadStore);
//									maskDuplicates(matches);	// overlapping parallelograms cause duplicates
							compactPairMatches(store, cnts, options, swiftPatternL, swiftPatternR);
							
							if (length(store.alignedReadStore) * 4 > oldSize)			// the threshold should not be raised
								options.compactThresh += (options.compactThresh >> 1);	// if too many matches were removed
							
							if (options._debugLevel >= 2)
								::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
						}
					}
					++options.countVerification;
				}
				++options.countFiltration;
			}
		}
	}
	
	if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
		if (orientation == 'R')	reverseComplement(genome);	// we have to restore original orientation

}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
	
template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
    typename TMatchNPolicy>
int _mapMatePairReadsParallel(
	FragmentStore<TFSSpec, TFSConfig>					& store,
	TCounts												& cnts,
	RazerSOptions<TSpec>								& options,
	TShape const										& shape,
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy>  const & mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	
	typedef typename Value<TReadSeqStore>::Type			TRead;
	typedef StringSet<TRead>							TReadSet;
#ifndef RAZERS_OPENADDRESSING
	typedef Index<TReadSet, IndexQGram<TShape> >	TIndex;			// q-gram index
#else
	typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >	TIndex;
#endif

	typedef typename If<
				IsSameType<TGapMode,RazerSGapped>::VALUE,
				SwiftSemiGlobal,
				SwiftSemiGlobalHamming>::Type			TSwiftSpec;
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

    // Save OpenMP maximal thread count so we can restore it below, then set
    // from options.
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Number of threads:               \t" << options.threadCount << std::endl;
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount);

//	std::cout << "SA-TYPE:" <<sizeof(typename SAValue<TIndex>::Type)<<std::endl;

	// split mate-pairs over two indices
	TReadSet readSetL, readSetR;
	unsigned pairCount = length(store.matePairStore);
	resize(readSetL, pairCount, Exact());
	resize(readSetR, pairCount, Exact());

	for (unsigned i = 0; i < pairCount; ++i)
	{
		assign(readSetL[i], store.readSeqStore[store.matePairStore[i].readId[0]]);
		assign(readSetR[i], store.readSeqStore[store.matePairStore[i].readId[1]]);
	}
	reverseComplement(readSetR);

	// configure q-gram index
	TIndex swiftIndexL(readSetL, shape);
	TIndex swiftIndexR(readSetR, shape);
#ifdef RAZERS_OPENADDRESSING
	swiftIndexL.alpha = options.loadFactor;
	swiftIndexR.alpha = options.loadFactor;
#endif
	reverse(indexShape(swiftIndexR));		// right mate qualities are reversed -> reverse right shape
	
	cargo(swiftIndexL).abundanceCut = options.abundanceCut;
	cargo(swiftIndexR).abundanceCut = options.abundanceCut;
	cargo(swiftIndexL)._debugLevel = options._debugLevel;
	cargo(swiftIndexR)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPatternL(swiftIndexL);
	TSwiftPattern swiftPatternR(swiftIndexR);
	swiftPatternL.params.minThreshold = options.threshold;
	swiftPatternR.params.minThreshold = options.threshold;
	swiftPatternL.params.tabooLength = options.tabooLength;
	swiftPatternR.params.tabooLength = options.tabooLength;
	swiftPatternL.params.printDots = false; // only one should print the dots
	swiftPatternR.params.printDots = options._debugLevel > 0;

#ifdef RAZERS_BANDED_MYERS
	typedef Nothing TPreprocessing;
	TPreprocessing forwardPatternsL;
	TPreprocessing forwardPatternsR;
#else  // #ifdef RAZERS_BANDED_MYERS
	// init edit distance verifiers
    typedef String<TMyersPattern> TPreprocessing;
	TPreprocessing forwardPatternsL;
	TPreprocessing forwardPatternsR;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (options.gapMode == RAZERS_GAPPED)
	{
		resize(forwardPatternsL, pairCount, Exact());
		resize(forwardPatternsR, pairCount, Exact());
		for(unsigned i = 0; i < pairCount; ++i)
		{
			setHost(forwardPatternsL[i], readSetL[i]);
			setHost(forwardPatternsR[i], readSetR[i]);
			_patternMatchNOfPattern(forwardPatternsL[i], options.matchN);
			_patternMatchNOfPattern(forwardPatternsR[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsL[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsR[i], options.matchN);
		}
	}
#endif  // #ifdef RAZERS_BANDED_MYERS

	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;
	SEQAN_PROTIMESTART(find_time);
	
	for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId) {
		// lock to prevent releasing and loading the same contig twice
		// (once per _mapSingleReadsToContig call)
		lockContig(store, contigId);
		
		if (options.forward)
			_mapMatePairReadsParallel(store, contigId, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'F', options, mode);
		
		if (options.reverse)
			_mapMatePairReadsParallel(store, contigId, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'R', options, mode);
		
		unlockAndFreeContig(store, contigId);
	}

	options.timeMapReads = SEQAN_PROTIMEDIFF(find_time);
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:  " << options.countFiltration << ::std::endl;
		::std::cerr << "Verfication counter: " << options.countVerification << ::std::endl;
	}

    // Restore global state.
    omp_set_num_threads(oldMaxThreads);

	return 0;
}


} // End namespace

#endif
