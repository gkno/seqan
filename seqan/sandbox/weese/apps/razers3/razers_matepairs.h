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

#ifndef SEQAN_HEADER_RAZERS_MATEPAIRS_H
#define SEQAN_HEADER_RAZERS_MATEPAIRS_H

#include <seqan/misc/misc_dequeue.h>

namespace SEQAN_NAMESPACE_MAIN
{

// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.


//////////////////////////////////////////////////////////////////////////////
// Definitions

#ifdef RAZERS_MEMOPT

    template <typename TMPReadSet, typename TShape>
	struct SAValue< Index<TMPReadSet, IndexQGram<TShape> > > {
		typedef Pair<
			unsigned,				
			unsigned,
			BitCompressed<24, 8>	// max. 16M reads of length < 256
		> Type;
	};
	
#else

	template <typename TMPReadSet, typename TShape>
	struct SAValue< Index<TMPReadSet, IndexQGram<TShape> > > {
		typedef Pair<
			unsigned,			// many reads
			unsigned,			// of arbitrary length
			Compressed
		> Type;
	};

#endif

	
// template <typename TMPReadSet, typename TShape, typename TSpec>
// struct Cargo< Index<TMPReadSet, IndexQGram<TShape, TSpec> > > {
// 	typedef struct {
// 		double		abundanceCut;
// 		int			_debugLevel;
// 	} Type;
// };

#ifdef RAZERS_PRUNE_QGRAM_INDEX

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TMPReadSet, typename TShape>
inline bool _qgramDisableBuckets(Index<TMPReadSet, IndexQGram<TShape> > &index) 
{
	typedef Index<TMPReadSet, IndexQGram<TShape>	>	TReadIndex;
	typedef typename Fibre<TReadIndex, QGramDir>::Type	TDir;
	typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
	typedef typename Value<TDir>::Type					TSize;
	
	TDir &dir    = indexDir(index);
	bool result  = false;
	unsigned counter = 0;
	TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
	if (thresh < 100) thresh = 100;
	
	TDirIterator it = begin(dir, Standard());
	TDirIterator itEnd = end(dir, Standard());
	for (; it != itEnd; ++it)
		if (*it > thresh) 
		{
			*it = (TSize)-1;
			result = true;
			++counter;
		}
	
	if (counter > 0 && cargo(index)._debugLevel >= 1)
		::std::cerr << "Removed " << counter << " k-mers" << ::std::endl;
	
	return result;
}

#endif

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TFSSpec, typename TFSConfig, typename TRazerSOptions>
bool loadReads(
	FragmentStore<TFSSpec, TFSConfig>	& store,
	const char							* fileNameL,		// left mates file
	const char							* fileNameR,		// right mates file
	TRazerSOptions						& options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);

	MultiFasta leftMates;
	MultiFasta rightMates;

	if (!open(leftMates.concat, fileNameL, OPEN_RDONLY)) return false;
	if (!open(rightMates.concat, fileNameR, OPEN_RDONLY)) return false;

	AutoSeqFormat formatL;
	guessFormat(leftMates.concat, formatL);
	split(leftMates, formatL);

	AutoSeqFormat formatR;
	guessFormat(rightMates.concat, formatR);
	split(rightMates, formatR);

	unsigned seqCount = length(leftMates);
	if (seqCount != length(rightMates))
	if (options._debugLevel > 1) 
	{
		::std::cerr << "Numbers of mates differ: " << seqCount << "(left) != " << length(rightMates) << "(right).\n";
		return false;
	}

	String<Dna5Q>	seq[2];
	CharString		qual[2];
	CharString		id[2];
	
	unsigned kickoutcount = 0;
	unsigned maxReadLength = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0 || options.readNaming == 3)
		{
			assignSeqId(id[0], leftMates[i], formatL);				// read left Fasta id
			assignSeqId(id[1], rightMates[i], formatR);				// read right Fasta id
			if (options.readNaming == 0)
			{
				append(id[0], "/L");
				append(id[1], "/R");
			}
		}
		
		assignSeq(seq[0], leftMates[i], formatL);					// read left Read sequence
		assignSeq(seq[1], rightMates[i], formatR);					// read right Read sequence
		assignQual(qual[0], leftMates[i], formatL);					// read left ascii quality values  
		assignQual(qual[1], rightMates[i], formatR);				// read right ascii quality values  
		
		if (countN)
		{
			for (int j = 0; j < 2; ++j)
			{
				int maxBase = (int)(0.8 * length(seq[j]));
				int allowed[5] = 
					{ maxBase, maxBase, maxBase, maxBase, (int)(options.errorRate * length(seq[j]))};
				for (unsigned k = 0; k < length(seq[j]); ++k)
					if (--allowed[ordValue(getValue(seq[j], k))] == 0)
					{
//						std::cout << "Ignoring mate-pair: " << seq[0] << " " << seq[1] << std::endl;
						clear(seq[0]);
						clear(seq[1]);
						clear(id[0]);
						clear(id[1]);
						++kickoutcount;
						break;
					}
			}
		}
		
		for (int j = 0; j < 2; ++j)
		{
			// store dna and quality together
			assignQualities(seq[j], qual[j]);
			
			if (options.trimLength > 0 && length(seq[j]) > (unsigned)options.trimLength)
				resize(seq[j], options.trimLength);
		}
		appendMatePair(store, seq[0], seq[1], id[0], id[1]);
		if (maxReadLength < length(seq[0]))
			maxReadLength = length(seq[0]);
		if (maxReadLength < length(seq[1]))
			maxReadLength = length(seq[1]);
	}
	// memory optimization
	reserve(store.readSeqStore.concat, length(store.readSeqStore.concat), Exact());
//	reserve(store.readNameStore.concat, length(store.readNameStore.concat), Exact());

	typedef Shape<Dna, SimpleShape> TShape;
	typedef typename SAValue< Index<StringSet<Dna5String>, IndexQGram<TShape, OpenAddressing> > >::Type TSAValue;
	TSAValue sa(0, 0);
	sa.i1 = ~sa.i1;
	sa.i2 = ~sa.i2;
	
	if ((unsigned)sa.i1 < length(store.readSeqStore) - 1)
	{
		::std::cerr << "Maximal read number of " << (unsigned)sa.i1 + 1 << " exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}
	if ((unsigned)sa.i2 < maxReadLength - 1)
	{
		::std::cerr << "Maximal read length of " << (unsigned)sa.i2 + 1 << " bps exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}

	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality mate-pairs.\n";
	return (seqCount > 0);
}

	template <typename TFragmentStore>
	struct LessPairScore : 
		public ::std::binary_function <
			typename Value<typename TFragmentStore::TAlignedReadStore>::Type, 
			typename Value<typename TFragmentStore::TAlignedReadStore>::Type,
			bool >
	{
		TFragmentStore &mainStore;
		TFragmentStore &threadStore;
		
		LessPairScore(TFragmentStore &_mainStore, TFragmentStore &_threadStore):
			mainStore(_mainStore), threadStore(_threadStore) {}
		
		inline bool operator() (
			typename Value<typename TFragmentStore::TAlignedReadStore>::Type const &a, 
			typename Value<typename TFragmentStore::TAlignedReadStore>::Type const &b) const 
		{
			typedef typename TFragmentStore::TReadStore			TReadStore;
			typedef typename TFragmentStore::TAlignedReadStore	TAlignedReadStore;
			typedef typename TFragmentStore::TAlignQualityStore	TAlignQualityStore;
			typedef typename Value<TReadStore>::Type			TRead;
			typedef typename Value<TAlignedReadStore>::Type		TAlignedRead;
			typedef typename Value<TAlignQualityStore>::Type	TQual;
			typedef typename Id<TRead>::Type					TId;

			// pair number
			if (b.readId == TAlignedRead::INVALID_ID) return false;
			if (a.readId == TAlignedRead::INVALID_ID) return true;
			TRead const &ra = mainStore.readStore[a.readId];
			TRead const &rb = mainStore.readStore[b.readId];
			if (ra.matePairId < rb.matePairId) return true;
			if (ra.matePairId > rb.matePairId) return false;

			// quality
			if (a.id == TAlignedRead::INVALID_ID) return false;
			if (b.id == TAlignedRead::INVALID_ID) return true;
			TQual const &qa = threadStore.alignQualityStore[a.id];
			TQual const &qb = threadStore.alignQualityStore[b.id];
			if (qa.pairScore > qb.pairScore) return true;
			if (qa.pairScore < qb.pairScore) return false;
			
			return a.pairMatchId < b.pairMatchId;
		}
	};
    
    template <typename TFragmentStore, typename TReadMatch>
	struct LessPairErrors : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
        typedef typename TFragmentStore::TReadStore			TReadStore;
        typedef typename Value<TReadStore>::Type			TRead;

        TFragmentStore const & mainStore;

        LessPairErrors(TFragmentStore const & mainStore_)
                : mainStore(mainStore_)
        {}

		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
            if (b.readId == TReadMatch::INVALID_ID) return false;
            if (a.readId == TReadMatch::INVALID_ID) return true;
			TRead const &ra = mainStore.readStore[a.readId];
			TRead const &rb = mainStore.readStore[b.readId];
			if (ra.matePairId < rb.matePairId) return true;
			if (ra.matePairId > rb.matePairId) return false;
            
			// quality
            if (a.orientation == '-') return false;
            if (b.orientation == '-') return true;
			if (a.pairScore > b.pairScore) return true;
			if (a.pairScore < b.pairScore) return false;

			if (a.pairMatchId < b.pairMatchId) return true;
			if (a.pairMatchId > b.pairMatchId) return false;

            return a.readId < b.readId;
		}
	};
	
    

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TFragmentStore, typename TMatches, typename TCounts, typename TSpec, typename TSwiftL, typename TSwiftR >
void compactPairMatches(
	TFragmentStore			& store,   // all but aligned reads up to the global writeback
	TMatches                & matches, // aligned read only
	TCounts					&, 
	RazerSOptions<TSpec>	& options, 
	TSwiftL					& swiftL, 
	TSwiftR					& swiftR)
{
	typedef typename Value<TMatches>::Type                          TMatch;
	typedef typename Iterator<TMatches, Standard>::Type             TIterator;
	
    // fprintf(stderr, "[pair-compact]");
    double beginTime = sysTime();
	unsigned matePairId = -2;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreDistCutOff = MinValue<int>::VALUE;

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;
    unsigned disabled = 0;
    
	// sort 
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
//	sortAlignedReads(threadStore, LessPairScore<TFragmentStore>(mainStore, threadStore));
    ::std::sort(it, itEnd, LessPairErrors<TFragmentStore, TMatch>(store));
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

	for (; it != itEnd; ++it) 
	{
        // std::cerr << *it << std::endl;
        // std::cerr << *(it + 1) << std::endl;
        // SEQAN_ASSERT(it->pairMatchId == (it + 1)->pairMatchId);
		if ((*it).orientation == '-' || (*it).readId == TMatch::INVALID_ID) continue;
		if (matePairId == store.readStore[(*it).readId].matePairId)
		{ 
			if (it->pairScore <= scoreDistCutOff) continue;
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					int maxErrors = -1 - (*it).pairScore;
					if (options.purgeAmbiguous)
						maxErrors = -1;

					setMaxErrors(swiftL, matePairId, maxErrors);
					setMaxErrors(swiftR, matePairId, maxErrors);

					if (maxErrors == -1 && options._debugLevel >= 2)
                        disabled += 1;
						// ::std::cerr << "(pair #" << matePairId << " disabled)";

					if (options.purgeAmbiguous)
						dit = ditBeg;
				}
#endif
				continue;
			}
		}
		else
		{
			matePairId = store.readStore[(*it).readId].matePairId;
			hitCount = 0;
			if (options.scoreDistanceRange > 0)
				scoreDistCutOff = (*it).pairScore - options.scoreDistanceRange;
			
			ditBeg = dit;
		}
		*dit = *it;	++dit; ++it;
		*dit = *it;	++dit;
	}
	unsigned origSize = length(matches);
	resize(matches, dit - begin(matches, Standard()));
//	compactAlignedReads(matches);

    options.timeCompactMatches += sysTime() - beginTime;

    fprintf(stderr, "[%u reads disabled]", disabled);
	unsigned newSize = length(matches);
	fprintf(stderr, "[%u of %u alignments removed]", unsigned(origSize - newSize), unsigned(origSize));
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
    typename TMatches,
	typename TFSSpec,
	typename TFSConfig,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TPreprocessing,
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode>
void _mapMatePairReads(
    TMatches                                & matches,
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

    typedef typename Value<TMatches>::Type TMatch;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinderL;
	typedef Finder<TGenomeInf, Swift<TSwiftSpec> >			TSwiftFinderR;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// MATE-PAIR FILTRATION
	typedef Pair<__int64, TMatch>	                        TDequeueValue;
	typedef Dequeue<TDequeueValue>							TDequeue;
	typedef typename TDequeue::TIter						TDequeueIterator;

	// VERIFICATION
	typedef MatchVerifier <
		TFragmentStore, 
        TMatches,
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
	TVerifier	verifierL(matches, options, swiftPatternL, cnts);
	TVerifier	verifierR(matches, options, swiftPatternR, cnts);

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

	resize(lastPotMatchNo, length(host(swiftPatternL)), (__int64)-2, Exact());

	TSize gLength = length(genome);

	TMatch mR;
	TDequeueValue fL(-2, mR);	// to supress uninitialized warnings

    // Iterate over all filtration results are returned by SWIFT.
	while (find(swiftFinderR, swiftPatternR, options.errorRate)) 
	{
		++options.countFiltration;

        CharString pref = prefix(store.readNameStore[2 * swiftPatternR.curSeqNo + 1], length("EAS20_8_6_1_248_1397"));
        CharString s = "EAS20_8_6_1_248_1397";
        if (pref == s)
            std::cerr << "GOTCHA" << std::endl;
        // #pragma omp critical
        // {
        //     std::cerr << tls.globalStore->readNameStore[2 * (threadIdOffset + itR->ndlSeqNo) + 1] << std::endl;
        // }

#ifdef RAZERS_DEBUG_MATEPAIRS
        std::cerr << "\nSWIFT\tR\t" << swiftPatternR.curSeqNo << "\t" << store.readNameStore[2 * swiftPatternR.curSeqNo + 1] << "\t" << scanShift + beginPosition(swiftFinderR) << "\t" << scanShift + endPosition(swiftFinderR) << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
        
		unsigned matePairId = swiftPatternR.curSeqNo;
		TGPos rEndPos = endPosition(swiftFinderR) + scanShift;
		TGPos doubleParWidth = 2 * (*swiftFinderR.curHit).bucketWidth;
		
        // (1) Remove out-of-window left mates from fifo.
		while (!empty(fifo) && (TSignedGPos)front(fifo).i2.endPos + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
		{
#ifdef RAZERS_DEBUG_MATEPAIRS
            if (front(fifo).i2.readId > length(store.readNameStore))
                std::cerr << "\nPOP\tL\t" << "[bad read]" << "\t" << front(fifo).i2.beginPos << "\t" << front(fifo).i2.endPos << std::endl;
            else
                std::cerr << "\nPOP\tL\t" << store.readNameStore[front(fifo).i2.readId & ~NOT_VERIFIED] << "\t" << front(fifo).i2.beginPos << "\t" << front(fifo).i2.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
			popFront(fifo);
			++firstNo;
		}

        // (2) Add within-window left mates to fifo.
		while (empty(fifo) || (TSignedGPos)back(fifo).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
		{
			if (find(swiftFinderL, swiftPatternL, options.errorRate))
			{
				++options.countFiltration;
#ifdef RAZERS_DEBUG_MATEPAIRS
                std::cerr << "\nSWIFT\tL\t" << swiftPatternL.curSeqNo << "\t" << store.readNameStore[2 * swiftPatternL.curSeqNo] << "\t" << beginPosition(swiftFinderL) << "\t" << endPosition(swiftFinderL) << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
				gPair = positionRange(swiftFinderL);
				if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
				{
					// link in
					fL.i1 = lastPotMatchNo[swiftPatternL.curSeqNo];
					lastPotMatchNo[swiftPatternL.curSeqNo] = lastNo++;
					
					fL.i2.readId = store.matePairStore[swiftPatternL.curSeqNo].readId[0] | NOT_VERIFIED;
					fL.i2.beginPos = beginPosition(swiftFinderL);
					fL.i2.endPos = gPair.i2;
					
					pushBack(fifo, fL);
				}
			} else {
				break;
            }
		}

		int	bestLeftScore = MinValue<int>::VALUE;
		int bestLibSizeError = MaxValue<int>::VALUE;
		TDequeueIterator bestLeft = TDequeueIterator();

        bool rightVerified = false;
		TDequeueIterator it;
		unsigned leftReadId = store.matePairStore[matePairId].readId[0];
		__int64 last = (__int64)-1;
		__int64 lastValid = (__int64)-1;
        __int64 i;
		for (i = lastPotMatchNo[matePairId]; firstNo <= i; last = i, i = (*it).i1)
		{
			it = &value(fifo, i - firstNo);

			// search left mate
//			if (((*it).i2.readId & ~NOT_VERIFIED) == leftReadId)	
//			        ^== we need not to test anymore, as only corr. left mates are traversed
//						via the linked list beginning from lastPotMatchNo[matePairId] 
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
#ifdef RAZERS_BANDED_MYERS
                        verifierL.patternState.leftClip = ((*it).i2.beginPos >= 0)? 0: -(*it).i2.beginPos;	// left clip if match begins left of the genome
#endif
#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "\nVERIFY\tL\t" << matePairId << "\t" << store.readNameStore[2 * matePairId] << "\t" << (TSignedGPos)(*it).i2.beginPos << "\t" << (*it).i2.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                        ++options.countVerification;
                        if (matchVerify(verifierL, infix(genome, ((*it).i2.beginPos >= 0)? (TSignedGPos)(*it).i2.beginPos: (TSignedGPos)0, (TSignedGPos)(*it).i2.endPos), 
                                        matePairId, readSetL, mode))
						{
#ifdef RAZERS_DEBUG_MATEPAIRS
                            std::cerr << "  YES: " << verifierL.m.beginPos << "\t" << verifierL.m.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS

							verifierL.m.readId = (*it).i2.readId & ~NOT_VERIFIED;		// has been verified positively
							(*it).i2 = verifierL.m;
						} else {
							(*it).i2.readId = ~NOT_VERIFIED;				// has been verified negatively
#ifdef RAZERS_DEBUG_MATEPAIRS
                            std::cerr << "  NO" << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
							continue;										// we intentionally do not set lastPositive to i
                        }													// to remove i from linked list
					} else {
						lastValid = i;
						continue;											// left pot. hit is out of tolerance window
                    }
				} //else {}													// left match is verified already

				// short-cut negative matches
				if (last != lastValid)
				{
					SEQAN_ASSERT_NEQ(lastValid, i);
					if (lastValid == (__int64)-1)
						lastPotMatchNo[matePairId] = i;
					else
						value(fifo, lastValid - firstNo).i1 = i;
				}
				lastValid = i;

				if (!rightVerified)											// here a verfied left match is available
				{
#ifdef RAZERS_DEBUG_MATEPAIRS
					std::cerr << "\nVERIFY\tR\t" << matePairId << "\t" << store.readNameStore[2 * matePairId + 1] << "\t" << beginPosition(swiftFinderR) << "\t" << endPosition(swiftFinderR) << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                    ++options.countVerification;
                    if (pref == s)
                        std::cerr << "\nVERIFY\tR\t" << matePairId << "\t" << store.readNameStore[2 * matePairId + 1] << "\t" << beginPosition(swiftFinderR) << "\t" << endPosition(swiftFinderR) << std::endl;
					if (matchVerify(verifierR, infix(swiftFinderR), matePairId, readSetR, mode)) {
#ifdef RAZERS_DEBUG_MATEPAIRS
						std::cerr << "  YES: " << verifierR.m.beginPos << "\t" << verifierR.m.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
						rightVerified = true;
						mR = verifierR.m;
                        if (pref == s)
                            std::cerr << "BREAK HERE" << std::endl;
					} else {
#ifdef RAZERS_DEBUG_MATEPAIRS
						std::cerr << "  NO" << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
						// Break out of lastPotMatch loop, rest of find(right SWIFT results loop will not
						// be executed since bestLeftScore remains untouched.
						i = (*it).i1;
						break;
					}
				}
				
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
					int score = (*it).i2.score;
					if (bestLeftScore <= score)
					{
                        // distance between left mate beginning and right mate end
                        __int64 dist = (__int64)verifierR.m.endPos - (__int64)(*it).i2.beginPos;
                        
						int libSizeError = options.libraryLength - dist;
#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "    libSizeError = " << libSizeError << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
						if (libSizeError < 0)
                            libSizeError = -libSizeError;
                        if (libSizeError > options.libraryError)
                            continue;
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
							// if (bestLeftScore == 0) break;	// TODO: replace if we have real qualities
						}
					}
				}
			}
		}

        // (3) Short-cut negative matches.
		if (last != lastValid)
		{
			SEQAN_ASSERT_NEQ(lastValid, i);
			if (lastValid == (__int64)-1)
				lastPotMatchNo[matePairId] = i;
			else
				value(fifo, lastValid - firstNo).i1 = i;
		}
		
		// verify right mate, if left mate matches
		if (bestLeftScore != MinValue<int>::VALUE)
		{
// //			if (matchVerify(
// //					mR, qR, infix(swiftFinderR),
// //					matePairId, readSetR, forwardPatternsR,
// //					options, TSwiftSpec()))
            // XXX
// 			if (matchVerify(verifierR, infix(swiftFinderR), 
// 					matePairId, readSetR, mode))
// 			{
            // XXX
                // XXX
				// distance between left mate beginning and right mate end
				// __int64 dist = (__int64)verifierR.m.endPos - (__int64)(*bestLeft).i2.beginPos;
                // std::cerr << "    COND = " << (dist <= options.libraryLength + options.libraryError) << " && " << (options.libraryLength <= dist + options.libraryError) << std::endl;
                // std::cerr << "    dist = " << dist << std::endl;
                // XXX
                // XXX
				// if (dist <= options.libraryLength + options.libraryError &&
				// 	options.libraryLength <= dist + options.libraryError)
				// {
                // XXX
					// mR = verifierR.m;
					// qR = verifierR.q;
					
					fL.i2 = (*bestLeft).i2;
					
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
						// dist = -dist;
					}
					
					// set a unique pair id
					fL.i2.pairMatchId = mR.pairMatchId = options.nextPairMatchId;
					if (++options.nextPairMatchId == TAlignedRead::INVALID_ID)
						options.nextPairMatchId = 0;

					// score the whole match pair
					fL.i2.pairScore = mR.pairScore = fL.i2.score + mR.score;

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
						appendValue(matches, fL.i2, Generous());
                        back(matches).orientation = orientation;
						appendValue(matches, mR, Generous());
                        back(matches).orientation = (orientation == 'F') ? 'R' : 'F';

                        if (pref == s) {
                            std::cerr << "ADDED" << std::endl;
                            std::cerr << "  (" << mR.beginPos << ", " << mR.endPos << ")" << std::endl;
                            std::cerr << "  (" << fL.i2.beginPos << ", " << fL.i2.endPos << ")" << std::endl;
                        }
#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "\nHIT\tL\t" << fL.i2.readId << "\t" << store.readNameStore[fL.i2.readId] << "\t" << fL.i2.beginPos << "\t" << fL.i2.endPos << std::endl;
                        std::cerr << "\nHIT\tR\t" << mR.readId << "\t" << store.readNameStore[mR.readId] << "\t" << mR.beginPos << "\t" << mR.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS

						if (length(store.alignedReadStore) > options.compactThresh)
						{
							typename Size<TAlignedReadStore>::Type oldSize = length(store.alignedReadStore);
                            // TODO(weese): Duplicates are hard to mask in paired-end mode.
                            // if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
                            //   maskDuplicates(matches);	// overlapping parallelograms cause duplicates
							compactPairMatches(store, matches, cnts, options, swiftPatternL, swiftPatternR);
							
							if (length(store.alignedReadStore) * 4 > oldSize)			// the threshold should not be raised
							  options.compactThresh *= options.compactMult;
								//options.compactThresh += (options.compactThresh >> 1);	// if too many matches were removed
							
							if (options._debugLevel >= 2)
								::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
						}
					}
                    // XXX
				// }
                // XXX
			}
            // XXX
		// }
        // XXX
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
int _mapMatePairReads(
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
	typedef Pattern<TRead const, MyersUkkonen>			TMyersPattern;	// verifier

    typedef typename TFragmentStore::TContigSeq TContigSeq;
    typedef typename Position<TContigSeq>::Type TContigPos;
    typedef MatchRecord<TContigPos> TMatchRecord;

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

    // Verifier preprocessing.
#ifdef RAZERS_BANDED_MYERS
	typedef Nothing TPreprocessing;
	TPreprocessing forwardPatternsL;
	TPreprocessing forwardPatternsR;
#else  // #ifdef RAZERS_BANDED_MYERS
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
	
	options.timeDumpResults = 0;
	SEQAN_PROTIMESTART(find_time);

    // We collect the matches in a more compact data structure than the
    // AlignedReadStoreElement from FragmentStore.
    String<TMatchRecord> matches;

	for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId) {
		// lock to prevent releasing and loading the same contig twice
		// (once per _mapSingleReadsToContig call)
		lockContig(store, contigId);
		
		if (options.forward)
			_mapMatePairReads(matches, store, contigId, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'F', options, mode);
		
		if (options.reverse)
			_mapMatePairReads(matches, store, contigId, swiftPatternL, swiftPatternR, forwardPatternsL, forwardPatternsR, cnts, 'R', options, mode);
		
		unlockAndFreeContig(store, contigId);
	}

    double beginCopyTime = sysTime();
    // Final compact matches.
    compactPairMatches(store, matches, cnts, options, swiftPatternL, swiftPatternR);
    // Write back to store.
    reserve(store.alignedReadStore, length(matches));
    reserve(store.alignQualityStore, length(matches));
    typedef typename Iterator<String<TMatchRecord>, Standard>::Type TIterator;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
	typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
    for (TIterator it = begin(matches), itEnd = end(matches); it != itEnd; ++it) {
        SEQAN_ASSERT_NEQ(it->orientation, '-');
        SEQAN_ASSERT(!(it->orientation == 'F') || (it->beginPos <= it->endPos));  // implication
        SEQAN_ASSERT(!(it->orientation == 'R') || (it->beginPos >= it->endPos));  // implication
        // if (it->orientation == 'R')
        //     ::std::swap(it->beginPos, it->endPos);
        appendValue(store.alignedReadStore, TAlignedReadStoreElem(length(store.alignQualityStore), it->readId, it->contigId, it->beginPos, it->endPos));
        back(store.alignedReadStore).pairMatchId = it->pairMatchId;
        appendValue(store.alignQualityStore, TAlignedQualStoreElem(it->pairScore, it->score, -it->score));
    }
    options.timeFsCopy = sysTime() - beginCopyTime;

	// restore original orientation (R-reads are infixes of ConcatDirect StringSet)
	reverseComplement(readSetR);

	options.timeMapReads = SEQAN_PROTIMEDIFF(find_time);
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << "Time for copying back            \t" << options.timeFsCopy << " seconds" << ::std::endl;
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:  " << options.countFiltration << ::std::endl;
		::std::cerr << "Verfication counter: " << options.countVerification << ::std::endl;
	}

	return 0;
}


} // End namespace

#endif