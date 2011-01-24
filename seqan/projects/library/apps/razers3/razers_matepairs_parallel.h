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
// We require mate-pairs to be stored together in one read string.
// Pair i has mates at positions 2*i and 2*i+1 in the read string.

#ifndef SEQAN_HEADER_RAZERS_MATEPAIRS_PARALLEL_H
#define SEQAN_HEADER_RAZERS_MATEPAIRS_PARALLEL_H

#include <numeric>

#include <seqan/misc/misc_dequeue.h>

#include "razers_parallel.h"
#include "razers_matepairs.h"

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// Stores the results of the verification.
//
// Put into its own class so it can be locked independently of other class
// members.
template <typename TFragmentStore>
class PairedVerificationResults
{
public:
    String<TFragmentStore *> localStores;
    Lock<Omp> * lock;

    PairedVerificationResults() : lock(new Lock<Omp>()) {}

    PairedVerificationResults(PairedVerificationResults const & other)
            : localStores(other.localStores), lock(new Lock<Omp>())
    {
        // Not thread-safe copying since this is only used at the beginning when resizing block local storages string.
    }

    PairedVerificationResults & operator=(PairedVerificationResults const & other)
    {
        if (this == &other)
            return *this;
        localStores = other.localStores;
    }

    ~PairedVerificationResults()
    {
        delete lock;
    }
};

template <typename TFragmentStore, typename TSwiftFinderL, typename TSwiftFinderR, typename TSwiftPattern, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename TPreprocessing>
struct MapPairedReads {};

// ThreadLocalStorage specialization for single-end read mapping in RazerS.
template <typename TFragmentStore, typename TSwiftFinderL_, typename TSwiftFinderR_, typename TSwiftPattern_, typename TShape_/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename TPreprocessing>
class ThreadLocalStorage<MapPairedReads<TFragmentStore, TSwiftFinderL_, TSwiftFinderR_, TSwiftPattern_, TShape_, TOptions, TCounts, TRazerSMode, TPreprocessing> >
{
public:
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type						TRead;
	typedef StringSet<TRead>										TReadSet;
    typedef TShape_ TShape;

    typedef TSwiftPattern_ TSwiftPattern;
    typedef TSwiftFinderL_ TSwiftFinderL;
    typedef TSwiftFinderR_ TSwiftFinderR;
    
	typedef typename TFragmentStore::TContigSeq				TGenome;
    typedef typename Infix<TGenome>::Type TGenomeInfix;

    // The id of this thread.
    unsigned threadId;

    // Each thread needs its local options since the compactionThreshold is changed.
    // TODO(holtgrew): Change overall program structure so this is factorized out of the options struct.
    TOptions options;
    TOptions /*const*/ * globalOptions;

    // Split the read seq store from fragment store into two parts.
    TReadSet readSetL, readSetR; 

    // Each thread has its own SWIFT finder and pattern object.
    TSwiftFinderL swiftFinderL;
    TSwiftFinderR swiftFinderR;
    TSwiftPattern swiftPatternL, swiftPatternR;

    TCounts counts;  // TODO(holtgrew): Artifact?

    TFragmentStore store;
    TFragmentStore /*const*/ * globalStore;

    TShape shape;

    typedef MatchVerifier<TFragmentStore, TOptions, TRazerSMode, TSwiftPattern, TCounts, TPreprocessing> TMatchVerifier;
    TMatchVerifier verifierL, verifierR;

    // Mailbox for the verification results.
    PairedVerificationResults<TFragmentStore> verificationResults;

    String<unsigned> splitters;

    unsigned completeWindows;

    TGenomeInfix genomeInf;

    // The last potential match fifo from the verification can be thread
    // local.
	typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
	typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
	typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
	typedef Triple<__int64, TAlignedRead, TAlignQuality>	TDequeueValue;
	typedef Dequeue<TDequeueValue>							TDequeue;
	typedef typename TDequeue::TIter						TDequeueIterator;
    TDequeue fifo;                      // stores left-mate potential matches
    String<__int64> fifoLastPotMatchNo; // last number of a left-mate potential
    __int64 fifoLastNo;                 // last number over all left-mate pot. matches in the queue
    __int64 fifoFirstNo;                // first number over all left-mate pot. match in the queue

    ThreadLocalStorage() : fifoLastNo(0), fifoFirstNo(0) {}
};

template <typename TThreadLocalStorage>
class SwiftPatternLSetMaxErrorsWrapper
{
public:
    TThreadLocalStorage & tls;
    SwiftPatternLSetMaxErrorsWrapper(TThreadLocalStorage & tls_) : tls(tls_) {}
};

template <typename TThreadLocalStorage>
class SwiftPatternRSetMaxErrorsWrapper
{
public:
    TThreadLocalStorage & tls;
    SwiftPatternRSetMaxErrorsWrapper(TThreadLocalStorage & tls_) : tls(tls_) {}
};

template <typename TFragmentStore, typename THitString, typename TOptions, typename TSwiftPattern>
struct PairedVerification;

template <typename THitString>
void
buildHitSplittersAndPartitionHits(String<size_t> & splitters, THitString & hitString, unsigned packageCount, unsigned matePairCount)
{
    typedef typename Iterator<THitString>::Type THitStringIterator;

    (void)matePairCount;

    // TODO(holtgrew): Optimize with packageCount power of 2 and/or use libdiv?
    // TODO(holtgrew): Or maybe logarithmic search faster than modulo?

    // Partition hitString into buckets by ndlSeqNo % packageCount.
    //
    // First, build counters.
    clear(splitters);
    resize(splitters, packageCount + 1, 0);
    splitters[0] = 0;
    for (THitStringIterator it = begin(hitString, Standard()), itEnd = end(hitString, Standard()); it != itEnd; ++it) {
        SEQAN_ASSERT_LEQ(it->ndlSeqNo, matePairCount);
        splitters[it->ndlSeqNo % packageCount + 1] += 1;
        SEQAN_ASSERT_LEQ(splitters[it->ndlSeqNo % packageCount + 1], length(hitString));
    }
    std::partial_sum(begin(splitters, Standard()), end(splitters, Standard()), begin(splitters, Standard()));
    // Second, copy into temporary buffer.
    String<size_t> offsets(splitters);
    THitString buffer;
    resize(buffer, length(hitString));
    for (THitStringIterator it = begin(hitString, Standard()), itEnd = end(hitString, Standard()); it != itEnd; ++it) {
        unsigned idx = it->ndlSeqNo % packageCount;
        buffer[offsets[idx]] = *it;
        offsets[idx] += 1;
        if (idx < length(offsets) - 1)
            SEQAN_ASSERT_LEQ(offsets[idx], offsets[idx + 1]);
        if (idx > 0u)
            SEQAN_ASSERT_LEQ(offsets[idx - 1], offsets[idx]);
    }

    // Finally, write out results.
    swap(hitString, buffer);

    // std::cout << "SPLITTERS: ";
    // for (unsigned i = 0; i < length(splitters); ++i) {
    //     std::cout << splitters[i] << " ";
    // }
    // std::cout << " last should be " << length(hitString) << std::endl;
}


template <typename TFragmentStore, typename THitString_, typename TOptions, typename TSwiftPattern>
class Job<PairedVerification<TFragmentStore, THitString_, TOptions, TSwiftPattern> >
{
public:
    typedef PairedVerificationResults<TFragmentStore> TVerificationResults;
    typedef THitString_ THitString;
    typedef std::tr1::shared_ptr<THitString> THitStringPtr;

    int threadId;
    TVerificationResults * verificationResults;
    TFragmentStore * globalStore;
    unsigned contigId;
    char orientation;
    // Hit strings of previous window for the left finder.
    THitStringPtr prevHitsPtrL;
    size_t prevHitsLBegin, prevHitsLEnd;
    // Current hit string of current finder.
    THitStringPtr hitsPtrL;
    size_t hitsLBegin, hitsLEnd;
    THitStringPtr hitsPtrR;
    size_t hitsRBegin, hitsREnd;
    __int64 rightWindowBegin;
    TOptions * options;
    TSwiftPattern * swiftPatternL;
    TSwiftPattern * swiftPatternR;

    Job() {}

    Job(int threadId_, TVerificationResults & verificationResults_, TFragmentStore & globalStore_, unsigned contigId_, char orientation_, THitStringPtr & prevHitsPtrL_, size_t prevHitsLBegin_, size_t prevHitsLEnd_, THitStringPtr & hitsPtrL_, size_t hitsLBegin_, size_t hitsLEnd_, THitStringPtr & hitsPtrR_, size_t hitsRBegin_, size_t hitsREnd_, __int64 rightWindowBegin_, TOptions & options_, TSwiftPattern & swiftPatternL_, TSwiftPattern & swiftPatternR_)
            : threadId(threadId_), verificationResults(&verificationResults_), globalStore(&globalStore_), contigId(contigId_), orientation(orientation_), prevHitsPtrL(prevHitsPtrL_), prevHitsLBegin(prevHitsLBegin_), prevHitsLEnd(prevHitsLEnd_), hitsPtrL(hitsPtrL_), hitsLBegin(hitsLBegin_), hitsLEnd(hitsLEnd_), hitsPtrR(hitsPtrR_), hitsRBegin(hitsRBegin_), hitsREnd(hitsREnd_), rightWindowBegin(rightWindowBegin_), options(&options_), swiftPatternL(&swiftPatternL_), swiftPatternR(&swiftPatternR_)
    {
        SEQAN_ASSERT_LEQ(hitsLBegin, length(*hitsPtrL));
        SEQAN_ASSERT_LEQ(hitsLEnd, length(*hitsPtrL));
        SEQAN_ASSERT_LEQ(hitsRBegin, length(*hitsPtrR));
        SEQAN_ASSERT_LEQ(hitsREnd, length(*hitsPtrR));
        if (prevHitsPtrL.get() != 0) {
            SEQAN_ASSERT_LEQ(prevHitsLBegin, length(*prevHitsPtrL));
            SEQAN_ASSERT_LEQ(prevHitsLEnd, length(*prevHitsPtrL));
        }
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================


// Allow disabling reads in compactPairMatches() for left-mate read set.
//
// We do not disable the read right 
template <typename TThreadLocalStorage, typename TReadNo, typename TMaxErrors>
void
setMaxErrors(SwiftPatternLSetMaxErrorsWrapper<TThreadLocalStorage> & wrapper,
             TReadNo readNo,
             TMaxErrors maxErrors)
{
    // std::cerr << std::endl << "SET MAX ERRORS LEFT" << std::endl;
    SEQAN_ASSERT_LEQ(readNo, wrapper.tls.splitters[wrapper.tls.threadId + 1] - 1);
    SEQAN_ASSERT_GEQ(wrapper.tls.splitters[wrapper.tls.threadId], readNo);
	int localReadNo = readNo - wrapper.tls.splitters[wrapper.tls.threadId];

	int minT = _qgramLemma(wrapper.tls.swiftPatternL, localReadNo, maxErrors);
	if (minT > 1){
		if (maxErrors < 0)
            minT = MaxValue<int>::VALUE;
		setMinThreshold(wrapper.tls.swiftPatternL, localReadNo, static_cast<unsigned>(minT));
	}
}

// Allow disabling reads in compactPairMatches() for left-mate read set.
//
// We do not disable the read right 
template <typename TThreadLocalStorage, typename TReadNo, typename TMaxErrors>
void
setMaxErrors(SwiftPatternRSetMaxErrorsWrapper<TThreadLocalStorage> & wrapper,
             TReadNo readNo,
             TMaxErrors maxErrors)
{
    // std::cerr << std::endl << "SET MAX ERRORS LEFT" << std::endl;
    SEQAN_ASSERT_LEQ(readNo, wrapper.tls.splitters[wrapper.tls.threadId + 1] - 1);
    SEQAN_ASSERT_GEQ(wrapper.tls.splitters[wrapper.tls.threadId], readNo);
	int localReadNo = readNo - wrapper.tls.splitters[wrapper.tls.threadId];

	int minT = _qgramLemma(wrapper.tls.swiftPatternR, localReadNo, maxErrors);
	if (minT > 1){
		if (maxErrors < 0)
            minT = MaxValue<int>::VALUE;
		setMinThreshold(wrapper.tls.swiftPatternR, localReadNo, static_cast<unsigned>(minT));
	}
}

template <typename TFragmentStore>
inline
void
appendToVerificationResults(PairedVerificationResults<TFragmentStore> & verificationResults, TFragmentStore * storePtr)
{
    omp_set_lock(&verificationResults.lock->lock_);
    appendValue(verificationResults.localStores, storePtr);
    omp_unset_lock(&verificationResults.lock->lock_);
}

template <typename TThreadLocalStorages, typename TFragmentStore, typename TSplitters, typename TShape, typename TOptions>
void initializeThreadLocalStoragesPaired(TThreadLocalStorages & threadLocalStorages,
                                         TFragmentStore /*const*/ & store,
                                         TSplitters const & splitters,
                                         TShape /*const*/ & shape,
                                         TOptions /*const*/ & options)
{
    SEQAN_ASSERT_GT(length(splitters), 1u);
    int threadCount = length(splitters) - 1;

    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    typedef typename TThreadLocalStorage::TSwiftPattern TSwiftPattern;
    typedef typename Host<TSwiftPattern>::Type TIndex;
    typedef typename Position<typename TFragmentStore::TContigStore>::Type TPosition;
    typedef typename TThreadLocalStorage::TReadSet TReadSet;

    resize(threadLocalStorages, threadCount);
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < threadCount; ++i) {
        TThreadLocalStorage & tls = threadLocalStorages[i];

        // Initialize properties that are simple to set.
        tls.threadId = i;
        tls.globalStore = &store;
        tls.shape = shape;
        tls.options = options;  // TODO(holtgrew): Copy for stats and threshold, really good?
        tls.globalOptions = &options;
        tls.splitters = splitters;

        // Split mate-pairs over two indices.
        unsigned pairCount = splitters[i + 1] - splitters[i];
        TReadSet & readSetL = tls.readSetL;
        TReadSet & readSetR = tls.readSetR;
        resize(readSetL, pairCount, Exact());
        resize(readSetR, pairCount, Exact());

        unsigned offset = splitters[i];
        for (unsigned j = 0; j < pairCount; ++j) {
            assign(readSetL[j], store.readSeqStore[store.matePairStore[offset + j].readId[0]]);
            assign(readSetR[j], store.readSeqStore[store.matePairStore[offset + j].readId[1]]);
        }
        reverseComplement(readSetR);

        // Clear patterns and set parameters.
        TSwiftPattern & swiftPatternL = tls.swiftPatternL;
        clear(swiftPatternL);
        swiftPatternL.params.minThreshold = options.threshold;
        swiftPatternL.params.tabooLength = options.tabooLength;
        swiftPatternL.params.printDots = false;
        TSwiftPattern & swiftPatternR = tls.swiftPatternR;
        clear(swiftPatternR);
        swiftPatternR.params.minThreshold = options.threshold;
        swiftPatternR.params.tabooLength = options.tabooLength;
        swiftPatternR.params.printDots = (tls.threadId == 0) && (tls.options._debugLevel > 0);

        // Initialize the indices.
        // TODO(holtgrew): Necessary to split into readSetL and readSetR if we assign separately anyway?
        TIndex & indexL = host(tls.swiftPatternL);
        clear(indexL);
        clear(indexText(indexL));
        reserve(indexText(indexL), length(readSetL), Exact());
        for (TPosition j = 0, jEnd = length(readSetL); j < jEnd; ++j)
            appendValue(indexText(indexL), readSetL[j]);
        indexL.shape = shape;
#ifdef RAZERS_OPENADDRESSING
        indexL.alpha = options.loadFactor;
#endif
        cargo(indexL).abundanceCut = options.abundanceCut;
        cargo(indexL)._debugLevel = options._debugLevel;
        indexRequire(indexL, QGramSADir());

        TIndex & indexR = host(tls.swiftPatternR);
        clear(indexR);
        clear(indexText(indexR));
        reserve(indexText(indexR), length(readSetR), Exact());
        for (TPosition j = 0, jEnd = length(readSetR); j < jEnd; ++j)
            appendValue(indexText(indexR), readSetR[j]);
        indexR.shape = shape;
#ifdef RAZERS_OPENADDRESSING
        indexR.alpha = options.loadFactor;
#endif
        cargo(indexR).abundanceCut = options.abundanceCut;
        cargo(indexR)._debugLevel = options._debugLevel;
        indexRequire(indexR, QGramSADir());
    }
}

template <typename TFragmentStore, typename TSwiftFinderL, typename TSwiftFinderR, typename TSwiftPattern, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename THitString, typename TPreprocessing>
void workVerification(ThreadLocalStorage<MapPairedReads<TFragmentStore, TSwiftFinderL, TSwiftFinderR, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > & tls,
                      Job<PairedVerification<TFragmentStore, THitString, TOptions, TSwiftPattern> > & job,
                      String<unsigned> const & splitters)
{
    typedef ThreadLocalStorage<MapPairedReads<TFragmentStore, TSwiftFinderL, TSwiftFinderR, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > TThreadLocalStorage;

    typedef typename TThreadLocalStorage::TReadSet TReadSet;

	typedef typename TFragmentStore::TContigSeq				TGenome;

	typedef typename TFragmentStore::TMatePairStore			TMatePairStore;
	typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
	typedef typename Value<TMatePairStore>::Type			TMatePair;
	typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
	typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
    typedef typename TThreadLocalStorage::TReadSet TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape>	>	TReadIndex;
	typedef typename Id<TAlignedRead>::Type					TId;

	typedef typename TFragmentStore::TContigSeq				TGenome;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Position<TGenome>::Type				TGPos;
	typedef typename MakeSigned_<TGPos>::Type				TSignedGPos;
	typedef typename Infix<TGenome>::Type					TGenomeInf;

    typedef typename Iterator<THitString, Standard>::Type THitStringIter;

	typedef MatchVerifier <
		TFragmentStore, 
		TOptions, 
		TRazerSMode, 
		TSwiftPattern,
		TCounts,
        TPreprocessing>											TVerifier;

	// MATE-PAIR FILTRATION
	typedef Triple<__int64, TAlignedRead, TAlignQuality>	TDequeueValue;
	typedef Dequeue<TDequeueValue>							TDequeue;
	typedef typename TDequeue::TIter						TDequeueIterator;

    // buffer variable to extend window in previous left hits string by.
    // TODO(holtgrew): DELTA has to be set to a better value, probably.
    const unsigned DELTA = 2000;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE

    // Allocate fragmentstore for job's results.
    TFragmentStore * localStore = new TFragmentStore();

    // Thread-wide offset for reads.
    unsigned threadIdOffset = splitters[job.threadId];

    // Initialize verifiers.
    tls.verifierL.store = localStore;
    tls.verifierL.options = job.options;
    tls.verifierL.swiftPattern = job.swiftPatternL;
    tls.verifierL.cnts = 0;

    tls.verifierR.store = localStore;
    tls.verifierR.options = job.options;
    tls.verifierR.swiftPattern = job.swiftPatternR;
    tls.verifierR.cnts = 0;

	const unsigned NOT_VERIFIED = 1u << (8 * sizeof(unsigned) - 1);

    TOptions & options = tls.options;

    TSwiftPattern & swiftPatternL = *job.swiftPatternL;
    TSwiftPattern & swiftPatternR = *job.swiftPatternR;
    (void)swiftPatternR;
    // TSwiftFinderL & swiftFinderL = *job.swiftFinderL;
    // TSwiftFinderR & swiftFinderR = *job.swiftFinderR;
    TReadSet	&readSetL = indexText(host(*job.swiftPatternL));
    TReadSet	&readSetR = indexText(host(*job.swiftPatternR));
    TVerifier	&verifierL = tls.verifierL;
    TVerifier	&verifierR = tls.verifierR;

	// distance <= libLen + libErr + 2*(parWidth-readLen) - shapeLen
	// distance >= libLen - libErr - 2*parWidth + shapeLen
	TSize readLength = length(readSetL[0]);
	TSignedGPos maxDistance = options.libraryLength + options.libraryError - 2 * (int)readLength - (int)length(indexShape(host(swiftPatternL)));
	TSignedGPos minDistance = options.libraryLength - options.libraryError + (int)length(indexShape(host(swiftPatternL)));
	TGPos scanShift = (minDistance < 0) ? 0 : minDistance;

    Pair<TGPos> gPair;

    // Make sure the queue is empty and all entries in the last potential
    // match no map point before the first no in this queue.
    clear(tls.fifo);
    tls.fifoFirstNo = 0;
    tls.fifoLastNo = 0;
    // TODO(holtgrew): Could do length(...)/job.stride but would have to do so everywhere else below, too.
    // TODO(holtgrew): Get around the clear() and resize-with-fill somehow.
    clear(tls.fifoLastPotMatchNo);
    resize(tls.fifoLastPotMatchNo, length(indexText(host(swiftPatternL))), (__int64)-2);
    
	TGenome & genome = tls.globalStore->contigStore[job.contigId].seq;

    TSize gLength = length(genome);
    
    TAlignedRead mR;
    TAlignQuality qR;
    TDequeueValue fL(-2, mR, qR);	// to supress uninitialized warnings

    //	unsigned const preFetchMatches = 2048;

    // -----------------------------------------------------------------------
    // Enqueue hits from the previous left window.
    // -----------------------------------------------------------------------
    // if (job.prevHitsPtrL.get() == 0)
    //     std::cerr << "\nLENGTH NO PREV HITS" << std::endl;
    // else
    //     std::cerr << "\nLENGTH length(*job.prevHitsPtrL) == " << length(*job.prevHitsPtrL) << std::endl;
    if (job.prevHitsPtrL.get() != 0 && job.prevHitsLEnd > job.prevHitsLBegin) {
        THitStringIter it = iter(*job.prevHitsPtrL, job.prevHitsLEnd, Standard());
        THitStringIter itBegin = iter(*job.prevHitsPtrL, job.prevHitsLBegin, Standard());

        do {
            --it;
            if (it->hstkPos + maxDistance + DELTA < job.rightWindowBegin)
                break;
        } while (it != itBegin);

        // it now points either to beginning or left of the actual entry we
        // want to point to.  Move it right if it is not actually left of the entry.
        if (it->hstkPos + maxDistance + DELTA == job.rightWindowBegin)  // TODO(holtgrew): Keep in sync with above expression, except for </==
            ++it;

        // Now enqueue all hits in the overlap that belong to this job.
        for (THitStringIter itEnd = iter(*job.prevHitsPtrL, job.prevHitsLEnd, Standard()); it != itEnd; ++it) {
            // std::cerr << " [from left window]" << std::flush;
            gPair = Pair<TGPos, TGPos>(_max(static_cast<TSignedGPos>(0), static_cast<TSignedGPos>(it->hstkPos)), _min(it->hstkPos + it->bucketWidth, static_cast<TSignedGPos>(length(genome))));

            fL.i1 = tls.fifoLastPotMatchNo[it->ndlSeqNo];
            tls.fifoLastPotMatchNo[it->ndlSeqNo] = tls.fifoLastNo++;
					
            fL.i2.readId = tls.globalStore->matePairStore[it->ndlSeqNo].readId[0] | NOT_VERIFIED;
            fL.i2.beginPos = gPair.i1;
            fL.i2.endPos = gPair.i2;

            // std::cerr << "\nPREV SWIFT\tL\t" << tls.globalStore->matePairStore[it->ndlSeqNo].readId[0] << "\t" << tls.globalStore->readNameStore[tls.globalStore->matePairStore[it->ndlSeqNo].readId[0]] << "\t" << fL.i2.beginPos << "\t" << fL.i2.endPos << std::endl;
					
            pushBack(tls.fifo, fL);
        }
    }

    // -----------------------------------------------------------------------
    // Perform paired-end verification.
    // -----------------------------------------------------------------------
    THitStringIter itL = iter(*job.hitsPtrL, job.hitsLBegin, Standard());
    THitStringIter itEndL = iter(*job.hitsPtrL, job.hitsLEnd, Standard());
    // Iterate over all filtration results are returned by SWIFT.
    for (THitStringIter itR = iter(*job.hitsPtrR, job.hitsRBegin, Standard()), itEndR = iter(*job.hitsPtrR, job.hitsREnd, Standard()); itR != itEndR; ++itR)
    {
        ++options.countFiltration;

#ifdef RAZERS_DEBUG_MATEPAIRS
        std::cerr << "\nSWIFT\tR\t" << itR->ndlSeqNo << "\t" << tls.globalStore->readNameStore[2 * itR->ndlSeqNo + 1] << "\t" << scanShift + itR->hstkPos << "\t" << scanShift + itR->hstkPos + itR->bucketWidth << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS

        unsigned matePairId = itR->ndlSeqNo;
        TGPos rEndPos = itR->hstkPos + itR->bucketWidth + scanShift;
        TGPos doubleParWidth = 2 * itR->bucketWidth;

        // (1) Remove out-of-window left mates from fifo.
        while (!empty(tls.fifo) && (TSignedGPos)front(tls.fifo).i2.endPos + maxDistance + (TSignedGPos)doubleParWidth < (TSignedGPos)rEndPos)
        {
#ifdef RAZERS_DEBUG_MATEPAIRS
            if (front(tls.fifo).i2.readId > length(tls.globalStore->readNameStore))
                std::cerr << "\nPOP\tL\t" << "[bad read]" << "\t" << front(tls.fifo).i2.beginPos << "\t" << front(tls.fifo).i2.endPos << std::endl;
            else
                std::cerr << "\nPOP\tL\t" << tls.globalStore->readNameStore[front(tls.fifo).i2.readId & ~NOT_VERIFIED] << "\t" << front(tls.fifo).i2.beginPos << "\t" << front(tls.fifo).i2.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
            popFront(tls.fifo);
            ++tls.fifoFirstNo;
        }

        // (2) Add within-window left mates to fifo.
        while (empty(tls.fifo) || (TSignedGPos)back(tls.fifo).i2.endPos + minDistance < (TSignedGPos)(rEndPos + doubleParWidth))
        {
            // Get next left-mate hits, if any and go to next.  This
            // corresponds to a find() in the sequential non-window case.
            if (itL != itEndL)
            {
				++options.countFiltration;
#ifdef RAZERS_DEBUG_MATEPAIRS
                std::cerr << "\nSWIFT\tL\t" << itL->ndlSeqNo << "\t" << tls.globalStore->readNameStore[2 * itL->ndlSeqNo] << "\t" << itL->hstkPos << "\t" << itL->hstkPos + itL->bucketWidth << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                gPair = Pair<TGPos, TGPos>(_max(static_cast<TSignedGPos>(0), static_cast<TSignedGPos>(itL->hstkPos)), _min(itL->hstkPos + itL->bucketWidth, static_cast<TSignedGPos>(length(genome))));
                if ((TSignedGPos)gPair.i2 + maxDistance + (TSignedGPos)doubleParWidth >= (TSignedGPos)rEndPos)
                {
                    // link in
                    fL.i1 = tls.fifoLastPotMatchNo[itL->ndlSeqNo];
                    tls.fifoLastPotMatchNo[itL->ndlSeqNo] = tls.fifoLastNo++;
					
                    // Translate from thread-local to global id.
                    fL.i2.readId = (tls.globalStore->matePairStore[threadIdOffset + itL->ndlSeqNo].readId[0]) | NOT_VERIFIED;
                    fL.i2.beginPos = gPair.i1;
                    fL.i2.endPos = gPair.i2;
					
                    pushBack(tls.fifo, fL);
                }

                ++itL;
            } else {
                break;
            }
        }

        int	bestLeftScore = MinValue<int>::VALUE;
        int bestLibSizeError = MaxValue<int>::VALUE;
        TDequeueIterator bestLeft = TDequeueIterator();

        bool rightVerified = false;
        TDequeueIterator it;
        unsigned leftReadId = tls.globalStore->matePairStore[threadIdOffset + matePairId].readId[0];
		__int64 last = (__int64)-1;
		__int64 lastValid = (__int64)-1;
        __int64 i;
        for (i = tls.fifoLastPotMatchNo[matePairId]; tls.fifoFirstNo <= i; i = (*it).i1)
        {
            it = &value(tls.fifo, i - tls.fifoFirstNo);

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
                        std::cerr << "\nVERIFY\tL\t" << matePairId << "\t" << tls.globalStore->readNameStore[2 * matePairId] << "\t" << (TSignedGPos)(*it).i2.beginPos << "\t" << (*it).i2.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                        ++options.countVerification;
						if (matchVerify(verifierL, infix(genome, ((*it).i2.beginPos >= 0)? (TSignedGPos)(*it).i2.beginPos: (TSignedGPos)0, (TSignedGPos)(*it).i2.endPos), 
										matePairId, readSetL, TRazerSMode()))
                        {
#ifdef RAZERS_DEBUG_MATEPAIRS
                            std::cerr << "  YES: " << verifierL.m.beginPos << "\t" << verifierL.m.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS

							verifierL.m.readId = (*it).i2.readId & ~NOT_VERIFIED;		// has been verified positively
							(*it).i2 = verifierL.m;
							(*it).i3 = verifierL.q;
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
						tls.fifoLastPotMatchNo[matePairId] = i;
					else
						value(tls.fifo, lastValid - tls.fifoFirstNo).i1 = i;
				}
				lastValid = i;

				// XXX
				if (!rightVerified)											// here a verfied left match is available
				{
#ifdef RAZERS_DEBUG_MATEPAIRS
					std::cerr << "\nVERIFY\tR\t" << matePairId << "\t" << tls.globalStore->readNameStore[2 * matePairId + 1] << "\t" << itR->hstkPos << "\t" << itR->hstkPos + itR->bucketWidth + scanShift << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                    ++options.countVerification;
					if (matchVerify(verifierR, swiftInfix(*itR, tls.genomeInf), matePairId, readSetR, TRazerSMode())) {
#ifdef RAZERS_DEBUG_MATEPAIRS
						std::cerr << "  YES: " << verifierR.m.beginPos << "\t" << verifierR.m.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
						rightVerified = true;
						mR = verifierR.m;
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
                    int score = (*it).i3.score;
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
				tls.fifoLastPotMatchNo[matePairId] = i;
			else
				value(tls.fifo, lastValid - tls.fifoFirstNo).i1 = i;
		}
		
        // (4) Verify right mate, if left mate matches.
        if (bestLeftScore != MinValue<int>::VALUE)
        {
            //			if (matchVerify(
            //					mR, qR, infix(swiftFinderR),
            //					matePairId, readSetR, forwardPatternsR,
            //					options, TSwiftSpec()))
            // std::cerr << " [verify]" << std::flush;

            // if (matchVerify(verifierR, swiftInfix(*itR, tls.genomeInf),
            //                 matePairId, readSetR, TRazerSMode()))
            // {
                // distance between left mate beginning and right mate end
                // __int64 dist = (__int64)verifierR.m.endPos - (__int64)(*bestLeft).i2.beginPos;
                // if (dist <= options.libraryLength + options.libraryError &&
                //     options.libraryLength <= dist + options.libraryError)
                // {
                    // mR = verifierR.m;
                    qR = verifierR.q;
					
                    fL.i2 = (*bestLeft).i2;
                    fL.i3 = (*bestLeft).i3;
					
                    // transform mate readNo to global readNo
                    TMatePair &mp = tls.globalStore->matePairStore[threadIdOffset + matePairId];
                    fL.i2.readId = mp.readId[0];
                    mR.readId    = mp.readId[1];
					
                    // transform coordinates to the forward strand
                    if (job.orientation == 'F')
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
                    fL.i2.pairMatchId = mR.pairMatchId = options.nextPairMatchId * options.threadCount + tls.threadId;
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
                        fL.i2.id = length(localStore->alignedReadStore);
                        appendValue(localStore->alignedReadStore, fL.i2, Generous());
                        appendValue(localStore->alignQualityStore, fL.i3, Generous());
                        mR.id = length(localStore->alignedReadStore);
                        appendValue(localStore->alignedReadStore, mR, Generous());
                        appendValue(localStore->alignQualityStore, qR, Generous());

#ifdef RAZERS_DEBUG_MATEPAIRS
                        std::cerr << "\nHIT\tL\t" << fL.i2.readId << "\t" << tls.globalStore->readNameStore[fL.i2.readId] << "\t" << fL.i2.beginPos << "\t" << fL.i2.endPos << std::endl;
                        std::cerr << "\nHIT\tR\t" << mR.readId << "\t" << tls.globalStore->readNameStore[mR.readId] << "\t" << mR.beginPos << "\t" << mR.endPos << std::endl;
#endif  // #ifdef RAZERS_DEBUG_MATEPAIRS
                    }
                // }
            }
        // }
    }

    appendToVerificationResults(*job.verificationResults, localStore);

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE
}

template <typename TFragmentStore, typename TSwiftFinderL, typename TSwiftFinderR, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode, typename TPreprocessing>
void
writeBackToLocal(ThreadLocalStorage<MapPairedReads<TFragmentStore, TSwiftFinderL, TSwiftFinderR, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > & tls, String<TFragmentStore *> & verificationHits)
{
    // TODO(holtgrew): The same as single read writeback! Really? Combine?
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_WRITEBACK);
#endif  // #ifdef RAZERS_PROFILE
    if (tls.threadId == 0u && tls.options._debugLevel >= 3)
        fprintf(stderr, "[writeback]");
	typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
    typedef ThreadLocalStorage<MapPairedReads<TFragmentStore, TSwiftFinderL, TSwiftFinderR, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > TThreadLocalStorage;

    // Update IDs and calculate new size so the prefix increment can be used in the loops.
    TAlignedReadStoreSize oldSize = length(tls.store.alignedReadStore);
    TAlignedReadStoreSize sizeSum = oldSize;

    for (unsigned i = 0; i < length(verificationHits); ++i) {
        TFragmentStore * bucket = verificationHits[i];
        for (unsigned j = 0; j < length(bucket->alignedReadStore); ++j) {
            bucket->alignedReadStore[j].id = sizeSum++;
        }
    }

    // Resize the local read stores appropriately.
	resize(tls.store.alignedReadStore, sizeSum, Generous());
	resize(tls.store.alignQualityStore, sizeSum, Generous());

    // Write back all matches from verification to the block local store.
    for (unsigned i = 0; i < length(verificationHits); ++i) {
        TFragmentStore * bucket = verificationHits[i];
        for (unsigned j = 0; j < length(bucket->alignedReadStore); ++j) {
            move(tls.store.alignedReadStore[oldSize], bucket->alignedReadStore[j]);
            move(tls.store.alignQualityStore[oldSize++], bucket->alignQualityStore[j]);
        }
    }

    SEQAN_ASSERT_EQ(length(tls.store.alignedReadStore), length(tls.store.alignQualityStore));
    if (!empty(tls.store.alignedReadStore))
      SEQAN_ASSERT_EQ(length(tls.store.alignedReadStore), back(tls.store.alignedReadStore).id + 1);

    // Possibly compact matches.
    if (length(tls.store.alignedReadStore) > tls.options.compactThresh)
    {
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_COMPACT);
#endif  // #ifdef RAZERS_PROFILE
        typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
        typename Size<TAlignedReadStore>::Type oldSize = length(tls.store.alignedReadStore);

        // TODO(weese): Duplicates are hard to mask in paired-end mode.
        // if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
        //   maskDuplicates(tls.store, TRazerSMode());  // overlapping parallelograms cause duplicates

        SwiftPatternLSetMaxErrorsWrapper<TThreadLocalStorage> wrapperL(tls);
        SwiftPatternRSetMaxErrorsWrapper<TThreadLocalStorage> wrapperR(tls);
        compactPairMatches(*tls.globalStore, tls.store, tls.counts, tls.options, wrapperL, wrapperR);

        if (length(tls.store.alignedReadStore) * 4 > oldSize) {     // the threshold should not be raised if too many matches were removed
            while (tls.options.compactThresh < oldSize)
                tls.options.compactThresh *= 2;  // Using * 2 in parallel version for scalability reasons.
                // tls.options.compactThresh += (tls.options.compactThresh >> 1);  // if too many matches were removed
            if (tls.threadId == 0u && tls.options._debugLevel >= 3)
                fprintf(stderr, "[raising threshold to %u]", unsigned(tls.options.compactThresh));
        }
#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_COMPACT);
#endif  // #ifdef RAZERS_PROFILE
    }

    SEQAN_ASSERT_EQ(length(tls.store.alignedReadStore), length(tls.store.alignQualityStore));
    if (!empty(tls.store.alignedReadStore))
      SEQAN_ASSERT_EQ(length(tls.store.alignedReadStore), back(tls.store.alignedReadStore).id + 1);

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_WRITEBACK);
#endif  // #ifdef RAZERS_PROFILE
}


// Find read matches in one genome sequence.
//
// The parallelization is similar to the single-end mode.  We perform the
// filtering window-wise.  Then, we generate verification jobs from the swift
// hits.  Here, we distribute them to the jobs by a hash function (num % len),
// constructable and reconstructable by the stride (result of modulo) and the
// module number.  These are then processed by all "leading" threads, where
// leading threads are those with the largest number of processed windows.
template <
	typename TFSSpec,
	typename TFSConfig,
    typename TThreadLocalStorages,
	typename TPreprocessing,
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode>
void _mapMatePairReadsParallel(
	FragmentStore<TFSSpec, TFSConfig>		& store,
	unsigned								  contigId,				// ... and its sequence number
    TThreadLocalStorages                    & threadLocalStorages,
    String<unsigned> const & splitters,
	TPreprocessing							& preprocessingL,
	TPreprocessing							& preprocessingR,
	TCounts									& /*cnts*/,  // TODO(holtgrew): What about this?
	char									  orientation,			// q-gram index of reads
	TRazerSOptions							& options,
	TRazerSMode						  const & /*mode*/)
{
	typedef FragmentStore<TFSSpec, TFSConfig>				TFragmentStore;
	typedef typename TFragmentStore::TMatePairStore			TMatePairStore;
	typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
	typedef typename Value<TMatePairStore>::Type			TMatePair;
	typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
	typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    typedef typename TThreadLocalStorage::TReadSet TReadSet;
    typedef typename TThreadLocalStorage::TShape TShape;
	typedef Index<TReadSet, IndexQGram<TShape>	>	TReadIndex;
	typedef typename Id<TAlignedRead>::Type					TId;

	typedef typename TFragmentStore::TContigSeq				TGenome;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Position<TGenome>::Type				TGPos;
	typedef typename MakeSigned_<TGPos>::Type				TSignedGPos;
	typedef typename Infix<TGenome>::Type					TGenomeInf;

    typedef TRazerSOptions TOptions;

	// FILTRATION
	typedef typename TThreadLocalStorage::TSwiftFinderL	TSwiftFinderL;
	typedef typename TThreadLocalStorage::TSwiftFinderR	TSwiftFinderR;
	typedef typename TThreadLocalStorage::TSwiftPattern	TSwiftPattern;

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

	typedef typename TSwiftFinderL::THitString THitString;
    typedef Job<PairedVerification<TFragmentStore, THitString, TOptions, TSwiftPattern> > TVerificationJob;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

    // -----------------------------------------------------------------------
    // Guard against too small contigs.
    // -----------------------------------------------------------------------

	// distance <= libLen + libErr + 2*(parWidth-readLen) - shapeLen
	// distance >= libLen - libErr - 2*parWidth + shapeLen
	// TSize readLength = length(threadLocalStorages[0].readSetL[0]); // XXX
	// TSignedGPos maxDistance = options.libraryLength + options.libraryError - 2 * (int)readLength - (int)length(indexShape(host(threadLocalStorages[0].swiftPatternL))); // XXX
	TSignedGPos minDistance = options.libraryLength - options.libraryError + (int)length(indexShape(host(threadLocalStorages[0].swiftPatternL)));
	TGPos scanShift = (minDistance < 0) ? 0: minDistance;
	
	// exit if contig is shorter than library size
	TGenome &genome = store.contigStore[contigId].seq;
	if (length(genome) <= scanShift)
		return;

    // -----------------------------------------------------------------------
    // Reverse-complement contig if necessary.
    // -----------------------------------------------------------------------
	// lockContig(store, contigId);
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE
	if (orientation == 'R')
        reverseComplement(genome);
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE

    // -----------------------------------------------------------------------
    // Per-contig initialization of thread local storage objects.
    // -----------------------------------------------------------------------
    // TODO(holtgrew): Maybe put into its own function?
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < static_cast<int>(options.threadCount); ++i) {
		// Initialize verifier objects.
		threadLocalStorages[i].verifierL.onReverseComplement = (orientation == 'R');
		threadLocalStorages[i].verifierL.genomeLength = length(genome);
		threadLocalStorages[i].verifierL.oneMatchPerBucket = true;
		threadLocalStorages[i].verifierL.m.contigId = contigId;
		threadLocalStorages[i].verifierL.preprocessing = &preprocessingL;

		threadLocalStorages[i].verifierR.onReverseComplement = (orientation == 'R');
		threadLocalStorages[i].verifierR.genomeLength = length(genome);
		threadLocalStorages[i].verifierR.oneMatchPerBucket = true;
		threadLocalStorages[i].verifierR.m.contigId = contigId;
		threadLocalStorages[i].verifierR.preprocessing = &preprocessingR;

        threadLocalStorages[i].swiftFinderL  = TSwiftFinderL(genome, options.repeatLength, 1);
        threadLocalStorages[i].genomeInf = infix(genome, scanShift, length(genome));
        threadLocalStorages[i].swiftFinderR  = TSwiftFinderR(threadLocalStorages[i].genomeInf, options.repeatLength, 1);
    }

    // -----------------------------------------------------------------------
    // Perform filtration.
    // -----------------------------------------------------------------------
    TaskQueue<TVerificationJob, OmpLock> taskQueue;
    volatile unsigned leaderWindowsDone = 0;  // Number of windows done in leaders.
    volatile unsigned threadsFiltering = options.threadCount;

    #pragma omp parallel
    {
        unsigned windowsDone = 0;

        // Initialization.
        TThreadLocalStorage & tls = threadLocalStorages[omp_get_thread_num()];
        TSwiftPattern & swiftPatternL = tls.swiftPatternL;
        TSwiftPattern & swiftPatternR = tls.swiftPatternR;
        TSwiftFinderL & swiftFinderL = tls.swiftFinderL;
        TSwiftFinderR & swiftFinderR = tls.swiftFinderR;
        TReadSet	&readSetL = tls.readSetL;
        // TReadSet	&readSetR = tls.readSetR;  // XXX
        // TVerifier	&verifierL = tls.verifierL;  // XXX
        // TVerifier	&verifierR = tls.verifierR;  // XXX

        (void)readSetL;
        SEQAN_ASSERT_NOT(empty(readSetL));
        
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE
        if (!windowFindBegin(swiftFinderL, swiftPatternL, tls.options.errorRate))
            std::cerr << "ERROR: windowFindBegin() failed for left reads in thread " << tls.threadId << std::endl;
        if (!windowFindBegin(swiftFinderR, swiftPatternR, tls.options.errorRate))
            std::cerr << "ERROR: windowFindBegin() failed for right reads in thread " << tls.threadId << std::endl;
#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

        // Previous left hits will be stored here.
        std::tr1::shared_ptr<THitString> previousLeftHits;
        std::tr1::shared_ptr<THitString> leftHits;
        std::tr1::shared_ptr<THitString> rightHits;
        // Declare hits splitters and initialize for left since we need it as "previous left" below.
        String<size_t> leftHitsSplitters;
        resize(leftHitsSplitters, options.maxVerificationPackageCount + 1, 0);
        String<size_t> rightHitsSplitters;
        String<size_t> previousLeftHitsSplitters;

        // For each filtration window...
        bool hasMore = false;
        do {
#ifdef RAZERS_PROFILE
            timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

            // TDequeue fifo;						// stores left-mate potential matches  // XXX
            // String<__int64> lastPotMatchNo;		// last number of a left-mate potential // XXX
            // __int64 lastNo = 0;					// last number over all left-mate pot. matches in the queue // XXX
            // __int64 firstNo = 0;				// first number over all left-mate pot. match in the queue // XXX
            // Pair<TGPos> gPair;

            // resize(lastPotMatchNo, length(host(swiftPatternL)), (__int64)-1, Exact());  // XXX

            // TSize gLength = length(genome);  // XXX

            // TAlignedRead mR;
            // TAlignQuality qR;
            // TDequeueValue fL(-1, mR, qR);	// to supress uninitialized warnings

            // Search for hits from next window.
            int delta = windowsDone == 0 ? scanShift : 0;  // First window of right finder is smaller.
            hasMore = windowFindNext(tls.swiftFinderR, tls.swiftPatternR, tls.options.windowSize - delta);
            bool ret = windowFindNext(tls.swiftFinderL, tls.swiftPatternL, tls.options.windowSize);
            (void) ret;
            SEQAN_ASSERT_TRUE(ret == hasMore);
            windowsDone += 1;  // Local windows done count.
            atomicMax(leaderWindowsDone, windowsDone);

            size_t rightWindowBegin = beginPosition(tls.genomeInf) + tls.options.windowSize * (windowsDone - 1);

            // Create verification jobs.
            // std::cerr << "\nSWIFT HIT COUNT\t" << length(getSwiftHits(tls.swiftFinderL)) << "\t" << length(getSwiftHits(tls.swiftFinderR)) << std::endl;
            // if (windowsDone > 1)
            //     std::cerr << "\nPOSITION       \t" << tls.swiftFinderL.curPos << "\t" << tls.swiftFinderR.curPos << std::endl;
            if (length(getSwiftHits(tls.swiftFinderL)) > 0u || length(getSwiftHits(tls.swiftFinderR)) > 0u) {
                String<TVerificationJob> jobs;

                // Update previous left hits and splitters.
                previousLeftHits = leftHits;
                swap(previousLeftHitsSplitters, leftHitsSplitters);
                // Update new left hits and splitters.
                leftHits.reset(new THitString());
                std::swap(*leftHits, getSwiftHits(tls.swiftFinderL));
                buildHitSplittersAndPartitionHits(leftHitsSplitters, *leftHits, options.maxVerificationPackageCount, length(indexText(host(swiftPatternL))));
                rightHits.reset(new THitString());
                std::swap(*rightHits, getSwiftHits(tls.swiftFinderR));
                buildHitSplittersAndPartitionHits(rightHitsSplitters, *rightHits, options.maxVerificationPackageCount, length(host(swiftPatternL)));

                for (unsigned i = 0; i < options.maxVerificationPackageCount; ++i) {
                    // std::cerr << "i == " << i << " length(leftHitsSplitters) == " << length(leftHitsSplitters) << ", length(rightHitsSplitters) == " << length(rightHitsSplitters) << ", length(previousLeftHitsSplitters) == " << length(previousLeftHitsSplitters) << ", verificationPackageCount == " << options.maxVerificationPackageCount << std::endl;
                    // std::cerr << "length(leftHits) == " << length(leftHits) << ", length(rightHits) == " << length(rightHits) << std::endl;
                    SEQAN_ASSERT_LEQ(previousLeftHitsSplitters[i], previousLeftHitsSplitters[i + 1]);
                    SEQAN_ASSERT_LEQ(rightHitsSplitters[i], rightHitsSplitters[i + 1]);
                    SEQAN_ASSERT_LEQ(leftHitsSplitters[i], leftHitsSplitters[i + 1]);
                    if (rightHitsSplitters[i] == rightHitsSplitters[i + 1] || (previousLeftHitsSplitters[i] == previousLeftHitsSplitters[i + 1] && leftHitsSplitters[i] == leftHitsSplitters[i + 1]))
                        continue;
                    appendValue(jobs, TVerificationJob(tls.threadId, tls.verificationResults, store, contigId, orientation, previousLeftHits, previousLeftHitsSplitters[i], previousLeftHitsSplitters[i + 1], leftHits, leftHitsSplitters[i], leftHitsSplitters[i + 1], rightHits, rightHitsSplitters[i], rightHitsSplitters[i + 1], rightWindowBegin, *tls.globalOptions, tls.swiftPatternL, tls.swiftPatternR));
                }

                pushFront(taskQueue, jobs);
            }
#ifdef RAZERS_PROFILE
            timelineEndTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

            // Perform verification as long as we are a leader and there are filtration jobs to perform.
            while (leaderWindowsDone == windowsDone) {
                TVerificationJob job;
                if (!popFront(job, taskQueue))
                    break;
                workVerification(tls, job, splitters);
            }

            // Write back verification results for this thread so far.
            //
            // First, swap out the current set of local stores from the verification results.
            omp_set_lock(&tls.verificationResults.lock->lock_);
            String<TFragmentStore *> localStores;
            std::swap(localStores, tls.verificationResults.localStores);
            omp_unset_lock(&tls.verificationResults.lock->lock_);
            // Write back the contents of these stores to the thread-local store.
            writeBackToLocal(tls, localStores);
            clearLocalStores(localStores);
        } while(hasMore);

        // Finalization
        windowFindEnd(swiftFinderL, swiftPatternL);
        windowFindEnd(swiftFinderR, swiftPatternR);

        #pragma omp atomic
        threadsFiltering -= 1;

        // Continue to try to help verify.
        while (threadsFiltering > 0u) {
            TVerificationJob job;
            if (popFront(job, taskQueue))
                workVerification(tls, job, splitters);
        }

        // After every thread is done with everything, write back once more.
        #pragma omp barrier
        writeBackToLocal(tls, tls.verificationResults.localStores);
        clearLocalStores(tls.verificationResults.localStores);
    }

    // NOTE: We never re-reverse complement since this function is only called
    // twice, in the right order regarding the orientation parameters.
	// if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
	// 	if (orientation == 'R')	reverseComplement(genome);	// we have to restore original orientation
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE
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

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
	typedef RazerSOptions<TSpec> TOptions;
	
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
	typedef Pattern<TRead const, MyersUkkonen>				TMyersPattern;	// verifier

	typedef typename TFragmentStore::TContigSeq						TContigSeq;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >					TSwiftFinderL;
	typedef typename Infix<TContigSeq>::Type					    TContigInf;
	typedef Finder<TContigInf, Swift<TSwiftSpec> >					TSwiftFinderR;

    // -----------------------------------------------------------------------
    // Initialize global information.
    // -----------------------------------------------------------------------

    // Save OpenMP maximal thread count so we can restore it below, then set
    // from options.
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Number of threads:               \t" << options.threadCount << std::endl;
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount);

    // Verifier preprocessing.
#ifdef RAZERS_BANDED_MYERS
	typedef Nothing TPreprocessing;
	TPreprocessing forwardPatternsL;
	TPreprocessing forwardPatternsR;
#else  // #ifdef RAZERS_BANDED_MYERS
    // TODO(holtgrew): Parallelize preprocessing?
    typedef String<TMyersPattern> TPreprocessing;
	TPreprocessing forwardPatternsL;
	TPreprocessing forwardPatternsR;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (options.gapMode == RAZERS_GAPPED)
	{
	    unsigned pairCount = length(store.matePairStore);
		resize(forwardPatternsL, pairCount, Exact());
		resize(forwardPatternsR, pairCount, Exact());
		Dna5String tmp;
		for(unsigned i = 0; i < pairCount; ++i)
		{
			setHost(forwardPatternsL[i], store.readSeqStore[store.matePairStore[i].readId[0]]);
			tmp = store.readSeqStore[store.matePairStore[i].readId[1]];
			reverseComplement(tmp);
			setHost(forwardPatternsR[i], tmp);
			_patternMatchNOfPattern(forwardPatternsL[i], options.matchN);
			_patternMatchNOfPattern(forwardPatternsR[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsL[i], options.matchN);
			_patternMatchNOfFinder(forwardPatternsR[i], options.matchN);
		}
	}
#endif  // #ifdef RAZERS_BANDED_MYERS

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_INIT);
    double beginInit = sysTime();
#endif  // #ifdef RAZERS_PROFILE

    // Clear/initialize global stats.
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;
	
    // -----------------------------------------------------------------------
    // Initialize thread local storages.
    // -----------------------------------------------------------------------
	SEQAN_PROTIMESTART(initTime);
    String<unsigned> splitters;
    computeSplittersBySlotCount(splitters, length(store.matePairStore), options.threadCount);
    typedef ThreadLocalStorage<MapPairedReads<TFragmentStore, TSwiftFinderL, TSwiftFinderR, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > TThreadLocalStorage;
    String<TThreadLocalStorage> threadLocalStorages;
    initializeThreadLocalStoragesPaired(threadLocalStorages, store, splitters, shape, options);

#ifdef RAZERS_PROFILE
    double endInit = sysTime();
    std::cerr << "TIME initialization: " << (endInit - beginInit) << " s";
    timelineEndTask(TASK_INIT);
#endif  // #ifdef RAZERS_PROFILE
	double timeInitialization = SEQAN_PROTIMEDIFF(initTime);
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Initialization took              \t" << timeInitialization << " seconds" << ::std::endl;

    // -----------------------------------------------------------------------
    // Perform parallel mapping.
    // -----------------------------------------------------------------------

    // Save compaction threshold and set global threshold to infinity, so matchVerify does not compact!
    int oldThreshold = options.compactThresh;
    options.compactThresh = MaxValue<unsigned>::VALUE;

	SEQAN_PROTIMESTART(findTime);
	for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId) {
		// lock to prevent releasing and loading the same contig twice
		// (once per _mapSingleReadsToContig call)
		lockContig(store, contigId);
		if (options.forward)
			_mapMatePairReadsParallel(store, contigId, threadLocalStorages, splitters, forwardPatternsL, forwardPatternsR, cnts, 'F', options, mode);
		if (options.reverse)
			_mapMatePairReadsParallel(store, contigId, threadLocalStorages, splitters, forwardPatternsL, forwardPatternsR, cnts, 'R', options, mode);
		unlockAndFreeContig(store, contigId);
	}
#ifdef RAZERS_PROFILE
    double endMapping = sysTime();
    std::cerr << std::endl << "TIME mapping: " << (endMapping - endInit) << " s" << std::endl;
#endif  // #ifdef RAZERS_PROFILE

    // Write back local stores to global stores.
    for (unsigned i = 0; i < length(threadLocalStorages); ++i)
        std::cerr << "thread " << i << " has " << length(threadLocalStorages[i].store.alignedReadStore) << " aligned reads." << std::endl;
    writeBackToGlobalStore(store, threadLocalStorages);
#ifdef RAZERS_PROFILE
    double endWriteback = sysTime();
    std::cerr << "TIME back to global: " << (endWriteback - endMapping) << " s" << std::endl;
#endif  // #ifdef RAZERS_PROFILE

    // -----------------------------------------------------------------------
    // Collect global statistics, cleanup.
    // -----------------------------------------------------------------------

    // Restore old compaction threshold.
    options.compactThresh = oldThreshold;

    // Add up thread-local filtration and verification counters and print totals.
    for (unsigned i = 0; i < length(threadLocalStorages); ++i) {
        options.countFiltration += threadLocalStorages[i].options.countFiltration;
        options.countVerification += threadLocalStorages[i].options.countVerification;
    }

	options.timeMapReads = SEQAN_PROTIMEDIFF(findTime);
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
