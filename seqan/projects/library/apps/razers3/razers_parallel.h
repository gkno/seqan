#ifndef RAZERS_RAZERS_PARALLEL_H_
#define RAZERS_RAZERS_PARALLEL_H_

// TODO(holtgrew): Ideally, we do not need any locks.

// #ifdef PLATFORM_WINDOWS
// #include <memory>
// #else  // #ifdef PLATFORM_WINDOWS
// #include <tr1/memory>
// #endif  // #ifdef PLATFORM_WINDOWS

#include <seqan/parallel.h>

#include "parallel_misc.h"
#include "parallel_job_queue.h"

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

template <typename TSpec>
class Lock;

struct Omp_;
typedef Tag<Omp_> Omp;

template <>
class Lock<Omp>
{
public:
    omp_lock_t lock_;

    Lock() { omp_init_lock(&lock_); }

    ~Lock() { omp_destroy_lock(&lock_); }
};

// Stores the results of the verification.
//
// Put into its own class so it can be locked independently of other class
// members.
template <typename TMatches>
class SingleVerificationResults
{
public:
    String<TMatches *> localMatches;
    Lock<Omp> * lock;

    SingleVerificationResults() : lock(new Lock<Omp>()) {}

    SingleVerificationResults(SingleVerificationResults const & other)
            : localMatches(other.localMatches), lock(new Lock<Omp>())
    {
        // Not thread-safe copying since this is only used at the beginning when resizing block local storages string.
    }

    SingleVerificationResults & operator=(SingleVerificationResults const & other)
    {
        if (this == &other)
            return *this;
        localMatches = other.localMatches;
    }

    ~SingleVerificationResults()
    {
        delete lock;
    }
};

template <typename TMatches, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename TPreprocessing>
struct MapSingleReads {};

/**
.Class.ThreadLocalStorage:Encapsulates the thread-local data.
..signature:ThreadLocalStorage<TSpec>
..param.TSpec:Specialization template argument.
 */
template <typename TSpec>
class ThreadLocalStorage;

// ThreadLocalStorage specialization for single-end read mapping in RazerS.
template <typename TMatches_, typename TFragmentStore, typename TSwiftFinder_, typename TSwiftPattern_, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename TPreprocessing>
class ThreadLocalStorage<MapSingleReads<TMatches_, TFragmentStore, TSwiftFinder_, TSwiftPattern_, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> >
{
public:
    typedef TSwiftPattern_ TSwiftPattern;
    typedef TSwiftFinder_ TSwiftFinder;
    typedef TMatches_ TMatches;
    
    // The id of this thread.
    unsigned threadId;

    // Each thread needs its local options since the compactionThreshold is changed.
    // TODO(holtgrew): Change overall program structure so this is factorized out of the options struct.
    TOptions options;
    TOptions /*const*/ * globalOptions;

    // Each thread has its own SWIFT finder and pattern object.
    TSwiftFinder swiftFinder;
    TSwiftPattern swiftPattern;

    TCounts counts;  // TODO(holtgrew): Artifact?

    TMatches matches;
    TFragmentStore /*const*/ * globalStore;

    TShape shape;

    typedef MatchVerifier<TFragmentStore, TMatches, TOptions, TRazerSMode, TSwiftPattern, TCounts, TPreprocessing> TMatchVerifier;
    TMatchVerifier verifier;

    // Mailbox for the verification results.
    SingleVerificationResults<TMatches> verificationResults;

    String<unsigned> splitters;

    ThreadLocalStorage() {}
};

template <typename TMatches, typename TFragmentStore, typename THitString, typename TOptions, typename TSwiftPattern>
struct SingleVerification;

template <typename TMatches, typename TFragmentStore, typename THitString_, typename TOptions, typename TSwiftPattern>
class Job<SingleVerification<TMatches, TFragmentStore, THitString_, TOptions, TSwiftPattern> >
{
public:
    typedef SingleVerificationResults<TMatches> TVerificationResults;
    typedef THitString_ THitString;

    int threadId;
    TVerificationResults * verificationResults;
    TFragmentStore * globalStore;
    unsigned contigId;
    std::tr1::shared_ptr<THitString> hitsPtr;
    unsigned hitBegin;
    unsigned hitEnd;
    TOptions * options;
    TSwiftPattern * swiftPattern;

    Job() {}
    
    Job(int threadId_, TVerificationResults & verificationResults_, TFragmentStore & globalStore_, unsigned contigId_, std::tr1::shared_ptr<THitString> & hitsPtr_, unsigned hitBegin_, unsigned hitEnd_, TOptions & options_, TSwiftPattern & swiftPattern_)
            : threadId(threadId_), verificationResults(&verificationResults_), globalStore(&globalStore_), contigId(contigId_), hitsPtr(hitsPtr_), hitBegin(hitBegin_), hitEnd(hitEnd_), options(&options_), swiftPattern(&swiftPattern_)
    {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TMatches>
inline
void
appendToVerificationResults(SingleVerificationResults<TMatches> & verificationResults, TMatches * matchesPtr)
{
    omp_set_lock(&verificationResults.lock->lock_);
    appendValue(verificationResults.localMatches, matchesPtr);
    omp_unset_lock(&verificationResults.lock->lock_);
}

// Uses the read ID to find the correct SWIFT pattern of the TLS, and
// the correct local ID within this SWIFT pattern to update the max
// errors.
//
// We do not disable the read right 
template <typename TMatches, typename TFragmentStore, typename TSwiftFinder_, typename TSwiftPattern_, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename TReadNo, typename TMaxErrors, typename TPreprocessing>
void
setMaxErrors(ThreadLocalStorage<MapSingleReads<TMatches, TFragmentStore, TSwiftFinder_, TSwiftPattern_, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > & tls,
             TReadNo readNo,
             TMaxErrors maxErrors)
{
	int localReadNo = readNo - tls.splitters[tls.threadId];

	int minT = _qgramLemma(tls.swiftPattern, localReadNo, maxErrors);
	if (minT > 1){
		if (maxErrors < 0)
            minT = MaxValue<int>::VALUE;
		setMinThreshold(tls.swiftPattern, localReadNo, static_cast<unsigned>(minT));
	}
}

template <typename TMatches, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename THitString, typename TPreprocessing>
void workVerification(ThreadLocalStorage<MapSingleReads<TMatches, TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > & tls,
                      Job<SingleVerification<TMatches, TFragmentStore, THitString, TOptions, TSwiftPattern> > & job,
                      String<unsigned> const & splitters)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE

    typedef typename Iterator<THitString, Standard>::Type THitStringIterator;

    TMatches * localMatches = new TMatches();

    // Initialize verifier.
    tls.verifier.matches = localMatches;
    tls.verifier.options = job.options;
    tls.verifier.swiftPattern = job.swiftPattern;
    tls.verifier.cnts = 0;

    unsigned offset = splitters[job.threadId];
// #ifdef RAZERS_DEBUG
//     std::cerr << "Verifying block " << jobData.blockId << std::endl;
//     std::cerr << "  offset = " << offset << std::endl;
//     std::cerr << "  jobData.matchBeginIndex == " << jobData.matchBeginIndex << std::endl;
// #endif  // #ifdef RAZERS_DEBUG
    for (THitStringIterator it = iter(*job.hitsPtr, job.hitBegin), itEnd = iter(*job.hitsPtr, job.hitEnd); it != itEnd; ++it) 
	{
        if (length(swiftInfix(value(it), job.globalStore->contigStore[job.contigId].seq)) < length(tls.globalStore->readSeqStore[value(it).ndlSeqNo]))
            continue;  // Skip if hit length < read length.  TODO(holtgrew): David has to fix something in banded myers to make this work.

        unsigned absReadId = offset + value(it).ndlSeqNo;
        tls.verifier.m.readId = absReadId;

#ifdef RAZERS_BANDED_MYERS
		tls.verifier.patternState.leftClip = (value(it).hstkPos >= 0)? 0: -value(it).hstkPos;	// left clip if match begins left of the genome
#endif
		matchVerify(tls.verifier, swiftInfix(value(it), job.globalStore->contigStore[job.contigId].seq), absReadId, tls.globalStore->readSeqStore, TRazerSMode());
    }

    appendToVerificationResults(*job.verificationResults, localMatches);

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE
}

template <typename TMatches, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode, typename TPreprocessing>
void
writeBackToLocal(ThreadLocalStorage<MapSingleReads<TMatches, TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > & tls, String<TMatches *> & verificationHits, bool dontCompact)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_WRITEBACK);
#endif  // #ifdef RAZERS_PROFILE
    if (tls.threadId == 0u && tls.options._debugLevel >= 3)
        fprintf(stderr, "[writeback]");
	typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;

    TAlignedReadStoreSize oldSize = length(tls.matches);
    TAlignedReadStoreSize newSize = oldSize;

    for (unsigned i = 0; i < length(verificationHits); ++i)
        newSize += length(*verificationHits[i]);

    reserve(tls.matches, newSize);

    // Write back all matches from verification to the block local store.
    for (unsigned i = 0; i < length(verificationHits); ++i)
        append(tls.matches, *verificationHits[i]);

    // Possibly compact matches.
    if (!dontCompact && length(tls.matches) > tls.options.compactThresh)
    {
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_COMPACT);
#endif  // #ifdef RAZERS_PROFILE
        typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
        typename Size<TAlignedReadStore>::Type oldSize = length(tls.matches);

        // if (tls.threadId == 0u && tls.options._debugLevel >= 3)
        //     fprintf(stderr, "[compact]");
        if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
            maskDuplicates(tls.matches, tls.options, TRazerSMode());  // overlapping parallelograms cause duplicates

        compactMatches(tls.matches, tls.counts, tls.options, TRazerSMode(), tls, COMPACT);

        if (length(tls.matches) * 4 > oldSize) {     // the threshold should not be raised if too many matches were removed
            while (tls.options.compactThresh < oldSize)
                tls.options.compactThresh *= tls.options.compactMult;
            // tls.options.compactThresh += (tls.options.compactThresh >> 1);  // if too many matches were removed
            if (tls.threadId == 0u && tls.options._debugLevel >= 3)
                fprintf(stderr, "[raising threshold to %u]", unsigned(tls.options.compactThresh));
        }
#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_COMPACT);
#endif  // #ifdef RAZERS_PROFILE
    }

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_WRITEBACK);
#endif  // #ifdef RAZERS_PROFILE
}

template <typename TMatches>
void clearLocalMatches(String<TMatches *> & localMatches)
{
  for (unsigned i = 0; i < length(localMatches); ++i)
    delete localMatches[i];
  clear(localMatches);
}

// Find read matches in one genome sequence.
//
// The parallelization is simple.  We perform the filtering window-wise.
// Then, we generate verification jobs from the swift hits.  The SWIFT hits
// are distributed by a simple static load balancing of adjacent hits.  These
// are then processed by all "leading" threads, where leading threads are
// those with the largest number of processed windows.
template <
	typename TFSSpec,
	typename TFSConfig,
	typename TThreadLocalStorages,
    typename TContigId,
	typename TCounts,
	typename TSpec,
	typename TShape,
    typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
    typename TMatchNPolicy,
    typename TPreprocessing >
void _mapSingleReadsParallelToContig(
	FragmentStore<TFSSpec, TFSConfig>						& store,
	TThreadLocalStorages & threadLocalStorages,
    String<unsigned> const & splitters,
    TContigId const                                         & contigId,
	TCounts													& /*cnts*/,
    char                                                      orientation,
	RazerSOptions<TSpec>									& options,
	TShape const											& /*shape*/,
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const      & /*mode*/,
    TPreprocessing                                          & preprocessing)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE
	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TContigSeq						TContigSeq;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type const				TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
    typedef RazerSOptions<TSpec> TOptions;

    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    typedef typename TThreadLocalStorage::TMatches TMatches;

    // TODO(holtgrew): What about cnts, mode?

	typedef typename If<IsSameType<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >					TSwiftFinder;
	typedef Pattern<TIndex, Swift<TSwiftSpec> >						TSwiftPattern;
    typedef RazerSOptions<TSpec> TOptions;

	typedef typename TSwiftFinder::THitString THitString;
    typedef Job<SingleVerification<TMatches, TFragmentStore, THitString, TOptions, TSwiftPattern> > TVerificationJob;

	// Debug output...
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F') ::std::cerr << "[fwd]";
		else                    ::std::cerr << "[rev]";
	}

    // -----------------------------------------------------------------------
    // Reverse-complement contig if necessary.
    // -----------------------------------------------------------------------
	TContigSeq & contigSeq = store.contigStore[contigId].seq;
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE
	if (orientation == 'R')
		reverseComplement(contigSeq);
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_REVCOMP);
#endif  // #ifdef RAZERS_PROFILE

    // -----------------------------------------------------------------------
    // Per-contig initialization of thread local storage objects.
    // -----------------------------------------------------------------------
    // TODO(holtgrew): Maybe put into its own function?
    for (unsigned i = 0; i < options.threadCount; ++i) {
		// Initialize verifier object.
		threadLocalStorages[i].verifier.onReverseComplement = (orientation == 'R');
		threadLocalStorages[i].verifier.genomeLength = length(contigSeq);
		threadLocalStorages[i].verifier.oneMatchPerBucket = false;
		threadLocalStorages[i].verifier.m.contigId = contigId;
		threadLocalStorages[i].verifier.preprocessing = &preprocessing;
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
        tls.swiftFinder = TSwiftFinder(store.contigStore[contigId].seq, tls.options.repeatLength, 1);
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE
        if (!windowFindBegin(tls.swiftFinder, tls.swiftPattern, tls.options.errorRate))
            std::cerr << "ERROR: windowFindBegin() failed in thread " << tls.threadId << std::endl;
#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

        // For each filtration window...
        bool hasMore = false;
        do {
#ifdef RAZERS_PROFILE
            timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE
            // fprintf(stderr, "[filter]");
            hasMore = windowFindNext(tls.swiftFinder, tls.swiftPattern, tls.options.windowSize);

            windowsDone += 1;  // Local windows done count.
            atomicMax(leaderWindowsDone, windowsDone);

            std::tr1::shared_ptr<THitString> hitsPtr(new THitString());	//TODO (weese:) Could we reuse memory here?
            std::swap(*hitsPtr, getSwiftHits(tls.swiftFinder));
            THitString & hits = *hitsPtr;
            tls.options.countFiltration += length(hits);

            // Enqueue verification jobs.
            if (length(hits) > 0u) {
                // Compute splitters, given a verification package size and a
                // bound on the package count.
                String<unsigned> splitters;
                unsigned packageCount = tls.options.maxVerificationPackageCount * omp_get_max_threads();
                computeSplittersBySlotSize(splitters, length(hits), tls.options.verificationPackageSize, packageCount);

                // Push verification jobs to the job queue.
                String<TVerificationJob> jobs;
                reserve(jobs, length(splitters) - 1);
                for (unsigned i = 1; i < length(splitters); ++i)
                    appendValue(jobs, TVerificationJob(tls.threadId, tls.verificationResults, store, contigId, hitsPtr, splitters[i - 1], splitters[i], *tls.globalOptions, tls.swiftPattern));
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
                // fprintf(stderr, "[verify]");
                workVerification(tls, job, splitters);
            }

            // Write back verification results for this thread so far.
            //
            // First, swap out the current set of local stores from the verification results.
            omp_set_lock(&tls.verificationResults.lock->lock_);
            String<TMatches *> localMatches;
            std::swap(localMatches, tls.verificationResults.localMatches);
            omp_unset_lock(&tls.verificationResults.lock->lock_);
            // Don't compact matches if in configured 'block fraction' of genome.
            size_t hstckLen = tls.swiftFinder.endPos - tls.swiftFinder.startPos;
            size_t hstckLeft = tls.swiftFinder.endPos - tls.swiftFinder.curPos;
            double fracTodo = 1.0  * hstckLeft / hstckLen;
            bool dontCompact = tls.options.noCompactFrac >= fracTodo;
            // Write back the contents of these stores to the thread-local store.
            writeBackToLocal(tls, localMatches, dontCompact);
            clearLocalMatches(localMatches);
        } while (hasMore);

        // Finalization
        windowFindEnd(tls.swiftFinder, tls.swiftPattern);

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
        writeBackToLocal(tls, tls.verificationResults.localMatches, true);
        clearLocalMatches(tls.verificationResults.localMatches);
    }

    // NOTE: We never undo the reverse-complementing!
	// if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
    //     if (orientation == 'R')	reverseComplement(contigSeq);	// we have to restore original orientation
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE
}

// Global initialization of block local storages.
template <typename TThreadLocalStorages, typename TFragmentStore, typename TSplitters, typename TShape, typename TOptions>
void initializeThreadLocalStoragesSingle(TThreadLocalStorages & threadLocalStorages,
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

    resize(threadLocalStorages, threadCount);
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < threadCount; ++i) {
        TThreadLocalStorage & tls = threadLocalStorages[i];

        tls.threadId = i;
        tls.globalStore = &store;
        tls.shape = shape;
        tls.options = options;  // TODO(holtgrew): Copy for stats and threshold, really good?
        tls.globalOptions = &options;
        tls.splitters = splitters;

        // Clear pattern and set parameters.
        TSwiftPattern & swiftPattern = tls.swiftPattern;
        clear(swiftPattern);
        swiftPattern.params.minThreshold = options.threshold;
        swiftPattern.params.tabooLength = options.tabooLength;

        // Initialize the index.
        TIndex & index = host(tls.swiftPattern);
        clear(index);
        clear(indexText(index));
        for (TPosition j = splitters[i]; j < splitters[i + 1]; ++j)
            appendValue(indexText(index), store.readSeqStore[j]);
        //unsigned x = length(indexText(index));
        //fprintf(stderr, "Index #%d has %u entries.\n", i, x);
        index.shape = shape;

#ifdef RAZERS_OPENADDRESSING
        index.alpha = options.loadFactor;
#endif
        cargo(index).abundanceCut = options.abundanceCut;
        cargo(index)._debugLevel = options._debugLevel;

        indexRequire(index, QGramSADir());

        tls.swiftPattern.params.printDots = (tls.threadId == 0) && (tls.options._debugLevel > 0);
    }
}

// Write back from thread local storages to global store.
template <typename TFragmentStore,
          typename TThreadLocalStorages>
void
writeBackToGlobalStore(
        TFragmentStore & target,
        TThreadLocalStorages /*const*/ & threadLocalStorages,
        bool isSingleEnd)  // begin/end already swapped for paired-end reads
{
	typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
	typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;
    typedef typename TThreadLocalStorage::TMatches TMatches;
    typedef typename Iterator<TMatches, Standard>::Type TMatchesIterator;

	// Update the IDs and calculate new size so the prefix increment can be
	// used in the loops.
	TAlignedReadStoreSize oldSize = length(target.alignedReadStore);
	TAlignedReadStoreSize newSize = oldSize;

	for (unsigned i = 0; i < length(threadLocalStorages); ++i)
        newSize += length(threadLocalStorages[i].matches);

	// Resize first so copying happens at most once and not every for each
	// block in the worst case
	resize(target.alignedReadStore, newSize, Exact());
	resize(target.alignQualityStore, newSize, Exact());

	// Append single block stores.
	// TODO(holtgrew): Do in parallel!
	for (unsigned i = 0; i < length(threadLocalStorages); ++i) {
        TMatchesIterator it = begin(threadLocalStorages[i].matches, Standard());
        TMatchesIterator itEnd = end(threadLocalStorages[i].matches, Standard());
        for (; it != itEnd; ++it, ++oldSize) {
            if (isSingleEnd && it->orientation == 'R')
                ::std::swap(it->beginPos, it->endPos);
	 		target.alignedReadStore[oldSize] = TAlignedReadStoreElem(oldSize, it->readId, it->contigId, it->beginPos, it->endPos);
            if (!isSingleEnd)
                SEQAN_ASSERT_NEQ(it->pairMatchId, Value<TMatches>::Type::INVALID_ID);
            target.alignedReadStore[oldSize].pairMatchId = it->pairMatchId;
	 		target.alignQualityStore[oldSize] = TAlignedQualStoreElem(it->pairScore, it->score, -it->score);
        }
	}
}

template <typename TFragmentStore,
          typename TThreadLocalStorages>
void
writeBackToGlobalStore(
        TFragmentStore & target,
        TThreadLocalStorages /*const*/ & threadLocalStorages)
{
    writeBackToGlobalStore(target, threadLocalStorages, true);
}

// Performs splitting of reads, initialization of OpenMP and the calls
// mapSingleReadsParallelToContig for each contig.
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
int _mapSingleReadsParallel(
	FragmentStore<TFSSpec, TFSConfig>						& store,
	TCounts													& cnts,
	RazerSOptions<TSpec>									& options,
	TShape const											& shape,
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const      & mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type	const				TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

	typedef typename If<IsSameType<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
	// typedef Pattern<TIndex, Swift<TSwiftSpec> >					TSwiftPattern;
	typedef Pattern<TIndex, Swift<TSwiftSpec> >						TSwiftPattern;
	typedef Pattern<TRead, MyersUkkonen>							TMyersPattern;	// verifier	
	typedef RazerSOptions<TSpec> TOptions;

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
	typedef typename TFragmentStore::TContigSeq						TContigSeq;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >					TSwiftFinder;

    typedef typename TFragmentStore::TContigSeq TContigSeq;
    typedef typename Position<TContigSeq>::Type TContigPos;
    typedef MatchRecord<TContigPos> TMatchRecord;
    typedef String<TMatchRecord> TMatches;

    // -----------------------------------------------------------------------
    // Initialize global information.
    // -----------------------------------------------------------------------

    // Save OpenMP maximal thread count so we can restore it below, then set
    // from options.
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Number of threads:               \t" << options.threadCount << std::endl;
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount);

    // Verifier preprocessing
#ifdef RAZERS_BANDED_MYERS
	typedef Nothing TPreprocessing;
	TPreprocessing preprocessing;
#else
    // TODO(holtgrew): Parallelize preprocessing?
	typedef String<TMyersPattern> TPreprocessing;
    TPreprocessing preprocessing;
	if (options.gapMode == RAZERS_GAPPED) 
	{
		unsigned readCount = countSequences(store.readSeqStore);
		resize(preprocessing, readCount, Exact()); 
		for(unsigned i = 0; i < readCount; ++i) 
		{ 
			setHost(preprocessing[i], store.readSeqStore[i]); 
			_patternMatchNOfPattern(preprocessing[i], options.matchN); 
			_patternMatchNOfFinder(preprocessing[i], options.matchN); 
		} 
	}
#endif

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
    computeSplittersBySlotCount(splitters, length(store.readNameStore), options.threadCount);
    typedef ThreadLocalStorage<MapSingleReads<TMatches, TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode, TPreprocessing> > TThreadLocalStorage;
    String<TThreadLocalStorage> threadLocalStorages;
    initializeThreadLocalStoragesSingle(threadLocalStorages, store, splitters, shape, options);

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

    // For each contig: Map reads in parallel.
	SEQAN_PROTIMESTART(findTime);
    for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId) {
		lockContig(store, contigId);
		if (options.forward)
			_mapSingleReadsParallelToContig(store, threadLocalStorages, splitters, contigId, cnts, 'F', options, shape, mode, preprocessing);
		if (options.reverse)
			_mapSingleReadsParallelToContig(store, threadLocalStorages, splitters, contigId, cnts, 'R', options, shape, mode, preprocessing);
		unlockAndFreeContig(store, contigId);
    }
#ifdef RAZERS_PROFILE
    double endMapping = sysTime();
    std::cerr << std::endl << "TIME mapping: " << (endMapping - endInit) << " s" << std::endl;
#endif  // #ifdef RAZERS_PROFILE

    #pragma omp parallel
    {
        if (IsSameType<TGapMode, RazerSGapped>::VALUE)
            maskDuplicates(threadLocalStorages[omp_get_thread_num()].matches, options, mode);
        Nothing nothing;
        compactMatches(threadLocalStorages[omp_get_thread_num()].matches, cnts, options, mode, nothing, COMPACT_FINAL);
    }
    #pragma omp barrier

    // Write back local stores to global stores.
    writeBackToGlobalStore(store, threadLocalStorages);
#ifdef RAZERS_PROFILE
    double endWriteback = sysTime();
    std::cerr << "TIME back to global: " << (endWriteback - endMapping) << " s" << std::endl;
#endif  // #ifdef RAZERS_PROFILE

    // TODO(holtgrew): Sum up cnts?!

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
		::std::cerr << "Filtration counter:      " << options.countFiltration << ::std::endl;
		::std::cerr << "Successful verfications: " << options.countVerification << ::std::endl;
	}

    // Restore global state.
    omp_set_num_threads(oldMaxThreads);

	return 0;
}

}  // namespace seqan

#endif  // #ifndef RAZERS_RAZERS_PARALLEL_H_
