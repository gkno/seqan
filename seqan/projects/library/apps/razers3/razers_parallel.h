#ifndef RAZERS_RAZERS_PARALLEL_H_
#define RAZERS_RAZERS_PARALLEL_H_

// TODO(holtgrew): Ideally, we do not need any locks.

#include <tr1/memory>

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
template <typename TFragmentStore>
class VerificationResults
{
public:
    String<TFragmentStore *> localStores;
    Lock<Omp> * lock;

    VerificationResults() : lock(new Lock<Omp>()) {}

    VerificationResults(VerificationResults const & other)
            : localStores(other.localStores), lock(new Lock<Omp>())
    {
        // Not thread-safe copying since this is only used at the beginning when resizing block local storages string.
    }

    VerificationResults & operator=(VerificationResults const & other)
    {
        if (this == &other)
            return *this;
        localStores = other.localStores;
    }

    ~VerificationResults()
    {
        delete lock;
    }
};

/**
.Class.ThreadLocalStorage:Encapsulates the thread-local data.
..signature:ThreadLocalStorage<TSpec>
..param.TSpec:Specialization template argument.
 */
template <typename TSpec>
class ThreadLocalStorage;

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode>
struct MapSingleReads {};

// ThreadLocalStorage specialization for single-end read mapping in RazerS.
template <typename TFragmentStore, typename TSwiftFinder_, typename TSwiftPattern_, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode>
class ThreadLocalStorage<MapSingleReads<TFragmentStore, TSwiftFinder_, TSwiftPattern_, TShape, TOptions, TCounts, TRazerSMode> >
{
public:
    typedef TSwiftPattern_ TSwiftPattern;
    typedef TSwiftFinder_ TSwiftFinder;
    
    // The id of this thread.
    unsigned threadId;

    // Each thread needs its local options since the compactionThreshold is changed.
    // TODO(holtgrew): Change overall program structure so this is factorized out of the options struct.
    TOptions options;

    // Each thread has its own SWIFT finder and pattern object.
    TSwiftFinder swiftFinder;
    TSwiftPattern swiftPattern;

    TCounts counts;  // TODO(holtgrew): Artifact?

    TFragmentStore store;
    TFragmentStore /*const*/ * globalStore;

    TShape shape;

    typedef MatchVerifier<TFragmentStore, TOptions, TRazerSMode, TSwiftPattern, TCounts> TMatchVerifier;
    TMatchVerifier verifier;

    // Mailbox for the verification results.
    VerificationResults<TFragmentStore> verificationResults;

    String<unsigned> splitters;

    unsigned completeWindows;

    ThreadLocalStorage() {}
};

template <typename TFragmentStore, typename THitString, typename TOptions, typename TSwiftPattern>
struct Verification;

template <typename TFragmentStore, typename THitString_, typename TOptions, typename TSwiftPattern>
class Job<Verification<TFragmentStore, THitString_, TOptions, TSwiftPattern> >
{
public:
    typedef VerificationResults<TFragmentStore> TVerificationResults;
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

template <typename TFragmentStore>
inline
void
appendToVerificationResults(VerificationResults<TFragmentStore> & verificationResults, TFragmentStore * storePtr)
{
    omp_set_lock(&verificationResults.lock->lock_);
    appendValue(verificationResults.localStores, storePtr);
    omp_unset_lock(&verificationResults.lock->lock_);
}

// Uses the read ID to find the correct SWIFT pattern in the string handled by
// the ParallelSwiftPatternHandler, and the correct local ID within this SWIFT
// pattern to update the max errors.
//
// We do not disable the read right 
template <typename TFragmentStore, typename TSwiftFinder_, typename TSwiftPattern_, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename TReadNo, typename TMaxErrors>
void
setMaxErrors(ThreadLocalStorage<MapSingleReads<TFragmentStore, TSwiftFinder_, TSwiftPattern_, TShape, TOptions, TCounts, TRazerSMode> > & tls,
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

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape/*TODO(holtgrew): Superflous.*/, typename TOptions, typename TCounts, typename TRazerSMode, typename THitString>
void workVerification(ThreadLocalStorage<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls,
                      Job<Verification<TFragmentStore, THitString, TOptions, TSwiftPattern> > & job,
                      String<unsigned> const & splitters)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE

    typedef typename Iterator<THitString, Standard>::Type THitStringIterator;

    TFragmentStore * localStore = new TFragmentStore();

    // Initialize verifier.
    tls.verifier.store = localStore;
    tls.verifier.options = job.options;
    tls.verifier.swiftPattern = job.swiftPattern;
    tls.verifier.cnts = 0;

    unsigned offset = splitters[job.threadId];
// #ifdef RAZERS_DEBUG
//     std::cerr << "Verifying block " << jobData.blockId << std::endl;
//     std::cerr << "  offset = " << offset << std::endl;
//     std::cerr << "  jobData.matchBeginIndex == " << jobData.matchBeginIndex << std::endl;
// #endif  // #ifdef RAZERS_DEBUG
    for (THitStringIterator it = iter(*job.hitsPtr, job.hitBegin), itEnd = iter(*job.hitsPtr, job.hitEnd); it != itEnd; ++it) {
        if (length(swiftInfix(value(it), job.globalStore->contigStore[job.contigId].seq)) < length(tls.globalStore->readSeqStore[value(it).ndlSeqNo]))
            continue;  // Skip if hit length < read length.  TODO(holtgrew): David has to fix something in banded myers to make this work.

        unsigned absReadId = offset + value(it).ndlSeqNo;
        tls.verifier.m.readId = absReadId;

        // TODO(holtgrew): Unsure about global/local indices and global/local store.
		matchVerify(tls.verifier, swiftInfix(value(it), job.globalStore->contigStore[job.contigId].seq), absReadId, tls.globalStore->readSeqStore, TRazerSMode());
    }

    appendToVerificationResults(*job.verificationResults, localStore);

#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_VERIFY);
#endif  // #ifdef RAZERS_PROFILE
}

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
void
writeBackToLocal(ThreadLocalStorage<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, String<TFragmentStore *> & verificationHits)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_WRITEBACK);
#endif  // #ifdef RAZERS_PROFILE
    if (tls.threadId == 0u && tls.options._debugLevel >= 3)
        fprintf(stderr, "[writeback]");
    // TODO(holtgrew): It would be possible to use local sorting and multiway merging, removes necessity to sort in compactMatches/maskDuplicates.
	typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;

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

        if (tls.threadId == 0u && tls.options._debugLevel >= 3)
            fprintf(stderr, "[compact]");
        if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
          maskDuplicates(tls.store, TRazerSMode());  // overlapping parallelograms cause duplicates

        compactMatches(tls.store, tls.counts, tls.options, TRazerSMode(), tls, COMPACT);

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

template <typename TFragmentStore>
void clearLocalStores(String<TFragmentStore *> & localStores)
{
  for (unsigned i = 0; i < length(localStores); ++i)
    delete localStores[i];
  clear(localStores);
}

// Create thread local storages for each thread, fill with filtration
// jobs and fire up using the queue system.
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
    typename TMatchNPolicy >
void _mapSingleReadsParallelToContig(
	FragmentStore<TFSSpec, TFSConfig>						& store,
	TThreadLocalStorages & threadLocalStorages,
    String<unsigned> const & splitters,
    TContigId const                                         & contigId,
	TCounts													& /*cnts*/,
    char                                                      orientation,
	RazerSOptions<TSpec>									& options,
	TShape const											& /*shape*/,
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const      & /*mode*/)
{
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_ON_CONTIG);
#endif  // #ifdef RAZERS_PROFILE
	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type						TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
    typedef RazerSOptions<TSpec> TOptions;

    typedef typename Value<TThreadLocalStorages>::Type TThreadLocalStorage;

    // TODO(holtgrew): What about cnts, mode?

	typedef typename If<IsSameType<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;
    typedef RazerSOptions<TSpec> TOptions;

	typedef typename TSwiftFinder::THitString THitString;
    typedef Job<Verification<TFragmentStore, THitString, TOptions, TSwiftPattern> > TVerificationJob;

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
        threadLocalStorages[i].completeWindows = 0;
		// Initialize verifier object.
		threadLocalStorages[i].verifier.onReverseComplement = (orientation == 'R');
		threadLocalStorages[i].verifier.genomeLength = length(contigSeq);
		threadLocalStorages[i].verifier.oneMatchPerBucket = false;
		threadLocalStorages[i].verifier.m.contigId = contigId;
    }

    // -----------------------------------------------------------------------
    // Perform filtration.
    // -----------------------------------------------------------------------
    TaskQueue<TVerificationJob, OmpLock> taskQueue;
    volatile unsigned leaderWindowsDone = 0;  // Number of windows done in leaders.
    volatile unsigned threadsFiltering = options.threadCount;

    #pragma omp parallel
    {
#ifdef RAZERS_PROFILE
        timelineBeginTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE
        unsigned windowsDone = 0;

        // Initialization.
        TThreadLocalStorage & tls = threadLocalStorages[omp_get_thread_num()];
        tls.swiftFinder = TSwiftFinder(store.contigStore[contigId].seq, tls.options.repeatLength, 1);
        if (!windowFindBegin(tls.swiftFinder, tls.swiftPattern, tls.options.errorRate))
            std::cerr << "ERROR: windowFindBegin() failed in thread " << tls.threadId << std::endl;

        // For each filtration window...
        bool hasMore = false;
        do {
            hasMore = windowFindNext(tls.swiftFinder, tls.swiftPattern, tls.options.windowSize);

            windowsDone += 1;  // Local windows done count.
            // Update global window count for leaders, basically this is a
            // leaderWindowsDone "windowsDone = max(windowsDone, leaderWindowsDone)".
            {
                #pragma omp flush(leaderWindowsDone)
                unsigned unsafeDone = leaderWindowsDone;
                while (unsafeDone < windowsDone) {
                    // TODO(holtgrew): GCC intrinsic is non-portable.
                    __sync_val_compare_and_swap(&leaderWindowsDone, unsafeDone, windowsDone);
                    unsafeDone = leaderWindowsDone;
                    #pragma omp flush(leaderWindowsDone)
                }
            }

            std::tr1::shared_ptr<THitString> hitsPtr(new THitString());
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
                    appendValue(jobs, TVerificationJob(tls.threadId, tls.verificationResults, store, contigId, hitsPtr, splitters[i - 1], splitters[i], tls.options, tls.swiftPattern));
                pushFront(taskQueue, jobs);
            }

            // Perform verification as long as we are a leader and there are filtration jobs to perform.
            #pragma omp flush(leaderWindowsDone)
            while (leaderWindowsDone == windowsDone) {
                TVerificationJob job;
                if (!popFront(job, taskQueue))
                    break;
                workVerification(tls, job, splitters);
                #pragma omp flush(leaderWindowsDone)
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
        } while (hasMore);

        // Finalization
        windowFindEnd(tls.swiftFinder, tls.swiftPattern);

#ifdef RAZERS_PROFILE
        timelineEndTask(TASK_FILTER);
#endif  // #ifdef RAZERS_PROFILE

        #pragma omp atomic
        threadsFiltering -= 1;

        // Continue to try to help verify.
        while (threadsFiltering > 0u) {
            TVerificationJob job;
            if (popFront(job, taskQueue))
                workVerification(tls, job, splitters);
            #pragma omp flush(threadsFiltering)  // TODO(holtgrew): Look at OpenMP memory model and understand whether we have to do this here.
        }

        // After every thread is done with everything, write back once more.
        #pragma omp barrier
        writeBackToLocal(tls, tls.verificationResults.localStores);
        clearLocalStores(tls.verificationResults.localStores);
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
void initializeThreadLocalStorages(TThreadLocalStorages & threadLocalStorages,
                                   TFragmentStore /*const*/ & store,
                                   TSplitters const & splitters,
                                   TShape /*const*/ & shape,
                                   TOptions const & options)
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
        tls.globalStore = & store;
        tls.shape = shape;
        tls.options = options;  // TODO(holtgrew): Copy for stats and threshold, really good?
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
        TPosition jBegin = splitters[i];
        TPosition jEnd = splitters[i + 1];
        for (TPosition j = jBegin; j < jEnd; ++j)
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
        TThreadLocalStorages /*const*/ & threadLocalStorages)
{
	typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
	typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;

	// Update the IDs and calculate new size so the prefix increment can be
	// used in the loops.
	TAlignedReadStoreSize oldSize = length(target.alignedReadStore);
	TAlignedReadStoreSize sizeSum = oldSize;

	for (unsigned i = 0; i < length(threadLocalStorages); ++i)
		for (unsigned j = 0; j < length(threadLocalStorages[i].store.alignedReadStore); ++j)
			threadLocalStorages[i].store.alignedReadStore[j].id = sizeSum++;

	// Resize first so copying happens at most once and not every for each
	// block in the worst case
	resize(target.alignedReadStore, sizeSum, Exact());
	resize(target.alignQualityStore, sizeSum, Exact());

	// Append single block stores.
	// TODO(holtgrew): Do in parallel!
	for (unsigned i = 0; i < length(threadLocalStorages); ++i) {
		for (unsigned j = 0; j < length(threadLocalStorages[i].store.alignedReadStore); ++j) {
			move(target.alignedReadStore[oldSize], threadLocalStorages[i].store.alignedReadStore[j]);
			move(target.alignQualityStore[oldSize++], threadLocalStorages[i].store.alignQualityStore[j]);
		}
	}
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
	typedef typename Value<TReadSeqStore>::Type						TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

	typedef typename If<IsSameType<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
	// typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;
	typedef RazerSOptions<TSpec> TOptions;

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;

    typedef ThreadLocalStorage<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > TThreadLocalStorage;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_INIT);
    double beginInit = sysTime();
#endif  // #ifdef RAZERS_PROFILE

    // Save OpenMP maximal thread count so we can restore it below, then set
    // from options.
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Number of threads:               \t" << options.threadCount << std::endl;
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount);

    // Clear/initialize global stats.
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

    // -----------------------------------------------------------------------
    // Initialize thread local storages.
    // -----------------------------------------------------------------------
    String<unsigned> splitters;
    computeSplittersBySlotCount(splitters, length(store.readNameStore), options.threadCount);
    String<TThreadLocalStorage> threadLocalStorages;
    initializeThreadLocalStorages(threadLocalStorages, store, splitters, shape, options);

    // Save compaction threshold and set global threshold to infinity, so matchVerify does not compact!
    int oldThreshold = options.compactThresh;
    options.compactThresh = MaxValue<unsigned>::VALUE;

#ifdef RAZERS_PROFILE
    double endInit = sysTime();
    std::cerr << "TIME initialization: " << (endInit - beginInit) << " s";
    timelineEndTask(TASK_INIT);
#endif  // #ifdef RAZERS_PROFILE

    // For each contig: Map reads in parallel.
    for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId) {
		lockContig(store, contigId);
		if (options.forward)
			_mapSingleReadsParallelToContig(store, threadLocalStorages, splitters, contigId, cnts, 'F', options, shape, mode);
		if (options.reverse)
			_mapSingleReadsParallelToContig(store, threadLocalStorages, splitters, contigId, cnts, 'R', options, shape, mode);
		unlockAndFreeContig(store, contigId);
    }
#ifdef RAZERS_PROFILE
    double endMapping = sysTime();
    std::cerr << std::endl << "TIME mapping: " << (endMapping - endInit) << " s" << std::endl;
#endif  // #ifdef RAZERS_PROFILE

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
