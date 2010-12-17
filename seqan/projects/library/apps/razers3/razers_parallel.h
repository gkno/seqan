#ifndef RAZERS_RAZERS_PARALLEL_H_
#define RAZERS_RAZERS_PARALLEL_H_

#include <tr1/memory>  // for shared_ptr

#include "parallel_misc.h"
#include "parallel_job_queue.h"

namespace SEQAN_NAMESPACE_MAIN {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

struct RazerSStats
{
    size_t countFiltration;
    size_t countVerification;

    RazerSStats() : countFiltration(0), countVerification(0) {}
};

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
struct MapSingleReads {};

// Global state goes here.
//
// At the moment, this is the free list of local stores, only.
template <typename TFragmentStore>
class GlobalState
{
public:
    std::vector<TFragmentStore *> localStores;
    omp_lock_t lock;

    GlobalState()
    {
        omp_init_lock(&lock);
    }

    ~GlobalState()
    {
        for (typename std::vector<TFragmentStore *>::iterator it = localStores.begin(); it != localStores.end(); ++it)
            delete *it;
        omp_destroy_lock(&lock);
    }

private:
    // No copy constructor, no assignment.
    GlobalState(GlobalState const &) {}
    GlobalState & operator=(GlobalState const &) {}
};

template <typename TFragmentStore>
void allocateStore(TFragmentStore * & result, GlobalState<TFragmentStore> & globalState)
{
    omp_set_lock(&globalState.lock);
    // if (length(globalState.localStores) == 0u) {
    if (globalState.localStores.size() == 0u) {
        omp_unset_lock(&globalState.lock);
        result = new TFragmentStore();
        return;
    }
    result = globalState.localStores.back();
    globalState.localStores.pop_back();
    // result = back(globalState.localStores);
    // eraseBack(globalState.localStores);
    omp_unset_lock(&globalState.lock);
}

template <typename TFragmentStore>
void deallocateStore(GlobalState<TFragmentStore> & globalState, TFragmentStore * const & store)
{
    omp_set_lock(&globalState.lock);
    appendValue(globalState.localStores, store);
    omp_unset_lock(&globalState.lock);
}

template <typename TFragmentStore>
void deallocateStores(GlobalState<TFragmentStore> & globalState, std::vector<TFragmentStore *> & stores)
{
    omp_set_lock(&globalState.lock);
    for (unsigned i = 0; i < length(stores); ++i)
        appendValue(globalState.localStores, stores[i]);
    omp_unset_lock(&globalState.lock);
}

// Simple thread local storage implementation.
template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
class ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> >
{
public:
    // The following members are only accessed through the TLS interface.
    TaskQueue<TJob, OmpLock> jobQueue_;
    volatile bool working_;

    unsigned threadId;
    
    int previousFiltrationJobId;
    TFragmentStore * globalStore;
    TFragmentStore localStore;
    GlobalState<TFragmentStore> * globalState;
    TShape const * shape;
    TOptions const * options;

    TSwiftFinder swiftFinder;
    TSwiftPattern swiftPattern;

    RazerSStats stats;

    ThreadLocalStorage() : previousFiltrationJobId(-1), globalState(0) {}
};

enum RazerSJobType
{
    JOB_NONE,
    JOB_FILTRATION,
    JOB_VERIFICATION
};

// Stores the results of the verification.
//
// Put into its own class so it can be put into a (shared) pointer and still be locked centrally.
template <typename TFragmentStore>
class VerificationResults
{
public:
    omp_lock_t lock_;
    GlobalState<TFragmentStore> * globalState;
    volatile unsigned blocksTotal;
    volatile unsigned blocksDone;
    std::vector<TFragmentStore *> localStores;

    VerificationResults(GlobalState<TFragmentStore> * globalState_)
            : globalState(globalState_), blocksTotal(0), blocksDone(0)
    {
        omp_init_lock(&lock_);
    }

    ~VerificationResults()
    {
        omp_destroy_lock(&lock_);
    }

private:
    VerificationResults(VerificationResults const &) {}
    VerificationResults operator=(VerificationResults const &) {}
};

template <typename TFragmentStore>
inline
void
allocateLocalStore(TFragmentStore * & localStorePtr, VerificationResults<TFragmentStore> & verificationResults)
{
    omp_set_lock(&verificationResults.lock_);
    verificationResults.blocksTotal += 1;
    allocateStore(localStorePtr, *verificationResults.globalState);
    verificationResults.localStores.push_back(localStorePtr);
    // appendValue(verificationResults.localStores, localStorePtr);
    omp_unset_lock(&verificationResults.lock_);
}

template <typename TFragmentStore>
inline
void clear(VerificationResults<TFragmentStore> & verificationResults)
{
    omp_set_lock(&verificationResults.lock_);
    SEQAN_ASSERT_EQ_MSG(verificationResults.blocksTotal, verificationResults.blocksDone, "Can only clear completed verification results.");

    verificationResults.blocksTotal = 0;
    verificationResults.blocksDone = 0;
    deallocateStores(*verificationResults.globalState, verificationResults.localStores);
    omp_unset_lock(&verificationResults.lock_);
}

template <typename TSpec>
class JobData;

template <typename TFragmentStore>
struct Filtration;

// Job data for filtration.
//
// SWIFT results are stored in the thread local storage's SWIFT
// filter's hit string.  The verification results are stored in the
// verification job, too, used for a barrier to wait for all
// verification jobs to complete before putting results into the
// central fragment store and possibly doing compaction there.
template <typename TFragmentStore>
class JobData<Filtration<TFragmentStore> >
{
public:
    typedef typename Size<typename TFragmentStore::TReadSeqStore>::Type TSize;

    // TODO(holtgrew): "active" flag to mark non-stealable
    int jobId;
    unsigned contigId;
    TFragmentStore * fragmentStore;  // global FragmentStore
    TSize readBeginIndex;
    TSize readEndIndex;

    VerificationResults<TFragmentStore> * verificationResults;

    JobData(int jobId_, unsigned contigId_, TFragmentStore *fragmentStore_, GlobalState<TFragmentStore> * globalState_, TSize readBeginIndex_, TSize readEndIndex_)
            : jobId(jobId_), contigId(contigId_), fragmentStore(fragmentStore_), readBeginIndex(readBeginIndex_),
              readEndIndex(readEndIndex_), verificationResults(new VerificationResults<TFragmentStore>(globalState_))
    {}
};

template <typename TFragmentStore, typename TSwiftFinder>
struct Verification;

template <typename TFragmentStore, typename TSwiftFinder>
class JobData<Verification<TFragmentStore, TSwiftFinder> >
{
public:
	typedef typename TSwiftFinder::THitString THitString;
    typedef typename Position<typename TSwiftFinder::THitString>::Type TPosition;

    VerificationResults<TFragmentStore> * verificationResults;
    TFragmentStore * localStore;

    std::tr1::shared_ptr<THitString> & hitString;
    TPosition matchBeginIndex;
    TPosition matchEndIndex;

    // first could be weak
    JobData(VerificationResults<TFragmentStore> * & verificationResults_, TFragmentStore * & localStore_, std::tr1::shared_ptr<THitString> & hitString_, TPosition matchBeginIndex_, TPosition matchEndIndex_)
            : verificationResults(verificationResults_), localStore(localStore_), hitString(hitString_), matchBeginIndex(matchBeginIndex_), matchEndIndex(matchEndIndex_)
    {
    }
};

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
class Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> >
{
public:
    typedef JobData<Filtration<TFragmentStore> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder> > TVerificationJobData;

    RazerSJobType jobType;
    void * jobData;

    Job(Job const & other)
            : jobType(other.jobType), jobData(other.jobData)
    {
        if (jobType == JOB_FILTRATION)
            jobData = new TFiltrationJobData(*static_cast<TFiltrationJobData *>(jobData));
        else if (jobType == JOB_VERIFICATION)
            jobData = new TVerificationJobData(*static_cast<TVerificationJobData *>(other.jobData));
    }
    
    Job() : jobType(JOB_NONE), jobData(0) {}

    Job(TFiltrationJobData * jobData_)
            : jobType(JOB_FILTRATION), jobData(jobData_)
    {}

    Job(TVerificationJobData * jobData_)
            : jobType(JOB_VERIFICATION), jobData(jobData_)
    {}

    Job & operator=(Job const & other)
    {
        if (this == &other)
            return *this;
        if (jobType == JOB_FILTRATION)
            delete static_cast<TFiltrationJobData *>(jobData);
        else if (jobType == JOB_VERIFICATION)
            delete static_cast<TVerificationJobData *>(jobData);

        jobType = other.jobType;
        if (jobType == JOB_FILTRATION)
            jobData = new TFiltrationJobData(*static_cast<TFiltrationJobData *>(other.jobData));
        else if (jobType == JOB_VERIFICATION)
            jobData = new TVerificationJobData(*static_cast<TVerificationJobData *>(other.jobData));
        return *this;
    }

    ~Job()
    {
        if (jobType == JOB_FILTRATION)
            delete static_cast<TFiltrationJobData *>(jobData);
        else if (jobType == JOB_VERIFICATION)
            delete static_cast<TVerificationJobData *>(jobData);
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
struct JobQueue<ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > >
{
    typedef TaskQueue<TJob, OmpLock> Type;
};

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
struct JobQueue<ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > const>
{
    typedef TaskQueue<TJob, OmpLock> const Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
inline
bool
onlyVerificationStealable(Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > const & job)
{
    return job.jobType == JOB_VERIFICATION;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
inline
typename JobQueue<ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > >::Type &
jobQueue(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > & tls)
{
    return tls.jobQueue_;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
inline
bool
isWorking(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > & tls)
{
    return tls.working_;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
inline
void
setWorking(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > & tls, bool b)
{
    tls.working_ = b;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
void
initializeFiltration(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > & tls, JobData<Filtration<TFragmentStore> > const & jobData)
{
    // No need to initialize if we are doing the same filtration job as previously.
    fprintf(stderr, "INITIALIZE: jobData.jobId == %d, tls.previousFiltrationJobId == %d\n", jobData.jobId, tls.previousFiltrationJobId);
    if (jobData.jobId == tls.previousFiltrationJobId)
        return;
    tls.previousFiltrationJobId = jobData.jobId;

    // 1) Initialize Finder.
	tls.swiftFinder = TSwiftFinder(tls.globalStore->contigStore[jobData.contigId].seq, tls.options->repeatLength, 1);

    // 2) Initialize Pattern.
    typedef typename Host<TSwiftPattern>::Type TIndex;
    typedef typename Position<typename TFragmentStore::TContigStore>::Type TPosition;

    // Clear pattern and set parameters.
    clear(tls.swiftPattern);
    tls.swiftPattern.params.minThreshold = tls.options->threshold;
    tls.swiftPattern.params.tabooLength = tls.options->tabooLength;
    tls.swiftPattern.params.printDots = (tls.threadId == 0) && (tls.options->_debugLevel > 0);

    // Re-initialize the index.
    TIndex & index = host(tls.swiftPattern);
    clear(index);
    for (TPosition i = jobData.readBeginIndex; i < jobData.readEndIndex; ++i)
        appendValue(value(index.text), tls.globalStore->readSeqStore[i]);
    index.shape = *tls.shape;
    // TODO(holtgrew): Rebuild index?!

    // 3) Call windowFindBegin.
    bool res = windowFindBegin(tls.swiftFinder, tls.swiftPattern, tls.options->errorRate);
    if (!res)
        std::cerr << "ERROR: windowFindBegin() failed in thread " << tls.threadId << std::endl;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions>
void
workFiltration(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > & tls, JobData<Filtration<TFragmentStore> > & jobData)
{
    fprintf(stderr, "Filtration\n");
    typedef JobData<Filtration<TFragmentStore> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder> > TVerificationJobData;

	typedef typename TSwiftFinder::THitString THitString;

    // First, wait for all verification jobs to complete.

    // Make sure that the SWIFT filter is correctly initialized.
    initializeFiltration(tls, jobData);

    // Filter reads in the next window.
    bool referenceLeft = windowFindNext(tls.swiftFinder, tls.swiftPattern, tls.options->windowSize);
    THitString & hits = getSwiftHits(tls.swiftFinder);
    tls.stats.countFiltration += length(hits);
    {
        int swiftHits = length(hits);
        fprintf(stderr, "Thread %ud Found %d SWIFT hits!\n", tls.threadId, swiftHits);
    }

    // Enqueue filtration job for remaining sequence, if any.  We simply copy over the data.
    if (referenceLeft)
        pushFront(jobQueue(tls), TJob(new TFiltrationJobData(jobData)));
    // TODO(holtgrew): Job with this data MUST NOT BE STOLEN since we update it below.

    // Enqueue verification jobs, if any.
    if (length(hits) == 0u) return;
    String<unsigned> splitters;
    computeSplittersBySlotSize(splitters, length(hits), tls.options->verificationPackageSize);

    String<TJob> jobs;
    reserve(jobs, length(splitters) - 1);
    std::tr1::shared_ptr<THitString> swiftHitsPtr(new THitString());
    std::swap(*swiftHitsPtr, hits);
    for (unsigned i = 1; i < length(splitters); ++i) {
        // Each verification job gets a local FragmentStore from the free list in the global state.
        TFragmentStore * localStorePtr;
        allocateLocalStore(localStorePtr, *jobData.verificationResults);
        appendValue(jobs, TJob(new TVerificationJobData(jobData.verificationResults, localStorePtr, swiftHitsPtr, splitters[i - 1], splitters[i])));
    }
    pushFront(jobQueue(tls), jobs);
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TThreadId>
void
workVerification(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > & tls, JobData<Verification<TFragmentStore, TSwiftFinder> > & jobData, TThreadId thisThreadId, TThreadId jobThreadId)
{
    fprintf(stderr, "Verification\n");
    (void) tls;
    (void) jobData;
    (void) thisThreadId;
    (void) jobThreadId;
    // while not done verifying
    //   filter up to X reads
    //   push filtration job for remaining reads
    //   create verification jobs for these X reads
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TThreadId>
void
work(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > & tls, TJob & job, TThreadId thisThreadId, TThreadId jobThreadId)
{
    typedef JobData<Filtration<TFragmentStore> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder> > TVerificationJobData;

    if (job.jobType == JOB_FILTRATION)
        workFiltration(tls, *static_cast<TFiltrationJobData *>(job.jobData));
    else if (job.jobType == JOB_VERIFICATION)
        workVerification(tls, *static_cast<TVerificationJobData *>(job.jobData), thisThreadId, jobThreadId);
    else
        fprintf(stderr, "INVALID JOB TYPE\n");
}

// Create thread local storages for each thread, fill with filtration
// jobs and fire up using the queue system.
template <
	typename TFSSpec, 
	typename TFSConfig,
    typename TContigId,
    typename TSplitter,
	typename TCounts,
	typename TSpec, 
	typename TShape,
    typename TAlignMode,
	typename TGapMode,
	typename TScoreMode >
void _mapSingleReadsParallelToContig(
	FragmentStore<TFSSpec, TFSConfig>						& store,
    GlobalState<FragmentStore<TFSSpec, TFSConfig> >         & globalState,
    TContigId const                                         & contigId,
    String<TSplitter> const                                 & splitters,
	TCounts													& /*cnts*/,
    char                                                      orientation,
	RazerSOptions<TSpec>									& options,
	TShape const											& shape,
	RazerSMode<TAlignMode, TGapMode, TScoreMode> const      & /*mode*/)
{
	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type						TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, Index_QGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

    // TODO(holtgrew): What about cnts, mode?

	typedef typename IF<TYPECMP<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;
    typedef RazerSOptions<TSpec> TOptions;

    typedef Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > TJob;
    typedef ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions> > TThreadLocalStorage;
    typedef String<TThreadLocalStorage> TThreadLocalStorages;
    typedef JobData<Filtration<TFragmentStore> > TFiltrationJobData;

    // Lock contig and possibly reverse-complement it.
	TContigSeq & contigSeq = store.contigStore[contigId].seq;
	if (orientation == 'R')
		reverseComplementInPlace(contigSeq);

    // Create thread local storages and fill with filtration jobs.
    TThreadLocalStorages threadLocalStorages;
    resize(threadLocalStorages, options.threadCount);
    for (unsigned i = 0, k = 0; i < options.threadCount; ++i) {
        threadLocalStorages[i].threadId = i;
        threadLocalStorages[i].globalStore = &store;
        threadLocalStorages[i].shape = &shape;
        threadLocalStorages[i].options = &options;
        threadLocalStorages[i].globalState = &globalState;
        for (unsigned j = 0; j < options.splitFactor; ++j, ++k) {
            TJob job(new TFiltrationJobData(k, contigId, &store, &globalState, splitters[0], splitters[k + 1]));
            pushBack(jobQueue(threadLocalStorages[i]), job);
        }
    }

    // Fire up filtration jobs.
    work(threadLocalStorages, onlyVerificationStealable<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions>, StealOne());

    // Restore orientation of contig
    // TODO(holtgrew): Always reverse-complemented again since outer caller locks!
	if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
        if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation
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
	typename TScoreMode >
int _mapSingleReadsParallel(
	FragmentStore<TFSSpec, TFSConfig>						& store,
	TCounts													& cnts,
	RazerSOptions<TSpec>									& options,
	TShape const											& shape,
	RazerSMode<TAlignMode, TGapMode, TScoreMode> const      & mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type						TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, Index_QGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

    // Store global state we change below so we are side-effect free.
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount);

    // Compute initial load balancing / number of blocks.
    SEQAN_ASSERT_EQ_MSG(options.splitFactor, 1u, "For now...");
    unsigned blockCount = options.threadCount * options.splitFactor;
    String<TSize> splitters;
    computeSplittersBySlotCount(splitters, length(store.readNameStore), blockCount);
	if (options._debugLevel >= 1) {
		::std::cerr << ::std::endl << "Number of threads:               \t" << options.threadCount << std::endl;
		::std::cerr << ::std::endl << "Number of blocks:                \t" << blockCount << std::endl;
	}

    GlobalState<TFragmentStore> globalState;  // Store free list of local stores, for example.
    
    for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId) {
		lockContig(store, contigId);
		if (options.forward)
			_mapSingleReadsParallelToContig(store, globalState, contigId, splitters, cnts, 'F', options, shape, mode);
		if (options.reverse)
			_mapSingleReadsParallelToContig(store, globalState, contigId, splitters, cnts, 'R', options, shape, mode);
		unlockAndFreeContig(store, contigId);
    }

    // Restore global state.
    omp_set_num_threads(oldMaxThreads);

	return 0;
}

}  // namespace seqan

#endif  // #ifndef RAZERS_RAZERS_PARALLEL_H_
