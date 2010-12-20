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

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
struct MapSingleReads {};

// Global state goes here.
//
// At the moment, this is the free list of local stores, only.
template <typename TFragmentStore, typename TOptions>
class GlobalState
{
public:
    String<unsigned> splitters;
    String<TFragmentStore *> localStores;
    TOptions * globalOptions;
    omp_lock_t lock;

    GlobalState()
      : globalOptions(0)
    {
        omp_init_lock(&lock);
    }

    ~GlobalState()
    {
        typedef typename Iterator<String<TFragmentStore *> >::Type TIterator;
        for (TIterator it = begin(localStores); it != end(localStores); ++it)
          delete *it;
        omp_destroy_lock(&lock);
    }

private:
    // No copy constructor, no assignment.
    GlobalState(GlobalState const &) {}
    GlobalState & operator=(GlobalState const &) {}
};

template <typename TFragmentStore, typename TOptions>
void allocateStore(TFragmentStore * & result, GlobalState<TFragmentStore, TOptions> & globalState)
{
    omp_set_lock(&globalState.lock);
    // if (length(globalState.localStores) == 0u) {
    if (length(globalState.localStores) == 0u) {
        omp_unset_lock(&globalState.lock);
        result = new TFragmentStore();
        return;
    }
    result = back(globalState.localStores);
    clear(result->alignedReadStore);
    clear(result->alignQualityStore);
    eraseBack(globalState.localStores);
    // result = back(globalState.localStores);
    // eraseBack(globalState.localStores);
    omp_unset_lock(&globalState.lock);
}

template <typename TFragmentStore, typename TOptions>
void deallocateStore(GlobalState<TFragmentStore, TOptions> & globalState, TFragmentStore * const & store)
{
    omp_set_lock(&globalState.lock);
    appendValue(globalState.localStores, store);
    omp_unset_lock(&globalState.lock);
}

template <typename TFragmentStore, typename TOptions>
void deallocateStores(GlobalState<TFragmentStore, TOptions> & globalState, String<TFragmentStore *> & stores)
{
    omp_set_lock(&globalState.lock);
    for (unsigned i = 0; i < length(stores); ++i)
        appendValue(globalState.localStores, stores[i]);
    omp_unset_lock(&globalState.lock);
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
struct ParallelSwiftPatternHandler;

// Simple thread local storage implementation.
template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
class ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> >
{
public:
    // The following two members are only accessed through the TLS interface.
    TaskQueue<TJob, OmpLock> jobQueue_;
    volatile bool working_;

    unsigned threadId;
    
    int previousFiltrationJobId;
    TFragmentStore * globalStore;
    TFragmentStore localStore;
    GlobalState<TFragmentStore, TOptions> * globalState;
    TShape const * shape;
    TOptions options;

    TSwiftFinder swiftFinder;
    TSwiftPattern swiftPattern;

    RazerSStats stats;
    TCounts counts;  // TODO(holtgrew): Artifact? ?!?!?

    typedef ParallelSwiftPatternHandler<TJob, TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode>  TParallelSwiftPatternHandler;
    TParallelSwiftPatternHandler swiftPatternHandler;

    typedef MatchVerifier<TFragmentStore, TOptions, TRazerSMode, TParallelSwiftPatternHandler, TCounts> TMatchVerifier;
    TMatchVerifier verifier;

    ThreadLocalStorage() : previousFiltrationJobId(-1), globalState(0), swiftPatternHandler(*this) {}
};

// Pattern handler: is necessary to be able to calculate the overall read ID based on the local 
// read ID within one thread. Used in matchVerify.
template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
struct ParallelSwiftPatternHandler
{
    typedef ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> >  TThreadLocalStorage;
    TThreadLocalStorage & threadLocalStorage;

	ParallelSwiftPatternHandler(TThreadLocalStorage & _threadLocalStorage)
	    : threadLocalStorage(_threadLocalStorage) {}
};

// Uses the read ID to find the correct SWIFT pattern in the string handled by the ParallelSwiftPatternHandler,
// and the correct local ID within this SWIFT pattern to update the max errors.
template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode, typename TReadNo, typename TMaxErrors>
inline 
void 
setMaxErrors(ParallelSwiftPatternHandler<TJob, TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> & handler, TReadNo readNo, TMaxErrors maxErrors)
{
    
	int blockSize = length(host(host(handler.threadLocalStorage.swiftPattern)));
	int localReadNo = readNo % blockSize;

	int minT = _qgramLemma(handler.threadLocalStorage.swiftPattern, localReadNo, maxErrors);
	if (minT > 1){
		if (maxErrors < 0) minT = SupremumValue<int>::VALUE;
		setMinThreshold(handler.threadLocalStorage.swiftPattern, localReadNo, (unsigned)minT);
	}
}

enum RazerSJobType
{
    JOB_NONE,
    JOB_FILTRATION,
    JOB_VERIFICATION
};

// Stores the results of the verification.
//
// Put into its own class so it can be put into a (shared) pointer and still be locked centrally.
template <typename TFragmentStore, typename TOptions>
class VerificationResults
{
public:
    omp_lock_t lock_;
    GlobalState<TFragmentStore, TOptions> * globalState;
    volatile unsigned blocksTotal;
    volatile unsigned blocksDone;
    String<TFragmentStore *> localStores;

    VerificationResults(GlobalState<TFragmentStore, TOptions> * globalState_)
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

template <typename TFragmentStore, typename TOptions>
inline
void
allocateLocalStore(TFragmentStore * & localStorePtr, VerificationResults<TFragmentStore, TOptions> & verificationResults)
{
    omp_set_lock(&verificationResults.lock_);
    verificationResults.blocksTotal += 1;
    allocateStore(localStorePtr, *verificationResults.globalState);
    appendValue(verificationResults.localStores, localStorePtr);
    // appendValue(verificationResults.localStores, localStorePtr);
    omp_unset_lock(&verificationResults.lock_);
}

template <typename TFragmentStore, typename TOptions>
inline
void clear(VerificationResults<TFragmentStore, TOptions> & verificationResults)
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

template <typename TFragmentStore, typename TOptions>
struct Filtration;

// Job data for filtration.
//
// SWIFT results are stored in the thread local storage's SWIFT
// filter's hit string.  The verification results are stored in the
// verification job, too, used for a barrier to wait for all
// verification jobs to complete before putting results into the
// central fragment store and possibly doing compaction there.
template <typename TFragmentStore, typename TOptions>
class JobData<Filtration<TFragmentStore, TOptions> >
{
public:
    typedef typename Size<typename TFragmentStore::TReadSeqStore>::Type TSize;

    // TODO(holtgrew): "active" flag to mark non-stealable
    int jobId;
    unsigned contigId;
    unsigned blockId;
    TFragmentStore * fragmentStore;  // global FragmentStore
    VerificationResults<TFragmentStore, TOptions> * verificationResults;

    JobData(int jobId_, unsigned contigId_, unsigned blockId_, TFragmentStore *fragmentStore_, GlobalState<TFragmentStore, TOptions> * globalState_)
            : jobId(jobId_), contigId(contigId_), blockId(blockId_), fragmentStore(fragmentStore_), verificationResults(new VerificationResults<TFragmentStore, TOptions>(globalState_))
    {}
};

template <typename TFragmentStore, typename TSwiftFinder, typename TOptions>
struct Verification;

template <typename TFragmentStore, typename TSwiftFinder, typename TOptions>
class JobData<Verification<TFragmentStore, TSwiftFinder, TOptions> >
{
public:
	typedef typename TSwiftFinder::THitString THitString;
    typedef typename Position<typename TSwiftFinder::THitString>::Type TPosition;

    TFragmentStore * fragmentStore;
    unsigned contigId;
    unsigned blockId;
    
    VerificationResults<TFragmentStore, TOptions> * verificationResults;
    TFragmentStore * localStore;

    std::tr1::shared_ptr<THitString> hitString;
    TPosition matchBeginIndex;
    TPosition matchEndIndex;

    JobData(TFragmentStore * fragmentStore_, unsigned contigId_, unsigned blockId_, VerificationResults<TFragmentStore, TOptions> * & verificationResults_, TFragmentStore * & localStore_, std::tr1::shared_ptr<THitString> & hitString_, TPosition matchBeginIndex_, TPosition matchEndIndex_)
            : fragmentStore(fragmentStore_), contigId(contigId_), blockId(blockId_), verificationResults(verificationResults_), localStore(localStore_), hitString(hitString_), matchBeginIndex(matchBeginIndex_), matchEndIndex(matchEndIndex_)
    {
    }
};

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
class Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> >
{
public:
    typedef JobData<Filtration<TFragmentStore, TOptions> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder, TOptions> > TVerificationJobData;

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

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
struct JobQueue<ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > >
{
    typedef TaskQueue<TJob, OmpLock> Type;
};

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
struct JobQueue<ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > const>
{
    typedef TaskQueue<TJob, OmpLock> const Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
inline
bool
onlyVerificationStealable(Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > const & job)
{
    return job.jobType == JOB_VERIFICATION;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
inline
typename JobQueue<ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > >::Type &
jobQueue(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls)
{
    return tls.jobQueue_;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
inline
bool
isWorking(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls)
{
    return tls.working_;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
inline
void
setWorking(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, bool b)
{
    tls.working_ = b;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
void
initializeFiltration(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, JobData<Filtration<TFragmentStore, TOptions> > const & jobData)
{
    // TODO(holtgrew): Keep one swift filter per block!
    // No need to initialize if we are doing the same filtration job as previously.
    // fprintf(stderr, "INITIALIZE: jobData.jobId == %d, tls.previousFiltrationJobId == %d\n", jobData.jobId, tls.previousFiltrationJobId);
    if (jobData.jobId == tls.previousFiltrationJobId)
        return;
    tls.previousFiltrationJobId = jobData.jobId;

    // 1) Initialize Finder.
	tls.swiftFinder = TSwiftFinder(tls.globalStore->contigStore[jobData.contigId].seq, tls.options.repeatLength, 1);

    // 2) Initialize Pattern.
    typedef typename Host<TSwiftPattern>::Type TIndex;
    typedef typename Position<typename TFragmentStore::TContigStore>::Type TPosition;

    // Clear pattern and set parameters.
    clear(tls.swiftPattern);
    tls.swiftPattern.params.minThreshold = tls.options.threshold;
    tls.swiftPattern.params.tabooLength = tls.options.tabooLength;
    tls.swiftPattern.params.printDots = (tls.threadId == 0) && (tls.options._debugLevel > 0);

    // Re-initialize the index.
    TIndex & index = host(tls.swiftPattern);
    clear(index);
    clear(indexText(index));
    TPosition iBegin = tls.globalState->splitters[jobData.blockId];
    TPosition iEnd = tls.globalState->splitters[jobData.blockId + 1];
    for (TPosition i = iBegin; i < iEnd; ++i)
        appendValue(value(index.text), tls.globalStore->readSeqStore[i]);
    index.shape = *tls.shape;

    // 3) Call windowFindBegin.
    bool res = windowFindBegin(tls.swiftFinder, tls.swiftPattern, tls.options.errorRate);
    if (!res)
        std::cerr << "ERROR: windowFindBegin() failed in thread " << tls.threadId << std::endl;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
void
writeBackToLocal(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, JobData<Filtration<TFragmentStore, TOptions> > & jobData)
{
  // Write back all matches from verification to the thread local store.
  for (unsigned i = 0; i < length(jobData.verificationResults->localStores); ++i) {
    TFragmentStore * bucket = jobData.verificationResults->localStores[i];
    for (unsigned j = 0; j < length(bucket->alignedReadStore); ++j) {
      bucket->alignedReadStore[j].id = length(tls.localStore.alignedReadStore);
      appendValue(tls.localStore.alignedReadStore, bucket->alignedReadStore[j], Generous());
      appendValue(tls.localStore.alignQualityStore, bucket->alignQualityStore[j], Generous());
    }
  }
  // Possibly compact matches.
  if (length(tls.localStore.alignedReadStore) > tls.options.compactThresh)
  {
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typename Size<TAlignedReadStore>::Type oldSize = length(tls.localStore.alignedReadStore);

    if (TYPECMP<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
      maskDuplicates(tls.localStore, TRazerSMode());  // overlapping parallelograms cause duplicates

    fprintf(stderr, "COMPACTING MATCHES\n");
    compactMatches(tls.localStore, tls.counts, tls.options, TRazerSMode(), tls.swiftPattern, COMPACT);
    
    if (length(tls.localStore.alignedReadStore) * 4 > oldSize)      // the threshold should not be raised
      tls.options.compactThresh += (tls.options.compactThresh >> 1);  // if too many matches were removed
  }
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
void
workFiltration(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, JobData<Filtration<TFragmentStore, TOptions> > & jobData)
{
    fprintf(stderr, "Filtration\n");
    typedef JobData<Filtration<TFragmentStore, TOptions> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder, TOptions> > TVerificationJobData;

	typedef typename TSwiftFinder::THitString THitString;

    // First, wait for all verification jobs to complete, write back hits to global store and compact this store.
    while (jobData.verificationResults->blocksTotal != jobData.verificationResults->blocksDone) {
        // Busy waiting.
        SEQAN_ASSERT_LEQ(jobData.verificationResults->blocksTotal, jobData.verificationResults->blocksDone);
    }
    // Write back hits into local store or global store.
    writeBackToLocal(tls, jobData);
    // Clear verification results struct for now.
    clear(*jobData.verificationResults);

    // Make sure that the SWIFT filter is correctly initialized.
    initializeFiltration(tls, jobData);

    // Filter reads in the next window.
    bool referenceLeft = windowFindNext(tls.swiftFinder, tls.swiftPattern, tls.options.windowSize);
    THitString & hits = getSwiftHits(tls.swiftFinder);
    tls.stats.countFiltration += length(hits);
    // {
    //     int swiftHits = length(hits);
    //     fprintf(stderr, "Thread %ud Found %d SWIFT hits!\n", tls.threadId, swiftHits);
    // }

    // Enqueue filtration job for remaining sequence, if any.  We simply copy over the data.
    // TODO(holtgrew): Job with this data MUST NOT BE STOLEN since we update it below. Add to queue at all?
    if (referenceLeft)
        pushFront(jobQueue(tls), TJob(new TFiltrationJobData(jobData)));
    else
        windowFindEnd(tls.swiftFinder, tls.swiftPattern);

    // Enqueue verification jobs, if any.
    if (length(hits) == 0u) return;
    String<unsigned> splitters;
    computeSplittersBySlotSize(splitters, length(hits), tls.options.verificationPackageSize);

    String<TJob> jobs;
    reserve(jobs, length(splitters) - 1);
    std::tr1::shared_ptr<THitString> swiftHitsPtr(new THitString());
    std::swap(*swiftHitsPtr, hits);
    for (unsigned i = 1; i < length(splitters); ++i) {
        // Each verification job gets a local FragmentStore from the free list in the global state.
        TFragmentStore * localStorePtr;
        allocateLocalStore(localStorePtr, *jobData.verificationResults);
        appendValue(jobs, TJob(new TVerificationJobData(jobData.fragmentStore, jobData.contigId, jobData.blockId, jobData.verificationResults, localStorePtr, swiftHitsPtr, splitters[i - 1], splitters[i])));
    }
    pushFront(jobQueue(tls), jobs);
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TThreadId, typename TCounts, typename TRazerSMode>
void
workVerification(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, JobData<Verification<TFragmentStore, TSwiftFinder, TOptions> > & jobData, TThreadId thisThreadId, TThreadId jobThreadId)
{
    fprintf(stderr, "Verification\n");
    typedef JobData<Verification<TFragmentStore, TSwiftFinder, TOptions> > TVerificationJobData;
    typedef typename TVerificationJobData::THitString THitString;
    typedef typename Iterator<THitString>::Type THitStringIterator;

    unsigned offset = tls.globalState->splitters[jobData.blockId];
    for (THitStringIterator it = iter(*jobData.hitString, jobData.matchBeginIndex), itEnd = iter(*jobData.hitString, jobData.matchEndIndex); it != itEnd; ++it) {
        unsigned absReadId = offset + value(it).ndlSeqNo;
        tls.verifier.m.readId = absReadId;
  
        tls.verifier.store = jobData.localStore;

        // TODO(holtgrew): Unsure about global/local indices and global/local store.
		matchVerify(tls.verifier, swiftInfix(value(it), jobData.fragmentStore->contigStore[jobData.contigId].seq), absReadId, tls.globalStore->readSeqStore, TRazerSMode());
    }
    // Mark one more job as done in verification results.
    #pragma omp atomic
    jobData.verificationResults->blocksDone += 1;

    // TODO(holtgrew): Unused?
    (void) thisThreadId;
    (void) jobThreadId;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TThreadId, typename TCounts, typename TRazerSMode>
void
work(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, TJob & job, TThreadId thisThreadId, TThreadId jobThreadId)
{
    typedef JobData<Filtration<TFragmentStore, TOptions> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder, TOptions> > TVerificationJobData;

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
	typename TCounts,
	typename TSpec, 
	typename TShape,
    typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
    typename TMatchNPolicy >
void _mapSingleReadsParallelToContig(
	FragmentStore<TFSSpec, TFSConfig>						& store,
    GlobalState<FragmentStore<TFSSpec, TFSConfig>, RazerSOptions<TSpec> >         & globalState,
    TContigId const                                         & contigId,
	TCounts													& cnts,
    char                                                      orientation,
	RazerSOptions<TSpec>									& options,
	TShape const											& shape,
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const      & /*mode*/)
{
	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type						TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, Index_QGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;

    // TODO(holtgrew): What about cnts, mode?

	typedef typename IF<TYPECMP<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;
    typedef RazerSOptions<TSpec> TOptions;

    typedef Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > TJob;
    typedef ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > TThreadLocalStorage;
    typedef String<TThreadLocalStorage> TThreadLocalStorages;
    typedef JobData<Filtration<TFragmentStore, TOptions> > TFiltrationJobData;

    // Lock contig and possibly reverse-complement it.
	lockContig(store, contigId);
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
        threadLocalStorages[i].options = options;  // TODO(holtgrew): Copy for stats, really good?
        threadLocalStorages[i].globalState = &globalState;

        // Construct match verifier object.
        typedef typename TThreadLocalStorage::TMatchVerifier TVerifier;
        TVerifier oneVerifier(threadLocalStorages[i].localStore, *globalState.globalOptions, threadLocalStorages[i].swiftPatternHandler, cnts);
		oneVerifier.onReverseComplement = (orientation == 'R');
		oneVerifier.genomeLength = length(contigSeq);
		oneVerifier.m.contigId = contigId;
		// Assign it to the string.
		threadLocalStorages[i].verifier = oneVerifier;

        for (unsigned j = 0; j < options.splitFactor; ++j, ++k) {
            TJob job(new TFiltrationJobData(k, contigId, k, &store, &globalState));
            pushBack(jobQueue(threadLocalStorages[i]), job);
        }
    }

    // Fire up filtration jobs.
    int oldThreshold = globalState.globalOptions->compactThresh;
    globalState.globalOptions->compactThresh = SupremumValue<unsigned>::VALUE;
    work(threadLocalStorages, onlyVerificationStealable<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode>, StealOne());
    globalState.globalOptions->compactThresh = oldThreshold;

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
	typedef Index<TReadSet, Index_QGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

	typedef RazerSOptions<TSpec> TOptions;

    // Store global state we change below so we are side-effect free.
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount);

    GlobalState<TFragmentStore, TOptions> globalState;  // Store free list of local stores, splitters, for example.
    globalState.globalOptions = &options;

    // Compute initial load balancing / number of blocks.
    SEQAN_ASSERT_EQ_MSG(options.splitFactor, 1u, "For now...");
    unsigned blockCount = options.threadCount * options.splitFactor;
    computeSplittersBySlotCount(globalState.splitters, length(store.readNameStore), blockCount);
	if (options._debugLevel >= 1) {
		::std::cerr << ::std::endl << "Number of threads:               \t" << options.threadCount << std::endl;
		::std::cerr << ::std::endl << "Number of blocks:                \t" << blockCount << std::endl;
	}

    for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId) {
		lockContig(store, contigId);
		if (options.forward)
			_mapSingleReadsParallelToContig(store, globalState, contigId, cnts, 'F', options, shape, mode);
		if (options.reverse)
			_mapSingleReadsParallelToContig(store, globalState, contigId, cnts, 'R', options, shape, mode);
		unlockAndFreeContig(store, contigId);
    }

    // TODO(holtgrew): Sum up cnts?!

    // Restore global state.
    omp_set_num_threads(oldMaxThreads);

	return 0;
}

}  // namespace seqan

#endif  // #ifndef RAZERS_RAZERS_PARALLEL_H_
