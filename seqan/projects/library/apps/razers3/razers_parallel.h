#ifndef RAZERS_RAZERS_PARALLEL_H_
#define RAZERS_RAZERS_PARALLEL_H_

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

// For each block, we have one BlockLocalStorage object.
//
// Also works as the "parallel SWIFT pattern handler", i.e. deactivates
// needles in the SWIFT pattern when threshold is set to infinity in
// matchVerify.
template <typename TFragmentStore, typename TSwiftPattern, typename TOptions>
class BlockLocalStorage
{
public:
    TFragmentStore store;
    TSwiftPattern swiftPattern;
    // Need options locally since compaction threshold is modified.
    TOptions options;
};

// Uses the read ID to find the correct SWIFT pattern in the string handled by
// the ParallelSwiftPatternHandler, and the correct local ID within this SWIFT
// pattern to update the max errors.
template <typename TFragmentStore, typename TSwiftPattern, typename TOptions, typename TReadNo, typename TMaxErrors>
inline 
void 
setMaxErrors(BlockLocalStorage<TFragmentStore, TSwiftPattern, TOptions> & bls,
             TReadNo readNo,
             TMaxErrors maxErrors)
{
	int blockSize = length(host(host(bls.swiftPattern)));
	int localReadNo = readNo % blockSize;

	int minT = _qgramLemma(bls.swiftPattern, localReadNo, maxErrors);
	if (minT > 1){
		if (maxErrors < 0)
            minT = MaxValue<int>::VALUE;
		setMinThreshold(bls.swiftPattern, localReadNo, static_cast<unsigned>(minT));
	}
}

template <typename TFragmentStore,
          typename TSwiftPattern,
          typename TShape,
          typename TOptions>
void initializeBlockLocalStorages(
        String<BlockLocalStorage<TFragmentStore, TSwiftPattern, TOptions> > & blockLocalStorages,
        TFragmentStore /*const*/ & store,
        TOptions const & options,
        TShape const & shape,
        String<unsigned> /*const*/ & splitters)
{
    SEQAN_ASSERT_GEQ(length(splitters), 2u);
    int blockCount = length(splitters) - 1;
    resize(blockLocalStorages, blockCount, Exact());

    // Initialize block local SWIFT pattern.
    typedef typename Host<TSwiftPattern>::Type TIndex;
    typedef typename Position<typename TFragmentStore::TContigStore>::Type TPosition;

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < blockCount; ++i) {
        blockLocalStorages[i].options = options;

        // Clear pattern and set parameters.
        TSwiftPattern & swiftPattern = blockLocalStorages[i].swiftPattern;
        clear(swiftPattern);
        swiftPattern.params.minThreshold = options.threshold;
        swiftPattern.params.tabooLength = options.tabooLength;
        
        // Initialize the index.
        TIndex & index = host(blockLocalStorages[i].swiftPattern);
        clear(index);
        clear(indexText(index));
        TPosition jBegin = splitters[i];
        TPosition jEnd = splitters[i + 1];
        for (TPosition j = jBegin; j < jEnd; ++j)
            appendValue(value(index.text), store.readSeqStore[j]);
        index.shape = shape;

#ifdef RAZERS_OPENADDRESSING
        index.alpha = options.loadFactor;
#endif
        cargo(index).abundanceCut = options.abundanceCut;
        cargo(index)._debugLevel = options._debugLevel;

        indexRequire(index, QGramSADir());
    }
}

// Global state goes here.
//
// At the moment, this is the free list of local stores, only.
template <typename TFragmentStore, typename TSwiftPattern, typename TOptions>
class GlobalState
{
public:
    typedef BlockLocalStorage<TFragmentStore, TSwiftPattern, TOptions> TBlockLocalStorage;

    String<TBlockLocalStorage> blockLocalStorages;

    String<unsigned> splitters;
    String<TFragmentStore *> localStores;
    TOptions * globalOptions;
    unsigned nextJobId;
    omp_lock_t lock;

    GlobalState()
      : globalOptions(0), nextJobId(0)
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

template <typename TFragmentStore, typename TSwiftPattern, typename TOptions>
void allocateStore(TFragmentStore * & result, GlobalState<TFragmentStore, TSwiftPattern, TOptions> & globalState)
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

template <typename TFragmentStore, typename TSwiftPattern, typename TOptions>
void deallocateStore(GlobalState<TFragmentStore, TSwiftPattern, TOptions> & globalState, TFragmentStore * const & store)
{
    omp_set_lock(&globalState.lock);
    appendValue(globalState.localStores, store);
    omp_unset_lock(&globalState.lock);
}

template <typename TFragmentStore, typename TSwiftPattern, typename TOptions>
void deallocateStores(GlobalState<TFragmentStore, TSwiftPattern, TOptions> & globalState, String<TFragmentStore *> & stores)
{
    omp_set_lock(&globalState.lock);
    for (unsigned i = 0; i < length(stores); ++i)
        appendValue(globalState.localStores, stores[i]);
    omp_unset_lock(&globalState.lock);
}

struct WorkRecord
{
  double begin;
  double end;
  int type;

  WorkRecord() : begin(0), end(0), type(0) {}
};

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
    GlobalState<TFragmentStore, TSwiftPattern, TOptions> * globalState;
    TShape const * shape;
    TOptions options;

    TSwiftFinder swiftFinder;
    TSwiftPattern * swiftPattern;

    RazerSStats stats;
    TCounts counts;  // TODO(holtgrew): Artifact? ?!?!?

    typedef BlockLocalStorage<TFragmentStore, TSwiftPattern, TOptions> TBlockLocalStorage;
    TBlockLocalStorage * currentBlockLocalStorage;
    typedef MatchVerifier<TFragmentStore, TOptions, TRazerSMode, TBlockLocalStorage, TCounts> TMatchVerifier;
    TMatchVerifier verifier;

    String<WorkRecord> workRecords;

    ThreadLocalStorage() : previousFiltrationJobId(-1), globalState(0), currentBlockLocalStorage(0) {}
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
template <typename TFragmentStore, typename TSwiftPattern, typename TOptions>
class VerificationResults
{
public:
    omp_lock_t lock_;
    GlobalState<TFragmentStore, TSwiftPattern, TOptions> * globalState;
    volatile unsigned blocksTotal;
    volatile unsigned blocksDone;
    String<TFragmentStore *> localStores;

    VerificationResults(GlobalState<TFragmentStore, TSwiftPattern, TOptions> * globalState_)
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

template <typename TFragmentStore, typename TSwiftPattern, typename TOptions>
inline
void
allocateLocalStore(TFragmentStore * & localStorePtr, VerificationResults<TFragmentStore, TSwiftPattern, TOptions> & verificationResults)
{
    omp_set_lock(&verificationResults.lock_);
    verificationResults.blocksTotal += 1;
    allocateStore(localStorePtr, *verificationResults.globalState);
    appendValue(verificationResults.localStores, localStorePtr);
    // appendValue(verificationResults.localStores, localStorePtr);
    omp_unset_lock(&verificationResults.lock_);
}

template <typename TFragmentStore, typename TSwiftPattern, typename TOptions>
inline
void clear(VerificationResults<TFragmentStore, TSwiftPattern, TOptions> & verificationResults)
{
    omp_set_lock(&verificationResults.lock_);
    SEQAN_ASSERT_EQ_MSG(verificationResults.blocksTotal, verificationResults.blocksDone, "Can only clear completed verification results.");

    verificationResults.blocksTotal = 0;
    verificationResults.blocksDone = 0;
    deallocateStores(*verificationResults.globalState, verificationResults.localStores);
    clear(verificationResults.localStores);
    omp_unset_lock(&verificationResults.lock_);
}

template <typename TSpec>
class JobData;

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TOptions>
struct Filtration;

// Job data for filtration.
//
// SWIFT results are stored in the thread local storage's SWIFT
// filter's hit string.  The verification results are stored in the
// verification job, too, used for a barrier to wait for all
// verification jobs to complete before putting results into the
// central fragment store and possibly doing compaction there.
template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TOptions>
class JobData<Filtration<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> >
{
public:
	typedef typename TSwiftFinder::THitString THitString;
    typedef typename Size<typename TFragmentStore::TReadSeqStore>::Type TSize;

    // TODO(holtgrew): "active" flag to mark non-stealable
    int jobId;
    unsigned contigId;
    unsigned blockId;
    TFragmentStore * fragmentStore;  // global FragmentStore
    VerificationResults<TFragmentStore, TSwiftPattern, TOptions> * verificationResults;
    bool writeBackOnly;
    //THitString hitString;

    JobData(int jobId_, unsigned contigId_, unsigned blockId_, TFragmentStore *fragmentStore_, GlobalState<TFragmentStore, TSwiftPattern, TOptions> * globalState_)
            : jobId(jobId_), contigId(contigId_), blockId(blockId_), fragmentStore(fragmentStore_), verificationResults(new VerificationResults<TFragmentStore, TSwiftPattern, TOptions>(globalState_)), writeBackOnly(false)
    {}
};

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TOptions>
struct Verification;

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TOptions>
class JobData<Verification<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> >
{
public:
	typedef typename TSwiftFinder::THitString THitString;
    typedef typename Position<typename TSwiftFinder::THitString>::Type TPosition;

    TFragmentStore * fragmentStore;
    unsigned contigId;
    unsigned blockId;
    
    VerificationResults<TFragmentStore, TSwiftPattern, TOptions> * verificationResults;
    TFragmentStore * localStore;

    THitString & hitString;
    TPosition matchBeginIndex;
    TPosition matchEndIndex;

    JobData(TFragmentStore * fragmentStore_, unsigned contigId_, unsigned blockId_, VerificationResults<TFragmentStore, TSwiftPattern, TOptions> * & verificationResults_, TFragmentStore * & localStore_, THitString & hitString_, TPosition matchBeginIndex_, TPosition matchEndIndex_)
            : fragmentStore(fragmentStore_), contigId(contigId_), blockId(blockId_), verificationResults(verificationResults_), localStore(localStore_), hitString(hitString_), matchBeginIndex(matchBeginIndex_), matchEndIndex(matchEndIndex_)
    {
    }
};

template <typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
class Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> >
{
public:
    typedef JobData<Filtration<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > TVerificationJobData;

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
initializeFiltration(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, JobData<Filtration<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > const & jobData)
{
    // No need to initialize if we are doing the same filtration job as previously.
    // fprintf(stderr, "INITIALIZE: jobData.jobId == %d, tls.previousFiltrationJobId == %d\n", jobData.jobId, tls.previousFiltrationJobId);
    if (jobData.jobId == tls.previousFiltrationJobId)
        return;
    tls.previousFiltrationJobId = jobData.jobId;

    // Initialize Finder, pattern is initialized block-wise.
	tls.swiftFinder = TSwiftFinder(tls.globalStore->contigStore[jobData.contigId].seq, tls.options.repeatLength, 1);

    tls.swiftPattern = &tls.globalState->blockLocalStorages[jobData.blockId].swiftPattern;
    tls.swiftPattern->params.printDots = (tls.threadId == 0) && (tls.options._debugLevel > 0);

    // 3) Call windowFindBegin.
    bool res = windowFindBegin(tls.swiftFinder, *tls.swiftPattern, tls.options.errorRate);
    if (!res)
        std::cerr << "ERROR: windowFindBegin() failed in thread " << tls.threadId << std::endl;
}

template <typename TFragmentStore,
          typename TSwiftPattern,
          typename TOptions>
void
writeBackToGlobalStore(
        TFragmentStore & target,
        String<BlockLocalStorage<TFragmentStore, TSwiftPattern, TOptions> > /*const*/ & blockLocalStorages)
{
	typedef typename Size<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreSize;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
	typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
	
	// Update the IDs and calculate new size so the prefix increment can be
	// used in the loops.
	TAlignedReadStoreSize oldSize = length(target.alignedReadStore);
	TAlignedReadStoreSize sizeSum = --oldSize;  // TODO(holtgrew): Maybe no prefix-- here, below suffix++ and below that no ++.
	
	for (unsigned i = 0; i < length(blockLocalStorages); ++i)
		for (unsigned j = 0; j < length(blockLocalStorages[i].store.alignedReadStore); ++j)
			blockLocalStorages[i].store.alignedReadStore[j].id = ++sizeSum;
	++sizeSum;
	
	// Resize first so copying happens at most once and not every for each
	// block in the worst case
	resize(target.alignedReadStore, sizeSum, Exact());
	resize(target.alignQualityStore, sizeSum, Exact());			
	
	// Append single block stores.
	for (unsigned i = 0; i < length(blockLocalStorages); ++i) {
		for (unsigned j = 0; j < length(blockLocalStorages[i].store.alignedReadStore); ++j){
			target.alignedReadStore[++oldSize] = blockLocalStorages[i].store.alignedReadStore[j];
			target.alignQualityStore[oldSize] = blockLocalStorages[i].store.alignQualityStore[j];
		}
	}
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
void
writeBackToLocal(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, JobData<Filtration<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > & jobData)
{
    TFragmentStore & localStore = tls.globalState->blockLocalStorages[jobData.blockId].store;

    // TODO(holtgrew): Reserve enough space in localStore.

    // Write back all matches from verification to the block local store.
    for (unsigned i = 0; i < length(jobData.verificationResults->localStores); ++i) {
        TFragmentStore * bucket = jobData.verificationResults->localStores[i];
        for (unsigned j = 0; j < length(bucket->alignedReadStore); ++j) {
            bucket->alignedReadStore[j].id = length(localStore.alignedReadStore);
            appendValue(localStore.alignedReadStore, bucket->alignedReadStore[j], Generous());
            appendValue(localStore.alignQualityStore, bucket->alignQualityStore[j], Generous());
        }
    }

    SEQAN_ASSERT_EQ(length(localStore.alignedReadStore), length(localStore.alignQualityStore));
    if (!empty(localStore.alignedReadStore))
      SEQAN_ASSERT_EQ(length(localStore.alignedReadStore), back(localStore.alignedReadStore).id + 1);

    // Possibly compact matches.
    if (length(localStore.alignedReadStore) > tls.options.compactThresh)
    {
        typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
        typename Size<TAlignedReadStore>::Type oldSize = length(localStore.alignedReadStore);

        //fprintf(stderr, "[compact]");
        if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
          maskDuplicates(localStore, TRazerSMode());  // overlapping parallelograms cause duplicates

        compactMatches(localStore, tls.counts, tls.options, TRazerSMode(), tls.globalState->blockLocalStorages[jobData.blockId], COMPACT);
    
        if (length(localStore.alignedReadStore) * 4 > oldSize)      // the threshold should not be raised
            tls.options.compactThresh += (tls.options.compactThresh >> 1);  // if too many matches were removed
    }

    SEQAN_ASSERT_EQ(length(localStore.alignedReadStore), length(localStore.alignQualityStore));
    if (!empty(localStore.alignedReadStore))
      SEQAN_ASSERT_EQ(length(localStore.alignedReadStore), back(localStore.alignedReadStore).id + 1);
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TCounts, typename TRazerSMode>
void
workFiltration(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, JobData<Filtration<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > & jobData)
{
    //fprintf(stderr, "[filtration]");
    typedef JobData<Filtration<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > TVerificationJobData;

	typedef typename TSwiftFinder::THitString THitString;

    // First, wait for all verification jobs to complete, write back hits to thread local store and compact this store.
    typedef VerificationResults<TFragmentStore, TSwiftPattern, TOptions> TVerificationResults;
    // TODO(holtgrew): Verification results should be block local.
    TVerificationResults * verificationResults = jobData.verificationResults;
    while (jobData.verificationResults->blocksTotal != jobData.verificationResults->blocksDone) {
        // Busy waiting.
        SEQAN_ASSERT_LEQ(jobData.verificationResults->blocksTotal, jobData.verificationResults->blocksDone);
    }
    SEQAN_ASSERT_EQ(jobData.verificationResults->blocksTotal, length(jobData.verificationResults->localStores));
    // Set current block local storage for disabling in compaction.
    tls.currentBlockLocalStorage = &tls.globalState->blockLocalStorages[jobData.blockId];
    // Write back hits into local store or global store.
    writeBackToLocal(tls, jobData);
    // Clear verification results struct for now.
    clear(*jobData.verificationResults);
    // Set current swift handler to none, should never be called in thread-local verification.
    tls.currentBlockLocalStorage = 0;

    if (jobData.writeBackOnly)
      return;

    // Make sure that the SWIFT finder is correctly initialized.
    initializeFiltration(tls, jobData);

    // Filter reads in the next window.
    bool referenceLeft = windowFindNext(tls.swiftFinder, *tls.swiftPattern, tls.options.windowSize);
    THitString & hits = getSwiftHits(tls.swiftFinder);
    tls.globalState->blockLocalStorages[jobData.blockId].options.countFiltration += length(hits);
    // {
    //     int swiftHits = length(hits);
    //     fprintf(stderr, "Thread %ud Found %d SWIFT hits!\n", tls.threadId, swiftHits);
    // }

    // Enqueue filtration job for remaining sequence, if any.  We simply copy over the data.
    // TODO(holtgrew): Job with this data MUST NOT BE STOLEN since we update it below. Add to queue at all?
    if (referenceLeft) {
        TFiltrationJobData * data = new TFiltrationJobData(jobData);
        verificationResults = data->verificationResults;
        pushFront(jobQueue(tls), TJob(data));
    } else {
        TFiltrationJobData * data = new TFiltrationJobData(jobData);
        data->writeBackOnly = true;
        verificationResults = data->verificationResults;
        pushFront(jobQueue(tls), TJob(data));
        tls.swiftPattern->params.printDots = false;
        windowFindEnd(tls.swiftFinder, *tls.swiftPattern);
    }

    // Enqueue verification jobs, if any.
    if (length(hits) == 0u) return;
    String<unsigned> splitters;
    computeSplittersBySlotSize(splitters, length(hits), tls.options.verificationPackageSize);

    String<TJob> jobs;
    reserve(jobs, length(splitters) - 1);
    //clear(jobData.hitString);
    //swap(jobData.hitString, hits);
    for (unsigned i = 1; i < length(splitters); ++i) {
        unsigned j = length(splitters) - i;
        // Each verification job gets a local FragmentStore from the free list in the global state.
        TFragmentStore * localStorePtr;
        allocateLocalStore(localStorePtr, *verificationResults);
        appendValue(jobs, TJob(new TVerificationJobData(jobData.fragmentStore, jobData.contigId, jobData.blockId, verificationResults, localStorePtr, hits, splitters[j - 1], splitters[j])));
    }
    pushFront(jobQueue(tls), jobs);
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TThreadId, typename TCounts, typename TRazerSMode>
void
workVerification(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls, JobData<Verification<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > & jobData, TThreadId thisThreadId, TThreadId jobThreadId)
{
    //fprintf(stderr, "[verification]");
    typedef JobData<Verification<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > TVerificationJobData;
    typedef typename TVerificationJobData::THitString THitString;
    typedef typename Iterator<THitString>::Type THitStringIterator;

    // Initialize verifier.
    tls.verifier.store = jobData.localStore;
    tls.verifier.options = &tls.globalState->blockLocalStorages[jobData.blockId].options;
    tls.verifier.swiftPattern = &tls.globalState->blockLocalStorages[jobData.blockId];
    tls.verifier.cnts = 0;

    unsigned offset = tls.globalState->splitters[jobData.blockId];
// #ifdef RAZERS_DEBUG
//     std::cerr << "Verifying block " << jobData.blockId << std::endl;
//     std::cerr << "  offset = " << offset << std::endl;
//     std::cerr << "  jobData.matchBeginIndex == " << jobData.matchBeginIndex << std::endl;
// #endif  // #ifdef RAZERS_DEBUG
    for (THitStringIterator it = iter(jobData.hitString, jobData.matchBeginIndex), itEnd = iter(jobData.hitString, jobData.matchEndIndex); it != itEnd; ++it) {
        // std::cerr << "value(it).ndlSeqNo == " << value(it).ndlSeqNo << std::endl;
        unsigned absReadId = offset + value(it).ndlSeqNo;
        tls.verifier.m.readId = absReadId;
  
        // TODO(holtgrew): Unsure about global/local indices and global/local store.
		matchVerify(tls.verifier, swiftInfix(value(it), jobData.fragmentStore->contigStore[jobData.contigId].seq), absReadId, tls.globalStore->readSeqStore, TRazerSMode());
    }
    // Mark one more job as done in verification results.
    #pragma omp atomic
    jobData.verificationResults->blocksDone += 1;

    // Reset verifier.
    tls.verifier.store = 0;
    tls.verifier.options = 0;
    tls.verifier.swiftPattern = 0;
    tls.verifier.cnts = 0;

    // TODO(holtgrew): Unused?
    (void) thisThreadId;
    (void) jobThreadId;
}

template <typename TJob, typename TFragmentStore, typename TSwiftFinder, typename TSwiftPattern, typename TShape, typename TOptions, typename TThreadId, typename TCounts, typename TRazerSMode>
void
work(ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > & tls,
     TJob & job,
     TThreadId thisThreadId,
     TThreadId jobThreadId)
{
    typedef JobData<Filtration<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > TFiltrationJobData;
    typedef JobData<Verification<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > TVerificationJobData;

    appendValue(tls.workRecords, WorkRecord(

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
	typename TThreadLocalStorages,
    typename TContigId,
	typename TCounts,
    typename TSwiftPattern,
	typename TSpec, 
	typename TShape,
    typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
    typename TMatchNPolicy >
void _mapSingleReadsParallelToContig(
	FragmentStore<TFSSpec, TFSConfig>						& store,
	TThreadLocalStorages & threadLocalStorages,
    GlobalState<FragmentStore<TFSSpec, TFSConfig>, TSwiftPattern, RazerSOptions<TSpec> >         & globalState,
    TContigId const                                         & contigId,
	TCounts													& /*cnts*/,
    char                                                      orientation,
	RazerSOptions<TSpec>									& options,
	TShape const											& /*shape*/,
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const      & /*mode*/)
{
	typedef FragmentStore<TFSSpec, TFSConfig>						TFragmentStore;
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef typename TFragmentStore::TReadSeqStore					TReadSeqStore;
	typedef typename Value<TReadSeqStore>::Type						TRead;
	typedef StringSet<TRead>										TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape, OpenAddressing> >	TIndex;			// q-gram index
	typedef typename Size<TReadSeqStore>::Type						TSize;

    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;

    // TODO(holtgrew): What about cnts, mode?

	typedef typename If<IsSameType<TGapMode,RazerSGapped>::VALUE, SwiftSemiGlobal, SwiftSemiGlobalHamming>::Type TSwiftSpec;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
	// typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;
    typedef RazerSOptions<TSpec> TOptions;

    typedef Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > TJob;
    typedef JobData<Filtration<TFragmentStore, TSwiftFinder, TSwiftPattern, TOptions> > TFiltrationJobData;
	
	// Debug output...
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F') ::std::cerr << "[fwd]";
		else                    ::std::cerr << "[rev]";
	}

    // Lock contig and possibly reverse-complement it.
	lockContig(store, contigId);
	TContigSeq & contigSeq = store.contigStore[contigId].seq;
	if (orientation == 'R')
		reverseComplement(contigSeq);

    // Fill thread local storages with filtration jobs.
    for (unsigned i = 0, k = 0; i < options.threadCount; ++i) {
		// Initialize verifier object.
		threadLocalStorages[i].verifier.onReverseComplement = (orientation == 'R');
		threadLocalStorages[i].verifier.genomeLength = length(contigSeq);
		threadLocalStorages[i].verifier.oneMatchPerBucket = false;
		threadLocalStorages[i].verifier.m.contigId = contigId;

        for (unsigned j = 0; j < options.splitFactor; ++j, ++k) {
            TJob job(new TFiltrationJobData(globalState.nextJobId, contigId, k, &store, &globalState));
            pushBack(jobQueue(threadLocalStorages[i]), job);
            globalState.nextJobId += 1;
        }
    }

    // Fire up filtration jobs.
    int oldThreshold = globalState.globalOptions->compactThresh;
    globalState.globalOptions->compactThresh = MaxValue<unsigned>::VALUE;
    work(threadLocalStorages, onlyVerificationStealable<TFragmentStore, TSwiftFinder,
         TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode>, StealOne());
    globalState.globalOptions->compactThresh = oldThreshold;

    // Restore orientation of contig
    // TODO(holtgrew): Always reverse-complemented again since outer caller locks!
	if (!unlockAndFreeContig(store, contigId))						// if the contig is still used
        if (orientation == 'R')	reverseComplement(contigSeq);	// we have to restore original orientation
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

    // Store global state we change below so we are side-effect free.
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount);

    // Global state stores free list of local stores, splitters, for example.
    GlobalState<TFragmentStore, TSwiftPattern, TOptions> globalState;
    globalState.globalOptions = &options;

    double beginInit = sysTime();
    // Compute initial load balancing / number of blocks and initialize the
    // block-specific data structures.
    SEQAN_ASSERT_GEQ(options.splitFactor, 1u);
    unsigned blockCount = options.threadCount * options.splitFactor;
    computeSplittersBySlotCount(globalState.splitters, length(store.readNameStore), blockCount);
	if (options._debugLevel >= 1) {
		::std::cerr << ::std::endl << "Number of threads:               \t" << options.threadCount << std::endl;
		::std::cerr << ::std::endl << "Number of blocks:                \t" << blockCount << std::endl;
	}
    initializeBlockLocalStorages(globalState.blockLocalStorages, store, options, shape, globalState.splitters);

    // Create thread local storages.
    typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
    typedef Job<MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > TJob;
    typedef ThreadLocalStorage<TJob, MapSingleReads<TFragmentStore, TSwiftFinder, TSwiftPattern, TShape, TOptions, TCounts, TRazerSMode> > TThreadLocalStorage;
    typedef String<TThreadLocalStorage> TThreadLocalStorages;

    TThreadLocalStorages threadLocalStorages;
    resize(threadLocalStorages, options.threadCount);
    for (unsigned i = 0; i < options.threadCount; ++i) {
        threadLocalStorages[i].threadId = i;
        threadLocalStorages[i].globalStore = &store;
        threadLocalStorages[i].shape = &shape;
        threadLocalStorages[i].options = options;  // TODO(holtgrew): Copy for stats and threshold, really good?
        threadLocalStorages[i].globalState = &globalState;
    }
    double endInit = sysTime();
    std::cerr << "TIME initialization: " << (endInit - beginInit) << " s" << std::endl;

    // For each contig: Map reads in parallel.
    for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId) {
		lockContig(store, contigId);
		double beginForward = sysTime();
		if (options.forward)
			_mapSingleReadsParallelToContig(store, threadLocalStorages, globalState, contigId, cnts, 'F', options, shape, mode);
		double beginReverse = sysTime();
        std::cerr << "  TIME forward mapping (" << contigId << ") " << (beginReverse - beginForward) << " s" << std::endl;
		if (options.reverse)
			_mapSingleReadsParallelToContig(store, threadLocalStorages, globalState, contigId, cnts, 'R', options, shape, mode);
		double endReverse = sysTime();
        std::cerr << "  TIME reverse mapping (" << contigId << ") " << (endReverse - beginReverse) << " s" << std::endl;
		unlockAndFreeContig(store, contigId);
    }
    double endMapping = sysTime();
    std::cerr << "TIME mapping: " << (endMapping - endInit) << " s" << std::endl;

    // Write back local stores to global stores.
    writeBackToGlobalStore(store, globalState.blockLocalStorages);
    double endWriteback = sysTime();
    std::cerr << "TIME back to global: " << (endWriteback - endMapping) << " s" << std::endl;

    // TODO(holtgrew): Sum up cnts?!

    // Restore global state.
    omp_set_num_threads(oldMaxThreads);

    for (unsigned i = 0; i < length(globalState.blockLocalStorages); ++i) {
        options.countFiltration += globalState.blockLocalStorages[i].options.countFiltration;
        options.countVerification += globalState.blockLocalStorages[i].options.countVerification;
    }

	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:      " << options.countFiltration << ::std::endl;
		::std::cerr << "Successful verfications: " << options.countVerification << ::std::endl;
	}
	return 0;
}

}  // namespace seqan

#endif  // #ifndef RAZERS_RAZERS_PARALLEL_H_
