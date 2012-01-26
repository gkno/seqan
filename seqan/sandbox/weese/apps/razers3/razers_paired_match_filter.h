#ifndef RAZERS_PAIRED_MATCH_FILTER_H_
#define RAZERS_PAIRED_MATCH_FILTER_H_

#include <tr1/unordered_map>

#include <seqan/graph_types/graph_idmanager.h>

#include "razers.h"

namespace seqan {

template <typename TThreadLocalStorage>
class FilterPatternLSetMaxErrorsWrapper;
template <typename TThreadLocalStorage>
class FilterPatternRSetMaxErrorsWrapper;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
class PairedMatchFilter
{
public:
    // Number of reads.
    unsigned readCount;
    // If a read has more than matchThreshold matches then it gets a histogram.
    unsigned matchThreshold;
    // Fraction of reads expect to have histograms for.
    double frac;
    // Count matches for each read
    String<unsigned> hitCount;
    // Map from read number to histogram id.
    std::tr1::unordered_map<unsigned, unsigned> pairIdToHistogramId;
    // Id manager for histogram.
    IdManager<unsigned> idManager;
    // Histograms.
    String<String<unsigned> > histograms;
    // Ids of the reads that are purged.
    String<unsigned> purgedPairIds;
    // The callback context object.
    Holder<TCallback> callback;
    // Read ID offset.
    unsigned readOffset;
    // The read sequences.
    TReadSeqSet const & readSeqs;
    // The options object.
    RazerSOptions<TOptionsSpec> const & options;

    PairedMatchFilter(unsigned readCount_, unsigned matchThreshold_, double frac_, TCallback & callback_, unsigned readOffset_, TReadSeqSet const & readSeqs_, RazerSOptions<TOptionsSpec> const & options_)
            : readCount(readCount_), matchThreshold(matchThreshold_), frac(frac_), callback(callback_), readOffset(readOffset_), readSeqs(readSeqs_), options(options_)
    {
        resize(hitCount, readCount, 0, Exact());
        reserve(histograms, unsigned(frac * readCount));
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline unsigned
_createHistogram(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned pairId)
{
    unsigned result = obtainId(filter.idManager);
    if (result >= length(filter.histograms))
        resize(filter.histograms, length(filter.histograms) + 1);
    resize(filter.histograms[result], ceil(filter.options.errorRate * (length(filter.readSeqs[2 * pairId]) + length(filter.readSeqs[2 * pairId + 1]))), 0, Exact());
    SEQAN_ASSERT_LT(result, length(filter.histograms));
    return result;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline void
_freeHistogram(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned histogramId)
{
    clear(filter.histograms[histogramId]);
    releaseId(filter.idManager, histogramId);
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline void
_incrementCount(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned histogramId, int score)
{
    SEQAN_ASSERT_LEQ(score, 0);
    filter.histograms[histogramId][-score] += 1;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline bool
_canBePurged(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> const & filter, unsigned histogramId)
{
    // TODO(holtgrew): This could be speeded up if using prefix sum data structures for histograms.
    if (!filter.options.purgeAmbiguous || filter.options.maxHits == 0)
        return false;
    unsigned sumUpTo = length(filter.histograms[histogramId]);
    if (filter.options.errorDistanceRange != 0)
        sumUpTo = filter.options.errorDistanceRange;
    // std::cerr << "sumUpTo == " << sumUpTo << std::endl;
    unsigned sum = 0;
    for (unsigned i = 0; i < sumUpTo; ++i) {
        sum += filter.histograms[histogramId][i];
        if (sum > filter.options.maxHits) {
            // std::cerr << "sum == " << sum << std::endl;
            return true;
        }
    }
    return false;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline bool
_canBeDisabled(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> const & filter, unsigned histogramId)
{
    if (filter.options.maxHits == 0)
        return false;
    // std::cerr << "histo " << histogramId << "\t" << filter.histograms[histogramId][0] << "\t" << filter.options.maxHits << std::endl;
    return filter.histograms[histogramId][0] >= filter.options.maxHits;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline int
_newLimit(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> const & filter, unsigned histogramId)
{
    // TODO(holtgrew): This could be speeded up if using prefix sum data structures for histograms.
    if (filter.options.maxHits == 0)
        return -1;
    unsigned sumUpTo = length(filter.histograms[histogramId]);
    unsigned sum = 0;
    for (unsigned i = 0; i < sumUpTo; ++i) {
        sum += filter.histograms[histogramId][i];
        if (sum > filter.options.maxHits)
            return i;
    }
    return -1;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
void
registerRead(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned pairId, int score)
{
    // std::cerr << "registering read " << pairId << std::endl;
    if (filter.hitCount[pairId - filter.readOffset] == MaxValue<unsigned>::VALUE)
        return;
    filter.hitCount[pairId - filter.readOffset] += 1;

    // TODO(holtgrew): Maybe global read to histogram map; faster?

    // Get histogram id, insert new histogram if necessary, exit if no histogram yet.
    unsigned histogramId = 0;
    if (filter.hitCount[pairId - filter.readOffset] == filter.matchThreshold) {
        // std::cerr << "new histogram for read " << pairId << std::endl;
        histogramId = _createHistogram(filter, pairId);
        filter.pairIdToHistogramId[pairId] = histogramId;
    } else if (filter.hitCount[pairId - filter.readOffset] > filter.matchThreshold) {
        // std::cerr << "updating histogram for read " << pairId << std::endl;
        typedef typename std::tr1::unordered_map<unsigned, unsigned>::iterator TIterator;
        TIterator it = filter.pairIdToHistogramId.find(pairId);
        SEQAN_ASSERT(it != filter.pairIdToHistogramId.end());
        histogramId = it->second;
    } else {
        return;
    }

    // Insert value into histogram.
    _incrementCount(filter, histogramId, score);
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
bool
processRead(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned pairId)
{
    typedef typename std::tr1::unordered_map<unsigned, unsigned>::iterator TIterator;

    if (filter.hitCount[pairId - filter.readOffset] < filter.matchThreshold)
        return false;

     // std::cerr << "processing read " << pairId << std::endl;
    // Get histogram id, insert new histogram if necessary, exit if no histogram yet.
    TIterator it = filter.pairIdToHistogramId.find(pairId);
    if (it == filter.pairIdToHistogramId.end())
        return false;  // Must have been disabled before.
    unsigned histogramId = it->second;

    // Perform actions.
    int newLimit;
    if (_canBePurged(filter, histogramId)) {
        // std::cerr << "PURGED " << pairId << std::endl;
        appendValue(filter.purgedPairIds, pairId);
        FilterPatternLSetMaxErrorsWrapper<TCallback> wrapperL(value(filter.callback));
        FilterPatternRSetMaxErrorsWrapper<TCallback> wrapperR(value(filter.callback));
        setMaxErrors(wrapperL, pairId, -1);
        setMaxErrors(wrapperR, pairId, -1);
        _freeHistogram(filter, histogramId);
        filter.pairIdToHistogramId.erase(pairId);
        filter.hitCount[pairId - filter.readOffset] = MaxValue<unsigned>::VALUE;
        return true;
    } else if (_canBeDisabled(filter, histogramId)) {
        // std::cerr << "DISABLED " << pairId << "\t" << filter.histograms[histogramId][0] << "\t" << filter.hitCount[pairId - filter.readOffset] << std::endl;
        FilterPatternLSetMaxErrorsWrapper<TCallback> wrapperL(value(filter.callback));
        FilterPatternRSetMaxErrorsWrapper<TCallback> wrapperR(value(filter.callback));
        setMaxErrors(wrapperL, pairId, -1);
        setMaxErrors(wrapperR, pairId, -1);
        _freeHistogram(filter, histogramId);
        filter.pairIdToHistogramId.erase(pairId);
        filter.hitCount[pairId - filter.readOffset] = MaxValue<unsigned>::VALUE;
        return true;
    } else if ((newLimit = _newLimit(filter, histogramId)) >= 0) {
        // std::cerr << "LIMITING " << pairId << "\t" << filter.histograms[histogramId][0] << "\t" << filter.hitCount[pairId - filter.readOffset] << "\t" << newLimit << std::endl;
        FilterPatternLSetMaxErrorsWrapper<TCallback> wrapperL(value(filter.callback));
        FilterPatternRSetMaxErrorsWrapper<TCallback> wrapperR(value(filter.callback));
        setMaxErrors(wrapperL, pairId, newLimit);
        setMaxErrors(wrapperR, pairId, newLimit);
    }
    return false;
}

}  // namespace seqan

#endif  // #ifndef RAZERS_PAIRED_MATCH_FILTER_H_
