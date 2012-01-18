
#ifndef RAZERS_MATCH_FILTER_H_
#define RAZERS_MATCH_FILTER_H_

#include <tr1/unordered_map>

#include <seqan/graph_types/graph_idmanager.h>

#include "razers.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
class MatchFilter
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
    std::tr1::unordered_map<unsigned, unsigned> readIdToHistogramId;
    // Id manager for histogram.
    IdManager<unsigned> idManager;
    // Histograms.
    String<String<unsigned> > histograms;
    // Ids of the reads that are purged.
    String<unsigned> purgedReadIds;
    // The callback context object.
    Holder<TCallback> callback;
    // Read ID offset.
    unsigned readOffset;
    // The read sequences.
    TReadSeqSet const & readSeqs;
    // The options object.
    RazerSOptions<TOptionsSpec> const & options;

    MatchFilter(unsigned readCount_, unsigned matchThreshold_, double frac_, TCallback & callback_, unsigned readOffset_, TReadSeqSet const & readSeqs_, RazerSOptions<TOptionsSpec> const & options_)
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
_createHistogram(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned readId)
{
    unsigned result = obtainId(filter.idManager);
    if (result >= length(filter.histograms))
        resize(filter.histograms, length(filter.histograms) + 1);
    resize(filter.histograms[result], ceil(filter.options.errorRate * length(filter.readSeqs[readId])), 0, Exact());
    SEQAN_ASSERT_LT(result, length(filter.histograms));
    return result;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline void
_freeHistogram(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned histogramId)
{
    clear(filter.histograms[histogramId]);
    releaseId(filter.idManager, histogramId);
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline void
_incrementCount(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned histogramId, int score)
{
    SEQAN_ASSERT_LEQ(score, 0);
    filter.histograms[histogramId][-score] += 1;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline bool
_canBePurged(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> const & filter, unsigned histogramId)
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
_canBeDisabled(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> const & filter, unsigned histogramId)
{
    if (filter.options.maxHits == 0)
        return false;
    // std::cerr << "histo " << histogramId << "\t" << filter.histograms[histogramId][0] << "\t" << filter.options.maxHits << std::endl;
    return filter.histograms[histogramId][0] >= filter.options.maxHits;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline int
_newLimit(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> const & filter, unsigned histogramId)
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
registerRead(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned readId, int score)
{
    // std::cerr << "registering read " << readId << std::endl;
    if (filter.hitCount[readId - filter.readOffset] == MaxValue<unsigned>::VALUE)
        return;
    filter.hitCount[readId - filter.readOffset] += 1;

    // TODO(holtgrew): Maybe global read to histogram map; faster?

    // Get histogram id, insert new histogram if necessary, exit if no histogram yet.
    unsigned histogramId = 0;
    if (filter.hitCount[readId - filter.readOffset] == filter.matchThreshold) {
        // std::cerr << "new histogram for read " << readId << std::endl;
        histogramId = _createHistogram(filter, readId);
        filter.readIdToHistogramId[readId] = histogramId;
    } else if (filter.hitCount[readId - filter.readOffset] > filter.matchThreshold) {
        // std::cerr << "updating histogram for read " << readId << std::endl;
        typedef typename std::tr1::unordered_map<unsigned, unsigned>::iterator TIterator;
        TIterator it = filter.readIdToHistogramId.find(readId);
        SEQAN_ASSERT(it != filter.readIdToHistogramId.end());
        histogramId = it->second;
    } else {
        return;
    }

    // Insert value into histogram.
    _incrementCount(filter, histogramId, score);
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
bool
processRead(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned readId)
{
    typedef typename std::tr1::unordered_map<unsigned, unsigned>::iterator TIterator;

    if (filter.hitCount[readId - filter.readOffset] < filter.matchThreshold)
        return false;

     // std::cerr << "processing read " << readId << std::endl;
    // Get histogram id, insert new histogram if necessary, exit if no histogram yet.
    TIterator it = filter.readIdToHistogramId.find(readId);
    if (it == filter.readIdToHistogramId.end())
        return false;  // Must have been disabled before.
    unsigned histogramId = it->second;

    // Perform actions.
    int newLimit;
    if (_canBePurged(filter, histogramId)) {
        // std::cerr << "PURGED " << readId << std::endl;
        appendValue(filter.purgedReadIds, readId);
        disableRead(value(filter.callback), readId);
        _freeHistogram(filter, histogramId);
        filter.readIdToHistogramId.erase(readId);
        filter.hitCount[readId - filter.readOffset] = MaxValue<unsigned>::VALUE;
        return true;
    } else if (_canBeDisabled(filter, histogramId)) {
        // std::cerr << "DISABLED " << readId << "\t" << filter.histograms[histogramId][0] << "\t" << filter.hitCount[readId - filter.readOffset] << std::endl;
        disableRead(value(filter.callback), readId);
        _freeHistogram(filter, histogramId);
        filter.readIdToHistogramId.erase(readId);
        filter.hitCount[readId - filter.readOffset] = MaxValue<unsigned>::VALUE;
        return true;
    } else if ((newLimit = _newLimit(filter, histogramId)) >= 0) {
        // std::cerr << "LIMITING " << readId << "\t" << filter.histograms[histogramId][0] << "\t" << filter.hitCount[readId - filter.readOffset] << "\t" << newLimit << std::endl;
        limitRead(value(filter.callback), readId, newLimit);
    }
    return false;
}

}  // namespace seqan

#endif  // #ifndef RAZERS_MATCH_FILTER_H_
