#ifndef BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_
#define BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/store.h>  // For Sort* tags.

#include "intervals.h"

using namespace seqan;

// Stores closed intervals [firstPos, lastPos].
struct IntervalOfReadOnContig {
    static size_t invalidId() { return ~0ul; }
    static size_t additionalId() { return ~1ul; }
    static size_t superflousId() { return ~2ul; }
    
    size_t id;
    size_t readId;
    unsigned distance;
    size_t contigId;
    bool isForward;
    size_t firstPos;
    size_t lastPos;

    IntervalOfReadOnContig()
            : id(invalidId()), readId(0), distance(0), contigId(0), isForward(0), firstPos(0),
              lastPos(0) {}

    IntervalOfReadOnContig(size_t _readId, unsigned _distance, size_t _contigId, bool _isForward, size_t _firstPos, size_t _lastPos)
            : id(invalidId()), readId(_readId), distance(_distance), contigId(_contigId),
              isForward(_isForward), firstPos(_firstPos), lastPos(_lastPos) {}
};


template <typename TStream>
inline
TStream & operator<<(TStream & stream, IntervalOfReadOnContig const & record) {
    stream << "(" << record.readId << ", " << record.contigId << ", "
           << record.distance << ", " << record.isForward << ", "
           << record.firstPos << ", " << record.lastPos << ")";
    return stream;
}


struct WitStore {
    typedef StringSet<CharString> TNameSet;
    typedef String<IntervalOfReadOnContig> TIntervalStore;

    Holder<TNameSet> readNames;
    Holder<TNameSet> contigNames;

    // TODO(holtgrew): Rename to witRecords.
    TIntervalStore intervals;
};


template <typename TStream>
inline
TStream & operator<<(TStream & stream, WitStore const & store) {
    typedef typename Iterator<typename WitStore::TIntervalStore, Standard>::Type TIterator;
    stream << ",-- WIT Store" << std::endl;
    for (TIterator it = begin(store.intervals, Standard()); it != end(store.intervals, Standard()); ++it) {
        stream << "| " << value(store.readNames)[value(it).readId] << " read id = " << value(it).readId << "\t" << value(it).distance << "\t" << value(store.contigNames)[value(it).contigId]
               << "\t" << (value(it).isForward ? "F" : "R") << "\t" << value(it).firstPos << "\t" << value(it).lastPos << std::endl;
    }
    stream << "`--" << std::endl;
    return stream;
}

void move(WitStore & target, WitStore & source) {
    move(target.readNames, source.readNames);
    move(target.contigNames, source.contigNames);
    move(target.intervals, source.intervals);
}


template <typename TSpec>
struct WitStoreLess;


template <>
struct WitStoreLess<SortReadId>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.readId < b.readId;
    }
};


template <>
struct WitStoreLess<SortId>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.id < b.id;
    }
};


template <>
struct WitStoreLess<SortContigId>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.contigId < b.contigId;
    }
};


struct _SortDistance {};
typedef Tag<_SortDistance> const SortDistance;

template <>
struct WitStoreLess<SortDistance>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.distance < b.distance;
    }
};


struct _SortFirstPos {};
typedef Tag<_SortFirstPos> const SortFirstPos;


template <>
struct WitStoreLess<SortFirstPos>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.firstPos < b.firstPos;
    }
};


struct _SortLastPos {};
typedef Tag<_SortLastPos> const SortLastPos;


template <>
struct WitStoreLess<SortLastPos>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.lastPos < b.lastPos;
    }
};


template <typename TSortTag>
void
sortWitRecords(WitStore const & store, TSortTag &) {
    std::stable_sort(begin(store.intervals, Standard()), end(store.intervals, Standard()), WitStoreLess<TSortTag>(store));
}


inline
void appendValue(WitStore & store, IntervalOfReadOnContig const & record) {
    IntervalOfReadOnContig tmp(record);
    tmp.id = length(store.intervals);
    appendValue(store.intervals, tmp);
}


void loadWitFile(WitStore & store,
                 StringSet<CharString> /*const*/ & readNames,
                 StringSet<CharString> /*const*/ & contigNames,
                 CharString const & fileName) {
    // Assign read and contig names into wit store members.
    set(store.readNames, readNames);
    set(store.contigNames, contigNames);

    // Build mapping from existing read and contig names to their indices
    // (ids) in readNames and contigNames.
    NameStoreCache<StringSet<CharString> > readNameCache(readNames);
    refresh(readNameCache);
    NameStoreCache<StringSet<CharString> > contigNameCache(contigNames);
    refresh(contigNameCache);

    // Then, load the file.
    std::fstream fstrm(toCString(fileName), std::ios_base::in | std::ios_base::binary);
    if (!fstrm.is_open()) {
        std::cerr << "Could not open WIT file " << fileName << std::endl;
        exit(1);
    }

    // Read header.
    char c = '\0';
    readWitHeader(fstrm, c);

    // Temporary data for one record.
    CharString readName;
    size_t distance;
    CharString contigName;
    bool isForward;
    size_t firstPos;
    size_t lastPos;

    // Read WIT file.
    while (!_streamEOF(fstrm)) {
        // Skip comments.
        if (c == '#') {
            _parse_skipLine(fstrm, c);
            continue;
        }

        clear(readName);
        clear(contigName);

        // Read line.
        _parse_readIdentifier(fstrm, readName, c);
        _parse_skipWhitespace(fstrm, c);
        distance = _parse_readNumber(fstrm, c);
        _parse_skipWhitespace(fstrm, c);
        _parse_readIdentifier(fstrm, contigName, c);
        _parse_skipWhitespace(fstrm, c);
        isForward = (_parse_readChar(fstrm, c) == 'F');
        _parse_skipWhitespace(fstrm, c);
        firstPos = _parse_readNumber(fstrm, c);
        _parse_skipWhitespace(fstrm, c);
        lastPos = _parse_readNumber(fstrm, c);
        _parse_skipLine(fstrm, c);

        // Insert record into read store.
        IntervalOfReadOnContig record;
        getIdByName(store.readNames, readName, record.readId, readNameCache);
        record.distance = distance;
        getIdByName(store.contigNames, contigName, record.contigId, contigNameCache);
        record.isForward = isForward;
        record.firstPos = firstPos;
        record.lastPos = lastPos;
        appendValue(store, record);
    }
}


template <typename TStream>
void writeWitFile(TStream & stream, WitStore const & witStore) {
    writeWitHeader(stream);
    writeWitComment(stream , WIT_COLUMN_NAMES);
    typedef typename WitStore::TIntervalStore TIntervalStore;
    typedef typename Iterator<TIntervalStore, Standard>::Type TIterator;
    for (TIterator it = begin(witStore.intervals, Standard()); it != end(witStore.intervals, Standard()); ++it) {
        stream << value(witStore.readNames)[it->readId] << '\t'
               << it->distance << '\t'
               << value(witStore.contigNames)[it->contigId] << '\t'
               << (it->isForward ? 'F' : 'R') << '\t'
               << it->firstPos << '\t'
               << it->lastPos << std::endl;
    }
}

#endif  // BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_
