#ifndef BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_
#define BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_

#include <algorithm>

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
    stream << "(id=" << record.id << ", readId=" << record.readId << ", contigId=" << record.contigId << ", distance="
           << record.distance << ", isForward=" << record.isForward << ", firstPos="
           << record.firstPos << ", lastPos=" << record.lastPos << ")";
    return stream;
}


struct WitStore {
    typedef StringSet<CharString, Owner<> > TNameSet;
    typedef String<IntervalOfReadOnContig> TIntervalStore;

    char mateSeparator;

    Holder<TNameSet> readNames;
    Holder<TNameSet> contigNames;

    // TODO(holtgrew): Rename to witRecords.
    TIntervalStore intervals;

    WitStore() {}

    WitStore(TNameSet & _readNames, TNameSet & _contigNames)
            : readNames(_readNames), contigNames(_contigNames) {}

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
struct WitStoreLess<seqan::SortId>
{
    WitStore const & store;

    WitStoreLess(WitStore const & _store) : store(_store) {};

    bool
    operator()(IntervalOfReadOnContig const & a, IntervalOfReadOnContig const & b)  const{
        return a.id < b.id;
    }
};


template <>
struct WitStoreLess<seqan::SortContigId>
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
sortWitRecords(WitStore const & store, TSortTag const &) {
    WitStoreLess<TSortTag const> less(store);
    std::stable_sort(begin(store.intervals, Standard()), end(store.intervals, Standard()), less);
}


inline size_t
appendValue(WitStore & store, IntervalOfReadOnContig const & record) {
    IntervalOfReadOnContig tmp(record);
    tmp.id = length(store.intervals);
    appendValue(store.intervals, tmp);
    return tmp.id;
}


template <typename TFragmentStore>
void loadWitFile(WitStore & store,
                 TFragmentStore /*const*/ & fragments,
                 CharString const & fileName) {
    typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;

    // Assign read and contig names into wit store members.
    setValue(store.readNames, fragments.readNameStore);
    setValue(store.contigNames, fragments.contigNameStore);

    // Build mapping from existing read and contig names to their indices
    // (ids) in readNames and contigNames.
    NameStoreCache<StringSet<CharString> > readNameCache(value(store.readNames));
    refresh(readNameCache);
    NameStoreCache<StringSet<CharString> > contigNameCache(value(store.contigNames));
    refresh(contigNameCache);

    // Then, load the file.
    std::fstream fstrm(toCString(fileName), std::ios_base::in | std::ios_base::binary);
    if (!fstrm.is_open()) {
        std::cerr << "Could not open WIT file " << fileName << std::endl;
        exit(1);
    }

    // Read header.
    char c = '\0';
    // Separator char of read name and mate identifier, index of first mate id.
    char mateSeparator = '/';  // TODO(holtgrew): Un-hardcode.
    int mateStart = 0;  // TODO(holtgrew): Un-hardcode.
    readWitHeader(fstrm, c);

    // Temporary data for one record.
    CharString readName;
    int mateNo;
    size_t distance;
    CharString contigName;
    bool isForward;
    size_t firstPos;
    size_t lastPos;
    bool wasInSam = true;

    // Read WIT file.
    while (!_streamEOF(fstrm)) {
        // Skip comments.
        if (c == '#') {
            _parse_skipLine(fstrm, c);
            continue;
        }

        clear(readName);
        clear(contigName);
        mateNo = -1;

        // Read line.
        _parse_readIdentifier(fstrm, readName, c);
        //std::cout << "readName " << readName << ", " << mateNo << std::endl;
        if (readName[length(readName) - 2] == mateSeparator) {
          mateNo = readName[length(readName) - 1] - '0' - mateStart;
          resize(readName, length(readName) - 2);
        }
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
        //
        // We also need to insert it into the fragment store if it does not
        // exist there yet.
        IntervalOfReadOnContig record;
        if (!getIdByName(value(store.readNames), readName, record.readId, readNameCache)) {
          wasInSam = false;
          record.readId = length(value(store.readNames));
          appendName(value(store.readNames), readName, readNameCache);

          if (mateNo != -1) {
            // If read is paired, create new entry in read name store.
            TMatePairElement mateElem;
            // set the first or second read ID in the mate pair element
            size_t matePairId = length(fragments.matePairStore);
            mateElem.readId[mateNo] = record.readId;
            // get a new mate pair ID and add the new mate pair element
            appendValue(fragments.matePairStore, mateElem);
            // set the new mate pair ID in the read element
            appendRead(fragments, "", matePairId);
            SEQAN_ASSERT_TRUE(getIdByName(value(store.readNames), readName, record.readId, readNameCache));
          }
        } else if (mateNo != -1) {
          // Handle case if we know this read's mate but not the read itself.
          size_t matePairId = fragments.readStore[record.readId].matePairId;
          SEQAN_ASSERT_NEQ(matePairId, TMatePairElement::INVALID_ID);
          record.readId = fragments.matePairStore[matePairId].readId[mateNo];
          if (record.readId == TMatePairElement::INVALID_ID) {
            // create new entry in read and read name store
            // set sequence and mate pair ID in new read store element
            record.readId = appendRead(fragments, "", matePairId);
            // add the identifier to the read name store
            appendName(fragments.readNameStore, readName, fragments.readNameStoreCache);
            // set the ID in the mate pair store
            fragments.matePairStore[matePairId].readId[mateNo] = record.readId;
          }
        }
        // Handle case of mate-paired read.
        //std::cerr << "-  " << readName << ", " << mateNo << std::endl;
        //std::cerr << "-  " << fragments.readNameStore[record.readId] << "/" << getMateNo(fragments, record.readId) << std::endl;
        if (mateNo != -1 && getMateNo(fragments, record.readId) != mateNo)
          record.readId = fragments.matePairStore[fragments.readStore[record.readId].matePairId].readId[mateNo];
        //std::cerr << "+  " << readName << ", " << mateNo << std::endl;
        //std::cerr << "+  " << fragments.readNameStore[record.readId] << "/" << getMateNo(fragments, record.readId) << std::endl;
        // end of "handle case of mate-paired read"
        record.distance = distance;
        if (!getIdByName((store.contigNames), contigName, record.contigId, contigNameCache)) {
          record.contigId = length(value(store.contigNames));
          appendName(value(store.contigNames), contigName, contigNameCache);
        }
        record.isForward = isForward;
        record.firstPos = firstPos;
        record.lastPos = lastPos;
        SEQAN_ASSERT_LEQ(record.firstPos, record.lastPos);
        /*int id = */appendValue(store, record);
        //std::cerr << "   record.id == " << id << std::endl;
    }
}


template <typename TStream>
void writeWitFile(TStream & stream, WitStore const & witStore) {
    writeWitHeader(stream);
    writeWitComment(stream , WIT_COLUMN_NAMES);
    typedef typename WitStore::TIntervalStore TIntervalStore;
    typedef typename Iterator<TIntervalStore, Standard>::Type TIterator;
    for (TIterator it = begin(witStore.intervals, Standard()); it != end(witStore.intervals, Standard()); ++it) {
        for (unsigned i = 0; i < length(value(witStore.readNames)[it->readId]); ++i) {
            char c = value(witStore.readNames)[it->readId][i];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
                break;
            stream << c;
        }
        stream << '\t'
               << it->distance << '\t'
               << value(witStore.contigNames)[it->contigId] << '\t'
               << (it->isForward ? 'F' : 'R') << '\t'
               << it->firstPos << '\t'
               << it->lastPos << std::endl;
    }
}

#endif  // BENCHMARKS_READ_MAPPERS_APPS_WIT_STORE_H_
