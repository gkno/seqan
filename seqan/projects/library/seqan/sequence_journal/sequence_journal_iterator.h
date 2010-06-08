/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
  ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ============================================================================
  Code for the Journaled string iterator.
  ==========================================================================*/

#ifndef SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_
#define SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

template <typename TJournalStringSpec>
struct JournalStringIterSpec;

template <typename TJournalString, typename TJournalSpec>
class Iter<TJournalString, JournalStringIterSpec<TJournalSpec> >
{
public:
    typedef Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > TIterator;
    typedef typename TJournalString::TValue TValue;
    typedef typename JournalType<TJournalString>::Type TJournalEntries;
    // We need a rooted iterator for iterating the journal tree since we need atEnd().
    typedef typename Iterator<TJournalEntries, Rooted>::Type TJournalEntriesIterator;
    typedef typename Host<TJournalString>::Type THost;
    typedef typename Iterator<THost, Standard>::Type THostIterator;
    typedef typename TJournalString::TInsertionBuffer TInsertionBuffer;
    typedef typename Iterator<TInsertionBuffer, Standard>::Type TInsertionBufferIterator;

    // The journal string we iterate over.
    TJournalString * _journalStringPtr;
    // Iterator over the segments in the journal tree.
    TJournalEntriesIterator _journalEntriesIterator;
    // Begin and end iterator in the host string of the journal string.
    THostIterator _hostSegmentBegin;
    THostIterator _hostSegmentEnd;
    // Current iterator in the host segment.
    THostIterator _currentHostIt;
    // Begin and end iterator in the insertion buffer of the journal string.
    TInsertionBufferIterator _insertionBufferSegmentBegin;
    TInsertionBufferIterator _insertionBufferSegmentEnd;
    // Current iterator in the insertion buffer.
    TInsertionBufferIterator _currentInsertionBufferIt;

    Iter() : _journalStringPtr(0) { SEQAN_CHECKPOINT; }

    Iter(TIterator const & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalEntriesIterator(other._journalEntriesIterator),
              _hostSegmentBegin(other._hostSegmentBegin),
              _hostSegmentEnd(other._hostSegmentEnd),
              _currentHostIt(other._currentHostIt),
              _insertionBufferSegmentBegin(other._insertionBufferSegmentBegin),
              _insertionBufferSegmentEnd(other._insertionBufferSegmentEnd),
              _currentInsertionBufferIt(other._currentInsertionBufferIt)
    {
        SEQAN_CHECKPOINT;
    }

    Iter(typename IterComplementConst<TIterator>::Type const & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalEntriesIterator(other._journalEntriesIterator),
              _hostSegmentBegin(other._hostSegmentBegin),
              _hostSegmentEnd(other._hostSegmentEnd),
              _currentHostIt(other._currentHostIt),
              _insertionBufferSegmentBegin(other._insertionBufferSegmentBegin),
              _insertionBufferSegmentEnd(other._insertionBufferSegmentEnd),
              _currentInsertionBufferIt(other._currentInsertionBufferIt)
    {
        SEQAN_CHECKPOINT;
    }

    explicit
    Iter(TJournalString & journalString)
    {
        SEQAN_CHECKPOINT;
        _initJournalStringIterator(*this, journalString);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// For SequenceJournal<TSequence, TJournalSpec>

/**
.Metafunction.Iterator:
..param.TValue:Spec.Journal String
 */
template <typename TSequence, typename TJournalSpec>
struct Iterator<SequenceJournal<TSequence, TJournalSpec>, Standard>
{
    typedef Iter<SequenceJournal<TSequence, TJournalSpec>, JournalStringIterSpec<TJournalSpec> > Type;
};

template <typename TSequence, typename TJournalSpec>
struct Iterator<SequenceJournal<TSequence, TJournalSpec> const, Standard>
{
    typedef Iter<SequenceJournal<TSequence, TJournalSpec> const, JournalStringIterSpec<TJournalSpec> > Type;
};

// For Iter<TJournalString, JournalStringIterSpec>

// ============================================================================
// Functions
// ============================================================================

// For SequenceJournal<TSequence, TJournalSpec>

template <typename TSequence, typename TJournalSpec>
inline
typename Iterator<SequenceJournal<TSequence, TJournalSpec> const, Standard>::Type
begin(SequenceJournal<TSequence, TJournalSpec> const & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<SequenceJournal<TSequence, TJournalSpec> const, Standard>::Type TResult;
    return TResult(journalString);
}

template <typename TSequence, typename TJournalSpec>
inline
typename Iterator<SequenceJournal<TSequence, TJournalSpec>, Standard>::Type
begin(SequenceJournal<TSequence, TJournalSpec> & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<SequenceJournal<TSequence, TJournalSpec>, Standard>::Type TResult;
    return TResult(journalString);
}

template <typename TSequence, typename TJournalSpec>
inline
typename Iterator<SequenceJournal<TSequence, TJournalSpec> const, Standard>::Type
end(SequenceJournal<TSequence, TJournalSpec> const & journalString, Standard)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<SequenceJournal<TSequence, TJournalSpec> const, Standard>::Type TResult;
    TResult result;
    _initJournalStringIteratorEnd(result, journalString);
    return result;
}

template <typename TSequence, typename TJournalSpec>
inline
typename Iterator<SequenceJournal<TSequence, TJournalSpec>, Standard>::Type
end(SequenceJournal<TSequence, TJournalSpec> & journalString, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<SequenceJournal<TSequence, TJournalSpec>, Standard>::Type TResult;
    TResult result;
    _initJournalStringIteratorEnd(result, journalString);
    return result;
}

// For Iter<TJournalString, JournalStringIterSpec>

template <typename TJournalString, typename TJournalSpec>
inline
void
_initJournalStringIterator(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > & iterator,
                           TJournalString & journalString)
{
    SEQAN_CHECKPOINT;
    iterator._journalStringPtr = &journalString;
    iterator._journalEntriesIterator = begin(journalString._journalEntries);
    // Update iterators on the segment.
    _updateSegmentIterators(iterator);
}

template <typename TJournalString, typename TJournalSpec>
inline
void
_initJournalStringIteratorEnd(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > & iterator,
                              TJournalString & journalString)
{
    SEQAN_CHECKPOINT;
    iterator._journalStringPtr = &journalString;
    iterator._journalEntriesIterator = end(journalString._journalEntries);
}

template <typename TJournalString, typename TJournalSpec>
inline
void
_updateSegmentIterators(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (atEnd(iterator._journalEntriesIterator))
        return;
    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
//             static_cast<int>(begin(host(*iterator._journalStringPtr), Standard()));
//             static_cast<int>(host(*iterator._journalStringPtr));
//             static_cast<int>(*iterator._journalStringPtr);
//             static_cast<int>(iterator._journalStringPtr);
//             static_cast<int>(iterator._hostSegmentBegin);
            iterator._hostSegmentBegin = begin(host(*iterator._journalStringPtr), Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
            iterator._hostSegmentEnd = iterator._hostSegmentBegin + value(iterator._journalEntriesIterator).length;
            iterator._currentHostIt = iterator._hostSegmentBegin;
            break;
        case SOURCE_PATCH:
            iterator._insertionBufferSegmentBegin = begin(iterator._journalStringPtr->_insertionBuffer, Standard()) + value(iterator._journalEntriesIterator).physicalPosition;
            iterator._insertionBufferSegmentEnd = iterator._insertionBufferSegmentBegin + value(iterator._journalEntriesIterator).length;
            iterator._currentInsertionBufferIt = iterator._insertionBufferSegmentBegin;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
}

template <typename TJournalString, typename TJournalSpec>
inline
typename Value<TJournalString>::Type
value(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        return value(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(value(iterator._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        return value(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournalString, typename TJournalSpec>
inline
typename Value<TJournalString>::Type
value(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & iterator)
{
    SEQAN_CHECKPOINT;
    if (value(iterator._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        return value(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(value(iterator._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        return value(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournalString, typename TJournalSpec>
inline
Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > &
operator++(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    switch (value(iterator._journalEntriesIterator).segmentSource) {
        case SOURCE_ORIGINAL:
            ++iterator._currentHostIt;
            if (iterator._currentHostIt == iterator._hostSegmentEnd) {
                ++iterator._journalEntriesIterator;
                _updateSegmentIterators(iterator);
            }
            break;
        case SOURCE_PATCH:
            ++iterator._currentInsertionBufferIt;
            if (iterator._currentInsertionBufferIt == iterator._insertionBufferSegmentEnd) {
                ++iterator._journalEntriesIterator;
                _updateSegmentIterators(iterator);
            }
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
    return iterator;
}

template <typename TJournalString, typename TJournalSpec>
inline
Iter<TJournalString, JournalStringIterSpec<TJournalSpec> >
operator++(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > & iterator, int /*postfix*/)
{
    SEQAN_CHECKPOINT;
    Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > temp(iterator);
    ++iterator;
    return temp;    
}

template <typename TJournalString, typename TJournalSpec>
inline
typename Value<TJournalString>::Type
operator*(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > & iterator)
{
    SEQAN_CHECKPOINT;
    return value(iterator);
}

template <typename TJournalString, typename TJournalSpec>
inline
typename Value<TJournalString>::Type
operator*(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & iterator)
{
    SEQAN_CHECKPOINT;
    return value(iterator);
}

template <typename TJournalString, typename TJournalSpec>
inline
Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > &
operator+=(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > & iterator,
          typename Size<TJournalString>::Type len)
{
    SEQAN_CHECKPOINT;
    typedef typename Size<TJournalString>::Type TSize;
    while (len > 0) {
        TSize remaining;
        switch (value(iterator._journalEntriesIterator).segmentSource) {
            case SOURCE_ORIGINAL:
                remaining = iterator._hostSegmentEnd - iterator._currentHostIt;
                SEQAN_ASSERT_GT(remaining, 0u);
                if (len >= remaining) {
                    len -= remaining;
                    ++iterator._journalEntriesIterator;
                    _updateSegmentIterators(iterator);
                } else {
                    iterator._currentHostIt += len;
                    len = 0;
                }
                break;
            case SOURCE_PATCH:
                remaining = iterator._insertionBufferSegmentEnd - iterator._currentInsertionBufferIt;
                SEQAN_ASSERT_GT(remaining, 0u);
                if (len >= remaining) {
                    len -= remaining;
                    ++iterator._journalEntriesIterator;
                    _updateSegmentIterators(iterator);
                } else {
                    iterator._currentInsertionBufferIt += len;
                    len = 0;
                }
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid segment source!");
        }
    }
    return iterator;
}

template <typename TJournalString, typename TJournalSpec>
inline
Iter<TJournalString, JournalStringIterSpec<TJournalSpec> >
operator+(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & iterator,
          typename Size<TJournalString>::Type const & len)
{
    SEQAN_CHECKPOINT;
    Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > temp(iterator);
    temp += len;
    return temp;
}

template <typename TJournalString, typename TJournalSpec>
inline
bool
operator==(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    if (atEnd(a._journalEntriesIterator) && atEnd(b._journalEntriesIterator))
        return true;
    if (a._journalEntriesIterator != b._journalEntriesIterator)
        return false;
    if (value(a._journalEntriesIterator).segmentSource == SOURCE_ORIGINAL) {
        if (a._currentHostIt != b._currentHostIt)
            return false;
    } else {
        SEQAN_ASSERT_EQ(value(a._journalEntriesIterator).segmentSource, SOURCE_PATCH);
        if (a._currentInsertionBufferIt != b._currentInsertionBufferIt)
            return false;
    }
    return true;
}

template <typename TJournalString, typename TJournalSpec>
inline
bool
operator==(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    typedef typename IterMakeConst<Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > >::Type TConstIter;
    return static_cast<TConstIter>(a) == static_cast<TConstIter>(b);
}

template <typename TJournalString, typename TJournalSpec>
inline
bool
operator!=(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & a,
           Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & b)
{
    SEQAN_CHECKPOINT;
    return !(a == b);
}

template <typename TJournalString, typename TJournalSpec>
inline
bool
operator!=(Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > const & a,
           typename IterComplementConst<Iter<TJournalString, JournalStringIterSpec<TJournalSpec> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return !(a == b);
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_
