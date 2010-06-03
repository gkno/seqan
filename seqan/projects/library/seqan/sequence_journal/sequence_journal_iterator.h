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
  ==========================================================================*/

#ifndef SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_
#define SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

struct JournalStringIterSpec;

template <typename TJournalString>
class Iter<TJournalString, JournalStringIterSpec>
{
public:
    typedef Iter<TJournalString, JournalStringIterSpec> TIterator;
    typedef typename TJournalString::TValue TValue;
    typedef typename TJournalString::TJournalSpec TJournalSpec;
    typedef typename TJournalString::TJournalTree TJournalTree;
    typedef typename Iterator<TJournalTree, Standard>::Type TJournalTreeIterator;
    typedef typename TJournalString::THost THost;
    typedef typename Iterator<THost, Standard>::Type THostIterator;
    typedef typename TJournalString::TInsertionBuffer TInsertionBuffer;
    typedef typename Iterator<TInsertionBuffer, Standard>::Type TInsertionBufferIterator;

    // The journal string we iterate over.
    TJournalString * _journalStringPtr;
    // Iterator over the segments in the journal tree.
    TJournalTreeIterator _journalTreeIterator;
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
              _journalTreeIterator(other._journalTreeIterator),
              _hostSegmentBegin(other._hostSegmentBegin),
              _hostSegmentEnd(other._hostSegmentEnd),
              _currentHostIt(other._currentHostIt),
              _insertionBufferSegmentBegin(other._insertionBufferSegmentBegin),
              _insertionBufferSegmentEnd(other._insertionBufferSegmentEnd),
              _currentInsertionBufferIt(other._currentInsertionBufferIt)
    { SEQAN_CHECKPOINT; }

    // Always allow conversion from non-const.
    Iter(typename _RemoveConst<TIterator>::Type & other)
            : _journalStringPtr(other._journalStringPtr),
              _journalTreeIterator(other._journalTreeIterator),
              _hostSegmentBegin(other._hostSegmentBegin),
              _hostSegmentEnd(other._hostSegmentEnd),
              _currentHostIt(other._currentHostIt),
              _insertionBufferSegmentBegin(other._insertionBufferSegmentBegin),
              _insertionBufferSegmentEnd(other._insertionBufferSegmentEnd),
              _currentInsertionBufferIt(other._currentInsertionBufferIt)
    { SEQAN_CHECKPOINT; }

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

// For SequenceJournal<TSequence, Journal<TJournalSpec> >

/**
.Metafunction.Iterator:
..param.TValue:Spec.Journal String
 */
template <typename TSequence, typename TJournalSpec>
struct Iterator<SequenceJournal<TSequence, TJournalSpec>, Standard>
{
    typedef Iter<SequenceJournal<TSequence, TJournalSpec>, JournalStringIterSpec> Type;
};

template <typename TSequence, typename TJournalSpec>
struct Iterator<SequenceJournal<TSequence, TJournalSpec> const, Standard>
{
    typedef Iter<SequenceJournal<TSequence, TJournalSpec> const, JournalStringIterSpec> Type;
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
end(SequenceJournal<TSequence, TJournalSpec> const & journalString, Standard const &)
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
    
template <typename TJournalString>
inline
void
_initJournalStringIterator(Iter<TJournalString, JournalStringIterSpec> & iterator,
                           TJournalString & journalString)
{
    SEQAN_CHECKPOINT;
    iterator._journalStringPtr = &journalString;
    // Initialize iterator on the journal tree, will point to leftmost child,
    // then.
    _initJournalTreeIterator(iterator._journalTreeIterator, journalString._journalTree);
    // Update iterators on the segment.
    _updateSegmentIterators(iterator);
}

template <typename TJournalString>
inline
void
_initJournalStringIteratorEnd(Iter<TJournalString, JournalStringIterSpec> & iterator,
                              TJournalString & journalString)
{
    SEQAN_CHECKPOINT;
    iterator._journalStringPtr = &journalString;
    // Initialize iterator on the journal tree to the end.
    _initJournalTreeIteratorEnd(iterator._journalTreeIterator, journalString._journalTree);
}

template <typename TJournalString>
inline
void
_updateSegmentIterators(Iter<TJournalString, JournalStringIterSpec> & iterator)
{
    SEQAN_CHECKPOINT;
    if (atEnd(iterator._journalTreeIterator))
        return;
    switch (iterator._journalTreeIterator._currentNode->segmentSource) {
        case SOURCE_ORIGINAL:
            iterator._hostSegmentBegin = begin(host(value(iterator._journalStringPtr)), Standard()) + iterator._journalTreeIterator._currentNode->physicalPosition;
            iterator._hostSegmentEnd = iterator._hostSegmentBegin + iterator._journalTreeIterator._currentNode->length;
            iterator._currentHostIt = iterator._hostSegmentBegin;
            break;
        case SOURCE_PATCH:
            iterator._insertionBufferSegmentBegin = begin(value(iterator._journalStringPtr)._insertionBuffer, Standard()) + iterator._journalTreeIterator._currentNode->physicalPosition;
            iterator._insertionBufferSegmentEnd = iterator._insertionBufferSegmentBegin + iterator._journalTreeIterator._currentNode->length;
            iterator._currentInsertionBufferIt = iterator._insertionBufferSegmentBegin;
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
}

template <typename TJournalString>
inline
typename Value<TJournalString>::Type
value(Iter<TJournalString, JournalStringIterSpec> & iterator)
{
    SEQAN_CHECKPOINT;
    if (iterator._journalTreeIterator._currentNode->segmentSource == SOURCE_ORIGINAL) {
        return value(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(iterator._journalTreeIterator._currentNode->segmentSource, SOURCE_PATCH);
        return value(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournalString>
inline
typename Value<TJournalString>::Type
value(Iter<TJournalString, JournalStringIterSpec> const & iterator)
{
    SEQAN_CHECKPOINT;
    if (iterator._journalTreeIterator._currentNode->segmentSource == SOURCE_ORIGINAL) {
        return value(iterator._currentHostIt);
    } else {
        SEQAN_ASSERT_EQ(iterator._journalTreeIterator._currentNode->segmentSource, SOURCE_PATCH);
        return value(iterator._currentInsertionBufferIt);
    }
}

template <typename TJournalString>
inline
Iter<TJournalString, JournalStringIterSpec> &
operator++(Iter<TJournalString, JournalStringIterSpec> & iterator)
{
    SEQAN_CHECKPOINT;
    switch (iterator._journalTreeIterator._currentNode->segmentSource) {
        case SOURCE_ORIGINAL:
            ++iterator._currentHostIt;
            if (iterator._currentHostIt == iterator._hostSegmentEnd) {
                ++iterator._journalTreeIterator;
                _updateSegmentIterators(iterator);
            }
            break;
        case SOURCE_PATCH:
            ++iterator._currentInsertionBufferIt;
            if (iterator._currentInsertionBufferIt == iterator._insertionBufferSegmentEnd) {
                ++iterator._journalTreeIterator;
                _updateSegmentIterators(iterator);
            }
            break;
        default:
            SEQAN_ASSERT_FAIL("Invalid segment source!");
    }
    return iterator;
}

template <typename TJournalString>
inline
Iter<TJournalString, JournalStringIterSpec>
operator++(Iter<TJournalString, JournalStringIterSpec> & iterator, int postfix)
{
    SEQAN_CHECKPOINT;
    typename Iterator<TJournalString, JournalStringIterSpec>::Type temp;
    ++iterator;
    return temp;    
}

template <typename TJournalString>
inline
typename Value<TJournalString>::Type
operator*(Iter<TJournalString, JournalStringIterSpec> & iterator)
{
    SEQAN_CHECKPOINT;
    return value(iterator);
}

template <typename TJournalString>
inline
typename Value<TJournalString>::Type
operator*(Iter<TJournalString, JournalStringIterSpec> const & iterator)
{
    SEQAN_CHECKPOINT;
    return value(iterator);
}

template <typename TJournalString>
inline
Iter<TJournalString, JournalStringIterSpec> &
operator+=(Iter<TJournalString, JournalStringIterSpec> & iterator,
          typename Size<TJournalString>::Type len)
{
    SEQAN_CHECKPOINT;
    typedef typename Size<TJournalString>::Type TSize;
    while (len > 0) {
        TSize remaining;
        switch (iterator._journalTreeIterator._currentNode->segmentSource) {
            case SOURCE_ORIGINAL:
                remaining = iterator._hostSegmentEnd - iterator._currentHostIt;
                SEQAN_ASSERT_GT(remaining, 0u);
                if (len >= remaining) {
                    len -= remaining;
                    ++iterator._journalTreeIterator;
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
                    ++iterator._journalTreeIterator;
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

template <typename TJournalString>
inline
Iter<TJournalString, JournalStringIterSpec>
operator+(Iter<TJournalString, JournalStringIterSpec> const & iterator,
          typename Size<TJournalString>::Type const & len)
{
    SEQAN_CHECKPOINT;
    Iter<TJournalString, JournalStringIterSpec> temp(iterator);
    temp += len;
    return temp;
}

template <typename TJournalString>
inline
bool
operator==(Iter<TJournalString, JournalStringIterSpec> const & a,
           Iter<TJournalString, JournalStringIterSpec> const & b)
{
    SEQAN_CHECKPOINT;
    if (atEnd(a._journalTreeIterator) && atEnd(b._journalTreeIterator))
        return true;
    if (a._journalTreeIterator._currentNode != b._journalTreeIterator._currentNode)
        return false;
    if (a._journalTreeIterator._iterationDirection != b._journalTreeIterator._iterationDirection)
        return false;
    if (a._journalTreeIterator._currentNode->segmentSource == SOURCE_ORIGINAL) {
        if (a._currentHostIt != b._currentHostIt)
            return false;
    } else {
        SEQAN_ASSERT_EQ(a._journalTreeIterator._currentNode->segmentSource, SOURCE_PATCH);
        if (a._currentInsertionBufferIt != b._currentInsertionBufferIt)
            return false;
    }
    return true;
}

template <typename TJournalString>
inline
bool
operator!=(Iter<TJournalString, JournalStringIterSpec> const & a,
           Iter<TJournalString, JournalStringIterSpec> const & b)
{
    SEQAN_CHECKPOINT;
    return !(a == b);
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_ITERATOR_H_
