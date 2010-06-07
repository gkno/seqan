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

#ifndef SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_
#define SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

/**
.Class.SequenceJournal:Journaled versions of arbitrary underlying sequences.
..signature:SequenceJournal<TSequence, Journal<TJournalSpec> >
 */
//template<typename TSequence, typename TJournalSpec>
//class SequenceJournal;


template<typename _TSequence, typename _TJournalSpec>
class SequenceJournal
{
public:
    typedef _TSequence TSequence;
    typedef _TJournalSpec TJournalSpec;
    typedef typename Value<TSequence>::Type TValue;
    typedef TSequence THost;
    typedef String<TValue, Alloc<> > TInsertionBuffer;
    typedef typename Size<THost>::Type TSize;
    typedef typename Position<THost>::Type TPosition;
    typedef JournalEntry<TSize, TPosition> TJournalEntry;
    typedef JournalTree<TJournalEntry, TJournalSpec> TJournalTree;

    // The underlying, hosting string.
    Holder<THost> _host;
    // A buffer for inserted strings.
    TInsertionBuffer _insertionBuffer;
    // The journal is a binary search tree.
    TJournalTree _journalTree;
    // The journal string's size.
    TSize _length;

	SequenceJournal() {}

	SequenceJournal(TSequence & host)
    {
        setHost(*this, host);
	}
};

// ============================================================================
// Metafunctions
// ============================================================================

/**
.Metafunction.Size:
..param.TValue:Spec.Journal String
 */
template<typename TSequence, typename TJournalSpec>
struct Host<SequenceJournal<TSequence, TJournalSpec> >
{
    typedef TSequence Type;
};

template<typename TSequence, typename TJournalSpec>
struct Host<SequenceJournal<TSequence, TJournalSpec> const>
{
    typedef TSequence const Type;
};

/**
.Metafunction.Size:
..param.TValue:Spec.Journal String
 */
template<typename TSequence, typename TJournalSpec>
struct Size<SequenceJournal<TSequence, TJournalSpec> >
        : public Size<TSequence> {};

template<typename TSequence, typename TJournalSpec>
struct Size<SequenceJournal<TSequence, TJournalSpec> const>
        : public Size<TSequence> {};

/**
.Metafunction.Position:
..param.TValue:Spec.Journal String
 */
template<typename TSequence, typename TJournalSpec>
struct Position<SequenceJournal<TSequence, TJournalSpec> const>
        : public Position<TSequence> {};

template<typename TSequence, typename TJournalSpec>
struct Position<SequenceJournal<TSequence, TJournalSpec> >
        : public Position<TSequence> {};

/**
.Metafunction.Value:
..param.TValue:Spec.Journal String
 */
template<typename TSequence, typename TJournalSpec>
struct Value<SequenceJournal<TSequence, TJournalSpec> const>
        : public Value<TSequence> {};

template<typename TSequence, typename TJournalSpec>
struct Value<SequenceJournal<TSequence, TJournalSpec> >
        : public Value<TSequence> {};

/**
.Metafunction.GetValue:
..param.TValue:Spec.Journal String
 */
template<typename TSequence, typename TJournalSpec>
struct GetValue<SequenceJournal<TSequence, TJournalSpec> const>
        : public GetValue<TSequence> {};
 
template<typename TSequence, typename TJournalSpec>
struct GetValue<SequenceJournal<TSequence, TJournalSpec> >
        : public GetValue<TSequence> {};

/**
.Metafunction.Reference:
..param.TValue:Spec.Journal String
 */
template<typename TSequence, typename TJournalSpec>
struct Reference<SequenceJournal<TSequence, TJournalSpec> const>
        : public Reference<TSequence> {};

template<typename TSequence, typename TJournalSpec>
struct Reference<SequenceJournal<TSequence, TJournalSpec> >
        : public Reference<TSequence> {};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TSequence, typename TJournalSpec>
inline
TStream &
operator<<(TStream & stream, SequenceJournal<TSequence, TJournalSpec> const & sequenceJournal)
{
    SEQAN_CHECKPOINT;
    typedef SequenceJournal<TSequence, TJournalSpec> TSequenceJournal;
    typedef typename TSequenceJournal::TJournalTree TJournalTree;
    typedef typename Iterator<TJournalTree const, Standard>::Type TIterator;

    for (TIterator it = begin(sequenceJournal._journalTree), itend = end(sequenceJournal._journalTree); it != itend; ++it) {
        if (value(it).segmentSource == SOURCE_ORIGINAL) {
            stream << infix(value(sequenceJournal._host), value(it).physicalPosition, value(it).physicalPosition + value(it).length);
        } else {
            SEQAN_ASSERT_EQ(value(it).segmentSource, SOURCE_PATCH);
            stream << infix(sequenceJournal._insertionBuffer, value(it).physicalPosition, value(it).physicalPosition + value(it).length);
        }
    }
    return stream;
}

/**
.Function.setHost:
..param.object.type:Spec.Journal String
..param.host.type:Class.String
*/
template<typename TSequence, typename TJournalSpec>
inline
void setHost(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal, TSequence & str)
{
    SEQAN_CHECKPOINT;
    setValue(sequenceJournal._host, str);
    sequenceJournal._length = length(str);
    reinit(sequenceJournal._journalTree, length(str));
}

/**
.Function.host:
..param.object.type:Spec.Journal String
*/
template<typename TSequence, typename TJournalSpec>
inline
typename Host<SequenceJournal<TSequence, TJournalSpec> >::Type &
host(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal)
{
    SEQAN_CHECKPOINT;
    return value(sequenceJournal._host);
}
template<typename TSequence, typename TJournalSpec>
inline
typename Host<SequenceJournal<TSequence, TJournalSpec> >::Type &
host(SequenceJournal<TSequence, TJournalSpec> const & sequenceJournal)
{
    SEQAN_CHECKPOINT;
    return value(sequenceJournal._host);
}

/**
.Function.clear:
..param.object.type:Spec.Journal String
 */
// TODO(holtgrew): Behaviour is to clear the journal, not the string!
template<typename TSequence, typename TJournalSpec>
inline
void
clear(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal)
{
    SEQAN_CHECKPOINT;
    reinit(sequenceJournal._journalTree, length(host(sequenceJournal)));
}

/**
.Function.flatten:
..summary:Apply the journal to the underlying string, destructively on the underlying string.
..signature:flatten(sequenceJournal)
..param.sequenceJournal:The journal string to flatten.
...type:Spec.Journal String
 */
// TODO(holtgrew): Write me! What about non-destructive version that creates a new copy and sets holder to it?

template<typename TSequence, typename TJournalSpec, typename TPos>
inline
void
erase(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal,
      TPos const & pos,
      TPos const & posEnd)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(static_cast<TPos>(sequenceJournal._length), posEnd - pos);
    sequenceJournal._length -= posEnd - pos;
    recordErase(sequenceJournal._journalTree, pos, posEnd);
    if (length(sequenceJournal._journalTree) == 0)
        clear(sequenceJournal._insertionBuffer);
}

template<typename TSequence, typename TJournalSpec, typename TPos>
inline
void
erase(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal,
      TPos const & pos)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(sequenceJournal._length, 1u);
    erase(sequenceJournal, pos, pos + 1);
}


template<typename TSequence, typename TJournalSpec, typename TPos, typename TString>
inline
void
insert(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal,
       TPos const & pos,
       TString const & seq)
{
    SEQAN_CHECKPOINT;
    sequenceJournal._length += length(seq);
    TPos beginPos = length(sequenceJournal._insertionBuffer);
    append(sequenceJournal._insertionBuffer, seq);
    recordInsertion(sequenceJournal._journalTree, pos, beginPos, length(seq));
}


template<typename TSequence, typename TJournalSpec, typename TPos, typename TValue>
inline
void
insertValue(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal,
            TPos const & pos,
            TValue const & value)
{
    SEQAN_CHECKPOINT;
    TPos beginPos = length(sequenceJournal._insertionBuffer);
    appendValue(sequenceJournal._insertionBuffer, value);
    recordInsertion(sequenceJournal._journalTree, pos, beginPos, 1u);
}


template <typename TSequence, typename TJournalSpec, typename TPos, typename TSequence2>
inline
void
assignInfix(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal,
            TPos const & beginPos,
            TPos const & endPos,
            TSequence2 const & valueString)
{
    SEQAN_CHECKPOINT;
    erase(sequenceJournal, beginPos, endPos);
    insert(sequenceJournal, beginPos, valueString);
}


template <typename TSequence, typename TJournalSpec, typename TPos, typename TValue>
inline
void
assignValue(SequenceJournal<TSequence, TJournalSpec> & sequenceJournal,
            TPos const & pos,
            TValue const & value)
{
    SEQAN_CHECKPOINT;
    erase(sequenceJournal, pos);
    insertValue(sequenceJournal, pos, value);
}


// TODO(holtgrew): Batch-Assignment of values through segments?

// TODO(holtgrew): begin
// TODO(holtgrew): empty
// TODO(holtgrew): end
// TODO(holtgrew): flatten
// TODO(holtgrew): fill
// TODO(holtgrew): getValue
// TODO(holtgrew): infix
// TODO(holtgrew): infixWithLength
// TODO(holtgrew): iter

// TODO(holtgrew): Unused, remove?
/*
template<typename TSequence, typename TJournalSpec>
inline
typename Value<TSequence>::Type const &
front(SequenceJournal<TSequence, TJournalSpec> const & sequenceJournal)
{
    SEQAN_XXXCHECKPOINT;
    typedef SequenceJournal<TSequence, TJournalSpec> TString;
    typedef typename TString::TNode TNode;
    TNode frontNode = front(sequenceJournal._journalTree);
    if (frontNode->segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(sequenceJournal._host), frontNode->virtualPosition + frontNode->length - 1);
    } else {
        SEQAN_ASSERT_EQ(frontNode->segmentSource, SOURCE_PATCH);
        return getValue(sequenceJournal._insertionBuffer, frontNode->virtualPosition + frontNode->length - 1);
    }
}

// front/back clash with general sequence definitions.
template<typename TSequence, typename TJournalSpec>
inline
TValue const &
back(SequenceJournal<TSequence, TJournalSpec> const & sequenceJournal)
{
    SEQAN_XXXCHECKPOINT;
    typedef SequenceJournal<TSequence, TJournalSpec> TString;
    typedef typename TString::TNode TNode;
    TNode backNode = back(sequenceJournal._journalTree);
    if (backNode->segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(sequenceJournal._host), backNode->virtualPosition + backNode->length - 1);
    } else {
        SEQAN_ASSERT_EQ(backNode->segmentSource, SOURCE_PATCH);
        return getValue(sequenceJournal._insertionBuffer, backNode->virtualPosition + backNode->length - 1);
    }
}
*/

template<typename TSequence, typename TJournalSpec>
inline
typename Size<SequenceJournal<TSequence, TJournalSpec> >::Type
length(SequenceJournal<TSequence, TJournalSpec> const & sequenceJournal)
{
    SEQAN_CHECKPOINT;
    return sequenceJournal._length;
}


// TODO(holtgrew): toCString
// TODO(holtgrew): value

// TOOD(holtgrew): operator<
// TOOD(holtgrew): operator>
// TOOD(holtgrew): operator<=
// TOOD(holtgrew): operator>=
// TOOD(holtgrew): operator==
// TOOD(holtgrew): operator!=

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_SEQUENCE_JOURNAL_H_
