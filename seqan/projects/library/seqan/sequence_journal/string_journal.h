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

#ifndef SEQAN_SEQUENCE_JOURNAL_STRING_JOURNAL_H_
#define SEQAN_SEQUENCE_JOURNAL_STRING_JOURNAL_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

/**
.Tag.Journal
..signature:Journal<TStringSpec, TJournalSpec>
 */
template <typename TStringSpec, typename TJournalSpec>
struct Journal;


/**
.Spec.Journal String:Journaled versions of arbitrary underlying strings.
..signature:String<TValue, Journal<TStringSpec, TJournalSpec>
 */
// TODO(holtgrew): What about SequenceJournal<TSequence, TJournalSpec>?
template<typename _TValue, typename _TStringSpec, typename _TJournalSpec>
class String<_TValue, Journal<_TStringSpec, _TJournalSpec > >
{
public:
    typedef _TValue TValue;
    typedef _TStringSpec TStringSpec;
    typedef _TJournalSpec TJournalSpec;
    typedef String<TValue, TStringSpec> THost;
    typedef String<TValue, Alloc<> > TInsertionBuffer;
    typedef typename Size<THost>::Type TSize;
    typedef typename Position<THost>::Type TPosition;
    typedef SegmentNode<TSize, TPosition> TNode;
    typedef JournalTree<TNode, TJournalSpec> TJournalTree;
    
	// The underlying, hosting string.
	Holder<THost> _host;
    // A buffer for inserted strings.
    TInsertionBuffer _insertionBuffer;
    // The journal is a binary search tree.
    TJournalTree _journalTree;
    // The journal string's size.
    TSize _length;
    
	String() {}

	String(String<TValue, TStringSpec> & host)
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
template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct Size<String<TValue, Journal<TStringSpec, TJournalSpec> > >
        : public Size<String<TValue, TStringSpec> > {};

template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct Size<String<TValue, Journal<TStringSpec, TJournalSpec> > const>
        : public Size<String<TValue, TStringSpec> const> {};

/**
.Metafunction.Position:
..param.TValue:Spec.Journal String
 */
template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct Position<String<TValue, Journal<TStringSpec, TJournalSpec> > >
        : public Position<String<TValue, TStringSpec> > {};

template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct Position<String<TValue, Journal<TStringSpec, TJournalSpec> > const>
        : public Position<String<TValue, TStringSpec> const> {};

/**
.Metafunction.Value:
..param.TValue:Spec.Journal String
 */
template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct Value<String<TValue, Journal<TStringSpec, TJournalSpec> > >
        : public Value<String<TValue, TStringSpec> > {};

template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct Value<String<TValue, Journal<TStringSpec, TJournalSpec> > const>
        : public Value<String<TValue, TStringSpec> const> {};

/**
.Metafunction.GetValue:
..param.TValue:Spec.Journal String
 */
template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct GetValue<String<TValue, Journal<TStringSpec, TJournalSpec> > >
        : public GetValue<String<TValue, TStringSpec> > {};
 
template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct GetValue<String<TValue, Journal<TStringSpec, TJournalSpec> > const>
        : public GetValue<String<TValue, TStringSpec> const> {};

/**
.Metafunction.Reference:
..param.TValue:Spec.Journal String
 */
template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct Reference<String<TValue, Journal<TStringSpec, TJournalSpec> > >
        : public Value<String<TValue, TStringSpec> > {};

template<typename TValue, typename TStringSpec, typename TJournalSpec>
struct Reference<String<TValue, Journal<TStringSpec, TJournalSpec> > const>
        : public Reference<String<TValue, TStringSpec> const> {};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TValue, typename TStringSpec, typename TJournalSpec>
inline
TStream &
operator<<(TStream & stream, String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString)
{
    typedef String<TValue, Journal<TStringSpec, TJournalSpec> > TSequenceJournal;
    typedef typename TSequenceJournal::TJournalTree TJournalTree;
    typedef typename Iterator<TJournalTree>::Type TIterator;

//     std::cout << __FILE__ << ":" << __LINE__ << " -- " << journalString._journalTree << std::endl;
    TIterator itend = end(journalString._journalTree, Standard());
//     std::cout << __FILE__ << ":" << __LINE__ << " -- itend = " << value(itend) << ", " << itend._iterationDirection << std::endl;
    for (TIterator it = begin(journalString._journalTree, Standard()); it != itend; ++it) {
//         std::cout << __FILE__ << ":" << __LINE__ << " --   " << value(it) << ", " << it._iterationDirection << std::endl;
        if (value(it)->segmentSource == SOURCE_ORIGINAL) {
            stream << infix(value(journalString._host), value(it)->physicalPosition, value(it)->physicalPosition + value(it)->length);
        } else {
            SEQAN_ASSERT_EQ(value(it)->segmentSource, SOURCE_PATCH);
            stream << infix(journalString._insertionBuffer, value(it)->physicalPosition, value(it)->physicalPosition + value(it)->length);
        }
    }
    return stream;
    /*

    // Depth-first, in-order traversal of the tree using an explicit stack.
    String<TNode*> nodePointers;
    TNode * current = journalString._journalTree._root;
    while (true) {
        while (current) {
            appendValue(nodePointers, current);
            current = current->left;
        }

        // POP - or break.
        if (empty(nodePointers))
            break;
        current = back(nodePointers);
        eraseBack(nodePointers);

        if (current->segmentSource == SOURCE_ORIGINAL) {
//             stream << "((" << infix(value(journalString._host), current->physicalPosition, current->physicalPosition + current->length) << "))";
            stream << infix(value(journalString._host), current->physicalPosition, current->physicalPosition + current->length);
        } else {
            SEQAN_ASSERT_EQ(current->segmentSource, SOURCE_PATCH);
//             stream << "((" << infix(journalString._insertionBuffer, current->physicalPosition, current->physicalPosition + current->length) << "))";
            stream << infix(journalString._insertionBuffer, current->physicalPosition, current->physicalPosition + current->length);
        }
        current = current->right;
    }

    return stream;
    */
}


/**
.Function.setHost:
..param.object.type:Spec.Journal String
..param.host.type:Class.String
*/
template<typename TValue, typename TStringSpec, typename TJournalSpec>
inline
void setHost(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString, String<TValue, TStringSpec> & str)
{
    SEQAN_CHECKPOINT;
    setValue(journalString._host, str);
    journalString._length = length(str);
    reinit(journalString._journalTree, length(str));
}

/**
.Function.host:
..param.object.type:Spec.Journal String
*/
template<typename TValue, typename TStringSpec, typename TJournalSpec>
inline
String<TValue, TStringSpec> &
host(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString)
{
    SEQAN_CHECKPOINT;
    return value(journalString._host);
}

/**
.Function.clear:
..param.object.type:Spec.Journal String
 */
// TODO(holtgrew): Behaviour is to clear the journal, not the string!
template<typename TValue, typename TStringSpec, typename TJournalSpec>
inline
void
clear(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString)
{
    SEQAN_CHECKPOINT;
    reinit(journalString._journalTree, length(host(journalString)));
}

/**
.Function.flatten:
..summary:Apply the journal to the underlying string, destructively on the underlying string.
..signature:flatten(journalString)
..param.journalString:The journal string to flatten.
...type:Spec.Journal String
 */
// TODO(holtgrew): Write me! What about non-destructive version that creates a new copy and sets holder to it?

template<typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos>
inline
void
erase(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString,
      TPos const & pos,
      TPos const & posEnd)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(static_cast<TPos>(journalString._length), posEnd - pos);
    journalString._length -= posEnd - pos;
    recordErase(journalString._journalTree, pos, posEnd);
    if (length(journalString) == 0)
        clear(journalString._insertionBuffer);
}

template<typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos>
inline
void
erase(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString,
      TPos const & pos)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ(journalString._length, 1u);
    erase(journalString, pos, pos + 1);
}


template<typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos, typename TString>
inline
void
insert(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString,
       TPos const & pos,
       TString const & seq)
{
    SEQAN_CHECKPOINT;
    journalString._length += length(seq);
    TPos beginPos = length(journalString._insertionBuffer);
    append(journalString._insertionBuffer, seq);
    recordInsertion(journalString._journalTree, pos, beginPos, length(seq));
}


template<typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos>
inline
void
insertValue(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString,
            TPos const & pos,
            TValue const & value)
{
    SEQAN_CHECKPOINT;
    TPos beginPos = length(journalString._insertionBuffer);
    appendValue(journalString._insertionBuffer, value);
    recordInsertion(journalString._journalTree, pos, beginPos, 1u);
}


template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos, typename TString>
inline
void
assignInfix(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString,
            TPos const & beginPos,
            TPos const & endPos,
            TString const & valueString)
{
    SEQAN_CHECKPOINT;
    erase(journalString, beginPos, endPos);
    insert(journalString, beginPos, valueString);
}


template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos>
inline
void
assignValue(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString,
            TPos const & pos,
            TValue const & value)
{
    erase(journalString, pos);
    insertValue(journalString, pos, value);
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

template<typename TValue, typename TStringSpec, typename TJournalSpec>
inline
TValue const &
front(String<TValue, Journal<TStringSpec, TJournalSpec> > const & journalString)
{
    typedef String<TValue, Journal<TStringSpec, TJournalSpec> > TString;
    typedef typename TString::TNode TNode;
    TNode frontNode = front(journalString._journalTree);
    if (frontNode->segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(journalString._host), frontNode->virtualPosition + frontNode->length - 1);
    } else {
        SEQAN_ASSERT_EQ(frontNode->segmentSource, SOURCE_PATCH);
        return getValue(journalString._insertionBuffer, frontNode->virtualPosition + frontNode->length - 1);
    }
}

/* front/back clash with general sequence definitions.
template<typename TValue, typename TStringSpec, typename TJournalSpec>
inline
TValue const &
back(String<TValue, Journal<TStringSpec, TJournalSpec> > const & journalString)
{
    typedef String<TValue, Journal<TStringSpec, TJournalSpec> > TString;
    typedef typename TString::TNode TNode;
    TNode backNode = back(journalString._journalTree);
    if (backNode->segmentSource == SOURCE_ORIGINAL) {
        return getValue(value(journalString._host), backNode->virtualPosition + backNode->length - 1);
    } else {
        SEQAN_ASSERT_EQ(backNode->segmentSource, SOURCE_PATCH);
        return getValue(journalString._insertionBuffer, backNode->virtualPosition + backNode->length - 1);
    }
}
*/

template<typename TValue, typename TStringSpec, typename TJournalSpec>
inline
typename Size<String<TValue, Journal<TStringSpec, TJournalSpec> > >::Type
length(String<TValue, Journal<TStringSpec, TJournalSpec> > const & journalString)
{
    return journalString._length;
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

#endif  // SEQAN_SEQUENCE_JOURNAL_STRING_JOURNAL_H_
