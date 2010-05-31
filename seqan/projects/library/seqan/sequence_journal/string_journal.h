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
  Author: Paul Theodor Pyl <paul.pyl@fu-berlin.de>
  ==========================================================================*/

#ifndef SEQAN_SEQUENCE_JOURNAL_STRING_JOURNAL_H_
#define SEQAN_SEQUENCE_JOURNAL_STRING_JOURNAL_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

/**
.Tag.Journal
..signature:Journal<TStringSpec, TJournalSpec>
 */
template <typename TStringSpec, typename TJournalSpec>
struct Journal;

// ============================================================================
// Specialization Journal String
// ============================================================================

/**
.Spec.Journal String:Journaled versions of arbitrary underlying strings.
..signature:String<TValue, Journal<TStringSpec, TJournalSpec>
 */
template<typename TValue, typename TStringSpec, typename TJournalSpec>
class String<TValue, Journal< TStringSpec, TJournalSpec > >
{
public:
    typedef String<TValue, TStringSpec> THost;
    typedef typename Size<THost>::Type TSize;
    typedef typename Position<THost>::Type TPosition;
    typedef SegmentNode<TSize, TPosition> TNode;
    
    // TODO(holtgrew): Why this?
    // typename DefaultStringSpec< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type m_insertion_string;

	// The underlying, hosting string.
	Holder<THost> _host;
    // A buffer for inserted strings.
    THost _insertionBuffer;
    // The journal is a binary search tree.
    JournalTree<TNode, TJournalSpec> _journalTree;
    
	String() {}

	String(String<TValue, TStringSpec> & host) {
        setHost(*this, host);
	}
};

// ============================================================================
// Metafunctions
// ============================================================================

// TODO(holtgrew): Iterator<>, Iterator<const>

/**
.Metafunction.Size:
..param.TValue:Spec.Journal String
 */
template<typename TValue, typename TStringSpec, typename TJournalSpec>
class Size<String<TValue, Journal< TStringSpec, TJournalSpec > > >
        : public Size<String<TValue, TStringSpec> > {};

template<typename TValue, typename TStringSpec, typename TJournalSpec>
class Size<String<TValue, Journal< TStringSpec, TJournalSpec > > const>
        : public Size<String<TValue, TStringSpec> const> {};

/**
.Metafunction.Position:
..param.TValue:Spec.Journal String
 */
template<typename TValue, typename TStringSpec, typename TJournalSpec>
class Position<String<TValue, Journal< TStringSpec, TJournalSpec > > >
        : public Position<String<TValue, TStringSpec> > {};

template<typename TValue, typename TStringSpec, typename TJournalSpec>
class Position<String<TValue, Journal< TStringSpec, TJournalSpec > > const>
        : public Position<String<TValue, TStringSpec> const> {};

// TODO(holtgrew): Value<>, Value<const>

/**
.Metafunction.GetValue:
..param.TValue:Spec.Journal String
 */
template<typename TValue, typename TStringSpec, typename TJournalSpec>
class GetValue<String<TValue, Journal< TStringSpec, TJournalSpec > > >
        : public GetValue<String<TValue, TStringSpec> > {};
 
template<typename TValue, typename TStringSpec, typename TJournalSpec>
class GetValue<String<TValue, Journal< TStringSpec, TJournalSpec > > const>
        : public GetValue<String<TValue, TStringSpec> const> {};

// TODO(holtgrew): Reference<>, Reference<const>

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TValue, typename TStringSpec, typename TJournalSpec>
inline
TStream &
operator<<(TStream & stream, String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString)
{
    typedef String<TValue, Journal<TStringSpec, TJournalSpec> > TString;
    typedef typename Position<TString>::Type TPos;
    typedef typename TString::TNode TNode;
    String<TNode*> nodePointers;
    appendValue(nodePointers, journalString._journalTree._root);

    // Depth-first, in-order traversal of the tree using an explicit stack.
    while (!empty(nodePointers)) {
        TNode * current = back(nodePointers);
        eraseBack(nodePointers);

        if (current->segmentSource == SOURCE_ORIGINAL) {
            for (TPos i = current->physicalPosition; i < current->physicalPosition + current->length; ++i) {
                stream << getValue(value(journalString._host), i);
            }
        } else {
            SEQAN_ASSERT_EQ(current->segmentSource, SOURCE_PATCH);
            for (TPos i = current->physicalPosition; i < current->physicalPosition + current->length; ++i) {
                stream << getValue(journalString._insertionBuffer, i);
            }
        }

        if (current->left != 0)
            appendValue(nodePointers, current->left);
        if (current->right != 0)
            appendValue(nodePointers, current->right);
    }
    return stream;
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
    recordErase(journalString._journalTree, pos, posEnd);
}

template<typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos>
inline
void
erase(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString,
      TPos const & pos)
{
    SEQAN_CHECKPOINT;
    erase(journalString, pos, pos + 1);
}

// TODO(holtgrew): erase
// TODO(holtgrew): insertValue
// TODO(holtgrew): insert
// TODO(holtgrew): assignValue
// TODO(holtgrew): Batch-Assignment of values through segments?

// TODO(holtgrew): back
// TODO(holtgrew): begin
// TODO(holtgrew): empty
// TODO(holtgrew): end
// TODO(holtgrew): fill
// TODO(holtgrew): front
// TODO(holtgrew): getValue
// TODO(holtgrew): infix
// TODO(holtgrew): infixWithLength
// TODO(holtgrew): iter
// TODO(holtgrew): length
// TODO(holtgrew): prefix
// TODO(holtgrew): resize
// TODO(holtgrew): suffix
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
