 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_JOURNAL_GENERATED_FORWARDS_H 
#define SEQAN_HEADER_SEQUENCE_JOURNAL_GENERATED_FORWARDS_H 

//////////////////////////////////////////////////////////////////////////////
// NOTE: This file is automatically generated by build_forwards.py
//       Do not edit this file manually!
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// Journal

template <typename TStringSpec, typename TJournalSpec> struct Journal;       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(34)

//____________________________________________________________________________
// JournalTree

template <typename TNode, typename TTreeSpec> class JournalTree;       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(15)

//____________________________________________________________________________
// SegmentNode

template <typename TPos, typename TSize> struct SegmentNode;       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree_node.h"(18)

//____________________________________________________________________________
// Unbalanced

struct Unbalanced;       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(11)

} //namespace seqan


//////////////////////////////////////////////////////////////////////////////
// TYPEDEFS


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// _addToVirtualPositionsRightOf

template <typename TNode> inline void _addToVirtualPositionsRightOf(TNode * node, typename Position<TNode>::Type const & pos, typename Position<TNode>::Type const & delta);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(109)

//____________________________________________________________________________
// _subtractFromVirtualPositionsRightOf

template <typename TNode> inline void _subtractFromVirtualPositionsRightOf(TNode * node, typename Position<TNode>::Type const & pos, typename Position<TNode>::Type const & delta);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(85)

//____________________________________________________________________________
// assignInfix

template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos, typename TString> inline void assignInfix(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString, TPos const & beginPos, TPos const & endPos, TString const & valueString);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(271)

//____________________________________________________________________________
// assignValue

template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos> inline void assignValue(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString, TPos const & pos, TValue const & value);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(284)

//____________________________________________________________________________
// back

template <typename TNode> inline TNode const * back(JournalTree<TNode, Unbalanced> const & tree);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(343)

//____________________________________________________________________________
// clear

template <typename TValue, typename TStringSpec, typename TJournalSpec> inline void clear(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(196)

//____________________________________________________________________________
// erase

template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos> inline void erase(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString, TPos const & pos, TPos const & posEnd);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(216)
template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos> inline void erase(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString, TPos const & pos);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(228)

//____________________________________________________________________________
// findNodeWithVirtualPos

template <typename TNode, typename TPos> inline void findNodeWithVirtualPos(TNode * & node, TNode * & parent, JournalTree<TNode, Unbalanced> const & tree, TPos const & pos);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(132)

//____________________________________________________________________________
// front

template <typename TNode> inline TNode const * front(JournalTree<TNode, Unbalanced> const & tree);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(355)
template <typename TValue, typename TStringSpec, typename TJournalSpec> inline TValue const & front(String<TValue, Journal<TStringSpec, TJournalSpec> > const & journalString);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(306)

//____________________________________________________________________________
// host

template <typename TValue, typename TStringSpec, typename TJournalSpec> inline String<TValue, TStringSpec> & host(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(182)

//____________________________________________________________________________
// insert

template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos, typename TString> inline void insert(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString, TPos const & pos, TString const & seq);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(241)

//____________________________________________________________________________
// insertValue

template <typename TValue, typename TStringSpec, typename TJournalSpec, typename TPos> inline void insertValue(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString, TPos const & pos, TValue const & value);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(256)

//____________________________________________________________________________
// length

template <typename TValue, typename TStringSpec, typename TJournalSpec> inline typename Size<String<TValue, Journal<TStringSpec, TJournalSpec> > >::Type length(String<TValue, Journal<TStringSpec, TJournalSpec> > const & journalString);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(340)

//____________________________________________________________________________
// operator<<

template <typename TStream, typename TNode, typename TTreeSpec> inline TStream & operator<<(TStream & stream, JournalTree<TNode, TTreeSpec> const & tree);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(58)
template <typename TStream, typename TPos, typename TSize> TStream & operator<<(TStream & stream, SegmentNode<TPos, TSize> const & node);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree_node.h"(78)
template <typename TStream, typename TValue, typename TStringSpec, typename TJournalSpec> inline TStream & operator<<(TStream & stream, String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(126)

//____________________________________________________________________________
// recordErase

template <typename TNode> inline void recordErase(JournalTree<TNode, Unbalanced> & tree, typename Position<TNode>::Type const & pos, typename Position<TNode>::Type const & posEnd);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(161)

//____________________________________________________________________________
// recordInsertion

template <typename TNode> inline void recordInsertion(JournalTree<TNode, Unbalanced> & tree, typename Position<TNode>::Type const & virtualPos, typename Position<TNode>::Type const & physicalBeginPos, typename Size<TNode>::Type const & length);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(280)

//____________________________________________________________________________
// reinit

template <typename TNode> inline void reinit(JournalTree<TNode, Unbalanced> & tree, typename Size<TNode>::Type originalStringLength);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/journal_tree.h"(67)

//____________________________________________________________________________
// setHost

template <typename TValue, typename TStringSpec, typename TJournalSpec> inline void setHost(String<TValue, Journal<TStringSpec, TJournalSpec> > & journalString, String<TValue, TStringSpec> & str);       	// "/Users/manuel/Development/seqan-trunk/projects/library/seqan/sequence_journal/string_journal.h"(167)

} //namespace seqan

#endif

