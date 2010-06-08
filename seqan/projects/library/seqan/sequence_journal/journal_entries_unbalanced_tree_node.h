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
  Code for the node of JournalEntries<UnorderedTree>, i.e. a node of an
  unbalanced binary search tree.  This could easily be extended to a Red Black
  Tree with a flag (which could actually be stored in the lowest bits of one
  of the pointers).  However, since std::set already implements balanced
  trees, it should be sufficient to use a STL implementation of std::set for
  balanced binary search trees.
  ==========================================================================*/

#ifndef SEQAN_SEQUENCE_JOURNAL_JOURNAL_ENTRIES_UNBALANCED_TREE_NODE_H_
#define SEQAN_SEQUENCE_JOURNAL_JOURNAL_ENTRIES_UNBALANCED_TREE_NODE_H_

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

template <typename TCargo>
struct JournalEntriesUnorderedTreeNode
{
    // Left child.
    JournalEntriesUnorderedTreeNode * left;
    // Right child.
    JournalEntriesUnorderedTreeNode * right;
    // Parent, 0 for root.
    JournalEntriesUnorderedTreeNode * parent;
    // The actual payload:  The journal entry.
    TCargo cargo;
    
    JournalEntriesUnorderedTreeNode() : left(0), right(0), parent(0) {}

    JournalEntriesUnorderedTreeNode(TCargo const & _cargo)
            : left(0),
              right(0),
              parent(0),
              cargo(_cargo) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// TODO(holtgrew): Rename value to cargo?

template <typename T>
struct Cargo;

template <typename TCargo>
struct Cargo<JournalEntriesUnorderedTreeNode<TCargo> >
{
    typedef TCargo Type;
};

template <typename TCargo>
struct Cargo<JournalEntriesUnorderedTreeNode<TCargo> const>
{
    typedef TCargo const Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TCargo>
TStream & operator<<(TStream & stream, JournalEntriesUnorderedTreeNode<TCargo> const & node)
{
    stream << "JournalEntriesUnorderedTreeNode(add=" << &node
           << ", cargo=" << node.cargo
           << ", parent=" << node.parent
           << ", left=";
    if (node.left)
        stream << *node.left;
    else
        stream << "NULL";
    stream << ", right=";
    if (node.right)
        stream << *node.right;
    else
        stream << "NULL";
    return stream << ")";
}

template <typename TCargo>
typename Cargo<JournalEntriesUnorderedTreeNode<TCargo> >::Type &
cargo(JournalEntriesUnorderedTreeNode<TCargo> & node)
{
    SEQAN_CHECKPOINT;
    return node.cargo;
}

// TODO(holtgrew): Unused, remove?
/*
template <typename TCargo>
typename Cargo<JournalEntriesUnorderedTreeNode<TCargo> const>::Type &
cargo(JournalEntriesUnorderedTreeNode<TCargo> const & node)
{
    SEQAN_XXXCHECKPOINT;
    return node.cargo;
}

template <typename TCargo>
typename Cargo<JournalEntriesUnorderedTreeNode<TCargo> const>::Type
getCargo(JournalEntriesUnorderedTreeNode<TCargo> const & node)
{
    SEQAN_XXXCHECKPOINT;
    return node.cargo;
}
*/

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_JOURNAL_ENTRIES_UNBALANCED_TREE_NODE_H_
