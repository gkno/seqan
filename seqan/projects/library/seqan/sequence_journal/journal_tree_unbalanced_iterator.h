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

#ifndef SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_UNBALANCED_ITERATOR_H_
#define SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_UNBALANCED_ITERATOR_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

template <typename TJournalTreeSpec>
struct JournalTreeIterSpec;


enum IterationDirection {
    DIRECTION_NULL,
    DIRECTION_DOWN,
    DIRECTION_UP_LEFT,
    DIRECTION_UP_RIGHT
};

template <typename TJournalTree>
class Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> >
{
public:
    typedef Iter<typename _RemoveConst<TJournalTree>::Type, JournalTreeIterSpec<Unbalanced> > TNonConstIterator;
    typedef Iter<typename _RemoveConst<TJournalTree>::Type const, JournalTreeIterSpec<Unbalanced> > TConstIterator;
    typedef typename TJournalTree::TNode TNode;

    // The current node.
    TNode * _currentNode;
    // The direction we arrived from at this node.  The end node is
    // encoded as the root and _iterationDirection == DIRECTION_UP_RIGHT.
    IterationDirection _iterationDirection;

    Iter() : _currentNode(0), _iterationDirection(DIRECTION_NULL) {}

    Iter(TConstIterator const & other)
            : _currentNode(other._currentNode),
              _iterationDirection(other._iterationDirection)
    { SEQAN_CHECKPOINT; }

    // Always allow conversion from non-const.
    Iter(TNonConstIterator const & other)
            : _currentNode(other._currentNode),
              _iterationDirection(other._iterationDirection)
    { SEQAN_CHECKPOINT; }

    explicit
    Iter(TJournalTree & tree)
    {
        SEQAN_CHECKPOINT;
        _initJournalTreeIterator(*this, tree);
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

template <typename TCargo>
struct Iterator<JournalTree<TCargo, Unbalanced>, Standard>
{
    typedef Iter<JournalTree<TCargo, Unbalanced>, JournalTreeIterSpec<Unbalanced> > Type;
};

template <typename TCargo>
struct Iterator<JournalTree<TCargo, Unbalanced> const, Standard>
{
    typedef Iter<JournalTree<TCargo, Unbalanced> const, JournalTreeIterSpec<Unbalanced> > Type;
};

template <typename TJournalTree>
struct Value<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >
{
    typedef typename TJournalTree::TCargo & Type;
};

template <typename TJournalTree>
struct Value<Iter<TJournalTree const, JournalTreeIterSpec<Unbalanced> > >
{
    typedef typename TJournalTree::TCargo const & Type;
};

template <typename TJournalTree>
struct GetValue<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >
{
    typedef typename TJournalTree::TCargo Type;
};

template <typename TJournalTree>
struct GetValue<Iter<TJournalTree const, JournalTreeIterSpec<Unbalanced> > >
{
    typedef typename TJournalTree::TCargo Type;
};

template <typename TJournalTree>
struct Reference<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >
{
    typedef typename TJournalTree::TCargo & Type;
};

template <typename TJournalTree>
struct Reference<Iter<TJournalTree const, JournalTreeIterSpec<Unbalanced> > >
{
    typedef typename TJournalTree::TCargo const & Type;
};

// ============================================================================
// Functions
// ============================================================================

// For JournalTree<TNode, Unbalanced>

template <typename TNode>
inline
typename Iterator<JournalTree<TNode, Unbalanced> const, Standard>::Type
begin(JournalTree<TNode, Unbalanced> const & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return Iter<JournalTree<TNode, Unbalanced> const, JournalTreeIterSpec<Unbalanced> >(journalTree);
}

template <typename TNode>
inline
typename Iterator<JournalTree<TNode, Unbalanced>, Standard>::Type
begin(JournalTree<TNode, Unbalanced> & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return Iter<JournalTree<TNode, Unbalanced>, JournalTreeIterSpec<Unbalanced> >(journalTree);
}

template <typename TNode>
inline
typename Iterator<JournalTree<TNode, Unbalanced> const, Standard>::Type
end(JournalTree<TNode, Unbalanced> const & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<JournalTree<TNode, Unbalanced> const, Standard>::Type TIterator;
    TIterator result;
    result._currentNode = journalTree._root;
    result._iterationDirection = DIRECTION_UP_RIGHT;
    return result;
}

template <typename TNode>
inline
typename Iterator<JournalTree<TNode, Unbalanced>, Standard>::Type
end(JournalTree<TNode, Unbalanced> & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    typedef Iter<JournalTree<TNode, Unbalanced>, JournalTreeIterSpec<Unbalanced> > TIterator;
    TIterator result;
    result._currentNode = journalTree._root;
    result._iterationDirection = DIRECTION_UP_RIGHT;
    return result;
}

// For JournalTreeIterator<TJournalTree, JournalTreeIterSpec<Unbalanced> >

template <typename TJournalTree>
inline
void
_initJournalTreeIterator(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator,
                         TJournalTree & tree)
{
    SEQAN_CHECKPOINT;
    iterator._currentNode = tree._root;
    if (tree._root == 0) {
        iterator._iterationDirection = DIRECTION_UP_RIGHT;
    } else {
        iterator._iterationDirection = DIRECTION_DOWN;
        while (goLeft(iterator))
            continue;  // Left-only traversal.
    }
}

// Initialize journal tree iterator to end, i.e. root & direction is up right.
template <typename TJournalTree>
inline
void
_initJournalTreeIteratorEnd(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator,
                            TJournalTree & tree)
{
    SEQAN_CHECKPOINT;
    iterator._currentNode = tree._root;
    iterator._iterationDirection = DIRECTION_UP_RIGHT;
}

// TODO(holtgrew): Unused, remove?
/*
template <typename TJournalTree>
inline
void
setValue(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator,
         typename Value<TJournalTree>::Type const & value)
{
    SEQAN_XXXCHECKPOINT;
    setValue(*iterator._currentNode, value);
}
*/

template <typename TJournalTree>
inline
typename Value<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >::Type
value(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
{
    SEQAN_CHECKPOINT;
    return cargo(*iterator._currentNode);
}

template <typename TJournalTree>
inline
typename Value<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >::Type
value(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & iterator)
{
    SEQAN_CHECKPOINT;
    return cargo(*iterator._currentNode);
}

// TODO(holtgrew): Unused, remove?
/*
template <typename TJournalTree>
inline
typename Value<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >::Type
operator*(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
{
    SEQAN_XXXCHECKPOINT;
    return value(iterator);
}
*/

// TODO(holtgrew): Unused, remove?
/*
template <typename TJournalTree>
inline
typename Value<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >::Type
operator*(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & iterator)
{
    SEQAN_XXXCHECKPOINT;
    return value(iterator);
}
*/

template <typename TJournalTree>
inline
bool
hasLeftChild(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
{
    SEQAN_CHECKPOINT;
    return iterator._currentNode->left != 0;
}

template <typename TJournalTree>
inline
bool
goLeft(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (!hasLeftChild(iterator))
        return false;
    iterator._iterationDirection = DIRECTION_DOWN;
    iterator._currentNode = iterator._currentNode->left;
    return true;
}

template <typename TJournalTree>
inline
bool
hasRightChild(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
{
    SEQAN_CHECKPOINT;
    return iterator._currentNode->right != 0;
}

template <typename TJournalTree>
inline
bool
goRight(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
{
    SEQAN_CHECKPOINT;
    if (!hasRightChild(iterator))
        return false;
    iterator._iterationDirection = DIRECTION_DOWN;
    iterator._currentNode = iterator._currentNode->right;
    return true;
}

template <typename TJournalTree>
inline
bool
hasParent(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & iterator)
{
    SEQAN_CHECKPOINT;
    return iterator._currentNode->parent != 0;
}
    
template <typename TJournalTree>
inline
bool
goUp(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
{
    if (!hasParent(iterator)) {
        // Up from the root means we go to end, set direction to up-right.
        iterator._iterationDirection = DIRECTION_UP_RIGHT;
        return false;
    }
    if (iterator._currentNode->parent->left == iterator._currentNode)
        iterator._iterationDirection = DIRECTION_UP_LEFT;
    else
        iterator._iterationDirection = DIRECTION_UP_RIGHT;
    iterator._currentNode = iterator._currentNode->parent;
    return true;
}

// template <typename TJournalTree>
// inline
// bool
// atEnd(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & iterator)
// {
//     SEQAN_CHECKPOINT;
//     return (iterator._currentNode == 0) || (!hasParent(iterator) && (iterator._iterationDirection == DIRECTION_UP_RIGHT));
// }

// template <typename TJournalTree>
// inline
// bool
// atEnd(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
// {
//     SEQAN_CHECKPOINT;
//     typedef Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > TIterator;
//     return (iterator._currentNode == 0) || (!hasParent(iterator) && (iterator._iterationDirection == DIRECTION_UP_RIGHT));
// }

template <typename TJournalTree>
inline
Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > &
operator++(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator)
{
    SEQAN_CHECKPOINT;
    switch (iterator._iterationDirection) {
        case DIRECTION_DOWN:
            // Arrived here by going down (either right or left).  Next is
            // either right-left traversal (if node has a child to the right)
            // or going up if we do not have a child on the right.
            if (goRight(iterator)) {
                while (goLeft(iterator))
                    continue;
            } else {
                while (goUp(iterator) && iterator._iterationDirection == DIRECTION_UP_RIGHT)
                    continue;
            }
            break;
        case DIRECTION_UP_LEFT:
            // We came up from our left child.  Next is either right-left
            // traversal (if node has a child to the right) or going up if
            // we do not have a child on the left.
            if (goRight(iterator)) {
                while (goLeft(iterator))
                    continue;
            } else {
                while (goUp(iterator) && iterator._iterationDirection == DIRECTION_UP_RIGHT)
                    continue;
            }
            break;
        case DIRECTION_UP_RIGHT:
            // We came here from our right child.  Next is up until we got
            // to the current node from it's left child.
        default:
            SEQAN_ASSERT_FAIL("Invalid iteration direction.");
    }
        
    return iterator;    
}

// TODO(holtgrew): Unused, remove?
/*
template <typename TJournalTree>
inline
Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> >
operator++(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > & iterator,
           int postfix)
{
    SEQAN_XXXCHECKPOINT;
    Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > temp(iterator);
    ++iterator;
    return temp;
}
*/

template <typename TJournalTree>
inline
bool
operator==(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & a,
           Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & b)
{
    SEQAN_CHECKPOINT;
    return (a._currentNode == b._currentNode) && (a._iterationDirection == b._iterationDirection);
}

template <typename TJournalTree>
inline
bool
operator==(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & a,
           typename IterComplementConst<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return (a._currentNode == b._currentNode) && (a._iterationDirection == b._iterationDirection);
}

template <typename TJournalTree>
inline
bool
operator!=(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & a,
           Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & b)
{
    SEQAN_CHECKPOINT;
    return (a._currentNode != b._currentNode) || (a._iterationDirection != b._iterationDirection);
}

template <typename TJournalTree>
inline
bool
operator!=(Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > const & a,
           typename IterComplementConst<Iter<TJournalTree, JournalTreeIterSpec<Unbalanced> > >::Type const & b)
{
    SEQAN_CHECKPOINT;
    return (a._currentNode != b._currentNode) || (a._iterationDirection != b._iterationDirection);
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_JOURNAL_UNBALANCED_TREE_ITERATOR_H_
