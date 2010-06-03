#ifndef SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_H_
#define SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

// Tag: Unbalanced tree.
struct Unbalanced;


template <typename TNode, typename TTreeSpec>
class JournalTree;


template <typename _TNode>
class JournalTree<_TNode, Unbalanced>
{
  public:
    typedef _TNode TNode;
    typedef typename Position<TNode>::Type TPos;
    typedef typename Size<TNode>::Type TSize;
    typedef Unbalanced TSpec;

    // Allocator for the nodes.
    Allocator<SinglePool<sizeof(TNode)> > _nodeAllocator;
    // Length of the underlying string.
    TSize _originalStringLength;
    // The root node.
    TNode * _root;
    
    JournalTree()
            : _root(0)
    {
        SEQAN_CHECKPOINT;
    }

    JournalTree(TSize const & originalStringLength)
    {
        SEQAN_CHECKPOINT;
        reinit(*this, originalStringLength);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TNode, typename TTreeSpec>
struct Iterator<JournalTree<TNode, TTreeSpec>, Standard>
{
    typedef Iter<JournalTree<TNode, TTreeSpec>, JournalTreeIterSpec<TTreeSpec> > Type;
};

template <typename TNode, typename TTreeSpec>
struct Iterator<JournalTree<TNode, TTreeSpec> const, Standard>
{
    typedef Iter<JournalTree<TNode, TTreeSpec> const, JournalTreeIterSpec<TTreeSpec> > Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TNode>
inline
bool checkVirtualPositionsRec(TNode * const & node, unsigned & virtualPosition)
{
    SEQAN_CHECKPOINT;
    if (node == 0)
        return true;
    bool res = true;
    if (node->left)
        res = res && checkVirtualPositionsRec(node->left, virtualPosition);
    if (node->virtualPosition != virtualPosition)
        res = false;
    virtualPosition += node->length;
    if (node->right)
        res = res && checkVirtualPositionsRec(node->right, virtualPosition);
    return res;
}

template <typename TNode>
inline
bool checkVirtualPositions(TNode * const & node)
{
    SEQAN_CHECKPOINT;
    unsigned virtualPosition = 0;
    return checkVirtualPositionsRec(node, virtualPosition);
}

template <typename TNode>
inline
bool checkOrderRec(TNode * const & node)
{
    SEQAN_CHECKPOINT;
    bool result = true;
    if (node->left != 0) {
        result = result && (node->virtualPosition > node->left->virtualPosition);
        result = result && checkOrderRec(node->left);
    }
    if (node->right != 0) {
        result = result && (node->virtualPosition < node->right->virtualPosition);
        result = result && checkOrderRec(node->right);
    }
    return result;
}

template <typename TNode>
inline
bool checkOrder(TNode * const & node)
{
    SEQAN_CHECKPOINT;
    if (node == 0)
        return true;
    return checkOrderRec(node);
}

template <typename TNode>
inline
bool checkStructureRec(TNode * const & node)
{
    SEQAN_CHECKPOINT;
    bool result = true;
    if (node->left != 0) {
        result = result && (node->left->parent == node);
        result = result && checkStructureRec(node->left);
    }
    if (node->right != 0) {
        result = result && (node->right->parent == node);
        result = result && checkStructureRec(node->right);
    }
    return result;
}

template <typename TNode>
inline
bool checkStructure(TNode * const & node)
{
    if (node == 0)
        return true;
    if (node->parent == 0)
        return checkStructureRec(node);
    if (node->parent->left != 0 && node->parent->right != 0) {
        if (!((node->parent->left == node) ^ (node->parent->right == node)))
            return false;
    } else if (node->parent->left == 0 && node->parent->right != 0) {
        if (node->parent->right != node)
            return false;
    } else if (node->parent->left != 0 && node->parent->right == 0) {
        if (node->parent->left != node)
            return false;
    } else {  // both == 0
        return false;
    }
    return checkStructureRec(node);
}

template <typename TStream, typename TNode, typename TTreeSpec>
inline
TStream &
operator<<(TStream & stream, JournalTree<TNode, TTreeSpec> const & tree)
{
    SEQAN_CHECKPOINT;
    if (tree._root != 0)
        return stream << "JournalTree(" << value(tree._root) << ")";
    else
        return stream << "JournalTree()";
}

template <typename TNode>
inline
void reinit(JournalTree<TNode, Unbalanced> & tree,
            typename Size<TNode>::Type originalStringLength)
{
    SEQAN_CHECKPOINT;
    clear(tree._nodeAllocator);
    tree._originalStringLength = originalStringLength;
    TNode *tmp;
    allocate(tree._nodeAllocator, tmp, 1);
    tree._root = new (tmp) TNode(SOURCE_ORIGINAL, 0, 0, originalStringLength);
}


// Subtract delta from all nodes with virtual positions right of,
// respectively >= pos.  Note that this must not make the tree invalid.
template <typename TNode>
inline
void
_subtractFromVirtualPositionsRightOf(TNode * node,
                                     typename Position<TNode>::Type const & pos,
                                     typename Position<TNode>::Type const & delta)
{
    SEQAN_CHECKPOINT;

    if (node->virtualPosition >= pos) {
        node->virtualPosition -= delta;
        if (node->left != 0)
            _subtractFromVirtualPositionsRightOf(node->left, pos, delta);
        if (node->right != 0)
            _subtractFromVirtualPositionsRightOf(node->right, pos, delta);
    } else {  // node->virtualPosition < pos
        if (node->right != 0)
            _subtractFromVirtualPositionsRightOf(node->right, pos, delta);
    }
}


// Add delta to all nodes with virtual positions right of,
// respectively >= pos.  Note that this must not make the tree invalid.
template <typename TNode>
inline
void
_addToVirtualPositionsRightOf(TNode * node,
                              typename Position<TNode>::Type const & pos,
                              typename Position<TNode>::Type const & delta)
{
    SEQAN_CHECKPOINT;

    if (node->virtualPosition >= pos) {
        node->virtualPosition += delta;
        if (node->left != 0)
            _addToVirtualPositionsRightOf(node->left, pos, delta);
        if (node->right != 0)
            _addToVirtualPositionsRightOf(node->right, pos, delta);
    } else {  // node->virtualPosition < pos
        if (node->right != 0)
            _addToVirtualPositionsRightOf(node->right, pos, delta);
    }
}


template <typename TNode, typename TPos>
inline
void
findNodeWithVirtualPos(TNode * & node,
                       TNode * & parent,
                       JournalTree<TNode, Unbalanced> const & tree,
                       TPos const & pos)
{
    SEQAN_CHECKPOINT;

    parent = 0;
    node = tree._root;
    while (true) {
        SEQAN_ASSERT_NEQ(node, static_cast<TNode *>(0));
        if ((node->virtualPosition <= pos) && (node->virtualPosition + node->length > pos)) {
            break;
        } else if ((node->virtualPosition + node->length) <= pos) {
            if (node->right == 0)
                break;
            TNode * tmp = node;
            node = node->right;
            parent = tmp;
            continue;
        } else {  // pos < node->virtualPosition
            TNode * tmp = node;
            node = node->left;
            parent = tmp;
            continue;
        }
    }
}


template <typename TNode>
inline
void recordErase(JournalTree<TNode, Unbalanced> & tree,
                 typename Position<TNode>::Type const & pos,
                 typename Position<TNode>::Type const & posEnd)
{
    SEQAN_CHECKPOINT;
    typedef typename Position<TNode>::Type TPos;
    typedef typename Size<TNode>::Type TSize;

    SEQAN_ASSERT_TRUE(checkStructure(tree._root));
    SEQAN_ASSERT_TRUE(checkOrder(tree._root));
    SEQAN_ASSERT_TRUE(checkVirtualPositions(tree._root));

    // Handle special case of removing all of the root node.
    if (tree._root->left == 0 && tree._root->right == 0 && pos == 0 && posEnd == tree._root->length) {
        tree._root = 0;
        clear(tree._nodeAllocator);
        return;
    }
    // Handle case of an empty journal tree.
    if (tree._root == 0) {
        SEQAN_ASSERT_EQ(pos, 0u);
        SEQAN_ASSERT_EQ(posEnd, 0u);
        return;
    }

    // Find node with virtual position pos.
    TNode * node = 0;
    TNode * parent = 0;
    findNodeWithVirtualPos(node, parent, tree, pos);

    // The position to subtract values from right of.
    TPos subtractRightOf = pos;
    // Virtual begin and end position of the node.
    TPos nodeBegin = node->virtualPosition;
    TPos nodeEnd = node->virtualPosition + node->length;

    // The simple cases are: A prefix or suffix but not all of the node are erased.
    // Otherwise, we have to delete or split the existing node.
    if (nodeBegin == pos && nodeEnd == posEnd) {
//         std::cout << "whole node" << std::endl;
        // The whole node is removed.  If there is <= one child, things are
        // simple, otherwise, we replace node with its right child and perform
        // a left-right traversal to find the leaf for which node's left child
        // can become the left child of.
        TNode * left = node->left;
        TNode * right = node->right;
        if (left == 0) {
            // Replace parent's pointer.
            if (parent == 0) {
                // Node was root.
                tree._root = node->right;
                tree._root->parent = 0;
            } else {
                // Node was not root.
                if (parent->left == node) {
                    parent->left = node->right;
                } else {
                    SEQAN_ASSERT_EQ(parent->right, node);
                    parent->right = node->right;
                }
                if (node->right != 0)
                    node->right->parent = parent;
            }
        } else {
            if (right == 0) { // left != 0 && right == 0
                if (parent == 0) {
                    // Node was root.
                    node->left->parent = 0;
                    tree._root = node->left;
                } else {
                    // Node was not root.
                    if (parent->left == node) {
                        node->left->parent = parent;
                        parent->left = node->left;
                    } else {
                        node->left->parent = parent;
                        parent->right = node->left;
                    }
                }
            } else {
                // left != 0 && right != 0
                if (parent == 0) {
                    // node is root
                    node->right->parent = 0;
                    tree._root = node->right;
                } else {
                    if (parent->left == node) {
                        node->right->parent = parent;
                        parent->left = node->right;
                    } else {
                        node->right->parent = parent;
                        parent->right = node->right;
                    }
                }
                TNode * tmp = node->right->left;
                node->right->left = node->left;
                node->left->parent = node->right;
                // Left-right traversal from node on.
                TNode * current = node->left;
                SEQAN_ASSERT_NEQ(current, static_cast<TNode *>(0));
                while (current->right != 0)
                    current = current->right;
                current->right = tmp;
                if (tmp != 0)
                    tmp->parent = current;
            }
        }
        // Actually deallocate node.
        deallocate(tree._nodeAllocator, node, 1);
        // Adjust virtual positions.
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
    } else if (nodeBegin == pos && nodeEnd > posEnd) {
//         std::cout << "true prefix" << std::endl;
        // A true prefix is removed, not the whole node.
        node->length -= posEnd - pos;
        node->physicalPosition += posEnd - pos;
        subtractRightOf = posEnd;  // No need to update this node!
        // Adjust virtual positions.
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
    } else if (nodeBegin < pos && nodeEnd == posEnd) {
//         std::cout << "true suffix" << std::endl;
        // A true suffix is removed, not the whole node.
        node->length -= posEnd - pos;
        subtractRightOf = posEnd;  // No need to update this node!
        // Adjust virtual positions.
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
    } else if (nodeBegin < pos && nodeEnd > posEnd) {
//         std::cout << "true infix" << std::endl;
        // A true infix of the node is removed.  This node stores the prefix
        // and we do a right-left traversal to place the newly created node.
        TNode * tmp;
        allocate(tree._nodeAllocator, tmp, 1);
        // Perform index calculations for prefix and suffix node.
        TSize prefixLength = pos - node->virtualPosition;
        TSize deletedInfixLength = posEnd - pos;
        TSize suffixLength = node->length - prefixLength - deletedInfixLength;
        // Update prefix node.
        node->length = prefixLength;
        // Construct suffix node.
        TNode * suffixNode = new (tmp) TNode(node->segmentSource,
                                             node->physicalPosition + prefixLength + deletedInfixLength,
                                             node->virtualPosition + prefixLength,
                                             suffixLength);
        // Insert node for suffix.
        if (node->right == 0) {
            node->right = suffixNode;
            suffixNode->parent = node;
        } else {
            TNode * current = node->right;
            while (current->left != 0) {
                current = current->left;
            }
            current->left = suffixNode;
            suffixNode->parent = current;
        }
        subtractRightOf = node->virtualPosition + prefixLength + deletedInfixLength + 1;
        // Adjust virtual positions.
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
    } else {
//         std::cout << "spans more than one" << std::endl;
        // The deletion range spans more than this node.  We solve this by
        // calling this function recursively.  First, we delete a suffix of
        // this node (possibly all of it) and then we recursively remove the
        // remaining suffix of [pos, posEnd).
        TSize len = node->length - (pos - node->virtualPosition);
        recordErase(tree, pos, pos + len);
        recordErase(tree, pos, posEnd - len);
    }
//     unsigned nextId = mtRand();
//     journalTreeToDot(std::cerr, nextId, tree);
//     std::cout << tree << std::endl;
//     std::cout << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;
    SEQAN_ASSERT_TRUE(checkStructure(tree._root));
    SEQAN_ASSERT_TRUE(checkOrder(tree._root));
    SEQAN_ASSERT_TRUE(checkVirtualPositions(tree._root));
}


template <typename TNode>
inline
void recordInsertion(JournalTree<TNode, Unbalanced> & tree,
                     typename Position<TNode>::Type const & virtualPos,
                     typename Position<TNode>::Type const & physicalBeginPos,
                     typename Size<TNode>::Type const & length)
{
    SEQAN_CHECKPOINT;
    typedef typename Position<TNode>::Type TPos;
    typedef typename Size<TNode>::Type TSize;

    SEQAN_ASSERT_TRUE(checkStructure(tree._root));
    SEQAN_ASSERT_TRUE(checkOrder(tree._root));
    SEQAN_ASSERT_TRUE(checkVirtualPositions(tree._root));

    // Handle special case of empty tree.
    if (tree._root == 0) {
        SEQAN_ASSERT_EQ(virtualPos, 0u);
        if (length == 0)
            return;
        TNode * tmp;
        allocate(tree._nodeAllocator, tmp, 1);
        tree._root = new (tmp) TNode(SOURCE_PATCH, physicalBeginPos, virtualPos, length);
        return;
    }
    
    // TODO(holtgrew): Does this really, really work?
    TNode * node;
    TNode * parent;
    findNodeWithVirtualPos(node, parent, tree, virtualPos);
    SEQAN_ASSERT_LEQ(node->virtualPosition, virtualPos);

    if (node->virtualPosition + node->length > virtualPos) {
        // Found node that contains virtualPos.
        SEQAN_ASSERT_LEQ(node->virtualPosition, virtualPos);
        if (node->virtualPosition == virtualPos) {
            // Simple case: Insert left of current.
            TNode * tmp;
            allocate(tree._nodeAllocator, tmp, 1);
            TNode * insertNode = new (tmp) TNode(SOURCE_PATCH, physicalBeginPos, virtualPos, length);
            _addToVirtualPositionsRightOf(tree._root, virtualPos, length);
            insertNode->left = node->left;
            if (insertNode->left != 0)
                insertNode->left->parent = insertNode;
            node->left = insertNode;
            insertNode->parent = node;
        } else {
            // Harder case: Split current and insert new node.
            _addToVirtualPositionsRightOf(tree._root, node->virtualPosition + node->length, length);
            TPos offset = virtualPos - node->virtualPosition;
            // Create right part of the node.
            TNode * tmp;
            allocate(tree._nodeAllocator, tmp, 1);
            TNode * splitNode = new (tmp) TNode(node->segmentSource, node->physicalPosition + offset, node->virtualPosition + offset + length, node->length - offset);
            // Current node becomes left part of current node.
            node->length = offset;
            // Create insertion node.
            allocate(tree._nodeAllocator, tmp, 1);
            TNode * insertNode = new (tmp) TNode(SOURCE_PATCH, physicalBeginPos, virtualPos, length);
            // Insert into tree...
            insertNode->left = node;
            node->parent = insertNode;
            insertNode->right = splitNode;
            splitNode->parent = insertNode;
            splitNode->right = node->right;
            if (node->right != 0)
                node->right->parent = splitNode;
            node->right = 0;
            if (parent == 0) {
                // current is the root node.
                tree._root = insertNode;
                insertNode->parent = 0;
            } else {
                if (parent->left == node)
                    parent->left = insertNode;
                else
                    parent->right = insertNode;
                insertNode->parent = parent;
            }
        }
    } else {
        // Returned node with highest virtualPosition but we need to insert
        // right of it.
        SEQAN_ASSERT_EQ(node->right, static_cast<TNode *>(0));
        SEQAN_ASSERT_EQ(node->virtualPosition + node->length, virtualPos);
        TNode * tmp;
        allocate(tree._nodeAllocator, tmp, 1);
        TNode * insertNode = new (tmp) TNode(SOURCE_PATCH, physicalBeginPos, virtualPos, length);
        node->right = insertNode;
        insertNode->parent = node;
    }
//     unsigned nextId = mtRand();
//     journalTreeToDot(std::cerr, nextId, tree);
//     std::cout << tree << std::endl;
//     std::cout << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;
    SEQAN_ASSERT_TRUE(checkStructure(tree._root));
    SEQAN_ASSERT_TRUE(checkOrder(tree._root));
    SEQAN_ASSERT_TRUE(checkVirtualPositions(tree._root));
}


template <typename TNode>
inline
TNode const * 
back(JournalTree<TNode, Unbalanced> const & tree)
{
    SEQAN_CHECKPOINT;
    TNode * current = tree._root;
    while (current->right != 0)
        current = current->right;
    return current;
}


template <typename TNode>
inline
TNode const * 
front(JournalTree<TNode, Unbalanced> const & tree)
{
    SEQAN_CHECKPOINT;
    TNode * current = tree._root;
    while (current->left != 0)
        current = current->left;
    return current;
}


template <typename TStream, typename TNode>
void journalTreeToDotRec(TStream & stream, unsigned & nextId, TNode const & node)
{
    SEQAN_CHECKPOINT;
    unsigned currentId = nextId;
    nextId += 1;
    stream << "  node_" << currentId << "[label=\"source=" << node.segmentSource << ", vpos=" << node.virtualPosition << ", phpos=" << node.physicalPosition << ", len=" << node.length << "\"]" << std::endl;
    if (node.left != 0) {
        stream << "  node_" << currentId << " -> node_" << nextId << "[label=\"L\"]" << std::endl;
        journalTreeToDotRec(stream, nextId, *(node.left));
        nextId += 1;
    }
    if (node.right != 0) {
        stream << "  node_" << currentId << " -> node_" << nextId << "[label=\"R\"]" << std::endl;
        journalTreeToDotRec(stream, nextId, *(node.right));
        nextId += 1;
    }
}


template <typename TStream, typename TNode>
void journalTreeToDot(TStream & stream, unsigned & nextId, JournalTree<TNode, Unbalanced> const & journalTree)
{
    SEQAN_CHECKPOINT;
    stream << "ROOTPTR" << nextId << " -> node_" << nextId << std::endl;
    journalTreeToDotRec(stream, nextId, *journalTree._root);
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_H_
