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

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TNode, typename TTreeSpec>
inline
TStream &
operator<<(TStream & stream, JournalTree<TNode, TTreeSpec> const & tree)
{
    SEQAN_CHECKPOINT;
    return stream << "JournalTree(" << value(tree._root) << ")";
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
    std::cout << "_subtractFromVirtualPositionsRightOf(" << *node << ", pos=" << pos << ", delta=" << delta << ")" << std::endl;
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


template <typename TNode, typename TPos>
inline
void
findNodeWithVirtualPos(TNode * & node, TNode * & parent, JournalTree<TNode, Unbalanced> const & tree, TPos const & pos)
{
    parent = 0;
    node = tree._root;
    while (true) {
        SEQAN_ASSERT_NEQ(node, static_cast<TNode *>(0));
        if ((node->virtualPosition <= pos) && (node->virtualPosition + node->length > pos)) {
            break;
        } else if ((node->virtualPosition + node->length) <= pos) {
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
    // TODO(holtgrew): Guard against removing the root node.

    // Find node with virtual position pos.
    TNode * node = 0;
    TNode * parent = 0;
    findNodeWithVirtualPos(node, parent, tree, pos);

    // The position to subtract values from right of.
    TPos subtractRightOf = pos;
    // Virtual begin and end position of the node.
    TPos nodeBegin = node->virtualPosition;
    TPos nodeEnd = node->virtualPosition + node->length;
    std::cout << "node begin = " << nodeBegin << ", nodeEnd = " << nodeEnd << std::endl;

    // The simple cases are: A prefix or suffix but not all of the node are erased.
    // Otherwise, we have to delete or split the existing node.
    if (nodeBegin == pos && nodeEnd == posEnd) {
        std::cout << "WHOLE NODE" << std::endl;
        // The whole node is removed.  If there is <= one child, things are
        // simple, otherwise, we replace node with its right child and perform
        // a left-right traversal to find the leaf for which node's left child
        // can become the left child of.
        TNode * left = node->left;
        TNode * right = node->right;
        if (left == 0) {
            // Replace parent's pointer.
            if (parent->left == node)
                parent->left = node->right;
            else
                parent->right = node->right;
        } else {
            if (right == 0) {
                // left != 0 && right == 0
                if (parent->left == node)
                    parent->left = node->left;
                else
                    parent->right = node->left;
            } else {
                // left != 0 && right != 0
                if (parent->left == node) {
                    parent->left = node->right;
                } else {
                    parent->right = node->right;
                }
                TNode * tmp = node->right->left;
                node->right->left = node->left;
                // Left-right traversal from node on.
                TNode * current = node->left;
                SEQAN_ASSERT_NEQ(current, static_cast<TNode *>(0));
                while (current->right != 0)
                    current = current->right;
                current->right = tmp;
            }
        }
        // Actually deallocate node.
        deallocate(tree._nodeAllocator, node, 1);
        // Adjust virtual positions.
        std::cout << "before subtraction <<<" << tree << ">>>" << std::endl;
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
        std::cout << "after subtraction <<<" << tree << ">>>" << std::endl;
    } else if (nodeBegin == pos && nodeEnd > posEnd) {
        // A true prefix is removed, not the whole node.
        std::cout << "TRUE PREFIX" << std::endl;
        node->length -= posEnd - pos;
        node->physicalPosition += posEnd - pos;
        subtractRightOf = posEnd;  // No need to update this node!
        // Adjust virtual positions.
        std::cout << "before subtraction <<<" << tree << ">>>" << std::endl;
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
        std::cout << "after subtraction <<<" << tree << ">>>" << std::endl;
    } else if (nodeBegin < pos && nodeEnd == posEnd) {
        // A true suffix is removed, not the whole node.
        std::cout << "TRUE SUFFIX" << std::endl;
        node->length -= posEnd - pos;
        subtractRightOf = posEnd;  // No need to update this node!
        // Adjust virtual positions.
        std::cout << "before subtraction <<<" << tree << ">>>" << std::endl;
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
        std::cout << "after subtraction <<<" << tree << ">>>" << std::endl;
    } else if (nodeBegin < pos && nodeEnd > posEnd) {
        std::cout << "TRUE INFIX" << std::endl;
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
        } else {
            TNode * current = node->right;
            while (current->left != 0)
                current = current->left;
            current->left = suffixNode;
        }
        subtractRightOf = node->physicalPosition + prefixLength + deletedInfixLength + 1;
        // Adjust virtual positions.
        std::cout << "before subtraction <<<" << tree << ">>>" << std::endl;
        _subtractFromVirtualPositionsRightOf(tree._root, subtractRightOf, posEnd - pos);
        std::cout << "after subtraction <<<" << tree << ">>>" << std::endl;
    } else {
        std::cout << "RANGE SPANS" << std::endl;
        // The deletion range spans more than this node.  We solve this by
        // calling this function recursively.  First, we delete a suffix of
        // this node (possibly all of it) and then we recursively remove the
        // remaining suffix of [pos, posEnd).
        TSize len = node->length - (pos - node->virtualPosition);
        std::cout << "recordErase(tree, " << pos << ", " << pos + len << ")" << std::endl;
        recordErase(tree, pos, pos + len);
        std::cout << "recordErase(tree, " << pos << ", " << posEnd - len << ")" << std::endl;
        recordErase(tree, pos, posEnd - len);
    }
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_H_
