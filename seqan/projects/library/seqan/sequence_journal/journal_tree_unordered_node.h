#ifndef SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_UNORDERED_NODE_H_
#define SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_UNORDERED_NODE_H_

namespace seqan {

// ============================================================================
// Enums, Classes
// ============================================================================

// TODO(holtgrew): s/SegmentNode/JournalTreeUnorderedNode/
template <typename TCargo>
struct SegmentNode
{
    // Left child.
    SegmentNode * left;
    // Right child.
    SegmentNode * right;
    // Parent, 0 for root.
    SegmentNode * parent;
    // The actual payload:  The journal entry.
    TCargo cargo;
    
    SegmentNode() : left(0), right(0), parent(0) {}

    SegmentNode(TCargo const & _cargo)
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
struct Cargo<SegmentNode<TCargo> >
{
    typedef TCargo Type;
};

template <typename TCargo>
struct Cargo<SegmentNode<TCargo> const>
{
    typedef TCargo const Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TCargo>
TStream & operator<<(TStream & stream, SegmentNode<TCargo> const & node)
{
    SEQAN_CHECKPOINT;
    stream << "SegmentNode(add=" << &node
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
typename Cargo<SegmentNode<TCargo> >::Type &
cargo(SegmentNode<TCargo> & node)
{
    SEQAN_CHECKPOINT;
    return node.cargo;
}

template <typename TCargo>
typename Cargo<SegmentNode<TCargo> const>::Type &
cargo(SegmentNode<TCargo> const & node)
{
    SEQAN_CHECKPOINT;
    return node.cargo;
}

template <typename TCargo>
typename Cargo<SegmentNode<TCargo> const>::Type
getCargo(SegmentNode<TCargo> const & node)
{
    SEQAN_CHECKPOINT;
    return node.cargo;
}

template <typename TCargo>
SegmentNode<TCargo> *
left(SegmentNode<TCargo> & node)
{
    SEQAN_CHECKPOINT;
    return node.left;
}

template <typename TCargo>
SegmentNode<TCargo> const *
left(SegmentNode<TCargo> const & node)
{
    SEQAN_CHECKPOINT;
    return node.left;
}

template <typename TCargo>
SegmentNode<TCargo> *
right(SegmentNode<TCargo> & node)
{
    SEQAN_CHECKPOINT;
    return node.right;
}

template <typename TCargo>
SegmentNode<TCargo> const *
right(SegmentNode<TCargo> const & node)
{
    SEQAN_CHECKPOINT;
    return node.right;
}

template <typename TCargo>
SegmentNode<TCargo> *
parent(SegmentNode<TCargo> & node)
{
    SEQAN_CHECKPOINT;
    return node.parent;
}

template <typename TCargo>
SegmentNode<TCargo> const *
parent(SegmentNode<TCargo> const & node)
{
    SEQAN_CHECKPOINT;
    return node.parent;
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_UNORDERED_NODE_H_
