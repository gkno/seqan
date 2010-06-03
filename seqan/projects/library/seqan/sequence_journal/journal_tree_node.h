#ifndef SEQAN_SEQUENCE_JOURNAL_SEGMENT_NODE_H_
#define SEQAN_SEQUENCE_JOURNAL_SEGMENT_NODE_H_

namespace seqan {

// ============================================================================
// Enums, Classes
// ============================================================================

enum SegmentSource {
    SOURCE_ORIGINAL,
    SOURCE_PATCH
};


template <typename TPos, typename TSize>
struct SegmentNode
{
    // Left child.
    SegmentNode *left;
    // Right child.
    SegmentNode *right;
    // Parent, 0 for root.
    SegmentNode *parent;
    // Flag for where the segment comes from.
    SegmentSource segmentSource;
    // Position in the original string or the insertion buffer,
    // depending on segmentSource.
    TPos virtualPosition;
    // Position in the virtual string.
    TPos physicalPosition;
    // Length of the segment.
    TSize length;
    
    SegmentNode() : left(0), right(0) {}

    SegmentNode(SegmentSource const & _segmentSource,
                TPos const & _physicalPosition,
                TPos const & _virtualPosition,
                TPos const & _length)
            : left(0),
              right(0),
              parent(0),
              segmentSource(_segmentSource),
              virtualPosition(_virtualPosition),
              physicalPosition(_physicalPosition),
              length(_length) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TPos, typename TSize>
struct Position<SegmentNode<TPos, TSize> >
{
    typedef TPos Type;
};

template <typename TPos, typename TSize>
struct Position<SegmentNode<TPos, TSize> const>
        : Position<SegmentNode<TPos, TSize> > {};


template <typename TPos, typename TSize>
struct Size<SegmentNode<TPos, TSize> >
{
    typedef TSize Type;
};

template <typename TPos, typename TSize>
struct Size<SegmentNode<TPos, TSize> const>
        : Size<SegmentNode<TPos, TSize> > {};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TPos, typename TSize>
TStream & operator<<(TStream & stream, SegmentNode<TPos, TSize> const & node)
{
    SEQAN_CHECKPOINT;
    stream << "SegmentNode(add=" << &node
           << ", segmentSource=" << node.segmentSource
           << ", virtualPosition=" << node.virtualPosition
           << ", physicalPosition=" << node.physicalPosition
           << ", length=" << node.length
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

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_SEGMENT_NODE_H_
