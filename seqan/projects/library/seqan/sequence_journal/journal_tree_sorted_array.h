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
  Journal implementation using a sorted string of journal elements.
  ==========================================================================*/

#ifndef SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_SORTED_ARRAY_H_
#define SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_SORTED_ARRAY_H_

namespace seqan {

// ============================================================================
// Tags, Classes
// ============================================================================

// Tag: SortedArray.
struct SortedArray {};


template <typename TNode, typename TTreeSpec>
class JournalTree;


template <typename _TCargo>
class JournalTree<_TCargo, SortedArray>
{
public:
    typedef _TCargo TCargo;
    typedef typename Size<TCargo>::Type TSize;

    String<TCargo> _journalNodes;
    TSize _originalStringLength;

    JournalTree() : _originalStringLength(0) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TCargo>
struct Iterator<JournalTree<TCargo, SortedArray>, Standard>
{
    typedef typename Iterator<String<TCargo>, Standard>::Type Type;
};


template <typename TCargo>
struct Iterator<JournalTree<TCargo, SortedArray> const, Standard>
{
    typedef typename Iterator<String<TCargo> const, Standard>::Type Type;
};

template <typename TCargo>
struct Reference<JournalTree<TCargo, SortedArray> >
{
    typedef TCargo & Type;
};

template <typename TCargo>
struct Reference<JournalTree<TCargo, SortedArray> const>
{
    typedef TCargo const & Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TNode>
inline
TStream &
operator<<(TStream & stream, JournalTree<TNode, SortedArray> const & tree)
{
    stream << "JournalTree(";
    for (unsigned i = 0; i < length(tree._journalNodes); ++i) {
        if (i > 0) stream << ", ";
        stream << tree._journalNodes[i];
    }
    stream << ")";
    return stream;
}

template <typename TCargo>
bool _checkSortedArrayTree(JournalTree<TCargo, SortedArray> const & tree)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TCargo> const, Standard>::Type TIterator;
    if (length(tree._journalNodes) == 0)
        return true;
    if (tree._journalNodes[0].virtualPosition != 0)
        return false;
    if (tree._journalNodes[0].length == 0)
        return false;
    for (TIterator it = begin(tree._journalNodes, Standard()) + 1, itend = end(tree._journalNodes, Standard()); it != itend; ++it) {
		if (it->length == 0)
			return false;
        if ((it - 1)->virtualPosition >= it->virtualPosition)
            return false;
        if ((it - 1)->virtualPosition + (it - 1)->length != it->virtualPosition)
            return false;
    }
	return true;
}

template <typename TCargo>
typename Iterator<JournalTree<TCargo, SortedArray>, Standard>::Type
begin(JournalTree<TCargo, SortedArray> & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return begin(journalTree._journalNodes, Standard());
}


template <typename TCargo>
typename Iterator<JournalTree<TCargo, SortedArray> const, Standard>::Type
begin(JournalTree<TCargo, SortedArray> const & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return begin(journalTree._journalNodes, Standard());
}


template <typename TCargo>
typename Iterator<JournalTree<TCargo, SortedArray>, Standard>::Type
end(JournalTree<TCargo, SortedArray> & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return end(journalTree._journalNodes, Standard());
}


template <typename TCargo>
typename Iterator<JournalTree<TCargo, SortedArray> const, Standard>::Type
end(JournalTree<TCargo, SortedArray> const & journalTree, Standard const &)
{
    SEQAN_CHECKPOINT;
    return end(journalTree._journalNodes, Standard());
}


template <typename TCargo>
inline
void reinit(JournalTree<TCargo, SortedArray> & tree,
            typename Size<TCargo>::Type originalStringLength)
{
    SEQAN_CHECKPOINT;
    clear(tree._journalNodes);
    appendValue(tree._journalNodes, TCargo(TCargo(SOURCE_ORIGINAL, 0, 0, originalStringLength)));
    tree._originalStringLength = originalStringLength;
}


template <typename TCargo>
inline
void recordInsertion(JournalTree<TCargo, SortedArray> & tree,
                     typename Position<TCargo>::Type const & virtualPosition,
                     typename Position<TCargo>::Type const & physicalBeginPos,
                     typename Size<TCargo>::Type const & len)
{
    SEQAN_CHECKPOINT;
    typedef typename Size<TCargo>::Type TSize;
    typedef typename Position<TCargo>::Type TPos;
    typedef typename Iterator<String<TCargo>, Standard>::Type TIterator;
    typedef JournalEntryLtByVirtualPos<TPos, TSize> TCmp;

	//std::cerr << __FILE__ << ":" << __LINE__ << " -- INSERT(" << virtualPosition << ", " << physicalBeginPos << ", " << len << ")" << std::endl;
    //std::cerr << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    // Handle special case that the entry list is empty.
    if (empty(tree._journalNodes)) {
        SEQAN_ASSERT_EQ(virtualPosition, 0u);
        if (len == 0)
            return;
        appendValue(tree._journalNodes, TCargo(SOURCE_PATCH, physicalBeginPos, virtualPosition, len));
        return;
    }

    // Find position in sorted array of nodes to insert in.
    TCargo refCargo;
    refCargo.virtualPosition = virtualPosition;
    TIterator iter = std::upper_bound(begin(tree._journalNodes, Standard()),
                                      end(tree._journalNodes, Standard()),
                                      refCargo,
                                      TCmp());
    // TODO(holtgrew): Maybe move and update entries right of pos at the same time?
    
	// MUST NOT find begin.
	SEQAN_ASSERT_TRUE(iter != begin(tree._journalNodes, Standard()));
	--iter;
	
    // Create new journal entries.
    String<TCargo> buffer;
    reserve(buffer, 3, Exact());
	if (iter->virtualPosition + iter->length > virtualPosition) {
		TPos pos = iter - begin(tree._journalNodes, Standard());
		TPos shiftRightOf = pos;
        // Found node that contains virtualPos.
        SEQAN_ASSERT_LEQ(iter->virtualPosition, virtualPosition);
        if (iter->virtualPosition == virtualPosition) {
            // Simple case:  Insert left of iter.
            insertValue(tree._journalNodes, pos, TCargo(SOURCE_PATCH, physicalBeginPos, virtualPosition, len));
			shiftRightOf += 1;
        } else {
            // Harder case:  Split current and insert new node.
            TPos offset = virtualPosition - iter->virtualPosition;
            appendValue(buffer, TCargo(iter->segmentSource, iter->physicalPosition, iter->virtualPosition, offset));
            appendValue(buffer, TCargo(SOURCE_PATCH, physicalBeginPos, virtualPosition, len));
            appendValue(buffer, TCargo(iter->segmentSource, iter->physicalPosition + offset, virtualPosition + len, iter->length - offset));
            // Insert new journal entries.
            infix(tree._journalNodes, pos, pos + 1) = buffer;
			shiftRightOf += 3;
        }
        // Update journal entries right of pos.
        for (TIterator it = begin(tree._journalNodes, Standard()) + shiftRightOf, itend = end(tree._journalNodes, Standard()); it != itend; ++it)
            it->virtualPosition += len;
    } else {
        // Insert at end.
		SEQAN_ASSERT_EQ(virtualPosition, iter->virtualPosition + iter->length);
        appendValue(tree._journalNodes, TCargo(SOURCE_PATCH, physicalBeginPos, virtualPosition, len));
    }
    //std::cerr << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    SEQAN_ASSERT_TRUE(_checkSortedArrayTree(tree));
}

template <typename TCargo>
inline
void recordErase(JournalTree<TCargo, SortedArray> & tree,
                 typename Position<TCargo>::Type const & pos,
                 typename Position<TCargo>::Type const & posEnd)
{
    SEQAN_CHECKPOINT;
    typedef typename Size<TCargo>::Type TSize;
    typedef typename Position<TCargo>::Type TPos;
    typedef typename Iterator<String<TCargo>, Standard>::Type TIter;
    typedef JournalEntryLtByVirtualPos<TPos, TSize> TCmp;
	//std::cerr << __FILE__ << ":" << __LINE__ << " -- ERASE(" << pos << ", " << posEnd << ")" << std::endl;
    //std::cerr << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    // Handle special case of removing all of the singleton existing entry.
    if (length(tree._journalNodes) == 1 && pos == 0 && front(tree._journalNodes).length == posEnd) {
        clear(tree._journalNodes);
        return;
    }
    // Handle case of an empty journal.
    if (length(tree._journalNodes) == 0) {
        SEQAN_ASSERT_EQ(pos, 0u);
        SEQAN_ASSERT_EQ(posEnd, 0u);
        return;
    }

    // Find node.
    TCargo refCargo;
    refCargo.virtualPosition = pos;
    TIter it = std::upper_bound(
            begin(tree._journalNodes, Standard()),
            end(tree._journalNodes, Standard()),
            refCargo,
            TCmp());

    // We will shift the virtual positions of all entries right of and
    // including beginShiftPos by delta positions to the left.
    TPos delta = 0;
    TPos beginShiftPos = 0;

	// MUST NOT find begin.
	SEQAN_ASSERT_TRUE(it != begin(tree._journalNodes, Standard()));
	--it;
	
	TPos itPos = it - begin(tree._journalNodes, Standard());
	if (it->virtualPosition == pos && it->length == posEnd - pos) {
		// Remove the whole entry.
		erase(tree._journalNodes, itPos);
		delta = posEnd - pos;
		beginShiftPos = itPos;
	} else if (it->virtualPosition == pos && it->length > posEnd - pos) {
		// Remove a prefix of the entry.
		SEQAN_ASSERT_LT(pos, it->virtualPosition + it->length);
		delta = posEnd - pos;
		it->physicalPosition += delta;
		it->length -= delta;
		beginShiftPos = itPos + 1;
	} else if (it->virtualPosition < pos && it->virtualPosition + it->length == posEnd) {
		// Remove a suffix of the entry.
		SEQAN_ASSERT_GT(pos, it->virtualPosition);
		delta = posEnd - pos;
		it->length -= delta;
		beginShiftPos = itPos + 1;
	} else if (it->virtualPosition < pos && it->virtualPosition + it->length > posEnd) {
		// Remove a true infix of the entry.
		TSize prefixLength = pos - it->virtualPosition;
		TSize suffixLength = it->length - prefixLength - (posEnd - pos);
		TSize removedInfixLength = posEnd - pos;
		// Insert a new entry for the right part.
		TCargo tmpEntry(it->segmentSource, it->physicalPosition + prefixLength + removedInfixLength, it->virtualPosition + prefixLength, suffixLength);
		insertValue(tree._journalNodes, itPos + 1, tmpEntry);
		// Update the left part.
		it->length -= removedInfixLength + suffixLength;
		// Set shift position and delta.
		delta = removedInfixLength;
		beginShiftPos = itPos + 2;
	} else {
		// Remove more than one entry.
		TPos rmBeginPos = itPos;
		TPos rmEndPos = itPos;
		if (it->virtualPosition != pos) {
			// Do not remove all of first.
			delta += it->length - (pos - it->virtualPosition);
			rmBeginPos += 1;
			rmEndPos += 1;
			it->length = (pos - it->virtualPosition);
		} else {
			// Remove all of first.
			delta = it->length;
			rmEndPos += 1;
		}
		it += 1;
		while (posEnd > it->virtualPosition + it->length) {
			rmEndPos += 1;
			delta += it->length;
			it += 1;
		}
		if (it->virtualPosition + it->length == posEnd) {
			// Remove all of last.
			rmEndPos += 1;
			delta += it->length;
			beginShiftPos = rmBeginPos;
		} else {
			// Do not remove all of last.
			SEQAN_ASSERT_GT(it->virtualPosition + it->length, posEnd);
			TSize tmpDelta = delta;
			delta += posEnd - it->virtualPosition;
			it->physicalPosition += posEnd - it->virtualPosition;
			it->length -= posEnd - it->virtualPosition;
			// We update this entry manually.
			it->virtualPosition -= tmpDelta;
			beginShiftPos = rmBeginPos + 1;
		}
		erase(tree._journalNodes, rmBeginPos, rmEndPos);
	}

    // Perform left-shift of the virtual positions.
    for (TIter it = begin(tree._journalNodes, Standard()) + beginShiftPos; it != end(tree._journalNodes, Standard()); ++it) {
        SEQAN_ASSERT_GEQ(it->virtualPosition, delta);
        it->virtualPosition -= delta;
    }
    //std::cerr << __FILE__ << ":" << __LINE__ << " -- " << tree << std::endl;

    SEQAN_ASSERT_TRUE(_checkSortedArrayTree(tree));
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNAL_JOURNAL_TREE_SORTED_ARRAY_H_
