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
  Code for journal entries.
  ==========================================================================*/

#ifndef SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRY_H_
#define SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRY_H_

namespace seqan {

// ============================================================================
// Enums, Classes
// ============================================================================

enum SegmentSource {
    SOURCE_NULL,
    SOURCE_ORIGINAL,
    SOURCE_PATCH
};


template <typename _TPos, typename _TSize>
struct JournalEntry
{
    typedef _TPos TPos;
    typedef _TSize TSize;

    // Flag for where the segment comes from.
    SegmentSource segmentSource;
    // Position in the original string or the insertion buffer,
    // depending on segmentSource.
    TPos physicalPosition;
    // Position in the virtual string.
    TPos virtualPosition;
    // Length of the segment.
    TSize length;

    JournalEntry()
            : segmentSource(SOURCE_NULL),
              physicalPosition(0),
              virtualPosition(0),
              length(0)
    {
		SEQAN_CHECKPOINT;
    }

    JournalEntry(SegmentSource const & _segmentSource,
                 TPos _physicalPosition,
                 TPos _virtualPosition,
                 TSize _length)
            : segmentSource(_segmentSource),
              physicalPosition(_physicalPosition),
              virtualPosition(_virtualPosition),
              length(_length)
    {
		SEQAN_CHECKPOINT;
    }
};


template <typename TPos, typename TSize>
struct JournalEntryLtByVirtualPos
{
    bool operator()(JournalEntry<TPos, TSize> const & a,
                    JournalEntry<TPos, TSize> const & b) const
    {
        SEQAN_CHECKPOINT;
        return a.virtualPosition < b.virtualPosition;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TPos, typename TSize>
struct Size<JournalEntry<TPos, TSize> >
{
    typedef TSize Type;
};

template <typename TPos, typename TSize>
struct Size<JournalEntry<TPos, TSize> const>
        : Size<JournalEntry<TPos, TSize> > {};

template <typename TPos, typename TSize>
struct Position<JournalEntry<TPos, TSize> >
{
    typedef TPos Type;
};

template <typename TPos, typename TSize>
struct Position<JournalEntry<TPos, TSize> const>
        : Position<JournalEntry<TPos, TSize> > {};

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TPos, typename TSize>
TStream & operator<<(TStream & stream, JournalEntry<TPos, TSize> const & entry)
{
    return stream << "{segmentSource=" << entry.segmentSource
                  << ", virtualPosition=" << entry.virtualPosition
                  << ", physicalPosition=" << entry.physicalPosition
                  << ", length=" << entry.length
                  << "}";
}

}  // namespace seqan

#endif  // SEQAN_SEQUENCE_JOURNALED_JOURNAL_ENTRY_H_

