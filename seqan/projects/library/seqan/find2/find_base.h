/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de
 ============================================================================
  Copyright (C) 2007-2010

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
  Basic definitions and types for the find module.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_BASE_H_
#define SEQAN_FIND2_FIND_BASE_H_

namespace seqan {

// Contains the state for finder and patterns.
struct _FindState {
    enum TState {
        STATE_EMPTY,           // Finder/pattern is empty.
        STATE_INITIAL,         // Finer/pattern has just been initialized.
        STATE_FOUND,           // Found the end position of a hit.
        STATE_NOTFOUND,        // No hit found, no more hits possible.
        STATE_BEGIN_FOUND,     // Found begin position.
        STATE_BEGIN_NOTFOUND,  // Found end but not begin, should not happen.
        STATE_NO_HIT           // Set manually to non-hit.
    };
};


template <typename TNeedle, typename TSpec>
struct Pattern;


template <typename THaystack, typename TSpec = Default>
struct Finder;


template <typename TPattern>
struct Needle;


template <typename TNeedle, typename TSpec>
struct Needle<Pattern<TNeedle, TSpec> > {
    typedef TNeedle Type;
};


// Metafunction for retrieving scoring schemes.
template <typename TPattern>
struct ScoringScheme;

}  // namespace seqan

#endif  // SEQAN_FIND2_FIND_BASE_H_
