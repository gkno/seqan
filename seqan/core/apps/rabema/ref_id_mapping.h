// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Consider for inclusion in SeqAn library.

#ifndef CORE_APPS_RABEMA_REF_ID_MAPPING_H_
#define CORE_APPS_RABEMA_REF_ID_MAPPING_H_

#include <seqan/sequence.h>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RefIdMapping
// ----------------------------------------------------------------------------

class RefIdMapping
{
public:
    String<unsigned> map;
    
    RefIdMapping() {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

inline unsigned length(RefIdMapping const & mapping)
{
    return length(mapping.map);
}

// ----------------------------------------------------------------------------
// Function rebuildMapping()
// ----------------------------------------------------------------------------

template <typename TTargetNameStore, typename TTargetNameStoreCache, typename TSourceNameStore>
void rebuildMapping(RefIdMapping & mapping,
                    TTargetNameStore const & targetNameStore,
                    TTargetNameStoreCache const & targetNameStoreCache,
                    TSourceNameStore const & sourceNameStore)
{
    clear(mapping.map);
    resize(mapping.map, length(sourceNameStore), maxValue<unsigned>());

    for (unsigned i = 0; i < length(sourceNameStore); ++i)
    {
        unsigned idx = 0;
        if (getIdByName(targetNameStore, sourceNameStore[i], idx, targetNameStoreCache))
            mapping.map[i] = idx;
    }
}

#endif  // #ifndef CORE_APPS_RABEMA_REF_ID_MAPPING_H_
