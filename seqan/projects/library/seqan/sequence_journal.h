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
  Umbrella Header for the sequence_journal module.
  ==========================================================================*/

#ifndef SEQAN_HEADER_SEQUENCE_JOURNAL_H
#define SEQAN_HEADER_SEQUENCE_JOURNAL_H

// ============================================================================
// Prerequisites.
// ============================================================================

#include <algorithm>
#include <map>
#include <ostream>
#include <string>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

// ============================================================================
// Forwards.
// ============================================================================

#include <seqan/sequence_journal/sequence_journal_forwards.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/sequence_journal/sequence_journal_generated_forwards.h>
#endif

// ============================================================================
// Journaled Sequences.
// ============================================================================

#include <seqan/sequence_journal/journal_entry.h>
#include <seqan/sequence_journal/journal_entries_unbalanced_tree_node.h>
#include <seqan/sequence_journal/journal_entries_unbalanced_tree.h>
#include <seqan/sequence_journal/journal_entries_unbalanced_tree_iterator.h>
#include <seqan/sequence_journal/journal_entries_sorted_array.h>
#include <seqan/sequence_journal/sequence_journal.h>
#include <seqan/sequence_journal/sequence_journal_iterator.h>

// ============================================================================
// Incremental Indices.
// ============================================================================

// TODO(holtgrew): Port back the incremental index stuff to the new code.

#endif  // SEQAN_HEADER_SEQUENCE_JOURNAL_H
