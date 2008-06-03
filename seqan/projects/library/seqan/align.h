 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_ALIGN_H
#define SEQAN_HEADER_ALIGN_H

//____________________________________________________________________________
// prerequisites

#include <cmath>
#include <stack>

#include <seqan/sequence.h>
#include <seqan/score.h>

#include <seqan/map.h>

//____________________________________________________________________________

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/align/align_generated_forwards.h>
#endif

#include <seqan/align/gaps_base.h>
#include <seqan/align/gaps_iterator_base.h>

#include <seqan/align/gaps_array.h>
#include <seqan/align/gaps_sumlist.h>
#include <seqan/align/gaps_simple.h>

#include <seqan/align/align_base.h>
#include <seqan/align/align_cols_base.h>
#include <seqan/align/align_iterator_base.h>

#include <seqan/align/matrix_base.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

#include <seqan/align/align_dynprog.h>
#include <seqan/align/align_local_dynprog.h>
#include <seqan/align/hirschberg_set.h>
#include <seqan/align/align_myers.h>
#include <seqan/align/align_hirschberg.h>


#endif //#ifndef SEQAN_HEADER_...
