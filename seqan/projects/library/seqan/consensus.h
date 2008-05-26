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
  $Id: consensus.h 1901 2008-04-28 13:07:56Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CONSENSUS_H
#define SEQAN_HEADER_CONSENSUS_H

// Seqan
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align/matrix_base.h>
#include <seqan/misc/misc_random.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
//#include <seqan/align.h>

// External / STL
#include <deque>
#include <map>
#include <set>
#include <queue>
#include <typeinfo>
#include <sstream> 
#include <fstream> 
#include <limits>



#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/refinement.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/consensus/consensus_generated_forwards.h>
#endif

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/refinement/refinement_generated_forwards.h>
#endif

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/graph_types/graph_types_generated_forwards.h>
#endif

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/graph_align/graph_align_generated_forwards.h>
#endif

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/graph_algorithms/graph_algorithms_generated_forwards.h>
#endif

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/graph_msa/graph_msa_generated_forwards.h>
#endif




// Consensus tool
#include <seqan/consensus/graph_consensus_base.h>
#include <seqan/consensus/graph_consensus_ctgstore.h>
#include <seqan/consensus/graph_consensus_library.h>
#include <seqan/consensus/graph_consensus_readstore.h>
#include <seqan/consensus/graph_consensus_frgstore.h>
#include <seqan/consensus/graph_consensus_libstore.h>

#endif //#ifndef SEQAN_HEADER_...
