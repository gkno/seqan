/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010 by Freie Universitaet Berlin
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  Author: Tobias Rausch <rausch@embl.de>
  ===========================================================================
  Tests for the Graph MSA module.
  ===========================================================================
*/

#define SEQAN_ENABLE_CHECKPOINTS 0

// External / STL
#include <iostream>
#include <fstream>
#include <string>

// SeqAn
#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/graph_msa.h>
#include <seqan/misc/misc_random.h>

// Test files
#include "test_graph_msa_guide_tree.h"
#include "test_graph_tcoffee.h"

SEQAN_BEGIN_TESTSUITE(test_graph_msa)
{
    // Call Tests.
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_neighbour_joining);
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_upgma_weight_avg);
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_upgma_avg);
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_upgma_min);
    SEQAN_CALL_TEST(test_graph_msa_guide_tree_upgma_max);

    SEQAN_CALL_TEST(test_distances);
	SEQAN_CALL_TEST(test_libraries);
	SEQAN_CALL_TEST(test_external_libraries);
	SEQAN_CALL_TEST(test_triplet_extension);
	SEQAN_CALL_TEST(test_sop);
	SEQAN_CALL_TEST(test_progressive);
	SEQAN_CALL_TEST(test_reversable_fragments);	
	
    // Verify Checkpoints.
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_kmer.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_distance.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_guidetree.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_library.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_io.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_progressive.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_msa/graph_align_tcoffee_msa.h");
}
SEQAN_END_TESTSUITE

