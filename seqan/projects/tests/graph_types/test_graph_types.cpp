/*==========================================================================
 SeqAn - The Library for Sequence Analysis
 http://www.seqan.de 
 ===========================================================================
 Copyright (C) 2010
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.
 
 ===========================================================================
 Author: @@Your Name@@ <@@Your Email@@>
 ===========================================================================
 @@Description of what is tested here@@
 ===========================================================================
*/


// External STL
#include <iostream>
#include <fstream>
#include <string>

// Seqan
#include <seqan/graph_types.h>

// Test files
#include "test_graph_basic.h"
#include "test_graph_types.h"
#include "test_graph_iterators.h"
#include "test_graph_properties.h"
#include "test_graph_derived.h"
#include "test_graph_utils.h"


SEQAN_BEGIN_TESTSUITE(test_graph_types) {
    // Call Tests.
    SEQAN_CALL_TEST(test_graph_basics);
    SEQAN_CALL_TEST(test_graph_types);	
    SEQAN_CALL_TEST(test_graph_iterators);
    SEQAN_CALL_TEST(test_graph_properties);
    SEQAN_CALL_TEST(test_graph_derived);
    SEQAN_CALL_TEST(test_graph_utils);
	
    // Verify Checkpoints.
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_interface.h");

	// basic
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_idmanager.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_edgestump.h");

	// types
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_directed.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_undirected.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_automaton.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_wordgraph.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_tree.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_fragment.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_hmm.h");
	
	// iterators
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_iterator.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_iterator_vertex.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_iterator_outedge.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_iterator_adjacency.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_iterator_edge.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_iterator_bfs.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_iterator_dfs.h");
	
	// properties
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_property.h");
	
	// derived
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_oracle.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_impl_trie.h");
	
	// utils	
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_drawing.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/graph_types/graph_utility_parsing.h");
}
SEQAN_END_TESTSUITE
