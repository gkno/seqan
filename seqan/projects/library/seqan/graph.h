#ifndef SEQAN_HEADER_GRAPH_H
#define SEQAN_HEADER_GRAPH_H

// Seqan
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align/matrix_base.h>
#include <seqan/misc/misc_random.h>

#include <seqan/score.h>
#include <seqan/align.h>

// External / STL
#include <deque>
#include <map>
#include <set>
#include <queue>
#include <typeinfo>
#include <sstream> 


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/graph/graph_generated_forwards.h>
#endif


// Basic graph stuff
#include <seqan/graph/graph_base.h>
#include <seqan/graph/graph_idmanager.h>	// Id manager
#include <seqan/graph/graph_edgestump.h>	// EdgeStumps
#include <seqan/graph/graph_interface.h>	// Graph metafunctions

// Graph types
#include <seqan/graph/graph_impl_directed.h>	// Directed Graph
#include <seqan/graph/graph_impl_undirected.h>	// Undirected graph
#include <seqan/graph/graph_impl_automaton.h>	// Automaton
#include <seqan/graph/graph_impl_wordgraph.h>	// Specialized automaton: Word graph
#include <seqan/graph/graph_impl_tree.h>		// Tree
#include <seqan/graph/graph_impl_fragment.h>    // Fragment
#include <seqan/graph/graph_impl_align.h>		// Alignment Graph

// Graph iterators
#include <seqan/graph/graph_iterator.h>
#include <seqan/graph/graph_iterator_vertex.h>
#include <seqan/graph/graph_iterator_outedge.h>
#include <seqan/graph/graph_iterator_adjacency.h>
#include <seqan/graph/graph_iterator_edge.h>

// Graph property maps
#include <seqan/graph/graph_property.h>

// Specializations
#include <seqan/graph/graph_impl_oracle.h>		// Oracle
#include <seqan/graph/graph_impl_trie.h>		// Trie

// Graph drawing
// #include <seqan/find.h>
#include <seqan/graph/graph_drawing.h>

// Specialized iterators
#include <seqan/graph/graph_iterator_bfs.h>
#include <seqan/graph/graph_iterator_dfs.h>

// Graph algorithms
#include <seqan/graph/graph_algorithm_lis_his.h>
#include <seqan/graph/graph_algorithm.h>

// Alignment
#include <seqan/graph/graph_align_base.h>
#include <seqan/graph/graph_align_config.h>
#include <seqan/graph/graph_align_needleman_wunsch.h>
#include <seqan/graph/graph_align_gotoh.h>
#include <seqan/graph/graph_align_myers.h>
#include <seqan/graph/graph_align_hirschberg.h>
#include <seqan/graph/graph_align_smith_waterman.h>
#include <seqan/graph/graph_align_smith_waterman_clump.h>
#include <seqan/graph/graph_align_interface.h>

// T-Coffee
#include <seqan/graph/graph_align_tcoffee_base.h>
#include <seqan/graph/graph_align_tcoffee_kmer.h>
#include <seqan/graph/graph_align_tcoffee_distance.h>
#include <seqan/graph/graph_align_tcoffee_guidetree.h>
#include <seqan/graph/graph_align_tcoffee_library.h>
#include <seqan/graph/graph_align_tcoffee_progressive.h>

// Folding
#include <seqan/graph/graph_fold_nussinov.h>

// Utilities
#include <seqan/graph/graph_utility_alphabets.h>
#include <seqan/graph/graph_utility_parsing.h>
#include <seqan/graph/graph_utility_match_parsing.h>

// Refinement
#include <seqan/graph/graph_impl_interval_types.h>
#include <seqan/graph/graph_impl_interval_tree.h>
#include <seqan/graph/graph_algorithm_refine_scoring.h>
#include <seqan/graph/graph_algorithm_refine_fragment.h>
//#include <seqan/graph/graph_algorithm_refine_aligngraph.h>
#include <seqan/graph/graph_algorithm_refine_align.h>
#include <seqan/graph/graph_algorithm_refine_exact.h>
#include <seqan/graph/graph_algorithm_refine_inexact.h>

#endif //#ifndef SEQAN_HEADER_...
