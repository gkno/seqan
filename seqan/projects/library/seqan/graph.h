#ifndef SEQAN_HEADER_GRAPH_H
#define SEQAN_HEADER_GRAPH_H

// External / STL
#include <deque>
#include <set>
#include <queue>
#include <typeinfo>

// Seqan
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align/matrix_base.h>

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

// Graph iterators
#include <seqan/graph/graph_iterator.h>
#include <seqan/graph/graph_vertexiterator.h>
#include <seqan/graph/graph_outedgeiterator.h>
#include <seqan/graph/graph_adjacencyiterator.h>
#include <seqan/graph/graph_edgeiterator.h>

// Graph property maps
#include <seqan/graph/graph_property.h>

// Special automatons: Just create functions
#include <seqan/graph/graph_impl_oracle.h>		// Oracle
#include <seqan/graph/graph_impl_trie.h>		// Trie

// Graph drawing
#include <seqan/graph/graph_drawing.h>

// Specialized iterators
#include <seqan/graph/graph_bfsiterator.h>
#include <seqan/graph/graph_dfsiterator.h>

// Graph algorithms
#include <seqan/graph/graph_algorithm.h>

#endif //#ifndef SEQAN_HEADER_...
