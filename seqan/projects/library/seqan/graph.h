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


// Basic metafunctions
#include <seqan/graph/graph_base.h>

// Id manager
#include <seqan/graph/graph_idmanager.h>

// Various EdgeStumps
#include <seqan/graph/graph_edgestump.h>	// Directed graph
#include <seqan/graph/graph_edgestumpu.h>	// Undirected graph
#include <seqan/graph/graph_edgestumpt.h>	// Tree
#include <seqan/graph/graph_edgestumpa.h>	// Automaton

// Graph interface and property maps
#include <seqan/graph/graph_interface.h>	// Basic graph stuff
#include <seqan/graph/graph_property.h>		// Property maps

// Graph types
#include <seqan/graph/graph_impl_edgelist.h>	// Directed Graph
#include <seqan/graph/graph_impl_edgelistu.h>	// Undirected graph
#include <seqan/graph/graph_impl_automaton.h>	// Automaton
#include <seqan/graph/graph_impl_tree.h>		// Tree
#include <seqan/graph/graph_impl_wordgraph.h>	// Specialized automaton: Word graph

// Graph iterators
#include <seqan/graph/graph_iterator.h>
#include <seqan/graph/graph_vertexiterator.h>
#include <seqan/graph/graph_outedgeiterator.h>
#include <seqan/graph/graph_adjacencyiterator.h>
#include <seqan/graph/graph_edgeiterator.h>
#include <seqan/graph/graph_bfsiterator.h>
#include <seqan/graph/graph_dfsiterator.h>

// Special automaton: Oracle
#include <seqan/graph/graph_impl_oracle.h>
// Special automaton: Trie
#include <seqan/graph/graph_impl_trie.h>

// Graph drawing
#include <seqan/graph/graph_drawing.h>

// Graph algorithms
#include <seqan/graph/graph_algorithm.h>



#endif //#ifndef SEQAN_HEADER_...
