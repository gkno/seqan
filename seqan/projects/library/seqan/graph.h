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


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/graph/graph_generated_forwards.h>
#endif


#include <seqan/graph/graph_stack.h>
#include <seqan/align/matrix_base.h>

// Basic metafunctions
#include <seqan/graph/graph_base.h>

// Property maps
#include <seqan/graph/graph_property.h>

// Id manager
#include <seqan/graph/graph_idmanager.h>

// Graph iterators
#include <seqan/graph/graph_iterator.h>
#include <seqan/graph/graph_vertexiterator.h>
#include <seqan/graph/graph_outedgeiterator.h>
#include <seqan/graph/graph_adjacencyiterator.h>
#include <seqan/graph/graph_edgeiterator.h>
#include <seqan/graph/graph_bfsiterator.h>
#include <seqan/graph/graph_dfsiterator.h>

// Graph interface
#include <seqan/graph/graph_interface.h>

// Graph implementations
// EdgeList
#include <seqan/graph/graph_edgestump.h>
#include <seqan/graph/graph_impl_edgelist.h>
// Automaton
#include <seqan/graph/graph_edgeautomaton.h>
#include <seqan/graph/graph_impl_automaton.h>
// Special automaton: Word graph
#include <seqan/graph/graph_impl_wordgraph.h>
// Special automaton: Oracle
#include <seqan/graph/graph_impl_oracle.h>
// Special automaton: Trie
#include <seqan/graph/graph_impl_trie.h>


// Graph drawing
#include <seqan/graph/graph_drawing.h>

// Graph algorithms
#include <seqan/graph/graph_algorithm.h>



#endif //#ifndef SEQAN_HEADER_...
