// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Implementation of the graph decomposition by Stiege.
//
// Relevant literature:
//
// Tarjan, R E. Depth-First Search and Linear Graph Algorithms. SIAM Journal
// on Computing 1, no. 2 (July 1972): 146.
// ==========================================================================

#ifndef SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_GRAPH_STIEGE_H_
#define SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_GRAPH_STIEGE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes, Structs, Enums, Tags
// ============================================================================

/**
.Enum.StandardDecompositionBlockType
..cat:Graph Decomposition
..summary:Types of building blocks in standard graph decomposition.
..signature:StandardDecompositionBlockType
..include:seqan/graph_decomposition.h
..value:BLOCK_NULL
..value:BLOCK_GRAPH
..value:BLOCK_IMPROPER_COMPONENT
..value:BLOCK_PROPER_COMPONENT
..value:BLOCK_PERIPHERAL_TREE
..value:BLOCK_STOPFREE_KERNEL
..value:BLOCK_INTERNAL_TREE
..value:BLOCK_SUBCOMPONENT
..value:BLOCK_BIBLOCK
 */

enum StandardDecompositionBlockType
{
    BLOCK_NULL,
    BLOCK_GRAPH,
    BLOCK_IMPROPER_COMPONENT,
    BLOCK_PROPER_COMPONENT,
    BLOCK_PERIPHERAL_TREE,
    BLOCK_STOPFREE_KERNEL,
    BLOCK_INTERNAL_TREE,
    BLOCK_SUBCOMPONENT,
    BLOCK_BIBLOCK
};

/**
.Enum.StandardDecompositionVertexType
..cat:Graph Decomposition
..summary:Types of vertices in standard graph decompositon.
..signature:StandardDecompositionVertexType
..include:seqan/graph_decomposition.h
..value:VERTEX_TOKEN
..value:VERTEX_MARK
..value:VERTEX_ISOLATED
..value:VERTEX_PERIPHERAL_TREE
..value:VERTEX_PERIPHERAL_TREE_ROOT
..value:VERTEX_PERIPHERAL_TREE_BORDER
..value:VERTEX_BIBLOCK
..value:VERTEX_INTERNAL_TREE
..value:VERTEX_INTERNAL_TREE_ROOT
 */

enum StandardDecompositionVertexType
{
    VERTEX_TOKEN,
    VERTEX_MARK,
    VERTEX_ISOLATED,
    VERTEX_PERIPHERAL_TREE,         // non-root of peripheral tree
    VERTEX_PERIPHERAL_TREE_ROOT,    // root of "isolated peripheral tree"
    VERTEX_PERIPHERAL_TREE_BORDER,  // border point/root of non-isolated ptree
    VERTEX_BIBLOCK,                 // in a biblock
    VERTEX_INTERNAL_TREE,           // in internal tree
    VERTEX_INTERNAL_TREE_ROOT       // root of internal tree
};

/**
.Enum.StandardDecompositionEdgeType
..cat:Graph Decomposition
..summary:Types of edges in standard graph decompositon.
..signature:StandardDecompositionEdgeType
..include:seqan/graph_decomposition.h
..value:EDGE_PERIPHERAL_TREE
..value:EDGE_INTERNAL_TREE
..value:EDGE_BIBLOCK
 */

enum StandardDecompositionEdgeType
{
    EDGE_PERIPHERAL_TREE,
    EDGE_INTERNAL_TREE,
    EDGE_BIBLOCK
};

/**
.Spec.Standard Graph Decomposition
..summary:Decomposition of a graph following Stiege's standard decomposition.
..general:Class.GraphDecomposition
..cat:Graph Decomposition
..signature:GraphDecomposition<TGraph, StandardDecomposition>
..param.TGraph:The type of the graph to decompose.
...type:Spec.Undirected Graph
..include:seqan/graph_decomposition.h

.Memfunc.Standard Graph Decomposition#GraphDecomposition
..class:Spec.Standard Graph Decomposition
..summary:Constructor
..signature:GraphDecomposition()
..signature:GraphDecomposition(g)
..param.g:The graph this decomposition is for.
...type:Spec.Undirected Graph
 */

struct StandardDecomposition_;
typedef Tag<StandardDecomposition_> StandardDecomposition;

template <typename TCargo, typename TSpec>
class GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition>
{
public:
    typedef Graph<Undirected<TCargo, TSpec> > TGraph_;
    typedef Graph<Tree<> > TDecompositionTree_;
    typedef String<StandardDecompositionBlockType> TClusterNodeCargo_;
    typedef typename VertexDescriptor<TDecompositionTree_>::Type TTreeVertexDescriptor_;
    typedef String<String<TTreeVertexDescriptor_> > TVertexToClusterMap_;
    typedef String<TTreeVertexDescriptor_> TEdgeToClusterMap_;

    // ------------------------------------------------------------------------
    // Member Variables
    // ------------------------------------------------------------------------

    // The Holder for the host, adhering to "Hosted Type" concept.
    Holder<TGraph_> data_host;

    TDecompositionTree_ _clusterTree;
    TClusterNodeCargo_ _clusterNodeCargoMap;
    TVertexToClusterMap_ _vertexToClusterMap;
    TEdgeToClusterMap_ _edgeToClusterMap;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    GraphDecomposition() {}
    
    GraphDecomposition(TGraph_ /*const*/ & graph)  // TODO(holtgrew): make const when holder is fixed
            : data_host(graph)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ClusterTree
// ----------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
struct ClusterTree<GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> >
{
    typedef GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> TDecomposition_;
    typedef typename TDecomposition_::TDecompositionTree_ Type;
};

template <typename TCargo, typename TSpec>
struct ClusterTree<GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> const>
{
    typedef GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> const TDecomposition_;
    typedef typename TDecomposition_::TDecompositionTree_ const Type;
};

// ----------------------------------------------------------------------------
// Metafunction ClusterNodeCargo
// ----------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
struct ClusterNodeCargo<GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> >
{
    typedef GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> TDecomposition_;
    typedef typename TDecomposition_::TClusterNodeCargo_ Type;
};

template <typename TCargo, typename TSpec>
struct ClusterNodeCargo<GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> const>
{
    typedef GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> const TDecomposition_;
    typedef typename TDecomposition_::TClusterNodeCargo_ const Type;
};

// ----------------------------------------------------------------------------
// Metafunction VertexToClusterMap
// ----------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
struct VertexToClusterMap<GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> >
{
    typedef GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> TDecomposition_;
    typedef typename TDecomposition_::TVertexToClusterMap_ Type;
};

template <typename TCargo, typename TSpec>
struct VertexToClusterMap<GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> const>
{
    typedef GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> const TDecomposition_;
    typedef typename TDecomposition_::TVertexToClusterMap_ const Type;
};

// ----------------------------------------------------------------------------
// Metafunction EdgeToClusterMap
// ----------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
struct EdgeToClusterMap<GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> >
{
    typedef GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> TDecomposition_;
    typedef typename TDecomposition_::TEdgeToClusterMap_ Type;
};

template <typename TCargo, typename TSpec>
struct EdgeToClusterMap<GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> const>
{
    typedef GraphDecomposition<Graph<Undirected<TCargo, TSpec> >, StandardDecomposition> const TDecomposition_;
    typedef typename TDecomposition_::TEdgeToClusterMap_ const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function decomposeGraph()
// ----------------------------------------------------------------------------

template <typename TGraphCargo, typename TGraphSpec, typename TVertexDescriptor, typename TComponents, typename TComponentToBlock>
void
_classifyAsIsolated(
        GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> & gd,
        TVertexDescriptor const & v,
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    typedef GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> TGraphDecomposition;
    typedef typename ClusterTree<TGraphDecomposition>::Type TClusterTree;
    typedef typename VertexDescriptor<TClusterTree>::Type TTreeVertexDescriptor;

    // Mark the connected component node in tree as improper component.
    TTreeVertexDescriptor x = componentToBlock[getProperty(components, v)];
    assignProperty(clusterNodeCargoMap(gd), x, BLOCK_IMPROPER_COMPONENT);
    appendValue(property(vertexToClusterMap(gd), v), x);
}

template <typename TVertexFlags, typename TEdgeFlags, typename TGraphCargo, typename TGraphSpec, typename TVertexDescriptor, typename TGraph, typename TComponents, typename TComponentToBlock>
void
_collectPeripheralTreeAndClassify(
        TVertexFlags /*const*/ & vertexFlags,  // we use the token flag for DFS, TODO(holtgrew): Better use external map?
        TEdgeFlags const & edgeFlags,
        GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> & gd,
        TVertexDescriptor const & v,
        TGraph /*const*/ & g,
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef String<TVertexDescriptor> TStack;
    typedef GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> TGraphDecomposition;
    typedef typename ClusterTree<TGraphDecomposition>::Type TTree;
    typedef typename VertexDescriptor<TTree>::Type TTreeVertexDescriptor;

    TStack stack;

    TTreeVertexDescriptor blockNode = componentToBlock[getProperty(components, v)];
    TTreeVertexDescriptor peripheralTreeNode = addChild(clusterTree(gd), blockNode);
    resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
    assignProperty(clusterNodeCargoMap(gd), peripheralTreeNode, BLOCK_PERIPHERAL_TREE);

    // Perform a DFS using only peripheral tree vertices.
    appendValue(stack, v);
    while (!empty(stack)) {
        TVertexDescriptor x = back(stack);
        SEQAN_ASSERT(isBitSet(getProperty(vertexFlags, x), VERTEX_PERIPHERAL_TREE));
        setBit(property(vertexFlags, x), VERTEX_MARK);
        appendValue(property(vertexToClusterMap(gd), x), peripheralTreeNode);
        eraseBack(stack);

        for (TOutEdgeIterator itE(g, x); !atEnd(itE); goNext(itE)) {
            TEdgeDescriptor e = *itE;
            TVertexDescriptor y = getTarget(e);
            if (y == x)
                y = getSource(e);
            if (isBitSet(getProperty(vertexFlags, y), VERTEX_MARK))
                continue;
            if (isBitSet(getProperty(vertexFlags, y), VERTEX_PERIPHERAL_TREE)) {
                SEQAN_ASSERT(isBitSet(getProperty(edgeFlags, e), EDGE_PERIPHERAL_TREE));
                assignProperty(edgeToClusterMap(gd), e, peripheralTreeNode);
                appendValue(stack, y);
            }
        }
    }

    // Release token and mark on root again so we can use the flag again later on.
    clearBit(property(vertexFlags, v), VERTEX_TOKEN);
    clearBit(property(vertexFlags, v), VERTEX_MARK);
}

template <typename TGraphCargo, typename TGraphSpec, typename TComponents, typename TComponentToBlock>
void
_decomposeGraphStiegeFindPeripheralTrees(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> & gd,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g, // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    // typename TTreeCargo, typename TTreeSpec, typename TBlockDescriptorsMap, typename TTreeVertexDescriptor, 

    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    // We need a vertex map with unmarked edge counts for the algorithm to be
    // O(m).  In star graphs, the complexity is O(n^2) = O(m^2) otherwise.
    //
    // TODO(holtgrew): Why is degree() not in O(1)?
    String<unsigned> unmarkedEdgeCounts;
    resizeVertexMap(g, unmarkedEdgeCounts);
    for (TVertexIterator it(g); !atEnd(it); ++it) {
        assignProperty(unmarkedEdgeCounts, *it, degree(g, *it));
        // std::cerr << "unmarkedEdgeCounts[" << *it << "] = " << degree(g, *it) << std::endl;
    }

    // Perform token-propagation algorithm.
    for (TVertexIterator it(g); !atEnd(it); ++it) {
        TVertexDescriptor v = value(it);

        // Easy case: v is isolated
        unsigned deg = degree(g, v);
        if (deg == 0u) {
            setBit(property(vertexFlags, v), VERTEX_ISOLATED);
            _classifyAsIsolated(gd, v, components, componentToBlock);
            continue;
        }

        // Ignore v if it has more than one edge or already has a token.
        if (deg > 1u)
            continue;
        if (isBitSet(getProperty(vertexFlags, v), VERTEX_TOKEN))
            continue;

        // Otherwise, it has degree 1 and no token.  It is the leaf of an
        // external tree, an end vertex.  Mark it so.
        setBit(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE);
        // std::cerr << "PERIPHERAL TREE " << v << std::endl;
        // Follow along this edge as long as the next vertex and continue as
        // long as the vertex is incident with only one unmarked edges,
        // marking the edges as peripheral tree edges.
        TVertexDescriptor w = v;
        // std::cerr << "w = " << w << std::endl;
        while (getProperty(unmarkedEdgeCounts, w) == 1u) {
#if SEQAN_ENABLE_DEBUG
            bool assignedE = false;
#endif  // #if SEQAN_ENABLE_DEBUG
            TEdgeDescriptor e;
            for (TOutEdgeIterator itE(g, w); !atEnd(itE); goNext(itE)) {
                if (getProperty(edgeFlags, value(itE)) == 0) {
#if SEQAN_ENABLE_DEBUG
                    assignedE = true;
#endif  // #if SEQAN_ENABLE_DEBUG
                    e = *itE;  // descriptor of one unmarked edge {w, x}
                    break;
                }
            }
            SEQAN_ASSERT(assignedE);
            TVertexDescriptor x = targetVertex(g, e);
            if (x == w)  // TODO(holtgrew): Necessary?
                x = sourceVertex(g, e);

            // Move token to x, joining any two conceptional tokens.
            // std::cerr << "MOVE TOKEN FROM " << w << " TO " << x << std::endl;
            setBit(property(vertexFlags, x), VERTEX_TOKEN);
            clearBit(property(vertexFlags, w), VERTEX_TOKEN);

            // Label e as peripheral tree edge.
            // std::cerr << "PERIPHERAL TREE {" << w << ", " << x << "}" << std::endl;
            setBit(property(edgeFlags, e), EDGE_PERIPHERAL_TREE);
            SEQAN_ASSERT_GT(unmarkedEdgeCounts[w], 0u);
            property(unmarkedEdgeCounts, w)  -= 1;
            SEQAN_ASSERT_GT(unmarkedEdgeCounts[x], 0u);
            property(unmarkedEdgeCounts, x) -= 1;

            // Set w to other end point of e.
            w = x;

            // Characterize w as peripheral tree vertex.
            // std::cerr << "PERIPHERAL TREE " << w << std::endl;
            setBit(property(vertexFlags, w), VERTEX_PERIPHERAL_TREE);
        }
    }

    // Characterize vertices as roots of peripheral trees that are proper
    // components and border point of peripheral trees.
    for (TVertexIterator it(g); !atEnd(it); ++it)
    {
        TVertexDescriptor v = *it;

        // std::cerr << "isBitSet(getProperty(vertexFlags, " << v << "), VERTEX_TOKEN) == " << isBitSet(getProperty(vertexFlags, v), VERTEX_TOKEN) << std::endl;
        if (!isBitSet(getProperty(vertexFlags, v), VERTEX_TOKEN))
            continue;
        clearBit(property(vertexFlags, v), VERTEX_TOKEN);

        if (getProperty(unmarkedEdgeCounts, v) == 0) {
            // std::cerr << "PERIPHERAL ROOT " << v << std::endl;
            setBit(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE);
            setBit(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE_ROOT);
        } else {
            // std::cerr << "PERIPHERAL BORDER " << v << std::endl;
            setBit(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE);
            setBit(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE_BORDER);
        }

        // Collect peripheral tree.
        _collectPeripheralTreeAndClassify(vertexFlags, edgeFlags, gd, v, g, components, componentToBlock);
    }
}

template <typename TGraphCargo, typename TGraphSpec, typename TStack, typename TVertexDescriptor, typename TComponents, typename TComponentToBlock>
void
_decomposeGraphStiegeBiconnect(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> & gd,
        String<unsigned> & lowPt,
        String<unsigned> & number,
        String<unsigned> & edgeComponentIds,
        unsigned & i,
        unsigned & j,  // component id
        TStack & stack,
        TVertexDescriptor const & v,
        TVertexDescriptor const & u,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g, // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    // std::cerr << "BIBLOCK(" << v << ", " << u << ")" << std::endl;

    typedef GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> TGraphDecomposition;
    typedef typename ClusterTree<TGraphDecomposition>::Type TClusterTree;
    typedef typename VertexDescriptor<TClusterTree>::Type TTreeVertexDescriptor;
    typedef typename Iterator<TClusterTree, OutEdgeIterator>::Type TTreeOutEdgeIterator;
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Value<TStack>::Type TTriple;

    i += 1;
    // std::cerr << "LOWPT[" << v << "]  := " << i << std::endl;
    // std::cerr << "NUMBER[" << v << "] := " << i << std::endl;
    assignProperty(number, v, i);
    assignProperty(lowPt, v, i);


    // The subcomponent vertex.  It is initialized with the tree root
    // since the tree root can never be a subcomponent vertex.  It is
    // really created assigned on the creation of the first biblock.
    TTreeVertexDescriptor subComponentVertex = root(clusterTree(gd));
            
    // For w in the adjacency list of v...
    for (TOutEdgeIterator itE(g, v); !atEnd(itE); goNext(itE)) {
        TEdgeDescriptor e = *itE;
        TVertexDescriptor w = getTarget(e);
        if (w == v)
            w = getSource(e);

        // Assign edge to the given component.
        assignProperty(edgeComponentIds, e, getProperty(components, v));

        // If w is not yet numbered then...
        // std::cerr << "getProperty(number, " << w << ") {== " << getProperty(number, w) << "} == 0u" << std::endl;
        if (getProperty(number, w) == 0u) {
            // Skip vertices which are on the peripheral tree and not border points.
            if (isBitSet(property(vertexFlags, w), VERTEX_PERIPHERAL_TREE) &&
                !isBitSet(property(vertexFlags, w), VERTEX_PERIPHERAL_TREE_BORDER))
                continue;

            // std::cerr << "PUSH(" << v << ", " << w << ")" << std::endl;
            appendValue(stack, TTriple(v, w, e));

            // Recursive call.
            _decomposeGraphStiegeBiconnect(vertexFlags, edgeFlags, gd, lowPt, number, edgeComponentIds, i, j, stack, w, v, g, components, componentToBlock);

            assignProperty(lowPt, v, _min(getProperty(lowPt, v), getProperty(lowPt, w)));

            // std::cerr << "LOWPT[" << w << "] == " << getProperty(lowPt, w) << std::endl;
            // std::cerr << "NUMBER[" << v << "] == " << getProperty(number, v) << std::endl;
            if (getProperty(lowPt, w) >= getProperty(number, v)) {
                // Unroll stack [(u1, u2, e)] until NUMBER[u1] <= NUMBER[v].
                SEQAN_ASSERT_NOT(empty(stack));
                // std::cerr << "stop condition: " << getProperty(number, back(stack).i1) << " > " << getProperty(number, v) << std::endl;
                bool wasComponent = false;  // True if top edges are biblock edges.
                TTreeVertexDescriptor biblockKernelVertex;
                if (!empty(stack) && getProperty(number, back(stack).i1) > getProperty(number, v)) {
                    // Start new biblock (and subcomponent if necessary).
                    wasComponent = true;
                    // Create biblock vertex and subcomponent vertex if necessary.
                    if (subComponentVertex == root(clusterTree(gd))) {
                        TTreeVertexDescriptor componentVertex = componentToBlock[getProperty(components, v)];
                        // The first child of the component vertex is the stopfree
                        // kernel vertex.  It is added after identifying the
                        // connected components. Search for it.
#if SEQAN_ENABLE_DEBUG
                        bool found = false;
#endif  // #if SEQAN_ENABLE_DEBUG
                        TTreeVertexDescriptor stopfreeKernelVertex;
                        for (TTreeOutEdgeIterator itF(clusterTree(gd), componentVertex); !atEnd(itF); goNext(itF)) {
                            if (childVertex(clusterTree(gd), *itF) == componentVertex)
                                continue;
                            if (getProperty(clusterNodeCargoMap(gd), childVertex(clusterTree(gd), *itF)) != BLOCK_STOPFREE_KERNEL)
                                continue;
                            stopfreeKernelVertex = childVertex(clusterTree(gd), *itF);
#if SEQAN_ENABLE_DEBUG
                            found = true;
#endif  // #if SEQAN_ENABLE_DEBUG
                            break;
                        }
                        SEQAN_ASSERT(found);
                        // Add a new vertex for the subcomponent.
                        // std::cerr << "STARTING SUBCOMPONENT FROM " << v << " stopfreeKernelVertex == " << stopfreeKernelVertex << std::endl;
                        subComponentVertex = addChild(clusterTree(gd), stopfreeKernelVertex);
                        // std::cerr << "ADD CHILD " << stopfreeKernelVertex << " == " << subComponentVertex << std::endl;
                        resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
                        assignProperty(clusterNodeCargoMap(gd), subComponentVertex, BLOCK_SUBCOMPONENT);
                    }
                    // std::cerr << "STARTING BIBLOCK FROM " << v << std::endl;
                    // Add a new vertex to the tree for the biblock we are about to add.
                    // std::cerr << "ADD CHILD " << subComponentVertex << std::endl;
                    biblockKernelVertex = addChild(clusterTree(gd), subComponentVertex);
                    resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
                    assignProperty(clusterNodeCargoMap(gd), biblockKernelVertex, BLOCK_BIBLOCK);
                }
                while (!empty(stack) && getProperty(number, back(stack).i1) > getProperty(number, v)) {
                    // if (!wasComponent) {
                    //     std::cerr << "Starting BIBLOCK " << j << std::endl;
                    // }
                    // std::cerr << "BIBLOCK " << back(stack).i1 << std::endl;
                    // std::cerr << "BIBLOCK " << back(stack).i2 << std::endl;
                    // std::cerr << "BIBLOCK {" << getSource(back(stack).i3) << ", " << getTarget(back(stack).i3) << "}" << std::endl;
                    // Update flags for edges and vertices.
                    setBit(property(edgeFlags, back(stack).i3), EDGE_BIBLOCK);
                    setBit(property(vertexFlags, back(stack).i1), VERTEX_BIBLOCK);
                    setBit(property(vertexFlags, back(stack).i2), VERTEX_BIBLOCK);
                    // Assign edge to biblock and add biblock to the blocks for the vertices.
                    assignProperty(edgeToClusterMap(gd), back(stack).i3, biblockKernelVertex);
                    appendValue(property(vertexToClusterMap(gd), back(stack).i1), biblockKernelVertex);
                    // Only adding vertex for the i1-end of the edge to avoid duplicates.
                    // std::cerr << "POP() == (" << back(stack).i1 << ", " << back(stack).i2 << ") " << __LINE__ << std::endl;
                    eraseBack(stack);
                    // if (!empty(stack))
                    //     std::cerr << "stop condition: " << getProperty(number, back(stack).i1) << " > " << getProperty(number, v) << std::endl;
                }
                // delete (v, w) from stack
                SEQAN_ASSERT_NOT(empty(stack));
                SEQAN_ASSERT_EQ(back(stack).i1, v);
                SEQAN_ASSERT_EQ(back(stack).i2, w);
                if (wasComponent) {
                    // std::cerr << "BIBLOCK " << back(stack).i1 << std::endl;
                    // std::cerr << "BIBLOCK " << back(stack).i2 << std::endl;
                    // std::cerr << "BIBLOCK {" << getSource(back(stack).i3) << ", " << getTarget(back(stack).i3) << "}" << std::endl;
                    // Update flags for edges and vertices.
                    setBit(property(edgeFlags, back(stack).i3), EDGE_BIBLOCK);
                    setBit(property(vertexFlags, back(stack).i1), VERTEX_BIBLOCK);
                    setBit(property(vertexFlags, back(stack).i2), VERTEX_BIBLOCK);
                    // Assign edge to biblock and add biblock to the blocks for the vertices.
                    assignProperty(edgeToClusterMap(gd), back(stack).i3, biblockKernelVertex);
                    appendValue(property(vertexToClusterMap(gd), back(stack).i1), biblockKernelVertex);
                    // Only adding vertex for the i1-end of the edge to avoid duplicates.
                    // std::cerr << "POP() == (" << back(stack).i1 << ", " << back(stack).i2 << ") " << __LINE__ << std::endl;
                    j += 1; // id of next component
                } else {
                    // Mark edge and incident vertices as belonging to the
                    // internal tree.  The node in clusterTree will be added
                    // when collecting internal trees.
                    setBit(property(vertexFlags, back(stack).i1), VERTEX_INTERNAL_TREE);
                    setBit(property(vertexFlags, back(stack).i2), VERTEX_INTERNAL_TREE);
                    setBit(property(edgeFlags, back(stack).i3), EDGE_INTERNAL_TREE);
                    // std::cerr << "INTERNAL TREE {" << back(stack).i1 << ", " << back(stack).i2 << "}" << std::endl;
                }
                eraseBack(stack);
            }
        } else if (getProperty(number, w) < getProperty(number, v) && w != u) {
            // std::cerr << "PUSH(" << v << ", " << w << ")" << std::endl;
            appendValue(stack, TTriple(v, w, e));
            assignProperty(lowPt, v, _min(getProperty(lowPt, v), getProperty(number, w)));
        }
    }
}

template <typename TGraphCargo, typename TGraphSpec, typename TVertexDescriptor, typename TComponents, typename TComponentToBlock>
void
_decomposeGraphStiegeCollectInternalTrees(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> & gd,
        TVertexDescriptor v,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g, // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    // std::cout << "COLLECT-INTERNAL-TREE(" << v << ")" << std::endl;

    // Skip vertex if not part of the internal tree.
    if (!isBitSet(getProperty(vertexFlags, v), VERTEX_INTERNAL_TREE))
        return;
    // std::cout << "  internal tree" << std::endl;
    if (isBitSet(property(vertexFlags, v), VERTEX_TOKEN))
        return;
    // std::cout << "  no token" << std::endl;

    // If we reach here, v is the root of an internal tree which we have not
    // yet collect.
    setBit(property(vertexFlags, v), VERTEX_INTERNAL_TREE_ROOT);

    // Create new node in cluster tree for internal tree.
    //
    // Start by getting the node vertex of the component.
    typedef GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> TGraphDecomposition;
    typedef typename ClusterTree<TGraphDecomposition>::Type TClusterTree;
    typedef typename VertexDescriptor<TClusterTree>::Type TTreeVertexDescriptor;
    typedef typename Iterator<TClusterTree, OutEdgeIterator>::Type TTreeOutEdgeIterator;
    TTreeVertexDescriptor componentVertex = componentToBlock[getProperty(components, v)];
    // The first child of the component vertex is the stopfree kernel.  Get it
#if SEQAN_ENABLE_DEBUG
    bool found = false;
#endif  // #if SEQAN_ENABLE_DEBUG
    TTreeVertexDescriptor stopfreeKernelVertex;
    for (TTreeOutEdgeIterator itF(clusterTree(gd), componentVertex); !atEnd(itF); goNext(itF)) {
        if (childVertex(clusterTree(gd), *itF) == componentVertex)
            continue;
        stopfreeKernelVertex = childVertex(clusterTree(gd), *itF);
#if SEQAN_ENABLE_DEBUG
        found = true;
#endif  // #if SEQAN_ENABLE_DEBUG
    }
    SEQAN_ASSERT(found);
    // Create new internal tree vertex below stopfree kernel vertex.
    TVertexDescriptor internalTreeVertex = addChild(clusterTree(gd), stopfreeKernelVertex);
    resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
    assignProperty(clusterNodeCargoMap(gd), internalTreeVertex, BLOCK_INTERNAL_TREE);

    // Assign root vertex to internal tree vertex.
    appendValue(property(vertexToClusterMap(gd), v), internalTreeVertex);

    // Stack for DFS search.
    typedef String<TVertexDescriptor> TStack;
    TStack stack;
    appendValue(stack, v);

    // Perform DFS on all internal tree edges to vertices that are not marked
    // with a token yet.
    while (!empty(stack)) {
        v = back(stack);
        setBit(property(vertexFlags, v), VERTEX_TOKEN);
        // std::cout << "putting token on " << v << std::endl;
        eraseBack(stack);
        
        for (TOutEdgeIterator itE(g, v); !atEnd(itE); goNext(itE)) {
            TEdgeDescriptor e = *itE;
            if (!isBitSet(getProperty(edgeFlags, e), EDGE_INTERNAL_TREE))
                continue;  // Skip non-internal-tree edges.
            
            TVertexDescriptor u = getTarget(e);
            if (u == v)
                u = getSource(e);
            SEQAN_ASSERT(isBitSet(getProperty(vertexFlags, u), VERTEX_INTERNAL_TREE));
            if (isBitSet(getProperty(vertexFlags, u), VERTEX_TOKEN))
                continue;  // Skip already marked vertices.

            // Add this edge and other end vertex to the internal tree.
            assignProperty(edgeToClusterMap(gd), e, internalTreeVertex);
            appendValue(property(vertexToClusterMap(gd), u), internalTreeVertex);
            
            appendValue(stack, u);
        }
    }
}

// This is a slightly extended version of Tarjan's classic biconnected
// component search.  The function _decomposeGraphStiegeBiconnect corresponds
// to the routine BICONNECT in (Tarjan, 1973).
template <typename TGraphCargo, typename TGraphSpec, typename TComponents, typename TComponentToBlock>
void
_decomposeGraphStiegeFindStopfreeKernels(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> & gd,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g, // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef String<Triple<TVertexDescriptor, TVertexDescriptor, TEdgeDescriptor> > TStack;

    unsigned i = 0;
    // lowPt[v] is the number of the smallest vertex in the DFS palmtree
    // reachable from v.
    String<unsigned> lowPt;
    resizeVertexMap(g, lowPt);
    // number[v] is the rank of v when visited in the DFS.
    String<unsigned> number;
    resizeVertexMap(g, number, 0);
    // The stack of vertex-pairs.
    TStack stack;
    // Edge map of component ids.
    unsigned j = 1;
    String<unsigned> edgeComponentIds;  // TODO(holtgrew): Not written out, though... why? Superflous?
    resizeEdgeMap(g, edgeComponentIds, maxValue<unsigned>());

    // For each node v which is not a border point or not a vertex of a
    // peripheral tree.
    for (TVertexIterator it(g); !atEnd(it); ++it) {
        TVertexDescriptor v = *it;

        // Skip vertices which are on the peripheral tree and not border points.
        if (isBitSet(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE) &&
            !isBitSet(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE_BORDER))
            continue;
        // Skip marked vertices.
        if (isBitSet(property(vertexFlags, v), VERTEX_TOKEN))
            continue;

        // Start biconnected components DFS.
        clear(stack);
        _decomposeGraphStiegeBiconnect(vertexFlags, edgeFlags, gd, lowPt, number, edgeComponentIds, i, j, stack, v, maxValue<TVertexDescriptor>(), g, components, componentToBlock);
    }

    // Finally, collect internal trees.  The edges have already been marked
    // appropriately.
    for (TVertexIterator it(g); !atEnd(it); ++it) {
        TVertexDescriptor v = *it;
        _decomposeGraphStiegeCollectInternalTrees(vertexFlags, edgeFlags, gd, v, g, components, componentToBlock);
    }
}

template <typename TGraphCargo, typename TGraphSpec>
void
decomposeGraph(GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> & gd,
               Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g) // TODO(holtgrew): Remove const after fixing iterator for const-graphs
{
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename VertexDescriptor<TGraph>::Type TGraphVertexDescriptor;
    typedef GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> TGraphDecomposition;
    typedef typename ClusterTree<TGraphDecomposition>::Type TTree;
    typedef typename VertexDescriptor<TTree>::Type TTreeVertexDescriptor;

    // Perform the standard decomposition of the given graph.  We build the
    // cluster tree at the same time.
    
    // Edge and vertex type markers.
    // TODO(holtgrew): These are written out nowhere!
    String<__uint32> vertexFlags;
    resizeVertexMap(g, vertexFlags, 0u);
    String<__uint32> edgeFlags;
    resizeEdgeMap(g, edgeFlags, 0u);

    // Initialize decomposition tree.
    //
    // Add root vertex which corresponds to the graph.
    clear(clusterTree(gd));
    TTreeVertexDescriptor r = addVertex(clusterTree(gd));
    assignRoot(clusterTree(gd), r);
    resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
    assignProperty(clusterNodeCargoMap(gd), r, BLOCK_GRAPH);
    // Initialize maps from vertices and edges into the decomposition tree.
    resizeVertexMap(g, vertexToClusterMap(gd));
    resizeEdgeMap(g, edgeToClusterMap(gd));

    // Collect preliminary stopfree kernels from cluster tree.
    String<TTreeVertexDescriptor> tentativeStopfreeKernels;

    // Find connected components.
    String<unsigned> componentIds;
    unsigned componentCount = connectedComponents(g, componentIds);
    // std::cerr << "vertex -> cc" << std::endl;
    // for (unsigned i = 0; i < numVertices(g); ++i) {
    //     std::cerr << i << " -> " << getProperty(componentIds, i) << std::endl;
    // }
    String<TTreeVertexDescriptor> componentToBlock;
    resize(componentToBlock, componentCount);
    for (unsigned i = 0; i < componentCount; ++i) {
        TTreeVertexDescriptor x = addChild(clusterTree(gd), r);
        resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
        assignProperty(clusterNodeCargoMap(gd), x, BLOCK_PROPER_COMPONENT);  // <-- will update later to improper if isolated
        componentToBlock[i] = x;

        // Each connected component contains at most one stopfree kernel.
        // This will always be the first child of that building block's tree
        // node.  If there is no stopfree kernel, it will be removed later on.
        TTreeVertexDescriptor y = addChild(clusterTree(gd), x);
        resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
        assignProperty(clusterNodeCargoMap(gd), y, BLOCK_STOPFREE_KERNEL);
        // std::cerr << "cc node y == " << y << " -> " << i << std::endl;
        appendValue(tentativeStopfreeKernels, y);
    }

    // Find isolated vertices and peripheral trees.
    _decomposeGraphStiegeFindPeripheralTrees(vertexFlags, edgeFlags, gd, g, componentIds, componentToBlock);

    // Find stopfree kernels and their internal structure.
    _decomposeGraphStiegeFindStopfreeKernels(vertexFlags, edgeFlags, gd, g, componentIds, componentToBlock);

    typedef typename Iterator<String<TTreeVertexDescriptor> >::Type TTentativeIterator;
    for (TTentativeIterator it = begin(tentativeStopfreeKernels); it != end(tentativeStopfreeKernels); ++it) {
        if (numChildren(clusterTree(gd), *it) > 0u)
            continue;
        removeVertex(clusterTree(gd), *it);
    }
}

// ----------------------------------------------------------------------------
// Function operator<<();  For stream output.
// ----------------------------------------------------------------------------

template <typename TStream, typename TGraphCargo, typename TGraphSpec>
void writeDecompositionTree(
        TStream & stream,
        GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> & gd)
{
    typedef GraphDecomposition<Graph<Undirected<TGraphCargo, TGraphSpec> >, StandardDecomposition> TGraphDecomposition;
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > TGraph;
    typedef typename ClusterTree<TGraphDecomposition>::Type TTree;
    typedef typename VertexDescriptor<TTree>::Type TTreeVertexDescriptor;

    // Collect vertices in elementary building blocks.
    String<String<TTreeVertexDescriptor> > blocksForVertex;
    resizeVertexMap(clusterTree(gd), blocksForVertex);

    typedef typename Iterator<TGraph, VertexIterator>::Type TGraphVertexIterator;
    TGraphVertexIterator itV(host(gd));
    for (; !atEnd(itV); goNext(itV)) {
        typedef typename Iterator<String<TTreeVertexDescriptor>, Rooted>::Type TIter;
        for (TIter it(begin(property(vertexToClusterMap(gd), *itV))); !atEnd(it); goNext(it)) {
            appendValue(property(blocksForVertex, *it), *itV);
        }
    }

    // DOT output.
    stream << "digraph G {" << std::endl;

    stream << std::endl;
    stream << "/* Nodes */" << std::endl;
    stream << std::endl;

    static char const * blockTypeNames[] = {
        "(null)", "graph", "improper component", "proper component",
        "peripheral tree", "stopfree kernel", "internal tree",
        "subcomponent", "biblock"
    };

    typedef typename Iterator<TTree, VertexIterator>::Type TTreeVertexIterator;
    for (TTreeVertexIterator itV(clusterTree(gd)); !atEnd(itV); goNext(itV)) {
        CharString label = blockTypeNames[getProperty(clusterNodeCargoMap(gd), *itV)];
        if (getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_BIBLOCK ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_PERIPHERAL_TREE ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_INTERNAL_TREE ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_IMPROPER_COMPONENT ) {
            append(label, "\\n(");
            typedef typename Iterator<String<TTreeVertexDescriptor>, Standard>::Type TIter;
            TIter itBegin = begin(property(blocksForVertex, *itV));
            TIter itEnd = end(property(blocksForVertex, *itV));
            for (TIter it = itBegin; it != itEnd; ++it) {
                if (it != itBegin)
                    append(label, ", ");
                std::stringstream ss;
                ss << *it;
                append(label, ss.str());
            }
            stream << *itV << " [label = \"" << label << ")\", shape = box];" << std::endl;
        } else {
            stream << *itV << " [label = \"" << label << "\"];" << std::endl;
        }
    }

    stream << std::endl;
    stream << "/* Edges */" << std::endl;
    stream << std::endl;

    typedef typename Iterator<TTree, DfsPreorder>::Type TTreeDfsIterator;
    typedef typename Iterator<TTree, OutEdgeIterator>::Type TOutEdgeIterator;
    for (TTreeDfsIterator itV(clusterTree(gd), root(clusterTree(gd))); !atEnd(itV); goNext(itV)) {
        for (TOutEdgeIterator itE(clusterTree(gd), *itV); !atEnd(itE); goNext(itE)) {
            if (childVertex(clusterTree(gd), *itE) == *itV)
                continue;
            stream << parentVertex(clusterTree(gd), *itE) << " -> " << childVertex(clusterTree(gd), *itE) << std::endl;
        }
    }

    stream << std::endl;
    stream << "}" << std::endl;
}

} // namespace seqan

#endif // #ifndef SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_GRAPH_STIEGE_H_