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

// TODO(holtgrew): How to document Enums?
enum StiegeUndirectedBlockType
{
    BLOCK_NULL,
    BLOCK_GRAPH,
    BLOCK_IMPROPER_COMPONENT,
    BLOCK_PROPER_COMPONENT,
    BLOCK_ACYCLIC_COMPONENT,
    BLOCK_CYCLIC_COMPONENT,
    BLOCK_PERIPHERAL_TREE,
    BLOCK_STOPFREE_KERNEL,
    BLOCK_INTERNAL_TREE,
    BLOCK_SUBCOMPONENT,
    BLOCK_BIBLOCK
};

enum StiegeVertexType
{
    VERTEX_TOKEN,
    VERTEX_MARK,
    VERTEX_ISOLATED,
    VERTEX_PERIPHERAL_TREE,         // non-root of peripheral tree
    VERTEX_PERIPHERAL_TREE_ROOT,    // root of "isolated peripheral tree"
    VERTEX_PERIPHERAL_TREE_BORDER   // border point/root of non-isolated ptree
};

enum StiegeEdgeType
{
    EDGE_PERIPHERAL_TREE,
    EDGE_INTERNAL_TREE,
    EDGE_SUBCOMPONENT,  // not in biblock!
    // Two b-markers, so we can also store processed/unprocessed state.
    EDGE_BIBLOCK,
    EDGE_PROCESSED,
    EDGE_LABEL_B
};

/**
.Class.UndirectedBuildingBlock
..cat:Graph Decomposition
..summary:Building block of an undirected graph after Stiege.
..signature:UndirectedBuildingBlock

.Memvar.UndirectedBuildingBlock#blockType
..summary:Type of an building block in the decomposition of undirected graphs.
..class:Class.UndirectedBuildingBlock
..remarks:The type is a value from the set {$BLOCK_NULL$, $BLOCK_GRAPH$, $BLOCK_IMPROPER_COMPONENT$, $BLOCK_PROPER_COMPONENT$, $BLOCK_ACYCLIC_COMPONENT$, $BLOCK_CYCLIC_COMPONENT$, $BLOCK_PERIPHERAL_TREE$, $BLOCK_STOPFREE_KERNEL$, $BLOCK_INTERNAL_TREE$, $BLOCK_SUBCOMPONENT$, $BLOCK_BIBLOCK$}.
 */
class UndirectedBuildingBlock
{
public:
    UndirectedBuildingBlock()
            : blockType(BLOCK_NULL)
    {}

    UndirectedBuildingBlock(StiegeUndirectedBlockType const & blockType_ )
            : blockType(blockType_)
    {}
    
    StiegeUndirectedBlockType blockType;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function decomposeGraphStiege()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Move these functions into misc_bitmasks.h?
template <typename TWord>
inline
void
setBitTo(TWord & word, unsigned index, bool value)
{
    // See http://www-graphics.stanford.edu/~seander/bithacks.html#ConditionalSetOrClearBitsWithoutBranching
    word = (word & ~(1 << index)) | (-value & (1<< index));
}

template <typename TWord>
inline
void
setBit(TWord & word, unsigned index)
{
    word |= (1 << index);
}

template <typename TWord>
inline
void
clearBit(TWord & word, unsigned index)
{
    word &= ~(1 << index);
}

template <typename TWord>
inline
void
clearBits(TWord & word)
{
    word = 0;
}

template <typename TWord>
inline
bool
isBitSet(TWord const & word, unsigned index)
{
    return (word & (1 << index)) != 0;
}

template <typename TGraphCargo, typename TGraphSpec>
void
_decomposeGraphStiegeFindPeripheralTrees(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g) // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
{
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    
    // We need a vertex map with unmarked edge counts for the algorithm to be
    // O(m).  In star graphs, the complexity is O(n^2) = O(m^2) otherwise.
    //
    // TODO(holtgrew): Wh is degree() not in O(1)?
    String<unsigned> unmarkedEdgeCounts;
    resizeVertexMap(g, unmarkedEdgeCounts);
    for (TVertexIterator it(g); !atEnd(it); ++it) {
        assignProperty(unmarkedEdgeCounts, *it, degree(g, *it));
        std::cerr << "unmarkedEdgeCounts[" << *it << "] = " << degree(g, *it) << std::endl;
    }

    // Perform token-propagation algorithm.
    for (TVertexIterator it(g); !atEnd(it); ++it)
    {
        TVertexDescriptor v = value(it);

        // Easy case: v is isolated
        unsigned deg = degree(g, v);
        if (deg == 0u) {
            setBit(property(vertexFlags, v), VERTEX_ISOLATED);
            std::cerr << "ISOLATED " << v << std::endl;
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
        std::cerr << "PERIPHERAL TREE " << v << std::endl;
        // Follow along this edge as long as the next vertex and continue as
        // long as the vertex is incident with only one unmarked edges,
        // marking the edges as peripheral tree edges.
        TVertexDescriptor w = v;
        std::cerr << "w = " << w << std::endl;
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
            SEQAN_ASSERT_TRUE(assignedE);
            TVertexDescriptor x = targetVertex(g, e);
            if (x == w)  // TODO(holtgrew): Necessary?
                x = sourceVertex(g, e);

            // Move token to x, joining any two conceptional tokens.
            setBit(property(vertexFlags, x), VERTEX_TOKEN);
            clearBit(property(vertexFlags, w), VERTEX_TOKEN);

            // Label e as peripheral tree edge.
            std::cerr << "PERIPHERAL TREE {" << w << ", " << x << "}" << std::endl;
            setBit(property(edgeFlags, e), EDGE_PERIPHERAL_TREE);
            SEQAN_ASSERT_GT(unmarkedEdgeCounts[w], 0u);
            property(unmarkedEdgeCounts, w)  -= 1;
            SEQAN_ASSERT_GT(unmarkedEdgeCounts[x], 0u);
            property(unmarkedEdgeCounts, x) -= 1;

            // Set w to other end point of e.
            w = x;

            // Characterize w as peripheral tree vertex.
            std::cerr << "PERIPHERAL TREE " << w << std::endl;
            setBit(property(vertexFlags, w), VERTEX_PERIPHERAL_TREE);
        }
    }

    // Characterize vertices as roots of peripheral trees that are proper
    // components and border point of peripheral trees.
    for (TVertexIterator it(g); !atEnd(it); ++it)
    {
        TVertexDescriptor v = *it;

        if (!isBitSet(getProperty(vertexFlags, v), VERTEX_TOKEN))
            continue;
        clearBit(property(vertexFlags, v), VERTEX_TOKEN);

        if (getProperty(unmarkedEdgeCounts, v) == 0) {
            std::cerr << "PERIPHERAL ROOT " << v << std::endl;
            setBit(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE_ROOT);
        } else {
            std::cerr << "PERIPHERAL BORDER " << v << std::endl;
            setBit(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE_BORDER);
        }
    }
}

template <typename TEdgeDescriptor, typename TVertexDescriptor, typename TGraphCargo, typename TGraphSpec>
void
_decomposeGraphStiegeInternalStructure(
        int & counter,
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        TEdgeDescriptor const & entryEdge,
        TVertexDescriptor const & vertex,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g) // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
{
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ TGraph;
    // typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    // typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    setBit(property(vertexFlags, vertex), VERTEX_MARK);  // mark vertex
    int oldCounter = counter;

    // For each non-labeled exit edge e of vertex.
    for (TOutEdgeIterator itE(g, vertex); !atEnd(itE); goNext(itE)) {
        TEdgeDescriptor e = *itE;
        if (isBitSet(getProperty(edgeFlags, e), EDGE_LABEL_B) || entryEdge == e)
            continue;

        // v is the other end vertex of e.
        TVertexDescriptor v = getTarget(e);
        if (vertex == v)
            v = getSource(e);

        // if the other end vertex v of e is unmarked...
        if (!isBitSet(getProperty(vertexFlags, v), VERTEX_MARK)) {
            _decomposeGraphStiegeInternalStructure(counter, vertexFlags, edgeFlags, e, v, g);
            if (counter > oldCounter) {
                // for each b-labeled exit edge f not yet processed
                // TODO(holtgrew): b-labeled edge incident to vertex?
                for (TOutEdgeIterator itF(g, v); !atEnd(itF); goNext(itF)) {
                    TEdgeDescriptor f = *itF;
                    // Skip b-labeled edges, processed edges and entry edge.
                    if (!isBitSet(property(edgeFlags, f), EDGE_LABEL_B) ||
                        !isBitSet(property(edgeFlags, f), EDGE_PROCESSED) ||
                        entryEdge != f)
                        continue;

                    counter -= 1;
                    setBit(property(edgeFlags, f), EDGE_PROCESSED);
                }
            }
        } else {
            setBit(property(edgeFlags, e), EDGE_LABEL_B);
            counter += 1;
        }
    }

    // TODO(holtgrew): Classify edges and vertex and integrate into new data structure.
}

template <typename TGraphCargo, typename TGraphSpec>
void
_decomposeGraphStiegeFindStopfreeKernels(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g) // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
{
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    // Fo reach node v which is not a border point or not a vertex of a
    // peripheral tree.
    for (TVertexIterator it(g); !atEnd(it); ++it)
    {
        TVertexDescriptor v = *it;

        // Skip vertices which are on the peripheral tree and not border points.
        if (isBitSet(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE) &&
            !isBitSet(property(vertexFlags, v), VERTEX_PERIPHERAL_TREE_BORDER))
            continue;
        // Skip marked vertices.
        if (isBitSet(property(vertexFlags, v), VERTEX_MARK))
            continue;

        // Mark v.
        setBit(property(vertexFlags, v), VERTEX_MARK);

        // For each unlabeled edge e incident with v.
        for (TOutEdgeIterator itE(g, v); !atEnd(itE); goNext(itE)) {
            TEdgeDescriptor e = *itE;
            if (isBitSet(getProperty(edgeFlags, e), EDGE_LABEL_B))
                continue;

            int counter = 0;  // TODO(holtgrew): __int64 would be safer but is probably unnecessary.
            // w = other end vertex of e
            TVertexDescriptor w = getTarget(e);
            if (w == v)
                w = getSource(e);

            // Call DFS recursion.
            _decomposeGraphStiegeInternalStructure(counter, vertexFlags, edgeFlags, e, w, g);

            if (counter > 0) {
                // for each b-labeled edge f not yet processed
                // TODO(holtgrew): b-labeled edge incident to w?
                for (TOutEdgeIterator itF(g, v); !atEnd(itF); goNext(itF)) {
                    TEdgeDescriptor f = *itF;
                    // The inversion of "is b-labeled and not processed" is
                    // "is not b-labeled or processed."
                    if (!isBitSet(getProperty(edgeFlags, f), EDGE_LABEL_B) ||
                        isBitSet(getProperty(edgeFlags, f), EDGE_PROCESSED))
                        continue;

                    counter -= 1;
                    // Mark edge as processed.
                    setBit(property(edgeFlags, f), EDGE_PROCESSED);
                }
            }

            // TODO(holtgrew): Classify edges and node v and integrate into new data structure.
        }
    }
}

/**
.Function.decomposeGraphStiege
..cat:Graph Decomposition
..summary:Decompose a graph using Stiege's digraph decomposition.
..signature:decomposeGraphStiege(clusterTree, blockDescriptors, blockMap, g)
..param.clusterTree:The resulting hierarchy stored in the tree.
...type:Spec.Tree
..param.blockDescriptors: An external property map assigning each vertex of the tree its @Class.UndirectedBuildingBlock@ description.
..param.blockMap:An edge property map assigning each edge to a vertex in $clusterTree$.
...remarks:Each vertex corresponds to a building block.
..param.g:The graph to decompose.
...type:Spec.Undirected Graph
..remarks:The graph will be decomposed into 1 improper connected components, 2 proper connected components, 2.1 acyclic component, 2.2 cyclic components, 2.2.1 peripheral trees, 2.2.2 stopfree kernels, 2.2.2.1 internal trees, 2.2.2.2 subcomponents, 2.2.2.2.1 biblocks (biconnected blocks)
..remarks:This algorithm follows "Günther Stiege. An Algorithm for Finding the Connectivity Structure of Undirected Graphs. Hildesheimer Informatik-Berichte, 13/97, May 1997."
 */
template <typename TTreeCargo, typename TTreeSpec, typename TBlockDescriptorMapValue, typename TBlockDescriptorMapSpec, typename TBlockMapValue, typename TBlockMapSpec, typename TGraphCargo, typename TGraphSpec>
void
decomposeGraphStiege(Graph<Tree<TTreeCargo, TTreeSpec> > & clusterTree,
                     String<TBlockDescriptorMapValue, TBlockDescriptorMapSpec> & blockDescriptors,
                     String<TBlockMapValue, TBlockMapSpec> & blockMap,
                     Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g) // TODO(holtgrew): Remove const after fixing iterator for const-graphs
{
    typedef Graph<Tree<TTreeCargo, TTreeSpec> > TTree;
    typedef typename VertexDescriptor<TTree>::Type TTreeVertexDescriptor;
    
    // Initialization.
    //
    // The root vertex corresponds to the graph.
    clear(clusterTree);
    TTreeVertexDescriptor root = addVertex(clusterTree);
    assignRoot(clusterTree, root);
    resizeVertexMap(clusterTree, blockMap);
    resizeVertexMap(clusterTree, blockDescriptors);
    assignProperty(blockDescriptors, root, UndirectedBuildingBlock(BLOCK_GRAPH));
    // Edge and vertex type markers.
    // TODO(holtgrew): resizeVertexMap and resizeEdgeMap should accept an optional prototype!
    String<unsigned> vertexTypes;
    resizeVertexMap(g, vertexTypes);
    String<unsigned> edgeTypes;
    resizeEdgeMap(g, edgeTypes);

    write(std::cout, g, DotDrawing());

    // Find isolated vertices and peripheral trees.
    _decomposeGraphStiegeFindPeripheralTrees(vertexTypes, edgeTypes, g);

    // Find stopfree kernels and their internal structure.
    _decomposeGraphStiegeFindStopfreeKernels(vertexTypes, edgeTypes, g);

    // Build clusterTree.
    // for all vertices
    //   if isolated: add as isolated component
    //   if peripheral tree root or border point: add tree component, traverse tree and map vertices to component
}

} // namespace seqan

#endif // #ifndef SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_GRAPH_STIEGE_H_
