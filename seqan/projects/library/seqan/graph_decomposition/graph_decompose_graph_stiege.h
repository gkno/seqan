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
    EDGE_PROCESSED
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
        String<unsigned> & vertexTypes, // TODO(holtgrew): s/vertexTypes/vertexFlags/g, s/edgeTypes/edgeFlags/g
        String<unsigned> & edgeTypes,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g) // TODO(holtgrew): Remove const after fixing iterator for const-graphs
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
            setBit(property(vertexTypes, v), VERTEX_ISOLATED);
            std::cerr << "ISOLATED " << v << std::endl;
            continue;
        }

        // Ignore v if it has more than one edge or already has a token.
        if (deg > 1u)
            continue;
        if (isBitSet(getProperty(vertexTypes, v), VERTEX_TOKEN))
            continue;

        // Otherwise, it has degree 1 and no token.  It is the leaf of an
        // external tree, an end vertex.  Mark it so.
        setBit(property(vertexTypes, v), VERTEX_PERIPHERAL_TREE);
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
            TOutEdgeIterator itE(g, w);
            for (; !atEnd(itE); goNext(itE)) {
                if (getProperty(edgeTypes, value(itE)) == 0) {
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
            setBit(property(vertexTypes, x), VERTEX_TOKEN);
            clearBit(property(vertexTypes, w), VERTEX_TOKEN);

            // Label e as peripheral tree edge.
            std::cerr << "PERIPHERAL TREE {" << w << ", " << x << "}" << std::endl;
            setBit(property(edgeTypes, e), EDGE_PERIPHERAL_TREE);
            SEQAN_ASSERT_GT(unmarkedEdgeCounts[w], 0u);
            property(unmarkedEdgeCounts, w)  -= 1;
            SEQAN_ASSERT_GT(unmarkedEdgeCounts[x], 0u);
            property(unmarkedEdgeCounts, x) -= 1;

            // Set w to other end point of e.
            w = x;

            // Characterize w as peripheral tree vertex.
            std::cerr << "PERIPHERAL TREE " << w << std::endl;
            setBit(property(vertexTypes, w), VERTEX_PERIPHERAL_TREE);
        }
    }

    // Characterize vertices as roots of peripheral trees that are proper
    // components and border point of peripheral trees.
    for (TVertexIterator it(g); !atEnd(it); ++it)
    {
        TVertexDescriptor v = *it;

        if (!isBitSet(getProperty(vertexTypes, v), VERTEX_TOKEN))
            continue;
        clearBit(property(vertexTypes, v), VERTEX_TOKEN);

        if (getProperty(unmarkedEdgeCounts, v) == 0) {
            std::cerr << "PERIPHERAL ROOT " << v << std::endl;
            setBit(property(vertexTypes, v), VERTEX_PERIPHERAL_TREE_ROOT);
        } else {
            std::cerr << "PERIPHERAL BORDER " << v << std::endl;
            setBit(property(vertexTypes, v), VERTEX_PERIPHERAL_TREE_BORDER);
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

    // Build clusterTree.
    // for all vertices
    //   if isolated: add as isolated component
    //   if peripheral tree root or border point: add tree component, traverse tree and map vertices to component
}

} // namespace seqan

#endif // #ifndef SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_GRAPH_STIEGE_H_
