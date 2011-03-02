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

// TODO(holtgrew): How to document Enums?
enum StiegeUndirectedBlockType
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

enum StiegeVertexType
{
    VERTEX_TOKEN,
    VERTEX_ISOLATED,
    VERTEX_PERIPHERAL_TREE,         // non-root of peripheral tree
    VERTEX_PERIPHERAL_TREE_ROOT,    // root of "isolated peripheral tree"
    VERTEX_PERIPHERAL_TREE_BORDER,  // border point/root of non-isolated ptree
    VERTEX_BIBLOCK,                 // in a biblock
    VERTEX_INTERNAL_TREE,           // in internal tree
    VERTEX_INTERNAL_TREE_ROOT       // root of internal tree
};

enum StiegeEdgeType
{
    EDGE_PERIPHERAL_TREE,
    EDGE_INTERNAL_TREE,
    EDGE_BIBLOCK
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

template <typename TBlockDescriptorsMap, typename TTreeVertexDescriptor, typename TVertexDescriptor, typename TComponents, typename TComponentToBlock>
void
_classifyAsIsolated(
        TBlockDescriptorsMap & blockDescriptor,
        String<String<TTreeVertexDescriptor> > & vertexBlocks,
        TVertexDescriptor const & v,
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    // Mark the connected component node in tree as improper component.
    TTreeVertexDescriptor x = componentToBlock[getProperty(components, v)];
    assignProperty(blockDescriptor, x, BLOCK_IMPROPER_COMPONENT);
    appendValue(property(vertexBlocks, v), x);
}

template <typename TTreeCargo, typename TTreeSpec, typename TBlockDescriptors, typename TTreeVertexDescriptor, typename TVertexFlags, typename TEdgeFlags, typename TVertexDescriptor, typename TGraph, typename TComponents, typename TComponentToBlock>
void
_collectPeripheralTreeAndClassify(
        Graph<Tree<TTreeCargo, TTreeSpec> > & clusterTree,
        TBlockDescriptors & blockDescriptor,
        String<String<TTreeVertexDescriptor> > & vertexBlocks,
        String<TTreeVertexDescriptor> & edgeBlock, 
        TVertexFlags /*const*/ & vertexFlags,  // we use the token flag for DFS, TODO(holtgrew): Better use external map?
        TEdgeFlags const & edgeFlags,
        TVertexDescriptor const & v,
        TGraph /*const*/ & g,
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef String<TVertexDescriptor> TStack;

    TStack stack;

    TTreeVertexDescriptor blockNode = componentToBlock[getProperty(components, v)];
    TTreeVertexDescriptor peripheralTreeNode = addChild(clusterTree, blockNode);
    resizeVertexMap(clusterTree, blockDescriptor);
    assignProperty(blockDescriptor, peripheralTreeNode, BLOCK_PERIPHERAL_TREE);

    // Perform a DFS using only peripheral tree vertices.
    appendValue(stack, v);
    while (!empty(stack)) {
        TVertexDescriptor x = back(stack);
        SEQAN_ASSERT_TRUE(isBitSet(getProperty(vertexFlags, x), VERTEX_PERIPHERAL_TREE));
        setBit(property(vertexFlags, x), VERTEX_TOKEN);
        appendValue(property(vertexBlocks, x), peripheralTreeNode);
        eraseBack(stack);

        for (TOutEdgeIterator itE(g, x); !atEnd(itE); goNext(itE)) {
            TEdgeDescriptor e = *itE;
            TVertexDescriptor y = getTarget(e);
            if (y == x)
                y = getSource(e);
            if (isBitSet(getProperty(vertexFlags, y), VERTEX_TOKEN))
                continue;
            if (isBitSet(getProperty(vertexFlags, y), VERTEX_PERIPHERAL_TREE)) {
                SEQAN_ASSERT_TRUE(isBitSet(getProperty(edgeFlags, e), EDGE_PERIPHERAL_TREE));
                assignProperty(edgeBlock, e, peripheralTreeNode);
                appendValue(stack, y);
            }
        }
    }

    // Release token on root again so we can use the flag again later on.
    clearBit(property(vertexFlags, v), VERTEX_TOKEN);
}

template <typename TTreeCargo, typename TTreeSpec, typename TBlockDescriptorsMap, typename TTreeVertexDescriptor, typename TGraphCargo, typename TGraphSpec, typename TComponents, typename TComponentToBlock>
void
_decomposeGraphStiegeFindPeripheralTrees(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        Graph<Tree<TTreeCargo, TTreeSpec> > & clusterTree,
        TBlockDescriptorsMap & blockDescriptor,
        String<String<TTreeVertexDescriptor> > & vertexBlocks,
        String<TTreeVertexDescriptor> & edgeBlock, 
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g, // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
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
        // std::cerr << "unmarkedEdgeCounts[" << *it << "] = " << degree(g, *it) << std::endl;
    }

    // Perform token-propagation algorithm.
    for (TVertexIterator it(g); !atEnd(it); ++it) {
        TVertexDescriptor v = value(it);

        // Easy case: v is isolated
        unsigned deg = degree(g, v);
        if (deg == 0u) {
            setBit(property(vertexFlags, v), VERTEX_ISOLATED);
            _classifyAsIsolated(blockDescriptor, vertexBlocks, v, components, componentToBlock);
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
            SEQAN_ASSERT_TRUE(assignedE);
            TVertexDescriptor x = targetVertex(g, e);
            if (x == w)  // TODO(holtgrew): Necessary?
                x = sourceVertex(g, e);

            // Move token to x, joining any two conceptional tokens.
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
        _collectPeripheralTreeAndClassify(clusterTree, blockDescriptor, vertexBlocks, edgeBlock, vertexFlags, edgeFlags, v, g, components, componentToBlock);
    }
}

template <typename TClusterTree, typename TBlockDescriptorsMap, typename TVertexBlocksMap, typename TEdgeBlockMap, typename TVertexDescriptor, typename TStack, typename TGraphCargo, typename TGraphSpec, typename TComponents, typename TComponentToBlock>
void
_decomposeGraphStiegeBiconnect(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        TClusterTree & clusterTree,
        TBlockDescriptorsMap & blockDescriptor,
        TVertexBlocksMap & vertexBlocks,
        TEdgeBlockMap & edgeBlock, 
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
    // CONTINUE HERE with structure extraction:
    // --> create one node in tree for stopfree kernel
    //     --> add one node per subcomponent
    //         --> add one per biconnected block

    typedef typename VertexDescriptor<TClusterTree>::Type TTreeVertexDescriptor;
    typedef typename Iterator<TClusterTree, OutEdgeIterator>::Type TTreeOutEdgeIterator;
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename Value<TStack>::Type TTriple;

    i += 1;
    assignProperty(number, v, i);
    assignProperty(lowPt, v, i);

    // For w in the adjacency list of v...
    for (TOutEdgeIterator itE(g, v); !atEnd(itE); goNext(itE)) {
        TEdgeDescriptor e = *itE;
        TVertexDescriptor w = getTarget(e);
        if (w == v)
            w = getSource(e);

        // Assign edge to the given component.
        assignProperty(edgeComponentIds, e, getProperty(components, v));

        // The subcomponent vertex.  It is initialized with the tree root
        // since the tree root can never be a subcomponent vertex.  It is
        // really created assigned on the creation of the first biblock.
        TTreeVertexDescriptor subComponentVertex = root(clusterTree);
            
        // If w is not yet numbered then...
        if (getProperty(number, w) == 0u) {
            // Skip vertices which are on the peripheral tree and not border points.
            if (isBitSet(property(vertexFlags, w), VERTEX_PERIPHERAL_TREE) &&
                !isBitSet(property(vertexFlags, w), VERTEX_PERIPHERAL_TREE_BORDER))
                continue;

            // std::cerr << "PUSH(" << v << ", " << w << ")" << std::endl;
            appendValue(stack, TTriple(v, w, e));

            _decomposeGraphStiegeBiconnect(vertexFlags, edgeFlags, clusterTree, blockDescriptor, vertexBlocks, edgeBlock, lowPt, number, edgeComponentIds, i, j, stack, w, v, g, components, componentToBlock);
            assignProperty(lowPt, v, _min(getProperty(lowPt, v), getProperty(lowPt, w)));
            if (getProperty(lowPt, w) >= getProperty(number, v)) {
                // Unroll stack [(u1, u2, e)] until NUMBER[u1] <= NUMBER[v].
                SEQAN_ASSERT_NOT(empty(stack));
                // std::cerr << "stop condition: " << getProperty(number, back(stack).i1) << " > " << getProperty(number, v) << std::endl;
                bool wasComponent = false;  // True if top edges are biblock edges.
                TTreeVertexDescriptor biblockKernelVertex;
                // TODO(holtgrew): What about doubly marked vertices? Only mark if source of edge?
                while (!empty(stack) && getProperty(number, back(stack).i1) > getProperty(number, v)) {
                    // Create biblock vertex and subcomponent vertex if necessary.
                    if (subComponentVertex == root(clusterTree)) {
                        TTreeVertexDescriptor componentVertex = componentToBlock[getProperty(components, v)];
                        // The first child of the component vertex is the stopfree
                        // kernel vertex.  It is added after identifying the
                        // connected components. Search for it.
#if SEQAN_ENABLE_DEBUG
                        bool found = false;
#endif  // #if SEQAN_ENABLE_DEBUG
                        TTreeVertexDescriptor stopfreeKernelVertex;
                        for (TTreeOutEdgeIterator itF(clusterTree, componentVertex); !atEnd(itF); goNext(itF)) {
                            if (childVertex(clusterTree, *itF) == componentVertex)
                                continue;
                            stopfreeKernelVertex = childVertex(clusterTree, *itF);
#if SEQAN_ENABLE_DEBUG
                            found = true;
#endif  // #if SEQAN_ENABLE_DEBUG
                            break;
                        }
                        SEQAN_ASSERT_TRUE(found);
                        // Add a new vertex for the subcomponent.
                        subComponentVertex = addChild(clusterTree, componentVertex);
                        resizeVertexMap(clusterTree, blockDescriptor);
                        assignProperty(blockDescriptor, subComponentVertex, BLOCK_SUBCOMPONENT);
                    }
                    // Add a new vertex to the tree for the biblock we are about to add.
                    biblockKernelVertex = addChild(clusterTree, subComponentVertex);
                    resizeVertexMap(clusterTree, blockDescriptor);
                    assignProperty(blockDescriptor, biblockKernelVertex, BLOCK_BIBLOCK);

                    // if (!wasComponent) {
                    //     std::cerr << "Starting BIBLOCK " << j << std::endl;
                    // }
                    wasComponent = true;
                    // std::cerr << "BIBLOCK " << back(stack).i1 << std::endl;
                    // std::cerr << "BIBLOCK " << back(stack).i2 << std::endl;
                    // std::cerr << "BIBLOCK {" << getSource(back(stack).i3) << ", " << getTarget(back(stack).i3) << "}" << std::endl;
                    // Update flags for edges and vertices.
                    setBit(property(edgeFlags, back(stack).i3), EDGE_BIBLOCK);
                    setBit(property(vertexFlags, back(stack).i1), VERTEX_BIBLOCK);
                    setBit(property(vertexFlags, back(stack).i2), VERTEX_BIBLOCK);
                    // Assign edge to biblock and add biblock to the blocks for the vertices.
                    assignProperty(edgeBlock, back(stack).i3, biblockKernelVertex);
                    appendValue(property(vertexBlocks, back(stack).i1), biblockKernelVertex);
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
                    assignProperty(edgeBlock, back(stack).i3, biblockKernelVertex);
                    appendValue(property(vertexBlocks, back(stack).i1), biblockKernelVertex);
                    // Only adding vertex for the i1-end of the edge to avoid duplicates.
                    // std::cerr << "POP() == (" << back(stack).i1 << ", " << back(stack).i2 << ") " << __LINE__ << std::endl;
                    j += 1; // id of next component
                } else {
                    // Mark edge and incident vertices as belonging to the
                    // internal tree.  The node in clusterTree will be added
                    // when collecting internal trees.
                    setBit(property(edgeFlags, back(stack).i1), VERTEX_INTERNAL_TREE);
                    setBit(property(edgeFlags, back(stack).i2), VERTEX_INTERNAL_TREE);
                    setBit(property(edgeFlags, back(stack).i3), EDGE_INTERNAL_TREE);
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

template <typename TClusterTree, typename TBlockDescriptorsMap, typename TVertexBlocksMap, typename TEdgeBlockMap, typename TVertexDescriptor, typename TGraphCargo, typename TGraphSpec, typename TComponents, typename TComponentToBlock>
void
_decomposeGraphStiegeCollectInternalTrees(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        TClusterTree & clusterTree,
        TBlockDescriptorsMap & blockDescriptor,
        TVertexBlocksMap & vertexBlocks,
        TEdgeBlockMap & edgeBlock,
        TVertexDescriptor v,
        Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g, // TODO(holtgrew): Uncomment const after fixing iterator for const-graphs
        TComponents const & components,
        TComponentToBlock const & componentToBlock)
{
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    // Skip vertex if not part of the internal tree.
    if (!isBitSet(getProperty(vertexFlags, v), VERTEX_INTERNAL_TREE))
        return;
    setBit(property(vertexFlags, v), VERTEX_TOKEN);

    // If we reach here, v is the root of an internal tree which we have not
    // yet collect.
    setBit(property(vertexFlags, v), VERTEX_INTERNAL_TREE_ROOT);

    // Create new node in cluster tree for internal tree.
    //
    // Start by getting the node vertex of the component.
    typedef typename VertexDescriptor<TClusterTree>::Type TTreeVertexDescriptor;
    typedef typename Iterator<TClusterTree, OutEdgeIterator>::Type TTreeOutEdgeIterator;
    TTreeVertexDescriptor componentVertex = componentToBlock[getProperty(components, v)];
    // The first child of the component vertex is the stopfree kernel.  Get it
#if SEQAN_ENABLE_DEBUG
    bool found = false;
#endif  // #if SEQAN_ENABLE_DEBUG
    TTreeVertexDescriptor stopfreeKernelVertex;
    for (TTreeOutEdgeIterator itF(clusterTree, componentVertex); !atEnd(itF); goNext(itF)) {
        if (childVertex(clusterTree, *itF) == componentVertex)
            continue;
        stopfreeKernelVertex = childVertex(clusterTree, *itF);
#if SEQAN_ENABLE_DEBUG
        found = true;
#endif  // #if SEQAN_ENABLE_DEBUG
    }
    SEQAN_ASSERT_TRUE(found);
    // Create new internal tree vertex below stopfree kernel vertex.
    TVertexDescriptor internalTreeVertex = addChild(clusterTree, stopfreeKernelVertex);
    resizeVertexMap(clusterTree, blockDescriptor);
    assignProperty(blockDescriptor, internalTreeVertex, BLOCK_INTERNAL_TREE);

    // Assign root vertex to internal tree vertex.
    appendValue(property(vertexBlocks, v), internalTreeVertex);

    // Stack for DFS search.
    typedef String<TVertexDescriptor> TStack;
    TStack stack;
    appendValue(stack, v);

    // Perform DFS on all internal tree edges to vertices that are not marked
    // with a token yet.
    while (!empty(stack)) {
        v = back(stack);
        eraseBack(stack);
        
        for (TOutEdgeIterator itE(g, v); !atEnd(itE); goNext(itE)) {
            TEdgeDescriptor e = *itE;
            if (!isBitSet(getProperty(edgeFlags, e), EDGE_INTERNAL_TREE))
                continue;  // Skip non-internal-tree edges.
            
            TVertexDescriptor u = getTarget(e);
            if (u == v)
                u = getSource(e);
            SEQAN_ASSERT_TRUE(isBitSet(getProperty(vertexFlags, u), VERTEX_INTERNAL_TREE));
            if (isBitSet(getProperty(vertexFlags, v), VERTEX_TOKEN))
                continue;  // Skip already marked vertices.

            // Add this edge and other end vertex to the internal tree.
            assignProperty(edgeBlock, e, internalTreeVertex);
            appendValue(property(vertexBlocks, u), internalTreeVertex);
            
            appendValue(stack, u);
        }
    }
}

// This is a slightly extended version of Tarjan's classic biconnected
// component search.  The function _decomposeGraphStiegeBiconnect corresponds
// to the routine BICONNECT in (Tarjan, 1973).
template <typename TClusterTree, typename TBlockDescriptorsMap, typename TVertexBlocksMap, typename TEdgeBlockMap, typename TGraphCargo, typename TGraphSpec, typename TComponents, typename TComponentToBlock>
void
_decomposeGraphStiegeFindStopfreeKernels(
        String<__uint32> & vertexFlags,
        String<__uint32> & edgeFlags,
        TClusterTree & clusterTree,
        TBlockDescriptorsMap & blockDescriptor,
        TVertexBlocksMap & vertexBlocks,
        TEdgeBlockMap & edgeBlock, 
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
        _decomposeGraphStiegeBiconnect(vertexFlags, edgeFlags, clusterTree, blockDescriptor, vertexBlocks, edgeBlock, lowPt, number, edgeComponentIds, i, j, stack, v, maxValue<TVertexDescriptor>(), g, components, componentToBlock);
    }

    // Finally, collect internal trees.  The edges have already been marked
    // appropriately.
    for (TVertexIterator it(g); !atEnd(it); ++it) {
        TVertexDescriptor v = *it;
        _decomposeGraphStiegeCollectInternalTrees(vertexFlags, edgeFlags, clusterTree, blockDescriptor, vertexBlocks, edgeBlock, v, g, components, componentToBlock);
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
                     String<String<TBlockMapValue>, TBlockMapSpec> & vertexBlocks,
                     String<TBlockMapValue, TBlockMapSpec> & edgeBlock,
                     Graph<Undirected<TGraphCargo, TGraphSpec> > /*const*/ & g) // TODO(holtgrew): Remove const after fixing iterator for const-graphs
{
    typedef Graph<Undirected<TGraphCargo, TGraphSpec> > TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename VertexDescriptor<TGraph>::Type TGraphVertexDescriptor;
    typedef Graph<Tree<TTreeCargo, TTreeSpec> > TTree;
    typedef typename VertexDescriptor<TTree>::Type TTreeVertexDescriptor;

    // Perform the standard decomposition of the given graph.  We build the
    // cluster tree at the same time.
    
    // Edge and vertex type markers.
    String<unsigned> vertexFlags;
    resizeVertexMap(g, vertexFlags, 0u);
    String<unsigned> edgeFlags;
    resizeEdgeMap(g, edgeFlags, 0u);

    // Initialize decomposition tree.
    //
    // Add root vertex which corresponds to the graph.
    clear(clusterTree);
    TTreeVertexDescriptor r = addVertex(clusterTree);
    assignRoot(clusterTree, r);
    resizeVertexMap(clusterTree, blockDescriptors);
    assignProperty(blockDescriptors, r, BLOCK_GRAPH);
    // Initialize maps from vertices and edges into the decomposition tree.
    resizeVertexMap(g, vertexBlocks);
    resizeEdgeMap(g, edgeBlock);

    // Find connected components.
    String<unsigned> componentIds;
    unsigned componentCount = connectedComponents(g, componentIds);
    String<TTreeVertexDescriptor> componentToBlock;
    resize(componentToBlock, componentCount);
    for (unsigned i = 0; i < componentCount; ++i) {
        TTreeVertexDescriptor x = addChild(clusterTree, r);
        resizeVertexMap(clusterTree, blockDescriptors);
        assignProperty(blockDescriptors, x, BLOCK_PROPER_COMPONENT);  // <-- will update later to improper if isolated
        componentToBlock[i] = x;

        // Each connected component contains exactly one stopfree kernel.
        // This will always be the first child of that building block's tree
        // node.
        TTreeVertexDescriptor y = addChild(clusterTree, x);
        resizeVertexMap(clusterTree, blockDescriptors);
        assignProperty(blockDescriptors, y, BLOCK_STOPFREE_KERNEL);
    }

    // Find isolated vertices and peripheral trees.
    _decomposeGraphStiegeFindPeripheralTrees(vertexFlags, edgeFlags, clusterTree, blockDescriptors, vertexBlocks, edgeBlock, g, componentIds, componentToBlock);

    // Find stopfree kernels and their internal structure.
    _decomposeGraphStiegeFindStopfreeKernels(vertexFlags, edgeFlags, clusterTree, blockDescriptors, vertexBlocks, edgeBlock, g, componentIds, componentToBlock);
}

} // namespace seqan

#endif // #ifndef SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_GRAPH_STIEGE_H_
