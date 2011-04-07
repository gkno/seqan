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
// Implementation of the Digraph decomposition by Stiege.
//
// References:
//
// * Stiege, G. Standard decomposition and periodicity of digraphs. Berichte
//   aus dem Fachbereich Informatik, Universit\"t Oldenburg, November 2001.
// ==========================================================================

// TODO(holtgrew): While all algorithms are linear-time, the constant factor could be improved by doing less steps, doing more at the same time, probably more DFS than iterator.

// TODO(holtgrew): Rather call "standard graph decomposition" and replace the word stiege by standard in file name and symbols.

// TODO(holtgrew): Note that the undirected decomposition is not the same as the one in the 2001 paper. This makes for some impedance mismatch between directed and undirected decomposition.

#ifndef SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_DIGRAPH_STIEGE_H_
#define SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_DIGRAPH_STIEGE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes, Structs, Enums, Tags
// ============================================================================

/**
.Enum.StandardDiDecompositionBlockType
..cat:Graph Decomposition
..summary:Types of building blocks in standard graph decomposition of directed graphs.
..signature:StandardDiDecompositionBlockType
..include:seqan/graph_decomposition.h
..value.BLOCK_DIGRAPH:Building block type digraph (root).
..value.BLOCK_IMPROPER_WEAK_COMPONENT:Building block type improper weak component, singleton.
..value.BLOCK_PROPER_WEAK_COMPONENT:Building block type proper weak component.
..value.BLOCK_ACYCLIC_WEAK_COMPONENT:Building block type acyclic weak component, i.e. underlying graph is connected and digraph itself is dag.
..value.BLOCK_CYCLIC_WEAK_COMPONENT:Building block type cyclick weak component, i.e. underlying graph is connected and digraph contains a cycle.
..value.BLOCK_EXTERNAL_DAG:Building block type external dag.
..value.BLOCK_STRONG_COMPONENT:Building block type strongly connected component.
..value.BLOCK_ACYCLIC_STRONG_COMPONENT:Building block type acyclic strong component, i.e. underlying graph has no cycle.
..value.BLOCK_CYCLIC_STRONG_COMPONENT:Building block type cylic strong component, i.e. underlying graph has cycle.
..value.BLOCK_PERIPHERAL_STREE:Building block type peripheral s-tree.
..value.BLOCK_STOPFREE_SKERNEL:Building block type stopfree s-kernel.
..value.BLOCK_INTERNAL_STREE:Building block type internal s-tree.
..value.BLOCK_SSUBCOMPONENT:Building block type s-subcomponent.
..value.BLOCK_SBIBLOCK:Building block type s-biblock.
 */

enum StandardDiDecompositionBlockType
{
    BLOCK_DIGRAPH,
    BLOCK_IMPROPER_WEAK_COMPONENT,
    BLOCK_PROPER_WEAK_COMPONENT,
    BLOCK_ACYCLIC_WEAK_COMPONENT,
    BLOCK_CYCLIC_WEAK_COMPONENT,
    BLOCK_EXTERNAL_DAG,
    BLOCK_STRONG_COMPONENT,
    BLOCK_ACYCLIC_STRONG_COMPONENT,
    BLOCK_CYCLIC_STRONG_COMPONENT,
    BLOCK_PERIPHERAL_STREE,
    BLOCK_STOPFREE_SKERNEL,
    BLOCK_INTERNAL_STREE,
    BLOCK_SSUBCOMPONENT,
    BLOCK_SBIBLOCK
};

/**
.Enum.StandardDiDecompositionVertexType
..cat:Graph Decomposition
..summary:Types of vertices in standard graph decompositon of directed graphs.
..signature.StandardDiDecompositionVertexType
..include.seqan/graph_decomposition.h
..value.VERTEX_DI_ISOLATED:Vertex is isolated.
..value.VERTEX_DI_ACYCLIC_WEAK_COMPONENT:Vertex is in acyclic weak component.
..value.VERTEX_DI_EXTERNAL_DAG:Vertex is in external dag.
..value.VERTEX_DI_ACYCLIC_STRONG_COMPONENT:Vertex is in acyclic strong component.
..value.VERTEX_DI_PERIPHERAL_STREE:Vertex is in peripheral s-tree.
..value.VERTEX_DI_PERIPHERAL_STREE_ROOT:Vertex is root of peripheral s-tree.
..value.VERTEX_DI_PERIPHERAL_TREE_BORDER:Vertex is border of peripheral s-tree.
..value.VERTEX_DI_SBIBLOCK:Vertex is in s-biblock.
..value.VERTEX_DI_INTERNAL_STREE:Vertex is in internal s-tree.
..value.VERTEX_DI_INTERNAL_STREE_ROOT:Vertex is root of internal s-tree.
 */

enum StandardDiDecompositionVertexType
{
    VERTEX_DI_ISOLATED,
    VERTEX_DI_ACYCLIC_WEAK_COMPONENT,
    VERTEX_DI_EXTERNAL_DAG,
    VERTEX_DI_ACYCLIC_STRONG_COMPONENT,
    VERTEX_DI_PERIPHERAL_STREE,        // non-root of peripheral tree
    VERTEX_DI_PERIPHERAL_STREE_ROOT,   // root of "isolated peripheral tree"
    VERTEX_DI_PERIPHERAL_TREE_BORDER,  // border point/root of non-isolated ptree
    VERTEX_DI_SBIBLOCK,                // in a biblock
    VERTEX_DI_INTERNAL_STREE,          // in internal tree
    VERTEX_DI_INTERNAL_STREE_ROOT      // root of internal tree
};

/**
.Enum.StandardDiDecompositionEdgeType
..cat:Graph Decomposition
..summary:Types of edges in standard graph decompositon of directed graphs.
..signature:StandardDiDecompositionEdgeType
..include:seqan/graph_decomposition.h
..value.EDGE_DI_ACYCLIC_WEAK_COMPONENT:Edge is in acyclic weak component.
..value.EDGE_DI_EXTERNAL_DAG:Edge is in external dag.
..value.EDGE_DI_ACYCLIC_STRONG_COMPONENT:Edge is in acyclic strong component.
..value.EDGE_DI_PERIPHERAL_STREE:Edge is in peripheral s-tree.
..value.EDGE_DI_INTERNAL_STREE:Edge is in internal s-tree.
..value.EDGE_DI_SBIBLOCK:Edge is in s-biblock.
 */

enum StandardDiDecompositionEdgeType
{
    EDGE_DI_ACYCLIC_WEAK_COMPONENT,
    EDGE_DI_EXTERNAL_DAG,
    EDGE_DI_ACYCLIC_STRONG_COMPONENT,
    EDGE_DI_PERIPHERAL_STREE,
    EDGE_DI_INTERNAL_STREE,
    EDGE_DI_SBIBLOCK
};

/**
.Spec.Standard Digraph Decomposition
..summary:Decomposition of a digraph following Stiege's standard decomposition.
..general:Class.GraphDecomposition
..cat:Graph Decomposition
..signature:GraphDecomposition<TGraph, StandardDigraphDecomposition>
..param.TGraph:The type of the digraph to decompose.
...type:Spec.Directed Graph
..include:seqan/graph_decomposition.h

.Memfunc.Standard Digraph Decomposition#GraphDecomposition
..class:Spec.Standard Digraph Decomposition
..summary:Constructor
..signature:GraphDecomposition()
..signature:GraphDecomposition(g)
..param.g:The graph this decomposition is for.
...type:Spec.Directed Graph
 */

// TODO(holtgrew): Maybe use StandardDecomposition tag?
struct StandardDigraphDecomposition_;
typedef Tag<StandardDigraphDecomposition_> StandardDigraphDecomposition;

template <typename TCargo, typename TSpec>
class GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition>
{
public:
    typedef Graph<Directed<TCargo, TSpec> > TGraph_;
    typedef Graph<Tree<> > TDecompositionTree_;
    typedef String<StandardDiDecompositionBlockType> TClusterNodeCargo_;
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
struct ClusterTree<GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> >
{
    typedef GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> TDecomposition_;
    typedef typename TDecomposition_::TDecompositionTree_ Type;
};

template <typename TCargo, typename TSpec>
struct ClusterTree<GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> const>
{
    typedef GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> const TDecomposition_;
    typedef typename TDecomposition_::TDecompositionTree_ const Type;
};

// ----------------------------------------------------------------------------
// Metafunction ClusterNodeCargo
// ----------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
struct ClusterNodeCargo<GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> >
{
    typedef GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> TDecomposition_;
    typedef typename TDecomposition_::TClusterNodeCargo_ Type;
};

template <typename TCargo, typename TSpec>
struct ClusterNodeCargo<GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> const>
{
    typedef GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> const TDecomposition_;
    typedef typename TDecomposition_::TClusterNodeCargo_ const Type;
};

// ----------------------------------------------------------------------------
// Metafunction VertexToClusterMap
// ----------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
struct VertexToClusterMap<GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> >
{
    typedef GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> TDecomposition_;
    typedef typename TDecomposition_::TVertexToClusterMap_ Type;
};

template <typename TCargo, typename TSpec>
struct VertexToClusterMap<GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> const>
{
    typedef GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> const TDecomposition_;
    typedef typename TDecomposition_::TVertexToClusterMap_ const Type;
};

// ----------------------------------------------------------------------------
// Metafunction EdgeToClusterMap
// ----------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
struct EdgeToClusterMap<GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> >
{
    typedef GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> TDecomposition_;
    typedef typename TDecomposition_::TEdgeToClusterMap_ Type;
};

template <typename TCargo, typename TSpec>
struct EdgeToClusterMap<GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> const>
{
    typedef GraphDecomposition<Graph<Directed<TCargo, TSpec> >, StandardDigraphDecomposition> const TDecomposition_;
    typedef typename TDecomposition_::TEdgeToClusterMap_ const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function decomposeGraph()
// ----------------------------------------------------------------------------

template <typename TGraphCargo, typename TGraphSpec>
void
decomposeGraph(GraphDecomposition<Graph<Directed<TGraphCargo, TGraphSpec> >, StandardDigraphDecomposition> & gd,
               Graph<Directed<TGraphCargo, TGraphSpec> > /*const*/ & g) // TODO(holtgrew): Remove const after fixing iterator for const-graphs
{
    typedef Graph<Directed<TGraphCargo, TGraphSpec> > TGraph;
    typedef GraphDecomposition<TGraph, StandardDigraphDecomposition> TGraphDecomposition;
    typedef typename ClusterTree<TGraphDecomposition>::Type TClusterTree;
    typedef typename VertexDescriptor<TClusterTree>::Type TTreeVertexDescriptor;
    typedef typename VertexDescriptor<TGraph>::Type TGraphVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TGraphEdgeDescriptor;
    typedef typename Iterator<TGraph, VertexIterator>::Type TGraphVertexIterator;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TGraphOutEdgeIterator;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TGraphEdgeIterator;
    typedef typename Iterator<String<unsigned>, Rooted>::Type TStringIterator;

    // ------------------------------------------------------------------------
    // Initialization.
    // ------------------------------------------------------------------------

    clear(clusterTree(gd));
    TTreeVertexDescriptor r = addVertex(clusterTree(gd));
    assignRoot(clusterTree(gd), r);
    resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
    assignProperty(clusterNodeCargoMap(gd), r, BLOCK_DIGRAPH);

    // Initialize maps from vertices and edges into the decomposition tree.
    //
    // We use r as a sentinel value, no vertex can be mapped there normally.
    resizeVertexMap(g, vertexToClusterMap(gd), r);
    resizeEdgeMap(g, edgeToClusterMap(gd), r);

    // Graph vertex and edge flags.
    String<__uint32> vertexFlags;
    resizeVertexMap(g, vertexFlags, 0);
    String<__uint32> edgeFlags;
    resizeEdgeMap(g, edgeFlags, 0);

    // ------------------------------------------------------------------------
    // Identify Weak Components, Label Improper Ones
    // ------------------------------------------------------------------------

    // Compute weakly connected components.
    String<unsigned> weakComponentMap;
    unsigned wccCount = weaklyConnectedComponents(g, weakComponentMap);

    // std::cerr << "wccCount == " << wccCount << std::endl;
    // std::cerr << "weak component map:" << std::endl;
    // for (TStringIterator it = begin(weakComponentMap); !atEnd(it); goNext(it)) {
    //     std::cerr << position(it) << " -> " << *it << std::endl;
    // }
    // std::cerr << "--" << std::endl;

    // Compute sizes of weakly connected components.
    String<unsigned> weakComponentSizes;
    resize(weakComponentSizes, wccCount, 0);
    for (TStringIterator it = begin(weakComponentMap); !atEnd(it); goNext(it)) 
        weakComponentSizes[*it] += 1;
    // std::cerr << "component sizes" << std::endl;
    // for (unsigned i = 0; i < length(weakComponentSizes); ++i)
    //     std::cerr << i << " -> " << weakComponentSizes[i] << std::endl;
    // std::cerr << "--" << std::endl;

    // Label improper components.
    unsigned improperWeakCount = 0;
    for (TGraphVertexIterator itV(g); !atEnd(itV); goNext(itV)) {
        if (weakComponentSizes[getProperty(weakComponentMap, *itV)] > 1u)
            continue;
        improperWeakCount += 1;
        // std::cerr << "IMPROPER COMPONENT" << std::endl;
        // Add new cluster.
        TTreeVertexDescriptor x = addChild(clusterTree(gd), r);
        resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
        assignProperty(clusterNodeCargoMap(gd), x, BLOCK_IMPROPER_WEAK_COMPONENT);
        // Map vertex to improper weak component cluster.
        appendValue(property(vertexToClusterMap(gd), *itV), x);
        setBit(property(vertexFlags, *itV), VERTEX_DI_ISOLATED);
    }

    // ------------------------------------------------------------------------
    // Compute Strong Components.
    //
    // This yields the cyclic and acyclic weak components and the external dag.
    // ------------------------------------------------------------------------

    // Compute strongly connected components.
    String<unsigned> strongComponentMap;
    unsigned sccCount = stronglyConnectedComponents(g, strongComponentMap);
    // std::cerr << "sccCount == " << sccCount << std::endl;
    // for (unsigned i = 0; i < length(strongComponentMap); ++i)
    //     std::cerr << i << " --> " << strongComponentMap[i] << std::endl;

    // Compute sizes of strongly connected components.
    String<unsigned> strongComponentSizes;
    resize(strongComponentSizes, sccCount, 0);
    for (TStringIterator it = begin(strongComponentMap); !atEnd(it); goNext(it))
        strongComponentSizes[*it] += 1;

    // Build SCC to WCC map and flags whether a WCC is cyclic.
    String<unsigned> sccToWeakComponent;
    resize(sccToWeakComponent, sccCount, wccCount);  // wccCount is sentinel
    String<bool> isCyclic;
    resize(isCyclic, length(weakComponentMap), false);
    for (TGraphVertexIterator itV(g); !atEnd(itV); goNext(itV)) {
        unsigned strongComponentId = getProperty(strongComponentMap, *itV);
        unsigned weakComponentId = getProperty(weakComponentMap, *itV);
        sccToWeakComponent[strongComponentId] = weakComponentId;
        bool b = strongComponentSizes[strongComponentId] > 1u;
        isCyclic[weakComponentId] = isCyclic[weakComponentId] || b;
        // std::cerr << "v == " << *itV << ", isCyclic[" << weakComponentId << "] == " << isCyclic[weakComponentId] << std::endl;
    }

    // Map vertices and edges from acyclic weak components to clusters.
    String<TTreeVertexDescriptor> weakComponentToClusterMap;
    resize(weakComponentToClusterMap, length(weakComponentMap), r);  // r is sentinel
    for (TGraphVertexIterator itV(g); !atEnd(itV); goNext(itV)) {
        unsigned c = getProperty(weakComponentMap, *itV);
        // Skip improper components and cyclic components.
        if (weakComponentSizes[getProperty(weakComponentMap, *itV)] == 1u)
            continue;
        if (isCyclic[c])
            continue;
        // std::cerr << c << " is not cyclic!" << std::endl;
        // Add a new acyclic cluster.
        TTreeVertexDescriptor x = weakComponentToClusterMap[c];
        // Add a new cluster if this vertex' weak component has none yet.
        if (x == r) {
            TTreeVertexDescriptor y = addChild(clusterTree(gd), r);
            x = addChild(clusterTree(gd), y);
            resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
            assignProperty(clusterNodeCargoMap(gd), y, BLOCK_PROPER_WEAK_COMPONENT);
            assignProperty(clusterNodeCargoMap(gd), x, BLOCK_ACYCLIC_WEAK_COMPONENT);
            weakComponentToClusterMap[c] = x;
        }
        // At this point, x describes the acyclic weak component cluster.
        appendValue(property(vertexToClusterMap(gd), *itV), x);
        // Mark vertex and all outgoing edges.
        TGraphVertexDescriptor v = *itV;
        setBit(property(vertexFlags, v), VERTEX_DI_ACYCLIC_WEAK_COMPONENT);
        for (TGraphOutEdgeIterator itE(g, v); !atEnd(itE); goNext(itE))
            setBit(property(edgeFlags, *itE), EDGE_DI_ACYCLIC_WEAK_COMPONENT);
    }

    // TODO(holtgrew): This loop could probablye be merged with the one above.

    // Maps SCC ids to "proper SCC ids", i.e. the SCCs with more than one
    // vertex.
    unsigned properSccCount = 0;
    String<unsigned> properSccIds;
    resize(properSccIds, sccCount, maxValue<unsigned>());

    // Now, create clusters for cyclic weak components and identify the external dags.
    String<TTreeVertexDescriptor> weakComponentToExternalDagMap;
    resize(weakComponentToExternalDagMap, length(weakComponentMap), r);  // r is sentinel
    for (TGraphVertexIterator itV(g); !atEnd(itV); goNext(itV)) {
        TGraphVertexDescriptor v = *itV;
        // Get id of WCC and skip if it is not cyclic.
        unsigned wccId = getProperty(weakComponentMap, v);
        if (!isCyclic[wccId])
            continue;

        // Get vertex for the WCC.  Add a new cluster if this vertex' weak
        // component has none yet.
        TTreeVertexDescriptor x = weakComponentToClusterMap[wccId];
        if (x == r) {
            TTreeVertexDescriptor y = addChild(clusterTree(gd), r);
            x = addChild(clusterTree(gd), y);
            TTreeVertexDescriptor z = addChild(clusterTree(gd), x);
            resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
            assignProperty(clusterNodeCargoMap(gd), y, BLOCK_PROPER_WEAK_COMPONENT);
            assignProperty(clusterNodeCargoMap(gd), x, BLOCK_CYCLIC_WEAK_COMPONENT);
            assignProperty(clusterNodeCargoMap(gd), z, BLOCK_EXTERNAL_DAG);
            weakComponentToClusterMap[wccId] = x;
            weakComponentToExternalDagMap[wccId] = z;
        }

        // If v is an external dag vertex then classify v and all outgoing edges of v.
        unsigned sccId = getProperty(strongComponentMap, v);
        unsigned sscSize = strongComponentSizes[sccId];
        if (sscSize == 1u) {
            // std::cerr << "v = " << v << " is in external DAG" << std::endl;
            // Mark and classify vertex and edges.
            unsigned c = weakComponentToExternalDagMap[getProperty(weakComponentMap, v)];
            appendValue(property(vertexToClusterMap(gd), v), c);
            setBit(property(vertexFlags, v), VERTEX_DI_EXTERNAL_DAG);
            for (TGraphOutEdgeIterator itE(g, v); !atEnd(itE); goNext(itE)) {
                setBit(property(edgeFlags, *itE), EDGE_DI_EXTERNAL_DAG);
                assignProperty(edgeToClusterMap(gd), *itE, c);
                appendValue(property(vertexToClusterMap(gd), targetVertex(itE)), c);
            }
        } else {
            if (properSccIds[sccId] == maxValue<unsigned>())
                properSccIds[sccId] = properSccCount++;
        }
    }

    // std::cerr << "proper SCC count " << properSccCount << std::endl;

    // TODO(holtgrew): Delete empty external dags.
    
    // ------------------------------------------------------------------------
    // For each SCC:  Decompose Underlying Graph with Standard Decomposition.
    //
    // Then, map the results back to the digraph and emboss the structure upon
    // it.
    // ------------------------------------------------------------------------

    // For each SCC with more than one vertex, we collect the edges
    // (min(u, v), max(u, v)).  From this, we build the subgraph for the new
    // component.  We also collect all vertices in one subgraph.
    typedef Triple<TGraphVertexDescriptor, TGraphVertexDescriptor, TGraphEdgeDescriptor> TTriple;
    String<String<TGraphVertexDescriptor> > subGraphVertices;
    resize(subGraphVertices, properSccCount);
    String<String<TTriple> > subGraphEdges;
    resize(subGraphEdges, properSccCount);
    
    for (TGraphEdgeIterator itE(g); !atEnd(itE); goNext(itE)) {
        TGraphVertexDescriptor u = sourceVertex(itE);
        TGraphVertexDescriptor v = targetVertex(itE);
        if (getProperty(strongComponentMap, u) != getProperty(strongComponentMap, v))
            continue;  // Skip edge between two different SCC.
        if (u > v)
            std::swap(u, v);
        SEQAN_ASSERT_NEQ(u, v);
        unsigned sccId = getProperty(strongComponentMap, u);
        unsigned properSccId = properSccIds[sccId];
        appendValue(subGraphEdges[properSccId], TTriple(u, v, *itE));
        // std::cerr << "PUSH(" << u << ", " << v << ")" << std::endl;
        appendValue(subGraphVertices[properSccId], u);
        appendValue(subGraphVertices[properSccId], v);
    }

    // Build the subgraphs.
    for (unsigned i = 0; i < properSccCount; ++i) {
        // Remove duplicates from edges.
        typedef typename Iterator<String<TGraphVertexDescriptor>, Standard>::Type TIterator;
        TIterator itBegin = begin(subGraphVertices[i], Standard());
        TIterator itEnd = end(subGraphVertices[i], Standard());
        // std::cerr << "SORT()" << std::endl;
        std::sort(itBegin, itEnd);
        TIterator itNewEnd = std::unique(itBegin, itEnd);
        if (itNewEnd != itEnd)
            SEQAN_FAIL("Mapping back edges does not work for multigraphs (yet)!");
        resize(subGraphVertices[i], itNewEnd - itBegin);
        itEnd = itNewEnd;
        // The position of a vertex descriptor in subGraphVertices[i] now
        // gives the id in the new graph.

        // Remove duplicates from edges.
        {
            typedef typename Iterator<String<TTriple>, Standard>::Type TTripleIterator;
            TTripleIterator itEBegin = begin(subGraphEdges[i], Standard());
            TTripleIterator itEEnd = end(subGraphEdges[i], Standard());
            std::sort(itEBegin, itEEnd);
            TTripleIterator itNewEEnd = std::unique(itEBegin, itEEnd);
            resize(subGraphEdges[i], itNewEEnd - itEBegin);
        }

        typedef Graph<Undirected<> > TSubGraph;
        typedef typename VertexDescriptor<TSubGraph>::Type TSubgraphVertexDescriptor;
        typedef typename EdgeDescriptor<TSubGraph>::Type TSubgraphEdgeDescriptor;
        TSubGraph subGraph;

        // Add vertices.
        for (unsigned j = 0; j < length(subGraphVertices[i]); ++j) {
            TSubgraphVertexDescriptor x = addVertex(subGraph);
            (void)x;  // In case we run without assertions.
            SEQAN_ASSERT_EQ(x, j);
        }

        // std::cerr << "subgraph" << std::endl;
        // for (unsigned j = 0; j < length(subGraphVertices[i]); ++j) {
        //     std::cerr << j << " --> " << subGraphVertices[i][j] << std::endl;
        // }
        // std::cerr << std::endl;
        // for (unsigned j = 0; j < length(subGraphEdges[i]); ++j) {
        //     std::cerr << subGraphEdges[i][j] << std::endl;
        // }

        // Add the edges to the undirected graph.
        String<TGraphEdgeDescriptor> subgraphEdgeToGraphEdges;
        reserve(subgraphEdgeToGraphEdges, length(subGraphEdges));
        typedef typename Iterator<String<TTriple>, Rooted>::Type TTripleIterator;
        for (TTripleIterator itP = begin(subGraphEdges[i]); !atEnd(itP); goNext(itP)) {
            // Get descriptors of vertices incient to the edge of subgraph by binary search.
            TGraphVertexDescriptor u = value(itP).i1;
            TIterator itU = std::lower_bound(itBegin, itEnd, u);
            SEQAN_ASSERT(itU != itEnd);
            SEQAN_ASSERT_EQ(*itU, u);
            TSubgraphVertexDescriptor us = itU - itBegin;
            TGraphVertexDescriptor v = value(itP).i2;
            TIterator itV = std::lower_bound(itBegin, itEnd, v);
            SEQAN_ASSERT(itV != itEnd);
            SEQAN_ASSERT_EQ(*itV, v);
            TSubgraphVertexDescriptor vs = itV - itBegin;
            SEQAN_ASSERT_NEQ(us, vs);
            // Add edge to subgraph.
            TSubgraphEdgeDescriptor e = addEdge(subGraph, us, vs);
            resizeEdgeMap(subGraph, subgraphEdgeToGraphEdges);
            // assignProperty(subgraphEdgeToGraphEdges, e, value(itP).i3);
            assignProperty(subgraphEdgeToGraphEdges, e, TGraphEdgeDescriptor());
        }
        // write(std::cerr, subGraph, DotDrawing());

        // Perform standard decomposition of the undirected graph.
        typedef GraphDecomposition<TSubGraph, StandardDecomposition> TSubGraphDecomposition;
        TSubGraphDecomposition subGraphDecomposition(subGraph);
        // std::cerr << "STANDARD DECOMPOSITION" << std::endl;
        decomposeGraph(subGraphDecomposition, subGraph);

        write(std::cout, subGraph, DotDrawing());
        writeDecompositionTree(std::cout, subGraphDecomposition);

        // ------------------------------------------------------------------------
        // Map back the graph standard decompositon decomposition to the
        // digraph and merge with digraph standard decomposition.
        // ------------------------------------------------------------------------

        // Add node to cluster tree for strong component.
        unsigned c = getProperty(weakComponentMap, front(subGraphVertices[i]));
        TTreeVertexDescriptor wcVertex = weakComponentToClusterMap[c];
        SEQAN_ASSERT_EQ(getProperty(clusterNodeCargoMap(gd), wcVertex), BLOCK_CYCLIC_WEAK_COMPONENT);
        TTreeVertexDescriptor sccVertex = addChild(clusterTree(gd), wcVertex);
        resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
        assignProperty(clusterNodeCargoMap(gd), sccVertex, BLOCK_STRONG_COMPONENT);

        // Traverse the cluster tree of the graph decomposition, copy it into
        // the corresponding structure in the digraph cluster tree and also
        // assign the vertices appropriately.
        typedef typename ClusterTree<TSubGraphDecomposition>::Type TSubGraphClusterTree;
        typedef typename VertexDescriptor<TSubGraphClusterTree>::Type TSubGraphTreeVertexDescriptor;
        typedef typename Iterator<TSubGraphClusterTree, DfsPreorder>::Type TDfsIterator;// TODO(holtgrew): Rename to DfsPreorderIterator to make things more consistent?
        typedef typename Iterator<TSubGraphClusterTree, OutEdgeIterator>::Type TOutEdgeIterator;
        typedef typename VertexDescriptor<TSubGraph>::Type TSubGraphVertexDescriptor;
        String<TTreeVertexDescriptor> subGraphTreeVertexToTreeVertex;
        resizeVertexMap(clusterTree(subGraphDecomposition), subGraphTreeVertexToTreeVertex, r);  // r is sentinel value
        assignProperty(subGraphTreeVertexToTreeVertex, root(clusterTree(subGraphDecomposition)), sccVertex);

        // We need to decide whether we have a cyclic or an acyclic SCC at
        // this point.  We can look at the decomposition tree of the
        // underlying graph.  If it only consist of an external tree then we
        // have an acyclic SCC.
        bool isAcyclic = false;
        TSubGraphTreeVertexDescriptor peripheralTreeVertex = root(clusterTree(subGraphDecomposition));  // sentinel
        if (numChildren(clusterTree(subGraphDecomposition), root(clusterTree(subGraphDecomposition))) == 1u) {
            // First child has to be proper connected component, since it is a
            // proper SCC.  Look into it.  If the only child is an external
            // tree then we are done.
            TOutEdgeIterator itE(clusterTree(subGraphDecomposition), root(clusterTree(subGraphDecomposition)));
            if (numChildren(clusterTree(subGraphDecomposition), targetVertex(itE)) == 1u) {
                TOutEdgeIterator itF(clusterTree(subGraphDecomposition), targetVertex(itE));
                if (getProperty(clusterNodeCargoMap(subGraphDecomposition), targetVertex(itF)) == BLOCK_PERIPHERAL_TREE) {
                    isAcyclic = true;
                    peripheralTreeVertex = targetVertex(itF);
                }
            }
        }

        if (isAcyclic) {
            // In case of acyclic components we are done after inserting one
            // vertex and updating the maps.
            TTreeVertexDescriptor aSccVertex = addChild(clusterTree(gd), sccVertex);
            assignProperty(subGraphTreeVertexToTreeVertex, peripheralTreeVertex, aSccVertex);
            assignProperty(clusterNodeCargoMap(gd), aSccVertex, BLOCK_ACYCLIC_STRONG_COMPONENT);
            // Update map for vertices.
            typedef typename Iterator<TSubGraph, VertexIterator>::Type TVertexIterator;
            typedef typename Iterator<TSubGraph, EdgeIterator>::Type TEdgeIterator;
            for (TVertexIterator itU(subGraph); !atEnd(itU); goNext(itU)) {
                TSubGraphVertexDescriptor u = *itU;
                TGraphVertexDescriptor uG = getProperty(subGraphVertices[i], u);
                appendValue(property(vertexToClusterMap(gd), uG), aSccVertex);
            }
            // Update map for edges.
            for (TEdgeIterator itE(subGraph); !atEnd(itE); goNext(itE))
                assignProperty(edgeToClusterMap(gd), *itE, aSccVertex);

            continue;
        }

        // If we reach here, we are not in an acyclic weak component but in a
        // cyclic weak component!
        //
        // Insert vertex into cluster tree for the cyclic strongly connected
        // component.
        typename Iterator<TSubGraphClusterTree, AdjacencyIterator>::Type adjIt(clusterTree(subGraphDecomposition), root(clusterTree(subGraphDecomposition)));
        TSubGraphTreeVertexDescriptor properComponentVertex = *adjIt;

        TTreeVertexDescriptor cSccVertex = addChild(clusterTree(gd), sccVertex);
        resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
        assignProperty(clusterNodeCargoMap(gd), cSccVertex, BLOCK_CYCLIC_STRONG_COMPONENT);
        assignProperty(subGraphTreeVertexToTreeVertex, properComponentVertex, cSccVertex);
        std::cerr << "sgtvttv[" << peripheralTreeVertex << "] := " << cSccVertex << std::endl;

        // Copy over the cluster tree from the undirected subgraph.
        // TDfsIterator itV(clusterTree(subGraphDecomposition), root(clusterTree(subGraphDecomposition)));
        for (TDfsIterator itV(clusterTree(subGraphDecomposition), properComponentVertex); !atEnd(itV); goNext(itV)) {
            // We are in/below the cyclic connected component vertex now.
            // std::cerr << "Working on " << *itV << std::endl;

            TSubGraphVertexDescriptor v = *itV;
            TTreeVertexDescriptor vG = getProperty(subGraphTreeVertexToTreeVertex, v);  // v in CT of original graph
            // std::cerr << "sgtvttv[" << v << "] == " << vG << std::endl;
            SEQAN_ASSERT_NEQ(vG, r);  // r was sentinel
            for (TOutEdgeIterator itE(clusterTree(subGraphDecomposition), v); !atEnd(itE); goNext(itE)) {  // edge (v, w)
                TSubGraphVertexDescriptor w = targetVertex(itE);
                // Insert correspondance of w into cluster tree for G, is wG,
                // assign map entry.
                TTreeVertexDescriptor wG = addChild(clusterTree(gd), vG);
                assignProperty(subGraphTreeVertexToTreeVertex, w, wG);
                // std::cerr << "sgtvttv[" << w << "] := " << wG << std::endl;
                resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
                // Set properties of vertex in new tree.
                StandardDecompositionBlockType blockType = getProperty(clusterNodeCargoMap(subGraphDecomposition), w);
                switch (blockType) {
                    case BLOCK_PERIPHERAL_TREE:
                        assignProperty(clusterNodeCargoMap(gd), wG, BLOCK_PERIPHERAL_STREE);
                        break;
                    case BLOCK_STOPFREE_KERNEL:
                        assignProperty(clusterNodeCargoMap(gd), wG, BLOCK_STOPFREE_SKERNEL);
                        break;
                    case BLOCK_INTERNAL_TREE:
                        assignProperty(clusterNodeCargoMap(gd), wG, BLOCK_INTERNAL_STREE);
                        break;
                    case BLOCK_SUBCOMPONENT:
                        assignProperty(clusterNodeCargoMap(gd), wG, BLOCK_SSUBCOMPONENT);
                        break;
                    case BLOCK_BIBLOCK:
                        assignProperty(clusterNodeCargoMap(gd), wG, BLOCK_SBIBLOCK);
                        break;
                    default:
                        SEQAN_FAIL("Wrong undirected graph block type (%d)!", int(blockType));
                }
            }
        }

        // Copy over the cluster-tree mapping for vertices in subgraph to those
        // in the original graph.
        typedef typename Iterator<TSubGraph, VertexIterator>::Type TVertexIterator;
        typedef typename Iterator<TSubGraph, EdgeIterator>::Type TEdgeIterator;
        for (TVertexIterator itU(subGraph); !atEnd(itU); goNext(itU)) {
            TSubGraphVertexDescriptor u = *itU;
            TGraphVertexDescriptor uG = getProperty(subGraphVertices[i], u);
            typedef typename VertexToClusterMap<TSubGraphDecomposition>::Type TVertexToClusterMap;
            typedef typename Value<TVertexToClusterMap>::Type TTreeVertices;
            typedef typename Iterator<TTreeVertices, Standard>::Type TIterator;
            TTreeVertices & map = property(vertexToClusterMap(subGraphDecomposition), u);
            for (TIterator it = begin(map); it != end(map); ++it)
                appendValue(property(vertexToClusterMap(gd), uG), getProperty(subGraphTreeVertexToTreeVertex, *it));
        }

        // Copy over the cluster-tree mapping for edges in subgraph to those
        // in the original graph.
        for (TEdgeIterator itE(subGraph); !atEnd(itE); goNext(itE)) {
            TGraphEdgeDescriptor e = getProperty(subgraphEdgeToGraphEdges, *itE);
            TSubGraphTreeVertexDescriptor v = getProperty(edgeToClusterMap(subGraphDecomposition), *itE);
            TTreeVertexDescriptor vG = getProperty(subGraphTreeVertexToTreeVertex, v);
            assignProperty(edgeToClusterMap(gd), e, vG);
        }
 

        /*
            for (TOutEdgeIterator itE(clusterTree(subGraphDecomposition), *itV); !atEnd(itE); goNext(itE)) {  // edge (v, w)
                // Add child from subgraph tree to graph tree.
                TSubGraphVertexDescriptor w = targetVertex(itE);
                SEQAN_ASSERT_EQ(getProperty(subGraphTreeVertexToTreeVertex, w), r);
                TTreeVertexDescriptor x = addChild(clusterTree(gd), getProperty(subGraphTreeVertexToTreeVertex, *itV));
                resizeVertexMap(clusterTree(gd), clusterNodeCargoMap(gd));
                // Map from subgraph cluster tree into graph cluster tree to
                // new vertex.
                assignProperty(subGraphTreeVertexToTreeVertex, w, x);
            }
        */

        // TODO(holtgrew): Edges are not properly marked and mapped yet!
    }

    // Collect vertices of cyclic weak components
    // for (TStringIterator it = begin(weakComponentMap); !atEnd(it); goNext(it)) {
    //     if (!isCyclic[c])
    //         continue;
    // }
}

// ----------------------------------------------------------------------------
// Function writeDecompositionTree();  For stream output.
// ----------------------------------------------------------------------------

template <typename TStream, typename TGraphCargo, typename TGraphSpec>
void writeDecompositionTree(
        TStream & stream,
        GraphDecomposition<Graph<Directed<TGraphCargo, TGraphSpec> >, StandardDigraphDecomposition> & gd)
{
    typedef GraphDecomposition<Graph<Directed<TGraphCargo, TGraphSpec> >, StandardDigraphDecomposition> TGraphDecomposition;
    typedef Graph<Directed<TGraphCargo, TGraphSpec> > TGraph;
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
    stream << "/* Decomposition tree for digraph */" << std::endl << std::endl;
    stream << "digraph G {" << std::endl;

    stream << std::endl;
    stream << "/* Nodes */" << std::endl;
    stream << std::endl;

    static char const * blockTypeNames[] = {
        "digraph", "improper weak component", "proper weak component",
        "acyclic weak component", "cyclic weak component", "external dag",
        "strong components", "acyclic strong component", "cyclic strong component",
        "peripheral s-tree", "stopfree s-kernel", "internal s-tree",
        "s-subcomponent", "s-biblock"
    };

    typedef typename Iterator<TTree, VertexIterator>::Type TTreeVertexIterator;
    for (TTreeVertexIterator itV(clusterTree(gd)); !atEnd(itV); goNext(itV)) {
        CharString label = blockTypeNames[getProperty(clusterNodeCargoMap(gd), *itV)];
        if (getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_SBIBLOCK ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_INTERNAL_STREE ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_PERIPHERAL_STREE ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_ACYCLIC_STRONG_COMPONENT ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_EXTERNAL_DAG ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_ACYCLIC_WEAK_COMPONENT ||
            getProperty(clusterNodeCargoMap(gd), *itV) == BLOCK_IMPROPER_WEAK_COMPONENT) {
            append(label, "\\n(");
            typedef typename Iterator<String<TTreeVertexDescriptor>, Standard>::Type TIter;
            TIter itBegin = begin(property(blocksForVertex, *itV));
            TIter itEnd = end(property(blocksForVertex, *itV));
            std::sort(itBegin, itEnd);
            TIter itNewEnd = std::unique(itBegin, itEnd);
            for (TIter it = itBegin; it != itNewEnd; ++it) {
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

#endif // #ifndef SEQAN_HEADER_GRAPH_DECOMPOSITON_GRAPH_DECOMPOSE_DIGRAPH_STIEGE_H_
