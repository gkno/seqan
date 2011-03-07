// ==========================================================================
//                             SeqAn app_template                            
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

#ifndef SEQAN_GRAPH_DECOMPOSITION_GRAPH_DECOMPOSITION_BASE_H_
#define SEQAN_GRAPH_DECOMPOSITION_GRAPH_DECOMPOSITION_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.GraphDecomposition
..summary:Decomposition of a graph.
..cat:Graph Decomposition
..signature:GraphDecomposition<TGraph, TSpec>
..param.TGraph:The type of the graph to decompose.
...type:Class.Graph
..param.TSpec:Specialization tag.
..include:seqan/graph_decomposition.h
 */

template <typename TGraph, typename TSpec>
class GraphDecomposition;

// ============================================================================
// Metafunctions
// ============================================================================

// TODO(holtgrew): Documentation for these metafunctions.

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TGraph, typename TSpec>
struct Host<GraphDecomposition<TGraph, TSpec> >
{
    typedef TGraph Type;
};

template <typename TGraph, typename TSpec>
struct Host<GraphDecomposition<TGraph, TSpec> const>
{
    typedef TGraph const Type;
};

// ----------------------------------------------------------------------------
// Metafunction ClusterTree
// ----------------------------------------------------------------------------

template <typename T>
struct ClusterTree;

// ----------------------------------------------------------------------------
// Metafunction ClusterNodeCargo
// ----------------------------------------------------------------------------

// TODO(holtgrew): RenameTo ClusterNodeCargoMap.

template <typename T>
struct ClusterNodeCargo;

// ----------------------------------------------------------------------------
// Metafunction VertexToClusterMap
// ----------------------------------------------------------------------------

template <typename T>
struct VertexToClusterMap;

// ----------------------------------------------------------------------------
// Metafunction EdgeToClusterMap
// ----------------------------------------------------------------------------

template <typename T>
struct EdgeToClusterMap;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function decomposeGraph()
// ----------------------------------------------------------------------------

/**
.Function.decomposeGraph
..cat:Graph Decomposition
..summary:Hierarchial decomposition / clustering of graphs.
..signature:decomposeGraph(graphDecomposition, g)
..param.graphDecomposition:Target object to store the decomposition / clustering of the graph. The algorithm is used by the type put here.
...type:Class.GraphDecomposition
..param.g:The graph to decompose.
...type:Class.Graph
..remarks:The type of $gd$ has to fit the type of $graphDecomposition$.
..example.code:
typedef Graph<Undirected<> > TGraph;
typedef GraphDecomposition<TGraph, StandardDecomposition> TGraphDecomposition;

TGraph g;
// ... construct graph

// Now, perform the graph decomposition.
TGraphDecomposition gd(g);
decomposeGraph(gd, g);
 */

// ----------------------------------------------------------------------------
// Function _dataHost()
// ----------------------------------------------------------------------------

///.Function.host.param.object.type:Class.GraphDecomposition

template <typename TGraph, typename TSpec>
inline
Holder<typename Host<GraphDecomposition<TGraph, TSpec> >::Type> &
_dataHost(GraphDecomposition<TGraph, TSpec> & graphDecomposition)
{
    return graphDecomposition.data_host;
}

template <typename TGraph, typename TSpec>
inline
Holder<typename Host<GraphDecomposition<TGraph, TSpec> const>::Type> &
_dataHost(GraphDecomposition<TGraph, TSpec> const & graphDecomposition)
{
    return graphDecomposition.data_host;
}

// ----------------------------------------------------------------------------
// Function clusterTree()
// ----------------------------------------------------------------------------

/**
.Function.clusterTree
..cat:Graph Decomposition
..summary:Return @Spec.Tree@ object used for the graph decomposition structure.
..signature:clusterTree(graphDecomposition)
..param.graphDecomposition:The @Class.GraphDecomposition@ object to query.
...type:Class.GraphDecomposition
..returns:Reference @Spec.Tree@ object that is used for storing the structure of the graph decomposition.
..see:Function.clusterNodeCargoMap
 */

template <typename TGraph, typename TSpec>
inline
typename ClusterTree<GraphDecomposition<TGraph, TSpec> >::Type &
clusterTree(GraphDecomposition<TGraph, TSpec> & graphDecomposition)
{
    return graphDecomposition._clusterTree;
}

template <typename TGraph, typename TSpec>
inline
typename ClusterTree<GraphDecomposition<TGraph, TSpec> const>::Type &
clusterTree(GraphDecomposition<TGraph, TSpec> const & graphDecomposition)
{
    return graphDecomposition._clusterTree;
}

// ----------------------------------------------------------------------------
// Function clusterNodeCargoMap()
// ----------------------------------------------------------------------------

/**
.Function.clusterNodeCargoMap
..cat:Graph Decomposition
..summary:Returns the external vertex map used for labeling the decomposition tree's nodes.
..signature:clusterNodeCargoMap(graphDecomposition)
..param.graphDecomposition:The @Class.GraphDecomposition@ object to query.
...type:Class.GraphDecomposition
..returns:Reference to external map for labeling the cluster tree's nodes.
..see:Function.clusterTree
 */

template <typename TGraph, typename TSpec>
inline
typename ClusterNodeCargo<GraphDecomposition<TGraph, TSpec> >::Type &
clusterNodeCargoMap(GraphDecomposition<TGraph, TSpec> & graphDecomposition)
{
    return graphDecomposition._clusterNodeCargoMap;
}

template <typename TGraph, typename TSpec>
inline
typename ClusterNodeCargo<GraphDecomposition<TGraph, TSpec> const>::Type &
clusterNodeCargoMap(GraphDecomposition<TGraph, TSpec> const & graphDecomposition)
{
    return graphDecomposition._clusterNodeCargoMap;
}

// ----------------------------------------------------------------------------
// Function vertexToClusterMap()
// ----------------------------------------------------------------------------

/**
.Function.vertexToClusterMap
..cat:Graph Decomposition
..summary:Returns the external vertex map used for assigning graph vertices to cluster tree node(s).
..signature:clusterNodeCargoMap(graphDecomposition)
..param.graphDecomposition:The @Class.GraphDecomposition@ object to query.
...type:Class.GraphDecomposition
..returns:Reference to external map that assign vertex identifiers of the cluster tree node(s) to graph vertices.
..see:Function.edgeToClusterMap
 */

template <typename TGraph, typename TSpec>
inline
typename VertexToClusterMap<GraphDecomposition<TGraph, TSpec> >::Type &
vertexToClusterMap(GraphDecomposition<TGraph, TSpec> & graphDecomposition)
{
    return graphDecomposition._vertexToClusterMap;
}

template <typename TGraph, typename TSpec>
inline
typename VertexToClusterMap<GraphDecomposition<TGraph, TSpec> const>::Type &
vertexToClusterMap(GraphDecomposition<TGraph, TSpec> const & graphDecomposition)
{
    return graphDecomposition._vertexToClusterMap;
}

// ----------------------------------------------------------------------------
// Function edgeToClusterMap()
// ----------------------------------------------------------------------------

/**
.Function.edgeToClusterMap
..cat:Graph Decomposition
..summary:Returns the external vertex map used for assigning graph edges to cluster tree node(s).
..signature:edgeToClusterMap(graphDecomposition)
..param.graphDecomposition:The @Class.GraphDecomposition@ object to query.
...type:Class.GraphDecomposition
..returns:Reference to external map that assign vertex identifiers of the cluster tree node(s) to graph edges.
..see:Function.vertexToClusterMap
 */

template <typename TGraph, typename TSpec>
inline
typename EdgeToClusterMap<GraphDecomposition<TGraph, TSpec> >::Type &
edgeToClusterMap(GraphDecomposition<TGraph, TSpec> & graphDecomposition)
{
    return graphDecomposition._edgeToClusterMap;
}

template <typename TGraph, typename TSpec>
inline
typename EdgeToClusterMap<GraphDecomposition<TGraph, TSpec> const>::Type &
edgeToClusterMap(GraphDecomposition<TGraph, TSpec> const & graphDecomposition)
{
    return graphDecomposition._edgeToClusterMap;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_GRAPH_DECOMPOSITION_GRAPH_DECOMPOSITION_BASE_H_
