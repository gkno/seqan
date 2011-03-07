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
// Tests for the graph decomposition algorithms.
// ==========================================================================

#ifndef TEST_GRAPH_DECOMPOSITION_DECOMPOSITION_H_
#define TEST_GRAPH_DECOMPOSITION_DECOMPOSITION_H_

SEQAN_DEFINE_TEST(test_graph_decomposition_digraph_condensed)
{
    // TODO(holtgrew): Write me!
}

SEQAN_DEFINE_TEST(test_graph_decomposition_digraph_stiege)
{
    // TODO(holtgrew): Write me!
}

SEQAN_DEFINE_TEST(test_graph_decomposition_graph_stiege)
{
    using namespace seqan;

    typedef Graph<Undirected<> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TGraphVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TGraphEdgeDescriptor;

    // ------------------------------------------------------------------------
    // Build Test Data
    // ------------------------------------------------------------------------

    // Build the same graph as in Stiege's example in the 1996 paper, Figure 1.
    TGraph graph;
    TGraphVertexDescriptor a = addVertex(graph);
    TGraphVertexDescriptor b = addVertex(graph);
    TGraphVertexDescriptor c1 = addVertex(graph);
    TGraphVertexDescriptor c2 = addVertex(graph);
    TGraphVertexDescriptor c3 = addVertex(graph);
    TGraphVertexDescriptor c4 = addVertex(graph);
    TGraphVertexDescriptor d1 = addVertex(graph);
    TGraphVertexDescriptor d2 = addVertex(graph);
    TGraphVertexDescriptor d3 = addVertex(graph);
    TGraphVertexDescriptor d4 = addVertex(graph);
    TGraphVertexDescriptor e1 = addVertex(graph);
    TGraphVertexDescriptor e2 = addVertex(graph);
    TGraphVertexDescriptor e3 = addVertex(graph);
    TGraphVertexDescriptor e4 = addVertex(graph);
    TGraphVertexDescriptor f1 = addVertex(graph);
    TGraphVertexDescriptor f2 = addVertex(graph);
    TGraphVertexDescriptor f3 = addVertex(graph);
    TGraphVertexDescriptor f4 = addVertex(graph);
    TGraphVertexDescriptor f5 = addVertex(graph);
    TGraphVertexDescriptor f6 = addVertex(graph);
    TGraphVertexDescriptor f7 = addVertex(graph);
    TGraphVertexDescriptor g = addVertex(graph);
    TGraphVertexDescriptor h1 = addVertex(graph);
    TGraphVertexDescriptor h2 = addVertex(graph);
    TGraphVertexDescriptor h3 = addVertex(graph);
    TGraphVertexDescriptor h4 = addVertex(graph);
    TGraphVertexDescriptor i = addVertex(graph);
    TGraphVertexDescriptor j1 = addVertex(graph);
    TGraphVertexDescriptor j2 = addVertex(graph);
    TGraphVertexDescriptor j3 = addVertex(graph);
    TGraphVertexDescriptor j4 = addVertex(graph);

    TGraphEdgeDescriptor a_c1 = addEdge(graph, a, c1);
    TGraphEdgeDescriptor b_c1 = addEdge(graph, b, c1);
    TGraphEdgeDescriptor c1_c2 = addEdge(graph, c1, c2);
    TGraphEdgeDescriptor c1_c4 = addEdge(graph, c1, c4);
    TGraphEdgeDescriptor c2_c3 = addEdge(graph, c2, c3);
    TGraphEdgeDescriptor c3_c4 = addEdge(graph, c3, c4);
    TGraphEdgeDescriptor c3_d1 = addEdge(graph, c3, d1);
    TGraphEdgeDescriptor d1_d2 = addEdge(graph, d1, d2);
    TGraphEdgeDescriptor d2_d3 = addEdge(graph, d2, d3);
    TGraphEdgeDescriptor d2_d4 = addEdge(graph, d2, d4);
    TGraphEdgeDescriptor d2_h1 = addEdge(graph, d2, h1);
    TGraphEdgeDescriptor d3_f1 = addEdge(graph, d3, f1);
    TGraphEdgeDescriptor d4_e1 = addEdge(graph, d4, e1);
    TGraphEdgeDescriptor e1_e2 = addEdge(graph, e1, e2);
    TGraphEdgeDescriptor e1_e3 = addEdge(graph, e1, e3);
    TGraphEdgeDescriptor e2_e4 = addEdge(graph, e2, e4);
    TGraphEdgeDescriptor e3_e4 = addEdge(graph, e3, e4);
    TGraphEdgeDescriptor f1_f2 = addEdge(graph, f1, f2);
    TGraphEdgeDescriptor f1_f4 = addEdge(graph, f1, f4);
    TGraphEdgeDescriptor f1_g = addEdge(graph, f1, g);
    TGraphEdgeDescriptor f2_f3 = addEdge(graph, f2, f3);
    TGraphEdgeDescriptor f3_f4 = addEdge(graph, f3, f4);
    TGraphEdgeDescriptor f1_f5 = addEdge(graph, f1, f5);
    TGraphEdgeDescriptor f1_f7 = addEdge(graph, f1, f7);
    TGraphEdgeDescriptor f5_f6 = addEdge(graph, f5, f6);
    TGraphEdgeDescriptor f6_f7 = addEdge(graph, f6, f7);
    TGraphEdgeDescriptor h1_h2 = addEdge(graph, h1, h2);
    TGraphEdgeDescriptor h1_h3 = addEdge(graph, h1, h3);
    TGraphEdgeDescriptor h1_h4 = addEdge(graph, h1, h4);
    TGraphEdgeDescriptor j1_j2 = addEdge(graph, j1, j2);
    TGraphEdgeDescriptor j2_j3 = addEdge(graph, j2, j3);
    TGraphEdgeDescriptor j3_j4 = addEdge(graph, j3, j4);

    // ------------------------------------------------------------------------
    // Call Function To Test
    // ------------------------------------------------------------------------

    typedef GraphDecomposition<TGraph, StandardDecomposition> TGraphDecomposition;

    TGraphDecomposition gd(graph);
    decomposeGraph(gd, graph);
    
    // write(std::cout, graph, DotDrawing());

    // ------------------------------------------------------------------------
    // Check Result
    // ------------------------------------------------------------------------

    // Check clusterTree(gd) and clusterNodeCargoMap(gd) first.
    SEQAN_ASSERT_EQ(length(clusterNodeCargoMap(gd)) - 2, numVertices(clusterTree(gd)));  // Two components without stopfree kernel.

    typedef typename ClusterTree<TGraphDecomposition>::Type TClusterTree;
    typedef typename Iterator<TClusterTree, AdjacencyIterator>::Type TAdjacencyIterator;
    typedef typename VertexDescriptor<TClusterTree>::Type TVertexDescriptor;

    TVertexDescriptor r = root(clusterTree(gd));
    TAdjacencyIterator itR(clusterTree(gd), r);
    SEQAN_ASSERT_NOT(atEnd(itR));
    TVertexDescriptor pc1 = *itR;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), pc1), BLOCK_PROPER_COMPONENT);
    // { pc1
    TAdjacencyIterator itPc1(clusterTree(gd), pc1);
    SEQAN_ASSERT_NOT(atEnd(itPc1));
    TVertexDescriptor pt = *itPc1;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), pt), BLOCK_PERIPHERAL_TREE);
    SEQAN_ASSERT_EQ(numChildren(clusterTree(gd), pt), 0u);
    goNext(itPc1);
    SEQAN_ASSERT(atEnd(itPc1));
    // } pc1
    goNext(itR);
    SEQAN_ASSERT_NOT(atEnd(itR));
    TVertexDescriptor ip1 = *itR;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), ip1), BLOCK_IMPROPER_COMPONENT);
    SEQAN_ASSERT_EQ(numChildren(clusterTree(gd), ip1), 0u);
    goNext(itR);
    SEQAN_ASSERT_NOT(atEnd(itR));
    TVertexDescriptor pc2 = *itR;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), pc2), BLOCK_PROPER_COMPONENT);
    // { pc2
    TAdjacencyIterator itPc2(clusterTree(gd), pc2);
    SEQAN_ASSERT_NOT(atEnd(itPc2));
    TVertexDescriptor pt1 = *itPc2;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), pt1), BLOCK_PERIPHERAL_TREE);
    goNext(itPc2);
    SEQAN_ASSERT_NOT(atEnd(itPc2));
    TVertexDescriptor pt2 = *itPc2;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), pt2), BLOCK_PERIPHERAL_TREE);
    goNext(itPc2);
    SEQAN_ASSERT_NOT(atEnd(itPc2));
    TVertexDescriptor pt3 = *itPc2;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), pt3), BLOCK_PERIPHERAL_TREE);
    goNext(itPc2);
    SEQAN_ASSERT_NOT(atEnd(itPc2));
    TVertexDescriptor sk = *itPc2;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), sk), BLOCK_STOPFREE_KERNEL);
    // { pc2 -> sk
    TAdjacencyIterator itSk(clusterTree(gd), sk);
    SEQAN_ASSERT_NOT(atEnd(itSk));
    TVertexDescriptor it = *itSk;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), it), BLOCK_INTERNAL_TREE);
    SEQAN_ASSERT_EQ(numChildren(clusterTree(gd), it), 0u);
    goNext(itSk);
    SEQAN_ASSERT_NOT(atEnd(itSk));
    TVertexDescriptor sc1 = *itSk;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), sc1), BLOCK_SUBCOMPONENT);
    // { pc2 -> sk -> sc1
    TAdjacencyIterator itSc1(clusterTree(gd), sc1);
    SEQAN_ASSERT_NOT(atEnd(itSc1));
    TVertexDescriptor bb1 = *itSc1;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), bb1), BLOCK_BIBLOCK);
    SEQAN_ASSERT_EQ(numChildren(clusterTree(gd), bb1), 0u);
    goNext(itSc1);
    SEQAN_ASSERT(atEnd(itSc1));
    // } pc2 -> sk -> sc1
    goNext(itSk);
    SEQAN_ASSERT_NOT(atEnd(itSk));
    TVertexDescriptor sc2 = *itSk;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), sc2), BLOCK_SUBCOMPONENT);
    // { pc2 -> sk -> sc2
    TAdjacencyIterator itSc2(clusterTree(gd), sc2);
    SEQAN_ASSERT_NOT(atEnd(itSc2));
    TVertexDescriptor bb2 = *itSc2;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), bb2), BLOCK_BIBLOCK);
    SEQAN_ASSERT_EQ(numChildren(clusterTree(gd), bb2), 0u);
    goNext(itSc2);
    SEQAN_ASSERT_NOT(atEnd(itSc2));
    TVertexDescriptor bb3 = *itSc2;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), bb3), BLOCK_BIBLOCK);
    SEQAN_ASSERT_EQ(numChildren(clusterTree(gd), bb3), 0u);
    goNext(itSc2);
    SEQAN_ASSERT(atEnd(itSc2));
    // } pc2 -> sk -> sc2
    goNext(itSk);
    SEQAN_ASSERT_NOT(atEnd(itSk));
    TVertexDescriptor sc3 = *itSk;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), sc3), BLOCK_SUBCOMPONENT);
    // { pc2 -> sk -> sc3
    TAdjacencyIterator itSc3(clusterTree(gd), sc3);
    SEQAN_ASSERT_NOT(atEnd(itSc3));
    TVertexDescriptor bb4 = *itSc3;
    SEQAN_ASSERT_EQ(property(clusterNodeCargoMap(gd), bb4), BLOCK_BIBLOCK);
    SEQAN_ASSERT_EQ(numChildren(clusterTree(gd), bb4), 0u);
    goNext(itSc3);
    SEQAN_ASSERT(atEnd(itSc3));
    // } pc2 -> sk -> sc3
    goNext(itSk);
    SEQAN_ASSERT(atEnd(itSk));
    // } pc2 -> sk
    goNext(itPc2);
    SEQAN_ASSERT(atEnd(itPc2));
    // } pc2
    goNext(itR);
    SEQAN_ASSERT(atEnd(itR));

    // Check vertex-to-building block assignment.
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), a)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), a)[0], pt3);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), b)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), b)[0], pt3);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), c1)), 2u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), c1)[0], pt3);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), c1)[1], bb1);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), c2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), c2)[0], bb1);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), c3)), 2u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), c3)[0], bb1);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), c3)[1], it);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), c4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), c4)[0], bb1);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), d1)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), d1)[0], it);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), d2)), 2u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), d2)[0], pt2);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), d2)[1], it);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), d3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), d3)[0], it);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), d4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), d4)[0], it);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), e1)), 2u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), e1)[0], bb4);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), e1)[1], it);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), e2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), e2)[0], bb4);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), e3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), e3)[0], bb4);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), e4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), e3)[0], bb4);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), f1)), 4u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f1)[0], pt1);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f1)[1], bb3);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f1)[2], bb2);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f1)[3], it);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), f2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f2)[0], bb2);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), f3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f3)[0], bb2);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), f4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f4)[0], bb2);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), f5)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f5)[0], bb3);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), f6)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f6)[0], bb3);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), f7)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), f7)[0], bb3);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), g)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), g)[0], pt1);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), h1)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), h1)[0], pt2);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), h2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), h2)[0], pt2);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), h3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), h3)[0], pt2);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), h4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), h4)[0], pt2);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), i)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), i)[0], ip1);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), j1)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), j1)[0], pt);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), j2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), j2)[0], pt);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), j3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), j3)[0], pt);
    SEQAN_ASSERT_EQ(length(property(vertexToClusterMap(gd), j4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexToClusterMap(gd), j4)[0], pt);

    // Check edge-to-building-block assignment.
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), a_c1), pt3);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), b_c1), pt3);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), c1_c2), bb1);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), c1_c4), bb1);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), c2_c3), bb1);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), c3_c4), bb1);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), c3_d1), it);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), d1_d2), it);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), d2_d3), it);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), d2_d4), it);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), d2_h1), pt2);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), d3_f1), it);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), d4_e1), it);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), e1_e2), bb4);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), e1_e3), bb4);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), e2_e4), bb4);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), e3_e4), bb4);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f1_f2), bb2);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f1_f4), bb2);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f1_g), pt1);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f2_f3), bb2);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f3_f4), bb2);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f1_f5), bb3);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f1_f7), bb3);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f5_f6), bb3);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), f6_f7), bb3);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), h1_h2), pt2);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), h1_h3), pt2);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), h1_h4), pt2);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), j1_j2), pt);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), j2_j3), pt);
    SEQAN_ASSERT_EQ(property(edgeToClusterMap(gd), j3_j4), pt);

    // TODO(holtgrew): Check flags on graph.
}

#endif  // TEST_GRAPH_DECOMPOSITION_DECOMPOSITION_H_
