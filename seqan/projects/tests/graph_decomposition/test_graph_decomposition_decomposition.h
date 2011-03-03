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

    typedef Graph<Tree<> > TTree;
    typedef typename VertexDescriptor<TTree>::Type TTreeVertexDescriptor;
    typedef String<UndirectedBuildingBlock> TBlockDescriptors;
    typedef Graph<Undirected<> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TGraphVertexDescriptor;
    typedef typename EdgeDescriptor<TGraph>::Type TGraphEdgeDescriptor;

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

    // write(std::cout, graph, DotDrawing());

    // Perform standard graph decomposition.
    typedef String<String<TTreeVertexDescriptor> > TVertexBlocksMap;
    typedef String<TTreeVertexDescriptor> TEdgeBlockMap;

    TTree clusterTree;
    TBlockDescriptors blockDescriptors;
    TVertexBlocksMap vertexBlocks;
    TEdgeBlockMap edgeBlock;
    decomposeGraphStiege(clusterTree, blockDescriptors, vertexBlocks, edgeBlock, graph);

    // Check result.
    //
    // Check clusterTree and blockDescriptors first.
    SEQAN_ASSERT_EQ(length(blockDescriptors) - 2, numVertices(clusterTree));  // Two components without stopfree kernel.

    typedef typename Iterator<TTree, AdjacencyIterator>::Type TAdjacencyIterator;
    typedef typename VertexDescriptor<TTree>::Type TVertexDescriptor;

    TVertexDescriptor r = root(clusterTree);
    TAdjacencyIterator itR(clusterTree, r);
    SEQAN_ASSERT_NOT(atEnd(itR));
    TVertexDescriptor pc1 = *itR;
    SEQAN_ASSERT_EQ(property(blockDescriptors, pc1).blockType, BLOCK_PROPER_COMPONENT);
    // { pc1
    TAdjacencyIterator itPc1(clusterTree, pc1);
    SEQAN_ASSERT_NOT(atEnd(itPc1));
    TVertexDescriptor pt = *itPc1;
    SEQAN_ASSERT_EQ(property(blockDescriptors, pt).blockType, BLOCK_PERIPHERAL_TREE);
    SEQAN_ASSERT_EQ(numChildren(clusterTree, pt), 0u);
    goNext(itPc1);
    SEQAN_ASSERT_TRUE(atEnd(itPc1));
    // } pc1
    goNext(itR);
    SEQAN_ASSERT_NOT(atEnd(itR));
    TVertexDescriptor ip1 = *itR;
    SEQAN_ASSERT_EQ(property(blockDescriptors, ip1).blockType, BLOCK_IMPROPER_COMPONENT);
    SEQAN_ASSERT_EQ(numChildren(clusterTree, ip1), 0u);
    goNext(itR);
    SEQAN_ASSERT_NOT(atEnd(itR));
    TVertexDescriptor pc2 = *itR;
    SEQAN_ASSERT_EQ(property(blockDescriptors, pc2).blockType, BLOCK_PROPER_COMPONENT);
    // { pc2
    TAdjacencyIterator itPc2(clusterTree, pc2);
    SEQAN_ASSERT_NOT(atEnd(itPc2));
    TVertexDescriptor pt1 = *itPc2;
    SEQAN_ASSERT_EQ(property(blockDescriptors, pt1).blockType, BLOCK_PERIPHERAL_TREE);
    goNext(itPc2);
    SEQAN_ASSERT_NOT(atEnd(itPc2));
    TVertexDescriptor pt2 = *itPc2;
    SEQAN_ASSERT_EQ(property(blockDescriptors, pt2).blockType, BLOCK_PERIPHERAL_TREE);
    goNext(itPc2);
    SEQAN_ASSERT_NOT(atEnd(itPc2));
    TVertexDescriptor pt3 = *itPc2;
    SEQAN_ASSERT_EQ(property(blockDescriptors, pt3).blockType, BLOCK_PERIPHERAL_TREE);
    goNext(itPc2);
    SEQAN_ASSERT_NOT(atEnd(itPc2));
    TVertexDescriptor sk = *itPc2;
    SEQAN_ASSERT_EQ(property(blockDescriptors, sk).blockType, BLOCK_STOPFREE_KERNEL);
    // { pc2 -> sk
    TAdjacencyIterator itSk(clusterTree, sk);
    SEQAN_ASSERT_NOT(atEnd(itSk));
    TVertexDescriptor it = *itSk;
    SEQAN_ASSERT_EQ(property(blockDescriptors, it).blockType, BLOCK_INTERNAL_TREE);
    SEQAN_ASSERT_EQ(numChildren(clusterTree, it), 0u);
    goNext(itSk);
    SEQAN_ASSERT_NOT(atEnd(itSk));
    TVertexDescriptor sc1 = *itSk;
    SEQAN_ASSERT_EQ(property(blockDescriptors, sc1).blockType, BLOCK_SUBCOMPONENT);
    // { pc2 -> sk -> sc1
    TAdjacencyIterator itSc1(clusterTree, sc1);
    SEQAN_ASSERT_NOT(atEnd(itSc1));
    TVertexDescriptor bb1 = *itSc1;
    SEQAN_ASSERT_EQ(property(blockDescriptors, bb1).blockType, BLOCK_BIBLOCK);
    SEQAN_ASSERT_EQ(numChildren(clusterTree, bb1), 0u);
    goNext(itSc1);
    SEQAN_ASSERT_TRUE(atEnd(itSc1));
    // } pc2 -> sk -> sc1
    goNext(itSk);
    SEQAN_ASSERT_NOT(atEnd(itSk));
    TVertexDescriptor sc2 = *itSk;
    SEQAN_ASSERT_EQ(property(blockDescriptors, sc2).blockType, BLOCK_SUBCOMPONENT);
    // { pc2 -> sk -> sc2
    TAdjacencyIterator itSc2(clusterTree, sc2);
    SEQAN_ASSERT_NOT(atEnd(itSc2));
    TVertexDescriptor bb2 = *itSc2;
    SEQAN_ASSERT_EQ(property(blockDescriptors, bb2).blockType, BLOCK_BIBLOCK);
    SEQAN_ASSERT_EQ(numChildren(clusterTree, bb2), 0u);
    goNext(itSc2);
    SEQAN_ASSERT_NOT(atEnd(itSc2));
    TVertexDescriptor bb3 = *itSc2;
    SEQAN_ASSERT_EQ(property(blockDescriptors, bb3).blockType, BLOCK_BIBLOCK);
    SEQAN_ASSERT_EQ(numChildren(clusterTree, bb3), 0u);
    goNext(itSc2);
    SEQAN_ASSERT_TRUE(atEnd(itSc2));
    // } pc2 -> sk -> sc2
    goNext(itSk);
    SEQAN_ASSERT_NOT(atEnd(itSk));
    TVertexDescriptor sc3 = *itSk;
    SEQAN_ASSERT_EQ(property(blockDescriptors, sc3).blockType, BLOCK_SUBCOMPONENT);
    // { pc2 -> sk -> sc3
    TAdjacencyIterator itSc3(clusterTree, sc3);
    SEQAN_ASSERT_NOT(atEnd(itSc3));
    TVertexDescriptor bb4 = *itSc3;
    SEQAN_ASSERT_EQ(property(blockDescriptors, bb4).blockType, BLOCK_BIBLOCK);
    SEQAN_ASSERT_EQ(numChildren(clusterTree, bb4), 0u);
    goNext(itSc3);
    SEQAN_ASSERT_TRUE(atEnd(itSc3));
    // } pc2 -> sk -> sc3
    goNext(itSk);
    SEQAN_ASSERT_TRUE(atEnd(itSk));
    // } pc2 -> sk
    goNext(itPc2);
    SEQAN_ASSERT_TRUE(atEnd(itPc2));
    // } pc2
    goNext(itR);
    SEQAN_ASSERT_TRUE(atEnd(itR));

    // Check vertex-to-building block assignment.
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, a)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, a)[0], pt3);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, b)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, b)[0], pt3);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, c1)), 2u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, c1)[0], pt3);
    SEQAN_ASSERT_EQ(property(vertexBlocks, c1)[1], bb1);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, c2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, c2)[0], bb1);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, c3)), 2u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, c3)[0], bb1);
    SEQAN_ASSERT_EQ(property(vertexBlocks, c3)[1], it);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, c4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, c4)[0], bb1);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, d1)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, d1)[0], it);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, d2)), 2u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, d2)[0], pt2);
    SEQAN_ASSERT_EQ(property(vertexBlocks, d2)[1], it);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, d3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, d3)[0], it);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, d4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, d4)[0], it);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, e1)), 2u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, e1)[0], bb4);
    SEQAN_ASSERT_EQ(property(vertexBlocks, e1)[1], it);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, e2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, e2)[0], bb4);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, e3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, e3)[0], bb4);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, e4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, e3)[0], bb4);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, f1)), 4u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f1)[0], pt1);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f1)[1], bb3);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f1)[2], bb2);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f1)[3], it);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, f2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f2)[0], bb2);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, f3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f3)[0], bb2);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, f4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f4)[0], bb2);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, f5)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f5)[0], bb3);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, f6)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f6)[0], bb3);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, f7)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, f7)[0], bb3);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, g)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, g)[0], pt1);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, h1)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, h1)[0], pt2);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, h2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, h2)[0], pt2);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, h3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, h3)[0], pt2);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, h4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, h4)[0], pt2);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, i)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, i)[0], ip1);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, j1)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, j1)[0], pt);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, j2)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, j2)[0], pt);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, j3)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, j3)[0], pt);
    SEQAN_ASSERT_EQ(length(property(vertexBlocks, j4)), 1u);
    SEQAN_ASSERT_EQ(property(vertexBlocks, j4)[0], pt);

    // Check edge-to-building-block assignment.
    SEQAN_ASSERT_EQ(property(edgeBlock, a_c1), pt3);
    SEQAN_ASSERT_EQ(property(edgeBlock, b_c1), pt3);
    SEQAN_ASSERT_EQ(property(edgeBlock, c1_c2), bb1);
    SEQAN_ASSERT_EQ(property(edgeBlock, c1_c4), bb1);
    SEQAN_ASSERT_EQ(property(edgeBlock, c2_c3), bb1);
    SEQAN_ASSERT_EQ(property(edgeBlock, c3_c4), bb1);
    SEQAN_ASSERT_EQ(property(edgeBlock, c3_d1), it);
    SEQAN_ASSERT_EQ(property(edgeBlock, d1_d2), it);
    SEQAN_ASSERT_EQ(property(edgeBlock, d2_d3), it);
    SEQAN_ASSERT_EQ(property(edgeBlock, d2_d4), it);
    SEQAN_ASSERT_EQ(property(edgeBlock, d2_h1), pt2);
    SEQAN_ASSERT_EQ(property(edgeBlock, d3_f1), it);
    SEQAN_ASSERT_EQ(property(edgeBlock, d4_e1), it);
    SEQAN_ASSERT_EQ(property(edgeBlock, e1_e2), bb4);
    SEQAN_ASSERT_EQ(property(edgeBlock, e1_e3), bb4);
    SEQAN_ASSERT_EQ(property(edgeBlock, e2_e4), bb4);
    SEQAN_ASSERT_EQ(property(edgeBlock, e3_e4), bb4);
    SEQAN_ASSERT_EQ(property(edgeBlock, f1_f2), bb2);
    SEQAN_ASSERT_EQ(property(edgeBlock, f1_f4), bb2);
    SEQAN_ASSERT_EQ(property(edgeBlock, f1_g), pt1);
    SEQAN_ASSERT_EQ(property(edgeBlock, f2_f3), bb2);
    SEQAN_ASSERT_EQ(property(edgeBlock, f3_f4), bb2);
    SEQAN_ASSERT_EQ(property(edgeBlock, f1_f5), bb3);
    SEQAN_ASSERT_EQ(property(edgeBlock, f1_f7), bb3);
    SEQAN_ASSERT_EQ(property(edgeBlock, f5_f6), bb3);
    SEQAN_ASSERT_EQ(property(edgeBlock, f6_f7), bb3);
    SEQAN_ASSERT_EQ(property(edgeBlock, h1_h2), pt2);
    SEQAN_ASSERT_EQ(property(edgeBlock, h1_h3), pt2);
    SEQAN_ASSERT_EQ(property(edgeBlock, h1_h4), pt2);
    SEQAN_ASSERT_EQ(property(edgeBlock, j1_j2), pt);
    SEQAN_ASSERT_EQ(property(edgeBlock, j2_j3), pt);
    SEQAN_ASSERT_EQ(property(edgeBlock, j3_j4), pt);

    // TODO(holtgrew): Check flags on graph.
}

#endif  // TEST_GRAPH_DECOMPOSITION_DECOMPOSITION_H_
