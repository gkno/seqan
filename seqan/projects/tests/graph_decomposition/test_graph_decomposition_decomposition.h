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

    // The same graph as in Stiege's example in the 1996 paper, Figure 1.
    TGraph graph;
    TTreeVertexDescriptor a = addVertex(graph);
    TTreeVertexDescriptor b = addVertex(graph);
    TTreeVertexDescriptor c1 = addVertex(graph);
    TTreeVertexDescriptor c2 = addVertex(graph);
    TTreeVertexDescriptor c3 = addVertex(graph);
    TTreeVertexDescriptor c4 = addVertex(graph);
    TTreeVertexDescriptor d1 = addVertex(graph);
    TTreeVertexDescriptor d2 = addVertex(graph);
    TTreeVertexDescriptor d3 = addVertex(graph);
    TTreeVertexDescriptor d4 = addVertex(graph);
    TTreeVertexDescriptor e1 = addVertex(graph);
    TTreeVertexDescriptor e2 = addVertex(graph);
    TTreeVertexDescriptor e3 = addVertex(graph);
    TTreeVertexDescriptor e4 = addVertex(graph);
    TTreeVertexDescriptor f1 = addVertex(graph);
    TTreeVertexDescriptor f2 = addVertex(graph);
    TTreeVertexDescriptor f3 = addVertex(graph);
    TTreeVertexDescriptor f4 = addVertex(graph);
    TTreeVertexDescriptor f5 = addVertex(graph);
    TTreeVertexDescriptor f6 = addVertex(graph);
    TTreeVertexDescriptor f7 = addVertex(graph);
    TTreeVertexDescriptor g = addVertex(graph);
    TTreeVertexDescriptor h1 = addVertex(graph);
    TTreeVertexDescriptor h2 = addVertex(graph);
    TTreeVertexDescriptor h3 = addVertex(graph);
    TTreeVertexDescriptor h4 = addVertex(graph);
    TTreeVertexDescriptor i = addVertex(graph);
    TTreeVertexDescriptor j1 = addVertex(graph);
    TTreeVertexDescriptor j2 = addVertex(graph);
    TTreeVertexDescriptor j3 = addVertex(graph);
    TTreeVertexDescriptor j4 = addVertex(graph);
    addEdge(graph, a, c1);
    addEdge(graph, b, c1);
    addEdge(graph, c1, c2);
    addEdge(graph, c1, c4);
    addEdge(graph, c2, c3);
    addEdge(graph, c3, c4);
    addEdge(graph, c3, d1);
    addEdge(graph, d1, d2);
    addEdge(graph, d2, d3);
    addEdge(graph, d2, d4);
    addEdge(graph, d2, h1);
    addEdge(graph, d3, f2);
    addEdge(graph, d4, e1);
    addEdge(graph, e1, e2);
    addEdge(graph, e1, e3);
    addEdge(graph, e2, e4);
    addEdge(graph, e3, e4);
    addEdge(graph, f1, f2);
    addEdge(graph, f1, f4);
    addEdge(graph, f1, g);
    addEdge(graph, f2, f3);
    addEdge(graph, f3, f4);
    addEdge(graph, f1, f5);
    addEdge(graph, f1, f7);
    addEdge(graph, f5, f6);
    addEdge(graph, f6, f7);
    addEdge(graph, h1, h2);
    addEdge(graph, h1, h3);
    addEdge(graph, h1, h4);
    addEdge(graph, j1, j2);
    addEdge(graph, j2, j3);
    addEdge(graph, j3, j4);

    (void)i; // no edge

    typedef String<String<TTreeVertexDescriptor> > TVertexBlockMap;
    typedef String<TTreeVertexDescriptor> TEdgeBlockMap;

    TTree clusterTree;
    TBlockDescriptors blockDescriptors;
    TVertexBlockMap vertexBlock;
    TEdgeBlockMap edgeBlock;
    decomposeGraphStiege(clusterTree, blockDescriptors, vertexBlock, edgeBlock, graph);
}

#endif  // TEST_GRAPH_DECOMPOSITION_DECOMPOSITION_H_
