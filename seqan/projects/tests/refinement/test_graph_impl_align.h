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

#include <seqan/basic.h>

#ifndef SEQAN_HEADER_TEST_GRAPH_IMPL_ALIGN_H
#define SEQAN_HEADER_TEST_GRAPH_IMPL_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN {


//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(AlignmentGraphFunctions)
{
//____________________________________________________________________________
// Alignment without edge weights
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Difference between Owner and Dependent
	typedef StringSet<TString, Owner<> > TOwnerStringSet;
	TStringSet str;
	TOwnerStringSet ownStr;
	TString str0("acaagtaacataaaaaaaaaaaaaaaacccccccccttttttttaaaaa");
	TId id0 = assignValueById(str, str0);
	appendValue(ownStr, str0);
	TString str1("cccaaagggtttttccccccccccccttttttttttaaaaaaagggggggg");
	TId id1 = assignValueById(str, str1);
	appendValue(ownStr, str1);
	TString str2("cacatgtaatcatgggggggggccccccttttaaaaaaaaaaatttt");
	TId id2 = assignValueById(str, str2);
	appendValue(ownStr, str2);

	// Check that the graph makes a dependent StringSet
	TGraph alwaysDependStringSet(ownStr);
	SEQAN_ASSERT_TRUE(length(value(stringSet(alwaysDependStringSet),0)) == length(str0));
	clear(alwaysDependStringSet);
	SEQAN_ASSERT_TRUE(length(value(ownStr,0)) == length(str0));
	assignStringSet(alwaysDependStringSet, ownStr);
	SEQAN_ASSERT_TRUE(length(value(stringSet(alwaysDependStringSet),0)) == length(str0));
	clear(ownStr);
	SEQAN_ASSERT_TRUE(empty(ownStr));


	TGraph g(str);
	SEQAN_ASSERT_TRUE(getStringSet(g)[0] == str0);
	SEQAN_ASSERT_TRUE(getStringSet(g)[1] == str1);
	SEQAN_ASSERT_TRUE(stringSet(g)[2] == str2);
	assignStringSet(g, str);
	SEQAN_ASSERT_TRUE(getStringSet(g)[0] == str0);
	SEQAN_ASSERT_TRUE(getStringSet(g)[1] == str1);
	SEQAN_ASSERT_TRUE(stringSet(g)[2] == str2);
	SEQAN_ASSERT_TRUE(numEdges(g) == 0);
	SEQAN_ASSERT_TRUE(numVertices(g) == 0);
	SEQAN_ASSERT_TRUE(empty(g) == true);

	TVertexDescriptor v0 = addVertex(g, id1,0,2);
	SEQAN_ASSERT_TRUE(v0 == 0);
	SEQAN_ASSERT_TRUE(outDegree(g, v0) == 0);	
	SEQAN_ASSERT_TRUE(inDegree(g, 0) == 0);
	SEQAN_ASSERT_TRUE(degree(g, 0) == 0);
	SEQAN_ASSERT_TRUE(numVertices(g) == 1);
	SEQAN_ASSERT_TRUE(empty(g) == false);
	SEQAN_ASSERT_TRUE(sequenceId(g, v0) == id1);
	SEQAN_ASSERT_TRUE(label(g, v0) == "cc");
	SEQAN_ASSERT_TRUE(fragmentBegin(g, v0) == 0);
	SEQAN_ASSERT_TRUE(fragmentLength(g, v0) == 2);
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 0) == v0);
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 1) == v0);
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 2) == nilVertex);

	TVertexDescriptor v1 = addVertex(g, id2, 0, 5);
	TEdgeDescriptor e =addEdge(g,0,1);
	SEQAN_ASSERT_TRUE(_getVertexString(g)[0] == e);
	SEQAN_ASSERT_TRUE(getIdUpperBound(_getVertexIdManager(g)) == 2);
	SEQAN_ASSERT_TRUE(getIdUpperBound(_getEdgeIdManager(g)) == 1);
	SEQAN_ASSERT_TRUE(v1 == 1);
	SEQAN_ASSERT_TRUE(numVertices(g) == 2);
	SEQAN_ASSERT_TRUE(targetVertex(g, e) == 1);
	SEQAN_ASSERT_TRUE(sourceVertex(g, e) == 0);
	SEQAN_ASSERT_TRUE(numEdges(g) == 1);
	SEQAN_ASSERT_TRUE(outDegree(g, v0) == 1);
	SEQAN_ASSERT_TRUE(inDegree(g, 1) == 1);
	SEQAN_ASSERT_TRUE(inDegree(g, 0) == 1);	
	SEQAN_ASSERT_TRUE(degree(g, 0) == 1);
	SEQAN_ASSERT_TRUE(findVertex(g, id2, 0) == v1);
	SEQAN_ASSERT_TRUE(findVertex(g, id2, 1) == v1);
	SEQAN_ASSERT_TRUE(findVertex(g, id2, 4) == v1);
	SEQAN_ASSERT_TRUE(findVertex(g, id2, 5) == nilVertex);

	// Add more vertices and edges
	addVertex(g, id1, 10, 20);  //2
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 0) == v0);
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 1) == v0);
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 2) == nilVertex);
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 10) == 2);
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 19) == 2);
	SEQAN_ASSERT_TRUE(findVertex(g, id1, 30) == nilVertex);
	TVertexDescriptor v3 = addVertex(g, id2, 5, 2);  //3
	addVertex(g, id1, 7, 3);  //4
	addEdge(g,3,4);
	TEdgeDescriptor my_edge = addEdge(g,3,1);
	addEdge(g,3,0);
	SEQAN_ASSERT_TRUE(v3 == 3);
	SEQAN_ASSERT_TRUE(numVertices(g) == 5);
	SEQAN_ASSERT_TRUE(targetVertex(g, my_edge) == 3);
	SEQAN_ASSERT_TRUE(sourceVertex(g, my_edge) == 1);
	SEQAN_ASSERT_TRUE(numEdges(g) == 4);
	SEQAN_ASSERT_TRUE(outDegree(g, v3) == 3);
	SEQAN_ASSERT_TRUE(inDegree(g, v3) == 3);
	SEQAN_ASSERT_TRUE(degree(g, v3) == 3);

	// Remove edges
	removeEdge(g,3,1);
	removeEdge(g,0,1);
	SEQAN_ASSERT_TRUE(numEdges(g) == 2);
	
	// Remove vertices 
	addVertex(g, id2, 14, 4);  //5
	addEdge(g,5,2);
	addEdge(g,2,3);
	addEdge(g,1,3);
	addEdge(g,1,4);
	SEQAN_ASSERT_TRUE(outDegree(g, 3) == 4);
	SEQAN_ASSERT_TRUE(outDegree(g, 4) == 2);
	removeVertex(g, v3);
	SEQAN_ASSERT_TRUE(outDegree(g, 4) == 1);
	SEQAN_ASSERT_TRUE(outDegree(g, 0) == 0);
	SEQAN_ASSERT_TRUE(numVertices(g) == 5);
	SEQAN_ASSERT_TRUE(numEdges(g) == 2);

	// Clear graph
	clearEdges(g);
	SEQAN_ASSERT_TRUE(numVertices(g) == 5);
	SEQAN_ASSERT_TRUE(numEdges(g) == 0);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clearVertices(g);
	SEQAN_ASSERT_TRUE(numVertices(g) == 0);
	SEQAN_ASSERT_TRUE(numEdges(g) == 0);
	addVertex(g,id1,0,1);addVertex(g,id1,1,1);addVertex(g,id1,2,1);
	addVertex(g,id1,3,1);addVertex(g,id1,4,1);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clear(g);
	assignStringSet(g, str);
	SEQAN_ASSERT_TRUE(numVertices(g) == 0);
	SEQAN_ASSERT_TRUE(numEdges(g) == 0);
	addVertex(g,id1,0,1);addVertex(g,id1,1,1);addVertex(g,id1,2,1);
	addVertex(g,id1,3,1);addVertex(g,id1,4,1);
	addEdge(g,2,0);
	addEdge(g,4,1);
	addEdge(g,4,2);
	removeVertex(g,3);
	SEQAN_ASSERT_TRUE(numVertices(g) == 4);
	SEQAN_ASSERT_TRUE(numEdges(g) == 3);
	SEQAN_ASSERT_TRUE(outDegree(g, 4) == 2);
	SEQAN_ASSERT_TRUE(inDegree(g, 4) == 2);

	// Transpose
	transpose(g); 
	SEQAN_ASSERT_TRUE(numVertices(g) == 4);
	SEQAN_ASSERT_TRUE(numEdges(g) == 3);
	SEQAN_ASSERT_TRUE(outDegree(g, 4) == 2);
	SEQAN_ASSERT_TRUE(inDegree(g, 4) == 2);
	TGraph g_copy(g);
	SEQAN_ASSERT_TRUE(numVertices(g_copy) == 4);
	SEQAN_ASSERT_TRUE(numEdges(g_copy) == 3);
	SEQAN_ASSERT_TRUE(outDegree(g_copy, 4) == 2);
	SEQAN_ASSERT_TRUE(inDegree(g_copy, 4) == 2);
	addVertex(g_copy, id0, 0, 3);
	addEdge(g_copy,3,0);
	g_copy = g;
	SEQAN_ASSERT_TRUE(numVertices(g_copy) == 4);
	SEQAN_ASSERT_TRUE(numEdges(g_copy) == 3);
	SEQAN_ASSERT_TRUE(outDegree(g_copy, 4) == 2);
	SEQAN_ASSERT_TRUE(inDegree(g_copy, 4) == 2);
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_ASSERT_TRUE(numVertices(g_copy) == 4);
	SEQAN_ASSERT_TRUE(numEdges(g_copy) == 3);
	SEQAN_ASSERT_TRUE(outDegree(g_copy, 4) == 2);
	SEQAN_ASSERT_TRUE(inDegree(g_copy, 4) == 2);

	// Adjacency matrix
	String<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = (unsigned int) std::sqrt((double) length(mat));
	SEQAN_ASSERT_TRUE(getValue(mat,0*len+2) == 1);
	SEQAN_ASSERT_TRUE(getValue(mat,3*len+2) == 0);
	SEQAN_ASSERT_TRUE(getValue(mat,0*len+2) == getValue(mat,2*len+0));
	SEQAN_ASSERT_TRUE(getValue(mat,1*len+4) == getValue(mat,4*len+1));
	SEQAN_ASSERT_TRUE(getValue(mat,2*len+4) == getValue(mat,4*len+2));

//____________________________________________________________________________
// Alignments with edge weights
	typedef Graph<Alignment<TStringSet> > TGraph2;
	typedef VertexDescriptor<TGraph2>::Type TVertexDescriptor2;
	typedef EdgeDescriptor<TGraph2>::Type TEdgeDescriptor2;


	TGraph2 g2(str);
	SEQAN_ASSERT_TRUE(numEdges(g2) == 0);
	SEQAN_ASSERT_TRUE(numVertices(g2) == 0);
	SEQAN_ASSERT_TRUE(empty(g2) == true);

	TVertexDescriptor2 v20 = addVertex(g2, id1,0, 2);
	TVertexDescriptor2 v21 = addVertex(g2, id2, 0, 5);
	SEQAN_ASSERT_TRUE(v20 == 0);
	SEQAN_ASSERT_TRUE(v21 == 1);
	TEdgeDescriptor2 e2 =addEdge(g2,0,1, 100);
	SEQAN_ASSERT_TRUE(getCargo(e2) == 100);
	SEQAN_ASSERT_TRUE(numEdges(g2) == 1);
	removeEdge(g2,e2);
	SEQAN_ASSERT_TRUE(numEdges(g2) == 0);
	e2 =addEdge(g2,0,1, 1005);
	SEQAN_ASSERT_TRUE(numEdges(g2) == 1);
	removeOutEdges(g2,0);
	SEQAN_ASSERT_TRUE(numEdges(g2) == 0);
	e2 =addEdge(g2,0,1, 1005);
	SEQAN_ASSERT_TRUE(findEdge(g2, 0, 1) == e2);
	SEQAN_ASSERT_TRUE(numEdges(g2) == 1);
	removeInEdges(g2,0);
	SEQAN_ASSERT_TRUE(numEdges(g2) == 0);



	// Vertex Iterator
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(g);
	SEQAN_ASSERT_TRUE(atBegin(itV)==true);
	SEQAN_ASSERT_TRUE(getValue(itV)==0);
	SEQAN_ASSERT_TRUE(value(itV)==0);
	SEQAN_ASSERT_TRUE(getValue(itV)==0);
	goNext(itV);
	SEQAN_ASSERT_TRUE(atBegin(itV)==false);
	SEQAN_ASSERT_TRUE(getValue(itV)==1);
	++itV;
	SEQAN_ASSERT_TRUE(getValue(itV)==2);
	SEQAN_ASSERT_TRUE(atEnd(itV)==false);
	goPrevious(itV);
	SEQAN_ASSERT_TRUE((*itV)==1);
	SEQAN_ASSERT_TRUE(atEnd(itV)==false);
	itV--;
	SEQAN_ASSERT_TRUE(getValue(itV)==0);
	SEQAN_ASSERT_TRUE(atBegin(itV)==true);

	// OutEdge Iterator
	typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator it(g, v0);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	SEQAN_ASSERT_TRUE(sourceVertex(g, value(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, *it)==2);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	goNext(it);
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goPrevious(it);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	--it;
	it--;
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	TOutEdgeIterator it2(g, v0);
	TOutEdgeIterator it3;
	it3 = it;
	SEQAN_ASSERT_TRUE(it == it2);
	SEQAN_ASSERT_TRUE(it2 == it3);
	goEnd(it);
	SEQAN_ASSERT_TRUE(it2 != it);
	goEnd(it2);
	SEQAN_ASSERT_TRUE(it2 == it);
	goBegin(it2);
	SEQAN_ASSERT_TRUE(it2 != it);
	SEQAN_ASSERT_TRUE(&g == &hostGraph(it));

	// EdgeIterator
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEdge(g);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(itEdge))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(itEdge))==2);
	SEQAN_ASSERT_TRUE(sourceVertex(g, value(itEdge))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, *itEdge)==2);
	SEQAN_ASSERT_TRUE(atEnd(itEdge)==false);
	SEQAN_ASSERT_TRUE(atBegin(itEdge)==true);
	goNext(itEdge);
	SEQAN_ASSERT_TRUE(sourceVertex(g, value(itEdge))==1);
	SEQAN_ASSERT_TRUE(targetVertex(g, *itEdge)==4);
	++itEdge;
	--itEdge;
	SEQAN_ASSERT_TRUE(sourceVertex(g, value(itEdge))==1);
	SEQAN_ASSERT_TRUE(targetVertex(g, *itEdge)==4);
	goEnd(itEdge);
	SEQAN_ASSERT_TRUE(atEnd(itEdge)==true);
	SEQAN_ASSERT_TRUE(atBegin(itEdge)==false);
	goBegin(itEdge);
	SEQAN_ASSERT_TRUE(atEnd(itEdge)==false);
	SEQAN_ASSERT_TRUE(atBegin(itEdge)==true);

	// Adjacency Iterator
	typedef Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
	TAdjacencyIterator itAdj(g,2);
	SEQAN_ASSERT_TRUE(getValue(itAdj) == 4);
	SEQAN_ASSERT_TRUE(&hostGraph(itAdj) == &g);
	SEQAN_ASSERT_TRUE(value(itAdj) == 4);
	SEQAN_ASSERT_TRUE(*itAdj == 4);
	SEQAN_ASSERT_TRUE(atEnd(itAdj)==false);
	SEQAN_ASSERT_TRUE(atBegin(itAdj)==true);
	goNext(itAdj);
	SEQAN_ASSERT_TRUE(*itAdj == 0);
	SEQAN_ASSERT_TRUE(atEnd(itAdj)==false);
	SEQAN_ASSERT_TRUE(atBegin(itAdj)==false);
	++itAdj;
	SEQAN_ASSERT_TRUE(atEnd(itAdj)==true);
	SEQAN_ASSERT_TRUE(atBegin(itAdj)==false);
	goPrevious(itAdj);--itAdj;
	SEQAN_ASSERT_TRUE(*itAdj == 4);
	goBegin(itAdj);
	SEQAN_ASSERT_TRUE(atBegin(itAdj)==true);
	goEnd(itAdj);
	SEQAN_ASSERT_TRUE(atEnd(itAdj)==true);

	// Bfs Iterator
	typedef Iterator<TGraph, BfsIterator>::Type TBfsIterator;
	TBfsIterator bfsIt(g,2);
	SEQAN_ASSERT_TRUE(atEnd(bfsIt)==false);
	SEQAN_ASSERT_TRUE(atBegin(bfsIt)==true);
	++bfsIt;
	SEQAN_ASSERT_TRUE(getValue(bfsIt) == 4);
	SEQAN_ASSERT_TRUE(&hostGraph(bfsIt) == &g);
	SEQAN_ASSERT_TRUE(value(bfsIt) == 4);
	SEQAN_ASSERT_TRUE(*bfsIt == 4);
	goNext(bfsIt);
	SEQAN_ASSERT_TRUE(value(bfsIt) == 0);

	// Dfs Iterator
	typedef Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder dfsIt(g,2);
	SEQAN_ASSERT_TRUE(atEnd(dfsIt)==false);
	SEQAN_ASSERT_TRUE(atBegin(dfsIt)==true);
	SEQAN_ASSERT_TRUE(*dfsIt == 2);
	++dfsIt;
	SEQAN_ASSERT_TRUE(getValue(dfsIt) == 0);
	SEQAN_ASSERT_TRUE(&hostGraph(dfsIt) == &g);
	SEQAN_ASSERT_TRUE(value(dfsIt) == 0);
	SEQAN_ASSERT_TRUE(*dfsIt == 0);
	goNext(dfsIt);

	// Alignments
	typedef String<char> TAlignString;
	typedef StringSet<TAlignString, Dependent<> > TAlignStringSet;
	typedef Graph<Alignment<TAlignStringSet, void> > TAlignmentGraph;
	typedef VertexDescriptor<TAlignmentGraph>::Type TVD;
	typedef EdgeDescriptor<TAlignmentGraph>::Type TED;
	
	TAlignStringSet al;
	TAlignString al0("Garfieldthelastfatcat");
	TId i0 = assignValueById(al, al0);
	TAlignString al1("Garfieldthefastcat");
	TId i1 = assignValueById(al, al1);
	TAlignString al2("Garfieldtheveryfastcat");
	TId i2 = assignValueById(al, al2);
	TAlignString al3("thefatcat");
	TId i3 = assignValueById(al, al3);

	TAlignmentGraph gAl(al);
	TVD vH = addVertex(gAl, i1, 8, 3);
	SEQAN_ASSERT_TRUE(getFirstCoveredPosition(gAl, i3) == 9);  // Not found, length of the sequence
	SEQAN_ASSERT_TRUE(getLastCoveredPosition(gAl, i3) == 0);  // Not found, length of the sequence
	SEQAN_ASSERT_TRUE(getFirstCoveredPosition(gAl, i1) == 8);
	SEQAN_ASSERT_TRUE(getLastCoveredPosition(gAl, i1) == 11);
	TVD vT = addVertex(gAl, i1, 13, 1);TVD vS = addVertex(gAl, i3, 6, 3);
	TVD vW = addVertex(gAl, i2, 18, 1);TVD vA = addVertex(gAl, i0, 0, 8);
	SEQAN_ASSERT_TRUE(getFirstCoveredPosition(gAl, i0) == 0);
	SEQAN_ASSERT_TRUE(getLastCoveredPosition(gAl, i0) == 8);
	TVD vM = addVertex(gAl, i2, 11, 4);
	TVD vK = addVertex(gAl, i2, 0, 8);TVD vC = addVertex(gAl, i0, 11, 4);TVD vD = addVertex(gAl, i0, 15, 2);
	TVD vF = addVertex(gAl, i0, 18, 3);TVD vG = addVertex(gAl, i1, 0, 8);addEdge(gAl, vA, vG);
	TVD vI = addVertex(gAl, i1, 11, 2);TVD vQ = addVertex(gAl, i3, 3, 2);TVD vB = addVertex(gAl, i0, 8, 3);
	TVD vU = addVertex(gAl, i1, 14, 1);TVD vE = addVertex(gAl, i0, 17, 1);TVD vJ = addVertex(gAl, i1, 15, 3);
	TVD vL = addVertex(gAl, i2, 8, 3);
	addEdge(gAl, vH, vL);
	TVD vN = addVertex(gAl, i2, 15, 2);TVD vV = addVertex(gAl, i2, 17, 1);
	TVD vO = addVertex(gAl, i2, 19, 3);TVD vP = addVertex(gAl, i3, 0, 3);TVD vR = addVertex(gAl, i3, 5, 1);
	addEdge(gAl, vA, vK);addEdge(gAl, vG, vK);addEdge(gAl, vB, vH);addEdge(gAl, vB, vL);
	addEdge(gAl, vB, vP);addEdge(gAl, vH, vP);addEdge(gAl, vL, vP);addEdge(gAl, vC, vM);
	addEdge(gAl, vD, vI);addEdge(gAl, vD, vQ);addEdge(gAl, vD, vN);addEdge(gAl, vI, vQ);
	addEdge(gAl, vI, vN);addEdge(gAl, vQ, vN);addEdge(gAl, vT, vV);addEdge(gAl, vE, vU);
	addEdge(gAl, vE, vW);addEdge(gAl, vE, vR);addEdge(gAl, vU, vW);addEdge(gAl, vU, vR);
	addEdge(gAl, vW, vR);addEdge(gAl, vF, vJ);addEdge(gAl, vF, vO);addEdge(gAl, vF, vS);
	addEdge(gAl, vJ, vO);addEdge(gAl, vJ, vS);addEdge(gAl, vO, vS);
	SEQAN_ASSERT_TRUE(getFirstCoveredPosition(gAl, i3) == 0);
	SEQAN_ASSERT_TRUE(getLastCoveredPosition(gAl, i3) == 9);

	// Output of the alignment graph
	std::cout << gAl << std::endl;
	StringSet<String<char> > seqs;
	appendValue(seqs, "seq1");
	appendValue(seqs, "seq2");
	appendValue(seqs, "seq3");
	appendValue(seqs, "seq4");
	write(std::cout,gAl,seqs,FastaFormat());
	write(std::cout,gAl,seqs,MsfFormat());
	write(std::cout,gAl,seqs,CgVizFormat());
}



//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(HeaviestCommonSubsequence)
{
	typedef String<AminoAcid> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, int> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	TString s1 = "aaa";
	TString s2 = "aa";
	TStringSet strSet;
	assignValueById(strSet, s1);
	assignValueById(strSet, s2);
	TGraph g(strSet);
	addVertex(g, 0, 0, 1);
	addVertex(g, 0, 1, 1);
	addVertex(g, 0, 2, 1);
	addVertex(g, 1, 0, 1);
	addVertex(g, 1, 1, 1);
	addEdge(g, 0, 3, 10); addEdge(g, 0, 4, 15);
	addEdge(g, 1, 3, 10); addEdge(g, 1, 4, 10);
	addEdge(g, 2, 3, 15); addEdge(g, 2, 4, 10);
	String<String<TVertexDescriptor> > str1;
	String<String<TVertexDescriptor> > str2;
	String<String<TVertexDescriptor> > align;
	String<TVertexDescriptor> tmp;
	clear(tmp); appendValue(tmp, 0); appendValue(str1, tmp);
	clear(tmp); appendValue(tmp, 1); appendValue(str1, tmp);
	clear(tmp); appendValue(tmp, 2); appendValue(str1, tmp);
	clear(tmp); appendValue(tmp, 3); appendValue(str2, tmp);
	clear(tmp); appendValue(tmp, 4); appendValue(str2, tmp);
	SEQAN_ASSERT_TRUE(heaviestCommonSubsequence(g, str1, str2, align) == 20);

	s1 = "aaaaa";
	s2 = "aaa";
	clear(strSet);
	assignValueById(strSet, s1);
	assignValueById(strSet, s2);
	assignStringSet(g, strSet);
	addVertex(g, 0, 0, 2);
	addVertex(g, 0, 2, 1);
	addVertex(g, 0, 3, 2);
	addVertex(g, 1, 0, 1);
	addVertex(g, 1, 1, 2);
	addEdge(g, 0, 4, 10); addEdge(g, 1, 3, 20);
	clear(align);
	clear(str1);
	clear(str2);
	clear(tmp); appendValue(tmp, 0); appendValue(str1, tmp);
	clear(tmp); appendValue(tmp, 1); appendValue(str1, tmp);
	clear(tmp); appendValue(tmp, 2); appendValue(str1, tmp);
	clear(tmp); appendValue(tmp, 3); appendValue(str2, tmp);
	clear(tmp); appendValue(tmp, 4); appendValue(str2, tmp);
	heaviestCommonSubsequence(g, str1, str2, align);
}




//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(OutEdgeIteratorAlignment)
{
//____________________________________________________________________________
// Graph AlignmentOutEdgeIterator
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	TString str1 = "aa";
	TString str2 = "ac";
	TStringSet strSet;
	appendValue(strSet, str1);
	appendValue(strSet, str2);
	TGraph g(strSet);
	TVertexDescriptor v0 = addVertex(g,0,0,1);
	addVertex(g,0,1,1);
	addVertex(g,1,0,1);
	addEdge(g,0,2,10);
	addVertex(g,1,1,1);
	
	typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator it(g, v0);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	SEQAN_ASSERT_TRUE(sourceVertex(g, value(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, *it)==2);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	// Slow
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	// Fast
	SEQAN_ASSERT_TRUE(sourceVertex(it)==0);
	SEQAN_ASSERT_TRUE(targetVertex(it)==2);
	SEQAN_ASSERT_TRUE(atEnd(it)==false);
	SEQAN_ASSERT_TRUE(atBegin(it)==true);
	++it;
	SEQAN_ASSERT_TRUE(atEnd(it)==true);
	SEQAN_ASSERT_TRUE(atBegin(it)==false);
	goPrevious(it);
	SEQAN_ASSERT_TRUE(sourceVertex(g, getValue(it))==0);
	SEQAN_ASSERT_TRUE(targetVertex(g, getValue(it))==2);
	goNext(it);
	--it;
	TOutEdgeIterator it2(g, v0);
	TOutEdgeIterator it3;
	it3 = it;
	SEQAN_ASSERT_TRUE(it == it2);
	SEQAN_ASSERT_TRUE(it2 == it3);
	goEnd(it);
	SEQAN_ASSERT_TRUE(it2 != it);
	goEnd(it2);
	SEQAN_ASSERT_TRUE(it2 == it);
	goBegin(it2);
	SEQAN_ASSERT_TRUE(it2 != it);
	SEQAN_ASSERT_TRUE(&g == &hostGraph(it));
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif

