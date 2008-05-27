#ifndef SEQAN_HEADER_TEST_GRAPH_ALIGNMENT_H
#define SEQAN_HEADER_TEST_GRAPH_ALIGNMENT_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

void Test_AlignmentGraph() {
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
	SEQAN_TASSERT(length(value(stringSet(alwaysDependStringSet),0)) == length(str0))
	clear(alwaysDependStringSet);
	SEQAN_TASSERT(length(value(ownStr,0)) == length(str0))
	assignStringSet(alwaysDependStringSet, ownStr);
	SEQAN_TASSERT(length(value(stringSet(alwaysDependStringSet),0)) == length(str0))
	clear(ownStr);
	SEQAN_TASSERT(empty(ownStr))


	TGraph g(str);
	SEQAN_TASSERT(getStringSet(g)[0] == str0)
	SEQAN_TASSERT(getStringSet(g)[1] == str1)
	SEQAN_TASSERT(stringSet(g)[2] == str2)
	assignStringSet(g, str);
	SEQAN_TASSERT(getStringSet(g)[0] == str0)
	SEQAN_TASSERT(getStringSet(g)[1] == str1)
	SEQAN_TASSERT(stringSet(g)[2] == str2)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	TVertexDescriptor v0 = addVertex(g, id1,0,2);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(sequenceId(g, v0) == id1)
	SEQAN_TASSERT(label(g, v0) == "cc")
	SEQAN_TASSERT(fragmentBegin(g, v0) == 0)
	SEQAN_TASSERT(fragmentLength(g, v0) == 2)
	SEQAN_TASSERT(findVertex(g, id1, 0) == v0)
	SEQAN_TASSERT(findVertex(g, id1, 1) == v0)
	SEQAN_TASSERT(findVertex(g, id1, 2) == nilVertex)

	TVertexDescriptor v1 = addVertex(g, id2, 0, 5);
	TEdgeDescriptor e =addEdge(g,0,1);
	SEQAN_TASSERT(_getVertexString(g)[0] == e)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e) == 1)
	SEQAN_TASSERT(sourceVertex(g, e) == 0)
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 1)
	SEQAN_TASSERT(findVertex(g, id2, 0) == v1)
	SEQAN_TASSERT(findVertex(g, id2, 1) == v1)
	SEQAN_TASSERT(findVertex(g, id2, 4) == v1)
	SEQAN_TASSERT(findVertex(g, id2, 5) == nilVertex)

	// Add more vertices and edges
	addVertex(g, id1, 10, 20);  //2
	SEQAN_TASSERT(findVertex(g, id1, 0) == v0)
	SEQAN_TASSERT(findVertex(g, id1, 1) == v0)
	SEQAN_TASSERT(findVertex(g, id1, 2) == nilVertex)
	SEQAN_TASSERT(findVertex(g, id1, 10) == 2)
	SEQAN_TASSERT(findVertex(g, id1, 19) == 2)
	SEQAN_TASSERT(findVertex(g, id1, 30) == nilVertex)
	TVertexDescriptor v3 = addVertex(g, id2, 5, 2);  //3
	addVertex(g, id1, 7, 3);  //4
	addEdge(g,3,4);
	TEdgeDescriptor my_edge = addEdge(g,3,1);
	addEdge(g,3,0);
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 3)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 1)
	SEQAN_TASSERT(numEdges(g) == 4)
	SEQAN_TASSERT(outDegree(g, v3) == 3)
	SEQAN_TASSERT(inDegree(g, v3) == 3)
	SEQAN_TASSERT(degree(g, v3) == 3)

	// Remove edges
	removeEdge(g,3,1);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 2)
	
	// Remove vertices 
	addVertex(g, id2, 14, 4);  //5
	addEdge(g,5,2);
	addEdge(g,2,3);
	addEdge(g,1,3);
	addEdge(g,1,4);
	SEQAN_TASSERT(outDegree(g, 3) == 4)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 4) == 1)
	SEQAN_TASSERT(outDegree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0);
	addEdge(g,4,1);
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g,id1,0,1);addVertex(g,id1,1,1);addVertex(g,id1,2,1);
	addVertex(g,id1,3,1);addVertex(g,id1,4,1);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clear(g);
	assignStringSet(g, str);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g,id1,0,1);addVertex(g,id1,1,1);addVertex(g,id1,2,1);
	addVertex(g,id1,3,1);addVertex(g,id1,4,1);
	addEdge(g,2,0);
	addEdge(g,4,1);
	addEdge(g,4,2);
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)

	// Transpose
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	TGraph g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy, id0, 0, 3);
	addEdge(g_copy,3,0);
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)

	// Adjacency matrix
	String<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = (unsigned int) std::sqrt((double) length(mat));
	SEQAN_TASSERT(getValue(mat,0*len+2) == 1)
	SEQAN_TASSERT(getValue(mat,3*len+2) == 0)
	SEQAN_TASSERT(getValue(mat,0*len+2) == getValue(mat,2*len+0))
	SEQAN_TASSERT(getValue(mat,1*len+4) == getValue(mat,4*len+1))
	SEQAN_TASSERT(getValue(mat,2*len+4) == getValue(mat,4*len+2))

//____________________________________________________________________________
// Alignments with edge weights
	typedef Graph<Alignment<TStringSet> > TGraph2;
	typedef VertexDescriptor<TGraph2>::Type TVertexDescriptor2;
	typedef EdgeDescriptor<TGraph2>::Type TEdgeDescriptor2;


	TGraph2 g2(str);
	SEQAN_TASSERT(numEdges(g2) == 0)
	SEQAN_TASSERT(numVertices(g2) == 0)
	SEQAN_TASSERT(empty(g2) == true)

	TVertexDescriptor2 v20 = addVertex(g2, id1,0, 2);
	TVertexDescriptor2 v21 = addVertex(g2, id2, 0, 5);
	SEQAN_TASSERT(v20 == 0)
	SEQAN_TASSERT(v21 == 1)
	TEdgeDescriptor2 e2 =addEdge(g2,0,1, 100);
	SEQAN_TASSERT(getCargo(e2) == 100)
	SEQAN_TASSERT(numEdges(g2) == 1)
	removeEdge(g2,e2);
	SEQAN_TASSERT(numEdges(g2) == 0)
	e2 =addEdge(g2,0,1, 1005);
	SEQAN_TASSERT(numEdges(g2) == 1)
	removeOutEdges(g2,0);
	SEQAN_TASSERT(numEdges(g2) == 0)
	e2 =addEdge(g2,0,1, 1005);
	SEQAN_TASSERT(findEdge(g2, 0, 1) == e2)
	SEQAN_TASSERT(numEdges(g2) == 1)
	removeInEdges(g2,0);
	SEQAN_TASSERT(numEdges(g2) == 0)



	// Vertex Iterator
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(g);
	SEQAN_TASSERT(atBegin(itV)==true)
	SEQAN_TASSERT(getValue(itV)==0)
	SEQAN_TASSERT(value(itV)==0)
	SEQAN_TASSERT(getValue(itV)==0)
	goNext(itV);
	SEQAN_TASSERT(atBegin(itV)==false)
	SEQAN_TASSERT(getValue(itV)==1)
	++itV;
	SEQAN_TASSERT(getValue(itV)==2)
	SEQAN_TASSERT(atEnd(itV)==false)
	goPrevious(itV);
	SEQAN_TASSERT((*itV)==1)
	SEQAN_TASSERT(atEnd(itV)==false)
	itV--;
	SEQAN_TASSERT(getValue(itV)==0)
	SEQAN_TASSERT(atBegin(itV)==true)

	// OutEdge Iterator
	typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator it(g, v0);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	SEQAN_TASSERT(sourceVertex(g, value(it))==0)
	SEQAN_TASSERT(targetVertex(g, *it)==2)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goPrevious(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	--it;
	it--;
	SEQAN_TASSERT(atBegin(it)==true)
	TOutEdgeIterator it2(g, v0);
	TOutEdgeIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
	SEQAN_TASSERT(&g == &hostGraph(it))

	// EdgeIterator
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEdge(g);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==2)
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==0)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==2)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==4)
	++itEdge;
	--itEdge;
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==4)
	goEnd(itEdge);
	SEQAN_TASSERT(atEnd(itEdge)==true)
	SEQAN_TASSERT(atBegin(itEdge)==false)
	goBegin(itEdge);
	SEQAN_TASSERT(atEnd(itEdge)==false)
	SEQAN_TASSERT(atBegin(itEdge)==true)

	// Adjacency Iterator
	typedef Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
	TAdjacencyIterator itAdj(g,2);
	SEQAN_TASSERT(getValue(itAdj) == 4)
	SEQAN_TASSERT(&hostGraph(itAdj) == &g)
	SEQAN_TASSERT(value(itAdj) == 4)
	SEQAN_TASSERT(*itAdj == 4)
	SEQAN_TASSERT(atEnd(itAdj)==false)
	SEQAN_TASSERT(atBegin(itAdj)==true)
	goNext(itAdj);
	SEQAN_TASSERT(*itAdj == 0)
	SEQAN_TASSERT(atEnd(itAdj)==false)
	SEQAN_TASSERT(atBegin(itAdj)==false)
	++itAdj;
	SEQAN_TASSERT(atEnd(itAdj)==true)
	SEQAN_TASSERT(atBegin(itAdj)==false)
	goPrevious(itAdj);--itAdj;
	SEQAN_TASSERT(*itAdj == 4)
	goBegin(itAdj);
	SEQAN_TASSERT(atBegin(itAdj)==true)
	goEnd(itAdj);
	SEQAN_TASSERT(atEnd(itAdj)==true)

	// Bfs Iterator
	typedef Iterator<TGraph, BfsIterator>::Type TBfsIterator;
	TBfsIterator bfsIt(g,2);
	SEQAN_TASSERT(atEnd(bfsIt)==false)
	SEQAN_TASSERT(atBegin(bfsIt)==true)
	++bfsIt;
	SEQAN_TASSERT(getValue(bfsIt) == 4)
	SEQAN_TASSERT(&hostGraph(bfsIt) == &g)
	SEQAN_TASSERT(value(bfsIt) == 4)
	SEQAN_TASSERT(*bfsIt == 4)
	goNext(bfsIt);
	SEQAN_TASSERT(value(bfsIt) == 0)

	// Dfs Iterator
	typedef Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder dfsIt(g,2);
	SEQAN_TASSERT(atEnd(dfsIt)==false)
	SEQAN_TASSERT(atBegin(dfsIt)==true)
	SEQAN_TASSERT(*dfsIt == 2)
	++dfsIt;
	SEQAN_TASSERT(getValue(dfsIt) == 0)
	SEQAN_TASSERT(&hostGraph(dfsIt) == &g)
	SEQAN_TASSERT(value(dfsIt) == 0)
	SEQAN_TASSERT(*dfsIt == 0)
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
	SEQAN_TASSERT(getFirstCoveredPosition(gAl, i3) == 9)  // Not found, length of the sequence
	SEQAN_TASSERT(getLastCoveredPosition(gAl, i3) == 0)  // Not found, length of the sequence
	SEQAN_TASSERT(getFirstCoveredPosition(gAl, i1) == 8)
	SEQAN_TASSERT(getLastCoveredPosition(gAl, i1) == 11)
	TVD vT = addVertex(gAl, i1, 13, 1);TVD vS = addVertex(gAl, i3, 6, 3);
	TVD vW = addVertex(gAl, i2, 18, 1);TVD vA = addVertex(gAl, i0, 0, 8);
	SEQAN_TASSERT(getFirstCoveredPosition(gAl, i0) == 0)
	SEQAN_TASSERT(getLastCoveredPosition(gAl, i0) == 8)
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
	SEQAN_TASSERT(getFirstCoveredPosition(gAl, i3) == 0)
	SEQAN_TASSERT(getLastCoveredPosition(gAl, i3) == 9) 
}


//////////////////////////////////////////////////////////////////////////////

void Test_OutEdgeIteratorAlignment() {
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
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	SEQAN_TASSERT(sourceVertex(g, value(it))==0)
	SEQAN_TASSERT(targetVertex(g, *it)==2)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	// Fast
	SEQAN_TASSERT(sourceVertex(it)==0)
	SEQAN_TASSERT(targetVertex(it)==2)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	++it;
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goPrevious(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	goNext(it);
	--it;
	TOutEdgeIterator it2(g, v0);
	TOutEdgeIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
	SEQAN_TASSERT(&g == &hostGraph(it))
}



//////////////////////////////////////////////////////////////////////////////

void  Test_NeedlemanWunsch() {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("annual");	assignValueById(str, str0);
	TString str1("annealing"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(0,-1,-1,0);
	int score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "annual")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "anneal")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "ing")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	int score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 3 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 3 == score2)
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;
	int score3 = globalAlignment(matches, stringSet(g), score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(length(matches) == 1)
	SEQAN_TASSERT(label(matches[0], stringSet(g), 0) == "annual")
	SEQAN_TASSERT(label(matches[0], stringSet(g), 1) == "anneal")
	SEQAN_TASSERT(score3 == score)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)


	str[0] = "annealing";
	str[1] = "annual";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "annual")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "anneal")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 7)) == "ing")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 3 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 3 == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ThisisGarfieldthecat";
	str[1] = "Garfield";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "Thisis")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 9)) == "Garfield")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "Garfield")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 17)) == "thecat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6  == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 12 == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "Garfield";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "Thisis")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 9)) == "Garfield")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "Garfield")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 17)) == "thecat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score  == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 12 == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cat";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ThisisGarfieldthe")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 18)) == "cat")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "cat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ThisisGarfieldthecat";
	str[1] = "cat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ThisisGarfieldthe")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 18)) == "cat")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "cat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Gotoh() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("ttagt");	assignValueById(str, str0);
	TString str1("ttgt"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(1,-1,-1,-2);
	int score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	int score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//double score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttgt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "tagt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttagt";
	str[1] = "tagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score + 4 == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "cttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	score2 = globalAlignment(g, score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,true,false>(), Gotoh() );
	SEQAN_TASSERT(score + 4 == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cttccagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "cc")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "ag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 7)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "ag")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 7)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "cttccagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "cc")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "ag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "ag")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 7)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttattaa";
	str[1] = "aaa";
	assignStringSet(g, str);
	Score<int> score_type2 = Score<int>(10,-1,-1,-2);
	score = globalAlignment(g, score_type2, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	score2 = globalAlignment(stringSet(g), score_type2, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 3 == score2)
	score2 = globalAlignment(stringSet(g), score_type2, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type2, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "aaa";
	str[1] = "ttattaa";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type2, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type2, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	// Note: Depending on the used recursion formula, the gotoh algorithms can differ !!!
	// Gotoh: Vertical gap and subsequent horizontal gap allowed
	// Gotoh3: Vertical gap and subsequent horizontal gap is not allowed
	//Score<double> score_scheme = Score<double>(5,-4,-0.5,-2);
	//str[0] = "tttggttt";
	//str[1] = "tttccttt";
	//assignStringSet(g, str);
	//score = globalAlignment(std::cout, str, score_scheme, Gotoh());
	//std::cout << std::endl;
	//score = globalAlignment(std::cout, str, score_scheme, Gotoh3());
}

//////////////////////////////////////////////////////////////////////////////

void Test_Hirschberg() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("ttagt");	assignValueById(str, str0);
	TString str1("ttgt"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(1,-1,-1,-2);
	int score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	int score2 = globalAlignment(stringSet(g), score_type, Hirschberg() );
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//double score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttgt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "tagt";
	str[1] = "attagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "at")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "attagt";
	str[1] = "tagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "at")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "cttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cttccagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "cc")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "ag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 7)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "ag")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 7)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "cttccagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "cc")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "ag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "ag")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 7)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttattaa";
	str[1] = "aaa";
	assignStringSet(g, str);
	Score<int> score_type2 = Score<int>(10,-1,-1,-2);
	score = globalAlignment(g, score_type2, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type2, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "aaa";
	str[1] = "ttattaa";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type2, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type2, Hirschberg() );
	//SEQAN_TASSERT(score == score2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_SmithWaterman() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("gctctgcgaata"); assignValueById(str, str0);
	TString str1("cgttgagatact"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(2,-1,-2,-2);
	int score = localAlignment(g, score_type, SmithWaterman());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "tgcg")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 8)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 9)) == "ata")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "tgag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "ata")
	SEQAN_TASSERT(numEdges(g) == 2)
	int score2 = localAlignment(stringSet(g), score_type, SmithWaterman());
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//score2 = localAlignment(std::cout, str, score_type, SmithWaterman() );
	//SEQAN_TASSERT(score == score2)
	
	str[0] = "TTGACACCCTCCCAATTGTA";
	str[1] = "ACCCCAGGCTTTACACAT";
	assignStringSet(g, str);
	score = localAlignment(g, score_type, SmithWaterman());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "TTGACAC")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 9)) == "TTTACAC")
	SEQAN_TASSERT(numEdges(g) == 1)
	score2 = localAlignment(stringSet(g), score_type, SmithWaterman());
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//score2 = localAlignment(std::cout, str, score_type, SmithWaterman() );
	//SEQAN_TASSERT(score == score2)
}
	

//////////////////////////////////////////////////////////////////////////////

void Test_SmithWatermanClump() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("gctctgcgaata"); assignValueById(str, str0);
	TString str1("cgttgagatact"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(2,-1,-2,-2);
	typedef String<Pair<int, int> > TAlignIdAndScore;
	TAlignIdAndScore alignIdScore;
	multiLocalAlignment(g, alignIdScore, score_type, 4, SmithWatermanClump() );
	
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "cgt")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "t")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 4)) == "ga")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 6)) == "g")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "at")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 9)) == "a")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 10)) == "c")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 11)) == "t")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "g")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "ct")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "c")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "t")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "gc")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 7)) == "g")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 8)) == "a")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 9)) == "at")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 11)) == "a")

	//SEQAN_TASSERT(numEdges(g) == 9)
	//SEQAN_TASSERT(numVertices(g) == 17)
		
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphAlignment() {
	// Alignment Graph
	Test_AlignmentGraph();
	Test_OutEdgeIteratorAlignment();

	// Global alignments
	Test_NeedlemanWunsch();
	Test_Gotoh();	
	Test_Hirschberg();

	// Local alignments
	Test_SmithWaterman();
	Test_SmithWatermanClump();

	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_impl_align.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_impl_align_adapt.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_config.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_needleman_wunsch.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_gotoh.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_hirschberg.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_smith_waterman.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_smith_waterman_clump.h");
	//debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_smith_waterman_island.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_interface.h");
}

}

#endif

