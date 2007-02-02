#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>
#include <string>

#define SEQAN_DEBUG
#define SEQAN_TEST

#define TEST_PATH "projects/tests/graph/"
#define LIB_PATH "projects/library/seqan/graph/"

#include <seqan/graph.h>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

void Test_IdManager() {
//____________________________________________________________________________
// IdManager
	typedef Id<IdManager<> >::Type TIdType;
	IdManager<> idm;
	
	// Obtain Ids
	TIdType id=obtainId(idm);
	SEQAN_TASSERT(id == 0)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 1)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 2)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 3)
	id=obtainId(idm);
	SEQAN_TASSERT(id == 4)
	SEQAN_TASSERT(getIdUpperBound(idm) == 5)
	SEQAN_TASSERT(idCount(idm) == 5)

	// Release Ids
	releaseId(idm,3);
	SEQAN_TASSERT(idInUse(idm,3) == false)
	releaseId(idm,1);
	SEQAN_TASSERT(idInUse(idm,1) == false)
	releaseId(idm,2);
	SEQAN_TASSERT(idInUse(idm,2) == false)
	SEQAN_TASSERT(idCount(idm) == 2)
	SEQAN_TASSERT(getIdUpperBound(idm) == 5)
	// Ensure a minimum id range
	releaseId(idm,4);
	SEQAN_TASSERT(idCount(idm) == 1)
	SEQAN_TASSERT(getIdUpperBound(idm) == 1)

	// Smallest Ids are reused first
	id=obtainId(idm);
	id=obtainId(idm);
	id=obtainId(idm);
	id=obtainId(idm);
	SEQAN_TASSERT(getIdUpperBound(idm) == 5)
	releaseId(idm,3);
	SEQAN_TASSERT(getIdLowerBound(idm) == 0)
	releaseId(idm,0);
	SEQAN_TASSERT(idInUse(idm,0) == false)
	SEQAN_TASSERT(getIdLowerBound(idm) == 1)
	releaseId(idm,2);
	id=obtainId(idm);
	SEQAN_TASSERT(id == 0)
	SEQAN_TASSERT(idInUse(idm,0) == true)
	
	// Check copy constructor and assignment operator
	IdManager<> idm2(idm);
	SEQAN_TASSERT(idCount(idm2) == 3)
	id=obtainId(idm2);
	SEQAN_TASSERT(id == 2)
	idm = idm2;
	SEQAN_TASSERT(idCount(idm) == 4)
	releaseAll(idm);
	SEQAN_TASSERT(idCount(idm) == 0)
	SEQAN_TASSERT(getIdUpperBound(idm) == 0)
	SEQAN_TASSERT(getIdLowerBound(idm) == 0)


//____________________________________________________________________________
// Dummy IdManager

	// Dummy IdManager
	typedef Id<IdManager<void> >::Type TIdType;
	IdManager<void> id_dummy;
	
	// Obtain Ids
	TIdType idd=obtainId(id_dummy);
	SEQAN_TASSERT(idd == 0)
	idd=obtainId(id_dummy); // Always zero
	SEQAN_TASSERT(idd == 0)
	SEQAN_TASSERT(idInUse(id_dummy, 1) == false) // Always false
	SEQAN_TASSERT(idInUse(id_dummy, 2) == false)
	idd=obtainId(id_dummy);
	SEQAN_TASSERT(idd == 0)
	idd=obtainId(id_dummy);
	SEQAN_TASSERT(idd == 0)
	idd=obtainId(id_dummy);
	SEQAN_TASSERT(idd == 0)
	SEQAN_TASSERT(getIdUpperBound(id_dummy) == 5) 
	SEQAN_TASSERT(getIdLowerBound(id_dummy) == 0)
	SEQAN_TASSERT(idCount(id_dummy) == 5) //But: Correct id count

	// Release Ids
	releaseId(id_dummy,3); 
	SEQAN_TASSERT(idCount(id_dummy) == 4)
	releaseId(id_dummy,1); 
	SEQAN_TASSERT(idCount(id_dummy) == 3)
	IdManager<void> id_dummy2(id_dummy);
	SEQAN_TASSERT(idCount(id_dummy2) == 3)
	idd=obtainId(id_dummy2);
	SEQAN_TASSERT(idd == 0)
	id_dummy = id_dummy2;
	SEQAN_TASSERT(idCount(id_dummy) == 4)
	releaseAll(id_dummy);
	SEQAN_TASSERT(idCount(idm) == 0)
	SEQAN_TASSERT(getIdUpperBound(idm) == 0)
	SEQAN_TASSERT(getIdLowerBound(idm) == 0)
}

//////////////////////////////////////////////////////////////////////////////

void Test_EdgeStump() {
//____________________________________________________________________________
// EdgeStump
	// Default: EdgeStump without edge cargo but with edge id
	EdgeStump<> es1;
	_assignId(&es1, 4);
	SEQAN_TASSERT(_getId(&es1) == 4)
	assignTarget(&es1, 5);
	SEQAN_TASSERT(getTarget(&es1) == 5)
	assignCargo(&es1, 15);  //Does nothing
	// No cargo -> pointer to 0
	SEQAN_TASSERT(getCargo(&es1) == (void*) 0)
	SEQAN_TASSERT(_getId(&es1) == 4)
	EdgeStump<> const es_const(es1);
	SEQAN_TASSERT(getCargo(&es1) == (void*) 0)
	SEQAN_TASSERT(getTarget(&es_const) == 5)
	
	// EdgeStump with edge cargo and with edge id
	EdgeStump<unsigned int> es2;
	assignTarget(&es2, 1);
	SEQAN_TASSERT(getTarget(&es2) == 1)
	_assignId(&es2, 1);
	SEQAN_TASSERT(_getId(&es2) == 1)
	assignCargo(&es2, 1);
	SEQAN_TASSERT(getCargo(&es2) == 1)
	cargo(&es2) = 3;
	SEQAN_TASSERT(getCargo(&es2) == 3)
	EdgeStump<unsigned int> const es3(es2);
	SEQAN_TASSERT(_getId(&es2) == 1)
	SEQAN_TASSERT(getCargo(&es3) == 3)
	SEQAN_TASSERT(cargo(&es3) == 3)
	SEQAN_TASSERT(target(&es3) == 1)

	// EdgeStump without edge cargo and without edge id
	EdgeStump<void, WithoutEdgeId> es4;
	assignTarget(&es4, 3);
	SEQAN_TASSERT(getTarget(&es4) == 3)
	SEQAN_TASSERT(target(&es4) == 3)
	target(&es4) = 1;
	SEQAN_TASSERT(getTarget(&es4) == 1)
	_assignId(&es4, 1); //Does nothing
	SEQAN_TASSERT(_getId(&es4) == 0)
	assignCargo(&es4, 1); //Does nothing
	SEQAN_TASSERT(getCargo(&es4) == (void*) 0)
	SEQAN_TASSERT(cargo(&es4) == (void*) 0); //Does nothing
	EdgeStump<void, WithoutEdgeId> const es5(es4);
	SEQAN_TASSERT(_getId(&es5) == 0)
	SEQAN_TASSERT(getCargo(&es5) == 0)
	SEQAN_TASSERT(cargo(&es5) == 0)

	// EdgeStump with edge cargo and without edge id
	EdgeStump<unsigned int, WithoutEdgeId> es6;
	assignTarget(&es6, 1);
	SEQAN_TASSERT(getTarget(&es6) == 1)
	SEQAN_TASSERT(target(&es6) == 1)
	_assignId(&es6, 1); //Does nothing
	SEQAN_TASSERT(_getId(&es6) == 0)
	assignCargo(&es6, 1); 
	SEQAN_TASSERT(getCargo(&es6) == 1)
	SEQAN_TASSERT(cargo(&es6) == 1);
	EdgeStump<unsigned int, WithoutEdgeId> const es7(es6);
	SEQAN_TASSERT(_getId(&es7) == 0)
	SEQAN_TASSERT(getCargo(&es7) == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_EdgeAutomaton() {
//____________________________________________________________________________
// EdgeAutomaton
	// Default edge automaton: No cargo
	EdgeAutomaton<> ea1;
	assignTarget(ea1, 5);
	SEQAN_TASSERT(getTarget(ea1) == 5)
	assignCargo(ea1, 15);  
	SEQAN_TASSERT(getCargo(ea1) == (void*) 0)
	SEQAN_TASSERT(cargo(ea1) == (void*) 0)
	EdgeAutomaton<> const ea_const(ea1);
	SEQAN_TASSERT(getCargo(ea1) == (void*) 0)
	SEQAN_TASSERT(getCargo(ea_const) == (void*) 0)
	SEQAN_TASSERT(cargo(ea_const) == (void*) 0)
	//SEQAN_TASSERT(cargo(ea_const) == (void*) 0)
	SEQAN_TASSERT(getTarget(ea_const) == 5)
	SEQAN_TASSERT(target(ea_const) == 5)


	// EdgeAutomaton with edge cargo
	EdgeAutomaton<unsigned int> ea2;
	assignTarget(ea2, 1);
	SEQAN_TASSERT(getTarget(ea2) == 1)
	target(ea2) = 4;
	SEQAN_TASSERT(getTarget(ea2) == 4)
	assignCargo(ea2, 1);
	SEQAN_TASSERT(getCargo(ea2) == 1)
	cargo(ea2) = 3;
	SEQAN_TASSERT(getCargo(ea2) == 3)
	EdgeAutomaton<unsigned int> const ea3(ea2);
	SEQAN_TASSERT(getCargo(ea3) == 3)
	SEQAN_TASSERT(cargo(ea3) == 3)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Automaton() {
//____________________________________________________________________________
// Standard automaton: No edge cargo

	typedef Graph<Automaton<Dna> > StandardAutomaton;
	typedef VertexDescriptor<StandardAutomaton>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardAutomaton>::Type TEdgeDescriptor;
	

	StandardAutomaton g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)

	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v0,'a');
	SEQAN_TASSERT(_getId(e1) == 0)
	SEQAN_TASSERT(targetVertex(g, e1) == 0)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0) 
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 1)
	SEQAN_TASSERT(degree(g, v0) == 2)
	
	// Add further edges and vertices
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1,'g');
	SEQAN_TASSERT(_getId(e2) == 2)
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 3)

	// Test a property map
	std::string nameEd[] = {"ar", "ae"};
	String<std::string> edMap;
	initEdgeMap(g,edMap, nameEd);
	SEQAN_TASSERT(getProperty(edMap, e1) == "ar")
	SEQAN_TASSERT(getProperty(edMap, e2) == "ae")

	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
	addEdge(g,3,4,'a');
	TEdgeDescriptor my_edge = addEdge(g,3,1,'c');
	SEQAN_TASSERT(_getId(my_edge) == 13)
	addEdge(g,3,0,'t');
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 1)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 3)
	SEQAN_TASSERT(numEdges(g) == 5)
	SEQAN_TASSERT(outDegree(g, v3) == 3)	

	// Remove edges
	removeEdge(g,my_edge);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 3)

	// Remove vertices 
	TEdgeDescriptor e3 = addEdge(g,3,3,'a');
	addEdge(g,1,3,'a');
	addEdge(g,0,3,'c');
	addEdge(g,0,4,'t');
	SEQAN_TASSERT(outDegree(g, 0) == 3)
	SEQAN_TASSERT(outDegree(g, 1) == 1)
	SEQAN_TASSERT(targetVertex(g, e3) == 3)
	SEQAN_TASSERT(sourceVertex(g, e3) == 3)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 0) == 2)
	SEQAN_TASSERT(outDegree(g, 1) == 0)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0,'a');
	addEdge(g,4,1,'c');
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0,'t');
	addEdge(g,4,1,'g');
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0,'c');
	addEdge(g,4,1,'g');
	addEdge(g,4,2,'t');
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 0)

	// Transpose
	StandardAutomaton g_tmp;
	//Transposes the graph in-place
	transpose(g, g_tmp); 
	g = g_tmp;
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 0)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardAutomaton g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
	addEdge(g_copy,3,0,'a');
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 0)
	removeVertex(g,0);

	// Adjacency matrix
	Matrix<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	typedef Iterator<Matrix<TVertexDescriptor>,PositionIterator>::Type TMatrixIterator;
	TMatrixIterator p = begin(mat);
	TMatrixIterator endP = end(mat);
	for(;p != endP;goNext(p)) {
		if (*p == 1) {
			// Both edges go to vertex 4
			SEQAN_TASSERT(coordinate(p,0) == 4)
		}
	}

	// Test iterators
	typedef Iterator<StandardAutomaton, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator itVert(g);
	SEQAN_TASSERT(getValue(itVert) == 1)
	++itVert;
	SEQAN_TASSERT(getValue(itVert) == 2)
	itVert++;
	SEQAN_TASSERT(getValue(itVert) == 4)
	goNext(itVert);
	SEQAN_TASSERT(atEnd(itVert) == true)

	addEdge(g,1,2,'T');
	typedef Iterator<StandardAutomaton, OutEdgeIterator<> >::Type TOutEdgeIterator;
	TOutEdgeIterator itEdge(g,1);
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==4)
	// Fast
	SEQAN_TASSERT(sourceVertex(itEdge)==1)
	SEQAN_TASSERT(targetVertex(itEdge)==4)
	
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==4)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(itEdge)==1)
	SEQAN_TASSERT(targetVertex(itEdge)==2)
	++itEdge;
	itEdge++;
	SEQAN_TASSERT(atEnd(itEdge)==true)
	SEQAN_TASSERT(atBegin(itEdge)==false)
	goPrevious(itEdge);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==2)
	--itEdge;
	itEdge++; 
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==2)
	itEdge--;
	itEdge--;
	SEQAN_TASSERT(atBegin(itEdge)==true)
	TOutEdgeIterator itEdge2(itEdge);
	TOutEdgeIterator itEdge3;
	itEdge3 = itEdge;
	SEQAN_TASSERT(itEdge == itEdge2)
	SEQAN_TASSERT(itEdge2 == itEdge3)
	goEnd(itEdge);
	SEQAN_TASSERT(itEdge2 != itEdge)
	goEnd(itEdge2);
	SEQAN_TASSERT(itEdge2 == itEdge)
	goBegin(itEdge2);
	SEQAN_TASSERT(atBegin(itEdge2)==true)
	SEQAN_TASSERT(&g == &hostGraph(itEdge))

	typedef Iterator<StandardAutomaton, EdgeIterator<> >::Type TEdgeIterator;
	TEdgeIterator itEd(g);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	SEQAN_TASSERT(sourceVertex(g, value(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, *itEd)==4)
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==true)
	goNext(itEd);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==false)
	++itEd;
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==false)
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	// Fast
	SEQAN_TASSERT(sourceVertex(itEd)==2)
	SEQAN_TASSERT(targetVertex(itEd)==4)
	itEd++;
	itEd++;
	SEQAN_TASSERT(atEnd(itEd)==true)
	SEQAN_TASSERT(atBegin(itEd)==false)
	goPrevious(itEd);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	--itEd;
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==2)
	TEdgeIterator itEd2(g);
	TEdgeIterator itEd3;
	goBegin(itEd);
	itEd3 = itEd;
	SEQAN_TASSERT(itEd == itEd2)
	SEQAN_TASSERT(itEd2 == itEd3)
	goEnd(itEd);
	SEQAN_TASSERT(itEd2 != itEd)
	goEnd(itEd2);
	SEQAN_TASSERT(itEd2 == itEd)
	goBegin(itEd2);
	SEQAN_TASSERT(itEd2 != itEd)
	SEQAN_TASSERT(&hostGraph(itEd) == &g)

	typedef Iterator<StandardAutomaton, AdjacencyIterator<> >::Type TAdjacencyIterator;
	TAdjacencyIterator itAd(g,1);
	SEQAN_TASSERT(getValue(itAd) == 4)
	SEQAN_TASSERT(&hostGraph(itAd) == &g)
	SEQAN_TASSERT(value(itAd) == 4)
	SEQAN_TASSERT(*itAd == 4)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==true)
	goNext(itAd);
	SEQAN_TASSERT(getValue(itAd)==2)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==false)
	++itAd;
	SEQAN_TASSERT(atEnd(itAd)==true)
	SEQAN_TASSERT(atBegin(itAd)==false)
	goBegin(itAd);
	itAd++;
	itAd++;
	itAd++;
	SEQAN_TASSERT(atEnd(itAd)==true)
	SEQAN_TASSERT(atBegin(itAd)==false)
	goPrevious(itAd);
	SEQAN_TASSERT(getValue(itAd)==2)
	--itAd;
	SEQAN_TASSERT(getValue(itAd)==4)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==true)
	goEnd(itAd);
	itAd--;
	SEQAN_TASSERT(getValue(itAd)==2)
	goBegin(itAd);
	TAdjacencyIterator itAd2(itAd);
	TAdjacencyIterator itAd3;
	itAd3 = itAd;
	SEQAN_TASSERT(itAd == itAd2)
	SEQAN_TASSERT(itAd2 == itAd3)
	goEnd(itAd);
	SEQAN_TASSERT(itAd2 != itAd)
	goEnd(itAd2);
	SEQAN_TASSERT(itAd2 == itAd)
	goBegin(itAd2);
	SEQAN_TASSERT(itAd2 != itAd)



//____________________________________________________________________________
// Automaton
	typedef VertexDescriptor<Graph<Automaton<char> > >::Type VertexDescriptorType;
	typedef EdgeDescriptor<Graph<Automaton<char> > >::Type EdgeDescriptorType;
	Graph<Automaton<char> > automaton;
	VertexDescriptorType rootVertex = addVertex(automaton); // A = 0
	addVertex(automaton); // B = 1
	addVertex(automaton); // C = 2
	addVertex(automaton); // D = 3
	addVertex(automaton); // E = 4
	addVertex(automaton); // F = 5
	addEdge(automaton,0,1,'2');
	addEdge(automaton,1,0,'1');
	addEdge(automaton,4,0,'6');
	addEdge(automaton,0,3,'7');
	addEdge(automaton,1,1,'3');
	addEdge(automaton,1,2,'4');
	addEdge(automaton,5,1,'8');
	addEdge(automaton,2,5,'5');
	addEdge(automaton,3,4,'2');
	addEdge(automaton,5,3,'7');

	VertexDescriptorType succ;
	succ = getSuccessor(automaton,rootVertex,'7');
	SEQAN_TASSERT(succ == 3)
	// Throws an error in debug mode because edge does not exist
	//succ = getSuccessor(automaton,rootVertex,'6');
	succ = getSuccessor(automaton,succ,'2');
	SEQAN_TASSERT(succ == 4)
	succ = getSuccessor(automaton,succ,'6');
	SEQAN_TASSERT(succ == 0)
	// If no map is specified it is assumed that an edge cargo exists!!!
	succ = getSuccessor(automaton,succ,'2');
	SEQAN_TASSERT(succ == 1)
	// Go backwards
	VertexDescriptorType pred;
	pred = getPredecessor(automaton,succ,'3');
	SEQAN_TASSERT(pred == 1)
	pred = getPredecessor(automaton,pred,'8');
	SEQAN_TASSERT(pred == 5)
	// If no map is specified it is assumed that an edge cargo exists!!!
	pred = getPredecessor(automaton,pred,'5');
	SEQAN_TASSERT(pred == 2)

	// Now using shortcuts
	succ = parseString(automaton,rootVertex,"7262");
	SEQAN_TASSERT(succ == 1)
	std::string str = "7262";
	succ = parseString(automaton,rootVertex, str.begin(), str.end());
	SEQAN_TASSERT(succ == 1)
	String<char> str2("7262");
	succ = parseString(automaton,rootVertex, begin(str2), end(str2));
	SEQAN_TASSERT(succ == 1)
	String<char> input("7262");
	succ = parseString(automaton,rootVertex, input);
	SEQAN_TASSERT(succ == 1)

	// Additional cargo
	typedef Graph<Automaton<Dna, short> > TGraph9;
	typedef VertexDescriptor<TGraph9>::Type TVertexDescriptor9;
	typedef EdgeDescriptor<TGraph9>::Type TEdgeDescriptor9;
	typedef Size<TGraph9>::Type TSize9;

	TGraph9 g9;
	TVertexDescriptor9 vert0 = addVertex(g9);
	TVertexDescriptor9 vert1 = addVertex(g9);
	Dna aDna('a');
	Dna gDna('g');
	TEdgeDescriptor9 edg1 = addEdge(g9,0,1,aDna,12);
	TEdgeDescriptor9 edg2 = addEdge(g9,1,0,gDna,21);
	TGraph9 g10;
	transpose(g9, g10);
	SEQAN_TASSERT(getCargo(g9.data_vertex[0].data_edge[(TSize9) aDna])==12)
	SEQAN_TASSERT(getCargo(g9.data_vertex[1].data_edge[(TSize9) gDna])==21)
	SEQAN_TASSERT(getCargo(g10.data_vertex[1].data_edge[(TSize9) aDna])==12)
	SEQAN_TASSERT(getCargo(g10.data_vertex[0].data_edge[(TSize9) gDna])==21)
}

//////////////////////////////////////////////////////////////////////////////

void Test_WordGraph() {
//____________________________________________________________________________
// Standard automaton: No edge cargo

	typedef Graph<Automaton<Dna, void, WordGraph<> > > TWordGraph;
	typedef VertexDescriptor<TWordGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TWordGraph>::Type TEdgeDescriptor;
	

	TWordGraph g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	addVertex(g);
	addVertex(g);
	TVertexDescriptor v3 = addVertex(g);


	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v3,"ag");
	SEQAN_TASSERT(_getId(e1) == 0)
	SEQAN_TASSERT(targetVertex(g, e1) == 3)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0) 
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 0)
	SEQAN_TASSERT(degree(g, v0) == 1)

	// Add further edges and vertices
	addVertex(g);
	TVertexDescriptor v5 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,5,"g");
	SEQAN_TASSERT(_getId(e2) == 2) //Third letter, first vertex --> 0*4 + 2 
	SEQAN_TASSERT(v5 == 5)
	SEQAN_TASSERT(numVertices(g) == 6)
	SEQAN_TASSERT(targetVertex(g, e2) == 5)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 5) == 1)
	SEQAN_TASSERT(degree(g, 0) == 2)
	SEQAN_TASSERT(getSuccessor(g, 0, "g") == 5)
	SEQAN_TASSERT(getSuccessor(g, 0, String<Dna>("ag")) == 3)  // The whole edge label or just the first letter
	SEQAN_TASSERT(getSuccessor(g, 0, "a") == _get_nil<TVertexDescriptor>())
	SEQAN_TASSERT(getPredecessor(g, 3, "a") == _get_nil<TVertexDescriptor>())
	SEQAN_TASSERT(getPredecessor(g, 3, "ag") == 0)
	SEQAN_TASSERT(getPredecessor(g, 5, 'g') == 0)
	addVertex(g);
	addVertex(g);
	addEdge(g,3,1,"aggg");
	addEdge(g,3,4,"gg");
	addEdge(g,5,2,"aggg");
	addEdge(g,5,7,"g");
	addEdge(g,7,6,"g");
	SEQAN_TASSERT(parseString(g, 0, "agaggg") == 1)
	SEQAN_TASSERT(parseString(g, 0, "aga") == 3)  // Does not reach 1
	SEQAN_TASSERT(parseString(g, 0, 'g') == 5)
	SEQAN_TASSERT(parseString(g, 0, "ggg") == 6)
	SEQAN_TASSERT(parseString(g, 0, "gaggg") == 2)
	SEQAN_TASSERT(parseString(g, 0, "gagggg") == 2)

	TWordGraph g_tmp(g);
	SEQAN_TASSERT(numVertices(g_tmp) == 8)
	SEQAN_TASSERT(parseString(g_tmp, 0, "agaggg") == 1)
	SEQAN_TASSERT(inDegree(g_tmp, 5) == 1)
	SEQAN_TASSERT(degree(g_tmp, 0) == 2)
	TWordGraph g_assign;
	g_assign = g;
	SEQAN_TASSERT(numVertices(g_assign) == 8)
	SEQAN_TASSERT(parseString(g_assign, 0, "agaggg") == 1)
	SEQAN_TASSERT(inDegree(g_assign, 5) == 1)
	SEQAN_TASSERT(degree(g_assign, 0) == 2)

	// Transpose
	transpose(g, g_tmp);
	SEQAN_TASSERT(numVertices(g_tmp) == 8)
	SEQAN_TASSERT(parseString(g_tmp, 2, "aggg") == 5)
	SEQAN_TASSERT(inDegree(g_tmp, 5) == 2)
	SEQAN_TASSERT(outDegree(g_tmp, 0) == 0)

	// Additional cargo
	typedef Graph<Automaton<Dna, short, WordGraph<> > > TWordGraph2;
	typedef VertexDescriptor<TWordGraph2>::Type TVertexDescriptor2;
	typedef EdgeDescriptor<TWordGraph2>::Type TEdgeDescriptor2;
	typedef Size<TWordGraph2>::Type TSize2;

	TWordGraph2 g2;
	TVertexDescriptor2 vert0 = addVertex(g2);
	TVertexDescriptor2 vert1 = addVertex(g2);
	TEdgeDescriptor2 ed0 = addEdge(g2,vert0,vert1,String<Dna>("acgt"),12);
	TEdgeDescriptor2 ed1 = addEdge(g2,vert1,vert0,"tgca",21);
	SEQAN_TASSERT(getCargo(g2.data_vertex[0].data_edge[(TSize2) Dna('a')])==12)
	SEQAN_TASSERT(getCargo(g2.data_vertex[1].data_edge[(TSize2) Dna('t')])==21)
	TWordGraph2 g2_transpose;
	transpose(g2,g2_transpose);
	SEQAN_TASSERT(getCargo(g2_transpose.data_vertex[1].data_edge[(TSize2) Dna('a')])==12)
	SEQAN_TASSERT(getCargo(g2_transpose.data_vertex[0].data_edge[(TSize2) Dna('t')])==21)
	clearEdges(g2_transpose);
	SEQAN_TASSERT(numEdges(g2_transpose) == 0)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Oracle() {
	Graph<Automaton<char> > g;
	createOracleOnReverse(g,"announce");
	SEQAN_TASSERT(parseString(g, 0, "e") == 1)
	SEQAN_TASSERT(parseString(g, 0, "ec") == 2)
	SEQAN_TASSERT(parseString(g, 0, "n") == 3)
	SEQAN_TASSERT(parseString(g, 0, "a") == 8)
	SEQAN_TASSERT(parseString(g, 0, "nn") == 7)

	Graph<Automaton<Dna> > g2;
	createOracle(g2,"ATATA");
	SEQAN_TASSERT(parseString(g2, 0, "A") == 1)
	SEQAN_TASSERT(parseString(g2, 0, "T") == 2)
	SEQAN_TASSERT(parseString(g2, 0, "AT") == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Graph() {
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

	typedef Graph<> StandardGraph;
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	
	StandardGraph g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)
	
	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	
	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v0);
	SEQAN_TASSERT(targetVertex(g, e1) == 0)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0)  //Expensive in standard graph!
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 1)
	SEQAN_TASSERT(degree(g, v0) == 2)
	
	// Add further edges and vertices
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1);
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 3)
	
	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
	addEdge(g,3,4);
	TEdgeDescriptor my_edge = addEdge(g,3,1);
	addEdge(g,3,0);
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 1)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 3)
	SEQAN_TASSERT(numEdges(g) == 5)
	SEQAN_TASSERT(outDegree(g, v3) == 3)	
	
	/*
	// Graph drawing
	// Raw output
	//std::cout << g << std::endl;

	// File output
	removeEdge(g,0,0); // ToDo: Self edges
	addEdge(g,4,3);  
	FILE * file = fopen(TEST_PATH "my_graph.tex", "w");
	write(file, g, PsTricks());
	fclose(file);
	removeEdge(g,4,3);
	addEdge(g,0,0);
*/

	// Remove edges
	removeEdge(g,my_edge);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 3)

	// Remove vertices 
	TEdgeDescriptor e3 = addEdge(g,3,3);
	addEdge(g,1,3);
	addEdge(g,0,3);
	addEdge(g,0,4);
	SEQAN_TASSERT(outDegree(g, 0) == 3)
	SEQAN_TASSERT(outDegree(g, 1) == 1)
	SEQAN_TASSERT(targetVertex(g, e3) == 3)
	SEQAN_TASSERT(sourceVertex(g, e3) == 3)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 0) == 2)
	SEQAN_TASSERT(outDegree(g, 1) == 0)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0);
	addEdge(g,4,1);
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	addEdge(g,4,2);
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 0)

	// Transpose
	StandardGraph g_tmp;
	//Transposes the graph in-place
	transpose(g, g_tmp); 
	g = g_tmp;
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 0)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardGraph g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
	addEdge(g_copy,3,0);
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 0)
	removeVertex(g,0);



	// Adjacency matrix
	Matrix<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	typedef Iterator<Matrix<TVertexDescriptor>,PositionIterator>::Type TMatrixIterator;
	TMatrixIterator p = begin(mat);
	TMatrixIterator endP = end(mat);
	for(;p != endP;goNext(p)) {
		if (*p == 1) {
			// Both edges go to vertex 4
			SEQAN_TASSERT(coordinate(p,0) == 4)
		}
	}


//____________________________________________________________________________
//Graph with edge cargo and edge ids
	typedef Pair<char, int> TPair;
	typedef EdgeList<TPair> TEdges;
	typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor2;
	typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor2;

	Graph<TEdges> g2;
	SEQAN_TASSERT(numVertices(g2) == 0)
	SEQAN_TASSERT(numEdges(g2) == 0)
	TVertexDescriptor2 ver0 = addVertex(g2);
	SEQAN_TASSERT(ver0 == 0)
	SEQAN_TASSERT(numVertices(g2) == 1)
	TVertexDescriptor2 ver1 = addVertex(g2);
	SEQAN_TASSERT(ver1 == 1)
	SEQAN_TASSERT(numVertices(g2) == 2)
	TEdgeDescriptor2 ed1 =addEdge(g2,ver0,ver0);
	TEdgeDescriptor2 ed2 =addEdge(g2,0,1);
	SEQAN_TASSERT(targetVertex(g2, ed1) == v0)
	SEQAN_TASSERT(targetVertex(g2, ed1) == 0)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == 0)
	SEQAN_TASSERT(targetVertex(g2, ed2) == 1)
	SEQAN_TASSERT(numEdges(g2) == 2)
	assignCargo(ed1, TPair('a',3));
	assignCargo(ed2, TPair('b',4));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT((getCargo(ed2)).i1 == 'b')
	SEQAN_TASSERT((getCargo(ed2)).i2 == 4)
	cargo(ed1) = TPair('c',1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'c')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 1)
	addVertex(g2);
	addVertex(g2);
	addVertex(g2);
	TEdgeDescriptor2 ed4 =addEdge(g2,1,4);
	cargo(ed4) = TPair('z',100);
	removeVertex(g2, 2);
	Graph<TEdges> g2_copy(g2);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 3)
	clearEdges(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 0)
	clearVertices(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	clear(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	SEQAN_TASSERT(numEdges(g2) == 3)
	SEQAN_TASSERT(outDegree(g2, 0) == 2)
	SEQAN_TASSERT(inDegree(g2, 0) == 1)
	transpose(g2, g2_copy);
	SEQAN_TASSERT(outDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(inDegree(g2_copy, 0) == 2)
	SEQAN_TASSERT(numEdges(g2_copy) == 3)
	TEdgeDescriptor2 edgCargo = addEdge(g2, 0, 0, TPair('m',3));
	SEQAN_TASSERT((getCargo(edgCargo)).i1 == 'm')
	SEQAN_TASSERT((getCargo(edgCargo)).i2 == 3)


//____________________________________________________________________________
//Graph without edge cargo and without edge ids
	typedef EdgeList<void, WithoutEdgeId> TEdges3;
	typedef VertexDescriptor<Graph<TEdges3> >::Type TVertexDescriptor3;
	typedef EdgeDescriptor<Graph<TEdges3> >::Type TEdgeDescriptor3;

	Graph<TEdges3> g3;
	addVertex(g3);addVertex(g3);addVertex(g3);
	addVertex(g3);addVertex(g3);
	addEdge(g3,1,4);
	SEQAN_TASSERT(numVertices(g3) == 5)
	SEQAN_TASSERT(numEdges(g3) == 1)
	TEdgeDescriptor3 edge3 = addEdge(g3,0,4);
	SEQAN_TASSERT(_getId(edge3) == 0);
	SEQAN_TASSERT(getCargo(edge3) == (void*) 0);
	addEdge(g3,0,2);
	addEdge(g3,0,0);
	removeEdge(g3,0,4);
	removeEdge(g3,0,2);
	SEQAN_TASSERT(numEdges(g3) == 2)

//____________________________________________________________________________
// Automaton
	typedef char CargoType;
	typedef EdgeList<CargoType> EdgeType;
	typedef VertexDescriptor<Graph<EdgeType> >::Type VertexDescriptorType;
	typedef EdgeDescriptor<Graph<EdgeType> >::Type EdgeDescriptorType;
	typedef Size<Graph<EdgeType> >::Type TSize;

	//Create the graph
	Graph<EdgeType> automaton;
	VertexDescriptorType rootVertex = addVertex(automaton); // A = 0
	addVertex(automaton); // B = 1
	addVertex(automaton); // C = 2
	addVertex(automaton); // D = 3
	addVertex(automaton); // E = 4
	addVertex(automaton); // F = 5
	addEdge(automaton,0,1,'2');
	addEdge(automaton,1,0,'1');
	addEdge(automaton,4,0,'6');
	addEdge(automaton,0,3,'7');
	addEdge(automaton,1,1,'3');
	addEdge(automaton,1,2,'4');
	addEdge(automaton,5,1,'8');
	addEdge(automaton,2,5,'5');
	addEdge(automaton,3,4,'2');
	addEdge(automaton,5,3,'7');

	InternalMap<char> eMap;
	initEdgeMap(automaton,eMap);
	VertexDescriptorType succ;
	succ = getSuccessor(automaton,rootVertex,eMap,'7');
	SEQAN_TASSERT(succ == 3)
	// Throws an error in debug mode because edge does not exist
	// succ = getSuccessor(automaton,rootVertex,eMap,'6');
	succ = getSuccessor(automaton,succ,eMap,'2');
	SEQAN_TASSERT(succ == 4)
	succ = getSuccessor(automaton,succ,eMap,'6');
	SEQAN_TASSERT(succ == 0)
	// If no map is specified it is assumed that an edge cargo exists!!!
	succ = getSuccessor(automaton,succ,'2');
	SEQAN_TASSERT(succ == 1)
	// Go backwards
	VertexDescriptorType pred;
	pred = getPredecessor(automaton,succ,eMap,'3');
	SEQAN_TASSERT(pred == 1)
	pred = getPredecessor(automaton,pred,eMap,'8');
	SEQAN_TASSERT(pred == 5)
	// If no map is specified it is assumed that an edge cargo exists!!!
	pred = getPredecessor(automaton,pred,'5');
	SEQAN_TASSERT(pred == 2)

	// Now using shortcuts
	succ = parseString(automaton,rootVertex,"7262");
	SEQAN_TASSERT(succ == 1)
	String<char> input("7262");
	succ = parseString(automaton,rootVertex, input);
	SEQAN_TASSERT(succ == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphExternalProperty() {
//____________________________________________________________________________
// Graph external property maps
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	
	Graph<> g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e1 =addEdge(g,v0,v0);
	TEdgeDescriptor e2 =addEdge(g,0,1);
	
	// Test external property maps
	String<int> dMap;
	initVertexMap(g,dMap);

	String<char> eMap;
	initEdgeMap(g,eMap);

	assignProperty(dMap, v0, 3);
	assignProperty(dMap, v1, 1);
	assignProperty(eMap, e1, 'a');
	assignProperty(eMap, e2, 'b');
	SEQAN_TASSERT(getProperty(dMap, v0) == 3)
	SEQAN_TASSERT(getProperty(dMap, v1) == 1)
	SEQAN_TASSERT(getProperty(eMap, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap, e1) == 'a')
	property(dMap, v1) = 2;
	property(eMap, e2) = 'c';
	SEQAN_TASSERT(getProperty(dMap, v1) == 2)
	SEQAN_TASSERT(getProperty(eMap, e2) == 'c')

	String<int> const dMap2(dMap);
	SEQAN_TASSERT(getProperty(dMap2, v0) == 3)
	SEQAN_TASSERT(getProperty(dMap2, v1) == 2)
	SEQAN_TASSERT(property(dMap2, v1) == 2)

	char names[] = {'r', 's'};
	String<char> nameMap;
	initVertexMap(g,nameMap, names);
	SEQAN_TASSERT(getProperty(nameMap, v0) == 'r')
	SEQAN_TASSERT(getProperty(nameMap, v1) == 's')

	std::string nameEd[] = {"ar", "ae"};
	String<std::string> edMap;
	initEdgeMap(g,edMap, nameEd);
	SEQAN_TASSERT(getProperty(edMap, e1) == "ar")
	SEQAN_TASSERT(getProperty(edMap, e2) == "ae")



	// Test prior resize before graph is created
	reserve(nameMap, 200);
	// Test prior resize before graph is created
	reserve(nameMap, 200);
}


//////////////////////////////////////////////////////////////////////////////

void Test_GraphInternalProperty() {
//____________________________________________________________________________
// Graph properties
	typedef Pair<char, int> TPair;
	typedef EdgeList<TPair> TEdges;
	typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor;

	// Create a simple graph
	Graph<TEdges> g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e1 =addEdge(g,v0,v0);
	TEdgeDescriptor e2 =addEdge(g,0,1);

	// First Variant: Explicit internal map with member Ids
	InternalMap<TPair, 1> eMap1; // This property map is used to access the first member
	InternalMap<TPair, 2> eMap2; // This property map is used to access the second member
	initEdgeMap(g,eMap1);
	initEdgeMap(g,eMap2);
	assignProperty(eMap1, e1, 'a');
	assignProperty(eMap2, e1, 20);
	assignProperty(eMap1, e2, 'b');
	assignProperty(eMap2, e2, 50);
	SEQAN_TASSERT(getProperty(eMap1, e1) == 'a')
	SEQAN_TASSERT(getProperty(eMap2, e1) == 20)
	SEQAN_TASSERT(getProperty(eMap1, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap2, e2) == 50)
	// Note: That these properties are stored inside the cargo of each edge
	SEQAN_TASSERT(getCargo(e1).i1 == 'a')
	SEQAN_TASSERT(getCargo(e1).i2 == 20)
	SEQAN_TASSERT(getCargo(e2).i1 == 'b')
	SEQAN_TASSERT(getCargo(e2).i2 == 50)
	assignProperty(eMap1, e1, 'c');
	assignProperty(eMap2, e1, 10);
	SEQAN_TASSERT(getProperty(eMap1, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap2, e1) == 10)
	SEQAN_TASSERT(property(eMap1, e1) == 'c')
	SEQAN_TASSERT(property(eMap2, e1) == 10)
	InternalMap<TPair, 1> const eMap3(eMap1);
	InternalMap<TPair, 2> const eMap31(eMap2);
	SEQAN_TASSERT(getProperty(eMap3, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap3, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap31, e1) == 10)
	SEQAN_TASSERT(property(eMap31, e1) == 10)
	SEQAN_TASSERT(property(eMap3, e2) == 'b')
	// Create a simple graph with unsigned int cargo
	typedef EdgeDescriptor<Graph<EdgeList<unsigned int> > >::Type TEdgeDescriptor2;
	Graph<EdgeList<unsigned int> > g2;
	addVertex(g2);
	addVertex(g2);
	TEdgeDescriptor2 edge1 =addEdge(g2,v0,v0);
	TEdgeDescriptor2 edge2 =addEdge(g2,0,1);
	InternalMap<unsigned int> edgeMap;
	initEdgeMap(g2,edgeMap);
	assignProperty(edgeMap, edge1 ,3);
	SEQAN_TASSERT(getProperty(edgeMap, edge1) == 3)
	SEQAN_TASSERT(property(edgeMap, edge1) == 3)
	InternalMap<unsigned int> const edgeMap2(edgeMap);
	SEQAN_TASSERT(getProperty(edgeMap2, edge1) == 3)
	SEQAN_TASSERT(property(edgeMap2, edge1) == 3)

	// Test dummy call
	reserve(edgeMap, 200);
	

	// Second Variant: Pointer to member using a class
	InternalPointerMap<char TPair:: *, &TPair::i1> eMap4;
	InternalPointerMap<int TPair:: *, &TPair::i2> eMap5;
	initEdgeMap(g,eMap4);
	initEdgeMap(g,eMap5);
	assignProperty(eMap4, e1, 'c');
	assignProperty(eMap5, e1, 10);
	assignProperty(eMap4, e2, 'd');
	assignProperty(eMap5, e2, 30);
	SEQAN_TASSERT(getProperty(eMap4, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap5, e1) == 10)
	SEQAN_TASSERT(getProperty(eMap4, e2) == 'd')
	SEQAN_TASSERT(getProperty(eMap5, e2) == 30)
	property(eMap4,e1)='z';
	property(eMap5,e1)=100;
	SEQAN_TASSERT(getProperty(eMap4, e1) == 'z')
	SEQAN_TASSERT(getProperty(eMap5, e1) == 100)
	InternalPointerMap<char TPair:: *, &TPair::i1> const eMap6(eMap4);
	SEQAN_TASSERT(getProperty(eMap6, e1) == 'z')
	SEQAN_TASSERT(getProperty(eMap6, e2) == 'd')
	SEQAN_TASSERT(property(eMap6, e2) == 'd')
	
	// Test dummy call
	reserve(eMap4, 200);

	// Third Variant: Raw pointer to member
	char TPair:: * pseudo_map = &TPair::i1;
	assignProperty(pseudo_map, e1, 'z');
	assignProperty(pseudo_map, e2, 'w');
	SEQAN_TASSERT(getProperty(pseudo_map, e1) == 'z')
	SEQAN_TASSERT(getProperty(pseudo_map, e2) == 'w')
	property(pseudo_map,e1)='k';
	SEQAN_TASSERT(getProperty(pseudo_map, e1) == 'k')

	// Test dummy call
	reserve(pseudo_map, 200);
}



void Test_GraphStack() {
//____________________________________________________________________________
// Graph stack

    String<int, BlockString<2> > stack;
	SEQAN_TASSERT(empty(stack) == true)
	push_back(stack,1);
	SEQAN_TASSERT(top(stack)==1)
	SEQAN_TASSERT(stack[0]==1)
	push_back(stack,20);
	SEQAN_TASSERT(top(stack)==20)
	SEQAN_TASSERT(stack[0]==1)
	SEQAN_TASSERT(stack[1]==20)
	push_back(stack,3);
	SEQAN_TASSERT(top(stack)==3)
	SEQAN_TASSERT(stack[2]==3)
	pop_back(stack);
	SEQAN_TASSERT(top(stack)==20)
	SEQAN_TASSERT(empty(stack) == false)
	pop_back(stack);
	pop_back(stack);
	pop_back(stack);
	pop_back(stack);
	SEQAN_TASSERT(empty(stack) == true)
	push_back(stack,20);
	SEQAN_TASSERT(top(stack)==20)
	String<char, BlockString<2> > wS("hallo");
	SEQAN_TASSERT(top(wS)=='o')
	SEQAN_TASSERT(wS[2]=='l')
	wS = "new_stack";
	SEQAN_TASSERT(top(wS)=='k')
	SEQAN_TASSERT(wS[2]=='w')
	SEQAN_TASSERT(length(wS)==9)
	typedef Iterator<String<char,BlockString<2> >, Standard>::Type TIter;
	TIter it = begin(wS, Standard());
	TIter itEnd = end(wS, Standard());
	SEQAN_TASSERT(*it=='n')
	++it;
	SEQAN_TASSERT(*it=='e')
	goNext(it);
	SEQAN_TASSERT(*it=='w')
	SEQAN_TASSERT(it != itEnd)
	SEQAN_TASSERT(atEnd(it, wS)==false)
	clear(wS);
	SEQAN_TASSERT(empty(wS) == true)
	append(wS, 'A');
	SEQAN_TASSERT(top(wS)=='A')
	append(wS, 'H');
	SEQAN_TASSERT(wS[1]=='H')
	appendValue(wS, 'A');
	SEQAN_TASSERT(wS[2]=='A')
	int a = 20;
	String<int, BlockString<2> > test(a);
	SEQAN_TASSERT(test[0]==20)
	String<int, BlockString<2> > test2(test);
	SEQAN_TASSERT(test2[0]==20)
	test[0]=10;
	SEQAN_TASSERT(test[0]==10)
	SEQAN_TASSERT(test2[0]==20)
	test2 = test;
	SEQAN_TASSERT(test[0]==10)
	SEQAN_TASSERT(test2[0]==10)
	String<char, BlockString<2> > const test3("Hallo");
	SEQAN_TASSERT(test3[0]=='H')
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphVertexIterator() {
//____________________________________________________________________________
// Graph VertexIterator
	// Works on all graph types
	//typedef Graph<Automaton<Dna> > TGraph;
	typedef Graph<EdgeList<char> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TEdgeDescriptor e1 =addEdge(g,v0,v0,'c');
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1,'t');
	addVertex(g); //2
	addVertex(g); //3
	addVertex(g); //4

	//Tricky case -> id 0 is released
	removeVertex(g,v0);
	removeVertex(g,3);
	typedef Iterator<TGraph, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	SEQAN_TASSERT(atBegin(it)==true)
	SEQAN_TASSERT(getValue(it)==1)
	SEQAN_TASSERT(value(it)==1)
	SEQAN_TASSERT(getValue(it)==1)
	goNext(it);
	SEQAN_TASSERT(atBegin(it)==false)
	SEQAN_TASSERT(getValue(it)==2)
	++it;
	SEQAN_TASSERT(getValue(it)==4)
	SEQAN_TASSERT(atEnd(it)==false)
	it++; 
	SEQAN_TASSERT(atEnd(it)==true)
	++it; 
	SEQAN_TASSERT(atEnd(it)==true)
	goPrevious(it);
	// No assignment to vertex iterators
	//*it = 3;
	SEQAN_TASSERT((*it)==4)
	SEQAN_TASSERT(atEnd(it)==false)
	it--;
	SEQAN_TASSERT(getValue(it)==2)
	SEQAN_TASSERT(atBegin(it)==false)
	--it;
	SEQAN_TASSERT(atBegin(it)==true)
	SEQAN_TASSERT(getValue(it)==1)
	--it;
	SEQAN_TASSERT(atBegin(it)==true)
	SEQAN_TASSERT(getValue(it)==1)
	TVertexIterator it2(g);
	TVertexIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
	SEQAN_TASSERT(&hostGraph(it) == &g)
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphOutEdgeIterator() {
//____________________________________________________________________________
// Graph EdgeIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	
	Graph<> g;
	TVertexDescriptor v0 = addVertex(g);
	TEdgeDescriptor e1 =addEdge(g,v0,v0);
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1);
	addVertex(g); //2
	addVertex(g); //3
	addVertex(g); //4

	typedef Iterator<Graph<>, OutEdgeIterator<> >::Type TOutEdgeIterator;
	TOutEdgeIterator it(g, v0);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==1)
	SEQAN_TASSERT(sourceVertex(g, value(it))==0)
	SEQAN_TASSERT(targetVertex(g, *it)==1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==0)
	// Fast
	SEQAN_TASSERT(sourceVertex(it)==0)
	SEQAN_TASSERT(targetVertex(it)==0)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	++it;
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goPrevious(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==0)
	--it;
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==1)
	// Fast
	SEQAN_TASSERT(sourceVertex(it)==0)
	SEQAN_TASSERT(targetVertex(it)==1)
	it++; 
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==0)
	it--;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==1)
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
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphEdgeIterator() {
//____________________________________________________________________________
// Graph EdgeIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	
	Graph<> g;
	TVertexDescriptor v0 = addVertex(g);
	TEdgeDescriptor e1 =addEdge(g,v0,v0);
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1);
	addVertex(g); //2
	addVertex(g); //3
	addVertex(g); //4
	addVertex(g); //5
	addEdge(g,3,2);
	addEdge(g,3,4);
	addEdge(g,4,4);

	typedef Iterator<Graph<>, EdgeIterator<> >::Type TEdgeIterator;
	TEdgeIterator it(g);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==1)
	SEQAN_TASSERT(sourceVertex(g, value(it))==0)
	SEQAN_TASSERT(targetVertex(g, *it)==1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==0)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==0)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	++it;
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==3)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==4)
	// Fast
	SEQAN_TASSERT(sourceVertex(it)==3)
	SEQAN_TASSERT(targetVertex(it)==4)
	it++;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==3)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	it++;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==4)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==4)
	++it;
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	++it;
	++it;
	goPrevious(it);
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==4)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==4)
	--it;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==3)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==2)
	it--;
	SEQAN_TASSERT(sourceVertex(g, getValue(it))==3)
	SEQAN_TASSERT(targetVertex(g, getValue(it))==4)
	TEdgeIterator it2(g);
	TEdgeIterator it3;
	goBegin(it);
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
	SEQAN_TASSERT(&hostGraph(it) == &g)
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphAdjacencyIterator() {
//____________________________________________________________________________
// Graph AdjacencyIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	
	Graph<> g;
	TVertexDescriptor v0 = addVertex(g);
	TEdgeDescriptor e1 =addEdge(g,v0,v0);
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1);
	addVertex(g); //2
	addVertex(g); //3
	addVertex(g); //4
	addVertex(g); //5
	addEdge(g,3,2);
	addEdge(g,3,4);
	addEdge(g,4,4);

	typedef Iterator<Graph<>, AdjacencyIterator<> >::Type TAdjacencyIterator;
	TAdjacencyIterator it(g,3);
	SEQAN_TASSERT(getValue(it) == 4)
	SEQAN_TASSERT(&hostGraph(it) == &g)
	SEQAN_TASSERT(value(it) == 4)
	SEQAN_TASSERT(*it == 4)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==2)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	++it;
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goBegin(it);
	it++;
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	it++;
	goPrevious(it);
	SEQAN_TASSERT(getValue(it)==2)
	--it;
	SEQAN_TASSERT(getValue(it)==4)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goEnd(it);
	it--;
	SEQAN_TASSERT(getValue(it)==2)
	goBegin(it);
	TAdjacencyIterator it2(g, 3);
	TAdjacencyIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphBfsIterator() {
//____________________________________________________________________________
// Graph BfsIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Create the graph
	Graph<> g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,1,0);
	addEdge(g,1,5);
	addEdge(g,0,4);
	addEdge(g,5,2);
	addEdge(g,5,6);
	addEdge(g,2,6);
	addEdge(g,2,3);
	addEdge(g,6,3);
	addEdge(g,6,7);
	addEdge(g,3,7);

	typedef Iterator<Graph<>, BfsIterator<> >::Type TBfsIterator;
	TBfsIterator it(g,1);
	SEQAN_TASSERT(getValue(it) == 1)
	SEQAN_TASSERT(&hostGraph(it) == &g)
	SEQAN_TASSERT(value(it) == 1)
	SEQAN_TASSERT(*it == 1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==5)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==0)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==6)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==2)
	++it;
	SEQAN_TASSERT(getValue(it)==4)
	it++;
	SEQAN_TASSERT(getValue(it)==7)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==3)
	goNext(it);
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goBegin(it);
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	TBfsIterator it2(g, 1);
	TBfsIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
}

//////////////////////////////////////////////////////////////////////////////
void Test_GraphDfsIterator() {
//____________________________________________________________________________
// Graph DfsIterator
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Create the graph
	Graph<> g;
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);addVertex(g);addVertex(g);
	addEdge(g,1,0);
	addEdge(g,1,5);
	addEdge(g,0,4);
	addEdge(g,5,2);
	addEdge(g,5,6);
	addEdge(g,2,6);
	addEdge(g,2,3);
	addEdge(g,6,3);
	addEdge(g,6,7);
	addEdge(g,3,7);

	typedef Iterator<Graph<>, DfsIterator<> >::Type TDfsIterator;
	TDfsIterator it(g,1);
	SEQAN_TASSERT(getValue(it) == 1)
	SEQAN_TASSERT(&hostGraph(it) == &g)
	SEQAN_TASSERT(value(it) == 1)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==0)
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==false)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==4)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==5)
	SEQAN_TASSERT(value(it)==5)
	SEQAN_TASSERT(*it == 5)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==2)
	++it;
	SEQAN_TASSERT(getValue(it)==3)
	it++;
	SEQAN_TASSERT(getValue(it)==7)
	goNext(it);
	SEQAN_TASSERT(getValue(it)==6)
	goNext(it);
	SEQAN_TASSERT(atEnd(it)==true)
	SEQAN_TASSERT(atBegin(it)==false)
	goBegin(it);
	SEQAN_TASSERT(atEnd(it)==false)
	SEQAN_TASSERT(atBegin(it)==true)
	TDfsIterator it2(g, 1);
	TDfsIterator it3;
	it3 = it;
	SEQAN_TASSERT(it == it2)
	SEQAN_TASSERT(it2 == it3)
	goEnd(it);
	SEQAN_TASSERT(it2 != it)
	goEnd(it2);
	SEQAN_TASSERT(it2 == it)
	goBegin(it2);
	SEQAN_TASSERT(it2 != it)
}

//////////////////////////////////////////////////////////////////////////////
void Test_BreadthFirstSearch() {
//____________________________________________________________________________
// Breadth-First Search
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 20;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 1,0, 0,4, 4,0, 1,5, 5,1, 2,5, 5,2, 2,6, 6,2, 2,3, 3,2, 3,6, 6,3, 3,7, 7,3, 5,6, 6,5, 6,7, 7,6};
	
	//Create the graph
	Graph<> g(edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bfs
	breadth_first_search(g, 1, predMap, distMap);
	
	SEQAN_TASSERT(getProperty(distMap, 0) == 1)
	SEQAN_TASSERT(getProperty(distMap, 1) == 0)
	SEQAN_TASSERT(getProperty(distMap, 2) == 2)
	SEQAN_TASSERT(getProperty(distMap, 3) == 3)
	SEQAN_TASSERT(getProperty(distMap, 4) == 2)
	SEQAN_TASSERT(getProperty(distMap, 5) == 1)
	SEQAN_TASSERT(getProperty(distMap, 6) == 2)
	SEQAN_TASSERT(getProperty(distMap, 7) == 3)
	SEQAN_TASSERT(getProperty(predMap, 0) == 1)
	SEQAN_TASSERT(getProperty(predMap, 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 2) == 5)
	SEQAN_TASSERT(getProperty(predMap, 3) == 6)
	SEQAN_TASSERT(getProperty(predMap, 4) == 0)
	SEQAN_TASSERT(getProperty(predMap, 5) == 1)
	SEQAN_TASSERT(getProperty(predMap, 6) == 5)
	SEQAN_TASSERT(getProperty(predMap, 7) == 6)
}

//////////////////////////////////////////////////////////////////////////////

void Test_DepthFirstSearch() {
//____________________________________________________________________________
// Depth-First Search
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 8;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,4, 2,4, 2,5, 3,1, 4,3, 5,5};

	//Create the graph
	Graph<> g(edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;

	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	SEQAN_TASSERT(getProperty(discoveryTimeMap, 0) == 1)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 1) == 2)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 2) == 9)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 3) == 4)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 4) == 3)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 5) == 10)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 0) == 8)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 1) == 7)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 2) == 12)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 3) == 5)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 4) == 6)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 5) == 11)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 0)
	SEQAN_TASSERT(getProperty(predMap, 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 3) == 4)
	SEQAN_TASSERT(getProperty(predMap, 4) == 1)
	SEQAN_TASSERT(getProperty(predMap, 5) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_TopologicalSort() {
//____________________________________________________________________________
// Topological Sort
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,2, 3,2, 5,7, 5,6, 6,7, 6,3, 8,7};
	
	//Create the graph
	Graph<> g(edges, numEdges);

	// Predecessor and distance map
	String<TVertexDescriptor> order;
	
	// Topological sort
	topological_sort(g, order);

	SEQAN_TASSERT(getValue(order, 0) == 8)
	SEQAN_TASSERT(getValue(order, 1) == 5)
	SEQAN_TASSERT(getValue(order, 2) == 6)
	SEQAN_TASSERT(getValue(order, 3) == 7)
	SEQAN_TASSERT(getValue(order, 4) == 4)
	SEQAN_TASSERT(getValue(order, 5) == 0)
	SEQAN_TASSERT(getValue(order, 6) == 3)
	SEQAN_TASSERT(getValue(order, 7) == 1)
	SEQAN_TASSERT(getValue(order, 8) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_StronglyConnectedComponents() {
//____________________________________________________________________________
// Strongly-Connected-Components

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};

	//Create the graph
	Graph<> g(edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> component;
	
	// Strongly Connected Components
	strongly_connected_components(g, component);

	SEQAN_TASSERT(getValue(component, 0) == 3)
	SEQAN_TASSERT(getValue(component, 1) == 3)
	SEQAN_TASSERT(getValue(component, 2) == 2)
	SEQAN_TASSERT(getValue(component, 3) == 2)
	SEQAN_TASSERT(getValue(component, 4) == 3)
	SEQAN_TASSERT(getValue(component, 5) == 1)
	SEQAN_TASSERT(getValue(component, 6) == 1)
	SEQAN_TASSERT(getValue(component, 7) == 0)
}


//////////////////////////////////////////////////////////////////////////////

void Test_DagShortestPath() {
//____________________________________________________________________________
// DAG-Shortest Paths
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,2, 0,1, 1,3, 1,2, 2,5, 2,4, 2,3, 3,5, 3,4, 4,5};
	int weights[] =             {3,   5,   6,   2,   2,   4,   7,   1,   -1,  -2};
	
	//Create the graph
	Graph<> g(edges, numEdges);
	String<int> weightMap;
	initEdgeMap(g,weightMap, weights);

	// Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// DAG-Shortest path(Graph, sourceVertex_vertex, weightMap, predMap, distMap)
	dag_shortest_path(g,1,weightMap,predMap,distMap);
	
	SEQAN_TASSERT(getProperty(distMap, 0) == getInfinityDistance(weightMap))
	SEQAN_TASSERT(getProperty(distMap, 1) == 0)
	SEQAN_TASSERT(getProperty(distMap, 2) == 2)
	SEQAN_TASSERT(getProperty(distMap, 3) == 6)
	SEQAN_TASSERT(getProperty(distMap, 4) == 5)
	SEQAN_TASSERT(getProperty(distMap, 5) == 3)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 1)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
	SEQAN_TASSERT(getProperty(predMap, 5) == 4)
}

//////////////////////////////////////////////////////////////////////////////

void Test_BellmanFord() {
//____________________________________________________________________________
// Bellman-Ford

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g(edges, numEdges);

	// In parameters of Bellman-Ford
	// Define the distances or weights of the edges
	String<unsigned int> weightMap;
	initEdgeMap(g,weightMap, weights);

	// Out parameters of Bellman-Ford: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bellman-Ford
	bool noNegativeCycle = bellman_ford_algorithm(g,0,weightMap,predMap,distMap);

	SEQAN_TASSERT(getProperty(distMap, 0) == 0)
	SEQAN_TASSERT(getProperty(distMap, 1) == 8)
	SEQAN_TASSERT(getProperty(distMap, 2) == 9)
	SEQAN_TASSERT(getProperty(distMap, 3) == 5)
	SEQAN_TASSERT(getProperty(distMap, 4) == 7)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 3)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 0)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
	SEQAN_TASSERT(noNegativeCycle == true)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Dijkstra() {
//____________________________________________________________________________
// Dijkstra

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g(edges, numEdges);

	// In parameters of Dijkstra
	// Define the distances or weights of the edges
	String<unsigned int> weightMap;
	initEdgeMap(g,weightMap, weights);

	// Out parameters of Dijkstra: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Dijkstra
	dijkstra(g,0,weightMap,predMap,distMap);

	SEQAN_TASSERT(getProperty(distMap, 0) == 0)
	SEQAN_TASSERT(getProperty(distMap, 1) == 8)
	SEQAN_TASSERT(getProperty(distMap, 2) == 9)
	SEQAN_TASSERT(getProperty(distMap, 3) == 5)
	SEQAN_TASSERT(getProperty(distMap, 4) == 7)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getProperty(predMap, 1) == 3)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 0)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
}

//////////////////////////////////////////////////////////////////////////////

void Test_AllPairsShortestPath() {
//____________________________________________________________________________
// All-Pairs Shortest Path

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	Graph<> g(edges, numEdges);

	// In parameters
	String<int> weightMap;
	initEdgeMap(g,weightMap, weights);
	// Out parameter
	Matrix<int> distMat;
	Matrix<TVertexDescriptor> predMat;

	// All-Pairs shortest path
	all_pairs_shortest_path(g,weightMap, distMat, predMat);

	unsigned int len = length(distMat, 0);
	SEQAN_TASSERT(getValue(distMat, 0*len + 0) == 0)
	SEQAN_TASSERT(getValue(distMat, 0*len + 1) == 1)
	SEQAN_TASSERT(getValue(distMat, 0*len + 2) == -3)
	SEQAN_TASSERT(getValue(distMat, 0*len + 3) == 2)
	SEQAN_TASSERT(getValue(distMat, 0*len + 4) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(distMat, 1*len + 1) == 0)
	SEQAN_TASSERT(getValue(distMat, 1*len + 2) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(distMat, 1*len + 4) == -1)
	SEQAN_TASSERT(getValue(distMat, 2*len + 0) == 7)
	SEQAN_TASSERT(getValue(distMat, 2*len + 1) == 4)
	SEQAN_TASSERT(getValue(distMat, 2*len + 2) == 0)
	SEQAN_TASSERT(getValue(distMat, 2*len + 3) == 5)
	SEQAN_TASSERT(getValue(distMat, 2*len + 4) == 3)
	SEQAN_TASSERT(getValue(distMat, 3*len + 0) == 2)
	SEQAN_TASSERT(getValue(distMat, 3*len + 1) == -1)
	SEQAN_TASSERT(getValue(distMat, 3*len + 2) == -5)
	SEQAN_TASSERT(getValue(distMat, 3*len + 3) == 0)
	SEQAN_TASSERT(getValue(distMat, 3*len + 4) == -2)
	SEQAN_TASSERT(getValue(distMat, 4*len + 0) == 8)
	SEQAN_TASSERT(getValue(distMat, 4*len + 1) == 5)
	SEQAN_TASSERT(getValue(distMat, 4*len + 2) == 1)
	SEQAN_TASSERT(getValue(distMat, 4*len + 3) == 6)
	SEQAN_TASSERT(getValue(distMat, 4*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 0*len + 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 0*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 0*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 0*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 0*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 1*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 1*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 2*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 2*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 2*len + 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 2*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 3*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 3*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 3) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 3*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 4*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 4*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 4*len + 4) == getNilPredecessor(g))
}

//////////////////////////////////////////////////////////////////////////////

void Test_FloydWarshall() {
//____________________________________________________________________________
// Floyd Warshall

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	Graph<> g(edges, numEdges);

	// In parameters
	String<int> weightMap;
	initEdgeMap(g,weightMap, weights);
	// Out parameter
	Matrix<int> distMat;
	Matrix<TVertexDescriptor> predMat;

	// Floyd-Warshall
	floyd_warshall(g,weightMap, distMat, predMat);

	unsigned int len = length(distMat, 0);
	SEQAN_TASSERT(getValue(distMat, 0*len + 0) == 0)
	SEQAN_TASSERT(getValue(distMat, 0*len + 1) == 1)
	SEQAN_TASSERT(getValue(distMat, 0*len + 2) == -3)
	SEQAN_TASSERT(getValue(distMat, 0*len + 3) == 2)
	SEQAN_TASSERT(getValue(distMat, 0*len + 4) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(distMat, 1*len + 1) == 0)
	SEQAN_TASSERT(getValue(distMat, 1*len + 2) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(distMat, 1*len + 4) == -1)
	SEQAN_TASSERT(getValue(distMat, 2*len + 0) == 7)
	SEQAN_TASSERT(getValue(distMat, 2*len + 1) == 4)
	SEQAN_TASSERT(getValue(distMat, 2*len + 2) == 0)
	SEQAN_TASSERT(getValue(distMat, 2*len + 3) == 5)
	SEQAN_TASSERT(getValue(distMat, 2*len + 4) == 3)
	SEQAN_TASSERT(getValue(distMat, 3*len + 0) == 2)
	SEQAN_TASSERT(getValue(distMat, 3*len + 1) == -1)
	SEQAN_TASSERT(getValue(distMat, 3*len + 2) == -5)
	SEQAN_TASSERT(getValue(distMat, 3*len + 3) == 0)
	SEQAN_TASSERT(getValue(distMat, 3*len + 4) == -2)
	SEQAN_TASSERT(getValue(distMat, 4*len + 0) == 8)
	SEQAN_TASSERT(getValue(distMat, 4*len + 1) == 5)
	SEQAN_TASSERT(getValue(distMat, 4*len + 2) == 1)
	SEQAN_TASSERT(getValue(distMat, 4*len + 3) == 6)
	SEQAN_TASSERT(getValue(distMat, 4*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 0*len + 0) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 0*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 0*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 0*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 0*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 1) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 1*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 1*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 2*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 2*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 2*len + 2) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 2*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 3*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 3*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 3) == getNilPredecessor(g))
	SEQAN_TASSERT(getValue(predMat, 3*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 4*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 4*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 4*len + 4) == getNilPredecessor(g))
}


//////////////////////////////////////////////////////////////////////////////

void Test_TransitiveClosure() {
//____________________________________________________________________________
// Transitive Closure

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 5;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {3,0, 1,2, 2,1, 1,3, 3,2};
	
	//Create the graph
	Graph<> g(edges, numEdges);

	// Transitive-Closure
	Matrix<bool> closure;
	transitive_closure(g,closure);
	
	unsigned int len = length(closure, 0);
	SEQAN_TASSERT(getValue(closure, 0*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 0*len + 1) == 0)
	SEQAN_TASSERT(getValue(closure, 0*len + 2) == 0)
	SEQAN_TASSERT(getValue(closure, 0*len + 3) == 0)
	SEQAN_TASSERT(getValue(closure, 1*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 3) == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphAlgorithms() {
//____________________________________________________________________________
// Graph Algorithms
	// Elementary graph algorithms
	Test_BreadthFirstSearch();
	Test_DepthFirstSearch();
	Test_TopologicalSort();
	Test_StronglyConnectedComponents();

	// Single-Source shortest paths
	Test_DagShortestPath();
	Test_BellmanFord();
	Test_Dijkstra();

	// All-Pairs Shortest paths
	Test_AllPairsShortestPath();
	Test_FloydWarshall();
	Test_TransitiveClosure();


	//Todo
	//Spanning Trees
	//Maximum Flow
	//Matching
}


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_IdManager();
	Test_EdgeStump();
	Test_EdgeAutomaton();
	Test_Automaton();
	Test_WordGraph();
	Test_Oracle();
	Test_Graph();
	Test_GraphExternalProperty();
	Test_GraphInternalProperty();
	Test_GraphStack();
	Test_GraphVertexIterator();
	Test_GraphOutEdgeIterator();
	Test_GraphEdgeIterator();
	Test_GraphAdjacencyIterator();
	Test_GraphBfsIterator();
	Test_GraphDfsIterator();
	Test_GraphAlgorithms();

//____________________________________________________________________________

	debug::verifyCheckpoints("projects/library/seqan/graph/graph_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_idmanager.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_property.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_iterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_vertexiterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_outedgeiterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_adjacencyiterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_edgeiterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_bfsiterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_dfsiterator.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_edgestump.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_edgelist.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_edgeautomaton.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_automaton.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_wordgraph.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_impl_oracle.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_algorithm.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_stack.h");
	
	SEQAN_TREPORT("TEST END")

	return 0;
}
