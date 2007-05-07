#ifndef SEQAN_HEADER_TEST_GRAPH_TYPES_H
#define SEQAN_HEADER_TEST_GRAPH_TYPES_H


using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

void Test_Directed() {
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

	typedef Graph<> StandardGraph;
	typedef VertexDescriptor<StandardGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardGraph>::Type TEdgeDescriptor;
	
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
	SEQAN_TASSERT(findEdge(g, v0, v0) == e1)
	SEQAN_TASSERT(_getVertexString(g)[0] == e1)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 1)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
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
	
	// Graph drawing
	removeEdge(g,0,0); // ToDo: Drawing of self edges
	addEdge(g,4,3);
	// Raw output
	// std::cout << g << ::std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_graph.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	StandardGraph gTmp;
	strm.open(TEST_PATH "my_graph.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	removeEdge(g,4,3);
	addEdge(g,0,0);

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
	transpose(g); 
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
	unsigned int len = length(mat, 0);
	SEQAN_TASSERT(getValue(mat, 1*len+4) == 1)
	SEQAN_TASSERT(getValue(mat, 2*len+4) == 1)
	SEQAN_TASSERT(getValue(mat, 2*len+2) == 0)

//____________________________________________________________________________
//Graph with edge cargo and edge ids
	typedef Pair<char, int> TPair;
	typedef Directed<TPair> TEdges;
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
	TEdgeDescriptor2 ed1 =addEdge(g2,ver0,ver0, TPair('a',3));
	TEdgeDescriptor2 ed2 =addEdge(g2,0,1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT(targetVertex(g2, ed1) == v0)
	SEQAN_TASSERT(targetVertex(g2, ed1) == 0)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == 0)
	SEQAN_TASSERT(targetVertex(g2, ed2) == 1)
	SEQAN_TASSERT(numEdges(g2) == 2)
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
	typedef Directed<void, WithoutEdgeId> TEdges3;
	typedef VertexDescriptor<Graph<TEdges3> >::Type TVertexDescriptor3;
	typedef EdgeDescriptor<Graph<TEdges3> >::Type TEdgeDescriptor3;

	Graph<TEdges3> g3;
	addVertex(g3);addVertex(g3);addVertex(g3);
	addVertex(g3);addVertex(g3);
	addEdge(g3,1,4);
	SEQAN_TASSERT(numVertices(g3) == 5)
	SEQAN_TASSERT(numEdges(g3) == 1)
	TEdgeDescriptor3 edge3 = addEdge(g3,0,4);
	//SEQAN_TASSERT(_getId(edge3) == 0);
	SEQAN_TASSERT(getCargo(edge3) == (void*) 0);
	addEdge(g3,0,2);
	addEdge(g3,0,0);
	removeEdge(g3,0,4);
	removeEdge(g3,0,2);
	SEQAN_TASSERT(numEdges(g3) == 2)

		
	clear(g);
	TVertexDescriptor edges[] = {0,1, 1,2};
	unsigned int numEdg = 2;
	std::string nameEd[] = {"ar", "ae"};
	addEdges(g, edges, numEdg);
	String<std::string> edMap;
	resizeEdgeMap(g, edMap, nameEd);
	SEQAN_TASSERT(getProperty(edMap, 0) == "ar")
	SEQAN_TASSERT(getProperty(edMap, 1) == "ae")
}

//////////////////////////////////////////////////////////////////////////////

void Test_Undirected() {
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

	typedef Graph<Undirected<void> > StandardGraph;
	typedef VertexDescriptor<StandardGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardGraph>::Type TEdgeDescriptor;
	
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
	// TEdgeDescriptor e1 =addEdge(g,v0,v0);  // Self edges are not allowed in undirected graphs
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e =addEdge(g,0,1);
	SEQAN_TASSERT(findEdge(g, 0, 1) == e)
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

	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
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

	// Graph drawing
	// Raw output
	// std::cout << g << ::std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_undirected_graph.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	StandardGraph gTmp;
	strm.open(TEST_PATH "my_undirected_graph.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	// Remove edges
	removeEdge(g,my_edge);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 2)

	
	// Remove vertices 
	addVertex(g);  //5
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
	SEQAN_TASSERT(inDegree(g, 4) == 2)

	// Transpose
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardGraph g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
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
	Matrix<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = getIdUpperBound(g.data_id_managerV);
	SEQAN_TASSERT(getValue(mat,0*len+2) == 1)
	SEQAN_TASSERT(getValue(mat,3*len+2) == 0)
	SEQAN_TASSERT(getValue(mat,0*len+2) == getValue(mat,2*len+0))
	SEQAN_TASSERT(getValue(mat,1*len+4) == getValue(mat,4*len+1))
	SEQAN_TASSERT(getValue(mat,2*len+4) == getValue(mat,4*len+2))

//____________________________________________________________________________
//Graph with edge cargo and edge ids
	typedef Pair<char, int> TPair;
	typedef Undirected<TPair> TEdges;
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
	TEdgeDescriptor2 ed1 =addEdge(g2,0,1);
	SEQAN_TASSERT(targetVertex(g2, ed1) == 1)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == 0)
	SEQAN_TASSERT(numEdges(g2) == 1)
	assignCargo(ed1, TPair('a',3));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
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
	SEQAN_TASSERT(numEdges(g2_copy) == 2)
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
	transpose(g2, g2_copy);
	SEQAN_TASSERT(outDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(inDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(numEdges(g2_copy) == 2)
	TEdgeDescriptor2 edgCargo = addEdge(g2, 0, 3, TPair('m',3));
	SEQAN_TASSERT((getCargo(edgCargo)).i1 == 'm')
	SEQAN_TASSERT((getCargo(edgCargo)).i2 == 3)

//____________________________________________________________________________
//Graph without edge cargo and without edge ids
	typedef Undirected<void, WithoutEdgeId> TEdges3;
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
	addEdge(g3,0,1);
	removeEdge(g3,0,4);
	removeEdge(g3,0,2);
	SEQAN_TASSERT(numEdges(g3) == 2)
	removeInEdges(g3,1);
	SEQAN_TASSERT(numEdges(g3) == 0)
	

//____________________________________________________________________________
// Undirected graph iterators
	typedef Graph<Undirected<> > TGraphIter;
	typedef VertexDescriptor<TGraphIter>::Type TVertexDescriptorIter;
	typedef EdgeDescriptor<TGraphIter>::Type TEdgeDescriptorIter;
	
	TGraphIter gIter;
	addVertex(gIter);addVertex(gIter);addVertex(gIter);addVertex(gIter);
	addVertex(gIter);addVertex(gIter);addVertex(gIter);addVertex(gIter);
	removeVertex(gIter,0);
	removeVertex(gIter,5);
	addEdge(gIter,2,7);
	addEdge(gIter,2,3);
	addEdge(gIter,2,4);
	addEdge(gIter,4,3);
	addEdge(gIter,3,6);
	addEdge(gIter,4,6);

	typedef Iterator<TGraphIter, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator itOutEdge(gIter,3);
	// Both ways are fast for undirected graphs
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==6)
	SEQAN_TASSERT(sourceVertex(itOutEdge)==3)
	SEQAN_TASSERT(targetVertex(itOutEdge)==6)
	SEQAN_TASSERT(sourceVertex(gIter, value(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, *itOutEdge)==6)
	SEQAN_TASSERT(atEnd(itOutEdge)==false)
	SEQAN_TASSERT(atBegin(itOutEdge)==true)
	goNext(itOutEdge);
	SEQAN_TASSERT(atEnd(itOutEdge)==false)
	SEQAN_TASSERT(atBegin(itOutEdge)==false)
	SEQAN_TASSERT(sourceVertex(itOutEdge)==3)
	SEQAN_TASSERT(targetVertex(itOutEdge)==4)
	++itOutEdge;
	itOutEdge++;
	SEQAN_TASSERT(atEnd(itOutEdge)==true)
	SEQAN_TASSERT(atBegin(itOutEdge)==false)
	goPrevious(itOutEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==3)
	--itOutEdge;
	itOutEdge--; 
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==6)
	itOutEdge--;
	itOutEdge--;
	SEQAN_TASSERT(atBegin(itOutEdge)==true)
	TOutEdgeIterator itEdge2(itOutEdge);
	TOutEdgeIterator itEdge3;
	itEdge3 = itOutEdge;
	SEQAN_TASSERT(itOutEdge == itEdge2)
	SEQAN_TASSERT(itEdge2 == itEdge3)
	goEnd(itOutEdge);
	SEQAN_TASSERT(itEdge2 != itOutEdge)
	goEnd(itEdge2);
	SEQAN_TASSERT(itEdge2 == itOutEdge)
	goBegin(itEdge2);
	SEQAN_TASSERT(atBegin(itEdge2)==true)
	SEQAN_TASSERT(&gIter == &hostGraph(itOutEdge))

	
	typedef Iterator<TGraphIter, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEdge(gIter);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(atBegin(itEdge)==false)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==7)
	++itEdge;
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	itEdge++;
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	goNext(itEdge);
	SEQAN_TASSERT(atBegin(itEdge)==false)
	SEQAN_TASSERT(atEnd(itEdge)==true)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	--itEdge;	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	itEdge--;	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==7)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==3)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	SEQAN_TASSERT(atEnd(itEdge)==false)
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
	createRoot(g);
	TVertexDescriptor v0 = getRoot(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)

	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v0,'a');
	SEQAN_TASSERT(findEdge(g, 0, 'a') == e1)
	SEQAN_TASSERT(&_getVertexString(g)[0].data_edge[0] == e1)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 1)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(_getId(e1) == 0)
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
	SEQAN_TASSERT(_getId(e2) == 1)
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
	addEdge(g,3,4,'g');
	TEdgeDescriptor my_edge = addEdge(g,3,1,'c');
	SEQAN_TASSERT(_getId(my_edge) == 3)
	addEdge(g,3,0,'t');
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 1)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 3)
	SEQAN_TASSERT(numEdges(g) == 5)
	SEQAN_TASSERT(outDegree(g, v3) == 3)	

	// Output
	// Raw output
	// std::cout << g << ::std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_automaton.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	StandardAutomaton gTmp;
	strm.open(TEST_PATH "my_automaton.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	// Remove edges
	removeEdge(g,3,1,'c');
	removeEdge(g,0,1,'g');
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
	SEQAN_TASSERT(numEdges(g) == 7)
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

	//Transposes the graph in-place
	transpose(g); 
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
	unsigned int len = length(mat,0);
	SEQAN_TASSERT(getValue(mat,1*len+4) == 1)
	SEQAN_TASSERT(getValue(mat,2*len+4) == 1)
	SEQAN_TASSERT(getValue(mat,0*len+2) == 0)

	// Test iterators
	typedef Iterator<StandardAutomaton, VertexIterator>::Type TVertexIterator;
	TVertexIterator itVert(g);
	SEQAN_TASSERT(getValue(itVert) == 1)
	++itVert;
	SEQAN_TASSERT(getValue(itVert) == 2)
	itVert++;
	SEQAN_TASSERT(getValue(itVert) == 4)
	goNext(itVert);
	SEQAN_TASSERT(atEnd(itVert) == true)

	addEdge(g,1,2,'T');
	typedef Iterator<StandardAutomaton, OutEdgeIterator>::Type TOutEdgeIterator;
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

	typedef Iterator<StandardAutomaton, EdgeIterator>::Type TEdgeIterator;
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

	typedef Iterator<StandardAutomaton, AdjacencyIterator>::Type TAdjacencyIterator;
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
// Automaton - Different alphabet
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
	addVertex(g9);
	addVertex(g9);
	Dna aDna('a');
	Dna gDna('g');
	TEdgeDescriptor9 edg1 = addEdge(g9,0,1,aDna,12);
	TEdgeDescriptor9 edg2 = addEdge(g9,1,0,gDna,21);
	TGraph9 g10;
	transpose(g9, g10);
	TEdgeDescriptor9 edg1_10 = findEdge(g10, 1, aDna);
	TEdgeDescriptor9 edg1_11 = findEdge(g10, 0, gDna);
	SEQAN_TASSERT(getCargo(edg1)==12)
	SEQAN_TASSERT(getCargo(edg2)==21)
	SEQAN_TASSERT(getCargo(edg1_10)==12)
	SEQAN_TASSERT(getCargo(edg1_11)==21)
}

//////////////////////////////////////////////////////////////////////////////

void Test_WordGraph() {
//____________________________________________________________________________
// Standard automaton: No edge cargo

	typedef Graph<Automaton<Dna, String<Dna>, WordGraph<> > > TWordGraph;
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
	SEQAN_TASSERT(isRoot(g, 0) == true)
	SEQAN_TASSERT(getRoot(g) == 0)
	assignRoot(g,3);
	SEQAN_TASSERT(getRoot(g) == 3)
	SEQAN_TASSERT(isRoot(g, 0) == false)
	SEQAN_TASSERT(isRoot(g, 3) == true)
	root(g) = 2;
	SEQAN_TASSERT(getRoot(g) == 2)
	SEQAN_TASSERT(isRoot(g, 3) == false)
	SEQAN_TASSERT(isRoot(g, 2) == true)

	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v3,"ag");
	SEQAN_TASSERT(findEdge(g,v0,'a') == e1)
	SEQAN_TASSERT(_getId(e1) == 0)
	// First letter -> edge label, all other letters into the cargo
	SEQAN_TASSERT(getCargo(e1) == "g")
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
	SEQAN_TASSERT(_getId(e2) == 1) 
	SEQAN_TASSERT(v5 == 5)
	SEQAN_TASSERT(numVertices(g) == 6)
	SEQAN_TASSERT(targetVertex(g, e2) == 5)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	removeEdge(g,0,5,String<Dna>("g"));
	SEQAN_TASSERT(numEdges(g) == 1)
	e2 =addEdge(g,0,5,"g");
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 5) == 1)
	SEQAN_TASSERT(degree(g, 0) == 2)
	SEQAN_TASSERT(getSuccessor(g, 0, "g") == 5)
	SEQAN_TASSERT(getSuccessor(g, 0, String<Dna>("ag")) == 3)  // The whole edge label or just the first letter
	SEQAN_TASSERT(getSuccessor(g, 0, "a") == getNil<TVertexDescriptor>())
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
	assignRoot(g,0);

	// Output
	// Raw output
	// std::cout << g << ::std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_wordgraph.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	TWordGraph gTmp;
	strm.open(TEST_PATH "my_wordgraph.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	assignRoot(g,2);
	TWordGraph g_tmp(g);
	SEQAN_TASSERT(numVertices(g_tmp) == 8)
	SEQAN_TASSERT(parseString(g_tmp, 0, "agaggg") == 1)
	SEQAN_TASSERT(inDegree(g_tmp, 5) == 1)
	SEQAN_TASSERT(degree(g_tmp, 0) == 2)
	SEQAN_TASSERT(isRoot(g_tmp, 2) == true)
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
	SEQAN_TASSERT(isRoot(g_tmp, 2) == true)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Tree() {
//____________________________________________________________________________
// Tree without edge cargo

	typedef Graph<Tree<void> > TTree;
	typedef VertexDescriptor<TTree>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TTree>::Type TEdgeDescriptor;
	
	TTree g;
	SEQAN_TASSERT(empty(g) == true)
	createRoot(g);
	TVertexDescriptor rootV = getRoot(g);
	SEQAN_TASSERT(rootV == 0)
	SEQAN_TASSERT(isRoot(g, rootV) == true)
	SEQAN_TASSERT(root(g) == rootV)
	SEQAN_TASSERT(empty(g) == false)
	TVertexDescriptor childC1 = addChild(g,rootV);
	TEdgeDescriptor childC1e = findEdge(g, rootV, childC1);
	SEQAN_TASSERT(_getVertexString(g)[0] == childC1e)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 2)
	SEQAN_TASSERT(targetVertex(g, childC1e) == childC1) // Target in a tree = child
	SEQAN_TASSERT(sourceVertex(g, childC1e) == rootV)  // Source in a tree = parent
	SEQAN_TASSERT(childVertex(g, childC1e) == childC1)  // Shortcuts
	SEQAN_TASSERT(parentVertex(g, childC1e) == rootV)
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(outDegree(g, rootV) == 1)
	TVertexDescriptor childC2 = addChild(g,rootV);
	TVertexDescriptor childC3 = addChild(g,rootV);
	SEQAN_TASSERT(outDegree(g, rootV) == 3)
	SEQAN_TASSERT(childC1 == 1)
	SEQAN_TASSERT(childC2 == 2)
	SEQAN_TASSERT(childC3 == 3)
	TVertexDescriptor childC2C1 = addChild(g,childC2);
	TVertexDescriptor childC2C1C1 = addChild(g,childC2C1);
	TVertexDescriptor childC2C1C1C1 = addChild(g,childC2C1C1);
	TVertexDescriptor childC2C1C1C2 = addChild(g,childC2C1C1);
	TVertexDescriptor childC4 = addChild(g,rootV);
	SEQAN_TASSERT(inDegree(g, childC2C1) == 1) 
	SEQAN_TASSERT(outDegree(g, childC2C1) == 1)
	SEQAN_TASSERT(degree(g, childC2C1) == 2)
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	SEQAN_TASSERT(numTreeEdges(g) == numVertices(g) - 1)
	TEdgeDescriptor childC2C1C1e = findEdge(g, childC2C1C1, childC2C1);
	
	// Raw output
	// std::cout << g << std::endl;
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_tree.dot", ios_base::out | ios_base::trunc);
	write(strm,g,DotDrawing());
	strm.close();
	// File read
	TTree gTmp;
	strm.open(TEST_PATH "my_tree.dot", ios_base::in);
	read(strm,gTmp,DotDrawing());
	strm.close();

	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==2)
	SEQAN_TASSERT(g.data_parent[5]==4)
	SEQAN_TASSERT(g.data_parent[6]==5)
	SEQAN_TASSERT(g.data_parent[7]==5)
	SEQAN_TASSERT(g.data_parent[8]==0)
	_rebuildParentMap(g);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==2)
	SEQAN_TASSERT(g.data_parent[5]==4)
	SEQAN_TASSERT(g.data_parent[6]==5)
	SEQAN_TASSERT(g.data_parent[7]==5)
	SEQAN_TASSERT(g.data_parent[8]==0)
	SEQAN_TASSERT(childVertex(g, childC2C1C1e) == childC2C1C1)  
	SEQAN_TASSERT(parentVertex(g, childC2C1C1e) == childC2C1)
	removeChild(g, rootV, childC2);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[5]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[6]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[7]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[8]==0)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(inDegree(g, rootV) == 0) 
	SEQAN_TASSERT(outDegree(g, rootV) == 3)
	SEQAN_TASSERT(degree(g, rootV) == 3)
	childC2 = addChild(g,rootV);
	childC2C1 = addChild(g,childC2);
	childC2C1C1 = addChild(g,childC2C1);
	childC2C1C1C1 = addChild(g,childC2C1C1);
	childC2C1C1C2 = addChild(g,childC2C1C1);
	removeAllChildren(g, rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(numTreeEdges(g) == 0)
	SEQAN_TASSERT(numVertices(g) == 1) // Just the root
	SEQAN_TASSERT(inDegree(g, rootV) == 0) 
	SEQAN_TASSERT(outDegree(g, rootV) == 0)
	SEQAN_TASSERT(degree(g, rootV) == 0)
	addChild(g,rootV);addChild(g,rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 2)
	clearEdges(g);
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(numVertices(g) == 3)
	SEQAN_TASSERT(empty(g) == false)
	addChild(g,rootV);addChild(g,rootV);
	clearVertices(g);
	SEQAN_TASSERT(empty(g) == true)
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	createRoot(g);
	childC1 = addChild(g,rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 1)
	childC3 = addChild(g,rootV);
	childC2 = addChild(g,rootV);
	childC2C1 = addChild(g,childC2);
	childC2C1C1 = addChild(g,childC2C1);
	childC2C1C1C1 = addChild(g,childC2C1C1);
	childC2C1C1C2 = addChild(g,childC2C1C1);
	childC4 = addChild(g,rootV);
	Matrix<unsigned int> mat; 	// Adjacency matrix
	getAdjacencyMatrix(g, mat);
	unsigned int len = length(mat, 0);
	SEQAN_TASSERT(getValue(mat, 0*len+8) == 1)
	SEQAN_TASSERT(getValue(mat, 8*len+0) == 0)
	SEQAN_TASSERT(getValue(mat, 3*len+0) == 0)
	SEQAN_TASSERT(getValue(mat, 0*len+3) == 1)
	SEQAN_TASSERT(getValue(mat, 0*len+4) == 0)
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	transpose(g); 
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	TTree g_copy(g);
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	clear(g_copy);
	g_copy = g;
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	transpose(g,g_copy);  
	g = g_copy;
	_rebuildParentMap(g);
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==3)
	SEQAN_TASSERT(g.data_parent[5]==4)
	SEQAN_TASSERT(g.data_parent[6]==5)
	SEQAN_TASSERT(g.data_parent[7]==5)
	SEQAN_TASSERT(g.data_parent[8]==0)
	removeOutEdges(g,childC2C1C1);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==3)
	SEQAN_TASSERT(g.data_parent[5]==4)
	SEQAN_TASSERT(g.data_parent[6]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[7]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[8]==0)
	SEQAN_TASSERT(numVertices(g) == 9)
	SEQAN_TASSERT(numEdges(g) == 6)
	removeVertex(g,childC2C1);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[5]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[6]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[7]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[8]==0)
	SEQAN_TASSERT(numEdges(g) == 4)

	SEQAN_TASSERT(numVertices(g) == 8)
	removeInEdges(g,childC2);
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(numVertices(g) == 8)
	removeOutEdges(g,rootV);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[2]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[3]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[5]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[6]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[7]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[8]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(numVertices(g) == 8)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == false) 
	addVertex(g);
	TEdgeDescriptor my_edge = addEdge(g,0,1);
	removeEdge(g,my_edge);

//____________________________________________________________________________
// Tree with cargo

	typedef Pair<char, int> TPair;
	typedef Tree<TPair> TEdges;
	typedef Graph<TEdges> TCargoGraph;
	typedef VertexDescriptor<TCargoGraph>::Type TVertexDescriptor2;
	typedef EdgeDescriptor<TCargoGraph>::Type TEdgeDescriptor2;

	TCargoGraph g2;
	createRoot(g2);
	SEQAN_TASSERT(numVertices(g2) == 1)
	SEQAN_TASSERT(numEdges(g2) == 0)
	TVertexDescriptor2 ver1 = addChild(g2, getRoot(g2), TPair('a',3));
	SEQAN_TASSERT(ver1 == 1)
	SEQAN_TASSERT(numVertices(g2) == 2)
	TVertexDescriptor2 ver2 = addChild(g2, getRoot(g2));
	SEQAN_TASSERT(ver2 == 2)
	SEQAN_TASSERT(numVertices(g2) == 3)
	TEdgeDescriptor2 ed1 =findEdge(g2,getRoot(g2),ver1);
	TEdgeDescriptor2 ed2 =findEdge(g2,getRoot(g2),ver2);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT(targetVertex(g2, ed1) == ver1)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == getRoot(g2))
	SEQAN_TASSERT(numEdges(g2) == 2)
	assignCargo(ed2, TPair('b',4));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT((getCargo(ed2)).i1 == 'b')
	SEQAN_TASSERT((getCargo(ed2)).i2 == 4)
	cargo(ed1) = TPair('c',1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'c')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 1)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Alignment() {
//____________________________________________________________________________
// Alignment without edge weights
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	
	TStringSet str;
	
	TString str0("acaagtaacataaaaaaaaaaaaaaaacccccccccttttttttaaaaa");
	TId id0 = assignValueById(str, str0);

	TString str1("cccaaagggtttttccccccccccccttttttttttaaaaaaagggggggg");
	TId id1 = assignValueById(str, str1);

	TString str2("cacatgtaatcatgggggggggccccccttttaaaaaaaaaaatttt");
	TId id2 = assignValueById(str, str2);


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

	//std::cout << g << std::endl;

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
	Matrix<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = getIdUpperBound(g.data_align.data_id_managerV);
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
	TVD vH = addVertex(gAl, i1, 8, 3);TVD vT = addVertex(gAl, i1, 13, 1);TVD vS = addVertex(gAl, i3, 6, 3);
	TVD vW = addVertex(gAl, i2, 18, 1);TVD vA = addVertex(gAl, i0, 0, 8);TVD vM = addVertex(gAl, i2, 11, 4);
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

	//std::cout << gAl << std::endl;

	// Graph drawing
	// Raw output
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_alignment.dot", ios_base::out | ios_base::trunc);
	write(strm,gAl,DotDrawing());
	strm.close();

}

}

#endif

