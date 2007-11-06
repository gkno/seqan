/// This code example illustrates a guide tree creation using neighbor-joing
#include <seqan/graph.h>
#include <iostream>

using namespace seqan;

int main() {
/// Create an artificial distance matrix
	String<double> mat;
	fill(mat, 8*8, 0);
	assignValue(mat, 0*8+0, 0);assignValue(mat, 0*8+1, 7);assignValue(mat, 0*8+2, 8);assignValue(mat, 0*8+3, 11);
	assignValue(mat, 0*8+4, 13);assignValue(mat, 0*8+5, 16);assignValue(mat, 0*8+6, 13);assignValue(mat, 0*8+7, 17);
	assignValue(mat, 1*8+0, 7);assignValue(mat, 1*8+1, 0);assignValue(mat, 1*8+2, 5);assignValue(mat, 1*8+3, 8);
	assignValue(mat, 1*8+4, 10);assignValue(mat, 1*8+5, 13);assignValue(mat, 1*8+6, 10);assignValue(mat, 1*8+7, 14);
	assignValue(mat, 2*8+0, 8);assignValue(mat, 2*8+1, 5);assignValue(mat, 2*8+2, 0);assignValue(mat, 2*8+3, 5);
	assignValue(mat, 2*8+4, 7);assignValue(mat, 2*8+5, 10);assignValue(mat, 2*8+6, 7);assignValue(mat, 2*8+7, 11);
	assignValue(mat, 3*8+0, 11);assignValue(mat, 3*8+1, 8);assignValue(mat, 3*8+2, 5);assignValue(mat, 3*8+3, 0);
	assignValue(mat, 3*8+4, 8);assignValue(mat, 3*8+5, 11);assignValue(mat, 3*8+6, 8);assignValue(mat, 3*8+7, 12);
	assignValue(mat, 4*8+0, 13);assignValue(mat, 4*8+1, 10);assignValue(mat, 4*8+2, 7);assignValue(mat, 4*8+3, 8);
	assignValue(mat, 4*8+4, 0);assignValue(mat, 4*8+5, 5);assignValue(mat, 4*8+6, 6);assignValue(mat, 4*8+7, 10);
	assignValue(mat, 5*8+0, 16);assignValue(mat, 5*8+1, 13);assignValue(mat, 5*8+2, 10);assignValue(mat, 5*8+3, 11);
	assignValue(mat, 5*8+4, 5);assignValue(mat, 5*8+5, 0);assignValue(mat, 5*8+6, 9);assignValue(mat, 5*8+7, 13);
	assignValue(mat, 6*8+0, 13);assignValue(mat, 6*8+1, 10);assignValue(mat, 6*8+2, 7);assignValue(mat, 6*8+3, 8);
	assignValue(mat, 6*8+4, 6);assignValue(mat, 6*8+5, 9);assignValue(mat, 6*8+6, 0);assignValue(mat, 6*8+7, 8);
	assignValue(mat, 7*8+0, 17);assignValue(mat, 7*8+1, 14);assignValue(mat, 7*8+2, 11);assignValue(mat, 7*8+3, 12);
	assignValue(mat, 7*8+4, 10);assignValue(mat, 7*8+5, 13);assignValue(mat, 7*8+6, 8);assignValue(mat, 7*8+7, 0);
/// Out-parameter: A guide tree based upon a distance matrix
	typedef Graph<Tree<double> > TGraph;
	TGraph njTreeOut;
/// Neighbor-joining
	slowNjTree(mat, njTreeOut);
/// Console output
	std::cout << njTreeOut << std::endl;
	typedef Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder dfsIt(njTreeOut,getRoot(njTreeOut));
	std::cout << "Dfs: ";
	for(;!atEnd(dfsIt);goNext(dfsIt)) {
		std::cout << *dfsIt << ", ";
	}
	std::cout << std::endl;
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIter;
	TEdgeIter edIt(njTreeOut);
	std::cout << "Edges: " << std::endl;
	for(;!atEnd(edIt);goNext(edIt)) {
		std::cout << sourceVertex(edIt) << " -- " << targetVertex(edIt) << ": Weight = " << getCargo(*edIt) << std::endl;
	}
	return 0;
}
