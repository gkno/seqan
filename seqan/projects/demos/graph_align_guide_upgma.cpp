/// This code example illustrates a guide tree creation using upgma
#include <seqan/graph.h>
#include <iostream>

using namespace seqan;

int main() {
/// Create an artificial distance matrix
	String<double> mat;
	fill(mat, 8*8, 0);
	assignValue(mat, 0*8+0, 0);assignValue(mat, 0*8+1, 32);assignValue(mat, 0*8+2, 48);assignValue(mat, 0*8+3, 51);
	assignValue(mat, 0*8+4, 50);assignValue(mat, 0*8+5, 48);assignValue(mat, 0*8+6, 98);assignValue(mat, 0*8+7, 148);
	assignValue(mat, 1*8+0, 32);assignValue(mat, 1*8+1, 0);assignValue(mat, 1*8+2, 26);assignValue(mat, 1*8+3, 34);
	assignValue(mat, 1*8+4, 29);assignValue(mat, 1*8+5, 33);assignValue(mat, 1*8+6, 84);assignValue(mat, 1*8+7, 136);
	assignValue(mat, 2*8+0, 48);assignValue(mat, 2*8+1, 26);assignValue(mat, 2*8+2, 0);assignValue(mat, 2*8+3, 42);
	assignValue(mat, 2*8+4, 44);assignValue(mat, 2*8+5, 44);assignValue(mat, 2*8+6, 92);assignValue(mat, 2*8+7, 152);
	assignValue(mat, 3*8+0, 51);assignValue(mat, 3*8+1, 34);assignValue(mat, 3*8+2, 42);assignValue(mat, 3*8+3, 0);
	assignValue(mat, 3*8+4, 44);assignValue(mat, 3*8+5, 38);assignValue(mat, 3*8+6, 86);assignValue(mat, 3*8+7, 142);
	assignValue(mat, 4*8+0, 50);assignValue(mat, 4*8+1, 29);assignValue(mat, 4*8+2, 44);assignValue(mat, 4*8+3, 44);
	assignValue(mat, 4*8+4, 0);assignValue(mat, 4*8+5, 24);assignValue(mat, 4*8+6, 89);assignValue(mat, 4*8+7, 142);
	assignValue(mat, 5*8+0, 48);assignValue(mat, 5*8+1, 33);assignValue(mat, 5*8+2, 44);assignValue(mat, 5*8+3, 38);
	assignValue(mat, 5*8+4, 24);assignValue(mat, 5*8+5, 0);assignValue(mat, 5*8+6, 90);assignValue(mat, 5*8+7, 142);
	assignValue(mat, 6*8+0, 98);assignValue(mat, 6*8+1, 84);assignValue(mat, 6*8+2, 92);assignValue(mat, 6*8+3, 86);
	assignValue(mat, 6*8+4, 89);assignValue(mat, 6*8+5, 90);assignValue(mat, 6*8+6, 0);assignValue(mat, 6*8+7, 148);
	assignValue(mat, 7*8+0, 148);assignValue(mat, 7*8+1, 136);assignValue(mat, 7*8+2, 152);assignValue(mat, 7*8+3, 142);
	assignValue(mat, 7*8+4, 142);assignValue(mat, 7*8+5, 142);assignValue(mat, 7*8+6, 148);assignValue(mat, 7*8+7, 0);
/// Out-parameter: A guide tree based upon a distance matrix
	typedef Graph<Tree<double> > TGraph;
	TGraph upgmaTreeOut;
/// UPGMA
	upgmaTree(mat, upgmaTreeOut);
/// Console output
	std::cout << upgmaTreeOut << std::endl;
	// Dfs Iterator
	typedef Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
	TDfsPreorder dfsIt(upgmaTreeOut,getRoot(upgmaTreeOut));
	std::cout << "Dfs: ";
	for(;!atEnd(dfsIt);goNext(dfsIt)) {
		std::cout << *dfsIt << ", ";
	}
	std::cout << std::endl;
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIter;
	TEdgeIter edIt(upgmaTreeOut);
	std::cout << "Edges: " << std::endl;
	for(;!atEnd(edIt);goNext(edIt)) {
		std::cout << sourceVertex(edIt) << " -- " << targetVertex(edIt) << ": Weight = " << getCargo(*edIt) << std::endl;
	}
	return 0;
}
