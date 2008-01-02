#include <iostream>
#include <fstream>
#include <seqan/graph.h>

using namespace seqan;

int main() {

	typedef double TProbability;
	typedef Dna TAlphabet;
	typedef Size<TAlphabet>::Type TSize;
	typedef Graph<Hmm<TAlphabet, TProbability> > THmm;
	typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
	typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;
	TSize alph_size = ValueSize<TAlphabet>::VALUE;
	
	Dna dnaA = Dna('A');
	Dna dnaC = Dna('C');
	Dna dnaG = Dna('G');
	Dna dnaT = Dna('T');

	// Create an empty HMM
	THmm hmm;

	// Add state1
	TVertexDescriptor state1 = addVertex(hmm);
	emissionProbability(hmm, state1, dnaA) = 0.2;
	emissionProbability(hmm, state1, dnaC) = 0.2;
	emissionProbability(hmm, state1, dnaG) = 0.3;
	emissionProbability(hmm, state1, dnaT) = 0.3;

	// Add state2
	String<TProbability> emis;
	resize(emis, alph_size);
	value(emis, (Byte) dnaA) = 0.5;
	value(emis, (Byte) dnaC) = 0.5;
	value(emis, (Byte) dnaG) = 0.0;
	value(emis, (Byte) dnaT) = 0.0;
	TVertexDescriptor state2 = addVertex(hmm, emis);

	// Add state3
	TVertexDescriptor state3 = addVertex(hmm, emis);
	assignEmissionProbability(hmm, state3, dnaA, 0.3);
	assignEmissionProbability(hmm, state3, dnaC, 0.3);
	assignEmissionProbability(hmm, state3, dnaG, 0.2);
	assignEmissionProbability(hmm, state3, dnaT, 0.2);

	// Add edges (transitions)
	addEdge(hmm, state1, state1, 0.95);
	TEdgeDescriptor e = addEdge(hmm, state1, state3);
	assignTransitionProbability(hmm, e, 0.05);
	e = addEdge(hmm, state3, state3);
	transitionProbability(hmm, e) = 0.4;
	e = addEdge(hmm, state3, state1);
	transitionProbability(hmm, state3, state1) = 0.1;
	e = addEdge(hmm, state2, state2);
	assignTransitionProbability(hmm, state2, state2, 1.0);
	
	// Add begin and end state
	TVertexDescriptor begState = addVertex(hmm);
	TVertexDescriptor eState = addVertex(hmm);
	addEdge(hmm, begState, state1, 1.0);
	addEdge(hmm, state3, eState, 0.5);
	addEdge(hmm, eState, eState, 1.0);
	beginState(hmm) = state3;
	assignBeginState(hmm, begState);
	endState(hmm) = state3;
	assignEndState(hmm, eState);

	// Print the whole model
	std::cout << hmm << std::endl;
	std::cout << "----" << std::endl;

	// Print parts of the model
	std::cout << getTransitionProbability(hmm, state1, state3) << std::endl;
	std::cout << getTransitionProbability(hmm, e) << std::endl;
	std::cout << getEmissionProbability(hmm, state1, dnaA) << std::endl;
	std::cout << "----" << std::endl;

	// Change model
	removeVertex(hmm, state2);
	std::cout << hmm << std::endl;
	std::cout << "----" << std::endl;

	// Viterbi algorithm
	String<Dna> sequence;
	appendValue(sequence, dnaA);
	appendValue(sequence, dnaC);
	String<TVertexDescriptor> path;
	TProbability p = viterbiAlgorithm(hmm, sequence, path);
	std::cout << "Viterbi algorithm" << std::endl;
	std::cout << "Probability of best path: " << p << std::endl;
	std::cout << "Best path: ";
	for(TSize i = 0; i<length(path); ++i) {
		std::cout << path[i] << ',';
	}
	std::cout << std::endl;

	// Forward algorithm
	std::cout << "Forward algorithm" << std::endl;
	p = forwardAlgorithm(hmm, sequence);
	std::cout << "Probability that the HMM generated the sequence: " << p << std::endl;

	// Backward algorithm
	std::cout << "Backward algorithm" << std::endl;
	p = backwardAlgorithm(hmm, sequence);
	std::cout << "Probability that the HMM generated the sequence: " << p << std::endl;

	return 0;
}
