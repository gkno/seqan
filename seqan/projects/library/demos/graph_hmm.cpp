#include <iostream>
#include <fstream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main() {

	typedef double TProbability;
	typedef Dna TAlphabet;
	typedef Size<TAlphabet>::Type TSize;
	typedef Graph<Hmm<TAlphabet, TProbability> > THmm;
	typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
	typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;
	
	Dna dnaA = Dna('A');
	Dna dnaC = Dna('C');
	Dna dnaG = Dna('G');
	Dna dnaT = Dna('T');

	// Create an empty HMM
	THmm hmm;

	// Add begin state
	TVertexDescriptor begState = addVertex(hmm);
	assignBeginState(hmm, begState);

	// Add exon state
	TVertexDescriptor exonState = addVertex(hmm);
	emissionProbability(hmm, exonState, dnaA) = 0.25;
	emissionProbability(hmm, exonState, dnaC) = 0.25;
	emissionProbability(hmm, exonState, dnaG) = 0.25;
	emissionProbability(hmm, exonState, dnaT) = 0.25;

	// Add 5' splice site
	TVertexDescriptor spliceState = addVertex(hmm);
	emissionProbability(hmm, spliceState, dnaA) = 0.05;
	emissionProbability(hmm, spliceState, dnaC) = 0.0;
	emissionProbability(hmm, spliceState, dnaG) = 0.95;
	emissionProbability(hmm, spliceState, dnaT) = 0;

	// Add intron state
	TVertexDescriptor intronState = addVertex(hmm);
	emissionProbability(hmm, intronState, dnaA) = 0.4;
	emissionProbability(hmm, intronState, dnaC) = 0.1;
	emissionProbability(hmm, intronState, dnaG) = 0.1;
	emissionProbability(hmm, intronState, dnaT) = 0.4;

	// Add end state
	TVertexDescriptor eState = addVertex(hmm);
	assignEndState(hmm, eState);

	// Add edges (transitions)
	addEdge(hmm, exonState, exonState, 0.9);
	addEdge(hmm, exonState, spliceState, 0.1);
	addEdge(hmm, spliceState, intronState, 1.0);
	addEdge(hmm, begState, exonState, 1.0);
	addEdge(hmm, intronState, intronState, 0.9);
	addEdge(hmm, intronState, eState, 0.1);

	// Print the whole model
	std::cout << hmm << std::endl;


	// Viterbi algorithm
	String<Dna> sequence = "CTTCATGTGAAAGCAGACGTAAGTCA";
	String<TVertexDescriptor> path;
	TProbability p = viterbiAlgorithm(hmm, sequence, path);
	std::cout << "Viterbi algorithm" << std::endl;
	std::cout << "Log probability of best path: " << p << std::endl;
	std::cout << "Best path: " << std::endl;
	for(TSize i = 0; i<length(sequence); ++i) std::cout << sequence[i] << ',';
	std::cout << std::endl;
	std::cout << "State path: " << std::endl;
	for(TSize i = 0; i<length(path); ++i) {
		std::cout << path[i];
		if (isSilent(hmm, path[i])) std::cout << " (Silent)";
		if (i < length(path) - 1) std::cout << ',';
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





	//// HMM with silent states
	//// Silent states must be numbered in increasing order from left to right!

	//typedef double TProbability;
	//typedef Dna TAlphabet;
	//typedef Size<TAlphabet>::Type TSize;
	//typedef Graph<Hmm<TAlphabet, TProbability> > THmm;
	//typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
	//typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;
	//
	//Dna dnaA = Dna('A');
	//Dna dnaC = Dna('C');
	//Dna dnaG = Dna('G');
	//Dna dnaT = Dna('T');

	//// Create an empty HMM
	//THmm hmm;

	//// Add begin state
	//TVertexDescriptor begState = addVertex(hmm);
	//assignBeginState(hmm, begState);

	//// Add emission state
	//TVertexDescriptor emitState1 = addVertex(hmm);
	//emissionProbability(hmm, emitState1, dnaA) = 0.0;
	//emissionProbability(hmm, emitState1, dnaC) = 0.8;
	//emissionProbability(hmm, emitState1, dnaG) = 0.1;
	//emissionProbability(hmm, emitState1, dnaT) = 0.1;

	//// Add emission state
	//TVertexDescriptor emitState2 = addVertex(hmm);
	//emissionProbability(hmm, emitState2, dnaA) = 0.0;
	//emissionProbability(hmm, emitState2, dnaC) = 0.2;
	//emissionProbability(hmm, emitState2, dnaG) = 0.2;
	//emissionProbability(hmm, emitState2, dnaT) = 0.6;

	//// Add emission state
	//TVertexDescriptor emitState3 = addVertex(hmm);
	//emissionProbability(hmm, emitState3, dnaA) = 1.0;
	//emissionProbability(hmm, emitState3, dnaC) = 0.0;
	//emissionProbability(hmm, emitState3, dnaG) = 0.0;
	//emissionProbability(hmm, emitState3, dnaT) = 0.0;

	//// Add delete states
	//TVertexDescriptor delState1 = addVertex(hmm, true);
	//TVertexDescriptor delState2 = addVertex(hmm, true);
	//TVertexDescriptor delState3 = addVertex(hmm, true);

	//// Add end state
	//TVertexDescriptor eState = addVertex(hmm);
	//assignEndState(hmm, eState);


	//// Add edges (transitions)
	//addEdge(hmm, begState, emitState1, 0.5);
	//addEdge(hmm, begState, delState1, 0.5);
	//addEdge(hmm, emitState1, emitState2, 0.5);
	//addEdge(hmm, emitState1, delState2, 0.5);
	//addEdge(hmm, delState1, emitState2, 0.5);
	//addEdge(hmm, delState1, delState2, 0.5);
	//addEdge(hmm, emitState2, emitState3, 0.5);
	//addEdge(hmm, emitState2, delState3, 0.5);
	//addEdge(hmm, delState2, emitState3, 0.5);
	//addEdge(hmm, delState2, delState3, 0.5);
	//addEdge(hmm, emitState3, eState, 1.0);
	//addEdge(hmm, delState3, eState, 1.0);

	//// Print the whole model
	//std::cout << hmm << std::endl;


	//// Viterbi algorithm
	//String<Dna> sequence = "CA";
	//String<TVertexDescriptor> path;
	//TProbability p = viterbiAlgorithm(hmm, sequence, path);
	//std::cout << "Viterbi algorithm" << std::endl;
	//std::cout << "Log probability of best path: " << p << std::endl;
	//std::cout << "Best path: " << std::endl;
	//for(TSize i = 0; i<length(sequence); ++i) std::cout << sequence[i] << ',';
	//std::cout << std::endl;
	//std::cout << "State path: " << std::endl;
	//for(TSize i = 0; i<length(path); ++i) {
	//	std::cout << path[i];
	//	if (isSilent(hmm, path[i])) std::cout << " (Silent)";
	//	if (i < length(path) - 1) std::cout << ',';
	//}
	//std::cout << std::endl;

	//// Forward algorithm
	//std::cout << "Forward algorithm" << std::endl;
	//p = forwardAlgorithm(hmm, sequence);
	//std::cout << "Probability that the HMM generated the sequence: " << p << std::endl;

	//// Backward algorithm
	//std::cout << "Backward algorithm" << std::endl;
	//p = backwardAlgorithm(hmm, sequence);
	//std::cout << "Probability that the HMM generated the sequence: " << p << std::endl;


	return 0;
}
