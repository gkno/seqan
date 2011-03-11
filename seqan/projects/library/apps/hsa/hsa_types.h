
#ifndef SEQAN_HEADER_HSA_TYPES_H
#define SEQAN_HEADER_HSA_TYPES_H

#include <seqan/sequence.h>

using namespace seqan;

struct MyOptions {
	typedef Graph<Tree<double> > TTree;

	CharString sequenceFiles;
	String<CharString> fileNames;
	CharString outputFile;

	unsigned recursions;
	unsigned initialMinLength;
	double initialEpsilon;
	unsigned deltaMinLength;
	double deltaEpsilon;

	// Should be set to true if the same guide tree computed from
	//   the full sequences shall be used in all recursion steps.
	// If set to false, a separate guide tree will be computed in
	//   each recursion step.
	bool globalGuideTree;
	TTree guideTree;
	
	// Matches from higher recursion levels will be carried to lower
	//   levels. Should be set to true if these matches shall obtain
	//   an infinite weight to guarantee that they appear in the final
	//   alignment and cannot be outweighted by lower level alignments.
	bool fixedHigherLevelMatches;
	
	// Should be set to true if not all sequence between cuts through
	//   the multiple alignment shall be pairwise compared during
	//   recursion but only those that can match according to anchors
	//   from higher recursions levels. This implies that the matches
	//   from higher recursion levels are fixed.
	bool anchoredPairwiseComparison;

	MyOptions() {
		outputFile = "ReSeAl.dot";

		recursions = 3;
		initialMinLength = 100;
		initialEpsilon = 0.1;

		deltaMinLength = 30;
		deltaEpsilon = 0.0;

		globalGuideTree = false;
		anchoredPairwiseComparison = false;
		fixedHigherLevelMatches = false;
	}

	MyOptions(bool tree, bool fixed, bool anchored) {
		MyOptions();
		SEQAN_ASSERT_EQ(anchored && fixed, fixed);
		globalGuideTree = tree;
		fixedHigherLevelMatches = fixed;
		anchoredPairwiseComparison = anchored;
	}
};

struct StellarParams {
	double epsilon;
	int minLength;
	unsigned xdrop;
	unsigned qgram;

	StellarParams(double eps, unsigned n0) {
		epsilon = eps;
		minLength = n0;
		xdrop = 5;

		unsigned n1 = ceil((floor(eps*n0) + 1) / eps);
		qgram = _min(ceil((n0 - floor(eps*n0)) / (floor(eps*n0) + 1)),
			         ceil((n1 - floor(eps*n1)) / (floor(eps*n1) + 1)));
	}

	StellarParams(double eps, unsigned n0, unsigned x, unsigned q) {
		epsilon = eps;
		minLength = n0;
		xdrop = x;
		qgram = q;
	}
};

#endif
