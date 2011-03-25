 /*==========================================================================
                  HSA - Hierarchical Segment-based Alignment

 ============================================================================
  Copyright (C) 2011 by Birte Kehr

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#ifndef SEQAN_HEADER_HSA_TYPES_H
#define SEQAN_HEADER_HSA_TYPES_H

#include <seqan/sequence.h>

using namespace seqan;

template<typename TSequence>
struct MyOptions {
	typedef Graph<Tree<double> > TTree;
	typedef std::map<typename Id<Graph<Alignment<TSequence> > >::Type, void const *> TIdMap;

	CharString sequenceFiles;
	String<CharString> fileNames;
	CharString outputFile;
	bool phylogeny;
	bool verbose;

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
	TTree guideTree; // guide tree from level 0  (on full sequences)
	TIdMap idMap; // map from string ids to stringset ids of level 0 guide tree
	
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
		phylogeny = false;
		verbose = false;

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
