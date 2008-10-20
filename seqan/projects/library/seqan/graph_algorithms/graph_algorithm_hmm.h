 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_algorithm_hmm.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/


#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_HMM_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_HMM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Basic HMM algorithms
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.viterbiAlgorithm:
..cat:Graph.Hmm
..summary:Implements the viterbi algorithm.
..signature:viterbiAlgorithm(hmm, seq, path)
..param.hmm:In-parameter:Input HMM.
...type:Spec.Hmm
..param.seq:In-parameter:Input sequence.
..param.path:Out-parameter:State path.
..returns:TCargo
...remarks:Probability of the path.
..see:Function.forwardAlgorithm
..see:Function.backwardAlgorithm
*/
template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequence, typename TPath>
inline TProbability
viterbiAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				 TSequence const& seq,
				 TPath& path)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	String<TProbability> vMat;
	String<TSize> traceback;
	TSize numCols = length(seq) + 2;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	fill(vMat, numCols * numRows, 0.0);
	resize(traceback, numCols * numRows);
	value(vMat, getBeginState(hmm)) = 1.0;
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Initialization for silent states connected to the begin state
	TVertexIterator itSilent(hmm);
	for(;!atEnd(itSilent);++itSilent) {
		if (!isSilent(hmm, value(itSilent))) continue;
		if ((value(itSilent) == bState) || (value(itSilent) == eState)) continue;
		TProbability maxValue = 0.0;
		TVertexDescriptor maxVertex = nilVertex;

		// Find maximum
		// A vertex iterator guarantees that the vertices are processend in increasing order
		// That is, the smallest silent states come first!!!
		TVertexIterator itMax(hmm);
		for(;!atEnd(itMax);++itMax) {
			if ((!isSilent(hmm, value(itMax))) ||
				((isSilent(hmm, value(itMax))) && (value(itMax) < value(itSilent)))) {
					TProbability local = value(vMat, value(itMax)) * getTransitionProbability(hmm, value(itMax), value(itSilent));
					if (local > maxValue) {
						maxValue = local;
						maxVertex = *itMax;
					}
			}
		}

		// Set traceback vertex
		if (maxVertex != nilVertex) {
			value(vMat, value(itSilent)) = maxValue;		
			value(traceback, value(itSilent)) = maxVertex;
		}
	}

	// Recurrence
	TSize len = length(seq);
	for(TSize i=1; i<=len; ++i) {

		// Iterate over real states
		TVertexIterator itV(hmm);
		for(;!atEnd(itV);++itV) {
			if ((value(itV) == bState) || (value(itV) == eState)) continue;
			if (isSilent(hmm, value(itV))) continue;
			TProbability maxValue = 0.0;
			TVertexDescriptor maxVertex = nilVertex;

			// Find maximum
			TVertexIterator itMax(hmm);
			for(;!atEnd(itMax);++itMax) {
				TProbability local = value(vMat, (i-1) * numRows + value(itMax)) * getTransitionProbability(hmm, value(itMax), value(itV));
				if (local > maxValue) {
					maxValue = local;
					maxVertex = *itMax;
				}
			}

			// Set traceback vertex
			if (maxVertex != nilVertex) {
				value(vMat, i * numRows + *itV) = maxValue * getEmissionProbability(hmm, value(itV), value(seq, i-1));;
				value(traceback, i * numRows + value(itV)) = maxVertex;
			}
		}

		// Iterate over silent states
		goBegin(itV);
		for(;!atEnd(itV);++itV) {
			if ((value(itV) == bState) || (value(itV) == eState)) continue;
			if (!isSilent(hmm, value(itV))) continue;
			TProbability maxValue = 0.0;
			TVertexDescriptor maxVertex = nilVertex;

			// Find maximum
			// A vertex iterator guarantees that the vertices are processend in increasing order
			// That is, the smallest silent states come first!!!
			TVertexIterator itMax(hmm);
			for(;!atEnd(itMax);++itMax) {
				if ((!isSilent(hmm, value(itMax))) ||
					((isSilent(hmm, value(itMax))) && (value(itMax) < value(itV)))) {
						TProbability local = value(vMat, i * numRows + value(itMax)) * getTransitionProbability(hmm, value(itMax), value(itV));
						if (local > maxValue) {
							maxValue = local;
							maxVertex = *itMax;
						}
				}
			}
			// Set traceback vertex
			if (maxVertex != nilVertex) {
				value(traceback, i * numRows + value(itV)) = maxVertex;
				value(vMat, i * numRows + *itV) = maxValue;		
			}
		}
	}

	// Termination
	TProbability maxValue = 0.0;
	TVertexDescriptor maxVertex = 0;
	TVertexIterator itMax(hmm);
	for(;!atEnd(itMax);++itMax) {
		TProbability local = value(vMat, len * numRows + *itMax) * getTransitionProbability(hmm, value(itMax), eState);
		if (local > maxValue) {
			maxValue = local;
			maxVertex = value(itMax);
		}
	}
	value(traceback, (len + 1) * numRows + eState) = maxVertex;
	if (maxVertex != nilVertex) value(vMat, (len+1) * numRows + eState) = maxValue;

	// Traceback
	clear(path);
	TVertexDescriptor oldState = eState;
	appendValue(path, oldState);
	for(TSize i = len + 1; i>=1; --i) {
		do {
			if ((!isSilent(hmm, oldState)) || (oldState == eState)) oldState = value(traceback, i * numRows + oldState);
			else oldState = value(traceback, (i - 1) * numRows + oldState);
			appendValue(path, oldState);
		} while ((isSilent(hmm, oldState)) && (oldState != bState));
	}
	std::reverse(begin(path), end(path));
	
	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(vMat, j*numRows + i) << ',';
	//		//std::cout << value(traceback, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//for(TSize i = 0; i<length(path); ++i) {
	//	std::cout << path[i] << ',';
	//}
	//std::cout << std::endl;

	return value(vMat, (len+1) * numRows + eState);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.forwardAlgorithm:
..cat:Graph.Hmm
..summary:Implements the forward algorithm.
..signature:forwardAlgorithm(hmm, seq)
..param.hmm:In-parameter:Input HMM.
...type:Spec.Hmm
..param.seq:In-parameter:Input sequence.
..returns:TProbability
...remarks:Probability of the sequence.
..see:Function.viterbiAlgorithm
..see:Function.backwardAlgorithm
*/
template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequence>
inline TProbability
forwardAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				 TSequence const& seq)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	String<TProbability> fMat;
	TSize numCols = length(seq) + 2;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	fill(fMat, numCols * numRows, 0.0);
	value(fMat, getBeginState(hmm)) = 1.0;
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);

	// Initialization for silent states connected to the begin state
	TVertexIterator itSilent(hmm);
	for(;!atEnd(itSilent);++itSilent) {
		if (!isSilent(hmm, value(itSilent))) continue;
		if ((value(itSilent) == bState) || (value(itSilent) == eState)) continue;
		TProbability sumValue = 0.0;
		
		// Add up all real states and silent states with lower number in increasing order		
		TVertexIterator itMax(hmm);
		for(;!atEnd(itMax);++itMax) {
			if ((!isSilent(hmm, value(itMax))) ||
				((isSilent(hmm, value(itMax))) && (value(itMax) < value(itSilent)))) {
					sumValue += value(fMat, value(itMax)) * getTransitionProbability(hmm, value(itMax), value(itSilent));
			}
		}
		value(fMat, value(itSilent)) = sumValue;		
	}

	// Recurrence
	TSize len = length(seq);
	for(TSize i=1; i<=len; ++i) {
		// Iterate over real states
		TVertexIterator itV(hmm);
		for(;!atEnd(itV);++itV) {
			if ((value(itV) == bState) || (value(itV) == eState)) continue;
			if (isSilent(hmm, value(itV))) continue;
			TProbability sum = 0.0;
			TVertexIterator itAll(hmm);
			for(;!atEnd(itAll);++itAll) sum += value(fMat, (i-1) * numRows + value(itAll)) * getTransitionProbability(hmm, value(itAll), value(itV));
			value(fMat, i * numRows + value(itV)) = getEmissionProbability(hmm, value(itV), seq[i-1]) * sum;
		}

		// Iterate over silent states
		goBegin(itV);
		for(;!atEnd(itV);++itV) {
			if ((value(itV) == bState) || (value(itV) == eState)) continue;
			if (!isSilent(hmm, value(itV))) continue;
			TProbability sumValue = 0.0;

			// Add up all real states and silent states with lower number in increasing order		
			TVertexIterator itMax(hmm);
			for(;!atEnd(itMax);++itMax) {
				if ((!isSilent(hmm, value(itMax))) ||
					((isSilent(hmm, value(itMax))) && (value(itMax) < value(itV)))) {
						sumValue += value(fMat, i * numRows + value(itMax)) * getTransitionProbability(hmm, value(itMax), value(itV));
				}
			}
			value(fMat, i * numRows + value(itV)) = sumValue;
		}
	}

	// Termination
	TProbability sum = 0.0;
	TVertexIterator itAll(hmm);
	for(;!atEnd(itAll);++itAll) {
		sum += value(fMat, len * numRows + value(itAll)) * getTransitionProbability(hmm, value(itAll), eState);
	}
	value(fMat, (len+1) * numRows + eState) = sum;

	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(fMat, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return value(fMat, (len+1) * numRows + eState);
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.backwardAlgorithm:
..cat:Graph.Hmm
..summary:Implements the backward algorithm.
..signature:backwardAlgorithm(hmm, seq)
..param.hmm:In-parameter:Input HMM.
...type:Spec.Hmm
..param.seq:In-parameter:Input sequence.
..returns:TProbability
...remarks:Probability of the sequence.
..see:Function.viterbiAlgorithm
..see:Function.forwardAlgorithm
*/
template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequence>
inline TProbability
backwardAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				  TSequence const& seq)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	String<TProbability> bMat;
	TSize numCols = length(seq) + 1;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	fill(bMat, numCols * numRows, 0.0);
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);
	TSize len = length(seq);
	TVertexIterator itAll(hmm);
	for(;!atEnd(itAll);++itAll) value(bMat, len * numRows + value(itAll)) = getTransitionProbability(hmm, value(itAll), eState);
	
	// Recurrence
	for(TSize i=len - 1; i>0; --i) {
		TVertexIterator itV(hmm);
		for(;!atEnd(itV);++itV) {
			TProbability sum = 0.0;
			TVertexIterator itAll(hmm);
			for(;!atEnd(itAll);++itAll) sum += value(bMat, (i+1) * numRows + value(itAll)) * getTransitionProbability(hmm, value(itV), value(itAll)) * getEmissionProbability(hmm, value(itAll), value(seq, i));
			value(bMat, i * numRows + value(itV)) =  sum;
		}
	}

	// Termination
	TProbability sum = 0.0;
	goBegin(itAll);
	for(;!atEnd(itAll);++itAll) {
		sum += value(bMat, 1 * numRows + value(itAll)) * getTransitionProbability(hmm, bState, value(itAll)) * getEmissionProbability(hmm, value(itAll), value(seq, 0));
	}
	value(bMat, bState) = sum;

	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(bMat, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return value(bMat, bState);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
