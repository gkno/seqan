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
template<typename TAlphabet, typename TCargo, typename TSpec, typename TSequence, typename TPath>
inline TCargo
viterbiAlgorithm(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& hmm,
				 TSequence const& seq,
				 TPath& path)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	String<TCargo> vMat;
	String<TSize> traceback;
	TSize numCols = length(seq) + 2;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	TCargo infVal = infimumValue<TCargo>();
	fill(vMat, numCols * numRows, infVal);
	reserve(traceback, numCols * numRows);
	value(vMat, getBeginState(hmm)) = (TCargo) 0;
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Initialization for silent states connected to the begin state
	TVertexIterator itSilent(hmm);
	for(;!atEnd(itSilent);++itSilent) {
		if (!isSilent(hmm, value(itSilent))) continue;
		if ((value(itSilent) == bState) || (value(itSilent) == eState)) continue;
		TCargo maxValue = infVal;
		TVertexDescriptor maxVertex = nilVertex;

		// Find maximum
		// A vertex iterator guarantees that the vertices are processend in increasing order
		// That is, the smallest silent states come first!!!
		TVertexIterator itMax(hmm);
		for(;!atEnd(itMax);++itMax) {
			if ((!isSilent(hmm, value(itMax))) ||
				((isSilent(hmm, value(itMax))) && (value(itMax) < value(itSilent)))) {
					TCargo local = value(vMat, value(itMax)) + std::log( (double) getTransitionProbability(hmm, *itMax, *itSilent));
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
			TCargo maxValue = infVal;
			TVertexDescriptor maxVertex = nilVertex;

			// Find maximum
			TVertexIterator itMax(hmm);
			for(;!atEnd(itMax);++itMax) {
				TCargo local = value(vMat, (i-1) * numRows + value(itMax)) + std::log( (double) getTransitionProbability(hmm, *itMax, *itV));
				if (local > maxValue) {
					maxValue = local;
					maxVertex = *itMax;
				}
			}

			// Set traceback vertex
			TCargo emis = std::log( (double) getEmissionProbability(hmm, *itV, value(seq, i-1)));
			if ((maxVertex != nilVertex) && (emis > infVal)) {
				value(vMat, i * numRows + *itV) = emis + maxValue;
				value(traceback, i * numRows + value(itV)) = maxVertex;
			}
		}

		// Iterate over silent states
		goBegin(itV);
		for(;!atEnd(itV);++itV) {
			if ((value(itV) == bState) || (value(itV) == eState)) continue;
			if (!isSilent(hmm, value(itV))) continue;
			TCargo maxValue = infVal;
			TVertexDescriptor maxVertex = nilVertex;

			// Find maximum
			// A vertex iterator guarantees that the vertices are processend in increasing order
			// That is, the smallest silent states come first!!!
			TVertexIterator itMax(hmm);
			for(;!atEnd(itMax);++itMax) {
				if ((!isSilent(hmm, value(itMax))) ||
					((isSilent(hmm, value(itMax))) && (value(itMax) < value(itV)))) {
						TCargo local = value(vMat, i * numRows + value(itMax)) + std::log( (double) getTransitionProbability(hmm, *itMax, *itV));
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
	TCargo maxValue = infVal;
	TVertexDescriptor maxVertex = 0;
	TVertexIterator itMax(hmm);
	for(;!atEnd(itMax);++itMax) {
		TCargo local = value(vMat, len * numRows + *itMax) + std::log( (double) getTransitionProbability(hmm, *itMax, eState));
		if (local > maxValue) {
			maxValue = local;
			maxVertex = *itMax;
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

	return (TCargo) value(vMat, (len+1) * numRows + eState);
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
..returns:TCargo
...remarks:Probability of the sequence.
..see:Function.viterbiAlgorithm
..see:Function.backwardAlgorithm
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TSequence>
inline TCargo
forwardAlgorithm(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& hmm,
				 TSequence const& seq)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	String<TCargo> fMat;
	TSize numCols = length(seq) + 2;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	fill(fMat, numCols * numRows, 0);
	value(fMat, getBeginState(hmm)) = (TCargo) 1;
	TVertexDescriptor eState = getEndState(hmm);
	TSize scaling = 10;

	// Recurrence
	TSize len = length(seq);
	for(TSize i=1; i<=len; ++i) {
		TVertexIterator itV(hmm);
		for(;!atEnd(itV);++itV) {
			TCargo sum = 0;
			TVertexIterator itAll(hmm);
			for(;!atEnd(itAll);++itAll) sum += value(fMat, (i-1) * numRows + *itAll) * getTransitionProbability(hmm, *itAll, *itV);
			value(fMat, i * numRows + *itV) = getEmissionProbability(hmm, *itV, seq[i-1]) * sum * scaling;
		}
	}

	// Termination
	TCargo sum = 0;
	TVertexIterator itAll(hmm);
	for(;!atEnd(itAll);++itAll) {
		sum += value(fMat, len * numRows + *itAll) * getTransitionProbability(hmm, *itAll, eState);
	}
	value(fMat, (len+1) * numRows + eState) = sum;

	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(fMat, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return (TCargo) (value(fMat, (len+1) * numRows + eState) / ( (double) std::pow( (double) scaling, (double) len)));
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
..returns:TCargo
...remarks:Probability of the sequence.
..see:Function.viterbiAlgorithm
..see:Function.forwardAlgorithm
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TSequence>
inline TCargo
backwardAlgorithm(Graph<Hmm<TAlphabet, TCargo, TSpec> > const& hmm,
				  TSequence const& seq)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	String<TCargo> bMat;
	TSize numCols = length(seq) + 1;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	fill(bMat, numCols * numRows, 0);
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);
	TSize len = length(seq);
	TVertexIterator itAll(hmm);
	TSize scaling = 10;
	for(;!atEnd(itAll);++itAll) value(bMat, len * numRows + *itAll) = getTransitionProbability(hmm, *itAll, eState) * scaling;
	
	// Recurrence
	for(TSize i=len - 1; i>0; --i) {
		TVertexIterator itV(hmm);
		for(;!atEnd(itV);++itV) {
			TCargo sum = 0;
			TVertexIterator itAll(hmm);
			for(;!atEnd(itAll);++itAll) sum += value(bMat, (i+1) * numRows + *itAll) * getTransitionProbability(hmm, *itV, *itAll) * getEmissionProbability(hmm, *itAll, seq[i]);
			value(bMat, i * numRows + *itV) =  sum * scaling;
		}
	}

	// Termination
	TCargo sum = 0;
	goBegin(itAll);
	for(;!atEnd(itAll);++itAll) {
		sum += value(bMat, 1 * numRows + *itAll) * getTransitionProbability(hmm, bState, *itAll) * getEmissionProbability(hmm, *itAll, seq[0]);
	}
	value(bMat, bState) = sum;

	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(bMat, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return (TCargo) (value(bMat, bState) / ( (double) std::pow( (double) scaling, (double) len)));
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
