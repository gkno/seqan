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

	// Distinct between silent states and real states
	typedef String<TVertexDescriptor> TStateSet;
	typedef typename Iterator<TStateSet>::Type TStateIter;
	TStateSet silentStates;
	TStateSet realStates;
	TVertexIterator itVertex(hmm);
	for(;!atEnd(itVertex);goNext(itVertex)) {
		if (isSilent(hmm, value(itVertex))) { 
			appendValue(silentStates, value(itVertex));
		} else {
			appendValue(realStates, value(itVertex));
		}
	}

	// Initialization for silent states connected to the begin state
	TStateIter itSilentStateEnd = end(silentStates);
	for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
		if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
		TProbability maxValue = 0.0;
		TVertexDescriptor maxVertex = nilVertex;

		for(TStateIter itBelow = begin(silentStates); itBelow != itSilentState; goNext(itBelow)) {
			TProbability local = value(vMat, value(itBelow)) * getTransitionProbability(hmm, value(itBelow), value(itSilentState));
			if (local > maxValue) {
				maxValue = local;
				maxVertex = value(itBelow);
			}
		}

		// Set traceback vertex
		if (maxVertex != nilVertex) {
			value(vMat, value(itSilentState)) = maxValue;		
			value(traceback, value(itSilentState)) = maxVertex;
		}
	}

	// Recurrence
	TSize len = length(seq);
	for(TSize i=1; i<=len; ++i) {
		// Iterate over real states
		TStateIter itRealStateEnd = end(realStates);
		for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
			// Find maximum
			TProbability maxValue = 0.0;
			TVertexDescriptor maxVertex = nilVertex;
			TVertexIterator itMax(hmm);
			for(;!atEnd(itMax);++itMax) {
				TProbability local = value(vMat, (i-1) * numRows + value(itMax)) * getTransitionProbability(hmm, value(itMax), value(itRealState));
				if (local > maxValue) {
					maxValue = local;
					maxVertex = *itMax;
				}
			}
			// Set traceback vertex
			if (maxVertex != nilVertex) {
				value(vMat, i * numRows + value(itRealState)) = maxValue * getEmissionProbability(hmm, value(itRealState), value(seq, i-1));;
				value(traceback, i * numRows + value(itRealState)) = maxVertex;
			}
		}

		// Iterate over silent states
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
			// Find maximum
			TProbability maxValue = 0.0;
			TVertexDescriptor maxVertex = nilVertex;
	
			// Iterate over real states
			for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
				TProbability local = value(vMat, i * numRows + value(itRealState)) * getTransitionProbability(hmm, value(itRealState), value(itSilentState));
				if (local > maxValue) {
					maxValue = local;
					maxVertex = value(itRealState);
				}
			}
			// Iterate over silent states in increasing order
			for(TStateIter itBelow = begin(silentStates); itBelow != itSilentState; goNext(itBelow)) {
				TProbability local = value(vMat, i * numRows + value(itBelow)) * getTransitionProbability(hmm, value(itBelow), value(itSilentState));
				if (local > maxValue) {
					maxValue = local;
					maxVertex = value(itBelow);
				}
			}
			// Set traceback vertex
			if (maxVertex != nilVertex) {
				value(traceback, i * numRows + value(itSilentState)) = maxVertex;
				value(vMat, i * numRows + value(itSilentState)) = maxValue;		
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
	if (maxValue > 0.0) {
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
	}
	
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

	// Distinct between silent states and real states
	typedef String<TVertexDescriptor> TStateSet;
	typedef typename Iterator<TStateSet>::Type TStateIter;
	TStateSet silentStates;
	TStateSet realStates;
	TVertexIterator itVertex(hmm);
	for(;!atEnd(itVertex);goNext(itVertex)) {
		if (isSilent(hmm, value(itVertex))) { 
			appendValue(silentStates, value(itVertex));
		} else {
			appendValue(realStates, value(itVertex));
		}
	}

	// Initialization for silent states connected to the begin state
	TStateIter itSilentStateEnd = end(silentStates);
	for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
		if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
		TProbability sumValue = 0.0;
		for(TStateIter itBelow = begin(silentStates); itBelow != itSilentState; goNext(itBelow)) {
			sumValue += value(fMat, value(itBelow)) * getTransitionProbability(hmm, value(itBelow), value(itSilentState));
		}
		value(fMat, value(itSilentState)) = sumValue;		
	}

	// Recurrence
	TSize len = length(seq);
	for(TSize i=1; i<=len; ++i) {
		// Iterate over real states
		TStateIter itRealStateEnd = end(realStates);
		for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
			TProbability sum = 0.0;
			TVertexIterator itAll(hmm);
			for(;!atEnd(itAll);++itAll) sum += value(fMat, (i-1) * numRows + value(itAll)) * getTransitionProbability(hmm, value(itAll), value(itRealState));
			value(fMat, i * numRows + value(itRealState)) = getEmissionProbability(hmm, value(itRealState), seq[i-1]) * sum;
		}

		// Iterate over silent states
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
			TProbability sumValue = 0.0;
			// Iterate over real states
			for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
				sumValue += value(fMat, i * numRows + value(itRealState)) * getTransitionProbability(hmm, value(itRealState), value(itSilentState));
			}
			// Iterate over silent states in increasing order
			for(TStateIter itBelow = begin(silentStates); itBelow != itSilentState; goNext(itBelow)) {
				sumValue += value(fMat, i * numRows + value(itBelow)) * getTransitionProbability(hmm, value(itBelow), value(itSilentState));
			}
			value(fMat, i * numRows + value(itSilentState)) = sumValue;
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
	TSize numCols = length(seq) + 2;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	fill(bMat, numCols * numRows, 0.0);
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);
	TSize len = length(seq);
	value(bMat, (len + 1) * numRows + eState) = 1.0;
	
	// Distinct between silent states and real states
	typedef String<TVertexDescriptor> TStateSet;
	typedef typename Iterator<TStateSet>::Type TStateIter;
	TStateSet silentStates;
	TStateSet realStates;
	TVertexIterator itVertex(hmm);
	for(;!atEnd(itVertex);goNext(itVertex)) {
		if (isSilent(hmm, value(itVertex))) { 
			appendValue(silentStates, value(itVertex));
		} else {
			appendValue(realStates, value(itVertex));
		}
	}
	// Reverse the silent states order
	std::reverse(begin(silentStates), end(silentStates));

	// Initialization for silent states connected to the end state
	TStateIter itSilentStateEnd = end(silentStates);
	for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
		if (value(itSilentState) == eState) continue;
		TProbability sumValue = getTransitionProbability(hmm, value(itSilentState), eState);
		for(TStateIter itAbove = begin(silentStates); itAbove != itSilentState; goNext(itAbove)) {
			sumValue += value(bMat, len * numRows + value(itAbove)) * getTransitionProbability(hmm, value(itSilentState), value(itAbove));
		}
		value(bMat, len * numRows + value(itSilentState)) = sumValue;
	}

	// Initialization for real states
	TStateIter itRealStateEnd = end(realStates);
	for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
		TProbability sumValue = getTransitionProbability(hmm, value(itRealState), eState);
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			sumValue += value(bMat, len * numRows + value(itSilentState)) * getTransitionProbability(hmm, value(itRealState), value(itSilentState));
		}
		value(bMat, len * numRows + value(itRealState)) = sumValue;
	}

	
	// Recurrence
	if (len > 0) {
		for(TSize i=len - 1; i>0; --i) {
			// Iterate over silent states
			for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
				if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
				TProbability sumValue = 0.0;
				// Iterate over real states
				for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
					sumValue += value(bMat, (i+1) * numRows + value(itRealState)) * getTransitionProbability(hmm, value(itSilentState), value(itRealState))* getEmissionProbability(hmm, value(itRealState), value(seq, i));
				}
				// Iterate over silent states in decreasing order
				for(TStateIter itAbove = begin(silentStates); itAbove != itSilentState; goNext(itAbove)) {
					if ((value(itAbove) == bState) || (value(itAbove) == eState)) continue;
					sumValue += value(bMat, i * numRows + value(itAbove)) * getTransitionProbability(hmm, value(itSilentState), value(itAbove));
				}
				value(bMat, i * numRows + value(itSilentState)) = sumValue;
			}

			// Iteration over real states
			for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
				TProbability sumValue = 0.0;
				// Iterate over silent states
				for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
					if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
					sumValue += value(bMat, i * numRows + value(itSilentState)) * getTransitionProbability(hmm, value(itRealState), value(itSilentState));
				}
				// Iterate over real states
				for(TStateIter itR = begin(realStates); itR != itRealStateEnd; goNext(itR)) {
					sumValue += value(bMat, (i+1) * numRows + value(itR)) * getTransitionProbability(hmm, value(itRealState), value(itR)) * getEmissionProbability(hmm, value(itR), value(seq, i));
				}
				value(bMat, i * numRows + value(itRealState)) =  sumValue;
			}
		}
	
		// Termination
		// Iterate over silent states
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
			TProbability sumValue = 0.0;
			// Iterate over real states
			for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
				sumValue += value(bMat, 1 * numRows + value(itRealState)) * getTransitionProbability(hmm, value(itSilentState), value(itRealState))* getEmissionProbability(hmm, value(itRealState), value(seq, 0));
			}
			// Iterate over silent states in decreasing order
			for(TStateIter itAbove = begin(silentStates); itAbove != itSilentState; goNext(itAbove)) {
				if ((value(itAbove) == bState) || (value(itAbove) == eState)) continue;
				sumValue += value(bMat, value(itAbove)) * getTransitionProbability(hmm, value(itSilentState), value(itAbove));
			}
			value(bMat, value(itSilentState)) = sumValue;
		}
		// Sum up all values
		TProbability sumValue = 0.0;
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
			sumValue += value(bMat, value(itSilentState)) * getTransitionProbability(hmm, bState, value(itSilentState));
		}
		for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
			sumValue += value(bMat, 1 * numRows + value(itRealState)) * getTransitionProbability(hmm, bState, value(itRealState)) * getEmissionProbability(hmm, value(itRealState), value(seq, 0));
		}
		value(bMat, bState) = sumValue;
	}

	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(bMat, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return value(bMat, bState);
}


///////////////////////////////////////////////////////////////////////////////////////

/**
.Function.generateSequence:
..cat:Graph.Hmm
..summary:Generates random state and alphabet sequences of a given HMM.
...remarks:Because of silent states, generated alphabet and state sequences might have different length.
..signature:generateSequence(hmm, sequences, states, numSeq, maxLength)
..param.hmm:In-parameter:Input HMM.
...type:Spec.Hmm
..param.sequences:The StringSet of alphabet sequences.
...type:Class.StringSet
..param.sequences:The StringSet of state sequences.
...type:Class.StringSet
..param.numSeq:The number of sequences to generate.
..param.maxLength:The maximum length of the sequences.
...remarks:Sequences might be shorter if the end state is reached prior to maxLength.
..returns:void
*/
template<typename TAlphabet, typename TProbability, typename TSpec,typename TSequenceSet, typename TStateSeqSet, typename TSize>
inline void
generateSequence(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				 TSequenceSet& sequences,
				 TStateSeqSet& states,
				 TSize numSeq,
				 TSize maxLength) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TSequenceSet>::Type TSequence;
	typedef typename Value<TStateSeqSet>::Type TStateSeq;
	
	// Initialization
	mtRandInit();
	clear(sequences);
	clear(states);
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	
	// Simulate sequences
	TVertexDescriptor currentState;
	TVertexDescriptor endState = getEndState(hmm);
	for(TSize i=0;i<numSeq;++i){
		currentState = getBeginState(hmm);
		TSequence seq;
		TStateSeq stat;
		appendValue(stat, getBeginState(hmm));
		bool stop = false;
		TSize pos = 0; 
		while (pos < maxLength) {
			TProbability prob = (double) (mtRand() % 100) / (double) (100);
			TProbability compareProb = 0.0;
			TOutEdgeIterator itState(hmm, currentState);
			// Determine the next state
			for(;!atEnd(itState);++itState) {
				// Probability of the next transition
				compareProb += getTransitionProbability(hmm, value(itState));
				// Compare with random probability
				if (prob <= compareProb){
					TVertexDescriptor nextState = targetVertex(hmm, value(itState));
					if (nextState == endState) {
						stop = true;
						break;
					}
					appendValue(stat, nextState);
					if (!isSilent(hmm, nextState)) {
						compareProb =0.0;
						prob = (double) (mtRand() % 100) / (double) (100);
						for (TSize c=0;c<alphSize;++c){
							compareProb += getEmissionProbability(hmm,targetVertex(hmm, value(itState)), TAlphabet(c));
							if (prob <= compareProb) {
								appendValue(seq, TAlphabet(c));
								++pos;
								break;
							}
						}
					}
					currentState = nextState;
					break;
				}
			}
			if (stop==true) break;
		}
		appendValue(stat, getEndState(hmm));
		appendValue(sequences, seq);
		appendValue(states, stat);
	}
}


template<typename TAlphabet, typename TProbability, typename TSpec,typename TEmissionCounter, typename TTransitionCounter>
inline void 
__parameterEstimator(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm,
					 TEmissionCounter const& emission,
					 TTransitionCounter const& transition)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	
	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	TSize pseudoCount = 1;

	// Estimate the parameters from the counter values
	TVertexIterator itAll(hmm);
	for(;!atEnd(itAll);++itAll){
		if (!isSilent(hmm, value(itAll))) {
			TSize summedCount =0;
			for(TSize i=0; i<alphSize;++i) summedCount += (value(value(emission, value(itAll)),i) + pseudoCount);
			for(TSize i=0; i<alphSize;++i) emissionProbability(hmm,value(itAll),TAlphabet(i)) = ((double) (value(value(emission, value(itAll)),i) + pseudoCount) / (double) summedCount);
		}
		TSize summedCount =0;
		TOutEdgeIterator itOutSum(hmm,value(itAll));
		for(;!atEnd(itOutSum);++itOutSum) summedCount += getProperty(transition, value(itOutSum)) + pseudoCount;
		TOutEdgeIterator itOut(hmm, value(itAll));
		for(;!atEnd(itOut);++itOut) transitionProbability(hmm, value(itOut)) = ((double) (getProperty(transition, value(itOut)) + pseudoCount) / (double) (summedCount));
	}
}



template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequenceSet, typename TStateSeqSet>
inline void 
estimationWithStates(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm,
					 TSequenceSet& sequences,
					 TStateSeqSet& states)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TStateSeqSet>::Type TStateSeq;
	typedef typename Value<TStateSeq>::Type TState;
	typedef typename Value<TSequenceSet>::Type TSequence;
	typedef typename Value<TSequence>::Type TChar;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	typedef String<TSize> TCountString;
	TCountString transitionCounter;
	fill(transitionCounter, getIdUpperBound(_getEdgeIdManager(hmm)), 0);
	StringSet<TCountString> emissionCounter;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	for(TSize i =0; i<numRows;++i) {
		TCountString emisCount;
		fill(emisCount,alphSize,0);
		appendValue(emissionCounter,emisCount);
	}
	
	// Iterate over all sequences
	for (TSize j=0;j<length(states);++j) {
		TSize posSeq = 0;
		for (TSize i = 0;i<length(value(states,j));++i) {
			TState s = value(value(states,j),i);
			if (!isSilent(hmm, s)) {
				TAlphabet c = value(value(sequences,j),posSeq);
				value(value(emissionCounter, s), ordValue(c)) += 1;
				++posSeq;
			}
			if (i<(length(value(states,j))-1)) {
				TEdgeDescriptor currEdge = findEdge(hmm, s, value(value(states,j), (i+1)));
				property(transitionCounter, currEdge) += 1;
			}
		}
	}

	//// Debug Code
	//typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	//TEdgeIterator itE(hmm);
	//for(;!atEnd(itE);goNext(itE)) std::cout << sourceVertex(itE) << ',' << targetVertex(itE) << ':' << property(transitionCounter, value(itE)) << std::endl;
	//for(TSize j=0; j<length(emissionCounter); ++j) {
	//	for(TSize i=0;i<length(value(emissionCounter, j));++i) {
	//		std::cout << j << ":" << TAlphabet(i) << '=' << value(value(emissionCounter, j), i) << std::endl;
	//	}
	//}

	// Estimate Parameters
	__parameterEstimator(hmm,emissionCounter, transitionCounter);
}


template<typename TAlphabet, typename TProbability, typename TSpec>
inline void
__fillHmmUniform(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	
	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	
	// Iterate over all states
	TVertexIterator itState(hmm);
	for(;!atEnd(itState);goNext(itState)) {		//pass through the states of the hmm	
		TSize oD = outDegree(hmm, value(itState));
		TOutEdgeIterator itOut(hmm, value(itState));
		for (;!atEnd(itOut);goNext(itOut)) transitionProbability(hmm, value(itOut)) = (TProbability) (1.0 / (double) (oD));
		if (!isSilent(hmm, value(itState))) {
			for(TSize i=0;i<alphSize;++i){
				emissionProbability(hmm,value(itState),TAlphabet(i)) = (TProbability) (1.0 / (double) (alphSize));
			}
		}
	}
}


template<typename TAlphabet, typename TProbability, typename TSpec>
inline void
__fillHmmRandom(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	
	// Initialization
	mtRandInit();
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	
	// Iterate over all states
	TVertexIterator itState(hmm);
	for(;!atEnd(itState);goNext(itState)) {		//pass through the states of the hmm	
		TSize oD = outDegree(hmm, value(itState));
		if (oD > 0) {
			String<TSize> counts;
			TSize sum = 0;
			for(TSize i = 0;i<oD;++i) {
				TSize rd = (mtRand() % 100) + 1;
				sum += rd;
				appendValue(counts, rd);
			}
			TSize pos = 0;
			TOutEdgeIterator itOut(hmm, value(itState));
			for (;!atEnd(itOut);goNext(itOut), ++pos) transitionProbability(hmm, value(itOut)) = (TProbability) ((double) (value(counts, pos)) / (double) (sum));
		}	
		if (!isSilent(hmm, value(itState))) {
			String<TSize> counts;
			TSize sum = 0;
			for(TSize i = 0;i<alphSize;++i) {
				TSize rd = (mtRand() % 100) + 1;
				sum += rd;
				appendValue(counts, rd);
			}
			for(TSize i=0;i<alphSize;++i) emissionProbability(hmm,value(itState),TAlphabet(i)) = (TProbability) ((double) (value(counts, i)) / (double) (sum));
		}
	}
}

template<typename TAlphabet, typename TProbability, typename TSpec>
inline void
randomizeHmm(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm)
{
	__fillHmmRandom(hmm);
	//__fillHmmUniform(hmm);
}


/*
template <typename TAlphabet, typename TCargo, typename TSpec, typename TSequence, typename TTag>
inline void 
__baumWelchAlgorithmMain(Graph<Hmm<TAlphabet, TCargo, TSpec > >& hmm,
				   StringSet<TSequence> const& seqSet,
				   TCargo border,
				   TCargo scalingFactor,
				   TTag const scalingTag)
{
	
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	
	std::fstream fstrm3;
	std::fstream fstrm4;
	std::fstream fstrm;
	std::fstream fstrm1;
	std::fstream fstrm2;
	fstrm.open("D:/Visio Profesional/Visio standard/Projects/Seqan_release_1.1/Test/hmms.txt", std::ios_base::out | std::ios_base::binary);
	fstrm1.open("D:/Visio Profesional/Visio standard/Projects/Seqan_release_1.1/Test/Emission.txt", std::ios_base::out | std::ios_base::binary);
	fstrm2.open("D:/Visio Profesional/Visio standard/Projects/Seqan_release_1.1/Test/Transition.txt", std::ios_base::out | std::ios_base::binary);
	fstrm3.open("D:/Visio Profesional/Visio standard/Projects/Seqan_release_1.1/Test/forw.txt", std::ios_base::out | std::ios_base::trunc);
	fstrm4.open("D:/Visio Profesional/Visio standard/Projects/Seqan_release_1.1/Test/backw.txt", std::ios_base::out | std::ios_base::trunc);
	//Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	String<TCargo> fMat;
	String<TCargo> bMat;
	TCargo fValue, em, tm;
	TCargo pseudo=1;
	TCargo localMax=0, oldLocalMax =0;
	String<TCargo> transProb;
	String<TCargo> estimatedEmission;
	StringSet<String<TCargo>> estimatedEmissionSet;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	for(TSize i =0; i<numRows;++i) appendValue(estimatedEmissionSet,estimatedEmission);
	TGraph backupHmm;

	for(TSize iter=0;iter<200;++iter){ //if no local maximum is found by the BW, it stops after it reached a certain iterationstep
		__logHmm(hmm, backupHmm, scalingFactor,TTag());
		fill(estimatedEmission,alphSize,0);
		for(TSize i =0; i<numRows;++i) value(estimatedEmissionSet,i)= estimatedEmission;	
		clear(transProb);
		fill(transProb, getIdUpperBound(_getEdgeIdManager(hmm)), 0);
				
		//MaximationStep
		//determine new contributions
		for (TSize i=0; i < length(seqSet);++i){		//sequences
			
			//compute forward und backward Algorithm
			fValue = __forwardAlgorithm(hmm, value(seqSet,i),fMat,border,scalingFactor, TTag());
			__backwardAlgorithm(hmm, value(seqSet,i),bMat,border ,scalingFactor, TTag());
		
			//compute logLikelihood of the current hmm
			__baumWelchTermination(hmm,fValue,scalingFactor,(TCargo)length(value(seqSet,i)), localMax, TTag());
			for (TSize j=0;j<length(value(seqSet,i));j++){  //passing through sequence
				TVertexIterator itAll(hmm);		
				for(;!atEnd(itAll);++itAll){	//all states
					
					//determine emission expactation values
					value(value(estimatedEmissionSet, *itAll),ordValue(value(value(seqSet,i),j)) ) += __baumWelchExpectedEmission(hmm,  fValue, scalingFactor, value(fMat, (j+1) * numRows + *itAll ), value(bMat, (j+1) * numRows + *itAll ),TTag());
					
					//determine transition expactation values
					TOutEdgeIterator itOut(hmm,*itAll);
					for(;!atEnd(itOut);++itOut){
						if (j==length(value(seqSet,i))-1){
							if (targetVertex(itOut)==endState(hmm)) {
								TCargo val = __baumWelchExpTransEState(hmm, fValue, scalingFactor, value(fMat, (j+1) * numRows + *itAll), getTransitionProbability(hmm,*itAll,endState(hmm)),TTag());
								assignProperty(transProb,*itOut,val);
							}
						}
						else { 
							TCargo val = getProperty(transProb,*itOut) + __baumWelchExpectedTransition(hmm, fValue, scalingFactor, value(fMat, j * numRows + *itAll), value(bMat, (j+1) * numRows + targetVertex(itOut) ), getTransitionProbability(hmm,*itOut), getEmissionProbability(hmm, targetVertex(itOut),value(value(seqSet,i),j) ),TTag());
							assignProperty(transProb,*itOut,val);
						}
					}
				}
			}	 
		}
		


		// Debug code

		//for(TSize i = 0; i<numRows; ++i) {
		//	for(TSize j=0; j<(4); ++j) {
		//		fstrm3 << value(fMat, j*numRows + i) << "	 ";
		//	}
		//	fstrm3<< std::endl;
		//}

		//for(TSize i = 0; i<numRows; ++i) {
		//	for(TSize j=0; j<3; ++j) {
		//		fstrm4 << value(bMat, j*numRows + i) << "	";
		//	}
		//	fstrm4 << std::endl;
		//}
		//
		////DEBUG
		//TVertexIterator itVert(hmm);
		//for(;!atEnd(itVert);++itVert){	
		//	for(TSize i = 0; i< alphSize;++i) fstrm1 << value(value(estimatedEmissionSet,*itVert),i)<<" ";
		//	fstrm1 <<std::endl;
		//	TOutEdgeIterator itVO(hmm,*itVert);
		//	for(;!atEnd(itVO);++itVO)fstrm2 << *itVert <<" -- " << targetVertex(itVO)<< " : " << getProperty(transProb, *itVO) << std::endl;
		//	fstrm2 << std::endl;
		//}
		//fstrm2 << std::endl;
		//fstrm1 << std::endl;
		//fstrm3 << std::endl;
		//fstrm4 << std::endl;
		//Compute new HMM parameters

		//retransfer the HMM
		hmm = backupHmm;

		//Estimation step 
		//calculate new model parameters
		parameterEstimator(hmm,estimatedEmissionSet, transProb);
		
		//DEBUG
		//std::cout << oldLocalMax  << std::endl;
		//std::cout << '(' <<localMax <<')';
		//std::cout << (localMax - oldLocalMax)  << std::endl;
		//std::cout << hmm << std::endl;
	
		write(fstrm,hmm);
		write(fstrm,"\n");
		
		std::cout<<".";
		//if ( (localMax - oldLocalMax)<0) std::cout << '(' << (localMax - oldLocalMax) << ") :"<<iter;
		
		
		//Termination by identical log likelihood of the models
		localMax /= length(seqSet);
		if ((iter != 0) && ((localMax - oldLocalMax) < 0.1) ) break; 
		else{
			oldLocalMax=localMax;
			localMax=0;
		}
		
	}
	fstrm.close();
	fstrm1.close();
	fstrm2.close();
	fstrm3.close();
	fstrm4.close();
	std::cout<<std::endl;
}
*/

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
