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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Library generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   LocalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TScore>::Type TScoreValue;

	// Clear graph
	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString, Rooted>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TScoreValue, Block<> > score_values;

	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			// Pairwise alignment graph
			TStringSet pairSet;
			assignValueById(pairSet, str, positionToId(str, i));
			assignValueById(pairSet, str, positionToId(str, j));
			typedef Graph<Alignment<TStringSet, unsigned int> > TPairGraph;
			typedef typename VertexDescriptor<TPairGraph>::Type TVD;
			typedef typename Iterator<TPairGraph, EdgeIterator>::Type TEI;
			TPairGraph pGraph(pairSet);
			localAlignment(pGraph, score_type, SmithWatermanClump() );
			
			// Remember the matches and their scores
			TEI it(pGraph);
			for(;!atEnd(it);++it) {
				TVD sV = sourceVertex(it);
				TVD tV = targetVertex(it);
				push_back(matches, TFragment( (unsigned int) sequenceId(pGraph, sV), (unsigned int) fragmentBegin(pGraph,sV), (unsigned int) sequenceId(pGraph, tV),  (unsigned int)  fragmentBegin(pGraph,tV),  (unsigned int)  fragmentLength(pGraph,tV)));
				push_back(score_values, cargo(*it));
			}
		}
	}

	// Refine all matches, rescore matches and create multiple alignment
	matchRefinement(matches,str,const_cast< TScore&>(score_type),g);

	// Adapt edge weights, scale weights by significance of local match
	TFragmentStringIter endIt = end(matches);
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			TVertexDescriptor p1 = findVertex(g, id1, pos1);
			TVertexDescriptor p2 = findVertex(g, id2, pos2);
			TEdgeDescriptor e = findEdge(g, p1, p2);
			if (e != 0) cargo(e) *= (TCargo) getValue(score_values, position(it));
			else addEdge(g, p1, p2, (TCargo) getValue(score_values, position(it)));
			pos1 += fragmentLength(g, p1);
			pos2 += fragmentLength(g, p2);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TDistanceMatrix, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TDistanceMatrix& dist,
					   TScore const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;

	// Clear graph
	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	resize(dist, nseq * nseq);
	
	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	// All pairwise alignments
	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			// Pairwise alignment, get the matches
			TStringSet pairSet;
			assignValueById(pairSet, str, positionToId(str, i));
			assignValueById(pairSet, str, positionToId(str, j));
			TSize from = length(matches);

			// Ends free-space alignment or not?
			//TSize len1 = length(pairSet[0]);
			//TSize len2 = length(pairSet[1]);
			//if (!((double) len1 / (double) len2 > 0.80)) {
			//	globalAlignment(matches, pairSet, score_type, AlignConfig<false,true,true,false>(), Gotoh() );
			//} else if (!((double) len1 / (double) len2 < 1.20)) {
			//	globalAlignment(matches, pairSet, score_type, AlignConfig<true,false,false,true>(), Gotoh() );
			//} else {
			//	globalAlignment(matches, pairSet, score_type, Gotoh() );
			//}
			globalAlignment(matches, pairSet, score_type, Gotoh() );
			
			// Determine a sequence weight
			TSize alignLen = 0;
			double seqSim = getSequenceSimilarity(matches, pairSet, from, length(matches), alignLen, typename Value<TStringSet>::Type() );

			// Normalize sequence similarity by alignment length
			double normalizedSimilarity = (double) seqSim / (double) alignLen;

			// Assign the values
			assignValue(dist, i*nseq+j, 1 - normalizedSimilarity);
			assignValue(dist, j*nseq+i, 1 - normalizedSimilarity);	
		}
	}

	// Refine all matches, rescore the matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;

	// Clear graph
	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	
	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString, Rooted>::Type TFragmentStringIter;
	TFragmentString matches;

	// All pairwise alignments
	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			// Pairwise alignment, get the matches
			TStringSet pairSet;
			assignValueById(pairSet, str, positionToId(str, i));
			assignValueById(pairSet, str, positionToId(str, j));

			// Ends free-space alignment or not?
			//TSize len1 = length(pairSet[0]);
			//TSize len2 = length(pairSet[1]);
			//if (!((double) len1 / (double) len2 > 0.80)) {
			//	globalAlignment(matches, pairSet, score_type, AlignConfig<false,true,true,false>(), Gotoh() );
			//} else if (!((double) len1 / (double) len2 < 1.20)) {
			//	globalAlignment(matches, pairSet, score_type, AlignConfig<true,false,false,true>(), Gotoh() );
			//} else {
			//	globalAlignment(matches, pairSet, score_type, Gotoh() );
			//}
			//globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
			globalAlignment(matches, pairSet, score_type, Gotoh() );
		}
	}

	// Refine all matches, rescore the matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList>
inline void 
selectPairsForLibraryGeneration(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
								TPairList& pList)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TPairList>::Type TPair;
	TStringSet& str = stringSet(g);
	TSize nseq = length(str);
	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			appendValue(pList, TPair(positionToId(str, i), positionToId(str, j)));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TScore, typename TDistanceMatrix>
inline void 
testGeneratePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						   TPairList const& pList,
						   TScore const& score_type,
						   TDistanceMatrix& dist,
						   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Value<TDistanceMatrix>::Type TValue;
	typedef typename Value<TScore>::Type TScoreValue;

	// Clear library graph and guide tree
	clearVertices(g);

	// Pairwise alignments for all pairs in the list
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	fill(dist, nseq * nseq, 0);
		
	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	// Distance graph
	Graph<Undirected<void> > topologicalGraph;
	for(TSize i=0;i<nseq; ++i) addVertex(topologicalGraph);

	// All pairwise alignments
	TSize amountOfPairs = length(pList);
	for(TSize k=0; k<amountOfPairs; ++k) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = (pList[k]).i1;
		TId id2 = (pList[k]).i2;
		assignValueById(pairSet, str, id1);
		assignValueById(pairSet, str, id2);
		
		// Overlap alignment
		TFragmentString pairwise_matches;
		globalAlignment(pairwise_matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
			
		// Determine a sequence weight
		TSize matchLen = 0;
		TSize overlapLen = getOverlapLength(pairwise_matches, pairSet, matchLen, typename Value<TStringSet>::Type() );
		double quality = (double) matchLen / (double) overlapLen;

		// Get only the good overlap alignments
		if (((quality >= 1) && (matchLen >= 8)) ||
			((quality >= 0.8) && (matchLen >= 15))) {
			TSize len = length(pairwise_matches);	
			for(TSize z = 0; z<len; ++z) {
				//TGraph tmp(pairSet);
				//globalAlignment(tmp, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
				//std::cout << "Match length: " << matchLen << std::endl;
				//std::cout << "Overlap length: " << overlapLen << std::endl;
				//std::cout << "Quality: " << quality << std::endl;
				//std::cout << tmp << std::endl;

				push_back(matches, pairwise_matches[z]);
			}
			// Create a corresponding edge
			TSize i = idToPosition(str, id1);
			TSize j = idToPosition(str, id2);
			if (i<j) assignValue(dist, i*nseq+j, matchLen * quality);
			else assignValue(dist, j*nseq+i, matchLen * quality);	
			addEdge(topologicalGraph, i, j);
		}
	}
	
	// Refine all matches, rescore the matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
	//matchRefinement(matches,str,g);

	// Eliminate non-cliques of 3
	for(TSize i=0;i<nseq;++i) {
		for(TSize j=0;j<i;++j) {
			if (getValue(dist, j*nseq+i) == 0) continue;
			bool foundTriplet = false;
			for(TSize k=0;k<j;++k) {
				if ((getValue(dist, k*nseq+j) != 0) &&
					(getValue(dist, k*nseq+i) != 0)) {
					foundTriplet = true;
					break;
				}
			}
			if (!foundTriplet) {
				for(TSize k=j+1;k<nseq;++k) {
					if (getValue(dist, j*nseq+k) != 0) {
						if (((i < k) && (getValue(dist, i*nseq+k) != 0)) ||
							((i > k) && (getValue(dist, k*nseq+i) != 0))) {
								foundTriplet = true;
								break;
						}
					}
				}
			}
			if (!foundTriplet) {
				//std::cout << "Eliminate pair: " << j << ',' << i << std::endl;
				value(dist, j*nseq+i) = 0;
			}
		}
		for(TSize j=i+1;j<nseq;++j) {
			if (getValue(dist, i*nseq+j) == 0) continue;
			bool foundTriplet = false;
			for(TSize k=0;k<j;++k) {
				if (getValue(dist, k*nseq+j) != 0) {
					if (((i < k) && (getValue(dist, i*nseq+k) != 0)) ||
						((i > k) && (getValue(dist, k*nseq+i) != 0))) {
							foundTriplet = true;
							break;
					}
				}
			}
			if (!foundTriplet) {
				for(TSize k=j+1;k<nseq;++k) {
					if ((getValue(dist, j*nseq+k) != 0) &&
						(getValue(dist, i*nseq+k) != 0)) {
						foundTriplet = true;
						break;
					}
				}
			}
			if (!foundTriplet) {
				//std::cout << "Eliminate pair: " << i << ',' << j << std::endl;
				value(dist, i*nseq+j) = 0;
			}
		}
	}

	// Prefer highly connected hubs
	TValue max_value = 0;
	for(TSize i=0;i<nseq;++i) {
		TSize degI = degree(topologicalGraph, i);
		for(TSize j=i+1;j<nseq;++j) {
			TSize degJ = degree(topologicalGraph, j);
			if ((degI == 0) || (degJ ==0)) {
				std::cout << "Not connected." << std::endl;
				//exit(0);
			}
			if (value(dist, i*nseq+j) != 0) {
				value(dist, i*nseq+j) += ((degI + degJ) / 2);
				if (value(dist, i*nseq+j) > max_value) max_value = value(dist, i*nseq+j);
			}
		}
	}

	// Make distances
	for(TSize i=0;i<nseq;++i) {
		for(TSize j=i+1;j<nseq;++j) {
			if (value(dist, i*nseq+j) > 0) {
				value(dist, i*nseq+j) = (max_value - getValue(dist, i*nseq+j));
				value(dist, j*nseq+i) = (max_value - getValue(dist, i*nseq+j));
			} else {
				value(dist, i*nseq+j) = 1000000;
				value(dist, j*nseq+i) = 1000000;
			}
		}
	}
}



/*
//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TPairList, typename TScore, typename TDistanceMatrix>
inline void 
testGeneratePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
						   TPairList const& pList,
						   TScore const& score_type,
						   TDistanceMatrix& dist,
						   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;

	// Clear library graph and guide tree
	clearVertices(g);

	// Pairwise alignments for all pairs in the list
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	fill(dist, nseq * nseq, 1000000);
		
	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	// Distance graph
	Graph<Directed<void> > topologicalGraph;
	for(TSize i=0;i<nseq; ++i) addVertex(topologicalGraph);

	// All pairwise alignments
	TSize amountOfPairs = length(pList);
	for(TSize k=0; k<amountOfPairs; ++k) {
		// Make a pairwise string-set
		TStringSet pairSet;
		TId id1 = (pList[k]).i1;
		TId id2 = (pList[k]).i2;
		assignValueById(pairSet, str, id1);
		assignValueById(pairSet, str, id2);
		
		// Overlap alignment
		TFragmentString pairwise_matches;
		globalAlignment(pairwise_matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
			
		// Determine a sequence weight
		TSize matchLen = 0;
		TSize overlapLen = getOverlapLength(pairwise_matches, pairSet, matchLen, typename Value<TStringSet>::Type() );
		double quality = (double) matchLen / (double) overlapLen;

		// Get only the good overlap alignments
		if (((quality >= 1) && (matchLen >= 5)) ||
			((quality >= 0.8) && (matchLen >= 10))) {
			TSize len = length(pairwise_matches);	
			for(TSize z = 0; z<len; ++z) {
				//TGraph tmp(pairSet);
				//globalAlignment(tmp, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
				//std::cout << "Match length: " << matchLen << std::endl;
				//std::cout << "Overlap length: " << overlapLen << std::endl;
				//std::cout << "Quality: " << quality << std::endl;
				//std::cout << tmp << std::endl;

				push_back(matches, pairwise_matches[z]);
			}
			// Create a corresponding edge
			TSize i = idToPosition(str, id1);
			TSize j = idToPosition(str, id2);
			double distance = 100 - matchLen * quality;
			if (distance < 0) distance = 0;
			assignValue(dist, i*nseq+j, distance);
			assignValue(dist, j*nseq+i, distance);	
			if (i < j) addEdge(topologicalGraph, i, j);
			else addEdge(topologicalGraph, j, i);
		}
	}
	
	// Refine all matches, rescore the matches and create multiple alignment
	matchRefinement(matches,str,g);

	typedef typename VertexDescriptor<Graph<Directed<void> > >::Type TVertexDescriptor;
	

	//typedef typename Iterator<Graph<Directed<void> >, VertexIterator>::Type TVertexIterator;
	//TVertexIterator itV(topologicalGraph);
	//for(;!atEnd(itV);goNext(itV)) {
	//	if (degree(topologicalGraph, *itV) == 0) {
	//		std::cout << "Not connected." << std::endl;
	//		exit(0);
	//	}
	//}

	// Topological sort the reads
	String<TVertexDescriptor> order;
	//std::cout << topologicalGraph << std::endl;
	topological_sort(topologicalGraph, order);

	// Adapt the distance matrix
	typedef Iterator<String<TVertexDescriptor> >::Type TStringIterator;
	TStringIterator itOrder = begin(order);
	TStringIterator itOrderEnd = end(order);
	TSize max_degree = 0;
	TSize max_index = 0;
	unsigned int count = 0;
	while(itOrder != itOrderEnd) {
		TSize current_degree = degree(topologicalGraph, *itOrder);
		if (current_degree == 0) {
			std::cout << "Not connected." << std::endl;
			exit(0);
		} else {
			if (current_degree > max_degree) {
				max_degree = current_degree;
				max_index = count;
			}
		}
		goNext(itOrder);
		++count;
	}
	std::cout << max_index << std::endl;
	std::cout << max_degree << std::endl;
	
	count = 0;
	itOrder = begin(order);
	while(itOrder != itOrderEnd) {
		if (max_index > count) {
			for(TSize i=0;i<nseq;++i) {
				value(dist, i*nseq+(*itOrder)) += (max_index-count) * 10;
				value(dist, (*itOrder)*nseq+i) += (max_index-count) * 10;
			}
		} else {
			for(TSize i=0;i<nseq;++i) {
				value(dist, i*nseq+(*itOrder)) += (count - max_index) * 10;
				value(dist, (*itOrder)*nseq+i) += (count - max_index) * 10;
			}
		}
		++count;
		goNext(itOrder);
	}
}
*/


//////////////////////////////////////////////////////////////////////////////
// T-Coffee Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar, typename TString>
inline void
_parse_readSequenceData(TFile & file,
						TChar & c,
						TString& str)
{
	SEQAN_CHECKPOINT

	append(str, c);

	// Read word
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(str, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
_readLibrary(TFile & file,
			 Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Value<TFile>::Type TValue;
	
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;

	TValue c;
	bool seq1ToN = false;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	typedef std::pair<unsigned int, unsigned int> TSeqRes;
	typedef std::map<TSeqRes, unsigned int> TNodeMap;
	TNodeMap node_map;
	TWord seq1 = 0;
	TWord seq2 = 0;
	while (!_streamEOF(file)) {
		_parse_skipWhitespace(file,c);
		if (_streamEOF(file)) break;
		if (c == '#') {
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			seq1 = _parse_readNumber(file, c);
			seq2 = _parse_readNumber(file, c);
			if (empty(g))
				if ((seq1 != 0) && (seq2 != 0)) seq1ToN = true;
			if (seq1ToN) {
				--seq1;
				--seq2;
			}
		} else if (c == '!') {
			_parse_skipLine(file, c);
		} else {
			unsigned int res1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			unsigned int res2 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			unsigned int weight = _parse_readNumber(file, c);
			_parse_skipLine(file,c);
		
			// Insert new vertex if necessary
			--res1;
			--res2;
			bool newEdge = false;
			TSeqRes key = std::make_pair(seq1, res1);
			TNodeMap::iterator nodePos = node_map.find(key);
			TId id1;
			if (nodePos == node_map.end()) {
				id1 = addVertex(g, seq1, res1, 1); 
				node_map.insert(std::make_pair(key, id1));
				newEdge = true;
			} else {
				id1 = nodePos->second;
			}

			key = std::make_pair(seq2, res2);
			nodePos = node_map.find(key);
			TId id2;
			if (nodePos == node_map.end()) {
				id2 = addVertex(g, seq2, res2, 1); 
				node_map.insert(std::make_pair(key, id2));
				newEdge = true;
			} else {
				id2 = nodePos->second;
			}

			// Insert a new edge or adapt the weight
			if (newEdge) addEdge(g,id1,id2,weight);
			else {
				TEdgeDescriptor e = findEdge(g,id1,id2);
				if( e == 0 ) addEdge(g,id1,id2,weight);
				else cargo(e) += weight;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
void 
read(TFile & file,
	 Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
	 TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Value<TStringSet>::Type TString;

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
	_parse_skipWhitespace(file, c);
	
	// Read number of sequences
	TWord nSeq = (TWord) _parse_readNumber(file, c);

	// Read sequences
	for(TWord i=0; i<nSeq; ++i) {
		_parse_skipWhitespace(file, c);
		std::cout << _parse_readIdentifier(file, c) << ", ";
		_parse_skipWhitespace(file, c);
		std::cout << _parse_readNumber(file, c) << ", ";
		_parse_skipWhitespace(file, c);
		TString sequence;
		_parse_readSequenceData(file,c,sequence);
		std::cout << sequence << std::endl;
	}
	// Reinitialize the graph, because we changed the sequences
	clearVertices(g);

	// Read library
	_readLibrary(file,g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TSpec, typename TNames>
void 
read(TFile & file,
	 StringSet<TString, TSpec>& oriStr,
	 TNames& names,
	 TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
	_parse_skipWhitespace(file, c);
	
	// Read number of sequences
	TWord nSeq = (TWord) _parse_readNumber(file, c);
	resize(oriStr, nSeq);

	// Read sequences
	for(TWord i=0; i<nSeq; ++i) {
		_parse_skipWhitespace(file, c);
		appendValue(names, _parse_readIdentifier(file, c));
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readSequenceData(file,c,oriStr[i]);
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;\
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TStringSet>::Type TSize;

	_streamWrite(file, "! TC_LIB_FORMAT_01\n");
	TSize len = length(getStringSet(g));
	_streamPutInt(file, len);
	_streamPut(file, '\n');
	for(TSize i=0;i<len;++i) {
		_streamWrite(file, "seq");
		_streamPutInt(file, i);
		_streamPut(file, ' ');
		TString str = value(getStringSet(g), i);
		_streamPutInt(file, length(str));
		_streamPut(file, ' ');
		_streamWrite(file, str);
		_streamPut(file, '\n');
	}


	typedef std::pair<unsigned int, unsigned int> TSeq;
	typedef Triple<unsigned int, unsigned int, unsigned int> TData;
	typedef std::multimap<TSeq, TData> TMap;
	TMap m;

	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		if (sequenceId(g,sV) > sequenceId(g,tV)) {
			TVertexDescriptor tmp = sV;
			sV = tV;
			tV = tmp;
		}
		m.insert(std::make_pair(TSeq(sequenceId(g,sV), sequenceId(g,tV)), 
								TData(fragmentBegin(g,sV) + 1, fragmentBegin(g,tV) + 1, getCargo(*it))));
	}
	TSeq old;
	for(TMap::iterator pos = m.begin();pos!=m.end();++pos) {
		if (old != pos->first) {
			old = pos->first; 
			_streamPut(file, '#');
			_streamPutInt(file, pos->first.first);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->first.second);
			_streamPut(file, '\n');		
		}
		_streamPutInt(file, pos->second.i1);
		_streamPut(file, ' ');
		_streamPutInt(file, pos->second.i2);
		_streamPut(file, ' ');
		_streamPutInt(file, pos->second.i3);
		_streamPut(file, '\n');	
	}
	_streamWrite(file, "! SEQ_0_TO_N-1");
	_streamPut(file, '\n');	
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
