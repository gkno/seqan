#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Library generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
	
template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   Lcs_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	//typedef typename Value<TStringSet>::Type TString;
	typedef String<AAGroupsDayhoff> TString;
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

	for(TSize i=0; i<nseq-2; ++i) {
		for(TSize j=i+1; j<nseq-1; ++j) {
			for(TSize k=j+1; k<nseq; ++k) {

				TString str1 = str[i];
				TString str2 = str[j];
				TString str3 = str[k];

				// Lcs between first and second string
				String<std::pair<unsigned int, unsigned int>, Block<> > pos1;
				longestCommonSubsequence(str1, str2, pos1);

				// Add the third string
				TString intermediateStr;
				TSize lastPos1 = length(pos1) - 1;
				reserve(intermediateStr, lastPos1 + 1);
				for(int z = lastPos1; z>=0; --z) {
					appendValue(intermediateStr, (str1)[pos1[z].first]);
				}
				String<std::pair<unsigned int, unsigned int>, Block<> > pos2;
				longestCommonSubsequence(intermediateStr, str3, pos2);

				// Get the significant matches in all 3 sequences
				bool firstRun = true;
				TSize lenMatch = 0;TSize kBegin = 0;
				TSize iBegin = 0;TSize jBegin = 0;
				for(int z = length(pos2)-1; z>=0; --z) {
					if (firstRun) {
						// Where do the matches start in each sequence?
						firstRun = false;
						lenMatch = 1;
						kBegin = pos2[z].second;
						iBegin = pos1[lastPos1 - (pos2[z].first)].first;
						jBegin = pos1[lastPos1 - (pos2[z].first)].second;
					} else {
						//Is it a consecutive run of characters in all 3 sequences?
						if ((pos2[z+1].second + 1 == pos2[z].second) &&
							(pos1[lastPos1 - (pos2[z+1].first)].first + 1 == pos1[lastPos1 - (pos2[z].first)].first) &&
							(pos1[lastPos1 - (pos2[z+1].first)].second + 1 == pos1[lastPos1 - (pos2[z].first)].second)) {
								++lenMatch;
						} else {
							// A new match started, what about the old one?
							if (lenMatch > 1) {
								//// Debug code
								//typedef typename Infix<TString>::Type TInfix;
								//TInfix inf1 = infix(str1,iBegin, iBegin + lenMatch);
								//TInfix inf2 = infix(str2,jBegin, jBegin + lenMatch);
								//TInfix inf3 = infix(str3,kBegin, kBegin + lenMatch);
								//std::cout << inf1 << std::endl;
								//std::cout << inf2 << std::endl;
								//std::cout << inf3 << std::endl;

								push_back(matches, TFragment(positionToId(str, i),iBegin,positionToId(str, j),jBegin,lenMatch));
								push_back(matches, TFragment(positionToId(str, j),jBegin,positionToId(str, k),kBegin,lenMatch));
								push_back(matches, TFragment(positionToId(str, i),iBegin,positionToId(str, k),kBegin,lenMatch));
							}
							lenMatch = 1;
							kBegin = pos2[z].second;
							iBegin = pos1[lastPos1 - (pos2[z].first)].first;
							jBegin = pos1[lastPos1 - (pos2[z].first)].second;
						}
					}
				}
				// Process last match
				if (lenMatch > 1) {
					push_back(matches, TFragment(positionToId(str, i),iBegin,positionToId(str, j),jBegin,lenMatch));
					push_back(matches, TFragment(positionToId(str, j),jBegin,positionToId(str, k),kBegin,lenMatch));
					push_back(matches, TFragment(positionToId(str, i),iBegin,positionToId(str, k),kBegin,lenMatch));
				}
			}
		}
	}

	// Clear graph
	clearVertices(g);

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}


//////////////////////////////////////////////////////////////////////////////
	
template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   LocalTriple_Library)
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

	for(TSize i=0; i<nseq-2; ++i) {
		for(TSize j=i+1; j<nseq-1; ++j) {
			TSize iterations = 1;
			if (nseq < 50) iterations = 50 - nseq;
			for(TSize k=j+1; (k - (j+1) < iterations) && (k<nseq);++k) {
				TFragmentString local_matches;
				
				//First pair
				TStringSet pairSet;
				assignValueById(pairSet, str, positionToId(str, i));
				assignValueById(pairSet, str, positionToId(str, j));
				localAlignment(local_matches, pairSet, score_type, SmithWaterman() );

				// Second pair
				clear(pairSet);
				assignValueById(pairSet, str, positionToId(str, i));
				assignValueById(pairSet, str, positionToId(str, k));
				localAlignment(local_matches, pairSet, score_type, SmithWaterman() );

				// Third pair
				clear(pairSet);
				assignValueById(pairSet, str, positionToId(str, j));
				assignValueById(pairSet, str, positionToId(str, k));
				localAlignment(local_matches, pairSet, score_type, SmithWaterman() );

				// Refine matches
				clearVertices(g);

				// Refine all matches and create multiple alignment
				matchRefinement(local_matches,str,g);

				// Get significant local matches
				TEdgeIterator edgeIt(g);
				for(;!atEnd(edgeIt);++edgeIt) {
					TVertexDescriptor p1 = sourceVertex(edgeIt);
					TVertexDescriptor p2 = targetVertex(edgeIt);
					//std::cout << sequenceId(g, p1) << ',' << fragmentBegin(g,p1) << ',' << sequenceId(g, p2) << ',' << fragmentBegin(g,p2) << ',' << fragmentLength(g,p2) << std::endl;
					if ((outDegree(g, p1) > 1) && 
						(outDegree(g, p2) > 1) &&
						(fragmentLength(g, p1) > 1)) 
					{
							push_back(matches, TFragment(sequenceId(g,p1),fragmentBegin(g,p1),sequenceId(g,p2),fragmentBegin(g,p2),fragmentLength(g,p1)));
					}
				}
			}
		}
	}

	// Clear graph
	clearVertices(g);

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,const_cast<TScore&>(score_type),g);
}


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
	typedef typename Iterator<TFragmentString, Rooted>::Type TFragmentStringIter;
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
// T-Coffee Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar, typename TStringSet, typename TCargo, typename TSpec, typename TString>
inline void
_parse_readSequenceData(TFile & file,
						TChar & c,
						Graph<Alignment<TStringSet, TCargo, TSpec> >&,
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
				if ((seq1 != 0) || (seq2 != 0)) seq1ToN = true;
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

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
	_parse_skipWhitespace(file, c);
	
	// Read number of sequences
	TWord nSeq = (TWord) _parse_readNumber(file, c);
	resize(stringSet(g), nSeq);

	// Read sequences
	for(TWord i=0; i<nSeq; ++i) {
		_parse_skipWhitespace(file, c);
		std::cout << _parse_readIdentifier(file, c) << ", ";
		_parse_skipWhitespace(file, c);
		std::cout << _parse_readNumber(file, c) << ", ";
		_parse_skipWhitespace(file, c);
		_parse_readSequenceData(file,c,g,stringSet(g)[i]);
		//std::cout << id << ", ";
		std::cout << getValueById(stringSet(g), i) << std::endl;
		//SEQAN_ASSERT(id < nSeq)		
	}
	// Reinitialize the graph, because we changed the sequences
	clearVertices(g);

	// Read library
	_readLibrary(file,g);
}

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
								TData(segmentBegin(g,sV) + 1, segmentBegin(g,tV) + 1, getCargo(*it))));
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
