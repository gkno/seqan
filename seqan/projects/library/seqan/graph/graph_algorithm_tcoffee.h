#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_TCOFFEE_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_TCOFFEE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph: T-Coffee - Distance Matrix
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
template<typename TValue1, typename TValue2, typename TSize>
void
_hitToScoreMatrix(Matrix<TValue1>& mat, Matrix<TValue2>& score, TSize const nseq) {
	SEQAN_CHECKPOINT

	typedef Matrix<TValue1> THitMatrix;
	typedef typename Size<THitMatrix>::Type TMatrixSize;

	// Initialize the score matrix
	setDimension(score, 2);setLength(score, 0, nseq);setLength(score, 1, nseq);
	fill(host(score), nseq*nseq, 0.0);

	// Calculate the score
	TValue2 score0;
	for (TMatrixSize row=0;row<nseq;++row) {
		score0 = getValue(mat, row*nseq+row);
		for(TMatrixSize col=0;col<nseq;++col) {		
			if (row == col) continue;
			// T-Coffee Score Function: Why score0 * 3 * 10.0 + 0.5 ???
			assignValue(score, row*nseq+col, (TValue2) ((score0 - getValue(mat, min(row,col)*nseq+max(row,col))) / score0 * 3 * 10.0 + 0.5));
		}
	}
	for (TMatrixSize row=0;row<nseq;++row) {
		for(TMatrixSize col=row+1;col<nseq;++col) {
			assignValue(score, row*nseq+col, 100 - min(getValue(score,row*nseq+col), getValue(score, col*nseq+row)));
			assignValue(score, col*nseq+row, getValue(score, row*nseq+col));
		}
	}

	
	//for (unsigned row=0;row<nseq;++row) {
	//	for(unsigned col=0;col<nseq;++col) {
	//		std::cout << getValue(score, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TValue>
void
getScoringMatrix(StringSet<String<Dna>, TSpec> const& strSet, Matrix<TValue>& score) {
	SEQAN_CHECKPOINT

	typedef unsigned int TWord;
	typedef StringSet<String<Dna>, TSpec> TStringSet;
	typedef typename Size<TStringSet>::Type TStringSetSize;
	typedef StringSet<String<Dna>, ConcatVirtual<> > TStringSetDna;
	typedef Index<TStringSetDna, Index_QGram<FixedShape<6> > > TIndex; // q-gram length = 6
	
	// Number of sequences
	TStringSetSize nseq = length(strSet);

	// Initialization
	Matrix<TWord> mat;   // Matrix for common k-tupels between sequence i and j
	setDimension(mat, 2);setLength(mat, 0, nseq);setLength(mat, 1, nseq);
	fill(host(mat), nseq*nseq, 0);

	// Index
	TIndex index;
	resize(indexText(index), nseq);
	for(TStringSetSize k=0;k<length(strSet);++k) indexText(index)[k] = strSet[k];
	indexCreate(index, QGram_SA());

	String<TWord> counter; // Counter for each sequence
	resize(counter, nseq);
	
	TWord nqgrams = length(indexDir(index)) - 1;
	if (nqgrams > 0) {
		TWord j1;
		TWord j2 = indexDir(index)[0];
		for(TWord i = 0; i < nqgrams; ++i) {  // Iterate over all q-grams
			j1 = j2;
			j2 = indexDir(index)[i+1];
			if (j1 == j2) continue;   //Empty hit-list for this q-gram

			// Clear the counters
			arrayFill(begin(counter, Standard()), end(counter, Standard()), 0);
			for(TWord j = j1; j < j2; ++j) ++counter[getValueI1(indexSA(index)[j])];

			// Add-up the values in the matrix
			for(TWord k = 0; k < nseq; ++k)
				for(TWord k2 = k; k2 < nseq; ++k2) {
					TWord minVal = counter[k];
					if (counter[k2] < minVal) minVal = counter[k2];
					assignValue(mat, k*nseq+k2, getValue(mat, k*nseq+k2) + minVal);
				}
		}
	}

	// Transform the hit matrix into a score matrix
	_hitToScoreMatrix(mat, score, nseq);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TValue>
void
getScoringMatrix(StringSet<TString, TSpec> const& strSet, Matrix<TValue>& score) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef StringSet<TString, TSpec> TStringSet;
	typedef typename Size<TStringSet>::Type TStringSetSize;
	typedef Matrix<TWord> THitMatrix;
	typedef typename Size<THitMatrix>::Type TMatrixSize;

	// Number of sequences
	TStringSetSize nseq = length(strSet);

	// Initialization
	THitMatrix mat;   // Matrix for common k-tupels between sequence i and j
	setDimension(mat, 2);setLength(mat, 0, nseq);setLength(mat, 1, nseq);
	fill(host(mat), nseq*nseq, 0);

	// StringSet where each string is a sequence of amino acid groups identifiers
	typedef StringSet<String<AAGroupsDayhoff>, ConcatVirtual<> > TStringSetAA;
	// q-gram length = 6
	typedef Index<TStringSetAA, Index_QGram<FixedShape<6> > > TIndex;
	TIndex index;
	resize(indexText(index), nseq);
	
	// Recode the strings into amino acid groups
	for(TStringSetSize k=0;k<length(strSet);++k) indexText(index)[k] = strSet[k];

	// Build index	
	indexCreate(index, QGram_SA());

	// Counter for each sequence
	String<TWord> counter;
	resize(counter, nseq);
	
	TWord nqgrams = length(indexDir(index)) - 1;
	if (nqgrams > 0) {
		TWord j1;
		TWord j2 = indexDir(index)[0];
		for(TWord i = 0; i < nqgrams; ++i) {  // Iterate over all q-grams
			j1 = j2;
			j2 = indexDir(index)[i+1];
			if (j1 == j2) continue;   //Empty hit-list for this q-gram

			// Clear the counters
			arrayFill(begin(counter, Standard()), end(counter, Standard()), 0);
			for(TWord j = j1; j < j2; ++j) ++counter[getValueI1(indexSA(index)[j])];

			// Add-up the values in the matrix
			for(TWord k = 0; k < nseq; ++k)
				for(TWord k2 = k; k2 < nseq; ++k2) {
					TWord minVal = counter[k];
					if (counter[k2] < minVal) minVal = counter[k2];
					assignValue(mat, k*nseq+k2, getValue(mat, k*nseq+k2) + minVal);
				}
		}
	}

	// Transform the hit matrix into a score matrix
	_hitToScoreMatrix(mat, score, nseq);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TScoreValue, typename TScoreSpec, typename TSimValue, typename TSimSpec, typename TVal>
void
scoreToSimilarityMatrix(Matrix<TScoreValue, TScoreSpec>& score, Matrix<TSimValue, TSimSpec>& sim, TVal maxSim) {
	SEQAN_CHECKPOINT

	typedef Matrix<TScoreValue, TScoreSpec> TScoreMatrix;
	typedef Matrix<TSimValue, TSimSpec> TSimMatrix;
	typedef typename Size<TScoreMatrix>::Type TSize;
	typedef typename Value<TSimMatrix>::Type TValue;

	TSize nseq = length(score, 0);
	setDimension(sim, 2);setLength(sim, 0, nseq);setLength(sim, 1, nseq);
	fill(host(sim), nseq*nseq, 0);

	for (TSize row=0;row<nseq;++row) {
		for(TSize col=row;col<nseq;++col) {
			if (row == col) assignValue(sim, row*nseq+col, maxSim);
			else {
				assignValue(sim, row*nseq+col, (TValue) getValue(score, row*nseq+col));
				assignValue(sim, col*nseq+row, (TValue) getValue(score, row*nseq+col));
			}
		}
	}
		
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(sim, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TScoreValue, typename TScoreSpec, typename TDistValue, typename TDistSpec, typename TVal>
void
scoreToDistanceMatrix(Matrix<TScoreValue, TScoreSpec>& score, Matrix<TDistValue, TDistSpec>& dist, TVal maxDist) {
	SEQAN_CHECKPOINT

	typedef Matrix<TScoreValue, TScoreSpec> TScoreMatrix;
	typedef Matrix<TDistValue, TDistSpec> TDistMatrix;
	typedef typename Size<TScoreMatrix>::Type TSize;
	typedef typename Value<TDistMatrix>::Type TValue;
	
	TSize nseq = length(score, 0);
	setDimension(dist, 2);setLength(dist, 0, nseq);setLength(dist, 1, nseq);
	fill(host(dist), nseq*nseq, 0);

	for (TSize row=0;row<nseq;++row) {
		for(TSize col=row;col<nseq;++col) {
			if (row == col) assignValue(dist, row*nseq+col, 0);
			else {
				assignValue(dist, row*nseq+col, (maxDist - (TValue) getValue(score, row*nseq+col)));
				assignValue(dist, col*nseq+row, (maxDist - (TValue) getValue(score, row*nseq+col)));
			}
		}
	}
	
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(dist, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix, typename TVal>
void
normalizeMatrix(TMatrix& mat, TVal const max, TVal const norm) {
	SEQAN_CHECKPOINT
	typedef typename Size<TMatrix>::Type TSize;
	
	TSize nseq = length(mat, 0);

	for (TSize row=0; row<nseq; ++row)
		for (TSize col=0; col<nseq; ++col)
	       assignValue(mat, row*nseq+col, (getValue(mat,row*nseq+col) * norm ) / max);

	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrixSpec, typename TCargo, typename TSpec>
void
slowNjTree(Matrix<double, TMatrixSpec>& mat, Graph<Tree<TCargo, TSpec> >& g) {
	SEQAN_CHECKPOINT
	
	typedef typename Size<Matrix<double, TMatrixSpec> >::Type TSize;
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	
	TSize nseq = length(mat, 0);

	//for(TSize i=0;i<nseq;++i) {
	//	for(TSize j=0;j<nseq;++j) {
	//		std::cout << getValue(mat, i*nseq+j) << ",";
	//	}
	//	std::cout << std::endl;
	//}

	// First initialization
	clear(g);
	String<double> av;    // Average branch length to a combined node
	fill(av,nseq,0.0);

	String<TVertexDescriptor> connector;   // Nodes that need to be connected
	resize(connector, nseq);

	for(TSize i=0;i<nseq;++i) {
		addVertex(g);  // Add all the nodes that correspond to sequences
		assignValue(connector, i, i);
		assignValue(mat, i*nseq+i, 0.0);
	}

	// Main cycle
	double fnseqs=(double) nseq;
	for(TSize nc=0; nc<(nseq-3); ++nc) {
		double sumOfBranches = 0.0;

		// Copy upper triangle matrix to lower triangle
		for(TSize col=1; col<nseq; ++col) {
			for(TSize row=0; row<col; ++row) {
				assignValue(mat, col*nseq+row, getValue(mat, row*nseq+col));
				// Determine the sum of all branches
				sumOfBranches = sumOfBranches + getValue(mat, row*nseq+col);
			}
		}

		// Compute the sum of branch lengths for all possible pairs
		double tmin = 0.0;	
		TSize mini = 0;  // Next pair of seq i and j to join
		TSize minj = 0;
		for(TSize col=1; col<nseq; ++col)  {
			if (getValue(connector,col) != nilVertex) {
				for(TSize row=0; row<col; ++row) {
					if (getValue(connector,row) != nilVertex) {
						double diToAllOthers = 0.0;
						double djToAllOthers = 0.0;
						
						for(TSize i=0; i<nseq; ++i) {
							diToAllOthers += getValue(mat, i*nseq+row);
							djToAllOthers += getValue(mat, i*nseq+col);
						}

						double dij = getValue(mat, row*nseq+col);
						double total = diToAllOthers + djToAllOthers + (fnseqs - 2.0)*dij +2.0*(sumOfBranches - diToAllOthers - djToAllOthers);
						total /= (2.0*(fnseqs - 2.0));

						if ((tmin == 0) || (total < tmin)) {
							tmin = total;
							mini = row;
							minj = col;
						}
					}
				}
			}
		}

		// Print nodes that are about to be joined
		//std::cout << mini << std::endl;
		//std::cout << minj << std::endl;
		//std::cout << tmin << std::endl;
		//std::cout << std::endl;
		
		// Compute branch lengths
		double dMinIToOthers = 0.0;
		double dMinJToOthers = 0.0;
		for(TSize i=0; i<nseq; ++i) {
			dMinIToOthers += getValue(mat, i*nseq + mini);
			dMinJToOthers += getValue(mat, i*nseq + minj);
		}
		double dmin = getValue(mat, mini*nseq + minj);
		dMinIToOthers = dMinIToOthers / (fnseqs - 2.0);
		dMinJToOthers = dMinJToOthers / (fnseqs - 2.0);
		double iBranch = (dmin + dMinIToOthers - dMinJToOthers) * 0.5;
		double jBranch = dmin - iBranch;
		iBranch -= av[mini];
		jBranch -= av[minj];
		
		// Set negative branch length to zero
		if( fabs(iBranch) < 0.0001) iBranch = 0.0;
		if( fabs(jBranch) < 0.0001) jBranch = 0.0;
	
		// Print branch lengths
		//std::cout << iBranch << std::endl;
		//std::cout << jBranch << std::endl;
		//std::cout << std::endl;
		
		// Build tree
		TVertexDescriptor internalVertex = addVertex(g);
		addEdge(g, internalVertex, getValue(connector, mini), (TCargo) iBranch);
		addEdge(g, internalVertex, getValue(connector, minj), (TCargo) jBranch);

		// Remember the average branch length for the new combined node
		// Must be subtracted from all branches that include this node
		if(dmin <= 0.0) dmin = 0.000001;
		av[mini] = dmin * 0.5;


		// Re-initialisation
		// mini becomes the new combined node, minj is killed
		fnseqs = fnseqs - 1.0;
		assignValue(connector, minj, nilVertex);
		assignValue(connector, mini, internalVertex);

		for(TSize j=0; j<nseq; ++j) {
			if( getValue(connector, j) != nilVertex ) {
				double minIminJToOther = ( getValue(mat, mini*nseq+j) + getValue(mat, minj*nseq+j)) * 0.5;
				// Use upper triangle
				if((TSize) mini < j) assignValue(mat, mini*nseq+j, minIminJToOther);
				if((TSize) mini > j) assignValue(mat, j*nseq+mini, minIminJToOther);
			}
		}
		for(TSize j=0; j<nseq; ++j) {
			assignValue(mat, minj*nseq+j, 0.0);
			assignValue(mat, j*nseq+minj, 0.0);
		}
	}

	// Only three nodes left

	// Find the remaining nodes
	String<TSize> l;
	fill(l,3,0);
	TSize count = 0;
	for(TSize i=0; i<nseq; ++i) {
		if(getValue(connector, i) != nilVertex) {
			l[count] = i;
			++count;
		}
	}

	// Remaining nodes
	//std::cout << l[0] << std::endl;
	//std::cout << l[1] << std::endl;
	//std::cout << l[2] << std::endl;
	//std::cout << std::endl;

	String<double> branch;
	resize(branch, 3);
	branch[0] = (getValue(mat, l[0]*nseq+l[1]) + getValue(mat, l[0]*nseq+l[2]) - getValue(mat, l[1]*nseq+l[2])) * 0.5;
	branch[1] = (getValue(mat, l[1]*nseq+l[2]) + getValue(mat, l[0]*nseq+l[1]) - getValue(mat, l[0]*nseq+l[2])) * 0.5;
	branch[2] =  (getValue(mat, l[0]*nseq+l[2]) + getValue(mat, l[0]*nseq+l[1]) - getValue(mat, l[1]*nseq+l[2])) * 0.5;
    
	branch[0] -= av[l[0]];
	branch[1] -= av[l[1]];
	branch[2] -= av[l[2]];

	// Print branch lengths
	//std::cout << branch[0] << std::endl;
	//std::cout << branch[1] << std::endl;
	//std::cout << branch[2] << std::endl;
	//std::cout << std::endl;
    
	// Reset tiny negative and positive branch lengths to zero
	if( fabs(branch[0]) < 0.0001) branch[0] = 0.0;
	if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
	if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;
    
	// Build tree
	TVertexDescriptor internalVertex = addVertex(g);
	addEdge(g, internalVertex, getValue(connector, l[0]), (TCargo) branch[0]);
	addEdge(g, internalVertex, getValue(connector, l[1]), (TCargo) branch[1]);
	addEdge(g, internalVertex, getValue(connector, l[2]), (TCargo) branch[2]);
	g.data_root = internalVertex;
}


//////////////////////////////////////////////////////////////////////////////
// T-Coffee Library Handling
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar, typename TStringSet, typename TCargo, typename TSpec>
inline typename Value<TStringSet>::Type&
_parse_readSequenceData(TFile & file,
						TChar & c,
						Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Value<TStringSet>::Type TString;

	// Read word
	TString* str = new TString(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(*str, c);
	}
	return *str;
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

	typedef std::pair<unsigned int, unsigned int> TSeqRes;
	typedef std::map<TSeqRes, unsigned int> TNodeMap;
	TNodeMap node_map;
	TWord seq1 = 0;
	TWord seq2 = 0;
	bool seq1ToN = false;
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		_parse_skipWhitespace(file,c);
		if (_streamEOF(file)) break;
		if (c == '#') {
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			seq1 = _parse_readNumber(file, c);
			seq2 = _parse_readNumber(file, c);
		} else if (c == '!') {
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			String<char> str = _parse_readIdentifier(file, c);
			_parse_skipLine(file, c);
			if (str == "SEQ_1_TO_N") seq1ToN = true;
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

	if (seq1ToN) {
		typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
		TConstIter it(g);
		for(;!atEnd(it);++it) {
			sequenceId(g, *it) -= 1;
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

	// Read sequences
	for(TWord i=0; i<nSeq; ++i) {
		_parse_skipWhitespace(file, c);
		std::cout << _parse_readIdentifier(file, c) << ", ";
		_parse_skipWhitespace(file, c);
		std::cout << _parse_readNumber(file, c) << ", ";
		_parse_skipWhitespace(file, c);
		TIdType id = assignValueById(stringSet(g), _parse_readSequenceData(file,c,g));
		std::cout << id << ", ";
		std::cout << getValueById(stringSet(g), id) << std::endl;
		SEQAN_ASSERT(id < nSeq)		
	}
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
