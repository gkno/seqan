#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_TCOFFEE_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_TCOFFEE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph: T-Coffee
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TCodedString, typename TTupelString, typename TAlphabetSize>
void
_getTupelString(TCodedString const& code, TTupelString& tupelString, TAlphabetSize const alphabet_size, unsigned int const ktup) {
	typedef typename Value<TCodedString>::Type TWord;
	typedef typename Size<TCodedString>::Type TSize;

	// Assign a unique number to each k-tupel
	String<TWord> prod;  // Scaling according to position in k-tupel
	resize(prod,ktup);
	for (TWord i=0; i<ktup;++i) prod[ktup-i-1]=(TWord)pow((double)alphabet_size,(double)i);

	TSize len = length(code);
	clear(tupelString);
	resize(tupelString, len-(ktup - 1)); 
	TSize tupelIndex = 0;
	TSize endTupel = 0;
	tupelString[tupelIndex] = 0;
	for(;endTupel<ktup;++endTupel) {
		tupelString[tupelIndex] += code[endTupel] * prod[endTupel];
	}
	++tupelIndex;
	for(;endTupel<len;++endTupel) {
		tupelString[tupelIndex] = tupelString[tupelIndex - 1];
		tupelString[tupelIndex] -= code[endTupel - ktup] * prod[0];
		tupelString[tupelIndex] *= alphabet_size;
		tupelString[tupelIndex] += code[endTupel];
		++tupelIndex;
	}
}


template<typename TInputString, typename TCodedString>
unsigned int 
_recodeSequence(TInputString const& input, TCodedString& out) {
	typedef typename Value<TCodedString>::Type TWord;
	typedef String<TInputString> TGroup;
	typedef typename Iterator<TInputString>::Type TStringIter;
	typedef typename Iterator<TGroup>::Type TGroupIter;
	
	// Predefined groups
	TGroup groups;
	appendValue(groups, "agjopstAGJOPST");
	appendValue(groups, "ilmvILMV");
	appendValue(groups, "bdenqzBDENQZ");
	appendValue(groups, "hkrHKR");
	appendValue(groups, "fwyFWY");
	appendValue(groups, "cC");

	// Assign to each letter the corresponding group
	String<TWord> lookupTable;
	for(TGroupIter gIt = begin(groups);!atEnd(gIt);++gIt) {
		for(TStringIter sIt = begin(*gIt);!atEnd(sIt);++sIt) {
			if ((TWord) *sIt >= length(lookupTable)) resize(lookupTable, (TWord) *sIt + 1, Generous());
			//std::cout << *sIt << "," << (TWord) *sIt << ":" << position(gIt) << std::endl;
			lookupTable[(TWord) *sIt] = position(gIt);
		}
	}

	// Recode the sequence
	typedef typename Size<TInputString>::Type TSize;
	TSize len = length(input);
	clear(out);
	resize(out, len, Exact());
	for(TSize s =0;s<len;++s) {
		out[s] = lookupTable[(TWord) input[s] ];
	}
	return length(groups);  // Return the new alphabet size
}

template<typename TString, typename TSpec, typename TValue>
void
getScoringMatrix(StringSet<TString, TSpec> const& strSet, Matrix<TValue>& score, unsigned int const ktup) {
	typedef unsigned int TWord;
	typedef StringSet<TString, TSpec> TStringSet;
	typedef typename Size<TStringSet>::Type TStringSetSize;
	typedef String<TWord> TTupelString;
	typedef typename Size<TTupelString>::Type TStringSize;
	typedef String<TTupelString> TTupelStringSet;
	typedef short TCountValue;
	typedef Matrix<TCountValue> THitMatrix;
	typedef typename Size<Matrix<short> >::Type TMatrixSize;

	// Initialization
	TStringSetSize nseq = length(strSet);
	THitMatrix mat;   // Matrix for common k-tupels between sequence i and j
	setDimension(mat, 2);setLength(mat, 0, nseq);setLength(mat, 1, nseq);
	resize(mat);
	setDimension(score, 2);setLength(score, 0, nseq);setLength(score, 1, nseq);
	resize(score);
	for (TMatrixSize row=0;row<nseq;++row) {
		for(TMatrixSize col=0;col<nseq;++col) {
			assignValue(mat, row*nseq+col, 0);
			assignValue(score, row*nseq+col, 0);
		}
	}

	// Transform the StringSet into a set of strings of k-tupels
	TWord alphabet_size;
	TTupelStringSet tupSet;
	resize(tupSet, length(strSet));
	for(TStringSetSize k=0;k<length(strSet);++k) {
		String<TWord> code;
		// Recode the sequence in a reduced alphabet	
		alphabet_size = _recodeSequence(strSet[k], code);
		
		// Break the sequence into k-tupels
		_getTupelString(code, tupSet[k], alphabet_size, ktup);

		//for(TStringSize i = 0;i<length(tupelString);++i) std::cout << tupelString[i] << ",";
		//std::cout << std::endl;
	}

	// Build for each sequence the q-gram Index and count common hits
	String<TCountValue> qIndex;
	String<TCountValue> compareIndex;
	for(TStringSetSize k=0;k<nseq;++k) {
		clear(qIndex);
		fill(qIndex, (unsigned int) pow((double)alphabet_size, (double)ktup), (TCountValue) 0, Exact());
		for(TStringSize i = 0;i<length(tupSet[k]);++i) ++qIndex[ tupSet[k][i] ];
		TWord value;
	    for (TStringSetSize k2=k; k2<nseq; ++k2) {
			clear(compareIndex);
			fill(compareIndex, (unsigned int) pow((double)alphabet_size, (double)ktup), (TCountValue) 0, Exact());
			value = 0;
			for(TStringSize i = 0;i<length(tupSet[k2]);++i) {
				//std::cout << tupSet[k2][i] << "," << compareIndex[ tupSet[k2][i] ] << "," << qIndex[ tupSet[k2][i] ]<< std::endl;
				if (compareIndex[ tupSet[k2][i] ] < qIndex[ tupSet[k2][i] ]) ++value;
				++compareIndex[ tupSet[k2][i] ];
			}
			assignValue(mat, k*nseq+k2, value);
		}
	}

	// Calculate the score
	TValue score0;
	for (TMatrixSize row=0;row<nseq;++row) {
		score0 = getValue(mat, row*nseq+row);
		for(TMatrixSize col=0;col<nseq;++col) {
			assignValue(score, row*nseq+col, (TValue) ((score0 - getValue(mat, min(row,col)*nseq+max(row,col))) / score0 * 3 * 10.0 + 0.5));
		}
	}
	for (TMatrixSize row=0;row<nseq;++row) {
		for(TMatrixSize col=row+1;col<nseq;++col) {
			assignValue(score, row*nseq+col, 100 - min(getValue(score,row*nseq+col), getValue(score, col*nseq+row)));
			assignValue(score, col*nseq+row, getValue(score, row*nseq+col));
		}
	}


	for (TMatrixSize row=0;row<nseq;++row) {
		for(TMatrixSize col=0;col<nseq;++col) {
			std::cout << getValue(score, row*nseq+col) << ",";
		}
		std::cout << std::endl;
	}
}





//////////////////////////////////////////////////////////////////////////////
// T-Coffee Library Handling
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile>
void 
_parse_skipLine(TFile& file)
{
	typedef typename Value<TFile>::Type TValue;
	
	TValue c;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == '\n') break;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
bool
_parse_isDigit(TChar& c)
{
	return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
			(c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
bool
_parse_isLetter(TChar& c)
{
	return ((c == 'a') || (c == 'b') || (c == 'c') || (c == 'd') || (c == 'e') || 
			(c == 'f') || (c == 'g') || (c == 'h') || (c == 'i') || (c == 'j') ||
			(c == 'k') || (c == 'l') || (c == 'm') || (c == 'n') || (c == 'o') || 
			(c == 'p') || (c == 'q') || (c == 'r') || (c == 's') || (c == 't') ||
			(c == 'u') || (c == 'v') || (c == 'w') || (c == 'x') || (c == 'y') || 
			(c == 'z') || (c == 'A') || (c == 'B') || (c == 'C') || (c == 'D') ||
			(c == 'E') || (c == 'F') || (c == 'G') || (c == 'H') || (c == 'I') || 
			(c == 'J') || (c == 'K') || (c == 'L') || (c == 'M') || (c == 'N') ||
			(c == 'O') || (c == 'P') || (c == 'Q') || (c == 'R') || (c == 'S') || 
			(c == 'T') || (c == 'U') || (c == 'V') || (c == 'W') || (c == 'X') ||
			(c == 'Y') || (c == 'Z'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
bool
_parse_isAlphanumericChar(TChar& c)
{
	return ((_parse_isDigit(c)) || (_parse_isLetter(c)) || (c == '_'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile>
int
_parse_readNumber(TFile & file)
{
	typedef typename Value<TFile>::Type TValue;

	// Move to first digit
	TValue c;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (_parse_isDigit(c)) break;
	}

	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isDigit(c)) break;
		append(str, c);
	}
	_streamSeek2G(file, -1);

 	return atoi(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile>
inline String<char>
_parse_readWord(TFile & file)
{
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TString;

	// Move to first letter
	TValue c;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (_parse_isLetter(c)) break;
	}

	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(str, c);
	}
	_streamSeek2G(file, -1);
	return str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile>
inline String<char>
_parse_readIdentifier(TFile & file)
{
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TString;

	// Move to first letter
	TValue c;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (_parse_isAlphanumericChar(c)) break;
	}

	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isAlphanumericChar(c)) break;
		append(str, c);
	}
	_streamSeek2G(file, -1);
	return str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline typename Value<TStringSet>::Type&
_parse_readSequenceData(TFile & file,
						Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Value<TFile>::Type TValue;
    typedef typename Value<TStringSet>::Type TString;

	// Move to first letter
	TValue c;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (_parse_isLetter(c)) break;
	}

	// Read word
	TString* str = new TString(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(*str, c);
	}

	_streamSeek2G(file, -1);
	return *str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
_readLibrary(TFile & file,
			 Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Value<TFile>::Type TValue;
	typedef std::pair<unsigned int, unsigned int> TSeqRes;
	typedef std::pair<TSeqRes, TSeqRes> TKey;
	typedef std::map<TKey, unsigned int> TMap;
	typedef std::map<TSeqRes, unsigned int> TNodeMap;

	TMap my_map;
	unsigned int seq1 = 0;
	unsigned int seq2 = 0;
	unsigned int seq1ToN = 0;
	TValue c;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == '#') {
			seq1 = _parse_readNumber(file);
			seq2 = _parse_readNumber(file);
			_parse_skipLine(file);
		} else if (c == '!') {
			String<char> str = _parse_readIdentifier(file);
			_parse_skipLine(file);
			if (str == "SEQ_1_TO_N") seq1ToN = 1;
		} else if (_parse_isDigit(c)) {
			_streamSeek2G(file, -1);
			unsigned int res1 = _parse_readNumber(file);
			unsigned int res2 = _parse_readNumber(file);
			unsigned int weight = _parse_readNumber(file);
			TKey k = std::make_pair(std::make_pair(seq1, res1 - 1), std::make_pair(seq2,res2 - 1));
			TMap::iterator pos = my_map.find(k);
			if (pos == my_map.end()) {
				my_map.insert(std::make_pair(k, weight));
			} else {
				pos->second += weight;
			}
			_parse_skipLine(file);
		}
	}
	
	// Create the graph

	TNodeMap node_map;
	for(TMap::iterator pos = my_map.begin(); pos!=my_map.end(); ++pos) {
		// Get the key
		TKey k = pos->first;
		//std::cout << k.first.first << "," << k.first.second << " / ";
		//std::cout << k.second.first << "," << k.second.second << ":" <<  pos->second << std::endl;

		// Insert new node_map if necessary
		TNodeMap::iterator nodePos = node_map.find(k.first);
		TId id1;
		if (nodePos == node_map.end()) {
			id1 = addVertex(g, k.first.first - seq1ToN, k.first.second, 1); 
			node_map.insert(std::make_pair(k.first, id1));
		} else {
			id1 = nodePos->second;
		}
		nodePos = node_map.find(k.second);
		TId id2;
		if (nodePos == node_map.end()) {
			id2 = addVertex(g, k.second.first - seq1ToN, k.second.second, 1); 
			node_map.insert(std::make_pair(k.second, id2));
		} else {
			id2 = nodePos->second;
		}

		// Insert a new edge
		addEdge(g,id1,id2,pos->second);
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
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Id<TGraph>::Type TIdType;

	unsigned int nSeq = 0;

	// Ignore first line
	_parse_skipLine(file);
	// Read number of sequences
	nSeq = (unsigned int) _parse_readNumber(file);
	_parse_skipLine(file);
	// Read sequences
	for(unsigned int i=0; i<nSeq; ++i) {
		std::cout << _parse_readIdentifier(file) << ", ";
		std::cout << _parse_readNumber(file) << ", ";
		TIdType id = assignValueById(stringSet(g), _parse_readSequenceData(file,g));
		std::cout << id << ", ";
		std::cout << getValueById(stringSet(g), id) << std::endl;
		SEQAN_ASSERT(id < nSeq)
		_parse_skipLine(file);
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
