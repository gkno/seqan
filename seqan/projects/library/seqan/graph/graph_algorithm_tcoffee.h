#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_TCOFFEE_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_TCOFFEE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph: T-Coffee
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////









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

	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		_streamPut(file, '#');
		_streamPutInt(file, sequenceId(g,sV));
		_streamPut(file, ' ');
		_streamPutInt(file, sequenceId(g,tV));
		_streamPut(file, '\n');		
		_streamPutInt(file, segmentBegin(g,sV) + 1);
		_streamPut(file, ' ');
		_streamPutInt(file, segmentBegin(g,tV) + 1);
		_streamPut(file, ' ');
		_streamPutInt(file, getCargo(*it));
		_streamPut(file, '\n');	
	}
	_streamWrite(file, "! SEQ_0_TO_N-1");
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
