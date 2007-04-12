#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_TCOFFEE_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_TCOFFEE_H

namespace SEQAN_NAMESPACE_MAIN
{

	
//////////////////////////////////////////////////////////////////////////////
// Graph: T-Coffee - New alphabet for amino acid groups
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template <typename T = void>
struct _Translate_Table_AA_2_AAGroups
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroups<T>::VALUE[24] = 
{
	0, // 0 Ala Alanine                 
	3, // 1 Arg Arginine                
	2, // 2 Asn Asparagine              
	2, // 3 Asp Aspartic Acid           
	5, // 4 Cys Cystine                 
	2, // 5 Gln Glutamine               
	2, // 6 Glu Glutamic Acid           
	0, // 7 Gly Glycine                 
	3, // 8 His Histidine               
	1, // 9 Ile Isoleucine              
	1, //10 Leu Leucine                 
	3, //11 Lys Lysine                  
	1, //12 Met Methionine              
	4, //13 Phe Phenylalanine           
	0, //14 Pro Proline                 
	0, //15 Ser Serine                  
	0, //16 Thr Threonine               
	4, //17 Trp Tryptophan              
	4, //18 Tyr Tyrosine                
	1, //19 Val Valine                  
	2, //20 Aspartic Acid, Asparagine   
	2, //21 Glutamic Acid, Glutamine    
	6, //22 Unknown                     
	6  //23 Terminator                  
};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroups:
..cat:Alphabets
..summary:Alphabet for Amino Acid Groups.
..general:Class.SimpleType
..signature:AAGroups
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroups$ is 7. 
The groups are defined in the following way: "agjopst"=0, "ilmv"=1, "bdenqz"=2, "hkr"=3, "fwy"=4, "c"=5, all others = 6
...text:Objects of type $AAGroups$ cannot be converted into other types.
...text:$AAGroups$ is typedef for $SimpleType<char,_AAGroups>$, while $_AAGroups$ is an auxilliary specialization tag class.
..see:Metafunction.ValueSize
*/
struct _AAGroups {};
typedef SimpleType<unsigned char,_AAGroups> AAGroups;

template <> struct ValueSize< AAGroups > { enum { VALUE = 7 }; };
template <> struct BitsPerValue< AAGroups > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroups assignment

template <>
struct CompareType<AAGroups, AminoAcid> { typedef AAGroups Type; };
inline void assign(AAGroups & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroups<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroups, Byte> { typedef AAGroups Type; };
inline void assign(AAGroups & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroups<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroups, Ascii> { typedef AAGroups Type; };
inline void assign(AAGroups & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroups<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroups, Unicode> { typedef AAGroups Type; };
inline void assign(AAGroups & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroups<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////




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
				for(TWord k2 = k; k2 < nseq; ++k2)
					assignValue(mat, k*nseq+k2, getValue(mat, k*nseq+k2) + min(counter[k], counter[k2]));
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
	typedef StringSet<String<AAGroups>, ConcatVirtual<> > TStringSetAA;
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
				for(TWord k2 = k; k2 < nseq; ++k2)
					assignValue(mat, k*nseq+k2, getValue(mat, k*nseq+k2) + min(counter[k], counter[k2]));
		}
	}

	// Transform the hit matrix into a score matrix
	_hitToScoreMatrix(mat, score, nseq);
}


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
				assignValue(dist, row*nseq+col, maxDist - (TValue) getValue(score, row*nseq+col));
				assignValue(dist, col*nseq+row, maxDist - (TValue) getValue(score, row*nseq+col));
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
