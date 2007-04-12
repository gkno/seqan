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

template<typename TFile, typename TChar>
void 
_parse_skipLine(TFile& file, TChar& c)
{
	SEQAN_CHECKPOINT
		
	if (c == '\n') return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == '\n') break;
	}
	c = _streamGet(file);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
void 
_parse_skipWhitespace(TFile& file, TChar& c)
{
	SEQAN_CHECKPOINT
	
	if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) break;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
bool
_parse_isDigit(TChar const c)
{
	SEQAN_CHECKPOINT
	//return (((unsigned) c >=  48) && ((unsigned) c <=  57));
	return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
		    (c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
bool
_parse_isLetter(TChar const c)
{
	SEQAN_CHECKPOINT
	//return ((((unsigned) c >=  97) && ((unsigned) c <=  122)) || (((unsigned) c >=  65) && ((unsigned) c <=  90)));
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
_parse_isAlphanumericChar(TChar const c)
{
	SEQAN_CHECKPOINT
	return ((_parse_isDigit(c)) || (_parse_isLetter(c)) || (c == '_'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
int
_parse_readNumber(TFile & file, TChar& c)
{
	SEQAN_CHECKPOINT
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isDigit(c)) break;
		append(str, c);
	}
 	return atoi(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline String<char>
_parse_readIdentifier(TFile & file, TChar& c)
{
	SEQAN_CHECKPOINT
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isAlphanumericChar(c)) break;
		append(str, c);
	}
	return str;
}

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
