#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_ALIGNMENT_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_ALIGNMENT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.NeedlemanWunsch
..summary:Switch to trigger Needleman Wunsch Alignments
..value.NeedlemanWunsch:Needleman Wunsch alignment
*/

struct NeedlemanWunsch_;
typedef Tag<NeedlemanWunsch_> const NeedlemanWunsch;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Gotoh
..summary:Switch to trigger Gotoh Alignments with affine gap costs
..value.Gotoh:Gotoh alignment
*/

struct Gotoh_;
typedef Tag<Gotoh_> const Gotoh;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Gotoh
..summary:Switch to trigger Gotoh Alignments with affine gap costs
..value.Gotoh:Gotoh alignment
*/

struct MyersBitVector_;
typedef Tag<MyersBitVector_> const MyersBitVector;



//////////////////////////////////////////////////////////////////////////////
// Alignment: Simple Traceback Alphabet
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_Byte_2_TraceBack
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_TraceBack<T>::VALUE[256] = 
{
	0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.TraceBack:
..cat:Alphabets
..summary: Trace back values.
..general:Class.SimpleType
..signature:TraceBack
..remarks:
...text:The @Metafunction.ValueSize@ of $TraceBack$ is 3. 
The values are defined in the following way: 0=Diagonal Move, 1=Horizontal Move, 2=Vertical Move
..see:Metafunction.ValueSize
*/
struct _TraceBack {};
typedef SimpleType<unsigned char, _TraceBack> TraceBack;

template <> struct ValueSize< TraceBack > { enum { VALUE = 3 }; };
template <> struct BitsPerValue< TraceBack > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<TraceBack, Byte> { typedef TraceBack Type; };
inline void assign(TraceBack & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_TraceBack<>::VALUE[c_source];
}


//////////////////////////////////////////////////////////////////////////////
// Alignment: Extended Traceback Alphabet (Gotoh)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template <typename T = void>
struct _Translate_Table_Byte_2_TraceBackGotoh
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_TraceBackGotoh<T>::VALUE[256] = 
{
	0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,   11,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.TraceBackGotoh:
..cat:Alphabets
..summary: Trace back values for gotoh.
..general:Class.SimpleType
..signature:TraceBackGotoh
..remarks:
...text:The @Metafunction.ValueSize@ of $TraceBackGotoh$ is 12. 
The values are defined in the following way: Move in Diagonal Matrix, Move in Horizontal Matrix, Move in Vertical Matrix
The values are: 
0=Diag, Diag, Diag; 1=Diag, Diag, Vert; 2=Diag, Hori, Diag; 3=Diag, Hori, Vert; 
4=Hori, Diag, Diag; 5=Hori, Diag, Vert; 6=Hori, Hori, Diag; 7=Hori, Hori, Vert; 
8=Vert, Diag, Diag; 9=Vert, Diag, Vert; 10=Vert, Hori, Diag; 11=Vert, Hori, Vert; 
..see:Metafunction.ValueSize
*/
struct _TraceBackGotoh {};
typedef SimpleType<unsigned char, _TraceBackGotoh> TraceBackGotoh;

template <> struct ValueSize< TraceBackGotoh > { enum { VALUE = 12 }; };
template <> struct BitsPerValue< TraceBackGotoh > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<TraceBackGotoh, Byte> { typedef TraceBackGotoh Type; };
inline void assign(TraceBackGotoh & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_TraceBackGotoh<>::VALUE[c_source];
}



//////////////////////////////////////////////////////////////////////////////
// Alignment: Meyer's Bit Vector algorithm
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TString>
unsigned int
_align_myers_bit_vector(TString const & str1,
						TString const & str2)
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Value<TString>::Type TAlphabet;
	typedef typename Size<TString>::Type TSize;
	typedef String<TWord> BitVector;
	typedef String<BitVector> TLookupTable;
	TSize alphLen = ValueSize<TAlphabet>::VALUE;

	// Preprocessing
	TLookupTable lT;
	resize(lT, alphLen);

	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TSize blockCount;
	if (len2<1) blockCount=1;
	else blockCount=((len2-1) / BitsPerValue<TWord>::VALUE)+1;
	for(TSize i = 0; i < alphLen; ++i) {
		fill(lT[i],blockCount, 0, Exact());
	}
	for(TSize j = 0; j < len2; ++j) {
		TSize pos = convert<TWord>(getValue(str2,j));
		TSize block = j / BitsPerValue<TWord>::VALUE;
		(lT[pos])[block] |= (1<<(j%BitsPerValue<TWord>::VALUE));
	}
	BitVector VP;
	BitVector VN;
	fill(VP, blockCount, ~0, Exact() );
	fill(VN, blockCount, 0, Exact() );
	TWord err = len2;
	
	/*
	// Debug code
	std::cout << "Alphabet size: " << alphLen << ::std::endl;
	std::cout << "Block count: " << blockCount << ::std::endl;
	for(unsigned int i=0;i<alphLen;++i) {
		if ((i<97) || (i>122)) continue;
		std::cout << static_cast<char>(i) << ": ";
		for(int j=0;j<(int)blockCount;++j) {
			for(int bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
			  std::cout << ((lT[i][j] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
	}
	*/

	BitVector X;	
	BitVector D0;
	BitVector HN;	
	BitVector HP;
	BitVector HNcopy;	
	BitVector HPcopy;
	resize(X, blockCount);
	resize(D0, blockCount);
	resize(HN, blockCount);
	resize(HP, blockCount);
	resize(HNcopy, blockCount);
	resize(HPcopy, blockCount);

	for(TSize col = 0; col<len1; ++col) {
		TWord pos = convert<TWord>(getValue(str1,col));	  
		// Addition might produce a carry
		bool carry = 0;
		bool newCarry;
		bool HPcarry = 1;  // We want edit distance!!!
		bool HPnewCarry;
		bool HNcarry = 0;
		bool HNnewCarry;
		for(TWord block=0;block<blockCount;++block) {
			X[block] = lT[pos][block] | VN[block];
			D0[block] = X[block] & VP[block];
			if (( (unsigned int) D0[block] + VP[block] < (unsigned int) D0[block] ) ||
				( (unsigned int) D0[block] + VP[block] < (unsigned int) VP[block] )) newCarry = 1;
			else newCarry = 0;
			D0[block] += VP[block];
			if ((carry) && ( (unsigned int) D0[block] == (unsigned int) ~0)) {
				if (newCarry) {
					std::cerr << "Two carries. Error !!!";
					exit(-1);
				} else {
					newCarry = 1;
				}
			}
			D0[block] += carry;
			carry = newCarry;
			D0[block] ^= VP[block];
			D0[block] |= X[block];
			HN[block] = VP[block] & D0[block];
			HP[block] = VN[block] | ~(VP[block] | D0[block]);
			HPcopy[block] = HP[block];
			HNcopy[block] = HN[block];
			HPnewCarry = ((HPcopy[block] & (1<< (BitsPerValue<TWord>::VALUE - 1)))!=0); 
			HPcopy[block] <<= 1;
			HPcopy[block] += HPcarry;
			HPcarry = HPnewCarry;
			X[block] = HPcopy[block];
			VN[block] = X[block] & D0[block];
			HNnewCarry = ((HNcopy[block] & (1<< (BitsPerValue<TWord>::VALUE - 1)))!=0); 
			HNcopy[block] <<= 1;
			HNcopy[block] += HNcarry;
			HNcarry = HNnewCarry;
			VP[block] = HNcopy[block] | ~(X[block] | D0[block]);
		}
		TSize finalBlock = (len2-1) / BitsPerValue<TWord>::VALUE;
		if ((HP[finalBlock] & (1<<((len2-1)%BitsPerValue<TWord>::VALUE))) != 0) ++err;
		else if ((HN[finalBlock] & (1<<((len2-1)%BitsPerValue<TWord>::VALUE))) != 0) --err;

		/*
		std::cout << err << std::endl;
		for(int block=0;block<(int)blockCount;++block) {
			for(int bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((D0[block] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		for(int block=0;block<(int)blockCount;++block) {
			for(int bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((HN[block] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		for(int block=0;block<(int)blockCount;++block) {
			for(int bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((HP[block] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		for(int block=0;block<(int)blockCount;++block) {
			for(int bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((VN[block] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		for(int block=0;block<(int)blockCount;++block) {
			for(int bit_pos=0;bit_pos<BitsPerValue<TWord>::VALUE;++bit_pos) {
				std::cout << ((VP[block] & (1<<(bit_pos % BitsPerValue<unsigned int>::VALUE))) !=0);
			}
		}
		std::cout << ::std::endl;
		std::cout << ::std::endl;
		*/

	}
	return err;
}


//////////////////////////////////////////////////////////////////////////////
// Alignment: Trace-back, internal functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringSet, typename TId, typename TPos, typename TTraceValue>
inline void
_align_trace_print(TFile& file,
				   TStringSet const& str,
				   TId const id1,
				   TPos const pos1,
				   TId const id2,
				   TPos const pos2,
				   TPos const segLen,
				   TTraceValue const tv)
{
	SEQAN_CHECKPOINT

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	if (tv == (Byte) Horizontal) {
		for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
			_streamPut(file, '(');
			_streamPut(file, (str[0])[i]);
			_streamPut(file, ',');
			_streamPut(file, gapValue<char>());
			_streamPut(file, ')');
			_streamPut(file, '\n');
		}
	}
	else if (tv == (Byte) Vertical) {
		for (int i = pos2 + segLen - 1; i>= (int) pos2;--i) {
			_streamPut(file, '(');
			_streamPut(file, gapValue<char>());
			_streamPut(file, ',');
			_streamPut(file, (str[1])[i]);
			_streamPut(file, ')');
			_streamPut(file, '\n');
		}
	}
	else if (tv == (Byte) Diagonal) {
		int j = pos2 + segLen - 1;
		for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
			_streamPut(file, '(');
			_streamPut(file, (str[0])[i]);
			_streamPut(file, ',');
			_streamPut(file, (str[1])[j]);
			_streamPut(file, ')');
			_streamPut(file, '\n');
			--j;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TStringSet2, typename TId, typename TPos, typename TTraceValue>
inline void
_align_trace_print(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				   TStringSet2 const& str,
				   TId const id1,
				   TPos const pos1,
				   TId const id2,
				   TPos const pos2,
				   TPos const segLen,
				   TTraceValue const tv)
{
	SEQAN_CHECKPOINT

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	if (tv == (Byte) Horizontal) addVertex(g, id1, pos1, segLen);
	else if (tv == (Byte) Vertical) addVertex(g, id2, pos2, segLen);
	else if (tv == (Byte) Diagonal) addEdge(g, addVertex(g, id1, pos1, segLen), addVertex(g, id2, pos2, segLen));
}


//////////////////////////////////////////////////////////////////////////////
// Alignment: Needleman-Wunsch Alignment, constant gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template <typename TAlign, typename TStringSet, typename TTrace>
void
_align_needleman_wunsch_trace(TAlign& align,
							  TStringSet const& str,
							  TTrace const& trace)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Value<TTrace>::Type TTraceValue;


	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};
	
	// Initialization
	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
	TSize numRows = len2;
		
	// Initialize everything
	TTraceValue tv = getValue(trace, (len1-1)*numRows + (len2-1));
	TTraceValue tvOld = tv;  // We need to know when the direction of the trace changes

	TSize segLen = 1;
	if (tv == (Byte) Diagonal) {
		--len1; --len2;
	}
	else if (tv == (Byte) Horizontal) --len1;
	else if (tv == (Byte) Vertical) --len2;

	// Now follow the trace
	do {
		tv = getValue(trace, (len1-1)*numRows + (len2-1));
		if (tv == (Byte) Diagonal) {
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				tvOld = tv; segLen = 1;
			} else ++segLen;
			--len1; --len2;
		} else if (tv == (Byte) Horizontal) {
			//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				tvOld = tv; segLen = 1;
			} else ++segLen;
			--len1;
		} else if (tv == (Byte) Vertical) {
			//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
			if (tv != tvOld) {
				_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
				tvOld = tv; segLen = 1;
			} else ++segLen;
			--len2;
		}
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);

	// Handle the remaining sequence
	if (len1 != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, (Byte) Horizontal);
	else if (len2 != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, (Byte) Vertical);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TStringSet, typename TScoreValue>
TScoreValue
_align_needleman_wunsch(String<TValue, TSpec>& trace,
						TStringSet const& str,
						Score<TScoreValue, Simple> const& sc) 
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TSpec> TTrace;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	// One DP Matrix column
	TColumn column;
		
	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	resize(column, len2 + 1);  
	resize(trace, len1*len2);
	for(TSize row = 0; row <= len2; ++row) assignValue(column, row, row * gap);

	// Classical DP
	TTraceIter it = begin(trace, Standard() );
	TTraceValue tv;
	TScoreValue maxVal = 0;
	TScoreValue tmp = 0;
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagVal = getValue(column, 0);
		assignValue(column, 0, col * gap);
		for(TSize row = 1; row <= len2; ++row) {
			// Get the new maximum	
			maxVal = diagVal + score(const_cast<Score<TScoreValue, Simple>&>(sc), str1[col-1], str2[row-1]);
			tv = (Byte) Diagonal;
			if ((tmp = getValue(column, row) + gap) > maxVal) {
				maxVal = tmp;
				tv = (Byte) Horizontal;
			}
			if ((tmp = getValue(column, (row - 1)) + gap) > maxVal) {
				maxVal = tmp;
				tv = (Byte) Vertical;
			}
			diagVal = getValue(column, row);
			// Assign the new value
			assignValue(column, row, maxVal);
			assignValue(it, tv);
			goNext(it);
		}
	}

	//for(unsigned int i= 0; i<len2;++i) {
	//	for(unsigned int j= 0; j<len1;++j) {
	//		std::cout << (unsigned int) getValue(trace, j*len2 + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return getValue(column, len2);
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Alignment: Gotoh Alignment, affine gap cost
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TStringSet, typename TTrace, typename TVal>
void
_align_gotoh_trace(TAlign& align,		 
				   TStringSet const& str,
				   TTrace const& trace,
				   TVal const initialDir)
{
	SEQAN_CHECKPOINT
	
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Id<TStringSet>::Type TId;

	// TraceBack values for Gotoh
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	TId id1 = positionToId(const_cast<TStringSet&>(str), 0);
	TId id2 = positionToId(const_cast<TStringSet&>(str), 1);	 
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
	TSize numRows = len2;

	// Initialize everything	
	TTraceValue tv = getValue(trace, (len1 - 1)*numRows + (len2 - 1));
	TTraceValue tvOld = initialDir;
	switch( (Byte) initialDir) {
		case Diagonal:
			if ((Byte) tv / (Byte) 4 == 0) tv = (Byte) Diagonal;
			else if ((Byte) tv / (Byte) 4 == 1) tv = (Byte) Horizontal;
			else tv = (Byte) Vertical;
			break;
		case Horizontal:
			if (((Byte) tv / (Byte) 2) % 2 == 0) tv = (Byte) Diagonal;
			else tv = (Byte) Horizontal;
			break;
		case Vertical:
			if ( (Byte) tv % 2 == 0) tv = (Byte) Diagonal;
			else tv = (Byte) Vertical;
			break;
	}
	TSize segLen = 0;

	// Now follow the trace
	do {
		switch( (Byte) tvOld) {
			case Diagonal: 
				//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
				--len1; --len2;
				++segLen;
				if (tv != tvOld) {
					_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
					tvOld = tv; 
					segLen = 0;
				}
				break;
			case Horizontal:
				//std::cout << '(' << ((*str)[0])[len1 - 1] << ',' << '-' << ')' << std::endl;
				--len1;
				++segLen;
				if (tv != tvOld) {
					_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
					tvOld = tv;
					segLen = 0;
				}
				break;
			case Vertical:
				//std::cout << '(' << '-' << ',' << ((*str)[1])[len2-1] << ')' << std::endl;
				--len2;
				++segLen;
				if (tv != tvOld) {
					_align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);
					tvOld = tv; 
					segLen = 0;
				}
				break;
		}
		if ((len1 != 0) && (len2 !=0)) {
			TTraceValue nextTraceValue = getValue(trace, (len1 - 1)*numRows + (len2 - 1));
			switch( (Byte) tv) {
				case Diagonal:
					if ( (Byte) nextTraceValue / (Byte) 4 == 0) tv = (Byte) Diagonal;
					else if ((Byte) nextTraceValue / (Byte) 4 == 1) tv = (Byte) Horizontal;
					else tv = (Byte) Vertical;
					break;
				case Horizontal:
					if (((Byte) nextTraceValue / (Byte) 2) % 2 == 0) tv = (Byte) Diagonal;
					else tv = (Byte) Horizontal;
					break;
				case Vertical:
					if ( (Byte) nextTraceValue % (Byte) 2 == 0) tv = (Byte) Diagonal;
					else tv = (Byte) Vertical;
					break;
			}
		}
	} while ((len1 != 0) && (len2 !=0));
	// Process left-overs
	if (segLen > 0) _align_trace_print(align, str, id1, len1, id2, len2, segLen, tvOld);

	// Handle the remaining sequence
	if (len1 != 0) _align_trace_print(align, str, (TId) id1, (TSize) 0, (TId) 0, (TSize) 0, (TSize) len1, (Byte) Horizontal);
	else if (len2 != 0) _align_trace_print(align, str, (TId) 0, (TSize) 0, (TId) id2, (TSize) 0, (TSize) len2, (Byte) Vertical);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TStringSet, typename TScoreValue>
TScoreValue
_align_gotoh(String<TValue, TSpec>& trace,
			 TStringSet const& str,
			 Score<TScoreValue, Simple> const & sc,
			 TValue& initialDir)
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TSpec> TTrace;
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Iterator<TTrace, Standard>::Type TTraceIter;
	typedef String<TScoreValue> TColumn;
	typedef typename Size<TColumn>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;

	// TraceBack values for Gotoh
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	// The DP Matrix for diagonal walks
	TColumn mat;
	// The DP Matrix for gaps from the left
	TColumn horizontal;
	// The DP Matrix for gaps from the top
	TColumn vertical;

	// Initialization
	TString const& str1 = str[0];
	TString const& str2 = str[1];		
	TSize len1 = length(str1);
	TSize len2 = length(str2);
	TScoreValue gap = scoreGapExtend(sc);
	TScoreValue gapOpen = scoreGapOpen(sc);
	TScoreValue maxVal = 0;
	TScoreValue tmp = 0;
	resize(mat, (len2+1));   // One column for the diagonal matrix
	resize(horizontal, (len2+1));   // One column for the horizontal matrix
	resize(vertical, (len2+1));   // One column for the vertical matrix
	resize(trace, len1*len2);

	TTraceValue tvMat, tvHorizontal, tvVertical;
	
	TScoreValue inf = infimumValue<TScoreValue>() / 2;

	// Classical DP
	TTraceIter it = begin(trace, Standard() );
	assignValue(mat, 0, 0);
	assignValue(horizontal, 0, inf);
	assignValue(vertical, 0, inf);
	for(TSize row = 1; row <= len2; ++row) {
		assignValue(mat, row, gapOpen + (row - 1) * gap);
		assignValue(horizontal, row, inf);
		assignValue(vertical, row, gapOpen + (row - 1) * gap);
	}
	for(TSize col = 1; col <= len1; ++col) {
		TScoreValue diagValMat = getValue(mat, 0);
		TScoreValue diagValHori = getValue(horizontal, 0);
		TScoreValue diagValVert = getValue(vertical, 0);
		TScoreValue diagValVertTmp;
		TScoreValue diagValHoriTmp;
		assignValue(mat, 0, gapOpen + (col - 1) * gap);
		assignValue(horizontal, 0, gapOpen + (col - 1) * gap);
		assignValue(vertical, 0, inf);
		for(TSize row = 1; row <= len2; ++row) {
			// Get the new maximum for vertical
			maxVal = getValue(mat, row - 1) + gapOpen;
			tvVertical = (Byte) Diagonal;
			if ((tmp = getValue(vertical, row - 1) + gap) > maxVal) {
				maxVal = tmp;
				tvVertical = (Byte) Vertical;
			}
			diagValVertTmp = getValue(vertical, row);
			assignValue(vertical, row, maxVal);

			// Get the new maximum for left
			maxVal = getValue(mat, row) + gapOpen;
			tvHorizontal = (Byte) Diagonal;
			if ((tmp = getValue(horizontal, row) + gap) > maxVal) {
				maxVal = tmp;
				tvHorizontal = (Byte) Horizontal;
			}
			diagValHoriTmp = getValue(horizontal, row);
			assignValue(horizontal, row, maxVal);

			// Get the new maximum for mat
			TScoreValue sc_ = score(const_cast<Score<TScoreValue, Simple>&>(sc), str1[col-1], str2[row-1]);
			maxVal = diagValMat + sc_;
			tvMat = (Byte) Diagonal;
			if ((tmp = diagValVert + sc_) > maxVal) {
				maxVal = tmp;
				tvMat = (Byte) Vertical;
			}
			if ((tmp = diagValHori + sc_) > maxVal) {
				maxVal = tmp;
				tvMat = (Byte) Horizontal;
			}

			// Assign the new diagonal values
			diagValMat = getValue(mat, row);
			diagValHori = diagValHoriTmp;
			diagValVert = diagValVertTmp;
			assignValue(mat, row, maxVal);

			// Assign the right trace value
			if (tvMat == (Byte) Diagonal) {
				if (tvHorizontal == (Byte) Diagonal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 0);
					else assignValue(it, 1);
				} else if (tvHorizontal == (Byte) Horizontal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 2);
					else assignValue(it, 3);
				}
			} else if (tvMat == (Byte) Horizontal) {
				if (tvHorizontal == (Byte) Diagonal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 4);
					else assignValue(it, 5);
				} else if (tvHorizontal == (Byte) Horizontal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 6);
					else assignValue(it, 7);
				}
			} else if (tvMat == (Byte) Vertical) {
				if (tvHorizontal == (Byte) Diagonal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 8);
					else assignValue(it, 9);
				} else if (tvHorizontal == (Byte) Horizontal) {
					if (tvVertical == (Byte) Diagonal) assignValue(it, 10);
					else assignValue(it, 11);
				}
			}
			goNext(it);
		}
	}

	maxVal = getValue(mat, len2);
	initialDir = (Byte) Diagonal;
	if ((tmp = getValue(horizontal, len2)) > maxVal) {
		maxVal = tmp;
		initialDir = (Byte) Horizontal;
	}
	if ((tmp = getValue(vertical, len2)) > maxVal) {
		maxVal = tmp;
		initialDir = (Byte) Vertical;
	}

	/*
	for(unsigned int i= 0; i<len2;++i) {
		for(unsigned int j= 0; j<len1;++j) {
			std::cout << (unsigned int) getValue(trace, j*len2 + i) << ',';
		}
		std::cout << std::endl;
	}
	*/

	return maxVal;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 Score<TScoreValue, Simple> const& sc,
				 NeedlemanWunsch)
{
	SEQAN_CHECKPOINT
	// Gap extension score is taken as the constant gap score!!!
	SEQAN_TASSERT(scoreGapOpen(sc) == 0)

	typedef typename Size<TStringSet>::Type TSize;
	  
	TScoreValue maxScore;
	TSize maxLen = length(str[0]);
	TSize tmp;
	if ((tmp = length(str[1])) > maxLen) maxLen = tmp;

	if (maxLen > 10000) {
		String<TraceBack, External<> > trace;
		//open(trace, "/media/sda5/seqan/version7/seqan.dat");

		// Create the trace
		maxScore = _align_needleman_wunsch(trace, str, sc);	

		// Follow the trace and create the graph
		_align_needleman_wunsch_trace(align, str, trace);	
	} else {
		String<TraceBack> trace;

		// Create the trace
		maxScore = _align_needleman_wunsch(trace, str, sc);	

		// Follow the trace and create the graph
		_align_needleman_wunsch_trace(align, str, trace);	
	}
	return maxScore;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue>
TScoreValue
_globalAlignment(TAlign& align,
				 TStringSet const& str,
				 Score<TScoreValue, Simple> const& sc,
				 Gotoh)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	  
	TScoreValue maxScore;
	TSize maxLen = length(str[0]);
	TSize tmp;
	if ((tmp = length(str[1])) > maxLen) maxLen = tmp;

	if (maxLen > 10000) {
		// Trace
		String<TraceBackGotoh, External<> > trace;
		TraceBackGotoh initialDir;

		// Create the trace
		maxScore = _align_gotoh(trace, str, sc, initialDir);	
		// Follow the trace and create the graph
		_align_gotoh_trace(align, str, trace, initialDir);
	} else {
		// Trace
		String<TraceBackGotoh> trace;
		TraceBackGotoh initialDir;

		// Create the trace
		maxScore = _align_gotoh(trace, str, sc, initialDir);	
		// Follow the trace and create the graph
		_align_gotoh_trace(align, str, trace, initialDir);
	}
	return maxScore;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.globalAlignment:
..summary:Computes the best global alignment of the two sequences.
..cat:Alignments
..signature:
globalAlignment(g, score, tag)
globalAlignment(file, str, score, tag)
..param.g:The alignment graph having 2 sequences.
...type:Class.Graph Alignment
..param.str:A string set with 2 sequences.
...type:Class.StringSet
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.tag:A tag indicating the alignment algorithm to use
...remarks:Either NeedlemanWunsch or Gotoh.
..returns:The score value of the best scoring global alignment.
*/
template<typename TFile, typename TStringSet, typename TScoreValue, typename TTag>
TScoreValue
globalAlignment(TFile& file,
				TStringSet const& str,
				Score<TScoreValue, Simple> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(file,str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TTag>
TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, Simple> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _globalAlignment(g,stringSet(g),sc,TTag());
}




//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet>
unsigned int
globalAlignment(TStringSet const& str,
				MyersBitVector)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TStringSet>::Type TSize;
	  
	TSize maxLen = length(str[0]);
	TSize tmp;
	if ((tmp = length(str[1])) > maxLen) maxLen = tmp;

	return _align_myers_bit_vector(str[0], str[1]);	
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
