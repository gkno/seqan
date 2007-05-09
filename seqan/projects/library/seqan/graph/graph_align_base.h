#ifndef SEQAN_HEADER_GRAPH_ALIGN_BASE_H
#define SEQAN_HEADER_GRAPH_ALIGN_BASE_H

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
.Tag.MyersBitVector
..summary:Switch to trigger a Myers bit vector alignment
..value.MyersBitVector:MyersBitVector alignment
*/

struct MyersBitVector_;
typedef Tag<MyersBitVector_> const MyersBitVector;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Hirschberg
..summary:Switch to trigger a Hirschberg algorithm that uses Gotoh
..value.Hirschberg:Hirschberg alignment
*/

struct Hirschberg_;
typedef Tag<Hirschberg_> const Hirschberg;

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




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
