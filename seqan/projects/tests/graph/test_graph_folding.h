#ifndef SEQAN_HEADER_TEST_GRAPH_FOLDING_H
#define SEQAN_HEADER_TEST_GRAPH_FOLDING_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrixAlign, typename TMatrix> 
inline void
_mutualInformationContent(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						  TMatrixAlign const& align,
						  TMatrix& mat)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Value<TMatrixAlign>::Type TChar;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;
				
	TSize nseq = length(stringSet(g));
	TSize colLen = length(align) / nseq;
	TSize value_size = ValueSize<TAlphabet>::VALUE;
	TChar gapChar = gapValue<TChar>();
	

	resize(mat, colLen * colLen);
	for(TSize i = 0; i < colLen; ++i) {
		for(TSize j = i + 1; j < colLen; ++j) {
			String<TValue> f_i;
			fill(f_i, value_size, 0);
			String<TValue> f_j;
			fill(f_j, value_size, 0);
			String<TValue> f_ij;
			fill(f_ij, value_size * value_size, 0);
			TSize count = 0;
			for(TSize seq = 0; seq < nseq; ++seq) {
				if ((getValue(align, seq * colLen + i) != gapChar) && (getValue(align, seq * colLen + j) != gapChar)) {
					TAlphabet c1 = (TAlphabet) getValue(align, seq * colLen + i);
					TAlphabet c2 = (TAlphabet) getValue(align, seq * colLen + j);
					value(f_i, (Byte) c1) += 1;
					value(f_j, (Byte) c2) += 1;
					value(f_ij, ((Byte) c1) * value_size + (Byte) c2) += 1;
					++count;
				}
			}
			// Calculate h_ij
			TValue h_ij = 0;
			if (count != 0) {
				for(TSize x = 0; x < value_size;++x) {
					value(f_i, x) /= count;
					//std::cout << getValue(f_i, x) << std::endl;
					for(TSize y = 0; y < value_size;++y) {
						if (x == 0) {
							value(f_j, y) /= count;
							//std::cout << getValue(f_j, y) << std::endl;
						}
						value(f_ij, x * value_size + y) /= count;
						//std::cout << getValue(f_ij, x * value_size + y) << std::endl;
						if (getValue(f_ij, x * value_size + y) != 0) {
							TValue ratio = getValue(f_ij, x * value_size + y) / (getValue(f_i, x) * getValue(f_j, y));
							h_ij += getValue(f_ij, x * value_size + y) * (log((TValue) ratio) / log((TValue) 2));						}
					}
				}
			}
			assignValue(mat, i * colLen + j, h_ij);
			assignValue(mat, j * colLen + i, h_ij);
			//std::cout << i << ',' << j << ':' << h_ij << std::endl;
		}
	}

	//// Debug code
	//for(TSize i = 0; i < colLen; ++i) {
	//	for(TSize j = i+1; j < colLen; ++j) {
	//		std::cout << i << ',' << j << ':' << getValue(mat, i * colLen + j) << std::endl;
	//	}
	//}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TMatrix> 
inline void
mutualInformationContent(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						 TMatrix& mat)
{
	SEQAN_CHECKPOINT	
	String<char> align;
	if (convertAlignment(g, align)) {
		_mutualInformationContent(g, align, mat);
	}
}


//////////////////////////////////////////////////////////////////////////////

void Test_Nussinov() {
	typedef String<char> TString;
	typedef Size<TString>::Type TSize;
	typedef Value<TString>::Type TCharacter;
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	
	typedef std::map<std::pair<TCharacter, TCharacter>, unsigned int> TBasePairMap;
	TBasePairMap pairMap;
	pairMap.insert(std::make_pair(std::make_pair('a','u'),1));
	pairMap.insert(std::make_pair(std::make_pair('u','a'),1));
	pairMap.insert(std::make_pair(std::make_pair('c','g'),1));
	pairMap.insert(std::make_pair(std::make_pair('g','c'),1));
	Score<TBasePairMap, ScoreNussinov> sc = Score<TBasePairMap, ScoreNussinov>(pairMap);
	TString str("gggaaaucc");
	TGraph g;
	unsigned int score = rnaFolding(g, str, sc, Nussinov() );
	SEQAN_TASSERT(numEdges(g) == 11)
	SEQAN_TASSERT(numVertices(g) == length(str))
	SEQAN_TASSERT(score == 3)

	pairMap.clear();
	pairMap.insert(std::make_pair(std::make_pair('g','u'),1));
	pairMap.insert(std::make_pair(std::make_pair('u','g'),1));
	pairMap.insert(std::make_pair(std::make_pair('a','u'),2));
	pairMap.insert(std::make_pair(std::make_pair('u','a'),2));
	pairMap.insert(std::make_pair(std::make_pair('c','g'),3));
	pairMap.insert(std::make_pair(std::make_pair('g','c'),3));
	sc = Score<TBasePairMap, ScoreNussinov>(pairMap);
	str = "gcagcacccaaagggaauaugggauacgcgua";
	clear(g);
	score = rnaFolding(g, str, sc, 3, Nussinov() );
	SEQAN_TASSERT(numEdges(g) == 41)
	SEQAN_TASSERT(numVertices(g) == length(str))
	SEQAN_TASSERT(score == 25)
	
	//String<char> names;
	//resize(names, length(str));
	//for(TSize i=0;i<length(str);++i) {
	//	assignValue(names,i,str[i]);
	//}
	//String<String<char> > nodeMap;
	//_createNodeAttributes(g,nodeMap,names);
	//String<String<char> > edgeMap;
	//_createEdgeAttributes(g,edgeMap);
	//fstream strm;
	//strm.open(TEST_PATH "my_rna_graph.dot", ios_base::out | ios_base::trunc);
	//write(strm,g,nodeMap,edgeMap,DotDrawing());
	//strm.close();
}

//////////////////////////////////////////////////////////////////////////////


void Test_MutualInformation() {
//____________________________________________________________________________
// Rna stuff

	// Mutual information content
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef	Id<TStringSet>::Type TId;
	typedef	Size<TStringSet>::Type TSize;

	// Create an alignment
	TStringSet str;
	TString str0("cgcgataa");
	assignValueById(str, str0);
	TString str1("cggccgcc");
	assignValueById(str, str1);
	TString str2("cgcggcgg");
	assignValueById(str, str2);
	TString str3("cggctatt");
	assignValueById(str, str3);

	TGraph g(str);
	TSize pos = 0;
	while(pos < length(getValue(str,0))) {
		TSize seq = 1;
		addVertex(g, 0, pos, 1);
		while (seq < length(str)) {
			TVertexDescriptor v2 = addVertex(g, seq, pos, 1);
			for(TSize i = seq; i>0;--i) {
				addEdge(g, (TVertexDescriptor) (v2 - i), (TVertexDescriptor) v2);
			}
			++seq;
		}
		++pos;
	}

	String<double> mat;
	mutualInformationContent(g, mat);
	TSize matrix_size = (TSize) sqrt((double) length(mat));

	SEQAN_TASSERT(getValue(mat, 0 * matrix_size + 1) == 0)
	SEQAN_TASSERT(getValue(mat, 2 * matrix_size + 3) == 1)
	SEQAN_TASSERT(getValue(mat, 4 * matrix_size + 5) == 2)
	SEQAN_TASSERT(getValue(mat, 6 * matrix_size + 7) == 2)
}


//////////////////////////////////////////////////////////////////////////////

void Test_GraphFolding() {
	// Nussinov
	Test_Nussinov();
	Test_MutualInformation();// Use of alignment graph for Rna stuff

	debug::verifyCheckpoints("projects/library/seqan/graph/graph_fold_nussinov.h");
}

}

#endif

