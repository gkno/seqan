#ifndef SEQAN_HEADER_TEST_GRAPH_FOLDING_H
#define SEQAN_HEADER_TEST_GRAPH_FOLDING_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

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

void Test_GraphFolding() {
	// Nussinov
	Test_Nussinov();

	debug::verifyCheckpoints("projects/library/seqan/graph/graph_fold_nussinov.h");
}

}

#endif

