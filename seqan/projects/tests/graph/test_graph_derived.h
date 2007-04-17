#ifndef SEQAN_HEADER_TEST_GRAPH_DERIVED_H
#define SEQAN_HEADER_TEST_GRAPH_DERIVED_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

void Test_Oracle() {
	Graph<Automaton<char> > g;
	createOracleOnReverse(g,"announce");
	SEQAN_TASSERT(parseString(g, 0, "e") == 1)
	SEQAN_TASSERT(parseString(g, 0, "ec") == 2)
	SEQAN_TASSERT(parseString(g, 0, "n") == 3)
	SEQAN_TASSERT(parseString(g, 0, "a") == 8)
	SEQAN_TASSERT(parseString(g, 0, "nn") == 7)

	Graph<Automaton<Dna> > g2;
	createOracle(g2,"ATATA");
	SEQAN_TASSERT(parseString(g2, 0, "A") == 1)
	SEQAN_TASSERT(parseString(g2, 0, "T") == 2)
	SEQAN_TASSERT(parseString(g2, 0, "AT") == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Trie() {
	Graph<Automaton<char> > g;
	String<String<unsigned int> > pos;
	String<String<char> > keywords;
	appendValue(keywords, String<char>("announce"));
	appendValue(keywords, String<char>("annual"));
	appendValue(keywords, String<char>("annually"));
	createTrie(g,pos,keywords);

	// Output
	// File output
	fstream strm;
	strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
	String<String<char> > nodeMap;
	_createTrieNodeAttributes(g, pos, nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(strm,g,nodeMap,edgeMap,DotDrawing());
	strm.close();

	Graph<Automaton<Dna> > gDna;
	clear(pos);
	String<String<Dna> > keyw;
	appendValue(keyw, String<Dna>("ATATATA"));
	appendValue(keyw, String<Dna>("TATAT"));
	appendValue(keyw, String<Dna>("ACGATAT"));
	createTrie(gDna,pos,keyw);

	// Output
	// File output
	fstream strm2;
	strm2.open(TEST_PATH "my_trie_dna.dot", ios_base::out | ios_base::trunc);
	clear(nodeMap);
	_createTrieNodeAttributes(gDna, pos, nodeMap);
	clear(edgeMap);
	_createEdgeAttributes(gDna,edgeMap);
	write(strm2,gDna,nodeMap,edgeMap,DotDrawing());
	strm2.close();
}

}

#endif

