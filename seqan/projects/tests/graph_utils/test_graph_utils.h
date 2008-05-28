#ifndef SEQAN_HEADER_TEST_GRAPH_UTILS_H
#define SEQAN_HEADER_TEST_GRAPH_UTILS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TGraph>
void  Test_GraphDrawing_Tmp() {
	// Create a dummy graph
	TGraph g;
	addVertex(g);addVertex(g);addEdge(g,0,1);
	addVertex(g);addEdge(g,0,2);addVertex(g); 
	addVertex(g);addVertex(g);addEdge(g,3,4);

	// Dot Drawing
	write(std::cout,g,DotDrawing());
}

//////////////////////////////////////////////////////////////////////////////

void  Test_GraphDrawing() {
	Test_GraphDrawing_Tmp<Graph<Directed<> > >();
	Test_GraphDrawing_Tmp<Graph<Undirected<> > >();
	Test_GraphDrawing_Tmp<Graph<Tree<> > >();

	// Automat
	Graph<Automaton<Dna> > automat;
	createRoot(automat);addVertex(automat);addEdge(automat,0,1,'a');
	addVertex(automat);addEdge(automat,0,2,'g');
	// Dot Drawing
	write(std::cout,automat,DotDrawing());

	// Trie
	Graph<Automaton<char> > trie;
	String<String<unsigned int> > pos;
	String<String<char> > keywords;
	appendValue(keywords, String<char>("announce"));
	appendValue(keywords, String<char>("annual"));
	appendValue(keywords, String<char>("annually"));
	createTrie(trie,pos,keywords);
	// Dot Drawing
	String<String<char> > nodeMap;
	_createTrieNodeAttributes(trie, pos, nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(trie,edgeMap);
	write(std::cout,trie,nodeMap, edgeMap, DotDrawing());

	// WordGraph
	typedef Graph<Automaton<Dna, String<Dna>, WordGraph<> > > TWordGraph;
	TWordGraph wordGr;
	addVertex(wordGr);addVertex(wordGr);addVertex(wordGr);addVertex(wordGr);
	assignRoot(wordGr,3);root(wordGr) = 2;addEdge(wordGr,0,3,"ag");
	addVertex(wordGr);addVertex(wordGr);addEdge(wordGr,0,5,"g");
	addVertex(wordGr);addVertex(wordGr);addEdge(wordGr,3,1,"aggg");
	addEdge(wordGr,3,4,"gg");addEdge(wordGr,5,2,"aggg");addEdge(wordGr,5,7,"g");
	addEdge(wordGr,7,6,"g");assignRoot(wordGr,0);
	write(std::cout,wordGr,DotDrawing());

	// Alignment Graph
	typedef String<char> TAlignString;
	typedef StringSet<TAlignString, Dependent<> > TAlignStringSet;
	typedef Graph<Alignment<TAlignStringSet, void> > TAlignmentGraph;
	typedef Id<TAlignmentGraph>::Type TId;
	typedef VertexDescriptor<TAlignmentGraph>::Type TVD;
	typedef EdgeDescriptor<TAlignmentGraph>::Type TED;	
	TAlignStringSet al;
	TAlignString al0("Garfieldthelastfatcat");
	TId i0 = assignValueById(al, al0);
	TAlignString al1("Garfieldthefastcat");
	TId i1 = assignValueById(al, al1);
	TAlignString al2("Garfieldtheveryfastcat");
	TId i2 = assignValueById(al, al2);
	TAlignString al3("thefatcat");
	TId i3 = assignValueById(al, al3);
	TAlignmentGraph gAl(al);
	TVD vH = addVertex(gAl, i1, 8, 3);
	TVD vT = addVertex(gAl, i1, 13, 1);TVD vS = addVertex(gAl, i3, 6, 3);
	TVD vW = addVertex(gAl, i2, 18, 1);TVD vA = addVertex(gAl, i0, 0, 8);
	TVD vM = addVertex(gAl, i2, 11, 4);
	TVD vK = addVertex(gAl, i2, 0, 8);TVD vC = addVertex(gAl, i0, 11, 4);TVD vD = addVertex(gAl, i0, 15, 2);
	TVD vF = addVertex(gAl, i0, 18, 3);TVD vG = addVertex(gAl, i1, 0, 8);addEdge(gAl, vA, vG);
	TVD vI = addVertex(gAl, i1, 11, 2);TVD vQ = addVertex(gAl, i3, 3, 2);TVD vB = addVertex(gAl, i0, 8, 3);
	TVD vU = addVertex(gAl, i1, 14, 1);TVD vE = addVertex(gAl, i0, 17, 1);TVD vJ = addVertex(gAl, i1, 15, 3);
	TVD vL = addVertex(gAl, i2, 8, 3);
	addEdge(gAl, vH, vL);
	TVD vN = addVertex(gAl, i2, 15, 2);TVD vV = addVertex(gAl, i2, 17, 1);
	TVD vO = addVertex(gAl, i2, 19, 3);TVD vP = addVertex(gAl, i3, 0, 3);TVD vR = addVertex(gAl, i3, 5, 1);
	addEdge(gAl, vA, vK);addEdge(gAl, vG, vK);addEdge(gAl, vB, vH);addEdge(gAl, vB, vL);
	addEdge(gAl, vB, vP);addEdge(gAl, vH, vP);addEdge(gAl, vL, vP);addEdge(gAl, vC, vM);
	addEdge(gAl, vD, vI);addEdge(gAl, vD, vQ);addEdge(gAl, vD, vN);addEdge(gAl, vI, vQ);
	addEdge(gAl, vI, vN);addEdge(gAl, vQ, vN);addEdge(gAl, vT, vV);addEdge(gAl, vE, vU);
	addEdge(gAl, vE, vW);addEdge(gAl, vE, vR);addEdge(gAl, vU, vW);addEdge(gAl, vU, vR);
	addEdge(gAl, vW, vR);addEdge(gAl, vF, vJ);addEdge(gAl, vF, vO);addEdge(gAl, vF, vS);
	addEdge(gAl, vJ, vO);addEdge(gAl, vJ, vS);addEdge(gAl, vO, vS);
	write(std::cout,gAl, DotDrawing());
}

	

//////////////////////////////////////////////////////////////////////////////

void  Test_Alphabets() {
	// Test Rna
	Rna5 r = 'U';
	SEQAN_TASSERT(r == 3)
	r = Byte(0);
	SEQAN_TASSERT(r == 'a')
	char c = 'c';
	r = c;
	c = r;
	SEQAN_TASSERT(r == Rna5('c'))
	SEQAN_TASSERT((c == 'C' || c == 'c'))
	Unicode u = 'g';
	r = u;
	SEQAN_TASSERT(r == Rna5('g'))
	

	// Test Dayhoff
	AAGroupsDayhoff gr;
	gr = AminoAcid('T');
	SEQAN_TASSERT(gr == 0)
	gr = Byte(3);
	SEQAN_TASSERT(gr == 2)
	gr = char('c');
	SEQAN_TASSERT(gr == 5)
	gr = Unicode('j');
	SEQAN_TASSERT(gr == 6)

	// Test SeB6
	AAGroupsSeB6 gr1;
	gr1 = AminoAcid('T');
	SEQAN_TASSERT(gr1 == 0)
	gr1 = Byte(3);
	SEQAN_TASSERT(gr1 == 2)
	gr1 = char('c');
	SEQAN_TASSERT(gr1 == 1)
	gr1 = Unicode('j');
	SEQAN_TASSERT(gr1 == 6)

	// Test SeB8
	AAGroupsSeB8 gr2;
	gr2 = AminoAcid('T');
	SEQAN_TASSERT(gr2 == 0)
	gr2 = Byte(3);
	SEQAN_TASSERT(gr2 == 2)
	gr2 = char('c');
	SEQAN_TASSERT(gr2 == 1)
	gr2 = Unicode('j');
	SEQAN_TASSERT(gr2 == 8)

	// Test Murphy
	AAGroupsMurphy gr3;
	gr3 = AminoAcid('T');
	SEQAN_TASSERT(gr3 == 9)
	gr3 = Byte(3);
	SEQAN_TASSERT(gr3 == 2)
	gr3 = char('c');
	SEQAN_TASSERT(gr3 == 1)
	gr3 = Unicode('j');
	SEQAN_TASSERT(gr3 == 10)

	// Test SolisG10
	AAGroupsSolisG10 gr4;
	gr4 = AminoAcid('T');
	SEQAN_TASSERT(gr4 == 8)
	gr4 = Byte(3);
	SEQAN_TASSERT(gr4 == 2)
	gr4 = char('c');
	SEQAN_TASSERT(gr4 == 1)
	gr4 = Unicode('j');
	SEQAN_TASSERT(gr4 == 10)

	// Test SolisD10
	AAGroupsSolisD10 gr5;
	gr5 = AminoAcid('T');
	SEQAN_TASSERT(gr5 == 6)
	gr5 = Byte(3);
	SEQAN_TASSERT(gr5 == 2)
	gr5 = char('c');
	SEQAN_TASSERT(gr5 == 1)
	gr5 = Unicode('j');
	SEQAN_TASSERT(gr5 == 10)

	// Test LiB10
	AAGroupsLiB10 gr6;
	gr6 = AminoAcid('T');
	SEQAN_TASSERT(gr6 == 0)
	gr6 = Byte(3);
	SEQAN_TASSERT(gr6 == 2)
	gr6 = char('c');
	SEQAN_TASSERT(gr6 == 1)
	gr6 = Unicode('j');
	SEQAN_TASSERT(gr6 == 10)

	// Test LiA10
	AAGroupsLiA10 gr7;
	gr7 = AminoAcid('T');
	SEQAN_TASSERT(gr7 == 9)
	gr7 = Byte(3);
	SEQAN_TASSERT(gr7 == 1)
	gr7 = char('c');
	SEQAN_TASSERT(gr7 == 0)
	gr7 = Unicode('j');
	SEQAN_TASSERT(gr7 == 10)

	// Test SeV10
	AAGroupsSeV10 gr8;
	gr8 = AminoAcid('T');
	SEQAN_TASSERT(gr8 == 0)
	gr8 = Byte(3);
	SEQAN_TASSERT(gr8 == 2)
	gr8 = char('c');
	SEQAN_TASSERT(gr8 == 1)
	gr8 = Unicode('j');
	SEQAN_TASSERT(gr8 == 10)

	// Test SeB10
	AAGroupsSeB10 gr9;
	gr9 = AminoAcid('T');
	SEQAN_TASSERT(gr9 == 0)
	gr9 = Byte(3);
	SEQAN_TASSERT(gr9 == 2)
	gr9 = char('c');
	SEQAN_TASSERT(gr9 == 1)
	gr9 = Unicode('j');
	SEQAN_TASSERT(gr9 == 10)

	// Test SeB14
	AAGroupsSeB14 gr10;
	gr10 = AminoAcid('T');
	SEQAN_TASSERT(gr10 == 12)
	gr10 = Byte(3);
	SEQAN_TASSERT(gr10 == 2)
	gr10 = char('c');
	SEQAN_TASSERT(gr10 == 1)
	gr10 = Unicode('j');
	SEQAN_TASSERT(gr10 == 14)
}


//////////////////////////////////////////////////////////////////////////////


void Test_GraphUtils() {
	Test_GraphDrawing();
	Test_Alphabets();

	debug::verifyCheckpoints("projects/library/seqan/graph_utils/graph_drawing.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_utils/graph_utility_alphabets.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_utils/graph_utility_parsing.h");
}


}

#endif

