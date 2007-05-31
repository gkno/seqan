#ifndef SEQAN_HEADER_TEST_GRAPH_ALIGNMENT_H
#define SEQAN_HEADER_TEST_GRAPH_ALIGNMENT_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

void  Test_NeedlemanWunsch() {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("annual");	assignValueById(str, str0);
	TString str1("annealing"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(0,-1,-1,0);
	int score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "annual")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "anneal")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "ing")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//int score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "annealing";
	str[1] = "annual";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "annual")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "anneal")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 7)) == "ing")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ThisisGarfieldthecat";
	str[1] = "Garfield";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "Thisis")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 9)) == "Garfield")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "Garfield")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 17)) == "thecat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "Garfield";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "Thisis")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 9)) == "Garfield")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "Garfield")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 17)) == "thecat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cat";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ThisisGarfieldthe")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 18)) == "cat")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "cat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ThisisGarfieldthecat";
	str[1] = "cat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ThisisGarfieldthe")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 18)) == "cat")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "cat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
	//SEQAN_TASSERT(score == score2)
}

//////////////////////////////////////////////////////////////////////////////

void  Test_Gotoh() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("ttagt");	assignValueById(str, str0);
	TString str1("ttgt"); assignValueById(str, str1);
	TGraph g(str);
	Score<double> score_type = Score<double>(1,-1,-1,-2);
	double score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	//std::cout << g << std::endl;
	//double score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttgt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "tagt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttagt";
	str[1] = "tagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "cttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cttccagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "cc")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "ag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 7)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "ag")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 7)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "cttccagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "cc")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "ag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "ag")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 7)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttattaa";
	str[1] = "aaa";
	assignStringSet(g, str);
	Score<double> score_type2 = Score<double>(10,-1,-1,-2);
	score = globalAlignment(g, score_type2, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type2, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "aaa";
	str[1] = "ttattaa";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type2, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type2, Gotoh() );
	//SEQAN_TASSERT(score == score2)

	// Note: Depending on the used recursion formula, the gotoh algorithms can differ !!!
	// Gotoh: Vertical gap and subsequent horizontal gap allowed
	// Gotoh3: Vertical gap and subsequent horizontal gap is not allowed
	//Score<double> score_scheme = Score<double>(5,-4,-0.5,-2);
	//str[0] = "tttggttt";
	//str[1] = "tttccttt";
	//assignStringSet(g, str);
	//score = globalAlignment(std::cout, str, score_scheme, Gotoh());
	//std::cout << std::endl;
	//score = globalAlignment(std::cout, str, score_scheme, Gotoh3());
}

//////////////////////////////////////////////////////////////////////////////

void  Test_MyersBitVector() {
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("annealing");	assignValueById(str, str0);
	TString str1("annual"); assignValueById(str, str1);
	int score1 = globalAlignment(str, MyersBitVector() );
	Score<int> score_type = Score<int>(0,-1,-1,0);
	int score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	Score<int> score_type2 = Score<int>(0,-1,-1,-1);
	int score3 = globalAlignment(str, score_type2, Gotoh() );
	int score4 = globalAlignment(str, score_type2, Hirschberg() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)
	SEQAN_TASSERT(score3 == score4)

	str[0] = "annual";
	str[1] = "annealing";
	score1 = globalAlignment(str, MyersBitVector() );
	score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	score3 = globalAlignment(str, score_type2, Gotoh() );
	score4 = globalAlignment(str, score_type2, Hirschberg() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)
	SEQAN_TASSERT(score3 == score4)

	str[0] = "cttagt";
	str[1] = "ttag";
	score1 = globalAlignment(str, MyersBitVector() );
	score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	score3 = globalAlignment(str, score_type2, Gotoh() );
	score4 = globalAlignment(str, score_type2, Hirschberg() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)
	SEQAN_TASSERT(score3 == score4)

	str[0] = "ttag";
	str[1] = "cttccagt";
	score1 = globalAlignment(str, MyersBitVector() );
	score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	score3 = globalAlignment(str, score_type2, Gotoh() );
	score4 = globalAlignment(str, score_type2, Hirschberg() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)
	SEQAN_TASSERT(score3 == score4)
	 
	str[0] = "ttagttagttagttagttagttagttagttagttagttagttagttagttagttagttag";
	str[1] = "cttccagtcttccagtcttccagtcttccagtcttccagtcttccagtcttccagtcttccagt";
	score1 = globalAlignment(str, MyersBitVector() );
	score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	score3 = globalAlignment(str, score_type2, Gotoh() );
	score4 = globalAlignment(str, score_type2, Hirschberg() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)
	SEQAN_TASSERT(score3 == score4)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Hirschberg() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("ttagt");	assignValueById(str, str0);
	TString str1("ttgt"); assignValueById(str, str1);
	TGraph g(str);
	Score<double> score_type = Score<double>(1,-1,-1,-2);
	double score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	//std::cout << g << std::endl;
	//double score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttgt";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "tagt";
	str[1] = "attagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "at")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "attagt";
	str[1] = "tagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "at")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "cttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 4)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "cttccagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "cc")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "ag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 7)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "ag")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 7)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttag";
	str[1] = "cttccagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "c")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "cc")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "ag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "t")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "ag")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 7)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "ttattaa";
	str[1] = "aaa";
	assignStringSet(g, str);
	Score<double> score_type2 = Score<double>(10,-1,-1,-2);
	score = globalAlignment(g, score_type2, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type2, Hirschberg() );
	//SEQAN_TASSERT(score == score2)

	str[0] = "aaa";
	str[1] = "ttattaa";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type2, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type2, Hirschberg() );
	//SEQAN_TASSERT(score == score2)
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

void  Test_Runtime() {
	typedef String<Dna5, External<> > TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Id<TStringSet>::Type TId;

	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\matches\\");
	String<char> out_path("Z:\\matches\\out\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/matches/");
	String<char> out_path("/home/takifugu/rausch/matches/out/");
#endif

	Score<double> score_type = Score<double>(5,-4,-0.5,-10);
	TStringSet str;
	clock_t startTime;
	clock_t duration;
	TString str0;
	TString str1;

	//unsigned int chrom = 24;
	std::stringstream s1;
	//s1 << out_path << "H.chr." << chrom - 1;
	s1 << out_path << "a.fasta";
	bool f = open(str0, s1.str().c_str());
	if (!f) {
		exit(-1);
	}
	std::stringstream s2;
	//s2 << out_path << "W.chr." << chrom - 1;
	s2 << out_path << "b.fasta";
	f = open(str1, s2.str().c_str());
	if (!f) {
		exit(-1);
	}
	
	//ifstream file;
	//std::stringstream input;
	//input << in_path << "a.fasta";
	//file.open(input.str().c_str(), ios_base::in | ios_base::binary);
	//std::stringstream output;
	//output << out_path << "a.fasta";
	//open(str0, output.str().c_str());
	//read(file, str0, Fasta());
	//close(str0);
	//file.close();
	//exit(0);

		
	assignValueById(str, str0);
	assignValueById(str, str1);
	std::cout << "Length Seq0: " << length(str0) << std::endl;
	std::cout << "Length Seq1: " << length(str1) << std::endl;
	TGraph g(str);

	startTime = clock();
	double score1 = globalAlignment(g, score_type, Hirschberg() );
	duration = (clock() - startTime) / CLOCKS_PER_SEC;
	std::cout << g << std::endl;
	std::cout << "Score: " << score1 << " (Runtime: " << duration << ")" << std::endl;
	std::cout << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

void  Test_Runtime2() {
	typedef String<Dna5> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Id<TStringSet>::Type TId;

	Score<double> score_type = Score<double>(10,-1,-0.5,-2);
	TStringSet str;
	clock_t startTime;
	clock_t duration;
	TString str0;
	TString str1;

	mtRandInit();
	while(true) {
		clear(str);
		clear(str0);
		clear(str1);

		unsigned int stop = (unsigned int) mtRand() % 500 + 1;
		for (unsigned int i=0; i<stop ; ++i) {
			append(str0, (Byte) mtRand() % 5 );
		}
		stop = (unsigned int) mtRand() % 500 + 1;
		for (unsigned int i=0; i<stop; ++i) {
			append(str1,  (Byte) mtRand() % 5);
		}
		std::cout << str0 << std::endl;
		std::cout << str1 << std::endl;


		assignValueById(str, str0);
		assignValueById(str, str1);
		std::cout << "Length Seq0: " << length(str0) << std::endl;
		std::cout << "Length Seq1: " << length(str1) << std::endl;
		TGraph g(str);

		startTime = clock();
		double score = globalAlignment(g, score_type, Gotoh() );
		duration = clock() - startTime;
		std::cout << g << std::endl;
		std::cout << "Score: " << score << " (Runtime: " << duration << ")" << std::endl;
		std::cout << std::endl;

		startTime = clock();
		double score1 = globalAlignment(g, score_type, Hirschberg() );
		duration = clock() - startTime;
		std::cout << g << std::endl;
		std::cout << "Score: " << score1 << " (Runtime: " << duration << ")" << std::endl;
		std::cout << std::endl;

		SEQAN_TASSERT(score == score1)
	}
}


//////////////////////////////////////////////////////////////////////////////

void  Test_CompressedAlphabets() {
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

void Test_TCoffeeTmp() {
//____________________________________________________________________________
// T-Coffee
	typedef String<AminoAcid> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
	
	TString str1 = "DDGGG";
	TString str2 = "NAE";
	TString str3 = "AANGEE";
	TString str4 = "NCCA";
	//TString str1 = "GARFIELDTHELASTFATCAT";
	//TString str2 = "GARFIELDTHEFASTCAT";
	//TString str3 = "GARFIELDTHEVERYFASTCAT";
	//TString str4 = "THEFATCAT";
	TStringSet strSet;
	assignValueById(strSet, str1);
	assignValueById(strSet, str2);
	assignValueById(strSet, str3);
	assignValueById(strSet, str4);
	TGraph lib1(strSet);
	TGraph lib2(strSet);
	TGraph g(strSet);

	// Generate a primary library, i.e., all pairwise alignments
	generatePrimaryLibrary(lib1, AAGroupsDayhoff(), GlobalPairwise_Library() );
	generatePrimaryLibrary(lib2, 2, MUM_Library() );

	// Weighting of libraries (Signal addition)
	combineGraphs(g, lib1, lib2);

	// Triplet library extension
	tripletLibraryExtension(g);

	//// Debug code
	//// Print all possible library matches, i.e., our scoring system
	//typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	//typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	//typedef Infix<Value<TStringSet>::Type>::Type TInfix;
	//TEdgeIterator itEdge(g);
	//for(;!atEnd(itEdge);++itEdge) {
	//	TVertexDescriptor sourceV = sourceVertex(itEdge);
	//	TVertexDescriptor targetV = targetVertex(itEdge);
	//	TInfix inf1 = infix(getValueById(stringSet(g), sequenceId(g, sourceV)),fragmentBegin(g, sourceV), fragmentBegin(g, sourceV) + fragmentLength(g, sourceV));
	//	TInfix inf2 = infix(getValueById(stringSet(g), sequenceId(g, targetV)),fragmentBegin(g, targetV), fragmentBegin(g, targetV) + fragmentLength(g, targetV));
	//	std::cout << "SeqId " << sequenceId(g, sourceV) << ':' << inf1 << " (VertexId: " << sourceV << ')' << std::endl;
	//	std::cout << "SeqId " << sequenceId(g, targetV) << ':' << inf2 << " (VertexId: " << targetV << ')' << std::endl;
	//	std::cout << "Weight " << ':' << getCargo(*itEdge) << std::endl;
	//	std::cout << std::endl;
	//}


	// Calculate a distance matrix using a compressed alphabet or not
	Matrix<double> distanceMatrix; 
	getCommonKmerMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
	kmerToDistanceMatrix(distanceMatrix, FractionalDistance() );

	// Build the neighbor joining tree
	Graph<Tree<double> > njTreeOut;
	slowNjTree(distanceMatrix, njTreeOut);

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void> > gOut(strSet);
	progressiveAlignment(g, njTreeOut, gOut, Hirschberg() );

	// Print the alignment
	std::cout << gOut << std::endl;

}


//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeGarfield() {
//____________________________________________________________________________
// T-Coffee
	typedef String<AminoAcid> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
	
	TString str1 = "GARFIELDTHELASTFATCAT";
	TString str2 = "GARFIELDTHEFASTCAT";
	TString str3 = "GARFIELDTHEVERYFASTCAT";
	TString str4 = "THEFATCAT";
	TStringSet strSet;
	assignValueById(strSet, str1);
	assignValueById(strSet, str2);
	assignValueById(strSet, str3);
	assignValueById(strSet, str4);
	TGraph lib1(strSet);
	TGraph lib2(strSet);
	TGraph g(strSet);

	// Generate a primary library, i.e., all pairwise alignments
	generatePrimaryLibrary(lib1, AAGroupsDayhoff(), GlobalPairwise_Library() );

	fstream strm01; // Alignment graph as dot
	strm01.open(TEST_PATH "my_tcoffee01.dot", ios_base::out | ios_base::trunc);
	write(strm01,lib1,DotDrawing());
	strm01.close();

	generatePrimaryLibrary(lib2, 2, MUM_Library() );

	fstream strm02; // Alignment graph as dot
	strm02.open(TEST_PATH "my_tcoffee02.dot", ios_base::out | ios_base::trunc);
	write(strm02,lib2,DotDrawing());
	strm02.close();

	// Weighting of libraries (Signal addition)
	combineGraphs(g, lib1, lib2);

	fstream strm1; // Alignment graph as dot
	strm1.open(TEST_PATH "my_tcoffee1.dot", ios_base::out | ios_base::trunc);
	write(strm1,g,DotDrawing());
	strm1.close();

	// Triplet library extension
	tripletLibraryExtension(g);

	fstream strm2; // Alignment graph as dot
	strm2.open(TEST_PATH "my_tcoffee2.dot", ios_base::out | ios_base::trunc);
	write(strm2,g,DotDrawing());
	strm2.close();

	// Calculate a distance matrix using a compressed alphabet or not
	Matrix<double> distanceMatrix; 
	getCommonKmerMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
	kmerToDistanceMatrix(distanceMatrix, FractionalDistance() );

	// Build the neighbor joining tree
	Graph<Tree<double> > njTreeOut;
	slowNjTree(distanceMatrix, njTreeOut);

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void> > gOut(strSet);
	progressiveAlignment(g, njTreeOut, gOut, Hirschberg() );

	// Print the alignment
	std::cout << gOut << std::endl;

	fstream strm3; // Alignment graph as dot
	strm3.open(TEST_PATH "my_tcoffee3.dot", ios_base::out | ios_base::trunc);
	write(strm3,gOut,DotDrawing());
	strm3.close();
}


//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeFromRandomSeq() {
//____________________________________________________________________________
// T-Coffee

	while (true) {
	
		// Sequences
		typedef String<AminoAcid> TString;
		StringSet<String<AminoAcid>, Owner<> > origStrSet;

		// Count sequences
		mtRandInit();
		unsigned int seqCount = ((Byte) mtRand() % 10) + 3;
		for(unsigned int i=0;i<seqCount;++i) {
			TString str;

			unsigned int stop = (unsigned int) mtRand() % 20 + 2;
			for (unsigned int i=0; i<stop ; ++i) {
				append(str, (Byte) mtRand() % 10 );
			}	
			assignValueById(origStrSet, str);
		}
		
		// Generate additional primary libraries, e.g., all pairwise alignments
		typedef StringSet<TString, Dependent<> > TStringSet;
		typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
		TStringSet strSet;
		for(unsigned int i = 0; i<length(origStrSet); ++i) {
			std::cout << i << ':' << origStrSet[i] << std::endl;
			appendValue(strSet, origStrSet[i]);
		}
		TGraph g(strSet);

		// Generate primary libraries
		TGraph lib1(strSet);
		TGraph lib2(strSet);
		generatePrimaryLibrary(lib1, AAGroupsDayhoff(), GlobalPairwise_Library() );
		//std::cout << "Pairwise alignments done" << std::endl;

		generatePrimaryLibrary(lib2, 2, MUM_Library() );
		//std::cout << "MUMs done" << std::endl;

		// Weighting of libraries (Signal addition)
		combineGraphs(g, lib1, lib2);
		//std::cout << "Combining graphs done" << std::endl;

		// Triplet library extension
		tripletLibraryExtension(g);
		//std::cout << "Triplet done" << std::endl;

		// Calculate a distance matrix using a compressed alphabet or not
		Matrix<double> distanceMatrix; 
		getCommonKmerMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
		//std::cout << "Kmer done" << std::endl;
		kmerToDistanceMatrix(distanceMatrix, FractionalDistance() );
		//std::cout << "Distance done" << std::endl;


		// Build the neighbor joining tree
		Graph<Tree<double> > njTreeOut;
		slowNjTree(distanceMatrix, njTreeOut);
		//std::cout << "NjTree done" << std::endl;
		
		// Perform a progressive alignment
		Graph<Alignment<TStringSet, void> > gOut(strSet);
		progressiveAlignment(g, njTreeOut, gOut, NeedlemanWunsch() );
		//std::cout << "Alignment done" << std::endl;

		
		// Print the alignment
		std::cout << gOut << std::endl;

		String<char> align;
		bool convAl = convertAlignment(gOut, align);
		SEQAN_TASSERT(convAl);
	}
}
//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeFromFile() {
//____________________________________________________________________________
// T-Coffee

	// Sequences
	typedef String<AminoAcid> TString;
	StringSet<String<AminoAcid>, Owner<> > origStrSet;

	// Count sequences
	unsigned seqCount = 0;
	ifstream file;
	std::stringstream input;
	input << TEST_PATH << "3d_sample2.fasta";
	//input << TEST_PATH << "sproteases_small.fasta";
	file.open(input.str().c_str(), ios_base::in | ios_base::binary);
	if (!file.is_open()) return;
	while (!_streamEOF(file)) {
		String<char> id;
		readID(file, id, Fasta());
		//std::cout << id << std::endl;
		goNext(file, Fasta());
		++seqCount;
	}
	//std::cout << "Number of sequences: " << seqCount << std::endl;

	// Import sequences
	file.clear();
	file.seekg(0, ios_base::beg);
	resize(origStrSet, seqCount);
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 	{
		read(file, origStrSet[i], Fasta());
	}
    file.close();
	std::cout << "Import sequences done" << std::endl;

	// Generate additional primary libraries, e.g., all pairwise alignments
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
	TStringSet strSet;
	for(unsigned int i = 0; i<length(origStrSet); ++i) {
		//assignValueById(strSet, origStrSet, positionToId(origStrSet, i));
		appendValue(strSet, origStrSet[i]);
	}
	TGraph g(strSet);

	// Generate primary libraries
	TGraph lib1(strSet);
	TGraph lib2(strSet);
	generatePrimaryLibrary(lib1, AAGroupsDayhoff(), GlobalPairwise_Library() );
	std::cout << "Pairwise alignments done" << std::endl;

	fstream strm01; // Alignment graph as dot
	strm01.open(TEST_PATH "my_tcoffee01.dot", ios_base::out | ios_base::trunc);
	write(strm01,lib1,DotDrawing());
	strm01.close();

	generatePrimaryLibrary(lib2, 2, MUM_Library() );
	std::cout << "MUMs done" << std::endl;

	fstream strm02; // Alignment graph as dot
	strm02.open(TEST_PATH "my_tcoffee02.dot", ios_base::out | ios_base::trunc);
	write(strm02,lib2,DotDrawing());
	strm02.close();

	// Weighting of libraries (Signal addition)
	combineGraphs(g, lib1, lib2);
	std::cout << "Combining graphs done" << std::endl;

	fstream strm1; // Alignment graph as dot
	strm1.open(TEST_PATH "my_tcoffee1.dot", ios_base::out | ios_base::trunc);
	write(strm1,g,DotDrawing());
	strm1.close();

	// Triplet library extension
	tripletLibraryExtension(g);
	std::cout << "Triplet done" << std::endl;

	fstream strm2; // Alignment graph as dot
	strm2.open(TEST_PATH "my_tcoffee2.dot", ios_base::out | ios_base::trunc);
	write(strm2,g,DotDrawing());
	strm2.close();


	// Calculate a distance matrix using a compressed alphabet or not
	Matrix<double> distanceMatrix; 
	getCommonKmerMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
	std::cout << "Kmer done" << std::endl;
	kmerToDistanceMatrix(distanceMatrix, FractionalDistance() );
	std::cout << "Distance done" << std::endl;

	// Build the neighbor joining tree
	Graph<Tree<double> > njTreeOut;
	slowNjTree(distanceMatrix, njTreeOut);
	std::cout << "NjTree done" << std::endl;

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void> > gOut(strSet);
	progressiveAlignment(g, njTreeOut, gOut, NeedlemanWunsch() );
	std::cout << "Alignment done" << std::endl;

	// Print the alignment
	std::cout << gOut << std::endl;
	
	fstream strm3; // Alignment graph as dot
	strm3.open(TEST_PATH "my_tcoffee3.dot", ios_base::out | ios_base::trunc);
	write(strm3,gOut,DotDrawing());
	strm3.close();
}

//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeFromLibrary() {
//____________________________________________________________________________
// Graph TCoffee

	// Read a t-coffee library: AminoAcid Alphabet
	typedef StringSet<String<AminoAcid>, Owner<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
	TStringSet strSet;
	TGraph g(strSet);

	fstream strm; // Read the library
	strm.open(TEST_PATH "garfield.lib", ios_base::in);
	read(strm,g,TCoffeeLib());
	strm.close();

	/*
	fstream strmW; // Write the library
	strmW.open(TEST_PATH "my_garfield.lib", ios_base::out | ios_base::trunc);
	write(strmW,g,TCoffeeLib());
	strmW.close();

	fstream strm2; // Alignment graph as dot
	strm2.open(TEST_PATH "my_tcoffee.dot", ios_base::out | ios_base::trunc);
	write(strm2,g,DotDrawing());
	strm2.close();
	*/	

	// Generate additional primary libraries, e.g., all pairwise alignments
	typedef StringSet<String<AminoAcid>, Dependent<> > TDependentStringSet;
	typedef Graph<Alignment<TDependentStringSet, unsigned int, Default> > TDependentGraph;
	TStringSet& ownerStrSet = stringSet(g);
	TDependentStringSet depStrSet;
	for(unsigned int i = 0; i<length(ownerStrSet); ++i) {
		assignValueById(depStrSet, ownerStrSet, positionToId(ownerStrSet, i));
	}
	TDependentGraph gAux(depStrSet);
	generatePrimaryLibrary(gAux, AAGroupsDayhoff(), GlobalPairwise_Library() );
	tripletLibraryExtension(gAux);

	// Debug code
	// Print all library matches
	// Note: We do not break it down to the base level
	typedef Iterator<TDependentGraph, EdgeIterator>::Type TEdgeIterator;
	typedef VertexDescriptor<TDependentGraph>::Type TVertexDescriptor;
	typedef Infix<Value<TDependentStringSet>::Type>::Type TInfix;
	TEdgeIterator itEdge(gAux);
	for(;!atEnd(itEdge);++itEdge) {
		TVertexDescriptor sourceV = sourceVertex(itEdge);
		TVertexDescriptor targetV = targetVertex(itEdge);
		TInfix inf1 = infix(getValueById(stringSet(gAux), sequenceId(gAux, sourceV)),fragmentBegin(gAux, sourceV), fragmentBegin(gAux, sourceV) + fragmentLength(gAux, sourceV));
		TInfix inf2 = infix(getValueById(stringSet(gAux), sequenceId(gAux, targetV)),fragmentBegin(gAux, targetV), fragmentBegin(gAux, targetV) + fragmentLength(gAux, targetV));
		std::cout << "SeqId " << sequenceId(gAux, sourceV) << ':' << inf1 << " (VertexId: " << sourceV << ')' << std::endl;
		std::cout << "SeqId " << sequenceId(gAux, targetV) << ':' << inf2 << " (VertexId: " << targetV << ')' << std::endl;
		std::cout << "Weight " << ':' << getCargo(*itEdge) << std::endl;
		std::cout << std::endl;
	}


	// Calculate a distance matrix
	Matrix<double> distanceMatrix; 
	getCommonKmerMatrix(stringSet(gAux), distanceMatrix, 6, AAGroupsDayhoff() );
	kmerToDistanceMatrix(distanceMatrix, FractionalDistance() );

	// Create neighbor joining tree
	Graph<Tree<double> > njTreeOut;
	slowNjTree(distanceMatrix, njTreeOut);


	fstream strm2; // Alignment graph as dot
	strm2.open(TEST_PATH "my_tcoffee.dot", ios_base::out | ios_base::trunc);
	write(strm2,gAux,DotDrawing());
	strm2.close();

	// Perform a progressive alignment
	Graph<Alignment<TDependentStringSet, void> > gOut(depStrSet);
	progressiveAlignment(gAux, njTreeOut, gOut, Hirschberg() );

	//fstream strm2; // Alignment graph as dot
	//strm2.open(TEST_PATH "my_tcoffee.dot", ios_base::out | ios_base::trunc);
	//write(strm2,gOut,DotDrawing());
	//strm2.close();

	/*
	// Read a t-coffee library: Dna Alphabet
	typedef StringSet<String<Dna>, Owner<> > TStringSetDna;
	typedef Graph<Alignment<TStringSetDna, unsigned int, Default> > TGraphDna;
	TStringSetDna strSetDna;
	TGraphDna gDna(strSetDna);

	fstream strmDna; // Read the library
	strmDna.open(TEST_PATH "dna_seq.lib", ios_base::in);
	read(strmDna,gDna,TCoffeeLib());
	strmDna.close();

	fstream strmWDna; // Write the library
	strmWDna.open(TEST_PATH "my_dna_seq.lib", ios_base::out | ios_base::trunc);
	write(strmWDna,gDna,TCoffeeLib());
	strmWDna.close();

	//std::cout << g << std::endl;


	// Calculate a distance matrix
	Matrix<double> distanceMatrixDna; 
	getCommonKmerMatrix(stringSet(gDna), distanceMatrixDna, 6);
	kmerToDistanceMatrix(distanceMatrixDna, TCoffeeDistance() );

	// Create neighbor joining tree
	Graph<Tree<double> > njTreeOutDna;
	slowNjTree(distanceMatrixDna, njTreeOutDna);
	std::cout << njTreeOutDna << std::endl;
	*/

}


//////////////////////////////////////////////////////////////////////////////

void Test_GraphAlignment() {
	Test_NeedlemanWunsch();
	Test_Gotoh();	
	Test_MyersBitVector();
	Test_Hirschberg();
	Test_Nussinov();
}

//////////////////////////////////////////////////////////////////////////////

void Test_TCoffee() {
	//Test_TCoffeeTmp();
	Test_TCoffeeGarfield();
	//Test_TCoffeeFromRandomSeq();
	//Test_TCoffeeFromLibrary(); 

}


}

#endif

