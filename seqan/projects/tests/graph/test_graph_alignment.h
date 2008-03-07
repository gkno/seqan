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
	int score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 3 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 3 == score2)
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;
	int score3 = globalAlignment(matches, stringSet(g), score_type, NeedlemanWunsch() );
	SEQAN_TASSERT(length(matches) == 1)
	SEQAN_TASSERT(label(matches[0], stringSet(g), 0) == "annual")
	SEQAN_TASSERT(label(matches[0], stringSet(g), 1) == "anneal")
	SEQAN_TASSERT(score3 == score)
	//std::cout << g << std::endl;
	//score2 = globalAlignment(std::cout, str, score_type, NeedlemanWunsch() );
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 3 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 3 == score2)
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6  == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 12 == score2)
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score  == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 6 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,true>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,true,false>(), NeedlemanWunsch() );
	SEQAN_TASSERT(score + 12 == score2)
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

void Test_Gotoh() {
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
	Score<int> score_type = Score<int>(1,-1,-1,-2);
	int score = globalAlignment(g, score_type, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	int score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
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
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score + 4 == score2)
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
	score2 = globalAlignment(g, score_type, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score + 2 == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,true,false>(), Gotoh() );
	SEQAN_TASSERT(score + 4 == score2)
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
	Score<int> score_type2 = Score<int>(10,-1,-1,-2);
	score = globalAlignment(g, score_type2, Gotoh());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "aa")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 1)) == "aa")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 6)
	score2 = globalAlignment(stringSet(g), score_type2, AlignConfig<true,false,false,false>(), Gotoh() );
	SEQAN_TASSERT(score + 3 == score2)
	score2 = globalAlignment(stringSet(g), score_type2, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
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
	Score<int> score_type = Score<int>(1,-1,-1,-2);
	int score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "gt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "gt")
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(numVertices(g) == 5)
	int score2 = globalAlignment(stringSet(g), score_type, Hirschberg() );
	SEQAN_TASSERT(score == score2)
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
	Score<int> score_type2 = Score<int>(10,-1,-1,-2);
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

void Test_SmithWaterman() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("gctctgcgaata"); assignValueById(str, str0);
	TString str1("cgttgagatact"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(2,-1,-2,-2);
	int score = localAlignment(g, score_type, SmithWaterman());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "tgcg")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 8)) == "a")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 9)) == "ata")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "tgag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "ata")
	SEQAN_TASSERT(numEdges(g) == 2)
	int score2 = localAlignment(stringSet(g), score_type, SmithWaterman());
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//score2 = localAlignment(std::cout, str, score_type, SmithWaterman() );
	//SEQAN_TASSERT(score == score2)
	
	str[0] = "TTGACACCCTCCCAATTGTA";
	str[1] = "ACCCCAGGCTTTACACAT";
	assignStringSet(g, str);
	score = localAlignment(g, score_type, SmithWaterman());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "TTGACAC")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 9)) == "TTTACAC")
	SEQAN_TASSERT(numEdges(g) == 1)
	score2 = localAlignment(stringSet(g), score_type, SmithWaterman());
	SEQAN_TASSERT(score == score2)
	//std::cout << g << std::endl;
	//score2 = localAlignment(std::cout, str, score_type, SmithWaterman() );
	//SEQAN_TASSERT(score == score2)
}
	

//////////////////////////////////////////////////////////////////////////////

void Test_SmithWatermanClump() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;

	TStringSet str;
	TString str0("gctctgcgaata"); assignValueById(str, str0);
	TString str1("cgttgagatact"); assignValueById(str, str1);
	TGraph g(str);
	Score<int> score_type = Score<int>(2,-1,-2,-2);
	localAlignment(g, score_type, SmithWatermanClump());
	
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "cgt")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 3)) == "t")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 4)) == "ga")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 6)) == "g")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 7)) == "at")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 9)) == "a")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 10)) == "c")
	//SEQAN_TASSERT(label(g, findVertex(g, 1, 11)) == "t")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "g")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 1)) == "ct")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 3)) == "c")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "t")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 5)) == "gc")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 7)) == "g")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 8)) == "a")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 9)) == "at")
	//SEQAN_TASSERT(label(g, findVertex(g, 0, 11)) == "a")

	//SEQAN_TASSERT(numEdges(g) == 9)
	//SEQAN_TASSERT(numVertices(g) == 17)
		
}

/*

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

*/
//////////////////////////////////////////////////////////////////////////////

void Test_GraphAlignment() {
	// Global alignments
	Test_NeedlemanWunsch();
	Test_Gotoh();	
	Test_MyersBitVector();
	Test_Hirschberg();

	// Local alignments
	Test_SmithWaterman();
	Test_SmithWatermanClump();

	// To test the runtime
	//Test_Runtime2();

	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_config.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_needleman_wunsch.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_gotoh.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_myers.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_hirschberg.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_smith_waterman.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_smith_waterman_clump.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_interface.h");
}

}

#endif

