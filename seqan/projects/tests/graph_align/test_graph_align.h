#ifndef SEQAN_HEADER_TEST_GRAPH_ALIGNMENT_H
#define SEQAN_HEADER_TEST_GRAPH_ALIGNMENT_H

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

	str[0] = "cat";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ThisisGarfieldthe")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 18)) == "cat")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "cat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	std::cout << g << std::endl;

	str[0] = "ThisisGarfieldthecat";
	str[1] = "cat";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, NeedlemanWunsch());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ThisisGarfieldthe")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 18)) == "cat")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "cat")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	std::cout << g << std::endl;
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
	int score3 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	std::cout << g << std::endl;

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
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<true,false,false,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,true,false,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,false>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,false,true>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	score2 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), Gotoh() );
	SEQAN_TASSERT(score == score2)
	score3 = globalAlignment(stringSet(g), score_type, AlignConfig<false,false,true,true>(), -1 * length(stringSet(g)[1]), length(stringSet(g)[0]), BandedGotoh());
	SEQAN_TASSERT(score2 == score3)
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;


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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

	str[0] = "tagt";
	str[1] = "attagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 2)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "at")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	std::cout << g << std::endl;

	str[0] = "attagt";
	str[1] = "tagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 2)) == "tagt")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "at")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	std::cout << g << std::endl;

	str[0] = "ttagt";
	str[1] = "ttag";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 0, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	std::cout << g << std::endl;

	str[0] = "ttag";
	str[1] = "ttagt";
	assignStringSet(g, str);
	score = globalAlignment(g, score_type, Hirschberg());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 0)) == "ttag")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 4)) == "t")
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(numVertices(g) == 3)
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;

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
	std::cout << g << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScore> 
inline typename Value<TScore>::Type
__sumOfPairsScore(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				TScore const& score_type)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;

	// Convert the graph
	String<char> mat;
	convertAlignment(g, mat);
	char gapChar = gapValue<char>();

	TScoreValue gap = scoreGapExtend(score_type);
	TScoreValue gapOpen = scoreGapOpen(score_type);
	TSize nseq = length(stringSet(g));
	TSize len = length(mat) / nseq;
	
	bool gapOpeni = false;
	bool gapOpenj = false;
	TScoreValue totalScore = 0;
	for(TSize i = 0; i<nseq-1; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			for(TSize k=0;k<len; ++k) {
				if (value(mat, i*len+k) != gapChar) {
					if (value(mat, j*len + k) != gapChar) {
						gapOpeni = false;
						gapOpenj = false;
						totalScore += score(const_cast<TScore&>(score_type), TAlphabet(value(mat, i*len+k)), TAlphabet(value(mat, j*len + k)));
					} else {
						if (gapOpenj) {
							totalScore += gap;
						} else {
							gapOpenj = true;
							totalScore += gapOpen;
						}
					}
				} else if (value(mat, j*len + k) != gapChar) {
						if (gapOpeni) {
							totalScore += gap;
						} else {
							gapOpeni = true;
							totalScore += gapOpen;
						}
				}
			}
		}
	}
	return totalScore;
}

//////////////////////////////////////////////////////////////////////////////

template<bool TTop, bool TLeft, bool TRight, bool TBottom>
void _Test_AllAgainstAll(AlignConfig<TTop, TLeft, TRight, TBottom> ac) {
	typedef unsigned int TSize;
	typedef int TScore;

	mtRandInit();
	for(TSize i = 0; i < 10; ++i) {
		typedef String<Dna> TDnaString;
		typedef StringSet<TDnaString, Dependent<> > TStringSet;
		typedef Graph<Alignment<TStringSet, void> > TGraph;
		TSize lenN = mtRand() % 20 + 1;
		TSize lenM = mtRand() % 20 + 1;
		TDnaString dna1;
		TDnaString dna2;
		for(TSize i = 0; i<lenN; ++i) appendValue(dna1, mtRand() % 4);
		for(TSize j = 0; j<lenM; ++j) appendValue(dna2, mtRand() % 4);
		//dna1 = "TCG";
		//dna2 = "GCGGAG";
		TStringSet str;
		appendValue(str, dna1);
		appendValue(str, dna2);
		Score<int> score_type = Score<int>(4,-3,-2,-2);
		TGraph g(str);
		globalAlignment(g, score_type, ac, Gotoh() );
		// sumOfPairsScore with AlignConfig is missing!!!
		TScore sc1 = __sumOfPairsScore(g, score_type);
		//std::cerr << sc1 << std::endl;
		//std::cerr << g << std::endl;
		clear(g);
		assignStringSet(g, str);
		globalAlignment(g, score_type, ac, NeedlemanWunsch() );
		// sumOfPairsScore with AlignConfig is missing!!!
		TScore sc2 = __sumOfPairsScore(g, score_type);
		//std::cerr << sc2 << std::endl;
		//std::cerr << g << std::endl;
		if (sc1 != sc2) {
			std::cerr << "Randomized test failed:" << std::endl;
			std::cerr << "Seq1: " << dna1 << std::endl;
			std::cerr << "Seq2: " << dna2 << std::endl;
			std::cerr << "AlignConfig: " << TTop << ',' << TLeft << ',' << TRight << ',' << TBottom << std::endl;
			exit(0);
		}
	}
}




//////////////////////////////////////////////////////////////////////////////////


template<bool TTop, bool TLeft, bool TRight, bool TBottom>
void __AllAgainstAll(AlignConfig<TTop, TLeft, TRight, TBottom> ac) {
	typedef unsigned int TSize;
	typedef int TScore;

	mtRandInit();
	for(TSize i = 0; i < 10; ++i) {
		typedef Dna5Q TAlphabet;
		typedef String<TAlphabet> TSequence;
		
		TSize lenN = mtRand() % 5 + 1;
		TSize lenM = mtRand() % 5 + 1;
		TSequence dna1;
		TSequence dna2;
		for(TSize i = 0; i<lenN; ++i) appendValue(dna1, TAlphabet(mtRand() % ValueSize<TAlphabet>::VALUE));
		for(TSize j = 0; j<lenM; ++j) appendValue(dna2, TAlphabet(mtRand() % ValueSize<TAlphabet>::VALUE));
		//dna1 = "T";
		//dna2 = "CAG";
		
		TSize len1 = length(dna1);
		TSize len2 = length(dna2);
		typedef StringSet<TSequence, Dependent<> > TStringSet;
		typedef Graph<Alignment<TStringSet, void> > TGraph;
		TStringSet str;
		appendValue(str, dna1);
		appendValue(str, dna2);

		int matchScore =  mtRand() % 5;
		int misMatchScore =  -1 * (int) (mtRand() % 5);
		int gapScore =  -1 * (int) (mtRand() % 5);
		if (gapScore > misMatchScore) gapScore = misMatchScore;
		//matchScore = 1;
		//misMatchScore = -4;
		//gapScore = -4;

		Score<int> score_type = Score<int>(matchScore,misMatchScore,gapScore,gapScore);
		typedef String<Fragment<> > TFragmentString;
		TFragmentString matches;
		globalAlignment(matches, str, score_type, ac, NeedlemanWunsch());
		int lowDiag = length(dna1);
		int highDiag = -1 * length(dna2);
		typedef typename Iterator<TFragmentString, Standard>::Type TFragIter;
		TFragIter itFrag = begin(matches, Standard());
		TFragIter itFragEnd = end(matches, Standard());
		for(;itFrag != itFragEnd; ++itFrag) {
			int newDiff = (int) itFrag->begin1 - (int) itFrag->begin2;
			if (newDiff > highDiag) highDiag = newDiff;
			if (newDiff < lowDiag) lowDiag = newDiff;
		}
		itFrag = begin(matches, Standard());
		if (itFrag != itFragEnd) {
			int pos1 = (int) itFrag->begin1 + (int) itFrag->len;
			int pos2 = (int) itFrag->begin2 + (int) itFrag->len;
			if (highDiag < ((int) len1 - pos2)) highDiag = ((int) len1 - pos2);
			if (lowDiag > (pos1 - (int) len2)) lowDiag = (pos1 - (int) len2);
			itFrag = itFragEnd;
			--itFrag;
			if (highDiag < (int) itFrag->begin1) highDiag = (int) itFrag->begin1;
			if (lowDiag > -1 * (int) itFrag->begin2) lowDiag = -1 * (int) itFrag->begin2;
		} else {
			lowDiag = -1 * length(dna2);
			highDiag = length(dna1);
		}
		TGraph g1(str);
		globalAlignment(g1, score_type, ac, NeedlemanWunsch());
		int sc1 = _pairWiseSumOfPairsScore(g1, score_type, ac);
		//std::cerr << g1 << std::endl;
		TGraph g2(str);
		globalAlignment(g2, score_type, ac, lowDiag, highDiag, BandedNeedlemanWunsch());
		int sc2 = _pairWiseSumOfPairsScore(g2, score_type, ac);
		TGraph g3(str);
		globalAlignment(g3, score_type, ac, -1 * length(dna2), length(dna1), BandedNeedlemanWunsch());
		int sc3 = _pairWiseSumOfPairsScore(g3, score_type, ac);
		TGraph g4(str);
		globalAlignment(g4, score_type, ac, Gotoh());
		int sc4 = _pairWiseSumOfPairsScore(g4, score_type, ac);
		TGraph g5(str);
		globalAlignment(g5, score_type, ac, lowDiag, highDiag, BandedGotoh());
		int sc5 = _pairWiseSumOfPairsScore(g5, score_type, ac);
		TGraph g6(str);
		globalAlignment(g6, score_type, ac, -1 * length(dna2), length(dna1), BandedGotoh());
		int sc6 = _pairWiseSumOfPairsScore(g6, score_type, ac);
		if ((sc1 != sc2) || (sc2 != sc3) || (sc3 != sc4) || (sc4 != sc5) || (sc5 != sc6)) {
			std::cerr << "Randomized test failed:" << std::endl;
			std::cerr << "Seq1: " << dna1 << std::endl;
			std::cerr << "Seq2: " << dna2 << std::endl;
			std::cerr << "AlignConfig: " << TTop << ',' << TLeft << ',' << TRight << ',' << TBottom << std::endl;
			std::cerr << "Scores: (Matchscore: " << matchScore << ", Mismatchscore: " << misMatchScore << ", Gapscore: " << gapScore << ')' << std::endl;
			std::cerr << g1 << std::endl;
			std::cerr << "Score: " << sc1 << std::endl;
			std::cerr << g2 << std::endl;
			std::cerr << "Score: " << sc2 << std::endl;
			std::cerr << "Diagonals: " << lowDiag << ',' << highDiag << std::endl;
			std::cerr << g3 << std::endl;
			std::cerr << "Score: " << sc3 << std::endl;
			std::cerr << g4 << std::endl;
			std::cerr << "Score: " << sc4 << std::endl;
			std::cerr << g5 << std::endl;
			std::cerr << "Score: " << sc5 << std::endl;
			std::cerr << "Diagonals: " << lowDiag << ',' << highDiag << std::endl;
			std::cerr << g6 << std::endl;
			std::cerr << "Score: " << sc6 << std::endl;
			exit(-1);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

void Test_AllAgainstAll() {
	//_Test_AllAgainstAll(AlignConfig<false,false,false,false>() );

	__AllAgainstAll(AlignConfig<false,false,false,false>() );
	__AllAgainstAll(AlignConfig<false,false,false,true>() );
	__AllAgainstAll(AlignConfig<false,false,true,false>() );
	__AllAgainstAll(AlignConfig<false,false,true,true>() );
	__AllAgainstAll(AlignConfig<false,true,false,false>() );
	__AllAgainstAll(AlignConfig<false,true,false,true>() );
	__AllAgainstAll(AlignConfig<false,true,true,false>() );
	__AllAgainstAll(AlignConfig<false,true,true,true>() );
	__AllAgainstAll(AlignConfig<true,false,false,false>() );
	__AllAgainstAll(AlignConfig<true,false,false,true>() );
	__AllAgainstAll(AlignConfig<true,false,true,false>() );
	__AllAgainstAll(AlignConfig<true,false,true,true>() );
	__AllAgainstAll(AlignConfig<true,true,false,false>() );
	__AllAgainstAll(AlignConfig<true,true,false,true>() );
	__AllAgainstAll(AlignConfig<true,true,true,false>() );
	__AllAgainstAll(AlignConfig<true,true,true,true>() );
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
	std::cout << g << std::endl;
	
	str[0] = "TTGACACCCTCCCAATTGTA";
	str[1] = "ACCCCAGGCTTTACACAT";
	assignStringSet(g, str);
	score = localAlignment(g, score_type, SmithWaterman());
	SEQAN_TASSERT(label(g, findVertex(g, 0, 0)) == "TTGACAC")
	SEQAN_TASSERT(label(g, findVertex(g, 1, 9)) == "TTTACAC")
	SEQAN_TASSERT(numEdges(g) == 1)
	score2 = localAlignment(stringSet(g), score_type, SmithWaterman());
	SEQAN_TASSERT(score == score2)
	std::cout << g << std::endl;
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
	Score<int> score_type = Score<int>(2,-1,-2,-2);
	String<Fragment<> > matches;
	String<int> scores;
	multiLocalAlignment(str, matches, scores, score_type, 4, SmithWatermanClump() );
	std::cout << length(matches) << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphAlignment() {
	// Global alignments
	//Test_NeedlemanWunsch();
	//Test_Gotoh();	
	//Test_Hirschberg();
	Test_AllAgainstAll();

	// Local alignments
	//Test_SmithWaterman();
	//Test_SmithWatermanClump();

	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_config.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_needleman_wunsch.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_gotoh.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_hirschberg.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_smith_waterman.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_smith_waterman_clump.h");
	//debug::verifyCheckpoints("projects/library/seqan/graph_align/graph_align_smith_waterman_island.h");
}

}

#endif

