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
	int score = needlemanWunsch(g, SimpleScore() );
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
	std::cout << std::endl;

	str[0] = "annealing";
	str[1] = "annual";
	assignStringSet(g, str);
	score = needlemanWunsch(g, SimpleScore() );
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
	std::cout << std::endl;

	str[0] = "ThisisGarfieldthecat";
	str[1] = "Garfield";
	assignStringSet(g, str);
	score = needlemanWunsch(g, SimpleScore() );
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
	std::cout << std::endl;

	str[0] = "Garfield";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = needlemanWunsch(g, SimpleScore() );
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
	std::cout << std::endl;

	str[0] = "cat";
	str[1] = "ThisisGarfieldthecat";
	assignStringSet(g, str);
	score = needlemanWunsch(g, SimpleScore() );
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
	std::cout << std::endl;

	str[0] = "ThisisGarfieldthecat";
	str[1] = "cat";
	assignStringSet(g, str);
	score = needlemanWunsch(g, SimpleScore() );
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
	std::cout << std::endl;

	str[0] = "ThisisGarfieldthecat";
	str[1] = "MyGarfieldcat";
	assignStringSet(g, str);
	score = needlemanWunsch(g, SimpleScore() );
	std::cout << g << std::endl;
	std::cout << "Score: " << score << std::endl;
	std::cout << std::endl;
}


}

#endif

