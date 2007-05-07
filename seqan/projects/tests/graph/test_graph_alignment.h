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
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)

	str[0] = "annual";
	str[1] = "annealing";
	score1 = globalAlignment(str, MyersBitVector() );
	score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	score3 = globalAlignment(str, score_type2, Gotoh() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)

	str[0] = "cttagt";
	str[1] = "ttag";
	score1 = globalAlignment(str, MyersBitVector() );
	score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	score3 = globalAlignment(str, score_type2, Gotoh() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)

	str[0] = "ttag";
	str[1] = "cttccagt";
	score1 = globalAlignment(str, MyersBitVector() );
	score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	score3 = globalAlignment(str, score_type2, Gotoh() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)
	 
	str[0] = "ttagttagttagttagttagttagttagttagttagttagttagttagttagttagttag";
	str[1] = "cttccagtcttccagtcttccagtcttccagtcttccagtcttccagtcttccagtcttccagt";
	score1 = globalAlignment(str, MyersBitVector() );
	score2 = globalAlignment(str, score_type, NeedlemanWunsch() );
	score3 = globalAlignment(str, score_type2, Gotoh() );
	SEQAN_TASSERT((-1) * score1 == score2)
	SEQAN_TASSERT(score2 == score3)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Hirschberg() {
	// ToDo!!!
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("annual");	assignValueById(str, str0);
	TString str1("annealing"); assignValueById(str, str1);
	//TGraph g(str);
	//Score<int> score_type = Score<int>(0,-1,-1,0);
	//int score1 = globalAlignment(g, Hirschberg_MyersBitVector() );
	//int score2 = globalAlignment(g, score_type, Hirschberg_NeedlemanWunsch() );
	//Score<int> score_type2 = Score<int>(0,-1,-1,-1);
	//int score3 = globalAlignment(g, score_type2, Hirschberg_Gotoh() );
	//SEQAN_TASSERT((-1) * score1 == score2)
	//SEQAN_TASSERT(score2 == score3)
}



/*
//////////////////////////////////////////////////////////////////////////////

void  Test_Runtime() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	Score<double> score_type = Score<double>(5,-4,-0.5,-10);
	Score<double> score_type2 = Score<double>(5,-4,-0.5,-9.5);
	double score;
	TStringSet str;
	clock_t startTime;
	clock_t duration;

	TString str0;
	fstream strm_in;
	strm_in.open(TEST_PATH "a.fasta", ios_base::in | ios_base::binary);
	read(strm_in, str0, Fasta());
	strm_in.close();
	assignValueById(str, str0);

	TString str1;
	fstream strm_in1;
	strm_in1.open(TEST_PATH "b.fasta", ios_base::in | ios_base::binary);
	read(strm_in1, str1, Fasta());
	strm_in1.close();
	assignValueById(str, str1);

	std::cout << "Length Seq0: " << length(str0) << std::endl;
	std::cout << "Length Seq1: " << length(str1) << std::endl;

	TGraph g(str);
	startTime = clock();
	score = globalAlignment(g, score_type, Gotoh() );
	duration = clock() - startTime;
	std::cout << g << std::endl;
	std::cout << "Score: " << score << " (Runtime: " << duration << ")" << std::endl;
	std::cout << std::endl;

	Align< String<Dna>, ArrayGaps> ali;
	resize(rows(ali), 2);
	assignSource(row(ali, 0), str[0]);
	assignSource(row(ali, 1), str[1]);
	startTime = clock();
	score = needlemanWunsch(ali,score_type2);
	duration = clock() - startTime;
	std::cout << ali << std::endl;
	std::cout << "Score: " << score << " (Runtime: " << duration << ")" << std::endl;
	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void  Test_LargeAlignment() {
	typedef String<Dna> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	Score<double> score_type = Score<double>(5,-4,-0.5,-10);
	double score;
	TStringSet str;
	clock_t startTime;
	clock_t duration;

	TString str0;
	fstream strm_in;
	strm_in.open(TEST_PATH "a.fasta", ios_base::in | ios_base::binary);
	read(strm_in, str0, Fasta());
	strm_in.close();
	assignValueById(str, str0);

	TString str1;
	fstream strm_in1;
	strm_in1.open(TEST_PATH "b.fasta", ios_base::in | ios_base::binary);
	read(strm_in1, str1, Fasta());
	strm_in1.close();
	assignValueById(str, str1);

	std::cout << "Length Seq0: " << length(str0) << std::endl;
	std::cout << "Length Seq1: " << length(str1) << std::endl;
	startTime = clock();
	score = globalAlignment(std::cout, str, score_type, Gotoh() );
	duration = clock() - startTime;
	std::cout << "Score: " << score << " (Runtime: " << duration << ")" << std::endl;
	std::cout << std::endl;
}
*/

//////////////////////////////////////////////////////////////////////////////

int Test_ConvertSequences(String<char> const in_path, String<char> const in_file, String<char> const path, String<char> const file_prefix) {
	typedef String<Dna5, External<ExternalConfig<File<>, 64*1024> > > TString;
	
	// count sequences
	unsigned seqCount = 0;

	ifstream file;
	std::stringstream input;
	input << in_path << in_file;
	file.open(input.str().c_str(), ios_base::in | ios_base::binary);
	if (!file.is_open()) return 0;
	while (!_streamEOF(file)) {
		String<char> id;
		readID(file, id, Fasta());
		std::cout << id << std::endl;
		goNext(file, Fasta());
		++seqCount;
	}
	std::cout << "Number of sequences: " << seqCount << std::endl;

	// import sequences
	file.clear();
	file.seekg(0, ios_base::beg);
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 
	{
		TString str;
		//String<TraceBack, External<> > trace;
		//open(trace, "D:\\seqan.dat");
		std::stringstream s;
		s << path << file_prefix << i;
		open(str, s.str().c_str());
		read(file, str, Fasta());
	}
    file.close();

	return seqCount;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TVal1, typename TVal2>
inline bool 
Test_ReadSequences(String<char> const path, String<char> const file_prefix, TStringSet& str, TVal1 const start, TVal2 const nseq) {
	for(unsigned i = start; i < start + nseq; ++i) {
		std::stringstream s;
		s << path << file_prefix << i - start;
		bool f = open(str[i], s.str().c_str());
		if (!f) return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

void Test_MatchRefinement() {
	typedef String<Dna5, External<ExternalConfig<File<>, 64*1024> > > TString;
	typedef StringSet<TString, Owner<> > TStringSet;
	typedef Id<TStringSet>::Type TId;
	typedef Size<TStringSet>::Type TSize;
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\matches\\");
	String<char> out_path("Z:\\matches\\out\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/matches/");
	String<char> out_path("/home/takifugu/rausch/matches/out/");
#endif
	

	TSize hSeq = 24;
	TSize wSeq = 24;
	TSize bSeq = 24;
	

	// Convert all sequences only once
	//TSize tmp = 0;
	//tmp = Test_ConvertSequences(in_path, "HUREF6CHROM.fasta",out_path,"H.chr.");
	//if (tmp != hSeq) {
	//	std::cerr << "Did not read all HUREF sequences." << std::endl;
	//	exit(1);
	//}
	//tmp = Test_ConvertSequences(in_path, "WGSACHROM.fasta",out_path,"W.chr.");
	//if (tmp != wSeq) {
	//	std::cerr << "Did not read all WGSA sequences." << std::endl;
	//	exit(1);
	//}
	//tmp = Test_ConvertSequences(in_path, "B36LCCHROM.fasta",out_path,"B.chr.");
	//if (tmp != bSeq) {
	//	std::cerr << "Did not read all B36LC sequences." << std::endl;
	//	exit(1);
	//}

	// Read all the sequences
	TStringSet str;
	resize(str, hSeq + wSeq + bSeq);
	bool f = Test_ReadSequences(out_path,"H.chr.", str, 0, hSeq);
	if (!f) {
		std::cerr << "Error importing HUREF sequences." << std::endl;
		exit(1);
	}
	f = Test_ReadSequences(out_path,"W.chr.", str, hSeq, wSeq);
	if (!f) {
		std::cerr << "Error importing WGSA sequences." << std::endl;
		exit(1);
	}
	f = Test_ReadSequences(out_path,"B.chr.", str, hSeq + wSeq, bSeq);
	if (!f) {
		std::cerr << "Error importing B36LC sequences." << std::endl;
		exit(1);
	}

	// Build a map:
	// SeqId -> Identifier
	typedef std::map<TId, String<char> > TIdToNameMap;
	TIdToNameMap idToName;
	for(TId i = 0; i < length(str); ++i) {
		TSize index = 0;
		std::stringstream s;
		if (i < 24) {
			s << "H";
			index = i;
		}
		else if (i < 48) {
			s << "W";
			index = i - hSeq;
		}
		else {
			s << "B";
			index = i - (hSeq + wSeq);
		}
		s << ":" << index;
		idToName.insert(std::make_pair(i, s.str().c_str()));
	}

	// Just to check that everything worked
	std::cout << "Number of sequences: " << length(str) << std::endl;
	for(TIdToNameMap::const_iterator pos =  idToName.begin(); pos != idToName.end(); ++pos) {
		std::cout << pos->second << ") ";
		for(TSize i=0; i<10;++i) {
			std::cout << str[pos->first][i];
		}
		std::cout << std::endl;
	}

	// Read the matches
	typedef StringSet<TString, Dependent<> > TMatchStringSet;
	String<Graph<Alignment<TMatchStringSet, void, Default()> > > matches;
	for(TSize i = 0; i<1; ++i) {
		fstream strm; 
		std::stringstream s;
		if (i==0 ) s << in_path << "TvsT.atac";
		else if (i==1 ) s << in_path << "BvsH.atac";
		else if (i==2 ) s << in_path << "BvsW.atac";
		else if (i==3 ) s << in_path << "WvsH.atac";
		strm.open(s.str().c_str(), ios_base::in);
		read(strm, matches, str, 0, AtacMatches());
		strm.close();
	}

	// Print all matches
	std::cout << "Number of matches: " << length(matches) << std::endl;
	typedef Infix<TString>::Type TInfix;
	for(TSize i = 0; i < length(matches); ++i) {
		TIdToNameMap::const_iterator pos1 =  idToName.find(sequenceId(matches[i],0));
		TIdToNameMap::const_iterator pos2 =  idToName.find(sequenceId(matches[i],1));
		std::cout << matches[i];
		std::cout << pos1->second << ") ";
		TInfix infix1 = label(matches[i],0);
		std::cout << infix1 << std::endl;
		std::cout << pos2->second << ") ";
		TInfix infix2 = label(matches[i],1);
		std::cout << infix2 << std::endl;
		std::cout << std::endl;
	}
	/*
	Score<int> score_type = Score<int>(1,-1,-2,0) ;
	typedef Graph<Alignment<TStringSet> > TAliGraph;
	
	TAliGraph ali_graph(str);

	matchRefinement(matches,str,score_type,ali_graph);//,StoreEdges());
	std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
	std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
	*/
}


//////////////////////////////////////////////////////////////////////////////

void Test_Fragment() {
	// Test Fragment
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef	Id<TStringSet>::Type TId;
	
	TStringSet str;
	TString str0("annual");	assignValueById(str, str0);
	TString str1("anneal"); assignValueById(str, str1);

	// Fragment: SeqId1, Begin1, SeqId2, Begin2, Length of Fragment
	Fragment<> f(0,4,1,4,2);
	SEQAN_TASSERT(f.seqId1 == 0)
	SEQAN_TASSERT(f.begin1 == 4)
	SEQAN_TASSERT(f.seqId2 == 1)
	SEQAN_TASSERT(f.begin2 == 4)
	SEQAN_TASSERT(f.len == 2)
	SEQAN_TASSERT(fragmentBegin(f, 0) == 4)
	SEQAN_TASSERT(fragmentBegin(f, 1) == 4)
	SEQAN_TASSERT(fragmentLength(f, 0) == 2)
	SEQAN_TASSERT(fragmentLength(f, 1) == 2)
	SEQAN_TASSERT(sequenceId(f, 0) == 0)
	SEQAN_TASSERT(sequenceId(f, 1) == 1)
	SEQAN_TASSERT(label(f, str, 0) == "al")
	SEQAN_TASSERT(label(f, str, 1) == "al")
	SEQAN_TASSERT(getProjectedPosition(f, 0, 5) == 5)
	SEQAN_TASSERT(getProjectedPosition(f, 1, 5) == 5)

	// Fragment: SeqId1, Begin1, SeqId2, Begin2, Length of Fragment
	Fragment<> f2(0,0,1,4,1);
	SEQAN_TASSERT(f2.seqId1 == 0)
	SEQAN_TASSERT(f2.begin1 == 0)
	SEQAN_TASSERT(f2.seqId2 == 1)
	SEQAN_TASSERT(f2.begin2 == 4)
	SEQAN_TASSERT(f2.len == 1)
	SEQAN_TASSERT(fragmentBegin(f2, 0) == 0)
	SEQAN_TASSERT(fragmentBegin(f2, 1) == 4)
	SEQAN_TASSERT(fragmentLength(f2, 0) == 1)
	SEQAN_TASSERT(fragmentLength(f2, 1) == 1)
	SEQAN_TASSERT(label(f2, str, 0) == "a")
	SEQAN_TASSERT(label(f2, str, 1) == "a")
	SEQAN_TASSERT(getProjectedPosition(f2, 0, 0) == 4)
	SEQAN_TASSERT(getProjectedPosition(f2, 1, 4) == 0)

	// The same stuff on an alignment graph
	TGraph g(str);
	TVertexDescriptor vert1 = addVertex(g,0,4,2);
	TVertexDescriptor vert2 = addVertex(g,1,4,2);
	addEdge(g, vert1, vert2);
	SEQAN_TASSERT(fragmentBegin(g, vert1) == 4)
	SEQAN_TASSERT(fragmentBegin(g, vert2) == 4)
	SEQAN_TASSERT(fragmentLength(g, vert1) == 2)
	SEQAN_TASSERT(fragmentLength(g, vert2) == 2)
	SEQAN_TASSERT(sequenceId(f, vert1) == 0)
	SEQAN_TASSERT(sequenceId(f, vert2) == 1)
	SEQAN_TASSERT(label(g, vert1) == "al")
	SEQAN_TASSERT(label(g, vert2) == "al")
	SEQAN_TASSERT(getProjectedPosition(g, 0, 5) == 5)
	SEQAN_TASSERT(getProjectedPosition(g, 1, 5) == 5)

	TGraph g2(str);
	TVertexDescriptor v1 = addVertex(g2,0,0,1);
	TVertexDescriptor v2 = addVertex(g2,1,4,1);
	addEdge(g2, v1, v2);
	SEQAN_TASSERT(fragmentBegin(g2, v1) == 0)
	SEQAN_TASSERT(fragmentBegin(g2, v2) == 4)
	SEQAN_TASSERT(fragmentLength(g2, v1) == 1)
	SEQAN_TASSERT(fragmentLength(g2, v2) == 1)
	SEQAN_TASSERT(label(g2, v1) == "a")
	SEQAN_TASSERT(label(g2, v2) == "a")
	SEQAN_TASSERT(getProjectedPosition(g2, 0, 0) == 4)
	SEQAN_TASSERT(getProjectedPosition(g2, 1, 4) == 0)
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

void Test_TCoffee() {
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

	// Generate additional primary libraries
	// Just slow-pair at the moment
	TGraph gAux(stringSet(g));
	generatePrimaryLibrary(gAux, AAGroupsDayhoff() );

	// Calculate a distance matrix
	Matrix<double> distanceMatrix; 
	getCommonKmerMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
	kmerToDistanceMatrix(distanceMatrix, FractionalDistance() );

	// Create neighbor joining tree
	Graph<Tree<double> > njTreeOut;
	slowNjTree(distanceMatrix, njTreeOut);
	std::cout << njTreeOut << std::endl;


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


void Test_GraphAlignment() {
	Test_NeedlemanWunsch();
	Test_Gotoh();	
	Test_MyersBitVector();
	Test_Hirschberg();
}

}

#endif

