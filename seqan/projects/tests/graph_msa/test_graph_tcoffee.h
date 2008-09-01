#ifndef SEQAN_HEADER_TEST_GRAPH_TCOFFEE_H
#define SEQAN_HEADER_TEST_GRAPH_TCOFFEE_H

using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template<typename TTag>
void 
Test_UpgmaGuideTree() {
	typedef unsigned int TSize;
	
	mtRandInit();
	for(TSize i = 0; i < 10; ++i) {
		// Set-up a sparse distance matrix
		TSize n = mtRand() % 10 + 2;
		Graph<Undirected<double> > distGraph;
		String<double> distMatrix;
		fill(distMatrix, n * n, 0);
		TSize all = (n * (n - 1)) / 2;
		typedef std::set<double> TDistanceSet;
		typedef TDistanceSet::iterator TSetIt;
		TDistanceSet distances;
		String<double> myDist;
		typedef Iterator<String<double> >::Type TStringIter;
		while (distances.size() < all) {
			double newVal = mtRand() % 1000000;
			//double newVal = mtRand() % 100;
			if (distances.insert(newVal).second) {
				appendValue(myDist, newVal);
			}
		}
		double infCargo = _getInfinity<double>();
		//clear(myDist); appendValue(myDist, infCargo); appendValue(myDist, infCargo); appendValue(myDist, 84);
		TStringIter strIt = begin(myDist);
		for(TSize row = 0; row < n; ++row) addVertex(distGraph);
		for(TSize row = 0; row < n; ++row) {
			for(TSize col = n - 1; col > row; --col) {
				addEdge(distGraph, row, col, value(strIt));
				value(distMatrix, row * n + col) = value(strIt);
				goNext(strIt);
			}
		}
		//removeEdge(distGraph, 0, 1);removeEdge(distGraph, 0, 2);
		for(TSize row = 0; row < n; ++row) {
			for(TSize col = n - 1; col > row; --col) {
				if (mtRand() % 100 < 50) {
					value(distMatrix, row * n + col) = infCargo;
					removeEdge(distGraph, row, col);
				}
			}
		}
		// Guide Tree
		Graph<Tree<double> > guideTreeGraph;
		upgmaTree(distGraph, guideTreeGraph, TTag() );
		Graph<Tree<double> > guideTreeMat;
		upgmaTree(distMatrix, guideTreeMat, TTag() );
		typedef Iterator<Graph<Tree<double> >, BfsIterator>::Type TBfsIterator;
		String<TSize> set1;
		TBfsIterator itBfs(guideTreeGraph, getRoot(guideTreeGraph));
		for(;!atEnd(itBfs);goNext(itBfs)) appendValue(set1, value(itBfs));
		String<TSize> set2;
		TBfsIterator itBfs2(guideTreeMat, getRoot(guideTreeMat));
		for(;!atEnd(itBfs2);goNext(itBfs2)) appendValue(set2, value(itBfs2));
		if (set1 != set2) {
			std::cout << "Randomized test failed:" << std::endl;
			std::cout << "Upgma Guide Trees:" << std::endl;
			std::cout << guideTreeMat << std::endl;
			std::cout << guideTreeGraph << std::endl;
			for(TSize i=0;i<n;++i) {
				for(TSize j=i+1;j<n;++j) {
					std::cout << value(distMatrix, i*n+j) << ",";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
			exit(0);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_GuideTree() {
//____________________________________________________________________________
// Neighbor Joining

	// Create a distance matrix
	String<double> mat;
	fill(mat, 8*8, 0);
	assignValue(mat, 0*8+1, 7);assignValue(mat, 0*8+2, 8);assignValue(mat, 0*8+3, 11);assignValue(mat, 0*8+4, 13);assignValue(mat, 0*8+5, 16);assignValue(mat, 0*8+6, 13);assignValue(mat, 0*8+7, 17);
	assignValue(mat, 1*8+2, 5);assignValue(mat, 1*8+3, 8);assignValue(mat, 1*8+4, 10);assignValue(mat, 1*8+5, 13);assignValue(mat, 1*8+6, 10);assignValue(mat, 1*8+7, 14);
	assignValue(mat, 2*8+3, 5);assignValue(mat, 2*8+4, 7);assignValue(mat, 2*8+5, 10);assignValue(mat, 2*8+6, 7);assignValue(mat, 2*8+7, 11);
	assignValue(mat, 3*8+4, 8);assignValue(mat, 3*8+5, 11);assignValue(mat, 3*8+6, 8);assignValue(mat, 3*8+7, 12);
	assignValue(mat, 4*8+5, 5);assignValue(mat, 4*8+6, 6);assignValue(mat, 4*8+7, 10);
	assignValue(mat, 5*8+6, 9);assignValue(mat, 5*8+7, 13);
	assignValue(mat, 6*8+7, 8);
	
	typedef Graph<Tree<double> > TGraph;
	TGraph guideTreeOut;
	slowNjTree(mat, guideTreeOut);
	//std::cout << guideTreeOut << std::endl;

	SEQAN_TASSERT(numVertices(guideTreeOut) == 15)
	SEQAN_TASSERT(findEdge(guideTreeOut, 8, 1) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 8, 0) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 9, 5) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 9, 4) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 10, 2) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 10, 8) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 10, 2) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 10, 8) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 11, 3) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 11, 10) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 12, 9) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 12, 11) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 13, 12) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 13, 6) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 14, 13) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 14, 7) != 0)
	SEQAN_TASSERT(getRoot(guideTreeOut) == 14)

//____________________________________________________________________________
// UPGMA
	Test_UpgmaGuideTree<UpgmaAvg>();
	Test_UpgmaGuideTree<UpgmaMin>();
	Test_UpgmaGuideTree<UpgmaMax>();
}


//////////////////////////////////////////////////////////////////////////////

void Test_Distances() {
	typedef String<AminoAcid> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	
	TString str1 = "GARFIELDTHELASTFATCAT";
	TString str2 = "GARFIELDTHEFASTCAT";
	TString str3 = "GARFIELDTHEVERYFASTCAT";
	TString str4 = "THEFATCAT";
	TStringSet strSet;
	assignValueById(strSet, str1);
	assignValueById(strSet, str2);
	assignValueById(strSet, str3);
	assignValueById(strSet, str4);
	TGraph g(strSet);

	String<double> distanceMatrix;
	getDistanceMatrix(g,distanceMatrix);
	SEQAN_TASSERT((unsigned int) getValue(distanceMatrix, 3) == (unsigned int) ((double) (1.0 - 5.0 / 7.0) * 100.0))
	SEQAN_TASSERT((unsigned int) getValue(distanceMatrix, 1 * length(strSet) + 3) == (unsigned int) ((double) (1.0 - 5.0 / 7.0) * 100.0))
	SEQAN_TASSERT((unsigned int) getValue(distanceMatrix, 2 * length(strSet) + 3) == (unsigned int) ((double) (1.0 - 3.0 / 7.0) * 100.0))

	clear(distanceMatrix);
	String<Pair<unsigned int, unsigned int> > pList;
	selectPairs(strSet, pList);
	Blosum62 score_type(-1,-11);
	String<Fragment<> > matches;
	String<int> scores;
	String<double> dist;
	appendSegmentMatches(strSet, pList, score_type, matches, scores, dist, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	getDistanceMatrix(g,distanceMatrix,LibraryDistance());
	SEQAN_TASSERT(getValue(distanceMatrix, 0 * length(strSet) + 1) < getValue(distanceMatrix, 2 * length(strSet) + 3))
	SEQAN_TASSERT(getValue(distanceMatrix, 1 * length(strSet) + 2) < getValue(distanceMatrix, 2 * length(strSet) + 3))
}


//////////////////////////////////////////////////////////////////////////////

void 
__testquickAlign(Graph<Alignment<StringSet<String<AminoAcid>, Dependent<> >, unsigned int> >& g) 
{
	Graph<Alignment<StringSet<String<AminoAcid>, Dependent<> >, void, WithoutEdgeId> > gOut(stringSet(g));
	tripletLibraryExtension(g);
	String<double> distForGuideTree;
	getDistanceMatrix(g,distForGuideTree,LibraryDistance());
	Graph<Tree<double> > guideTree;
	slowNjTree(distForGuideTree, guideTree);
	progressiveAlignment(g, guideTree, gOut);
	//std::cout << gOut << std::endl;
	String<char> alignMat;	
	convertAlignment(gOut,alignMat);
	unsigned int len = length(alignMat) / 4;
	SEQAN_TASSERT(String<char>(infix(alignMat, 0, 8)) == "GARFIELD");
	SEQAN_TASSERT(String<char>(infix(alignMat, 1*len + 0, 1*len+8)) == "GARFIELD");
	SEQAN_TASSERT(String<char>(infix(alignMat, 2*len + 0, 2*len+8)) == "GARFIELD");
	SEQAN_TASSERT(String<char>(infix(alignMat, 3*len + 0, 3*len+8)) == "--------");

	//std::cout << gOut << std::endl;
}

void Test_Libraries() {
	typedef String<AminoAcid> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	
	TString str1 = "GARFIELDTHELASTFATCAT";
	TString str2 = "GARFIELDTHEFASTCAT";
	TString str3 = "GARFIELDTHEVERYFASTCAT";
	TString str4 = "THEFATCAT";
	TStringSet strSet;
	assignValueById(strSet, str1);
	assignValueById(strSet, str2);
	assignValueById(strSet, str3);
	assignValueById(strSet, str4);
	TGraph g(strSet);
	Blosum62 score_type(-1,-11);
	String<Fragment<> > matches;
	String<int> scores;
	appendSegmentMatches(strSet, matches, scores, Lcs_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	__testquickAlign(g);
	clear(matches);
	clear(scores);
	appendSegmentMatches(strSet, matches, scores, Kmer_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	__testquickAlign(g);
	clear(matches);
	clear(scores);
	appendSegmentMatches(strSet, score_type, matches, scores, LocalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	__testquickAlign(g);
	clear(matches);
	clear(scores);
	appendSegmentMatches(strSet, score_type, matches, scores, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	__testquickAlign(g);
	clear(matches);
	clear(scores);
	String<double> distanceMatrix;
	String<Pair<unsigned int, unsigned int> > pList;
	selectPairs(strSet, pList);
	appendSegmentMatches(strSet, pList, score_type, matches, scores, distanceMatrix, AlignConfig<false,false,false,false>(), GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	__testquickAlign(g);
	Graph<Undirected<double> > distGraph;
	clear(matches);
	clear(scores);
	appendSegmentMatches(strSet, pList, score_type, matches, scores, distanceMatrix, AlignConfig<false,false,false,false>(), GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	__testquickAlign(g);
}

void Test_ExternalLibraries() {
	typedef String<char> TName;
	typedef StringSet<TName, Owner<> > TNameSet;
	typedef String<AminoAcid> TSequence;
	typedef StringSet<TSequence, Owner<> > TSequenceSet;
	
	TSequenceSet seqSet;
	appendValue(seqSet, "GARFIELDTHELASTFATCAT");
	appendValue(seqSet, "GARFIELDTHEFASTCAT");
	appendValue(seqSet, "GARFIELDTHEVERYFASTCAT");
	appendValue(seqSet, "THEFATCAT");
	TNameSet nameSet;
	appendValue(nameSet, "seq1");
	appendValue(nameSet, "seq2");
	appendValue(nameSet, "seq3");
	appendValue(nameSet, "seq4");

	typedef Graph<Alignment<StringSet<TSequence, Dependent<> >, unsigned int> > TGraph;
	TGraph g(seqSet);
	Blosum62 score_type(-1,-11);
	String<Fragment<> > matches;
	String<int> scores;
	appendSegmentMatches(seqSet, score_type, matches, scores, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	std::cout << g << std::endl;

	// T-Coffee lib format

	// Writing
	std::fstream strm;
	strm.open(TEST_PATH "test.lib", std::ios_base::out | std::ios_base::trunc);
	write(strm, g, nameSet, TCoffeeLib());
	strm.close();
	//_debugRefinedMatches(g);
	__testquickAlign(g);

	// Reading
	clear(seqSet);
	clear(nameSet);
	std::fstream strmRead;
	strmRead.open(TEST_PATH "test.lib", std::ios_base::in | std::ios_base::binary);
	read(strmRead,seqSet,nameSet,TCoffeeLib());
	strmRead.close();
	std::cout << g << std::endl;


	// Blast format

	// Writing
	std::fstream strm2;
	strm2.open(TEST_PATH "test.blast", std::ios_base::out | std::ios_base::trunc);
	write(strm2, g, nameSet, BlastLib());
	strm2.close();

	// Reading
	std::fstream strmRead2;
	strmRead2.open(TEST_PATH "test.blast", std::ios_base::in | std::ios_base::binary);
	clearVertices(g);
	clear(matches);
	clear(scores);
	read(strmRead2, matches, scores, nameSet,BlastLib());
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	strmRead2.close();
	std::cout << g << std::endl;


	// Re-read everything
	clear(g);
	assignStringSet(g, seqSet);
	std::fstream strm_lib;
	strm_lib.open(TEST_PATH "test.lib", std::ios_base::in | std::ios_base::binary);
	read(strm_lib,g,TCoffeeLib());	// Read library
	strm_lib.close();
	//_debugRefinedMatches(g);
	__testquickAlign(g);

	clear(g);
	assignStringSet(g, seqSet);
	clear(matches);
	clear(scores);
	appendSegmentMatches(seqSet, score_type, matches, scores, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	// Blast lib format

	// Writing
	std::fstream strmBlast;
	strmBlast.open(TEST_PATH "test.lib", std::ios_base::out | std::ios_base::trunc);
	write(strmBlast, g, nameSet, BlastLib());
	strmBlast.close();
	//_debugRefinedMatches(g);
	__testquickAlign(g);
	
	// Reading
	clear(g);
	assignStringSet(g, seqSet);
	std::fstream strmBlastLib;
	strmBlastLib.open(TEST_PATH "test.lib", std::ios_base::in | std::ios_base::binary);
	clear(matches);
	clear(scores);
	read(strmBlastLib, matches, scores, nameSet,BlastLib());
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	strmBlastLib.close();
	//_debugRefinedMatches(g);
	__testquickAlign(g);

}

void Test_GraphCombination() {
	typedef String<char> TName;
	typedef StringSet<TName, Owner<> > TNameSet;
	typedef String<AminoAcid> TSequence;
	typedef StringSet<TSequence, Owner<> > TSequenceSet;
	typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;
	
	TSequenceSet seqSet;
	appendValue(seqSet, "GARFIELD");
	appendValue(seqSet, "GARFIELDIELD");
	TNameSet nameSet;
	appendValue(nameSet, "seq1");
	appendValue(nameSet, "seq2");
	
	typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;
	TGraph lib1(seqSet);
	addEdge(lib1, addVertex(lib1, 0, 0, 8), addVertex(lib1, 1, 0, 8), 40);
	addVertex(lib1, 1, 8, 4);
	//_debugRefinedMatches(lib1);
	TGraph lib2(seqSet);
	addEdge(lib2, addVertex(lib2, 0, 4, 4), addVertex(lib2, 1, 4, 4), 20);
	addEdge(lib2, findVertex(lib2, 0, 4), addVertex(lib2, 1, 8, 4), 40);
	addVertex(lib2, 0, 0, 4);
	addVertex(lib2, 1, 0, 4);
	//_debugRefinedMatches(lib2);
	TGraph lib3(seqSet);
	addEdge(lib3, addVertex(lib3, 0, 2, 6), addVertex(lib3, 1, 2, 6), 30);
	addVertex(lib3, 0, 0, 2);
	addVertex(lib3, 1, 0, 2);
	addVertex(lib3, 1, 8, 4);
	//_debugRefinedMatches(lib3);
	TGraph g(seqSet);
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	appendValue(libs, &lib3);
	combineGraphs(g, libs);
	//_debugRefinedMatches(g);

	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 0), findVertex(g, 1, 0))) < cargo(findEdge(g, findVertex(g, 0, 2), findVertex(g, 1, 2))));
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 2), findVertex(g, 1, 2))) < cargo(findEdge(g, findVertex(g, 0, 4), findVertex(g, 1, 4))));
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 2), findVertex(g, 1, 2))) < cargo(findEdge(g, findVertex(g, 0, 4), findVertex(g, 1, 8))));
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 4), findVertex(g, 1, 8))) < cargo(findEdge(g, findVertex(g, 0, 4), findVertex(g, 1, 4))));

	clearVertices(g);
	combineGraphs(g, libs, FrequencyCounting() );
	//_debugRefinedMatches(g);
	
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 0), findVertex(g, 1, 0))) == 1);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 2), findVertex(g, 1, 2))) == 2);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 4), findVertex(g, 1, 4))) == 3);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 4), findVertex(g, 1, 8))) == 1);
}


void Test_TripletExtension() {
	typedef String<char> TName;
	typedef StringSet<TName, Owner<> > TNameSet;
	typedef String<AminoAcid> TSequence;
	typedef StringSet<TSequence, Owner<> > TSequenceSet;
	typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;

	TSequenceSet seqSet;//1234567890
	appendValue(seqSet, "GARFIELDTHE");
	appendValue(seqSet, "GARFIELDTHE");
	appendValue(seqSet, "GARFIELDTHE");
	appendValue(seqSet, "THE");
	TNameSet nameSet;
	appendValue(nameSet, "seq1");
	appendValue(nameSet, "seq2");
	appendValue(nameSet, "seq3");
	appendValue(nameSet, "seq4");

	// Triplet extension
	typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;
	TGraph g(seqSet);
	addEdge(g, addVertex(g, 0, 0, 8), addVertex(g, 1, 0, 8), 40);
	addEdge(g, findVertex(g, 0, 0), addVertex(g, 2, 0, 8), 30);
	addEdge(g, addVertex(g, 1, 8, 3), addVertex(g, 0, 8, 3), 40);
	addEdge(g, findVertex(g, 1, 8), addVertex(g, 2, 8, 3), 30);
	addEdge(g, findVertex(g, 1, 8), addVertex(g, 3, 0, 3), 20);
	addEdge(g, findVertex(g, 2, 8), findVertex(g, 3, 0), 40);
	tripletLibraryExtension(g);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 1, 0), findVertex(g, 2, 0))) == 30);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 8), findVertex(g, 3, 0))) == 20);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 0, 8), findVertex(g, 2, 8))) == 30);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 2, 8), findVertex(g, 3, 0))) == 60);

	// Triplet extension confined on a set of sequences
	clearVertices(g);
	addEdge(g, addVertex(g, 0, 0, 8), addVertex(g, 1, 0, 8), 40);
	addEdge(g, findVertex(g, 0, 0), addVertex(g, 2, 0, 8), 30);
	addEdge(g, addVertex(g, 1, 8, 3), addVertex(g, 0, 8, 3), 40);
	addEdge(g, findVertex(g, 1, 8), addVertex(g, 2, 8, 3), 30);
	addEdge(g, findVertex(g, 1, 8), addVertex(g, 3, 0, 3), 20);
	addEdge(g, findVertex(g, 2, 8), findVertex(g, 3, 0), 40);
	std::set<VertexDescriptor<TGraph>::Type> seqs;
	seqs.insert(1);seqs.insert(2);seqs.insert(3);
	tripletLibraryExtension(g, seqs);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 1, 0), findVertex(g, 2, 0))) == 30);
	SEQAN_TASSERT(findEdge(g, findVertex(g, 0, 8), findVertex(g, 3, 0)) == 0);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 2, 8), findVertex(g, 3, 0))) == 60);

	// Triplet extension between 2 sets of sequences
	clearVertices(g);
	addEdge(g, addVertex(g, 0, 0, 8), addVertex(g, 1, 0, 8), 40);
	addEdge(g, findVertex(g, 0, 0), addVertex(g, 2, 0, 8), 30);
	addEdge(g, addVertex(g, 1, 8, 3), addVertex(g, 0, 8, 3), 40);
	addEdge(g, findVertex(g, 1, 8), addVertex(g, 2, 8, 3), 30);
	addEdge(g, findVertex(g, 1, 8), addVertex(g, 3, 0, 3), 20);
	addEdge(g, findVertex(g, 2, 8), findVertex(g, 3, 0), 40);
	std::set<VertexDescriptor<TGraph>::Type> seqs1;
	std::set<VertexDescriptor<TGraph>::Type> seqs2;
	seqs1.insert(0);seqs1.insert(1);seqs1.insert(3);
	seqs2.insert(2);
	//_debugRefinedMatches(g);
	tripletLibraryExtension(g, seqs1, seqs2);
	//_debugRefinedMatches(g);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 1, 0), findVertex(g, 2, 0))) == 30);
	SEQAN_TASSERT(findEdge(g, findVertex(g, 0, 8), findVertex(g, 3, 0)) == 0);
	SEQAN_TASSERT(cargo(findEdge(g, findVertex(g, 1, 8), findVertex(g, 3, 0))) == 20);
}

void Test_SumOfPairsScore() {
	typedef String<char> TName;
	typedef StringSet<TName, Owner<> > TNameSet;
	typedef String<Dna> TSequence;
	typedef StringSet<TSequence, Owner<> > TSequenceSet;
	typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;
	typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;

	TSequenceSet seqSet;
	appendValue(seqSet, "ACAAGTA");
	appendValue(seqSet, "AA");
	appendValue(seqSet, "ACCTA");
	TNameSet nameSet;
	appendValue(nameSet, "seq1");
	appendValue(nameSet, "seq2");
	appendValue(nameSet, "seq3");

	TGraph g(seqSet);
	Score<int> score_type = Score<int>(5,-4,-2,-10);
	String<double> distanceMatrix;
	String<Fragment<> > matches;
	String<int> scores;
	appendSegmentMatches(seqSet, score_type, matches, scores, distanceMatrix, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	tripletLibraryExtension(g);
	Graph<Tree<double> > guideTree;
	slowNjTree(distanceMatrix, guideTree);
	TGraph gOut(seqSet);
	progressiveAlignment(g, guideTree, gOut);
	SEQAN_TASSERT(sumOfPairsScore(gOut, score_type) == -8)
	SEQAN_TASSERT(sumOfPairsScoreInd(gOut, score_type) == 16)

	seqSet[1] = "AAG";
	Score<int> scType = Score<int>(5,-4,-1,-2);
	clear(distanceMatrix);
	clearVertices(g);
	clear(matches);
	clear(scores);
	appendSegmentMatches(seqSet, scType, matches, scores, distanceMatrix, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	tripletLibraryExtension(g);
	clear(guideTree);
	slowNjTree(distanceMatrix, guideTree);
	clearVertices(gOut);
	progressiveAlignment(g, guideTree, gOut);
	SEQAN_TASSERT(sumOfPairsScore(gOut, scType) == 20)

	resize(seqSet, 2);
	seqSet[0] = "TTT";
	seqSet[1] = "AAA";
	clear(distanceMatrix);
	clear(g);
	assignStringSet(g, seqSet);
	clear(matches);
	clear(scores);
	appendSegmentMatches(seqSet, scType, matches, scores, distanceMatrix, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	tripletLibraryExtension(g);
	clear(guideTree);
	slowNjTree(distanceMatrix, guideTree);
	clear(gOut);
	assignStringSet(gOut, seqSet);
	progressiveAlignment(g, guideTree, gOut);
	SEQAN_TASSERT(sumOfPairsScore(gOut, scType) == -8)

	seqSet[0] = "TTTAAATTT";
	seqSet[1] = "AAA";
	clear(distanceMatrix);
	clearVertices(g);
	clear(matches);
	clear(scores);
	appendSegmentMatches(seqSet, scType, matches, scores, distanceMatrix, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	tripletLibraryExtension(g);
	clear(guideTree);
	slowNjTree(distanceMatrix, guideTree);
	clearVertices(gOut);
	progressiveAlignment(g, guideTree, gOut);
	SEQAN_TASSERT(sumOfPairsScore(gOut, scType) == 7)

	seqSet[0] = "AAAAAA";
	seqSet[1] = "TTTAAATTTAAATTT";
	clear(distanceMatrix);
	clearVertices(g);
	clear(matches);
	clear(scores);
	appendSegmentMatches(seqSet, scType, matches, scores, distanceMatrix, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	tripletLibraryExtension(g);
	clear(guideTree);
	slowNjTree(distanceMatrix, guideTree);
	clearVertices(gOut);
	progressiveAlignment(g, guideTree, gOut);
	SEQAN_TASSERT(sumOfPairsScore(gOut, scType) == 18)
}


void Test_Progressive() {
	typedef String<char> TName;
	typedef StringSet<TName, Owner<> > TNameSet;
	typedef String<AminoAcid> TSequence;
	typedef StringSet<TSequence, Owner<> > TSequenceSet;
	typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;
	
	TSequenceSet seqSet;
	appendValue(seqSet, "GARFIELDTHELASTFATCAT");
	appendValue(seqSet, "GARFIELDTHEFASTCAT");
	appendValue(seqSet, "GARFIELDTHEVERYFASTCAT");
	appendValue(seqSet, "THEFATCAT");
	appendValue(seqSet, "GARFIELDTHELASTFATCAT");
	appendValue(seqSet, "GARFIELDTHEFASTCAT");
	appendValue(seqSet, "GARFIELDTHEVERYFASTCAT");
	appendValue(seqSet, "THEFATCAT");
	TNameSet nameSet;
	appendValue(nameSet, "seq1");
	appendValue(nameSet, "seq2");
	appendValue(nameSet, "seq3");
	appendValue(nameSet, "seq4");
	appendValue(nameSet, "seq5");
	appendValue(nameSet, "seq6");
	appendValue(nameSet, "seq7");
	appendValue(nameSet, "seq8");

	typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;
	TGraph g(seqSet);
	Blosum62 score_type(-1,-11);
	String<double> distanceMatrix;
	String<Fragment<> > matches;
	String<int> scores;
	appendSegmentMatches(seqSet, score_type, matches, scores, distanceMatrix, GlobalPairwise_Library() );
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	tripletLibraryExtension(g);
	Graph<Tree<double> > guideTree;
	slowNjTree(distanceMatrix, guideTree);
	TGraph gOut(seqSet);
	progressiveAlignment(g, guideTree, gOut);
	int score = (int) sumOfPairsScore(gOut, score_type);
	clearVertices(gOut);
	progressiveAlignment(g, guideTree, gOut, 4);
	int score2 = (int) sumOfPairsScore(gOut, score_type);
	SEQAN_TASSERT(score == score2)
}


void Test_Alignments1() {
	typedef unsigned int TSize;
	typedef String<AminoAcid> TSequence;
	typedef StringSet<TSequence, Owner<> > TSequenceSet;
	typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;
	
	TSequenceSet seqSet;
	appendValue(seqSet, "GARFIELDTHELASTFATCAT");
	appendValue(seqSet, "GARFIELDTHEFASTCAT");
	appendValue(seqSet, "GARFIELDTHEVERYFASTCAT");
	appendValue(seqSet, "THEFATCAT");

	typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;
	TGraph gOut(seqSet);
	globalAlignment(seqSet, gOut, MSA_Protein());
	Blosum62 score_type(-2,-8);
	SEQAN_TASSERT(sumOfPairsScore(gOut, score_type) == 273)


	TSize gapExCount;
	TSize gapCount;
	TSize pairCount;
	String<TSize> numPairs;
	TSize len;
	alignmentEvaluation(gOut, score_type, gapExCount, gapCount, pairCount, numPairs, len);
	std::cout << gOut << std::endl;
	SEQAN_TASSERT(gapExCount == 33)
	SEQAN_TASSERT(gapCount == 11)
	SEQAN_TASSERT(pairCount == 83)
	SEQAN_TASSERT(len == 22)
}

void Test_Alignments2() {
	typedef String<Dna> TSequence;
	typedef StringSet<String<char>, Owner<> > TNameSet;
	typedef StringSet<TSequence, Owner<> > TSequenceSet;
	typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;
	
	TSequenceSet seqSet;
	TNameSet nameSet;
	appendValue(seqSet, "AACCTTGGAA");
	appendValue(nameSet, "seq1");
	appendValue(seqSet, "ACCTGAA");
	appendValue(nameSet, "seq2");
	appendValue(seqSet, "TTGG");
	appendValue(nameSet, "seq3");
	appendValue(seqSet, "AACCTT");
	appendValue(nameSet, "seq4");

	typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;
	TGraph gOut(seqSet);
	globalAlignment(seqSet, gOut, MSA_Dna());
	Score<int> score_type = Score<int>(5,-4,-2,-6);
	SEQAN_TASSERT(sumOfPairsScore(gOut, score_type) == 5)
	clearVertices(gOut);
	globalAlignment(seqSet, gOut, MSA_Genome());
	SEQAN_TASSERT(sumOfPairsScore(gOut, score_type) == 11)
	std::cout << gOut << std::endl;
	
	// Fasta format

	// Writing
	std::fstream strm3;
	strm3.open(TEST_PATH "align.fasta", std::ios_base::out | std::ios_base::trunc);
	write(strm3,gOut,nameSet,FastaFormat());
	strm3.close();

	// Reading
	clear(seqSet);
	clear(nameSet);
	std::fstream strmRead3a;
	strmRead3a.open(TEST_PATH "align.fasta", std::ios_base::in | std::ios_base::binary);
	read(strmRead3a,seqSet,nameSet,FastaAlign()); 
	strmRead3a.close();
	std::fstream strmRead3;
	strmRead3.open(TEST_PATH "align.fasta", std::ios_base::in | std::ios_base::binary);
	clearVertices(gOut);
	read(strmRead3, gOut, nameSet, score_type, FastaAlign());
	strmRead3.close();
	std::cout << gOut << std::endl;
}



void Test_ReversableFragments() {
	typedef unsigned int TSize;
	typedef String<Dna> TSequence;
	TSequence seq1 = "AACGTT";
	TSequence seq2 = "AACGTTC";
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet;
	appendValue(strSet, seq1);
	appendValue(strSet, seq2);
	typedef Fragment<TSize, ExactReversableFragment<> > TFragment;
	typedef String<TFragment> TFragmentString;
	TFragmentString matches;
	appendValue(matches, TFragment(0,0,1,0,2));
	appendValue(matches, TFragment(0,3,1,3,3));
	appendValue(matches, TFragment(0,1,1,2,3, true));
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	TGraph g(strSet);
	matchRefinement(matches,strSet,g);
	SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 1), findVertex(g, 1, 2)) == 0)
	SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 1), findVertex(g, 1, 4)) != 0)
	SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 2), findVertex(g, 1, 3)) != 0)
	SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 3), findVertex(g, 1, 2)) != 0)
	SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 3), findVertex(g, 1, 4)) == 0)
	SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 0), findVertex(g, 1, 0)) != 0)
	SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 1), findVertex(g, 1, 1)) != 0)
}

//////////////////////////////////////////////////////////////////////////////


void Test_GraphTCoffee() {
	Test_GuideTree();
	Test_Distances();
	Test_Libraries();
	Test_ExternalLibraries();
	Test_GraphCombination();
	Test_TripletExtension();
	Test_SumOfPairsScore();
	Test_Progressive();
	Test_Alignments1();
	Test_Alignments2();
	Test_ReversableFragments();
	
	
	debug::verifyCheckpoints("projects/library/seqan/graph_msa/graph_align_tcoffee_kmer.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_msa/graph_align_tcoffee_distance.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_msa/graph_align_tcoffee_guidetree.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_msa/graph_align_tcoffee_library.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_msa/graph_align_tcoffee_io.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_msa/graph_align_tcoffee_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_msa/graph_align_tcoffee_progressive.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_msa/graph_align_tcoffee_msa.h");
}


}

#endif

