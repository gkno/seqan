#ifndef SEQAN_HEADER_TEST_GRAPH_TCOFFEE_H
#define SEQAN_HEADER_TEST_GRAPH_TCOFFEE_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

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
	clear(mat);
	fill(mat, 8*8, 0);
	assignValue(mat, 0*8+1, 32);assignValue(mat, 0*8+2, 48);assignValue(mat, 0*8+3, 51);assignValue(mat, 0*8+4, 50);assignValue(mat, 0*8+5, 48);assignValue(mat, 0*8+6, 98);assignValue(mat, 0*8+7, 148);
	assignValue(mat, 1*8+2, 26);assignValue(mat, 1*8+3, 34);assignValue(mat, 1*8+4, 29);assignValue(mat, 1*8+5, 33);assignValue(mat, 1*8+6, 84);assignValue(mat, 1*8+7, 136);
	assignValue(mat, 2*8+3, 42);assignValue(mat, 2*8+4, 44);assignValue(mat, 2*8+5, 44);assignValue(mat, 2*8+6, 92);assignValue(mat, 2*8+7, 152);
	assignValue(mat, 3*8+4, 44);assignValue(mat, 3*8+5, 38);assignValue(mat, 3*8+6, 86);assignValue(mat, 3*8+7, 142);
	assignValue(mat, 4*8+5, 24);assignValue(mat, 4*8+6, 89);assignValue(mat, 4*8+7, 142);
	assignValue(mat, 5*8+6, 90);assignValue(mat, 5*8+7, 142);
	assignValue(mat, 6*8+7, 148);
	
	clear(guideTreeOut);
	upgmaTree(mat, guideTreeOut);

	//std::cout << guideTreeOut << std::endl;
	SEQAN_TASSERT(numVertices(guideTreeOut) == 15)
	SEQAN_TASSERT(findEdge(guideTreeOut, 8, 4) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 8, 5) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 9, 1) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 9, 2) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 10, 8) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 10, 9) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 11, 3) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 11, 10) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 12, 11) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 12, 0) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 13, 6) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 13, 12) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 14, 7) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 14, 13) != 0)
	SEQAN_TASSERT(getRoot(guideTreeOut) == 14)

//____________________________________________________________________________
// UpgmaMin vs. UpgmaAvg
	clear(mat);
	fill(mat, 4*4, 0);
	assignValue(mat, 0*4+1, 30);assignValue(mat, 0*4+2, 80);assignValue(mat, 0*4+3, 60);
	assignValue(mat, 1*4+2, 20);assignValue(mat, 1*4+3, 40);
	assignValue(mat, 2*4+3, 40);
	clear(guideTreeOut);
	upgmaTree(mat, guideTreeOut, UpgmaAvg());
	//std::cout << guideTreeOut << std::endl;
	SEQAN_TASSERT(findEdge(guideTreeOut, 1, 4) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 2, 4) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 3, 5) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 0, 6) != 0)
	clear(mat);
	fill(mat, 4*4, 0);
	assignValue(mat, 0*4+1, 30);assignValue(mat, 0*4+2, 80);assignValue(mat, 0*4+3, 60);
	assignValue(mat, 1*4+2, 20);assignValue(mat, 1*4+3, 40);
	assignValue(mat, 2*4+3, 40);
	clear(guideTreeOut);
	upgmaTree(mat, guideTreeOut, UpgmaMin());
	//std::cout << guideTreeOut << std::endl;
	SEQAN_TASSERT(findEdge(guideTreeOut, 1, 4) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 2, 4) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 0, 5) != 0)
	SEQAN_TASSERT(findEdge(guideTreeOut, 3, 6) != 0)
}

//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeGarfield() {
//____________________________________________________________________________
// T-Coffee
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

	// Score object
	Blosum62 score_type_global(-1,-11);
	Blosum62 score_type_local(-2,-8);
	
	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(strSet);
	generatePrimaryLibrary(lib1, score_type_local, GlobalPairwise_Library() );
	//generatePrimaryLibrary(lib1, score_type_local, MUMPairwise_Library() );
	//generatePrimaryLibrary(lib1, score_type_global, 5, false, Kmer_Library() );
	
	fstream strm01; // Alignment graph as dot
	strm01.open(TEST_PATH "my_tcoffee01.dot", ios_base::out | ios_base::trunc);
	write(strm01,lib1,DotDrawing());
	strm01.close();

	// Generate a primary library, i.e., all local pairwise alignments
	TGraph lib2(strSet);
	generatePrimaryLibrary(lib2, score_type_local, LocalPairwise_Library() );

	fstream strm02; // Alignment graph as dot
	strm02.open(TEST_PATH "my_tcoffee02.dot", ios_base::out | ios_base::trunc);
	write(strm02,lib2,DotDrawing());
	strm02.close();

	// Weighting of libraries (Signal addition)
	TGraph g(strSet);
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	combineGraphs(g, libs);

	// Clear the old libraries
	clear(lib1);
	clear(lib2);

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

	// Pairwise Distances
	String<double> distanceMatrix;
	getDistanceMatrix(g, distanceMatrix, LibraryDistance() ); 
	//getDistanceMatrix(g, distanceMatrix, KmerDistance() );
	//getDistanceMatrix(g, distanceMatrix, 6, AAGroupsDayhoff(), KmerDistance());

	// Build the neighbor joining tree
	Graph<Tree<double> > guideTree;
	upgmaTree(distanceMatrix, guideTree, UpgmaMin() );

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void> > gOut(strSet);
	progressiveAlignment(g, guideTree, gOut);
	//iterativeProgressiveAlignment(g, guideTree, score_type_local, gOut);
	
	// Print the alignment
	std::cout << gOut << std::endl;

	fstream strm3; // Alignment graph as dot
	strm3.open(TEST_PATH "my_tcoffee3.dot", ios_base::out | ios_base::trunc);
	write(strm3,gOut,DotDrawing());
	strm3.close();

	String<String<char> > names;
	appendValue(names, "seq1");	appendValue(names, "seq2");
	appendValue(names, "seq3");	appendValue(names, "seq4");
	fstream strm4; // Alignment graph as msf
	strm4.open(TEST_PATH "my_alignment.fasta", ios_base::out | ios_base::trunc);
	write(strm4,gOut,names, FastaFormat());
	strm4.close();
}

//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeFromFile(String<char> const in_path, String<char> const file_prefix, String<char> const file_suffix) {
//____________________________________________________________________________
// T-Coffee
	// Timing variables
	clock_t bigbang, startTime;

	// Import Sequences
	typedef String<AminoAcid> TString;
	StringSet<String<AminoAcid>, Owner<> > origStrSet;
	String<String<char> > names;
	startTime = clock();
	bigbang = startTime;
	unsigned int nucCount = _alignImportSequences(in_path, file_prefix, file_suffix, origStrSet, names);
	unsigned int seqCount = length(origStrSet);
	if (length(file_suffix)) std::cout << in_path << file_prefix << '.' << file_suffix << std::endl;
	else std::cout << in_path << file_prefix << std::endl;
	std::cout << "Total number of bp: " << nucCount << ", Number of sequences: " << seqCount << std::endl;
	_alignTiming(startTime, "Import sequences done: ");

	// Make dependent string set
	typedef StringSet<TString, Dependent<> > TStringSet;
	TStringSet strSet;
	for(unsigned int i = 0; i<seqCount; ++i) appendValue(strSet, origStrSet[i]);

	// Align the sequences
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(strSet);
	tCoffeeProteinAlignment(strSet, gOut);
	_alignTiming(startTime, "T-Coffee alignment done: ");
	
	// Output alignment
	std::stringstream output2;
	if (length(file_suffix)) output2 << in_path << file_prefix << "." << file_suffix << "_my";
	else output2 << in_path << file_prefix << "." << "my";
	fstream strm5; // Alignment graph as fasta
	strm5.open(output2.str().c_str(), ios_base::out | ios_base::trunc);
	write(strm5,gOut,names, MsfFormat());
	strm5.close();
	_alignTiming(startTime, "Alignment output done: ");

	// Finished
	clear(gOut);
	_alignTiming(bigbang, "Total time: ");
	std::cout << "==============================" << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeAssemblyFromFile(String<char> const readsPath, String<char> const sequencePath) {
//____________________________________________________________________________
// T-Coffee
	// Timing variables
	clock_t bigbang, startTime;

	// Import Sequences
	typedef String<Dna> TString;
	StringSet<TString, Owner<> > origStrSet;
	String<String<char> > names;
	startTime = clock();
	bigbang = startTime;
	unsigned int nucCount = _alignImportSequences(readsPath, "", "", origStrSet, names);
	unsigned int seqCount = length(origStrSet);
	std::cout << readsPath << std::endl;
	std::cout << "Total number of bp: " << nucCount << ", Number of sequences: " << seqCount << std::endl;
	_alignTiming(startTime, "Import sequences done: ");

	// Make dependent string set
	typedef StringSet<TString, Dependent<> > TStringSet;
	TStringSet strSet;
	for(unsigned int i = 0; i<seqCount; ++i) appendValue(strSet, origStrSet[i]);

	// Score objects
	Score<int> score_type_global = Score<int>(5,-4,-6,-14);

	// Align the sequences
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef Size<TGraph>::Type TSize;
	typedef Id<TGraph>::Type TId;

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph g(strSet);
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, names, pList);
	String<double> distanceMatrix;
	generatePrimaryLibrary(g, pList, score_type_global, distanceMatrix, Overlap_Library() );
	_alignTiming(startTime, "Overlap done: ");
	
	// Triplet library extension
	tripletLibraryExtension(g);
	std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;
	_alignTiming(startTime, "Triplet done: ");

	// Build the guide tree
	Graph<Tree<double> > guideTree;
	upgmaTree(distanceMatrix, guideTree);
	_alignTiming(startTime, "Guide tree done: ");

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void> > gOut(strSet);
	progressiveAlignment(g, guideTree, gOut);
	clear(guideTree);
	clear(g);
	_alignTiming(startTime, "Progressive alignment done: ");

	// Create the consensus alignment
	std::cout << "Consensus sequence: " << std::endl;
	TString consens;
	consensusSequence(gOut, consens);
	std::cout << consens << std::endl;
	_alignTiming(startTime, "Consensus done: ");
	

	exit(0);


	// Make quick alignment
	std::cout << "Alignment of consensus to reference: " << std::endl;
	StringSet<TString, Owner<> > seqSet;
	_alignImportSequences(sequencePath, "", "", seqSet, names);
	TStringSet str;
	assignValueById(str, consens);
	assignValueById(str, seqSet[0]);
	TGraph ali(str);
	Score<int> score_type = Score<int>(0,-1,-1,0);
	globalAlignment(ali, score_type, NeedlemanWunsch() );
	std::cout << ali << std::endl;

	// Output alignment
	std::stringstream output2;
	output2 << sequencePath << "." << "my";
	fstream strm5; // Alignment graph as fasta
	strm5.open(output2.str().c_str(), ios_base::out | ios_base::trunc);
	write(strm5,gOut,names, FastaReadFormat());
	strm5.close();
	_alignTiming(startTime, "Alignment output done: ");

	// Finished
	clear(gOut);
	_alignTiming(bigbang, "Total time: ");
	std::cout << "==============================" << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeFromLibrary(String<char> const in_path) {
//____________________________________________________________________________
// Graph TCoffee
	clock_t bigbang, startTime;
	startTime = clock();
	bigbang = startTime;

	// Read a t-coffee library
	typedef String<Rna5> TString;
	fstream strm;
	std::stringstream input;
	input << in_path;
	std::cout << input.str().c_str() << std::endl;
	strm.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
	String<String<char> > names;
	typedef StringSet<TString, Owner<> > TOwnerStringSet;
	TOwnerStringSet oriStr;
	read(strm,oriStr,names,TCoffeeLib());	// Read identifiers and strings
	strm.close();
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	TStringSet strSet;
	for(unsigned int i = 0; i<length(oriStr); ++i) {
		assignValueById(strSet, oriStr, i);
	}
	TGraph lib1(strSet);
	fstream strm2;
	strm2.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
	read(strm2,lib1,TCoffeeLib());	// Read sequences
	strm2.close();
	_alignTiming(startTime, "Library reading done: ");

	//// Write the library
	//fstream strmW; 
	//strmW.open(TEST_PATH "my_garfield.lib", ios_base::out | ios_base::trunc);
	//write(strmW,g,TCoffeeLib());
	//strmW.close();

	//Score<int> score_type = Score<int>(5,-4,-1,-10);
	Score<int> score_type = Score<int>(5,-4,-4,-14);
	
	//// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib2(strSet);
	String<double> distanceMatrix;
	generatePrimaryLibrary(lib2, distanceMatrix, score_type, GlobalPairwise_Library() );
	_alignTiming(startTime, "Global Pairwise alignments done: ");
	std::cout << "Library size: " << numVertices(lib2) << " Vertices, " << numEdges(lib2) << " Edges" << std::endl;

	//// Generate a primary library, i.e., all local pairwise alignments
	//TGraph lib3(strSet);
	//generatePrimaryLibrary(lib3, score_type, LocalPairwise_Library() );
	//_alignTiming(startTime, "Local Pairwise alignments done: "); 
	//std::cout << "Library size: " << numVertices(lib3) << " Vertices, " << numEdges(lib3) << " Edges" << std::endl;

	//// Weighting of libraries (Signal addition)
	TGraph g(strSet);
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	combineGraphs(g, libs);
	_alignTiming(startTime, "Combining graphs / libraries done: "); 
	std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;

	// Clear the old libraries
	clear(lib1);
	clear(lib2);
	
	// Triplet library extension
	tripletLibraryExtension(g);
	_alignTiming(startTime, "Triplet done: ");
	std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;

	// Build the guide tree
	Graph<Tree<double> > guideTree;
	upgmaTree(distanceMatrix, guideTree);
	_alignTiming(startTime, "Guide tree done: ");

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(strSet);
	//progressiveAlignment(g, guideTree, gOut, Gotoh() );
	progressiveAlignment(g, guideTree, gOut);
	_alignTiming(startTime, "Progressive alignment done: ");
	std::cout << gOut << std::endl;
	clear(guideTree);

	// Output alignment
	std::stringstream output2;
	output2 << in_path << "." << "fasta";
	fstream strm6; // Alignment graph as msf
	strm6.open(output2.str().c_str(), ios_base::out | ios_base::trunc);
	write(strm6,gOut,names, FastaFormat());
	strm6.close();
	_alignTiming(startTime, "Alignment output done: ");

	// Clean-up
	clear(distanceMatrix);
	clear(gOut);
	clear(g);
	_alignTiming(startTime, "Clean-up done: ");
	_alignTiming(bigbang, "Total time: ");
	std::cout << "==============================" << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeScratch(String<char> const in_path) {
//____________________________________________________________________________
// T-Coffee
	clock_t bigbang, startTime;
	startTime = clock();
	bigbang = startTime;

	// Import Sequences
	typedef String<Dna> TString;
	StringSet<String<Dna>, Owner<> > origStrSet;
	String<String<char> > names;
	unsigned int nucCount = _alignImportSequences(in_path, "", "", origStrSet, names);
	unsigned int seqCount = length(origStrSet);
	_alignTiming(startTime, "Library reading done: ");
	std::cout << "Total number of bp: " << nucCount << ", Number of sequences: " << seqCount << std::endl;

	// Align the sequences
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	typedef Size<TGraph>::Type TSize;
	typedef Id<TGraph>::Type TId;

	// Make dependent string set
	TStringSet strSet;
	for(unsigned int i = 0; i<seqCount; ++i) appendValue(strSet, origStrSet[i]);
	
	//Score<int> score_type = Score<int>(5,-4,-1,-10);
	Score<int> score_type = Score<int>(5,-4,-4,-14);

	// Generate kmer primary library
	TGraph g(strSet);
	//generatePrimaryLibrary(g, score_type, MUMPairwise_Library() );
	generatePrimaryLibrary(g, score_type, 15, false, Kmer_Library() );
	_alignTiming(startTime, "MUM library done: ");
	std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;
	
	// Triplet library extension
	tripletLibraryExtension(g);
	_alignTiming(startTime, "Triplet done: ");
	std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;

	// Pairwise Distances
	String<double> distanceMatrix;
	getDistanceMatrix(g, distanceMatrix, LibraryDistance() ); 

	// Build the guide tree
	Graph<Tree<double> > guideTree;
	upgmaTree(distanceMatrix, guideTree);
	_alignTiming(startTime, "Guide tree done: ");

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(strSet);
	progressiveAlignment(g, guideTree, gOut);
	_alignTiming(startTime, "Progressive alignment done: ");
	//std::cout << gOut << std::endl;
	clear(guideTree);

	// Output alignment
	std::stringstream output3;
	output3 << in_path << "." << "cgviz";
	fstream strm7; // Alignment graph as cgviz
	strm7.open(output3.str().c_str(), ios_base::out | ios_base::trunc);
	write(strm7,gOut,names, CgVizFormat());
	strm7.close();
	_alignTiming(startTime, "Alignment output done: ");

	// Clean-up
	clear(distanceMatrix);
	clear(gOut);
	clear(g);
	_alignTiming(startTime, "Clean-up done: ");
	_alignTiming(bigbang, "Total time: ");
	std::cout << "==============================" << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

void Test_BaliBaseRef11() {
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\Balibase\\bb3_release\\RV11\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/Balibase/bb3_release/RV11/");
#endif

	//for(int i = 1; i<39; ++i) {
	//	std::stringstream s;
	//	s << "BBS110";
	//	if (i < 10) s << '0';
	//	s << i;
	//	Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	//}

	for(int i = 1; i<39; ++i) {
		std::stringstream s;
		s << "BB110";
		if (i < 10) s << '0';
		s << i;
		Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_BaliBaseRef12() {
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\Balibase\\bb3_release\\RV12\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/Balibase/bb3_release/RV12/");
#endif

	//for(int i = 1; i<45; ++i) {
	//	std::stringstream s;
	//	s << "BBS120";
	//	if (i < 10) s << '0';
	//	s << i;
	//	Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	//}

	for(int i = 1; i<45; ++i) {
		std::stringstream s;
		s << "BB120";
		if (i < 10) s << '0';
		s << i;
		Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	}
}


//////////////////////////////////////////////////////////////////////////////

void Test_BaliBaseRef20() {
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\Balibase\\bb3_release\\RV20\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/Balibase/bb3_release/RV20/");
#endif

	//for(int i = 1; i<42; ++i) {
	//	std::stringstream s;
	//	s << "BBS200";
	//	if (i < 10) s << '0';
	//	s << i;
	//	Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	//}

	for(int i = 1; i<42; ++i) {
		std::stringstream s;
		s << "BB200";
		if (i < 10) s << '0';
		s << i;
		Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_BaliBaseRef30() {
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\Balibase\\bb3_release\\RV30\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/Balibase/bb3_release/RV30/");
#endif

	//for(int i = 1; i<31; ++i) {
	//	std::stringstream s;
	//	s << "BBS300";
	//	if (i < 10) s << '0';
	//	s << i;
	//	Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	//}

	for(int i = 1; i<31; ++i) {
		std::stringstream s;
		s << "BB300";
		if (i < 10) s << '0';
		s << i;
		Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_BaliBaseRef40() {
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\Balibase\\bb3_release\\RV40\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/Balibase/bb3_release/RV40/");
#endif

	for(int i = 1; i<50; ++i) {
		std::stringstream s;
		s << "BB400";
		if (i < 10) s << '0';
		s << i;
		Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_BaliBaseRef50() {
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\Balibase\\bb3_release\\RV50\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/Balibase/bb3_release/RV50/");
#endif

	//for(int i = 1; i<17; ++i) {
	//	std::stringstream s;
	//	s << "BBS500";
	//	if (i < 10) s << '0';
	//	s << i;
	//	Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	//}

	for(int i = 1; i<17; ++i) {
		std::stringstream s;
		s << "BB500";
		if (i < 10) s << '0';
		s << i;
		Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	}
}

//////////////////////////////////////////////////////////////////////////////

void _findFilesPrefab(char const * path, char const * index_file_name)
{
	char filepath[1024];
	strcpy(filepath, path);
	char * filename = filepath + strlen(path);

	strcpy(filename, index_file_name);

	fstream strm;
	strm.open(filepath, ios_base::in);
	while (!strm.eof())
	{
		strm.getline(filename, 512);
		//std::cout << filename << "\n";
		Test_TCoffeeFromFile(filepath,"", "");
	}
	//std::cout << "---------------" << "\n";
	strm.close();
}

//////////////////////////////////////////////////////////////////////////////

void Test_Prefab() {
	// Windows
#ifdef PLATFORM_WINDOWS
	_findFilesPrefab("Z:\\Prefab4.0\\in\\", "..\\dir.txt"); 
#else
	// Linux
	_findFilesPrefab("/home/takifugu/rausch/Prefab4.0/in/", "..//dir.txt"); 
#endif
}

//////////////////////////////////////////////////////////////////////////////

void _findFilesRna(char const * path, char const * index_file_name)
{
	char filepath[1024];
	strcpy(filepath, path);
	char * filename = filepath + strlen(path);

	strcpy(filename, index_file_name);

	fstream strm;
	strm.open(filepath, ios_base::in);
	while (!strm.eof())
	{
		strm.getline(filename, 512);
		std::cout << filename << "\n";
		Test_TCoffeeFromLibrary(filepath);
	}
	std::cout << "---------------" << "\n";
	strm.close();
}


//////////////////////////////////////////////////////////////////////////////

void Test_RnaLibraries() {
	// Windows
#ifdef PLATFORM_WINDOWS
	_findFilesRna("Z:\\Bralibase\\libs\\k3\\stacked_lara\\", "dir.txt"); 
	_findFilesRna("Z:\\Bralibase\\libs\\k5\\stacked_lara\\", "dir.txt"); 
	_findFilesRna("Z:\\Bralibase\\libs\\k7\\stacked_lara\\", "dir.txt"); 
	_findFilesRna("Z:\\Bralibase\\libs\\k10\\stacked_lara\\", "dir.txt"); 
	_findFilesRna("Z:\\Bralibase\\libs\\k15\\stacked_lara\\", "dir.txt"); 
	//_findFilesRna("Z:\\Bralibase\\libs\\k3\\lara\\", "dir.txt"); 
#else
	// Linux
	_findFilesRna("/home/takifugu/rausch/Bralibase/libs/k3/stacked_lara/", "dir.txt"); 
	_findFilesRna("/home/takifugu/rausch/Bralibase/libs/k5/stacked_lara/", "dir.txt"); 
	_findFilesRna("/home/takifugu/rausch/Bralibase/libs/k7/stacked_lara/", "dir.txt"); 
	_findFilesRna("/home/takifugu/rausch/Bralibase/libs/k10/stacked_lara/", "dir.txt"); 
	_findFilesRna("/home/takifugu/rausch/Bralibase/libs/k15/stacked_lara/", "dir.txt"); 
	//_findFilesRna("/home/takifugu/rausch/Bralibase/libs/k3/lara/", "dir.txt"); 
#endif
}

//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeAssembly() {
	// Windows
#ifdef PLATFORM_WINDOWS
	Test_TCoffeeAssemblyFromFile("D:\\matches\\reads\\readsim.out", "D:\\matches\\reads\\sequence.fasta"); 
#else
	// Linux
	Test_TCoffeeAssemblyFromFile("/home/takifugu/rausch/matches/reads/readsim.out", "/home/takifugu/rausch/matches/reads/sequence.fasta"); 
#endif
}

//////////////////////////////////////////////////////////////////////////////


void Test_GraphTCoffee() {
	Test_CompressedAlphabets();
	Test_GuideTree();
	
	// Toy example
	Test_TCoffeeGarfield();
	
	// Assembly
	//Test_TCoffeeAssembly();

	//// Balibase
	//Test_BaliBaseRef11();
	//Test_BaliBaseRef12();
	//Test_BaliBaseRef20();
	//Test_BaliBaseRef30();
	//Test_BaliBaseRef40();
	//Test_BaliBaseRef50();

	// Prefab
	//Test_Prefab();

	// RnaLibraries
	//Test_RnaLibraries();
	//Test_TCoffeeFromLibrary("/home/takifugu/rausch/Bralibase/small.lib");
	//Test_TCoffeeScratch("Z://matches//test//pox_subset.fasta");
	//Test_TCoffeeScratch("/home/takifugu/rausch/matches/test/adeno.fasta");
	//Test_TCoffeeScratch("Z://matches//test//adeno.fasta");

	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee_msa.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee_kmer.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee_distance.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee_guidetree.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee_library.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee_progressive.h");
}


}

#endif

