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

	// Score object
	Score<int> score_type = Score<int>(4,-1,-1,-10);

	// Generate a primary library, i.e., all pairwise alignments
	generatePrimaryLibrary(lib1, score_type, GlobalPairwise_Library() );
	generatePrimaryLibrary(lib2, score_type, LocalPairwise_Library() );

	// Weighting of libraries (Signal addition)
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	combineGraphs(g, libs);

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
	String<double> distanceMatrix; 
	getKmerSimilarityMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
	similarityToDistanceMatrix(distanceMatrix, KimuraDistance() );

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
	typedef Graph<Alignment<TStringSet, int, Default> > TGraph;
	
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
	Score<int, Pam<> > score_type_global(250, -1, -2);
	Score<int, Pam<> > score_type_local(250, -1, -2);

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(strSet);
	generatePrimaryLibrary(lib1, score_type_global, GlobalPairwise_Library() );

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

	// Calculate a distance matrix using a compressed alphabet or not
	String<double> distanceMatrix; 
	getKmerSimilarityMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
	similarityToDistanceMatrix(distanceMatrix, KimuraDistance() );

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

	String<String<char> > names;
	appendValue(names, "seq1");	appendValue(names, "seq2");
	appendValue(names, "seq3");	appendValue(names, "seq4");
	fstream strm4; // Alignment graph as msf
	strm4.open(TEST_PATH "my_alignment.msf", ios_base::out | ios_base::trunc);
	write(strm4,gOut,names, MsfFormat());
	strm4.close();
}

//////////////////////////////////////////////////////////////////////////////

void Test_TCoffeeFromFile(String<char> const in_path, String<char> const file_prefix, String<char> const file_suffix) {
//____________________________________________________________________________
// T-Coffee
	clock_t bigbang, startTime,endTime;
	double time;

	// Sequences
	typedef String<AminoAcid> TString;
	StringSet<String<AminoAcid>, Owner<> > origStrSet;

	// Count sequences and read names
	startTime = clock();
	bigbang = startTime;
	String<String<char> > names;
	unsigned seqCount = 0;
	ifstream file;
	std::stringstream input;
	input << in_path << file_prefix << '.' << file_suffix;
	std::cout << input.str().c_str() << std::endl;
	file.open(input.str().c_str(), ios_base::in | ios_base::binary);
	if (!file.is_open()) return;
	while (!_streamEOF(file)) {
		String<char> id;
		readID(file, id, Fasta());
		appendValue(names, id);
		//std::cout << id << std::endl;
		goNext(file, Fasta());
		++seqCount;
	}
	std::cout << "Number of sequences: " << seqCount << std::endl;

	// Import sequences
	file.clear();
	file.seekg(0, ios_base::beg);
	resize(origStrSet, seqCount);
	unsigned int count = 0;
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 	{
		read(file, origStrSet[i], Fasta());
		count += length(origStrSet[i]);
		std::cout << "Length of Sequence " << i << ": " << length(origStrSet[i]);
		if (i<seqCount-1) std::cout << std::endl;
	}
    file.close();
	std::cout << std::endl << "Total number of bp: " << count << std::endl;
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Import sequences done: " << time << " sec" << std::endl;


	// Make dependent string set
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
	TStringSet strSet;
	for(unsigned int i = 0; i<length(origStrSet); ++i) {
		appendValue(strSet, origStrSet[i]);
	}

	// Score objects
	Score<int, Pam<> > score_type_global(250, -1, -17);
	Score<int, Pam<> > score_type_local(250, -2, -13);

	// Generate a primary library, i.e., all global pairwise alignments
	TGraph lib1(strSet);
	generatePrimaryLibrary(lib1, score_type_global, GlobalPairwise_Library() );
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Stage 1)" << std::endl;
	std::cout << "Global Pairwise alignments done: " << time << " sec" << std::endl;
	std::cout << "Library size: " << numVertices(lib1) << " Vertices, " << numEdges(lib1) << " Edges" << std::endl;

	// Generate a primary library, i.e., all local pairwise alignments
	TGraph lib2(strSet);
	generatePrimaryLibrary(lib2, score_type_local, LocalPairwise_Library() );
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Local Pairwise alignments done: " << time << " sec" << std::endl;
	std::cout << "Library size: " << numVertices(lib2) << " Vertices, " << numEdges(lib2) << " Edges" << std::endl;

	// Weighting of libraries (Signal addition)
	TGraph g(strSet);
	String<TGraph*> libs;
	appendValue(libs, &lib1);
	appendValue(libs, &lib2);
	combineGraphs(g, libs);
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Combining graphs / libraries done: " << time << " sec" << std::endl;
	std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;


	// Clear the old libraries
	clear(lib1);
	clear(lib2);

	// Triplet library extension
	tripletLibraryExtension(g);
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Triplet done: " << time << " sec" << std::endl;
	std::cout << "Library size: " << numVertices(g) << " Vertices, " << numEdges(g) << " Edges" << std::endl;

	// Calculate a distance matrix using a compressed alphabet or not
	String<double> distanceMatrix; 
	getKmerSimilarityMatrix(stringSet(g), distanceMatrix, 4, AAGroupsDayhoff() );
	similarityToDistanceMatrix(distanceMatrix, KimuraDistance() );
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Kmer done: " << time << " sec" << std::endl;

	// Build the neighbor joining tree
	Graph<Tree<double> > njTreeOut;
	slowNjTree(distanceMatrix, njTreeOut);
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "NjTree done: " << time << " sec" << std::endl;

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void> > gOut(strSet);
	progressiveAlignment(g, njTreeOut, gOut, Hirschberg() );
	//std::cout << gOut << std::endl;
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Progressive alignment done: " << time << " sec" << std::endl;
	
	//fstream strm3; // Alignment graph as dot
	//strm3.open(TEST_PATH "my_tcoffee3.dot", ios_base::out | ios_base::trunc);
	//write(strm3,gOut,DotDrawing());
	//strm3.close();

	//std::stringstream output3;
	//output3 << in_path << "my_" << file_prefix << ".msf_middle";
	//fstream strm9; // Alignment graph as msf
	//strm9.open(output3.str().c_str(), ios_base::out | ios_base::trunc);
	//write(strm9,gOut,names, MsfFormat());
	//strm9.close();

	// Calculate a distance matrix using a compressed alphabet or not
	clear(distanceMatrix);
	getSequenceSimilarity(gOut, distanceMatrix, AminoAcid() );
	for(unsigned int row = 0; row<seqCount; ++row) {
		for(unsigned int col = row + 1; col<seqCount; ++col) {
			double val = getValue(distanceMatrix, row * seqCount + col);
			if (length((stringSet(g))[row]) > length((stringSet(g))[col])) val /= length((stringSet(g))[col]);
			else val /= length((stringSet(g))[row]);
			assignValue(distanceMatrix, row * seqCount + col, val);
			assignValue(distanceMatrix, col * seqCount + row, val);
		}
	}
	similarityToDistanceMatrix(distanceMatrix, KimuraDistance() );
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Stage 2)" << std::endl;
	std::cout << "Kmer done: " << time << " sec" << std::endl;

	// Build the neighbor joining tree
	clearVertices(njTreeOut);
	slowNjTree(distanceMatrix, njTreeOut);
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "NjTree done: " << time << " sec" << std::endl;

	// Perform a progressive alignment
	clearVertices(gOut);
	progressiveAlignment(g, njTreeOut, gOut, Hirschberg() );
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Progressive alignment done: " << time << " sec" << std::endl;

	// Print the alignment
	// std::cout << gOut << std::endl;

	std::stringstream output2;
	output2 << in_path << file_prefix << ".tfa_my";
	fstream strm5; // Alignment graph as msf
	strm5.open(output2.str().c_str(), ios_base::out | ios_base::trunc);
	write(strm5,gOut,names, MsfFormat());
	strm5.close();

	// Clean-up
	clear(g);
	clear(gOut);
	endTime=clock();
	time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << "Clean-up done: " << time << " sec" << std::endl;

	endTime=clock();
	time=((float)(endTime-bigbang)/CLOCKS_PER_SEC);
	std::cout << "Total time: " << time << " sec" << std::endl;
	std::cout << "==============================" << std::endl;
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
		unsigned int seqCount = ((Byte) mtRand() % 5) + 3;
		for(unsigned int i=0;i<seqCount;++i) {
			TString str;

			unsigned int stop = ((unsigned int) mtRand() % 200) + 6;
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

		// Score object
		Score<int> score_type = Score<int>(4,-1,-1,-10);

		// Generate primary libraries
		TGraph lib1(strSet);
		TGraph lib2(strSet);
		generatePrimaryLibrary(lib1, score_type, GlobalPairwise_Library() );
		std::cout << "Global Pairwise alignments done" << std::endl;

		generatePrimaryLibrary(lib2, score_type, LocalPairwise_Library());
		std::cout << "Local Pairwise alignments done" << std::endl;

		// Weighting of libraries (Signal addition)
		String<TGraph*> libs;
		appendValue(libs, &lib1);
		appendValue(libs, &lib2);
		combineGraphs(g, libs);
		std::cout << "Combining graphs done" << std::endl;

		// Triplet library extension
		tripletLibraryExtension(g);
		std::cout << "Triplet done" << std::endl;

		// Calculate a distance matrix using a compressed alphabet or not
		String<double> distanceMatrix; 
		getKmerSimilarityMatrix(stringSet(g), distanceMatrix, 6, AAGroupsDayhoff() );
		std::cout << "Kmer done" << std::endl;
		similarityToDistanceMatrix(distanceMatrix, KimuraDistance() );
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

		String<char> align;
		SEQAN_TASSERT(convertAlignment(gOut, align));
	}
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

	// Score object
	Score<int> score_type = Score<int>(4,-1,-1,-10);

	generatePrimaryLibrary(gAux, score_type, GlobalPairwise_Library() );
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
	String<double> distanceMatrix; 
	getKmerSimilarityMatrix(stringSet(gAux), distanceMatrix, 6, AAGroupsDayhoff() );
	similarityToDistanceMatrix(distanceMatrix, KimuraDistance() );

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

void Test_BaliBaseRef11() {
	
	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\Balibase\\bb3_release\\RV11\\");
#else
	// Linux
	String<char> in_path("/home/takifugu/rausch/Balibase/bb3_release/RV11/");
#endif

	for(int i = 1; i<39; ++i) {
		std::stringstream s;
		s << "BBS110";
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

	for(int i = 1; i<45; ++i) {
		std::stringstream s;
		s << "BBS120";
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

	for(int i = 1; i<42; ++i) {
		std::stringstream s;
		s << "BBS200";
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

	for(int i = 1; i<31; ++i) {
		std::stringstream s;
		s << "BBS300";
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

	for(int i = 1; i<17; ++i) {
		std::stringstream s;
		s << "BBS500";
		if (i < 10) s << '0';
		s << i;
		Test_TCoffeeFromFile(in_path,s.str().c_str(), "tfa");
	}
}



//////////////////////////////////////////////////////////////////////////////

void Test_GraphTCoffee() {
	Test_CompressedAlphabets();


	//Test_TCoffeeTmp();
	Test_TCoffeeGarfield();
	//Test_TCoffeeFromRandomSeq();
	//Test_BaliBaseRef11();
	//Test_BaliBaseRef12();
	//Test_BaliBaseRef20();
	//Test_BaliBaseRef30();
	//Test_BaliBaseRef40();
	//Test_BaliBaseRef50();
	
	//Test_TCoffeeFromLibrary(); 

	debug::verifyCheckpoints("projects/library/seqan/graph/graph_align_tcoffee.h");

}


}

#endif

