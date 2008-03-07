#include <seqan/graph.h>
#include <seqan/modifier.h>

#include <iostream>
#include <fstream>


using namespace seqan;


template <typename TStringSet, typename TCargo, typename TSpec, typename TAlignmentMatrix, typename TConsensus, typename TCoverage, typename TOptions>
inline void 
evaluationOfReadAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const&,
						  TAlignmentMatrix const&, 
						  TConsensus const& gappedConsensus,
						  TCoverage const& coverage,
						  TOptions& cfgOpt) 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TAlignmentMatrix>::Type TValue;
	typedef typename Value<TString>::Type TAlphabet;
	typedef typename Id<TString>::Type TId;
	TValue gapChar = gapValue<TValue>();

	// Calculate ungapped consensus
	TString ungappedConsensus;
	for(TSize col = 0; col<length(gappedConsensus); ++col) {
		if (gappedConsensus[col] != gapChar) appendValue(ungappedConsensus, gappedConsensus[col]);
	}

	
	// Make an alignment of haplotypes and consensus
	StringSet<TString, Owner<> > seqSet;
	String<String<char> > names;
	_alignImportSequences(value(cfgOpt, "haplotypes"), seqSet, names);
	
	typedef StringSet<TString, Dependent<> > TDepStringSet;
	typedef Graph<Alignment<TDepStringSet, unsigned int> > TAliGraph;
	TStringSet strSetAli;
	appendValue(strSetAli, seqSet[0]);
	appendValue(strSetAli, seqSet[1]);
	appendValue(strSetAli, ungappedConsensus);
	Score<int> score_type = Score<int>(5,-4,-4,-14);
	TAliGraph lib1(strSetAli);
	String<Pair<TId, TId> > pList2;
	selectPairsForLibraryGeneration(lib1, pList2);
	String<double> distanceMatrix;
	generatePrimaryLibrary(lib1, pList2, distanceMatrix, score_type, GlobalPairwise_Library() );
	tripletLibraryExtension(lib1);
	Graph<Tree<double> > aliGuideTree;
	upgmaTree(distanceMatrix, aliGuideTree);
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > aliOut(strSetAli);
	progressiveAlignment(lib1, aliGuideTree, aliOut);
	std::cout << "Alignment of Haplotype 1, Haplotype 2, and Consensus: " << std::endl;
	std::cout << aliOut << std::endl;


	// Where are the differences?
	String<char> align;
	char gap = gapValue<char>();
	convertAlignment(aliOut, align);
	TSize nseq = length(stringSet(aliOut));
	TSize colLen = length(align) / nseq;
	typedef Triple<unsigned int, unsigned int, unsigned int> TTriple;
	String<TTriple> posit;
	unsigned int posGappedConsensus = 0;
	unsigned int posUngappedConsensus = 0;
	unsigned int recoveredMatch = 0;
	unsigned int totalMatchCount = 0;
	unsigned int coverageThreshold = 2;
	unsigned int recoveredMatchAndCovered = 0;
	unsigned int totalCoveredMatchCount = 0;
	String<unsigned int> missedMatch;
	std::set<unsigned int> covertAlignmentPositions;
	for(unsigned int col = 0; col < colLen; ++col) {
		char curC = value(align, 2*colLen+col);
		if (curC != gap) {
			unsigned int oldValue = posGappedConsensus;
			while(gappedConsensus[posGappedConsensus] == gap) ++posGappedConsensus;
			for(unsigned int i = col - (posGappedConsensus - oldValue); i<=col;++i) covertAlignmentPositions.insert(i);
			appendValue(posit, TTriple(col, posUngappedConsensus, posGappedConsensus));
			if ((value(align, 0*colLen+col) == value(align, 1*colLen+col)) &&
				(value(align, 0*colLen+col) == value(align, 2*colLen+col))) {
					++recoveredMatch;
					if (coverage[posGappedConsensus] > coverageThreshold) ++recoveredMatchAndCovered;
			}
			++posUngappedConsensus;
			++posGappedConsensus;
		}
		if ((value(align, 0*colLen+col) == value(align, 1*colLen+col)) &&
			(covertAlignmentPositions.find(col) != covertAlignmentPositions.end())) {
			if (value(align, 0*colLen+col) != value(align, 2*colLen+col)) appendValue(missedMatch, col);
			++totalMatchCount;
			if (coverage[posGappedConsensus - 1] > coverageThreshold) ++totalCoveredMatchCount;
		}
	}
	std::cout << std::endl;
	std::cout << "#Matches Consensus to Hap1 and Hap2 / #Matches between Hap1 and Hap2 = " << recoveredMatch << "/" << totalMatchCount << "  (with Coverage > 0)" << std::endl;
	std::cout << "Missed Match positions (in final Alignment): ";
	for(unsigned int i = 0; i<length(missedMatch);++i) {
		std::cout << missedMatch[i] << ' ';
	}
	std::cout << std::endl;
	std::cout << "#Matches Consensus to Hap1 and Hap2 / #Matches between Hap1 and Hap2 = " << recoveredMatchAndCovered << "/" << totalCoveredMatchCount << "  (with Coverage > 2)" << std::endl;
	std::cout << std::endl;

	std::cout << "Alignment position, Position in ungapped consensus, Position in gapped consensus (Read alignment)" << std::endl;
	for(unsigned int i = 0; i<length(posit);++i) {
		std::cout << posit[i].i1 << ',' << posit[i].i2 << ',' << posit[i].i3 << std::endl;
	}
}


int main(int argc, const char *argv[]) {
	bool deltas = true;


	//////////////////////////////////////////////////////////////////////////////
	// Command line parsing
	//////////////////////////////////////////////////////////////////////////////
	
	// Set the keys
	typedef String<char> TKey;
	typedef String<char> TValue;
	ConfigOptions<TKey, TValue> cfgOpt;
	TKey keys[] = {"reads", "haplotypes","outfile"};
	assignKeys(cfgOpt, keys, 3);
	assign(cfgOpt, "outfile", "readAlign.txt");
	// Help Message
	String<char> helpMsg;
	append(helpMsg, "Usage: graph_consensus -reads <FASTA File with Reads> [ARGUMENTS]\n");
	assignHelp(cfgOpt, helpMsg);
	if (argc < 2) {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}
	if (!parseCmdLine(argc, argv, cfgOpt)) return -1;

	// Timing variables
	clock_t bigbang, startTime;
	startTime = clock();
	bigbang = startTime;

	//////////////////////////////////////////////////////////////////////////////
	// Import sequences
	//////////////////////////////////////////////////////////////////////////////

	typedef String<Dna> TSequence;
	StringSet<TSequence, Owner<> > origStrSet;
	typedef String<char> TName;
	String<TName> names;
	unsigned int nucCount = _alignImportSequences(value(cfgOpt, "reads"), origStrSet, names);
	unsigned int nseq = length(origStrSet);
	unsigned int avgReadLength = nucCount / nseq;
	std::cout << "Number of reads: " << nseq << ", Total number of nucleotides: " << nucCount << std::endl;
	std::cout << "Average read length: " << avgReadLength << std::endl;
	_alignTiming(startTime, "Import sequences done: ");

	//// Debug code
	//for(unsigned int i = 0; i<nseq; ++i) {
	//	std::cout << '>' << names[i] << std::endl;
	//	std::cout << origStrSet[i] << std::endl;
	//}

	// Extract positions and reverse complement reads
	String<Pair<unsigned int, unsigned int> > begEndPos;
	layoutReads(names, origStrSet, begEndPos);
	// Make dependent string set
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);
	unsigned int sequenceThreshold = 100;

	//// Debug code
	//for(unsigned int i = 0; i<nseq; ++i) std::cout << begEndPos[i].i1 << ',' << begEndPos[i].i2 << std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Align the sequences
	//////////////////////////////////////////////////////////////////////////////

	typedef Graph<Alignment<TDepSequenceSet, unsigned int> > TGraph;
	typedef Id<TGraph>::Type TId;
	//Score<int> score_type = Score<int>(5,-4,-4,-14);
	Score<int> score_type = Score<int>(2,-6,-4,-9);
	
	// Generate a primary library, i.e., all global pairwise alignments
	TGraph g(strSet);
	String<Pair<TId, TId> > pList;
	selectPairsForLibraryGeneration(g, begEndPos, pList, avgReadLength);
	Graph<Undirected<double> > pairGraph;
	String<double> distanceMatrix;
	if (nseq < sequenceThreshold) generatePrimaryLibrary(g, pList, distanceMatrix, score_type, Overlap_Library() );
	else generatePrimaryLibrary(g, pList, pairGraph, score_type, Overlap_Library() );
	_alignTiming(startTime, "Overlap done: ");

	// Triplet library extension
	tripletLibraryExtension(g);
	_alignTiming(startTime, "Triplet done: ");

	// Guide Tree
	Graph<Tree<double> > guideTree;
	if (nseq < sequenceThreshold) slowNjTree(distanceMatrix, guideTree);
	else upgmaTree(pairGraph, guideTree);
	_alignTiming(startTime, "Guide tree done: ");
	clear(distanceMatrix);
	clear(pairGraph);

	// Perform a progressive alignment
	Graph<Alignment<TDepSequenceSet, void, WithoutEdgeId> > gOut(strSet);
	progressiveAlignment(g, guideTree, gOut);
	clearVertices(g);
	_alignTiming(startTime, "Progressive alignment done: ");

	// Build the read alignment matrix
	String<char> alignmentMatrix;
	String<unsigned int> coverage;
	String<char> gappedConsensus;
	String<Triple<unsigned int, unsigned int, unsigned int> > readBegEndRowPos;
	consensusAlignment(gOut, alignmentMatrix, readBegEndRowPos, coverage, gappedConsensus);
	_alignTiming(startTime, "Consensus done: ");

	//// Debug code
	//for(unsigned int i = 0; i<nseq; ++i) std::cout << readBegEndRowPos[i].i1 << ',' << readBegEndRowPos[i].i2 << ',' << readBegEndRowPos[i].i3 << std::endl;

	// Realign disrupted reads
	unsigned int numUnalignedReads = realignLowQualityReads(gOut, pList, readBegEndRowPos, g);
	if (numUnalignedReads > 0) {
		std::cout << "Disrupted reads: " << numUnalignedReads << std::endl;
		clearVertices(gOut);
		progressiveAlignment(g, guideTree, gOut);
		clear(alignmentMatrix);
		clear(coverage);
		clear(gappedConsensus);
		clear(readBegEndRowPos);
		consensusAlignment(gOut, alignmentMatrix, readBegEndRowPos, coverage, gappedConsensus);
		_alignTiming(startTime, "Realignment done: ");
	}
	clear(g);
	clear(guideTree);

	//// Debug code
	//TSequence consensus;
	//char gapChar = gapValue<char>();
	//for(unsigned int i = 0; i<length(gappedConsensus); ++i) {
	//	if (gappedConsensus[i] != gapChar) appendValue(consensus, gappedConsensus[i]);
	//}
	//std::cout << consensus << std::endl;

	
	//////////////////////////////////////////////////////////////////////////////
	// Output of aligned reads
	//////////////////////////////////////////////////////////////////////////////

	std::fstream strm;
	strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
	write(strm, gOut, alignmentMatrix, begEndPos, readBegEndRowPos, gappedConsensus, deltas, FastaReadFormat());
	strm.close();
	_alignTiming(startTime, "Output done: ");

	std::cout << "==============================" << std::endl;
	_alignTiming(bigbang, "Total time: ");
	std::cout << "==============================" << std::endl;

	if (length(value(cfgOpt, "haplotypes"))) evaluationOfReadAlignment(gOut, alignmentMatrix, gappedConsensus, coverage, cfgOpt);

	return 0;
}
