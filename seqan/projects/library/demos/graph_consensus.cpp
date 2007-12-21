#include <iostream>
#include <fstream>
#include <seqan/graph.h>

using namespace seqan;

int main(int argc, const char *argv[]) {
	if (argc < 4) return 0;

	String<char> readsPath = argv[1];
	String<char> hap1Path = argv[2];
	String<char> hap2Path = argv[3];

	std::cout << "Reads file: " << readsPath << std::endl;
	std::cout << "Haplotype1 file: " << hap1Path << std::endl;
	std::cout << "Haplotype2 file: " << hap2Path << std::endl;
	
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
	unsigned int avgReadLength = nucCount / seqCount;
	std::cout << "Number of reads: " << seqCount << ", Total number of nucleotides: " << nucCount << std::endl;
	std::cout << "Average read length: " << avgReadLength << std::endl;
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
	selectPairsForLibraryGeneration(g, names, pList, avgReadLength);
	Graph<Undirected<double> > pairGraph;
	//String<double> pairGraph;
	generatePrimaryLibrary(g, pList, pairGraph, score_type_global, Overlap_Library() );
	_alignTiming(startTime, "Overlap done: ");
	
	// Triplet library extension
	tripletLibraryExtension(g);
	_alignTiming(startTime, "Triplet done: ");

	// Guide Tree
	Graph<Tree<double> > guideTree;
	upgmaTree(pairGraph, guideTree);
	_alignTiming(startTime, "Guide tree done: ");

	// Perform a progressive alignment
	Graph<Alignment<TStringSet, void> > gOut(strSet);
	progressiveAlignment(g, guideTree, gOut);
	clear(guideTree);
	clear(g);
	_alignTiming(startTime, "Progressive alignment done: ");

	String<char> alignmentMatrix;
	String<unsigned int> coverage;
	String<char> gappedConsensus;
	TString ungappedConsensus;
	typedef std::map<unsigned int, Dna> TPolymorphismMap;
	TPolymorphismMap polyMap;
	unsigned int maxCoverage = 0;
	consensusAlignment(gOut, alignmentMatrix, coverage, gappedConsensus, ungappedConsensus, polyMap, maxCoverage);
	_alignTiming(startTime, "Consensus done: ");
	std::cout << "==============================" << std::endl;
	_alignTiming(bigbang, "Total time: ");
	std::cout << "==============================" << std::endl;

	// Make an alignment of haplotypes and consensus
	StringSet<TString, Owner<> > seqSet1;
	StringSet<TString, Owner<> > seqSet2;
	_alignImportSequences(hap1Path, "", "", seqSet1, names);
	_alignImportSequences(hap2Path, "", "", seqSet2, names);
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;
	TStringSet strSetAli;
	appendValue(strSetAli, seqSet1[0]);
	appendValue(strSetAli, seqSet2[0]);
	appendValue(strSetAli, ungappedConsensus);
	Score<int> score_type = Score<int>(5,-4,-4,-14);
	TGraph lib1(strSetAli);
	String<double> distanceMatrix;
	generatePrimaryLibrary(lib1, distanceMatrix, score_type, GlobalPairwise_Library() );
	tripletLibraryExtension(lib1);
	Graph<Tree<double> > aliGuideTree;
	upgmaTree(distanceMatrix, aliGuideTree);
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > aliOut(strSetAli);
	progressiveAlignment(lib1, aliGuideTree, aliOut);
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
	unsigned int recoveredSnp = 0;
	unsigned int totalSnpCount = 0;
	unsigned int coverageThreshold = 2;
	unsigned int recoveredMatchAndCovered = 0;
	unsigned int totalCoveredMatchCount = 0;
	TPolymorphismMap leftOver(polyMap);
	String<unsigned int> missedSNP;
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
			if ((value(align, 0*colLen+col) != value(align, 1*colLen+col)) &&
				(value(align, 0*colLen+col) != gap) &&
				(value(align, 1*colLen+col) != gap)) {
					TPolymorphismMap::const_iterator pMap = polyMap.find(posGappedConsensus);
					if (pMap != polyMap.end()) {
						Dna a = pMap->second;
						Dna b = gappedConsensus[posGappedConsensus];
						if (((value(align, 0*colLen+col) == a) &&
							 (value(align, 1*colLen+col) == b)) ||
							((value(align, 0*colLen+col) == b) &&
							(value(align, 1*colLen+col) == a))) {
								leftOver.erase(pMap->first);
								++recoveredSnp;
						} 
					} else {
						appendValue(missedSNP, col);
					}
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
		if ((value(align, 0*colLen+col) != value(align, 1*colLen+col)) &&
			(value(align, 0*colLen+col) != gap) &&
			(value(align, 1*colLen+col) != gap)) {
			++totalSnpCount;
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
	std::cout << "#SNPs in the consensus = " << polyMap.size() << std::endl;
	std::cout << "#True SNPs in the consensus / #SNPs between Hap1 and Hap2 = " << recoveredSnp << "/" << totalSnpCount << std::endl;
	std::cout << "Missed SNP positions (in final Alignment): ";
	for(unsigned int i = 0; i<length(missedSNP);++i) {
		std::cout << missedSNP[i] << ' ';
	}
	std::cout << std::endl;
	std::cout << "False SNP positions (in gapped Consensus): ";
	for(TPolymorphismMap::const_iterator leftOverP = leftOver.begin();leftOverP != leftOver.end(); ++leftOverP) {
		std::cout << leftOverP->first << ' ';
	}
	std::cout << std::endl;
	std::cout << std::endl;

	// Output alignment of reads
	write(std::cout,alignmentMatrix, gappedConsensus, ungappedConsensus, polyMap, maxCoverage, FastaReadFormat());
	std::cout << std::endl;
	std::cout << "Alignment of Haplotype 1, Haplotype 2, and Consensus: " << std::endl;
	std::cout << aliOut << std::endl;

	std::cout << "Alignment position, Position in ungapped consensus, Position in gapped consensus (Read alignment)" << std::endl;
	for(unsigned int i = 0; i<length(posit);++i) {
		std::cout << posit[i].i1 << ',' << posit[i].i2 << ',' << posit[i].i3 << std::endl;
	}
	return 0;
}
