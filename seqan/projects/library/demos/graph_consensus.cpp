#include <seqan/graph.h>
#include <seqan/modifier.h>

#include <iostream>
#include <fstream>


using namespace seqan;


template<typename TName>
class ConfigOptions {
public:
	TName readsPath;
	TName hapPath;
	TName outfile;
	bool evaluationMode;

	ConfigOptions() : outfile("readAlign.txt"), evaluationMode(false) {}
};

inline bool
printErrorForCmd()
{
	std::cerr << "Usage: graph_consensus -reads <read file> [OPTIONS]" << std::endl;
	return false;
}

template<typename TConfigOptions>
inline bool
parseCommandLine(int argc, const char *argv[], TConfigOptions& cfgOpt) {
	if (argc <= 2) return printErrorForCmd();
	for(int i = 1; i < argc; ++i) {
		if (argv[i][0] == '-') {
			// parse option
			if (strcmp(argv[i], "-reads")==0) {
				if (i + 1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.readsPath = argv[i];
				continue;
			}
			if (strcmp(argv[i], "-haplotypes")==0) {
				if (i + 1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.hapPath = argv[i];
				cfgOpt.evaluationMode = true;
				continue;
			}
			if (strcmp(argv[i], "-outfile")==0) {
				if (i + 1 == argc) return printErrorForCmd();
				++i;
				cfgOpt.outfile = argv[i];
				continue;
			}
		}
	}
	return true;
}

template <typename TStringSet, typename TCargo, typename TSpec, typename TAlignmentMatrix, typename TConsensus, typename TOptions>
inline void 
evaluationOfReadAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						  TAlignmentMatrix const& alignmentMatrix, 
						  TConsensus const& gappedConsensus,
						  TOptions const& cfgOpt) 
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

	typedef std::map<unsigned int, Dna> TPolymorphismMap;
	TPolymorphismMap polyMap;
	
	// Make an alignment of haplotypes and consensus
	StringSet<TString, Owner<> > seqSet;
	String<String<char> > names;
	_alignImportSequences(cfgOpt.hapPath, "", "", seqSet, names);
	
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int> > TAliGraph;
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

	//String<char> align;
	//char gap = gapValue<char>();
	//convertAlignment(aliOut, align);
	//TSize nseq = length(stringSet(aliOut));
	//TSize colLen = length(align) / nseq;
	//typedef Triple<unsigned int, unsigned int, unsigned int> TTriple;
	//String<TTriple> posit;
	//unsigned int posGappedConsensus = 0;
	//unsigned int posUngappedConsensus = 0;
	//unsigned int recoveredMatch = 0;
	//unsigned int totalMatchCount = 0;
	//unsigned int recoveredSnp = 0;
	//unsigned int totalSnpCount = 0;
	//unsigned int coverageThreshold = 2;
	//unsigned int recoveredMatchAndCovered = 0;
	//unsigned int totalCoveredMatchCount = 0;
	//TPolymorphismMap leftOver(polyMap);
	//String<unsigned int> missedSNP;
	//String<unsigned int> missedMatch;
	//std::set<unsigned int> covertAlignmentPositions;
	//for(unsigned int col = 0; col < colLen; ++col) {
	//	char curC = value(align, 2*colLen+col);
	//	if (curC != gap) {
	//		unsigned int oldValue = posGappedConsensus;
	//		while(gappedConsensus[posGappedConsensus] == gap) ++posGappedConsensus;
	//		for(unsigned int i = col - (posGappedConsensus - oldValue); i<=col;++i) covertAlignmentPositions.insert(i);
	//		appendValue(posit, TTriple(col, posUngappedConsensus, posGappedConsensus));
	//		if ((value(align, 0*colLen+col) == value(align, 1*colLen+col)) &&
	//			(value(align, 0*colLen+col) == value(align, 2*colLen+col))) {
	//				++recoveredMatch;
	//				if (coverage[posGappedConsensus] > coverageThreshold) ++recoveredMatchAndCovered;
	//		}
	//		if ((value(align, 0*colLen+col) != value(align, 1*colLen+col)) &&
	//			(value(align, 0*colLen+col) != gap) &&
	//			(value(align, 1*colLen+col) != gap)) {
	//				TPolymorphismMap::const_iterator pMap = polyMap.find(posGappedConsensus);
	//				if (pMap != polyMap.end()) {
	//					Dna a = pMap->second;
	//					Dna b = gappedConsensus[posGappedConsensus];
	//					if (((value(align, 0*colLen+col) == a) &&
	//						 (value(align, 1*colLen+col) == b)) ||
	//						((value(align, 0*colLen+col) == b) &&
	//						(value(align, 1*colLen+col) == a))) {
	//							leftOver.erase(pMap->first);
	//							++recoveredSnp;
	//					} 
	//				} else {
	//					appendValue(missedSNP, col);
	//				}
	//		}
	//		++posUngappedConsensus;
	//		++posGappedConsensus;
	//	}
	//	if ((value(align, 0*colLen+col) == value(align, 1*colLen+col)) &&
	//		(covertAlignmentPositions.find(col) != covertAlignmentPositions.end())) {
	//		if (value(align, 0*colLen+col) != value(align, 2*colLen+col)) appendValue(missedMatch, col);
	//		++totalMatchCount;
	//		if (coverage[posGappedConsensus - 1] > coverageThreshold) ++totalCoveredMatchCount;
	//	}
	//	if ((value(align, 0*colLen+col) != value(align, 1*colLen+col)) &&
	//		(value(align, 0*colLen+col) != gap) &&
	//		(value(align, 1*colLen+col) != gap)) {
	//		++totalSnpCount;
	//	}
	//}
	//std::cout << std::endl;
	//std::cout << "#Matches Consensus to Hap1 and Hap2 / #Matches between Hap1 and Hap2 = " << recoveredMatch << "/" << totalMatchCount << "  (with Coverage > 0)" << std::endl;
	//std::cout << "Missed Match positions (in final Alignment): ";
	//for(unsigned int i = 0; i<length(missedMatch);++i) {
	//	std::cout << missedMatch[i] << ' ';
	//}
	//std::cout << std::endl;
	//std::cout << "#Matches Consensus to Hap1 and Hap2 / #Matches between Hap1 and Hap2 = " << recoveredMatchAndCovered << "/" << totalCoveredMatchCount << "  (with Coverage > 2)" << std::endl;
	//std::cout << "#SNPs in the consensus = " << polyMap.size() << std::endl;
	//std::cout << "#True SNPs in the consensus / #SNPs between Hap1 and Hap2 = " << recoveredSnp << "/" << totalSnpCount << std::endl;
	//std::cout << "Missed SNP positions (in final Alignment): ";
	//for(unsigned int i = 0; i<length(missedSNP);++i) {
	//	std::cout << missedSNP[i] << ' ';
	//}
	//std::cout << std::endl;
	//std::cout << "False SNP positions (in gapped Consensus): ";
	//for(TPolymorphismMap::const_iterator leftOverP = leftOver.begin();leftOverP != leftOver.end(); ++leftOverP) {
	//	std::cout << leftOverP->first << ' ';
	//}
	//std::cout << std::endl;
	//std::cout << std::endl;

	//// Output alignment of reads
	////write(std::cout,alignmentMatrix, gappedConsensus, ungappedConsensus, polyMap, maxCoverage, FastaReadFormat());
	//std::cout << std::endl;
	//std::cout << "Alignment of Haplotype 1, Haplotype 2, and Consensus: " << std::endl;
	//std::cout << aliOut << std::endl;

	//std::cout << "Alignment position, Position in ungapped consensus, Position in gapped consensus (Read alignment)" << std::endl;
	//for(unsigned int i = 0; i<length(posit);++i) {
	//	std::cout << posit[i].i1 << ',' << posit[i].i2 << ',' << posit[i].i3 << std::endl;
	//}
}


int main(int argc, const char *argv[]) {
	ConfigOptions<String<char> > cfgOpt;
	if (!parseCommandLine(argc, argv, cfgOpt)) return -1;

	// Timing variables
	clock_t bigbang, startTime;

	// Import Sequences
	typedef String<Dna> TString;
	StringSet<TString, Owner<> > origStrSet;
	String<String<char> > names;
	startTime = clock();
	bigbang = startTime;
	unsigned int nucCount = _alignImportSequences(cfgOpt.readsPath, "", "", origStrSet, names);
	unsigned int seqCount = length(origStrSet);
	unsigned int avgReadLength = nucCount / seqCount;
	std::cout << "Number of reads: " << seqCount << ", Total number of nucleotides: " << nucCount << std::endl;
	std::cout << "Average read length: " << avgReadLength << std::endl;
	_alignTiming(startTime, "Import sequences done: ");

	// Extract positions and reverseComplement reads
	String<Pair<unsigned int, unsigned int> > begEndPos;
	layoutReads(names, origStrSet, begEndPos);

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
	selectPairsForLibraryGeneration(g, begEndPos, pList, avgReadLength);
	Graph<Undirected<double> > pairGraph;
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

	// Build the read alignment matrix
	String<char> alignmentMatrix;
	String<unsigned int> coverage;
	String<char> gappedConsensus;
	String<Triple<unsigned int, unsigned int, unsigned int> > readBegEndRowPos;
	consensusAlignment(gOut, alignmentMatrix, readBegEndRowPos, coverage, gappedConsensus);
	_alignTiming(startTime, "Consensus done: ");
	
	// Write to file
	std::fstream strm;
	strm.open(toCString(cfgOpt.outfile), std::ios_base::out | std::ios_base::trunc);
	write(strm, gOut, alignmentMatrix, begEndPos, readBegEndRowPos, gappedConsensus, FastaReadFormat());
	strm.close();
	_alignTiming(startTime, "Output done: ");

	std::cout << "==============================" << std::endl;
	_alignTiming(bigbang, "Total time: ");
	std::cout << "==============================" << std::endl;

	if (cfgOpt.evaluationMode) evaluationOfReadAlignment(gOut, alignmentMatrix, gappedConsensus, cfgOpt);

	return 0;
}
