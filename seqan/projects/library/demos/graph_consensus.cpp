#include <seqan/consensus.h>


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
	_loadSequences(value(cfgOpt, "haplotypes"), seqSet, names);
	
	typedef StringSet<TString, Dependent<> > TDepStringSet;
	typedef Graph<Alignment<TDepStringSet, unsigned int> > TAliGraph;
	TStringSet strSetAli;
	for(TSize i = 0; i < length(seqSet); ++i) appendValue(strSetAli, seqSet[i]);
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
	std::cout << "Alignment of Haplotypes and Consensus: " << std::endl;
	std::cout << aliOut << std::endl;


	// Where are the differences?
	String<char> align;
	char gap = gapValue<char>();
	convertAlignment(aliOut, align);
	TSize nseq = length(stringSet(aliOut));
	TSize colLen = length(align) / nseq;
	unsigned int posGappedConsensus = 0;
	unsigned int recoveredMatch = 0;
	unsigned int totalMatchCount = 0;
	unsigned int coverageThreshold = 2;
	unsigned int recoveredMatchAndCovered = 0;
	unsigned int totalCoveredMatchCount = 0;
	String<unsigned int> missedMatch;
	std::set<unsigned int> covertAlignmentPositions;

	// Walk through the alignment
	for(unsigned int col = 0; col < colLen; ++col) {
		char curC = value(align, (length(strSetAli) - 1)*colLen+col);
		if (curC != gap) {
			// The consensus was gapped in the alignment
			unsigned int oldValue = posGappedConsensus;
			while(gappedConsensus[posGappedConsensus] == gap) ++posGappedConsensus;
			for(unsigned int i = col - (posGappedConsensus - oldValue); i<=col;++i) covertAlignmentPositions.insert(i);
			bool recov = true;
			for(TSize i = 1; i<length(strSetAli); ++i) {
				if (value(align, 0*colLen+col) != value(align, i*colLen+col)) {
					recov = false;
					break;
				}
			}
			if (recov) {
				++recoveredMatch;
				if (coverage[posGappedConsensus] > coverageThreshold) ++recoveredMatchAndCovered;
			}
			++posGappedConsensus;
		}
		bool allMatch = true;
		for(TSize i = 1; i<length(strSetAli) - 1; ++i) {
			if (value(align, 0*colLen+col) != value(align, 1*colLen+col)) {
				allMatch = false;
				break;
			}
		}
		if ((allMatch) &&
			(covertAlignmentPositions.find(col) != covertAlignmentPositions.end())) {
			if (value(align, 0*colLen+col) != value(align, (length(strSetAli) - 1)*colLen+col)) appendValue(missedMatch, col);
			++totalMatchCount;
			if (coverage[posGappedConsensus - 1] > coverageThreshold) ++totalCoveredMatchCount;
		}
	}
	std::cout << std::endl;
	std::cout << "#Matches Consensus to Haplotypes / #Matches between Haplotypes = " << recoveredMatch << "/" << totalMatchCount << "  (with Coverage > 0)" << std::endl;
	std::cout << "Missed Match positions (in final Alignment): ";
	for(unsigned int i = 0; i<length(missedMatch);++i) {
		std::cout << missedMatch[i] << ' ';
	}
	std::cout << std::endl;
	std::cout << "#Matches Consensus to Haplotypes / #Matches between Haplotypes = " << recoveredMatchAndCovered << "/" << totalCoveredMatchCount << "  (with Coverage > 2)" << std::endl;
	std::cout << std::endl;
}



template <typename TOptions, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
convertSimulationFile(ReadStore<TAlphabet, TSpec>& readSt,
					  TFragmentStore& frgSt,
					  TLibraryStore& libSt,
					  TContigStore& ctgSt,
					  TOptions& cfgOpt)
{
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef String<char> TName;

	// Convert Libraries
	String<char> filePath = value(cfgOpt, "reads");
	appendValue(filePath, 'L');
	StringSet<String<char>, Owner<> > libSet;
	String<TName> libIds;
	_loadSequences(filePath, libSet, libIds);
	for(TSize i = 0; i<length(libSet); ++i) {
		std::stringstream input;
		input << "L" << i;
		String<char> tmp(input.str().c_str());
		appendValue(libSt.data_names, tmp);
		TSize brPoint = 0;
		TSize mean;
		TSize std;
		for(TSize k = 0; k<length(value(libSet, i)); ++k) {
			if (value(value(libSet,i),k) == ',') brPoint = k;
		}
		String<char> inf1 = infix(value(libSet,i), 0, brPoint);
		String<char> inf2 = infix(value(libSet,i), brPoint+1, length(value(libSet,i)));
		std::stringstream ssStream1(toCString(inf1));
		ssStream1 >> mean; 
		std::stringstream ssStream2(toCString(inf2));
		ssStream2 >> std;
		appendValue(libSt.data_mean, mean);
		appendValue(libSt.data_std, std);
		++libSt.data_pos_count;
	}

	// Convert Fragments
	filePath = value(cfgOpt, "reads");
	appendValue(filePath, 'F');
	StringSet<String<char>, Owner<> > fragSet;
	String<TName> fragIds;
	_loadSequences(filePath, fragSet, fragIds);
	for(TSize i = 0; i<length(fragSet); ++i) {
		std::stringstream input;
		input << "F" << i;
		String<char> tmp(input.str().c_str());
		appendValue(frgSt.data_names, tmp);
		TSize brPoint1 = 0;
		TSize brPoint2 = 0;
		TSize libId;
		for(TSize k = 0; k<length(value(fragIds, i)); ++k) {
			if (value(value(fragIds,i),k) == '=') brPoint1 = k + 1;
			if (value(value(fragIds,i),k) == ']') brPoint2 = k;
		}
		String<char> inf = infix(value(fragIds,i), brPoint1, brPoint2);
		std::stringstream ssStream(toCString(inf));
		ssStream >> libId; 
		appendValue(frgSt.data_lib_id, libId - 1);
		TSize brPoint = 0;
		TSize rds1;
		TSize rds2;
		for(TSize k = 0; k<length(value(fragSet, i)); ++k) {
			if (value(value(fragSet,i),k) == ',') brPoint = k;
		}
		String<char> inf1 = infix(value(fragSet,i), 0, brPoint);
		String<char> inf2 = infix(value(fragSet,i), brPoint+1, length(value(fragSet,i)));
		std::stringstream ssStream1(toCString(inf1));
		ssStream1 >> rds1; 
		std::stringstream ssStream2(toCString(inf2));
		ssStream2 >> rds2;
		if (rds1 < rds2) appendValue(frgSt.data_rds, Pair<TSize, TSize>(rds1 - 1, rds2 - 1));
		else appendValue(frgSt.data_rds, Pair<TSize, TSize>(rds1 - 1, rds2 - 1));
		++frgSt.data_pos_count;
	}

	// Convert reads
	typedef String<TAlphabet> TSequence;
	StringSet<String<TAlphabet>, Owner<> > origStrSet;
	String<TName> names;
	clear(filePath);
	filePath = value(cfgOpt, "reads");
	_loadSequences(filePath, origStrSet, names);
	for(TSize i = 0; i<length(origStrSet); ++i) {
		TSize begRead;
		if (readSt.data_pos_count == 0) begRead = 0;
		else begRead = (value(readSt.data_begin_end, readSt.data_pos_count - 1)).i2;
		TSize endRead = begRead + length(value(origStrSet, i));
		typedef typename Iterator<String<TAlphabet> >::Type TSeqIter;
		TSeqIter itSeq = begin(value(origStrSet, i));
		TSeqIter itSeqEnd = end(value(origStrSet, i));
		for(;itSeq != itSeqEnd; ++itSeq) appendValue(readSt.data_reads, *itSeq);
		itSeq = begin(value(origStrSet, i));
		for(;itSeq != itSeqEnd; ++itSeq) appendValue(readSt.data_qualities, char(48+60));
		appendValue(readSt.data_begin_end, Pair<TSize,TSize>(begRead, endRead));
		TSize brPoint1 = 0;
		TSize brPoint2 = 0;
		TSize fragId = 0;
		for(TSize k = 0; k<length(value(names, i)); ++k) {
			if (value(value(names,i),k) == '=') brPoint1 = k + 1;
			if (value(value(names,i),k) == ']') brPoint2 = k;
		}
		if (brPoint1 != brPoint2) {
			String<char> inf = infix(value(names,i), brPoint1, brPoint2);
			std::stringstream ssStream(toCString(inf));
			ssStream >> fragId; 
			appendValue(readSt.data_frg_id, fragId - 1);
		}
		std::stringstream input;
		input << "R" << i;
		String<char> tmp(input.str().c_str());
		appendValue(readSt.data_names, tmp);
		appendValue(readSt.data_clr, Pair<TSize, TSize>(0, length(value(origStrSet, i))));
		++readSt.data_pos_count;
	}

	// Write minimal contig
	appendValue(ctgSt.data_reads, String<GappedRead<> >());
	for(TSize i = 0; i<length(names); ++i) {
		GappedRead<> gapRead;
		gapRead.data_source = i;
		TSize brPoint1 = 0;
		TSize brPoint2 = length(value(names, i));
		for(TSize k = 0; k<length(value(names, i)); ++k) {
			if (value(value(names,i),k) == ',') brPoint1 = k;
			if (value(value(names,i),k) == '[') {
				brPoint2 = k;
				break;
			}
		}
		TSize posI = 0;
		TSize posJ = 0;
		String<char> inf1 = infix(value(names,i), 0, brPoint1);
		String<char> inf2 = infix(value(names,i), brPoint1+1, brPoint2);
		std::stringstream ssStream1(toCString(inf1));
		ssStream1 >> posI; 
		std::stringstream ssStream2(toCString(inf2));
		ssStream2 >> posJ;
		if (posI < posJ) {
			gapRead.data_offset = posI;
			gapRead.data_clr = Pair<TSize, TSize>(0, length(value(origStrSet, i)));
		} else {
			gapRead.data_offset = posJ;
			gapRead.data_clr = Pair<TSize, TSize>(length(value(origStrSet, i)), 0);
		}
		appendValue(value(ctgSt.data_reads, ctgSt.data_pos_count), gapRead);
	}
	String<char> tmp = "C0";
	appendValue(ctgSt.data_names, tmp);
	StringSet<String<char>, Owner<> > seqSet;
	String<TName> seqNames;
	filePath = value(cfgOpt, "reads");
	appendValue(filePath, 'S');
	_loadSequences(filePath, seqSet, seqNames);
	if (!empty(seqSet)) {
		for(TSize l = 0;l < length(seqSet[0]); ++l) appendValue(ctgSt.data_contig, value(value(seqSet, 0), l));
		for(TSize l = 0;l < length(seqSet[0]); ++l) appendValue(ctgSt.data_quality, 'D');
		appendValue(ctgSt.data_begin_end, Pair<TSize,TSize>(0, length(seqSet[0])));
		++ctgSt.data_pos_count;
	}
}


int main(int argc, const char *argv[]) {
	//////////////////////////////////////////////////////////////////////////////
	// Command line parsing
	//////////////////////////////////////////////////////////////////////////////
	
	// Set the keys
	typedef String<char> TKey;
	typedef String<char> TValue;
	typedef Size<TValue>::Type TSize;
	ConfigOptions<TKey, TValue> cfgOpt;
	TKey keys[] = {"afg", "reads", "haplotypes","outfile", "output", "convert"};
	assignKeys(cfgOpt, keys, 6);
	assign(cfgOpt, "output", "seqan");
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

	//////////////////////////////////////////////////////////////////////////////
	// Read simulation file, afg file or celera consensus file
	//////////////////////////////////////////////////////////////////////////////
	ReadStore<> readSt;
	FrgStore<> frgSt;
	LibStore<> libSt;
	CtgStore<> ctgSt;

	if (!empty(value(cfgOpt, "reads"))) {
		convertSimulationFile(readSt, frgSt, libSt, ctgSt, cfgOpt);
	} else if (!empty(value(cfgOpt, "afg"))) {
		// Read Amos
		std::fstream strmReads;
		strmReads.open(toCString(value(cfgOpt, "afg")), std::ios_base::in | std::ios_base::binary);
		read(strmReads,readSt,frgSt,libSt,ctgSt,Amos());	
		strmReads.close();
	} else {
		return -1;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Just convert the input file
	//////////////////////////////////////////////////////////////////////////////
	if (!empty(value(cfgOpt, "convert"))) {
		if (value(cfgOpt, "convert") == "afg") {
			std::fstream strmWrite;
			strmWrite.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
			write(strmWrite,readSt,frgSt,libSt,ctgSt,Amos());	
			strmWrite.close();
		} else if (value(cfgOpt, "convert") == "frg") {
			std::fstream strmWrite;
			strmWrite.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
			write(strmWrite,readSt,frgSt,libSt,ctgSt,CeleraFrg());	
			strmWrite.close();
		} else if (value(cfgOpt, "convert") == "cgb") {
			std::fstream strmWrite;
			strmWrite.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
			write(strmWrite,readSt,frgSt,libSt,ctgSt,CeleraCgb());	
			strmWrite.close();
		}
		return 0;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Import sequences
	//////////////////////////////////////////////////////////////////////////////

	// Timing variables
	clock_t bigbang, startTime;
	startTime = clock();
	bigbang = startTime;

	typedef String<Dna> TSequence;
	StringSet<TSequence, Owner<> > origStrSet;
	String<Pair<unsigned int, unsigned int> > begEndPos;
	unsigned int currentContig = 0;
	loadReadsClr(readSt, ctgSt, currentContig, origStrSet, begEndPos);
	unsigned int nseq = length(origStrSet);
	unsigned int avgReadLength = 0;
	for(unsigned int i = 0; i<length(origStrSet);++i) avgReadLength += length(value(origStrSet,i));
	avgReadLength /= nseq;
	std::cout << "Number of reads: " << nseq << std::endl;
	std::cout << "Average read length: " << avgReadLength << std::endl;
	_alignTiming(startTime, "Import sequences done: ");
		
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

	//// Realign disrupted reads
	//unsigned int numUnalignedReads = realignLowQualityReads(gOut, pList, readBegEndRowPos, g);
	//if (numUnalignedReads > 0) {
	//	std::cout << "Disrupted reads: " << numUnalignedReads << std::endl;
	//	_alignTiming(startTime, "Disrupted read scan done: ");
	//	clearVertices(gOut);
	//	progressiveAlignment(g, guideTree, gOut);
	//	clear(alignmentMatrix);
	//	clear(coverage);
	//	clear(gappedConsensus);
	//	clear(readBegEndRowPos);
	//	consensusAlignment(gOut, alignmentMatrix, readBegEndRowPos, coverage, gappedConsensus);
	//	_alignTiming(startTime, "Realignment done: ");
	//}
	//clear(g);
	//clear(guideTree);

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

	if (value(cfgOpt, "output") == "seqan") {
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm, gOut, alignmentMatrix, begEndPos, readBegEndRowPos, gappedConsensus, FastaReadFormat());
		strm.close();
		_alignTiming(startTime, "Output done: ");
	} else if (value(cfgOpt, "output") == "afg") {
		TSize len = length(gappedConsensus);
		TSize begContig = (value(ctgSt.data_begin_end, ctgSt.data_pos_count - 1)).i2;
		TSize endContig = begContig + length(gappedConsensus);
		for(TSize i = 0; i<len; ++i) appendValue(ctgSt.data_contig, gappedConsensus[i]);
		for(TSize i = 0; i<len; ++i) appendValue(ctgSt.data_quality, 'D');
		value(ctgSt.data_begin_end, currentContig) = Pair<TSize, TSize>(begContig, endContig);

		for(TSize i = 0; i < length(readBegEndRowPos); ++i) {
			String<TSize> gaps;
			TSize letterCount = 0;
			TSize gapCount = 0;
			for(TSize column = (readBegEndRowPos[i]).i1; column<(readBegEndRowPos[i]).i2; ++column) {
				if (value(alignmentMatrix, (readBegEndRowPos[i]).i3 * len + column) == '-') {
					++gapCount;
					appendValue(gaps, letterCount);
				} else ++letterCount;
			}
			value(value(ctgSt.data_reads, currentContig), i).data_gap = gaps;
			value(value(ctgSt.data_reads, currentContig), i).data_offset = (readBegEndRowPos[i]).i1;
		}

		// Write Amos
		std::fstream strmWrite;
		strmWrite.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strmWrite,readSt,frgSt,libSt,ctgSt,Amos());	
		strmWrite.close();
	}

	std::cout << "==============================" << std::endl;
	_alignTiming(bigbang, "Total time: ");
	std::cout << "==============================" << std::endl;

	if (length(value(cfgOpt, "haplotypes"))) evaluationOfReadAlignment(gOut, alignmentMatrix, gappedConsensus, coverage, cfgOpt);

	return 0;
}
