#define SEQAN_PROFILE

#include <seqan/consensus.h>
//#include <seqan/misc/misc_random.h>

#include <iostream>
#include <fstream>


using namespace seqan;

/*

template <typename TOptions>
inline void 
evaluationOfConsensusAlignment(TOptions& cfgOpt) 
{
	typedef typename Size<TOptions>::Type TSize;
	typedef String<Dna> TSequence;
	typedef String<char> TName;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<TName> names;

	// Read the sequences
	std::fstream strm;
	strm.open(toCString(value(cfgOpt, "evaluation")), std::ios_base::in | std::ios_base::binary);
	read(strm,origStrSet,names,FastaAlign());	
	strm.close();

	// Make a dependent StringSet
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);
	
	// Score the alignment
	typedef Score<int> TScore;
	typedef Value<TScore>::Type TScoreValue;
	TScore scType = Score<int>(2,-6,-4,-9);

	// Read the alignment
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	TGraph g(strSet);
	std::fstream strm_lib;
	strm_lib.open(toCString(value(cfgOpt, "evaluation")), std::ios_base::in | std::ios_base::binary);
	read(strm_lib,g,names,scType,FastaAlign());	
	strm_lib.close();

	// Get the alignment matrix
	typedef char TValue;
	typedef String<TValue> TAlignmentMatrix;
	TAlignmentMatrix alignmentMatrix;
	if (convertAlignment(g, alignmentMatrix)) {

		// Find all differences
		TSize nseq = length(stringSet(g));
		TSize len = length(alignmentMatrix) / nseq;
		String<bool> diff;
		fill(diff, len, false);
		TSize countDiff = 0;
		for(TSize col = 0; col<len; ++col) {
			TValue c = value(alignmentMatrix, col);
			for(TSize row = 1; row<nseq; ++row) {
				if (value(alignmentMatrix, row * len + col) != c) {
					value(diff, col) = true;
					break;
				}
			}
		}

		// Output differences
		TSize col = 0;
		while(col < len) {
			// New diff?
			if (value(diff, col)) {
				// Extend diff
				TSize stopCol = col + 1;
				while (((stopCol < len) && (value(diff, stopCol)))) ++stopCol;

				// Output diff
				std::cout << "Difference: " << std::endl;
				TSize from = 0;
				TSize windowSize = 10;
				if (col > windowSize) from = col - windowSize;
				TSize to = len;
				if (stopCol + windowSize < len) to = stopCol + windowSize;
				for(TSize localRow = 0; localRow<nseq; ++localRow) {
					for(TSize localCol = from; localCol<to; ++localCol) {
						std::cout << value(alignmentMatrix, localRow * len + localCol);
					}
					std::cout << "  <<" << value(names, localRow) << std::endl;
				}
				col = stopCol;
			} else {
				++col;
			}
		}
	}
}
*/

template <typename TName, typename TKey>
inline void
__findThis(TName const& searchText,
		   TKey const& keyText,
		   TName& result) 
{
	typedef typename Iterator<TName>::Type TNameIter;
	TName key = keyText;
	TNameIter nameIt = begin(searchText);
	TNameIter nameItEnd = end(searchText);
	for(;nameIt != nameItEnd;++nameIt) {
		TNameIter keyIt = begin(key);
		TNameIter keyItEnd = end(key);
		TNameIter localNameIt = nameIt;
		while((localNameIt != nameItEnd) && (keyIt != keyItEnd) && (value(localNameIt) == value(keyIt))) {
			++keyIt;
			++localNameIt;
		}
		if (keyIt == keyItEnd) {
			while ((value(localNameIt) != ']') && (value(localNameIt) != ',')) {
				appendValue(result, value(localNameIt));
				++localNameIt;
			}
			break;
		}
	}
}

template <typename TOptions, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline bool
convertSimpleReadFile(ReadStore<TAlphabet, TSpec>& readSt,
					  TFragmentStore& frgSt,
					  TLibraryStore& libSt,
					  TContigStore& ctgSt,
					  TOptions& cfgOpt,
					  bool moveToFront)
{
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef String<char> TName;



	// Convert reads
	typedef String<TAlphabet> TSequence;
	StringSet<String<TAlphabet>, Owner<> > origStrSet;
	String<TName> names;
	String<char>  filePath;
	filePath = value(cfgOpt, "reads");
	TSize count = _loadSequences(filePath, origStrSet, names);
	if (count == 0) return false;
	else {
		for(TSize i = 0; i<length(origStrSet); ++i) {
			// Get the begin and end position
			TSize begRead;
			if (readSt.data_pos_count == 0) begRead = 0;
			else begRead = (value(readSt.data_begin_end, readSt.data_pos_count - 1)).i2;
			TSize endRead = begRead + length(value(origStrSet, i));
			appendValue(readSt.data_begin_end, Pair<TSize,TSize>(begRead, endRead));

			// Save the read
			typedef typename Iterator<String<TAlphabet> >::Type TSeqIter;
			TSeqIter itSeq = begin(value(origStrSet, i));
			TSeqIter itSeqEnd = end(value(origStrSet, i));
			for(;itSeq != itSeqEnd; ++itSeq) {
				appendValue(readSt.data_reads, *itSeq);
				appendValue(readSt.data_qualities, char(48+60)); // Dummy quality value
			}
			
			// Get the fragment id
			TSize fragId = i;
			TName valKey;
			__findThis(value(names, i), "fragId=", valKey);
			if (!empty(valKey)) fragId = _stringToNumber<TSize>(valKey);
			appendValue(readSt.data_frg_id, fragId);
			
			// Build eid string
			std::stringstream input;
			input << "R" << i;

			// Possibly append repeatId
			clear(valKey);
			__findThis(value(names, i), "repeatId=", valKey);
			if (empty(valKey)) valKey = "0";
			input << "-" << valKey;
			String<char> tmp(input.str().c_str());
			appendValue(readSt.data_names, tmp);

			// Append a clear range
			appendValue(readSt.data_clr, Pair<TSize, TSize>(0, length(value(origStrSet, i))));
			++readSt.data_pos_count;
		}
	}




	// Write minimal contig

	// Create a contig name
	String<char> tmp = "C0";
	appendValue(ctgSt.data_names, tmp);

	// Create a string of gapped reads
	appendValue(ctgSt.data_reads, String<GappedRead<> >());
	TSize lastLetter = 0;
	for(TSize i = 0; i<length(names); ++i) {
		// Add the original read id
		GappedRead<> gapRead;
		gapRead.data_source = i;

		// Add the offset and clear range
		typedef typename Iterator<TName>::Type TNameIter;
		TNameIter nameIt = begin(value(names, i));
		TNameIter nameItEnd = end(value(names, i));
		bool before = true;
		String<char> inf1; String<char> inf2;
		for(;nameIt != nameItEnd; goNext(nameIt)) {
			if (value(nameIt) == ',') before = false;
			else {
				if (before) appendValue(inf1, value(nameIt));
				else appendValue(inf2, value(nameIt));
			}
		}
		TSize posI = _stringToNumber<TSize>(inf1);
		TSize posJ = _stringToNumber<TSize>(inf2);
		if (posI < posJ) {
			if (posJ > lastLetter) lastLetter = posJ;
			gapRead.data_offset = posI;
			gapRead.data_clr = Pair<TSize, TSize>(0, length(value(origStrSet, i)));
		} else {
			if (posI > lastLetter) lastLetter = posI;
			gapRead.data_offset = posJ;
			gapRead.data_clr = Pair<TSize, TSize>(length(value(origStrSet, i)), 0);
		}
		appendValue(value(ctgSt.data_reads, ctgSt.data_pos_count), gapRead);
	}

	// Create a contig sequence	
	StringSet<String<char>, Owner<> > seqSet;
	String<TName> seqNames;
	clear(filePath);
	if (length(value(cfgOpt, "source"))) {
		filePath = value(cfgOpt, "source");
	} else {
		filePath = value(cfgOpt, "reads");
		appendValue(filePath, 'S');
	}
	count = _loadSequences(filePath, seqSet, seqNames);
	if (count != 0) {
		for(TSize l = 0;l < length(value(seqSet,0)); ++l) {
			appendValue(ctgSt.data_contig, value(value(seqSet, 0), l));
			appendValue(ctgSt.data_quality, char(48+60)); // Dummy quality value
		}
		appendValue(ctgSt.data_begin_end, Pair<TSize,TSize>(0, length(value(seqSet,0))));
	} else {
		for(TSize l = 0;l < lastLetter; ++l) {
			appendValue(ctgSt.data_contig, '-'); // Dummy letter
			appendValue(ctgSt.data_quality, char(48+60)); // Dummy quality value
		}
		appendValue(ctgSt.data_begin_end, Pair<TSize, TSize>(0, lastLetter));
	}
	++ctgSt.data_pos_count; // Added just one contig



	// Convert Libraries
	clear(filePath);
	filePath = value(cfgOpt, "reads");
	appendValue(filePath, 'L');
	StringSet<String<char>, Owner<> > libSet;
	String<TName> libIds;
	count = _loadSequences(filePath, libSet, libIds);
	if (count != 0) {
		for(TSize i = 0; i<length(libSet); ++i) {
			// Create the library name
			std::stringstream input;
			input << "L" << i;
			String<char> tmp(input.str().c_str());
			appendValue(libSt.data_names, tmp);

			// Get the mean and std
			typedef typename Iterator<String<char> >::Type TStrIter;
			TStrIter nameIt = begin(value(libSet, i));
			TStrIter nameItEnd = end(value(libSet, i));
			bool before = true;
			String<char> inf1; String<char> inf2;
			for(;nameIt != nameItEnd; goNext(nameIt)) {
				if (value(nameIt) == ',') before = false;
				else {
					if (before) appendValue(inf1, value(nameIt));
					else appendValue(inf2, value(nameIt));
				}
			}
			TSize mean = _stringToNumber<TSize>(inf1);
			TSize std = _stringToNumber<TSize>(inf2);
			appendValue(libSt.data_mean, mean);
			appendValue(libSt.data_std, std);
			++libSt.data_pos_count;
		}
	} else {
		appendValue(libSt.data_names, "L0");
		appendValue(libSt.data_mean, 0);
		appendValue(libSt.data_std, 0);
		++libSt.data_pos_count;
	}



	// Convert Fragments
	filePath = value(cfgOpt, "reads");
	appendValue(filePath, 'F');
	StringSet<String<char>, Owner<> > fragSet;
	String<TName> fragIds;
	count = _loadSequences(filePath, fragSet, fragIds);
	if (count != 0) {
		for(TSize i = 0;i<length(fragSet); ++i) {
			// Create the fragment name
			std::stringstream input;
			input << "F" << i;
			String<char> tmp(input.str().c_str());
			appendValue(frgSt.data_names, tmp);
			
			// Find the library id
			TSize libId = i;
			TName valKey;
			__findThis(value(fragIds, i), "libId=", valKey);
			if (!empty(valKey)) libId = _stringToNumber<TSize>(valKey);
			appendValue(frgSt.data_lib_id, libId);
			
			// Get rds1 and rds2
			typedef typename Iterator<String<char> >::Type TStrIter;
			TStrIter nameIt = begin(value(fragSet, i));
			TStrIter nameItEnd = end(value(fragSet, i));
			bool before = true;
			String<char> inf1; String<char> inf2;
			for(;nameIt != nameItEnd; goNext(nameIt)) {
				if (value(nameIt) == ',') before = false;
				else {
					if (before) appendValue(inf1, value(nameIt));
					else appendValue(inf2, value(nameIt));
				}
			}
			TSize rds1 = _stringToNumber<TSize>(inf1);
			TSize rds2 = _stringToNumber<TSize>(inf2);
			if (rds1 < rds2) appendValue(frgSt.data_rds, Pair<TSize, TSize>(rds1, rds2));
			else appendValue(frgSt.data_rds, Pair<TSize, TSize>(rds2, rds1));
			++frgSt.data_pos_count;
		}
	} else {
		for(TSize i = 0; i<readSt.data_pos_count; ++i) {
			std::stringstream input;
			input << "F" << i;
			String<char> tmp(input.str().c_str());
			appendValue(frgSt.data_names, tmp);
			appendValue(frgSt.data_lib_id, 0);
			appendValue(frgSt.data_rds, Pair<TSize, TSize>(0, 0));
			++frgSt.data_pos_count;
		}
	}



	// Find smallest offset and reset all contig offsets
	if (moveToFront) {
		for(TSize i = 0; i<length(ctgSt); ++i) {
			TSize smallestOffset = length(ctgSt.data_contig);
			String<GappedRead<> >& gappedReads = value(ctgSt.data_reads, i);
			for(TSize k = 0; k<length(gappedReads); ++k) {
				if ((value(gappedReads, k)).data_offset < smallestOffset) smallestOffset = (value(gappedReads, k)).data_offset;
			}
			if (smallestOffset > 0) {
				for(TSize k = 0; k<length(gappedReads); ++k) {
					(value(gappedReads, k)).data_offset -= smallestOffset;
				}
			}
		}
	}
	return true;
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
	TKey keys[] = {"afg", "reads", "bandwidth", "matchlength", "quality", "overlaps", "window", "snp", "evaluation","outfile", "output", "convert", "moveToFront", "source"};
	assignKeys(cfgOpt, keys, 14);
	assign(cfgOpt, "moveToFront", "false");
	assign(cfgOpt, "matchlength", "15");
	assign(cfgOpt, "quality", "80");
#ifdef CELERA_OFFSET
	assign(cfgOpt, "bandwidth", "15");
	assign(cfgOpt, "overlaps", "5");
#else
	assign(cfgOpt, "bandwidth", "8");
	assign(cfgOpt, "overlaps", "3");
#endif
	assign(cfgOpt, "window", "0");
	assign(cfgOpt, "call", "majority");
	assign(cfgOpt, "output", "seqan");
	assign(cfgOpt, "outfile", "readAlign.txt");

	// Help Message
	String<char> helpMsg;
	append(helpMsg, "Usage: ./seqcons -reads <FASTA File with Reads> [Options]\n");
	append(helpMsg, "\nOptions\n");

	// Main options
	append(helpMsg, "\nMain Options\n------------\n");
	append(helpMsg, "-reads <FASTA File with Reads>\n");
	append(helpMsg, "\tReads in FASTA format and approximate layout positions.\n\n");
	append(helpMsg, "-afg <AMOS message file>\n");
	append(helpMsg, "\tInput is read from an AMOS message file instead of a fasta file.\n\n");
	append(helpMsg, "-matchlength <Length>\n");
	append(helpMsg, "\tMinimum match-length for an overlap, default is 15.\n\n");
	append(helpMsg, "-quality <Number>\n");
	append(helpMsg, "\tMinimum quality of an overlap, default is 80 (for 80% identity).\n\n");
	append(helpMsg, "-overlaps <Number>\n");
	append(helpMsg, "\tNumber of overlaps that are computed per read, default is 3.\n\n");
	append(helpMsg, "-bandwidth <Number>\n");
	append(helpMsg, "\tSpecifies the bandwidth, default is 8.\n\n");
	append(helpMsg, "-call [majority | bayesian]\n");
	append(helpMsg, "\tHow to call consensus bases, default is majority.\n\n");
	append(helpMsg, "-outfile <Alignment Filename>\n");
	append(helpMsg, "\tThe name of the output file, default is readAlign.txt.\n\n");
	append(helpMsg, "-output [seqan | afg]\n");
	append(helpMsg, "\tThe output format, default is seqan.\n\n");

	// Insert Sequencing
	append(helpMsg, "\nInsert Sequencing\n------------\n");	
	append(helpMsg, "-window <Window-Size>\n");
	append(helpMsg, "\tIf this parameter is given all overlaps within a given window are computed (no banded alignment).\n\n");

	// Examples
	append(helpMsg, "\n\n\nExamples\n");
	append(helpMsg, "\nMulti-Read Alignment:\n");
	append(helpMsg, "\t./seqcons -reads reads.fasta\n");
	append(helpMsg, "\t./seqcons -reads reads.fasta -matchlength 15 -quality 80 -outfile align.txt\n");
	append(helpMsg, "\nInsert Sequencing (Reads of length 35):\n");
	append(helpMsg, "\t./seqcons -reads reads.fasta -window 200 -overlaps 100\n");
	append(helpMsg, "\nInsert Sequencing (Reads of length 200):\n");
	append(helpMsg, "\t./seqcons -reads reads.fasta -window 300 -overlaps 100\n");
	append(helpMsg, "\nInsert Sequencing (Reads of length 800):\n");
	append(helpMsg, "\t./seqcons -reads reads.fasta -window 1000 -overlaps 100\n");
	assignHelp(cfgOpt, helpMsg);

	std::cout << "**************************************************" << std::endl;
	std::cout << "* Consensus Computation                          *" << std::endl;
	std::cout << "*                                                *" << std::endl;
	std::cout << "* SeqCons                                        *" << std::endl;
	std::cout << "* Version: 0.202 (17. October 2008)              *" << std::endl;
	std::cout << "**************************************************" << std::endl;
	std::cout << std::endl;

	if (argc < 2) {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}
	if (!parseCmdLine(argc, argv, cfgOpt)) return -1;
	if (!empty(value(cfgOpt, "evaluation"))) {
		//evaluationOfConsensusAlignment(cfgOpt);
		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Read simulation file, afg file or celera consensus file
	//////////////////////////////////////////////////////////////////////////////
	ReadStore<> readSt;
	FrgStore<> frgSt;
	LibStore<> libSt;
	CtgStore<> ctgSt;

	TSize numberOfContigs = 0;
	if (!empty(value(cfgOpt, "reads"))) {
		bool success;
		if (value(cfgOpt, "moveToFront") == "false") success = convertSimpleReadFile(readSt, frgSt, libSt, ctgSt, cfgOpt, false);
		else success = convertSimpleReadFile(readSt, frgSt, libSt, ctgSt, cfgOpt, true);
		if (!success) { std::cerr << valueHelp(cfgOpt) << std::endl; return -1; }
		numberOfContigs = 1;
	} else if (!empty(value(cfgOpt, "afg"))) {
		// Read Amos
		std::fstream strmReads;
		strmReads.open(toCString(value(cfgOpt, "afg")), std::ios_base::in | std::ios_base::binary);
		read(strmReads,readSt,frgSt,libSt,ctgSt,Amos());	
		strmReads.close();
		numberOfContigs = length(ctgSt);
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

	// Initialize alignment parameters
	TSize matchlength = _stringToNumber<TSize>(value(cfgOpt, "matchlength"));
	TSize quality = _stringToNumber<TSize>(value(cfgOpt, "quality"));
	TSize bandwidth = _stringToNumber<TSize>(value(cfgOpt, "bandwidth"));
	TSize overlaps = _stringToNumber<TSize>(value(cfgOpt, "overlaps"));
	TSize window = _stringToNumber<TSize>(value(cfgOpt, "window"));

	//////////////////////////////////////////////////////////////////////////////
	// Iterate over all contigs
	//////////////////////////////////////////////////////////////////////////////	

	for(TSize currentContig = 0; currentContig < numberOfContigs; ++currentContig) {
	
		// Profiling
		SEQAN_PROTIMESTART(profileTime);


		//////////////////////////////////////////////////////////////////////////////
		// Import all reads of the given contig
		//////////////////////////////////////////////////////////////////////////////
		
		typedef String<Dna> TSequence;
		typedef Id<TSequence>::Type TId;
		StringSet<TSequence, Owner<> > origStrSet;
		String<Pair<TSize, TSize> > begEndPos;
		loadReadsClr(readSt, ctgSt, currentContig, origStrSet, begEndPos);
		TSize nseq = length(origStrSet);
		if (nseq == 0) continue;
		std::cout << "Import sequences done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;

		// Make a dependent StringSet
		typedef StringSet<TSequence, Dependent<> > TStringSet;
		TStringSet seqSet = origStrSet;

		//////////////////////////////////////////////////////////////////////////////
		// Align the reads
		//////////////////////////////////////////////////////////////////////////////

		// Get the average read length and estimate a bandwidth
		TSize avgReadLength = 0;
		for(TSize i = 0; i<length(origStrSet);++i) avgReadLength += length(value(origStrSet,i));
		avgReadLength /= nseq;
		//if (avgReadLength < 50) overlaps += 3;
						
		// Select all overlapping reads and record the diagonals of the band
		String<Pair<TId, TId> > pList;
		String<Pair<int, int> > diagList;
		if (window == 0) selectPairs(seqSet, begEndPos, bandwidth, pList, diagList);
		else selectPairsIndel(seqSet, begEndPos, window, pList, diagList);

		// Estimate the number of overlaps we want to compute
		if (window == 0) std::cout << "Matchlength: " << matchlength << ", " << "Quality: " << quality << ", " << "Bandwidth: " << bandwidth << ", " << "Overlaps: " << overlaps << std::endl;
		else std::cout << "Matchlength: " << matchlength << ", " << "Quality: " << quality << ", " << "Window: " << window << ", " << "Overlaps: " << overlaps << std::endl;
		std::cout << "Number of reads: " << nseq << std::endl;
		std::cout << "Average read length: " << avgReadLength << std::endl;
		if (window == 0) {
			TSize covEstim = length(pList) / nseq;
			std::cout << "Estimated coverage: " << covEstim << std::endl;
		}
		std::cout << "Pair selection done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;

		// Set-up a sparse distance matrix
		Graph<Undirected<double> > pairGraph;
		
		// Set-up alignment scoring matrices
		typedef Score<int> TScore;
		typedef Value<TScore>::Type TScoreValue;
		TScore scType = Score<int>(2,-6,-4,-9);
		//TScore scType = Score<int>(1,-2,-1,-2);

		// Containers for segment matches and corresponding scores 
		typedef String<Fragment<> > TFragmentString;
		TFragmentString matches;
		typedef String<TScoreValue> TScoreValues;
		TScoreValues scores;

		// Compute segment matches from global pairwise alignments
		appendSegmentMatches(seqSet, pList, diagList, begEndPos, scType, matchlength, quality, overlaps, matches, scores, pairGraph, Overlap_Library() );
		clear(pList);
		clear(diagList);
		std::cout << "Overlap done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;
	
		// Re-Score the matches
		scoreMatches(seqSet, scType, matches, scores);
		std::cout << "Re-scoring done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;

		// Use these segment matches for the initial alignment graph
		typedef Graph<Alignment<TStringSet, TSize> > TGraph;
		TGraph g(seqSet);
		buildAlignmentGraph(matches, scores, g, FractionalScore() );
		clear(matches);
		clear(scores);
		std::cout << "Construction of Alignment Graph done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;
		
		// Triplet library extension
		if ( ((2 * numEdges(g)) / numVertices(g) ) < 50 ) tripletLibraryExtension(g);
		else reducedTripletLibraryExtension(g);
		std::cout << "Triplet done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;
	
		// Guide Tree
		Graph<Tree<double> > guideTree;
		upgmaTree(pairGraph, guideTree);
		std::cout << "Guide tree done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;
		clear(pairGraph);

		// Perform a progressive alignment
		Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(seqSet);
		progressiveAlignment(g, guideTree, gOut);
		clear(g);
		clear(guideTree);
		std::cout << "Progressive alignment done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;
	
		// Build the read alignment matrix
		TSize alignDepth;
		String<char> alignmentMatrix;
		String<Triple<unsigned int, unsigned int, unsigned int> > readBegEndRowPos;
		multireadAlignment(gOut, alignmentMatrix, readBegEndRowPos, alignDepth);
		clear(gOut);
		std::cout << "Multi-read Alignment done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;

		// Find disrupted reads
		//TSize badReads = fixDisruptedReads(alignmentMatrix, seqSet, scType, readBegEndRowPos, alignDepth);
		//std::cout << "Bad reads: " << badReads << std::endl;
		
		// Call the consensus
		String<unsigned int> coverage;
		String<char> gappedConsensus;
		String<Dna> consensusSequence;
		if (value(cfgOpt, "call") == "majority") consensusCalling(alignmentMatrix, consensusSequence, gappedConsensus, coverage, alignDepth, Majority_Vote() );
		else consensusCalling(alignmentMatrix, consensusSequence, gappedConsensus, coverage, alignDepth, Bayesian() );
		std::cout << "Consensus done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;

		//////////////////////////////////////////////////////////////////////////////
		// Output of aligned reads
		//////////////////////////////////////////////////////////////////////////////

		if (value(cfgOpt, "output") == "seqan") {
			std::fstream strm;
			if (currentContig == 0) strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
			else strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::app);
			write(strm, seqSet, alignmentMatrix, begEndPos, readBegEndRowPos, gappedConsensus, coverage, FastaReadFormat());
			strm.close();

			//// Debug code for CA
			//mtRandInit();
			//String<char> fileTmp1 = "tmp1";
			//String<char> fileTmp2 = "tmp2";
			//for(int i = 0; i<10; ++i) {
			//	int file = (mtRand() % 20) + 65;
			//	appendValue(fileTmp1, char(file));
			//	appendValue(fileTmp2, char(file));
			//}
			//std::fstream strm3;
			//strm3.open(toCString(fileTmp2), std::ios_base::out | std::ios_base::trunc);
			//for(int i = 0;i<length(origStrSet); ++i) {
			//	std::stringstream name;
			//	name << value(begEndPos, i).i1 << "," << value(begEndPos, i).i2;
			//	String<char> myTitle = name.str();
			//	write(strm3, origStrSet[i], myTitle, Fasta());			
			//	if (value(begEndPos, i).i1 > value(begEndPos, i).i2) reverseComplementInPlace(origStrSet[i]);
			//}
			//strm3.close();
			//std::fstream strm2;
			//strm2.open(toCString(fileTmp1), std::ios_base::out | std::ios_base::trunc);
			//write(strm2, seqSet, alignmentMatrix, begEndPos, readBegEndRowPos, gappedConsensus, coverage, FastaReadFormat());
			//strm2.close();
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
		}

		std::cout << "Output done: " << SEQAN_PROTIMEUPDATE(profileTime) << " seconds" << std::endl;
	}
	
	// Write the AMOS message file
	if (value(cfgOpt, "output") == "afg") {
		// Write Amos
		std::fstream strmWrite;
		strmWrite.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strmWrite,readSt,frgSt,libSt,ctgSt,Amos());	
		strmWrite.close();
	}

	return 0;
}
