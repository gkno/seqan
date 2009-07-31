/*==========================================================================
               SeqAn - The Library for Sequence Analysis
                         http://www.seqan.de 
============================================================================
Copyright (C) 2007

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
==========================================================================*/

#define SEQAN_PROFILE

#include <seqan/basic.h>

// Profiling
#ifdef SEQAN_PROFILE
		SEQAN_PROTIMESTART(__myProfileTime); 
#endif

#include <seqan/consensus.h>
#include <seqan/modifier.h>
#include <seqan/misc/misc_cmdparser.h>

#include <iostream>
#include <fstream>


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

inline void
_addVersion(CommandLineParser& parser) {
	::std::string rev = "$Revision: 4637 $";
	addVersionLine(parser, "Version 0.21 (31. July 2009) Revision: " + rev.substr(11, 4) + "");
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TStrSpec, typename TPosPair, typename TStringSpec, typename TSpec, typename TConfig, typename TId>
inline void 
getContigReads(StringSet<TValue, Owner<TStrSpec> >& strSet,
			   String<TPosPair, TStringSpec>& startEndPos,
			   FragmentStore<TSpec, TConfig> const& fragStore,
			   TId const contigId)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// Sort aligned reads according to contig id
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	resize(strSet, length(fragStore.alignedReadStore));

	// Retrieve all reads, limit them to the clear range and if required reverse complement them
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TSize numRead = 0;
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		TSize offset = _min(alignIt->beginPos, alignIt->endPos);
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		value(strSet, numRead) = infix(fragStore.readSeqStore[alignIt->readId], begClr, endClr);
		TSize lenRead = length(value(strSet, numRead));
		if (alignIt->beginPos < alignIt->endPos) {
			appendValue(startEndPos, TPosPair(offset, offset + lenRead), Generous());
		} else {
			reverseComplementInPlace(value(strSet, numRead));
			appendValue(startEndPos, TPosPair(offset + lenRead, offset), Generous());
		}
		++numRead;
	}
	resize(strSet, numRead, Exact());
}



//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {
	// Command line parsing
	CommandLineParser parser;
	_addVersion(parser);

	addTitleLine(parser, "***************************************");
	addTitleLine(parser, "* Multi-read alignment - SeqCons      *");
	addTitleLine(parser, "* (c) Copyright 2009 by Tobias Rausch *");
	addTitleLine(parser, "***************************************");

	addUsageLine(parser, "-r <FASTA file with reads> [Options]");
	addUsageLine(parser, "-a <AMOS message file> [Options]");

	addSection(parser, "Main Options:");
	addOption(parser, addArgumentText(CommandLineOption("r", "reads", "file with reads", OptionType::String), "<FASTA reads file>"));
	addOption(parser, addArgumentText(CommandLineOption("a", "afg", "message file", OptionType::String), "<AMOS afg file>"));
	addOption(parser, addArgumentText(CommandLineOption("o", "outfile", "output filename", OptionType::String, "align.txt"), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("f", "format", "output format", OptionType::String, "seqan"), "[seqan | afg]"));
	addOption(parser, addArgumentText(CommandLineOption("m", "method", "multi-read alignment method", OptionType::String, "realign"), "[realign | msa]"));
	addOption(parser, addArgumentText(CommandLineOption("b", "bandwidth", "bandwidth", OptionType::Int, 8), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("c", "consensus", "consensus calling", OptionType::String, "majority"), "[majority | bayesian]"));
	addOption(parser, CommandLineOption("n", "noalign", "no align, only convert input", OptionType::Boolean));

	addSection(parser, "MSA Method Options:");
	addOption(parser, addArgumentText(CommandLineOption("ma", "matchlength", "minimum overlap length", OptionType::Int, 15), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("qu", "quality", "minimum overlap precent identity", OptionType::Int, 80), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("ov", "overlaps", "minimum number of overlaps per read", OptionType::Int, 3), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("wi", "window", "window size", OptionType::Int, 0), "<Int>"));
	addHelpLine(parser, "/*If this parameter is > 0 then all");
	addHelpLine(parser, "  overlaps within a given window");
	addHelpLine(parser, "  are computed.*/");

	addSection(parser, "ReAlign Method Options:");
	addOption(parser, CommandLineOption("in", "include", "include contig sequence", OptionType::Boolean));
	addOption(parser, addArgumentText(CommandLineOption("rm", "rmethod", "realign method", OptionType::String, "nw"), "[nw | gotoh]"));

	if (argc == 1)
	{
		shortHelp(parser, std::cerr);	// print short help and exit
		return 0;
	}

	if (!parse(parser, argc, argv, ::std::cerr)) return 1;
	if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit


	// Get all command line options
	ConsensusOptions consOpt;

	// Main options
	getOptionValueLong(parser, "reads", consOpt.readsfile);
	getOptionValueLong(parser, "afg", consOpt.afgfile);
	getOptionValueLong(parser, "outfile", consOpt.outfile);
	String<char> optionVal;
	getOptionValueLong(parser, "format", optionVal);
	if (optionVal == "seqan") consOpt.output = 0;
	else if (optionVal == "afg") consOpt.output = 1;
	else if (optionVal == "frg") consOpt.output = 2;
	else if (optionVal == "cgb") consOpt.output = 3;
	getOptionValueLong(parser, "method", optionVal);
	if (optionVal == "realign") consOpt.method = 0;
	else if (optionVal == "msa") consOpt.method = 1;
	getOptionValueLong(parser, "bandwidth", consOpt.bandwidth);
#ifdef CELERA_OFFSET
	if (!isSetLong(parser, "bandwidth") consOpt.bandwidth = 15;	
#endif
	getOptionValueLong(parser, "consensus", optionVal);
	if (optionVal == "majority") consOpt.consensus = 0;
	else if (optionVal == "bayesian") consOpt.consensus = 1;
	getOptionValueLong(parser, "noalign", consOpt.noalign);

	// Msa options
	getOptionValueLong(parser, "matchlength", consOpt.matchlength);
	getOptionValueLong(parser, "quality", consOpt.quality);
	getOptionValueLong(parser, "overlaps", consOpt.overlaps);
#ifdef CELERA_OFFSET
	if (!isSetLong(parser, "overlaps") consOpt.overlaps = 5;	
#endif
	getOptionValueLong(parser, "window", consOpt.window);
	
	// ReAlign options
	getOptionValueLong(parser, "include", consOpt.include);
	getOptionValueLong(parser, "rmethod", optionVal);
	if (optionVal == "nw") consOpt.rmethod = 0;
	else if (optionVal == "gotoh") consOpt.rmethod = 1;


	// Create a new fragment store
	typedef FragmentStore<> TFragmentStore;
	typedef Size<TFragmentStore>::Type TSize;
	TFragmentStore fragStore;

	// Load the reads and layout positions
	TSize numberOfContigs = 0;
	if (!empty(consOpt.readsfile)) {
		// Load simple read file
		FILE* strmReads = fopen(consOpt.readsfile.c_str(), "rb");
		bool success = _convertSimpleReadFile(strmReads, fragStore, consOpt.readsfile, true);
		fclose(strmReads);
		if (!success) { 
			shortHelp(parser, std::cerr);
			return 0;
		}
		numberOfContigs = 1;
	} else if (!empty(consOpt.afgfile)) {
		// Load Amos message file
		FILE* strmReads = fopen(consOpt.afgfile.c_str(), "rb");
		read(strmReads, fragStore, Amos());	
		fclose(strmReads);
		numberOfContigs = length(fragStore.contigStore);
	} else {
		shortHelp(parser, std::cerr);
		return 0;
	}



	// Multi-realignment desired or just conversion of the input
	if (!consOpt.noalign) {
	
		// Profiling
#ifdef SEQAN_PROFILE
				SEQAN_PROTIMEUPDATE(__myProfileTime); 
#endif
	
		// Iterate over all contigs
		for(TSize currentContig = 0; currentContig < numberOfContigs; ++currentContig) {

			if (consOpt.method == 0) {
#ifdef SEQAN_PROFILE
				::std::cout << "ReAlign method" << ::std::endl;
				if (consOpt.rmethod == 0) ::std::cout << "Realign algorithm: Needleman-Wunsch" << ::std::endl;
				else if (consOpt.rmethod == 1) ::std::cout << "Realign algorithm: Gotoh" << ::std::endl;
				::std::cout << "Bandwidth: " << consOpt.bandwidth << ::std::endl;
				::std::cout << "Include reference: " << consOpt.include << ::std::endl;
#endif
				Score<int, WeightedConsensusScore<Score<int, FractionalScore>, Score<int, ConsensusScore> > > combinedScore;
				reAlign(fragStore, combinedScore, currentContig, consOpt.rmethod, consOpt.bandwidth, consOpt.include);

#ifdef SEQAN_PROFILE
				::std::cout << "ReAlignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << ::std::endl;
#endif
			} else {

#ifdef SEQAN_PROFILE
				::std::cout << "MSA method" << ::std::endl;
				::std::cout << "Bandwidth: " << consOpt.bandwidth << ::std::endl;
				::std::cout << "Matchlength: " << consOpt.matchlength << ::std::endl;
				::std::cout << "Quality: " << consOpt.quality << ::std::endl;
				::std::cout << "Window: " << consOpt.window << ::std::endl;
#endif

				// Import all reads of the given contig
				typedef TFragmentStore::TReadSeq TReadSeq;
				typedef Id<TFragmentStore>::Type TId;
				StringSet<TReadSeq, Owner<> > readSet;
				String<Pair<TSize, TSize> > begEndPos;
	
				getContigReads(readSet, begEndPos, fragStore, currentContig);
				TSize nseq = length(readSet);
				if (nseq == 0) continue;

#ifdef SEQAN_PROFILE
				::std::cout << "Import sequences done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << ::std::endl;
#endif

				// Align the reads
				Graph<Alignment<StringSet<TReadSeq, Dependent<> >, void, WithoutEdgeId> > gOut(readSet);
				consensusAlignment(gOut, begEndPos, consOpt);
#ifdef SEQAN_PROFILE
				std::cout << "Multi-read Alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

			
				// Update the contig in the fragment store
				updateContig(fragStore, gOut, currentContig);
				clear(gOut);
#ifdef SEQAN_PROFILE
				std::cout << "Update contig done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

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
				//for(int i = 0;i<(int) length(origStrSet); ++i) {
				//	std::stringstream name;
				//	name << value(begEndPos, i).i1 << "," << value(begEndPos, i).i2;
				//	String<char> myTitle = name.str();
				//	write(strm3, origStrSet[i], myTitle, Fasta());			
				//	if (value(begEndPos, i).i1 > value(begEndPos, i).i2) reverseComplementInPlace(origStrSet[i]);
				//}
				//strm3.close();
			}

			// ToDo: Consensus calling methods go here
			//String<unsigned int> coverage;
			//		String<char> gappedConsensus;
			//		String<Dna> consensusSequence;
			//		if (consOpt.snp == 0) consensusCalling(alignmentMatrix, consensusSequence, gappedConsensus, coverage, alignDepth, Majority_Vote() );
			//		else consensusCalling(alignmentMatrix, consensusSequence, gappedConsensus, coverage, alignDepth, Bayesian() );
		
		} // end loop over all contigs
	}
	
	// Output
	if (consOpt.output == 0) {
		// Write old SeqAn multi-read alignment format
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		write(strmWrite, fragStore, FastaReadFormat());	
		fclose(strmWrite);
	} else if (consOpt.output == 1) {
		// Write Amos
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		write(strmWrite, fragStore, Amos());	
		fclose(strmWrite);
	} else if (consOpt.output == 2) {
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		_writeCeleraFrg(strmWrite, fragStore);	
		fclose(strmWrite);
	} else if (consOpt.output == 3) {
		FILE* strmWrite = fopen(consOpt.outfile.c_str(), "w");
		_writeCeleraCgb(strmWrite, fragStore);	
		fclose(strmWrite);
	}

	return 0;
}
