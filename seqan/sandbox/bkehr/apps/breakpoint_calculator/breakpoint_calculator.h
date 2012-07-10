// ==========================================================================
//                           breakpoint_calculator
// ==========================================================================
// Copyright (C) 2012 by Birte Kehr
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================


#ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_CALCULATOR_H_
#define SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_CALCULATOR_H_

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>

#include "parse_alignment.h"
#include "breakpoint_counts.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

enum AlignmentFormat
{
	XMFA,
	MAF
};

struct Options
{
	bool verbose;
	bool detailed;

	bool pairwiseCount;
	bool tripletCount;

	CharString inputFile;
	AlignmentFormat inputFormat;
	bool swapPositionsXmfa;

    Options()
    {
        // Set defaults.
		verbose = false;
		detailed = false;

		pairwiseCount = false;
		tripletCount = false;

		inputFormat = XMFA;
		swapPositionsXmfa = false;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

void
setupCommandLineParser(ArgumentParser & parser)
{
    setShortDescription(parser, "for multiple genome alignments");
    setVersion(parser, "0.2");
    setDate(parser, "Jul 2012");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIALIGNMENTFILE\\fP");
    
    addDescription(parser, "Can calculate pairwise and hidden threeway breakpoints in multiple genome alignments. The alignment must be given in MAF or XMFA format. If the format is not directly specified, it is guessed from the file extension.");
    addDescription(parser, "(c) 2012 by Birte Kehr");
    
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE, false, "IALIGNMENTFILE"));
    //setValidValues(parser, 0, "xmfa maf"); // allow only *.xmfa and *.maf files as input

	addSection(parser, "Main Options");
	addOption(parser, ArgParseOption("d2", "pairwiseCount", "Compute pairwise breakpoint counts."));
	addOption(parser, ArgParseOption("d3", "tripletCount", "Compute triplet breakpoint counts."));
	addOption(parser, ArgParseOption("d", "detailed", "Print breakpoint counts of all pairs/triplets."));

    addSection(parser, "Miscellaneous");
	addOption(parser, ArgParseOption("f", "inFormat", "Format of input file (xmfa or maf).", ArgParseArgument::STRING));
    setValidValues(parser, "f", "xmfa maf");
    addOption(parser, ArgParseOption("v", "verbose", "Turn on verbose output."));
    addOption(parser, ArgParseOption("x", "swapPositions", "Turn on swapping of start and end position for reverse orientation in XMFA format (necessary for sgEvolver output)."));
	hideOption(parser, "x");

    addTextSection(parser, "References");
    addText(parser, "Kehr, B., Reinert, K., Darling, A.: Hidden breakpoints in genome alignments. WABI 2012. To be presented.");
}

ArgumentParser::ParseResult
parseArgumentsAndCheck(Options & options,
                       ArgumentParser & parser,
                       int argc,
                       char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    
	if (res == ArgumentParser::PARSE_OK)
	{
		getArgumentValue(options.inputFile, parser, 0);
		getOptionValue(options.verbose, parser, "verbose");
		getOptionValue(options.detailed, parser, "detailed");
		getOptionValue(options.pairwiseCount, parser, "pairwiseCount");
		getOptionValue(options.tripletCount, parser, "tripletCount");
		getOptionValue(options.swapPositionsXmfa, parser, "swapPositions");

		CharString format;
		if (isSet(parser, "f"))
			getOptionValue(format, parser, "inFormat");
		else
			format = suffix(options.inputFile, length(options.inputFile) - 3);

		toUpper(format);
		if (format == "MAF")
			options.inputFormat = MAF;
		else
			options.inputFormat = XMFA;
	}

	return res;
}

int mainWithOptions(Options & options)
{
	typedef Dna5 TAlphabet;

	typedef String<TAlphabet> TSequence;
	typedef Size<TSequence>::Type TSize;

	typedef StringSet<TSequence> TStringSet;

	typedef Align<TSequence, ArrayGaps> TAlign;
	typedef std::map<CharString, AlignmentBlockRow<TSize, TSize> > TIdRowMap;

	if (options.verbose)
	{
		std::cout << "Alignment file: " << options.inputFile << std::endl;
		std::cout << std::endl;
	}

	// Open input file.
	std::fstream inStream(toCString(options.inputFile), std::ios::in | std::ios::binary);
	if (!inStream.good())
	{
		std::cerr << "ERROR: Could not open " << options.inputFile << "!\n";
		return 1;
	}

	String<TAlign> aligns;
	String<TIdRowMap> idToRowMaps;
	TStringSet seqs;

	// Parse the input alignment file
	if (options.inputFormat == MAF && parseAlignment(aligns, idToRowMaps, seqs, inStream, options.verbose, Maf())) 
	{
		std::cerr << "ERROR: Parsing file in MAF format failed for " << options.inputFile << "!" << std::endl;
		return 1;
	}
	else if (options.inputFormat == XMFA && options.swapPositionsXmfa && parseAlignment(aligns, idToRowMaps, seqs, inStream, options.verbose, XmfaSwap())) 
	{
		std::cerr << "ERROR: Parsing file in XMFA format failed for " << options.inputFile << "!" << std::endl;
		return 1;
	}
	else if (options.inputFormat == XMFA && !options.swapPositionsXmfa && parseAlignment(aligns, idToRowMaps, seqs, inStream, options.verbose, Xmfa())) 
	{
		std::cerr << "ERROR: Parsing file in XMFA format failed for " << options.inputFile << "!" << std::endl;
		return 1;
	}
	inStream.close();
	
	// Compute the breakpoint counts
	if (options.pairwiseCount || options.tripletCount)
	{
		typedef int TBlockId;

		std::map<CharString, String<TBlockId> > blockSeqs;
		sequencesOfBlocks(blockSeqs, idToRowMaps);

		// // Debug code:
		// for (std::map<CharString, String<TBlockId> >::const_iterator it = blockSeqs.begin(); it != blockSeqs.end(); ++it)
		// {
		// 	std::cout << it->first << ": ";
		// 	for (Iterator<String<TBlockId> >::Type itit = begin(it->second); itit != end(it->second); ++itit)
		// 		std::cout << *itit << "  ";
		// 	std::cout << std::endl;
		//}

		std::map<CharString, StringSet<String<TBlockId>, Dependent<> > > blockSeqSets;
		collateChromosomes(blockSeqSets, blockSeqs);

		if (options.pairwiseCount) {
			if (options.verbose) std::cout << "Computing pairwise count..." << std::endl;
			std::cout << pairwiseCounts(blockSeqSets, options.detailed) << " pairwise breakpoints" << std::endl;
		}
		if (options.tripletCount)
		{
			if (options.verbose) std::cout << "Computing triplet count..." << std::endl;
			std::cout << tripletCounts(blockSeqSets, options.detailed) << " 3-way breakpoints" << std::endl;
		}
	}

	return 0;
}

#endif  // #ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_CALCULATOR_H_
