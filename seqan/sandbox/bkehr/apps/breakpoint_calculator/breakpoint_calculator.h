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
#include <seqan/refinement.h>
#include <seqan/stream.h>
#include <seqan/score.h>

#include <seqan/misc/misc_cmdparser.h>

#include "parse_alignment.h"
#include "score_blocks.h"
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
    bool showHelp;
    bool showVersion;
	bool verbose;
	bool detailed;

	bool pairwiseCount;
	bool tripletCount;
	bool sopScore;

	CharString inputFile;
	AlignmentFormat inputFormat;
	bool swapPositionsXmfa;
    
	CharString matrixFile;
	bool hoxdMatrix;
	bool computeMatrix;
	int gapOpen;
	int gapExtend;


    Options()
    {
        // Set defaults.
        showHelp = false;
        showVersion = false;
		verbose = false;
		detailed = false;

		pairwiseCount = false;
		tripletCount = false;
		sopScore = false;

		hoxdMatrix = false;
		computeMatrix = true;
		gapOpen = -4;
		gapExtend = -1;

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
setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
    addVersionLine(parser, "0.1");
    
    addTitleLine(parser, "****************************");
    addTitleLine(parser, "*  breakpoint_calculator   *");
    addTitleLine(parser, "****************************");
    addTitleLine(parser, "");
    addTitleLine(parser, "(c) 2012 by Birte Kehr");

    addUsageLine(parser, "[OPTIONS] <ALIGNMENT FILE>");
    
	addSection(parser, "Main Options");
	addOption(parser, CommandLineOption("f", "inFormat", "input file format: XMFA, MAF", OptionType::String | OptionType::Label, "XMFA"));
	addOption(parser, CommandLineOption("x", "swapPositions", "swap start and end position for reverse orientation", OptionType::Bool | OptionType::Label));
	addHelpLine(parser, "in XMFA format (necessary for sgEvolver output) (default false)");
	addOption(parser, CommandLineOption("v", "verbose", "verbosity mode", OptionType::Bool | OptionType::Label, "false"));
	addOption(parser, CommandLineOption("d", "detailed", "print breakpoint counts of all pairs and triplets", OptionType::Bool | OptionType::Label, "false"));

	addOption(parser, CommandLineOption("d2", "pairwiseCount", "compute pairwise breakpoint counts", OptionType::Bool | OptionType::Label, "false"));
	addOption(parser, CommandLineOption("d3", "tripletCount", "compute triplet breakpoint counts", OptionType::Bool | OptionType::Label, "false"));
	addOption(parser, CommandLineOption("s", "score", "compute sum-of-pairs score for alignment blocks", OptionType::Bool | OptionType::Label, "false"));

	addSection(parser, "Scoring Options");
	addOption(parser, CommandLineOption("m", "scoreMatrix", "file containing a DNA scoring matrix OR 'compute'", OptionType::String | OptionType::Label, "compute"));
	addOption(parser, CommandLineOption("go", "gapOpen", "gap open penalty", OptionType::Int | OptionType::Label, options.gapOpen));
	addOption(parser, CommandLineOption("ge", "gapExtend", "gap extension penalty", OptionType::Int | OptionType::Label, options.gapExtend));
    
    requiredArguments(parser, 1);
}

int parseCommandLineAndCheck(Options & options,
                             CommandLineParser & parser,
                             int argc,
                             char const ** argv)
{
    bool stop = !parse(parser, argc, argv);
    if (stop)
        return 1;

    options.inputFile = getArgumentValue(parser, 0);

    if (isSetLong(parser, "help")) {
        options.showHelp = true;
        return 0;
    }
    if (isSetLong(parser, "version")) {
        options.showVersion = true;
        return 0;
    }
	if (isSetLong(parser, "verbose")) {
		options.verbose = true;
	}
	if (isSetLong(parser, "detailed")) {
		options.detailed = true;
	}
	if (isSetLong(parser, "pairwiseCount")) {
		options.pairwiseCount = true;
	}
	if (isSetLong(parser, "tripletCount")) {
		options.tripletCount = true;
	}
	if (isSetLong(parser, "score")) {
		options.sopScore = true;
	}
	if (isSetLong(parser, "scoreMatrix")) {
		CharString matrix;
		getOptionValueLong(parser, "scoreMatrix", matrix);
		options.matrixFile = matrix;
		toUpper(matrix);
		if (matrix == "HOXD") options.hoxdMatrix = true;
		else if (matrix == "COMPUTE") options.computeMatrix = true;
	}
	if (isSetLong(parser, "gapOpen")) {
		getOptionValueLong(parser, "gapOpen", options.gapOpen);
	}
	if (isSetLong(parser, "gapExtend")) {
		getOptionValueLong(parser, "gapExtend", options.gapExtend);
	}
	if (isSetLong(parser, "inFormat")) {
		CharString format;
		getOptionValueLong(parser, "inFormat", format);
		toUpper(format);
		if (format == "XMFA") options.inputFormat = XMFA;
		else if (format == "MAF") options.inputFormat = MAF;
		else { 
			std::cerr << "ERROR: Invalid input format." << std::endl;
			return 1;
		}
	}
	if (isSetShort(parser, "x")) {
		options.swapPositionsXmfa = true;
	}

	if (options.gapOpen > 0) std::cerr << "WARNING: Positive gap open penalty." << std::endl;
	if (options.gapExtend > 0) std::cerr << "WARNING: Positive gap extension penalty." << std::endl;

	return 0;
}

int mainWithOptions(Options & options)
{
	typedef Dna5 TAlphabet;
	typedef double TScoreValue;

	typedef String<TAlphabet> TSequence;
	typedef Size<TSequence>::Type TSize;

	typedef StringSet<TSequence> TStringSet;

	typedef Align<TSequence, ArrayGaps> TAlign;
	typedef std::map<CharString, AlignmentBlockRow<TSize, TSize> > TIdRowMap;

	if (options.verbose)
	{
		std::cout << "Alignment file: " << options.inputFile << std::endl;
		if (!options.computeMatrix && !options.hoxdMatrix)
			std::cout << "Scoring matrix file: " << options.matrixFile << std::endl;

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
		std::cerr << "ERROR: Parsing file " << options.inputFile << " in MAF format failed!" << std::endl;
		return 1;
	}
	else if (options.inputFormat == XMFA && options.swapPositionsXmfa && parseAlignment(aligns, idToRowMaps, seqs, inStream, options.verbose, XmfaSwap())) 
	{
		std::cerr << "ERROR: Parsing file " << options.inputFile << " in XMFA format failed!" << std::endl;
		return 1;
	}
	else if (options.inputFormat == XMFA && !options.swapPositionsXmfa && parseAlignment(aligns, idToRowMaps, seqs, inStream, options.verbose, Xmfa())) 
	{
		std::cerr << "ERROR: Parsing file " << options.inputFile << " in XMFA format failed!" << std::endl;
		return 1;
	}
	inStream.close();

	// Compute the sum-of-pairs alignment score
	if (options.sopScore)
	{
		if (options.verbose) std::cout << "Computing sum-of-pairs score..." << std::endl;

		Score<TScoreValue, ScoreMatrix<TAlphabet, Default> > score(static_cast<TScoreValue>(options.gapExtend), static_cast<TScoreValue>(options.gapOpen));

		if (options.computeMatrix)
		{
			CharString matrixFile = options.inputFile;
			append(matrixFile, ".matrix");
			if (sopScore(aligns, score, matrixFile))
			{
				std::cerr << "ERROR: Computation of sum-of-pairs score failed!" << std::endl;
				return 1;
			}
		}
		//else if (options.hoxdMatrix)
		//{
		//	if (sopScore(aligns, Score<TSCoreValue, ScoreMatrix<Dna, Hoxd> >()))
		//	{
		//		std::cerr << "ERROR: Computation of sum-of-pairs score failed!" << std::endl;
		//		return 1;
		//	}
		//}
		else // read score matrix from file
		{
			char * file = toCString(options.matrixFile);
			loadScoreMatrix(score, file);

			if (sopScore(aligns, score))
			{
				std::cerr << "ERROR: Computation of sum-of-pairs score failed!" << std::endl;
				return 1;
			}
		}
	}
	
	// Compute the breakpoint counts
	if (options.pairwiseCount || options.tripletCount)
	{
		typedef int TBlockId;

		std::map<CharString, String<TBlockId> > blockSeqs;
		sequencesOfBlocks(blockSeqs, idToRowMaps);

		//for (std::map<CharString, String<TBlockId> >::const_iterator it = blockSeqs.begin(); it != blockSeqs.end(); ++it)
		//{
		//	std::cout << it->first << ": ";
		//	for (Iterator<String<TBlockId> >::Type itit = begin(it->second); itit != end(it->second); ++itit)
		//		std::cout << *itit << "  ";
		//	std::cout << std::endl;
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
