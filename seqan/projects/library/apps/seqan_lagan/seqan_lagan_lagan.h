/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
 ===========================================================================
  Copyright (C) 2010
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ===========================================================================
  Definitions and so on for the Lagan implementation in Seqan::LAGAN.
 ===========================================================================*/

#ifndef SEQAN_LAGAN_LAGAN_H_
#define SEQAN_LAGAN_LAGAN_H_

#include <iostream>
#include <list>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/seeds2.h>

#include "seqan_lagan.h"

using namespace seqan;

// ===========================================================================
// Tags, Enums, Classes, Specializations
// ===========================================================================

enum OutputFormat
{
    OUTPUT_SCREEN,
    OUTPUT_FASTA_ALIGN
};

struct _AlignmentScores {};
typedef Tag<_AlignmentScores> AlignmentScores;

struct _PairwiseAlignment {};
typedef Tag<_PairwiseAlignment> PairwiseAlignment;

struct _Lagan : PairwiseAlignment {};
typedef Tag<_Lagan> Lagan;

struct _ClassicDP : PairwiseAlignment {};
typedef Tag<_ClassicDP> ClassicDP;

template <>
struct Options<AlignmentScores>
{
    // The score value for a match.
    int scoreMatch;
    // The score value for a mismatch.
    int scoreMismatch;
    // The score value for opening a gap.
    int scoreGapOpen;
    // The score value for extending gaps.
    int scoreGapExtend;

    Options()
            : scoreMatch(3),
              scoreMismatch(-2),
              scoreGapOpen(-3),
              scoreGapExtend(-1)
    {}
};


template <>
struct Options<PairwiseAlignment> : Options<AlignmentScores>
{
    // Whether or not to show help and exit.
    bool showHelp;

    // I/O related.

    // The path to the output file.  Use "-" to write to stdout.
    CharString outputFilename;
    // The output format.
    OutputFormat outputFormat;

    Options()
            : Options<AlignmentScores>(),
              showHelp(false),
              // General options
              outputFilename("-"),
              outputFormat(OUTPUT_SCREEN)
    {}
};


template <>
struct Options<ClassicDP> : Options<PairwiseAlignment>
{
    Options()
            : Options<PairwiseAlignment>()
    {}
};


template <>
struct Options<Lagan> : Options<PairwiseAlignment>
{
    // LAGAN specific options

    // Minimal length at least one sequence must have for recursive
    // chaining to continue.
    unsigned sequenceLengthRecursionThreshold;
    // The start value for q-gram length q.
    unsigned qMax;
    // The end value for q-gram length q.
    unsigned qMin;
    // The minimum score a seed must have to be of high quality.
    int seedScoreThreshold;
    // Maximal horizontal/veritcal distance when chaining.  This is
    // the largest number of gaps to bridge with CHAOS chaining.
    unsigned chainingMaxDistance;
    // Maximal diagonal distance when chaining.
    unsigned chainingMaxDiagonalDistance;
    // The number to extend the upper/lower diagonals of seeds with
    // for banded chain alignment.
    unsigned chainAlignmentBandwidthDelta;

    Options()
            : Options<PairwiseAlignment>(),
              // LAGAN specific options
              sequenceLengthRecursionThreshold(10),  // was 200
              qMax(6),
              qMin(3),
              seedScoreThreshold(5), // was 30
              chainingMaxDistance(200),
              chainingMaxDiagonalDistance(5),
              chainAlignmentBandwidthDelta(5)
    {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<AlignmentScores> const & options)
{
    std::cerr << "Alignment Scores options {" << std::endl
              << "  scoreMatch: " << options.scoreMatch << ", " << std::endl
              << "  scoreMismatch: " << options.scoreMismatch << ", " << std::endl
              << "  scoreGapOpen: " << options.scoreGapOpen << ", " << std::endl
              << "  scoreGapExtend: " << options.scoreGapExtend << std::endl
              << "}" << std::endl;
    return stream;
}

template <typename TStream>
TStream & operator<<(TStream & stream, Options<PairwiseAlignment> const & options)
{
    stream << static_cast<Options<AlignmentScores> >(options);
    std::cerr << "Pairwise Alignment options {" << std::endl
              << "  outputFilename: " << options.outputFilename << std::endl
              << "  outputFormat: " << options.outputFormat << std::endl
              << "}" << std::endl;
    return stream;
}

template <typename TStream>
TStream & operator<<(TStream & stream, Options<ClassicDP> const & options)
{
    stream << static_cast<Options<PairwiseAlignment> >(options);
    return stream;
}

template <typename TStream>
TStream & operator<<(TStream & stream, Options<Lagan> const & options)
{
    stream << static_cast<Options<PairwiseAlignment> >(options);
    stream << "LAGAN options {" << std::endl
           << "  sequenceLengthRecursionThreshold: " << options.sequenceLengthRecursionThreshold << ", " << std::endl
           << "  qMax: " << options.qMax << ", " << std::endl
           << "  qMin: " << options.qMin << ", " << std::endl
           << "  seedScoreThreshold: " << options.seedScoreThreshold << ", " << std::endl
           << "  chainingMaxDistance: " << options.chainingMaxDistance << ", " << std::endl
           << "  chainingMaxDiagonalDistance: " << options.chainingMaxDiagonalDistance << ", " << std::endl
           << "  chainAlignmentBandwidthDelta: " << options.chainAlignmentBandwidthDelta << ", " << std::endl
           << "}" << std::endl;
    return stream;
}

void setUpCommandLineParser(CommandLineParser & parser,
                            AlignmentScores const &)
{
    addSection(parser, "Alignment Scores");

    addOption(parser, CommandLineOption("sm", "score-match", "Score value for a match.", OptionType::Int));
    addOption(parser, CommandLineOption("smm", "score-mismatch", "Score value for a mismatch.", OptionType::Int));
    addOption(parser, CommandLineOption("sg", "score-gap", "Score value for both opening and extending gaps.", OptionType::Int));
    addOption(parser, CommandLineOption("sgo", "score-gap-open", "Score value for opening.", OptionType::Int));
    addOption(parser, CommandLineOption("sge", "score-gap-extend", "Score value for extending gap.", OptionType::Int));
    addHelpLine(parser, "Score values are applied in the order of the");
    addHelpLine(parser, "command line options above.  Thus, gap open/");
    addHelpLine(parser, "extend overwrite gap.");
}


void setUpCommandLineParser(CommandLineParser & parser,
                            PairwiseAlignment const &)
{
    addSection(parser, "Input / Output");

    addOption(parser, CommandLineOption("o", "output-filename", "Path to the output file, '-' for stdout.", OptionType::String));
    addOption(parser, CommandLineOption("of", "output-format", "Alignment format to write. Default is 'screen'.", OptionType::String));
    addHelpLine(parser, "  screen   Suitable for visual inspection on a screen.");
    addHelpLine(parser, "  fasta    FASTA Alignment Format.");

    setUpCommandLineParser(parser, AlignmentScores());
}

void setUpCommandLineParser(CommandLineParser & parser,
                            ClassicDP const &)
{
    addVersionLine(parser, "SeqAn::LAGAN");

    addTitleLine(parser, "Classic pairwise DP alignment program.");
    addUsageLine(parser, "dp LEFT.fasta RIGHT.fasta");

    addLine(parser, "");
    addLine(parser, "At the moment, only the first sequences in each FASTA file is interpreted.");

    setUpCommandLineParser(parser, PairwiseAlignment());

    requiredArguments(parser, 2);
}

void setUpCommandLineParser(CommandLineParser & parser,
                            Lagan const &)
{
    addVersionLine(parser, "SeqAn::LAGAN");

    addTitleLine(parser, "An implementation of LAGAN in SeqAn.");
    addUsageLine(parser, "lagan LEFT.fasta RIGHT.fasta");
    addLine(parser, "");
    addLine(parser, "At the moment, only the first sequences in each FASTA file is interpreted.");

    setUpCommandLineParser(parser, PairwiseAlignment());

    requiredArguments(parser, 2);
}

int parseCommandLineAndCheck(Options<AlignmentScores> & options,
                             CommandLineParser & parser,
                             const int /*argc*/,
                             const char ** /*argv[]*/)
{
    // Apply alignment score related options.
    if (isSetLong(parser, "score-match"))
        getOptionValueLong(parser, "score-match", options.scoreMatch);
    if (isSetLong(parser, "score-mismatch"))
        getOptionValueLong(parser, "score-mismatch", options.scoreMatch);
    if (isSetLong(parser, "score-gap")) {
        getOptionValueLong(parser, "score-gap", options.scoreGapOpen);
        getOptionValueLong(parser, "score-gap", options.scoreGapExtend);
    }
    if (isSetLong(parser, "score-gap-open"))
        getOptionValueLong(parser, "score-gap-open", options.scoreGapOpen);
    if (isSetLong(parser, "score-gap-extend"))
        getOptionValueLong(parser, "score-gap-extend", options.scoreGapExtend);

    return 0;
}

int parseCommandLineAndCheck(Options<PairwiseAlignment> & options,
                             CommandLineParser & parser,
                             const int argc,
                             const char * argv[])
{
    if (!parse(parser, argc, argv)) {
        if (!isSetShort(parser, "h"))
            shortHelp(parser, std::cerr);
        return 1;
    }
    if (isSetShort(parser, "h")) {
        options.showHelp = true;
        return 0;
    }

    // Apply I/O options.
    if (isSetLong(parser, "output-filename"))
        getOptionValueLong(parser, "output-filename", options.outputFilename);
    if (isSetLong(parser, "output-format")) {
        CharString outputFormat;
        getOptionValueLong(parser, "output-format", outputFormat);
        if (outputFormat == "screen") {
            options.outputFormat = OUTPUT_SCREEN;
        } else if (outputFormat == "fasta") {
            options.outputFormat = OUTPUT_FASTA_ALIGN;
        } else {
            std::cerr << "Invalid output format!" << std::endl;
            return 1;
        }
    }

    int ret = parseCommandLineAndCheck(static_cast<Options<AlignmentScores> &>(options),
                                       parser, argc, argv);
    if (ret != 0)
        return ret;

    return 0;
}


int parseCommandLineAndCheck(Options<ClassicDP> & options,
                             CharString & leftFilename,
                             CharString & rightFilename,
                             CommandLineParser & parser,
                             const int argc,
                             const char * argv[])
{
    int ret = parseCommandLineAndCheck(static_cast<Options<PairwiseAlignment> &>(options),
                                       parser,
                                       argc,
                                       argv);
    if (ret != 0)
        return ret;

    // First argument is "lagan", second and third ones are file name.
    leftFilename = getArgumentValue(parser, 1);
    rightFilename = getArgumentValue(parser, 2);

    return 0;
}


int parseCommandLineAndCheck(Options<Lagan> & options,
                             CharString & leftFilename,
                             CharString & rightFilename,
                             CommandLineParser & parser,
                             const int argc,
                             const char * argv[])
{
    int ret = parseCommandLineAndCheck(static_cast<Options<PairwiseAlignment> &>(options),
                                       parser,
                                       argc,
                                       argv);
    if (ret != 0)
        return ret;

    // First argument is "lagan", second and third ones are file name.
    leftFilename = getArgumentValue(parser, 1);
    rightFilename = getArgumentValue(parser, 2);
    
    return 0;
}


template <typename TSeed, typename TSequence>
void constructLaganChain(
        std::list<TSeed> & chain,
        TSequence const & sequence0,
        TSequence const & sequence1,
        Options<Lagan> const & options)
{
    if (length(sequence0) <= options.sequenceLengthRecursionThreshold &&
        length(sequence1) <= options.sequenceLengthRecursionThreshold)
        return;  // Break recursion.

    Score<int, Simple> scoringScheme(options.scoreMatch, options.scoreMismatch, options.scoreGapOpen, options.scoreGapExtend);

    // -----------------------------------------------------------------------
    // Find seeds
    // -----------------------------------------------------------------------
    typedef typename Position<TSeed>::Type TPosition;
	typedef SeedSet<Simple, Unordered, DefaultSeedSetConfigScore> TSeedSet;
    // typedef typename Value<TSeedSet>::Type TSeed;
	typedef Index<TSequence, Index_QGram<SimpleShape > > TQGramIndex;
	typedef Finder<TQGramIndex> TFinder;

	TSeedSet seedSet;
    setMinScoreThreshold(seedSet, options.seedScoreThreshold);

    // std::cerr << "Creating index of sequence1;  length(sequence1) == " << length(sequence1) << std::endl;
	TQGramIndex qGramIndex(sequence1);
	TFinder finder(qGramIndex);

    std::cerr << "Iterating over qgrams" << std::endl;
    unsigned q = options.qMax;
	while (length(seedSet) == 0) {
        // std::cerr << "Iteration..." << std::endl;
        // Do not recurse if we have reached the minimal q-gram size.
		if (q < options.qMin)
            return;

        // Search all q-grams from sequence0 in sequence1 with the
        // q-gram index.
		resize(indexShape(qGramIndex), q);
		for (unsigned i = 0; i < length(sequence0) - q + 1; ++i) {
			while (find(finder, TSequence(infix(sequence0, i, i + q)))) {
				TPosition pos0 = beginPosition(sequence0) + i;
				TPosition pos1 = beginPosition(sequence1) + position(finder);
                TSeed seed(pos0, pos1, q);
                setScore(seed, q * scoreMatch(scoringScheme));
                // std::cerr << "Adding " << seed << ", score == " << getScore(seed) << std::endl;
                if (addSeed(seedSet, seed, options.chainingMaxDistance, Nothing(), scoringScheme/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge())) {
                    // std::cerr << "  by merging" << std::endl;
                    continue;
                }
				if (addSeed(seedSet, seed, options.chainingMaxDistance, options.chainingMaxDiagonalDistance, scoringScheme, sequence0, sequence1, Chaos())) {
                    // std::cerr << "  by CHAOS chaining" << std::endl;
                    continue;
                }
                // std::cerr << "  as single seed" << std::endl;
				addSeed(seedSet, seed, Single());
			}
			clear(finder);
		}
		--q;
	}

    // -----------------------------------------------------------------------
    // Perform global chaining of these seeds
    // -----------------------------------------------------------------------
    std::cerr << "Global chaining..." << std::endl;
    // std::cerr << "length(seedSet) == " << length(seedSet) << std::endl;
    chainSeedsGlobally(chain, seedSet, SparseChaining());

    // -----------------------------------------------------------------------
    // Recursively fill gaps
    // -----------------------------------------------------------------------
	// if (q > q_min)
	// {
	// 	std::list<TSeed> subchain;
	// 	typedef typename std::list<TSeed>::iterator TIterator;

	// 	TIterator it = chain.begin();
	// 	TIterator it2 = it; 
	// 	++it2;

	// 	laganChaining(subchain, 
	// 		infix(host(a), beginPosition(a), leftDim0(*it)), 
	// 		infix(host(b), beginPosition(b), leftDim1(*it)), q);
	// 	chain.splice(it, subchain);

	// 	while(it2 != chain.end())
	// 	{
	// 		laganChaining(subchain, 
	// 			infix(host(a), rightDim0(*it), leftDim0(*it2)), 
	// 			infix(host(b), rightDim1(*it), leftDim1(*it2)), q);
	// 		chain.splice(it2, subchain);

	// 		it = it2;
	// 		++it2;
	// 	}

	// 	laganChaining(subchain, 
	// 		infix(host(a), rightDim0(*it), endPosition(a)), 
	// 		infix(host(b), rightDim1(*it), endPosition(b)), q);
	// 	chain.splice(it2, subchain);
	// }
}

// Helper for output writing
template <typename TStream, typename TSequenceIds, typename TAlignment, typename TAlgorithm>
int writeOutputPairwiseToStream(TStream & stream,
                                TAlignment const & alignment,
                                TSequenceIds const & sequenceIds,
                                Options<TAlgorithm> const & options)
{
    if (options.outputFormat == OUTPUT_SCREEN) {
        stream << "Alignment of sequences " << sequenceIds[0] << " and " << sequenceIds[1] << std::endl;
        stream << alignment;
        return 0;
    } else if (options.outputFormat == OUTPUT_FASTA_ALIGN) {
        write(stream, alignment, sequenceIds, FastaAlign());
        return 0;
    } else {
        SEQAN_ASSERT_FAIL("Invalid output format!");
        return 1;
    }
}

// Write the given alignment to a file or stdout, depending on options.
template <typename TAlignment, typename TSequenceIds, typename TAlgorithm>
int writeOutputPairwise(TAlignment const & alignment,
                        TSequenceIds const & sequenceIds,
                        Options<TAlgorithm> const & options)
{
    // Dynamically dispatch to the correct stream type.
    if (options.outputFilename == CharString("-")) {
        return writeOutputPairwiseToStream(std::cout, alignment, sequenceIds, options);
    } else {
        std::fstream outFile(toCString(options.outputFilename), std::fstream::out | std::fstream::binary);
        return writeOutputPairwiseToStream(outFile, alignment, sequenceIds, options);
    }
}


template <typename TAlignment>
int performPairwiseAlignment(
        TAlignment & alignment,
        Options<Global> const & /*globalOptions*/,
        Options<Lagan> const & options,
        Lagan const &)
{
    clearGaps(row(alignment, 0));
    clearGaps(row(alignment, 1));

    typedef typename Row<TAlignment>::Type TRow;
    typedef typename Source<TRow>::Type TSequence;
    TSequence & sequence0 = source(row(alignment, 0));
    TSequence & sequence1 = source(row(alignment, 1));
    
    // Execute LAGAN algorithm which constructs a chain of seeds.
    typedef Seed<Simple, DefaultSeedConfigScore> TSeed;
    typedef std::list<TSeed> TSeedChain;
    TSeedChain seedChain;
    constructLaganChain(seedChain, sequence0, sequence1, options);
    if (length(seedChain) == 0) {
        std::cerr << "ERROR: No similarity found!" << std::endl;
        return 1;
    }

    // Perform banded alignment around this seed.
    Score<int, Simple> scoringScheme(options.scoreMatch, options.scoreMismatch, options.scoreGapExtend, options.scoreGapOpen);
    int score = bandedChainAlignment(alignment, seedChain, options.chainAlignmentBandwidthDelta, scoringScheme, AlignConfig<false, false, false, false>());
    std::cerr << "Alignment score is " << score << std::endl;
    
    return 0;
}


template <typename TAlignment>
int performPairwiseAlignment(
        TAlignment & alignment,
        Options<Global> const & /*globalOptions*/,
        Options<ClassicDP> const & options,
        ClassicDP const &)
{
    clearGaps(row(alignment, 0));
    clearGaps(row(alignment, 1));

    Score<int, Simple> scoringScheme(options.scoreMatch, options.scoreMismatch, options.scoreGapExtend, options.scoreGapOpen);

    int score = 0;
    if (options.scoreGapOpen == options.scoreGapExtend) {
        score = globalAlignment(alignment, scoringScheme, NeedlemanWunsch());
    } else {
        score = globalAlignment(alignment, scoringScheme, Gotoh());
    }
    std::cerr << "Alignment score is " << score << std::endl;
    return 0;
}


template <typename TAlgorithm>
int executeParwiseAlignmentCommand(
        Options<Global> const & globalOptions,
        Options<TAlgorithm> const & options,
        CharString const & leftFilename,
        CharString const & rightFilename,
        TAlgorithm const & tag)
{
    typedef Dna5String TSequence;
    typedef CharString TSequenceIdentifier;
    typedef ArrayGaps TAlignSpec;

    std::cerr << "Pairwise LAGAN alignment of " << leftFilename << " and " << rightFilename << std::endl << std::endl;
    std::cerr << std::endl;
    std::cerr << options;
    std::cerr << std::endl;

    // Declare variables for sequences and sequence identifiers.
    TSequence sequence0;
    TSequence sequence1;
    StringSet<TSequenceIdentifier> sequenceIds;

    // Read sequence files.
    {
        TSequenceIdentifier sequenceId;
        std::fstream f0(toCString(leftFilename));
        readID(f0, sequenceId, Fasta());
        appendValue(sequenceIds, sequenceId);
        read(f0, sequence0, Fasta());
        std::fstream f1(toCString(rightFilename));
        readID(f1, sequenceId, Fasta());
        appendValue(sequenceIds, sequenceId);
        read(f1, sequence1, Fasta());
    }

    std::cerr << "Sequence Information" << std::endl;
    std::cerr << "  sequence 0 is " << sequenceIds[0] << std::endl;
    std::cerr << "  sequence 1 is " << sequenceIds[1] << std::endl;

    // Initialize Align object.
    Align<TSequence, TAlignSpec> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequence0);
    assignSource(row(alignment, 1), sequence1);
    
    // Perform pairwise alignment.
    int ret = performPairwiseAlignment(alignment, globalOptions, options, tag);
    if (ret != 0)
        return ret;
    // Write output.
    ret = writeOutputPairwise(alignment, sequenceIds, options);
    if (ret != 0)
        return ret;

    return 0;
}

#endif  // SEQAN_LAGAN_LAGAN_H_

