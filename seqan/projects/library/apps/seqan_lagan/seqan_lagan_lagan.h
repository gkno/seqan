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

struct _Lagan;
typedef Tag<_Lagan> Lagan;

template <>
struct Options<Lagan>
{
    // Options common to all commands, could to into a Common
    // specialization of options later.

    // Whether or not to show help and exit.
    bool showHelp;

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
    // The score value for a match.
    int scoreMatch;
    // The score value for a mismatch.
    int scoreMismatch;
    // The score value for opening a gap.
    int scoreGapOpen;
    // The score value for extending gaps.
    int scoreGapExtend;

    Options()
            : showHelp(false),
              // LAGAN specific options
              sequenceLengthRecursionThreshold(10),  // was 200
              qMax(6),
              qMin(3),
              seedScoreThreshold(5), // was 30
              chainingMaxDistance(200),
              chainingMaxDiagonalDistance(5),
              chainAlignmentBandwidthDelta(5),
              scoreMatch(3),
              scoreMismatch(-2),
              scoreGapOpen(-1),
              scoreGapExtend(-1)  // TODO(holtgrew): Change to -3 when affine banded chain alignment is in place.
    {}
};

// ===========================================================================
// Tags, Enums, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<Lagan> const & options)
{
    stream << "LAGAN options {" << std::endl
           << "  sequenceLengthRecursionThreshold: " << options.sequenceLengthRecursionThreshold << ", " << std::endl
           << "  qMax: " << options.qMax << ", " << std::endl
           << "  qMin: " << options.qMin << ", " << std::endl
           << "  seedScoreThreshold: " << options.seedScoreThreshold << ", " << std::endl
           << "  chainingMaxDistance: " << options.chainingMaxDistance << ", " << std::endl
           << "  chainingMaxDiagonalDistance: " << options.chainingMaxDiagonalDistance << ", " << std::endl
           << "  chainAlignmentBandwidthDelta: " << options.chainAlignmentBandwidthDelta << ", " << std::endl
           << "  scoreMatch: " << options.scoreMatch << ", " << std::endl
           << "  scoreMismatch: " << options.scoreMismatch << ", " << std::endl
           << "  scoreGapOpen: " << options.scoreGapOpen << ", " << std::endl
           << "  scoreGapExtend: " << options.scoreGapExtend << std::endl
           << "}" << std::endl;
    return stream;
}

void setUpCommandLineParser(CommandLineParser & parser,
                            Lagan const &)
{
    addVersionLine(parser, "SeqAn::LAGAN");

    addTitleLine(parser, "An implementation of LAGAN in SeqAn.");
    addUsageLine(parser, "lagan LEFT.fasta RIGHT.fasta");
    addLine(parser, "");
    addLine(parser, "At the moment, only the first sequences in each FASTA file is interpreted.");

    requiredArguments(parser, 2);
}


int parseCommandLineAndCheck(Options<Lagan> & options,
                             CharString & leftFilename,
                             CharString & rightFilename,
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

    // std::cout << "Creating index of sequence1;  length(sequence1) == " << length(sequence1) << std::endl;
	TQGramIndex qGramIndex(sequence1);
	TFinder finder(qGramIndex);

    std::cout << "Iterating over qgrams" << std::endl;
    unsigned q = options.qMax;
	while (length(seedSet) == 0) {
        // std::cout << "Iteration..." << std::endl;
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
                std::cout << "Adding " << seed << ", score == " << getScore(seed) << std::endl;
                if (addSeed(seedSet, seed, options.chainingMaxDistance, Nothing(), scoringScheme/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge())) {
                    std::cout << "  by merging" << std::endl;
                    continue;
                }
				if (addSeed(seedSet, seed, options.chainingMaxDistance, options.chainingMaxDiagonalDistance, scoringScheme, sequence0, sequence1, Chaos())) {
                    std::cout << "  by CHAOS chaining" << std::endl;
                    continue;
                }
                std::cout << "  as single seed" << std::endl;
				addSeed(seedSet, seed, Single());
			}
			clear(finder);
		}
		--q;
	}

    // -----------------------------------------------------------------------
    // Perform global chaining of these seeds
    // -----------------------------------------------------------------------
    std::cout << "Global chaining..." << std::endl;
    std::cout << "length(seedSet) == " << length(seedSet) << std::endl;
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


int executeCommand(Options<Global> const & /*globalOptions*/,
                   Options<Lagan> const & options,
                   CharString const & leftFilename,
                   CharString const & rightFilename,
                   Lagan const &)
{
    typedef Dna5String TSequence;
    typedef CharString TSequenceIdentifier;
    typedef ArrayGaps TAlignSpec;

    std::cout << "Pairwise LAGAN alignment of " << leftFilename << " and " << rightFilename << std::endl;
    std::cout << std::endl;
    std::cout << options;
    std::cout << std::endl;

    // Declare variables for sequences and sequence identifiers.
    TSequence sequence0;
    TSequenceIdentifier sequence0Id;
    TSequence sequence1;
    TSequenceIdentifier sequence1Id;

    // Read sequence files.
    {
        std::fstream f0(toCString(leftFilename));
        readID(f0, sequence0Id, Fasta());
        read(f0, sequence0, Fasta());
        std::fstream f1(toCString(rightFilename));
        readID(f1, sequence1Id, Fasta());
        read(f1, sequence1, Fasta());
    }

    std::cout << "Sequence Information" << std::endl;
    std::cout << "  sequence 0 is " << sequence0Id << std::endl;
    std::cout << "  sequence 1 is " << sequence1Id << std::endl;

    // Execute LAGAN algorithm which constructs a chain of seeds.
    typedef Seed<Simple, DefaultSeedConfigScore> TSeed;
    typedef std::list<TSeed> TSeedChain;
    TSeedChain seedChain;
    constructLaganChain(seedChain, sequence0, sequence1, options);
    if (length(seedChain) == 0) {
        std::cout << "ERROR: No similarity found!" << std::endl;
        return 1;
    }

    // Perform banded alignment around this seed.
    Align<TSequence, TAlignSpec> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequence0);
    assignSource(row(alignment, 1), sequence1);
    // TODO(holtgrew): Parameter order in bandedChainAlignment is subject to change, as is the name.
    Score<int, Simple> scoringScheme(options.scoreMatch, options.scoreMismatch, options.scoreGapOpen, options.scoreGapExtend);
    int score = bandedChainAlignment(alignment, seedChain, options.chainAlignmentBandwidthDelta, scoringScheme, AlignConfig<false, false, false, false>());

    // Write result to output file.
    std::cout << "LAGAN Alignment (score == " << score << ")" << std::endl
              << alignment;

    // Perform NW alignment.
    clearGaps(row(alignment, 0));
    clearGaps(row(alignment, 1));
    score = globalAlignment(alignment, scoringScheme, AlignConfig<false, false, false, false>(), NeedlemanWunsch());

    std::cout << "NW (score == " << score << ")" << std::endl
              << alignment;
    
    return 0;
}

#endif  // SEQAN_LAGAN_LAGAN_H_

