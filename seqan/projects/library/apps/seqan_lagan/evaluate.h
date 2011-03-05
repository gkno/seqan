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
  Code for evaluating an alignment.
 ===========================================================================*/

#ifndef SEQAN_LAGAN_EVALUATE_H_
#define SEQAN_LAGAN_EVALUATE_H_

#include <fstream>
#include <iostream>

using namespace seqan;

// ===========================================================================
// Tags, Enums, Classes, Specializations
// ===========================================================================

struct Evaluate_ {};
typedef Tag<Evaluate_> Evaluate;

template <>
struct Options<Evaluate> : Options<AlignmentScores>
{
    // Whether or not to show help and exit.
    bool showHelp;

    Options() : Options<AlignmentScores>(),
                showHelp(false) {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

void setUpCommandLineParser(CommandLineParser & parser,
                            Evaluate const &)
{
    addVersionLine(parser, "SeqAn::LAGAN");

    addTitleLine(parser, "Evaluation of alignments.");
    addUsageLine(parser, "evaluate ALIGNMENT.fasta");
    addLine(parser, "");
    addLine(parser, "At the moment, only the FASTA alignment format is supported.");

    setUpCommandLineParser(parser, AlignmentScores());

    requiredArguments(parser, 1);
}

int parseCommandLineAndCheck(Options<Evaluate> & options,
                             CharString & alignmentFilename,
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

    int ret = parseCommandLineAndCheck(static_cast<Options<AlignmentScores> &>(options),
                                       parser, argc, argv);
    if (ret != 0)
        return ret;

    // First argument is "evaluate", second one is file name.
    alignmentFilename = getArgumentValue(parser, 1);

    return 0;
}

template <typename TAlignment>
int evaluateAlignment(TAlignment const & alignment,
                      Options<Evaluate> const & options)
{
    typedef typename Row<TAlignment>::Type TAlignmentRow;
    typedef typename Iterator<TAlignmentRow>::Type TRowIterator;
    typedef typename Position<TAlignmentRow>::Type TPosition;

    // Get first and last column of alignment.
    // TODO(holtgrew): Why do iterators not work?
    TPosition beginPos = beginPosition(cols(alignment));
    TPosition endPos = endPosition(cols(alignment));
    // State whether we already are in a gap in the given sequence.
    bool inGap0 = false;
    bool inGap1 = false;

    int totalScore = 0;
    for (unsigned i = beginPos; i < endPos; ++i) {
        SEQAN_ASSERT_MSG(!isGap(row(alignment, 0), i) || !isGap(row(alignment, 1), i), "i = %u", i);
        if (isGap(row(alignment, 0), i)) {
            if (inGap0)
                totalScore += options.scoreGapExtend;
            else
                totalScore += options.scoreGapOpen;
            inGap0 = true;
        } else {
            inGap0 = false;
        }
        if (isGap(row(alignment, 1), i)) {
            if (inGap1)
                totalScore += options.scoreGapExtend;
            else
                totalScore += options.scoreGapOpen;
            inGap1 = true;
        } else {
            inGap1 = false;
        }
        if (!inGap0 && !inGap1) {
            if (row(alignment, 0)[i] == row(alignment, 1)[i]) {
                totalScore += options.scoreMatch;
            } else {
                totalScore += options.scoreMismatch;
            }
        }
    }

    return totalScore;
}

int executeEvaluation(
        Options<Global> const & /*globalOptions*/,
        Options<Evaluate> const & options,
        CharString const & alignmentFilename,
        Evaluate const &)
{
    typedef Dna5String TString;
    typedef Align<Dna5String> TAlignment;

    // Load alignment from file.
    TAlignment alignment;
    std::fstream f(toCString(alignmentFilename));
    read(f, alignment, FastaAlign());

    std::cerr << "Evaluation of pairwise alignment." << std::endl;
    std::cerr << options;
    std::cerr << std::endl;

    // Compute alignment score.
    int alignmentScore = evaluateAlignment(alignment, options);
    std::cout << "Alignment Score: " << alignmentScore << std::endl;
    
    return 0;
}

#endif  // #ifndef SEQAN_LAGAN_EVALUATE_H_
