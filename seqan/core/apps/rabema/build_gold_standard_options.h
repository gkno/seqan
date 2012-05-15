// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef APPS_RABEMA_BUILD_GOLD_STANDARD_OPTIONS_H_
#define APPS_RABEMA_BUILD_GOLD_STANDARD_OPTIONS_H_

#include <iostream>

#include <seqan/misc/misc_cmdparser.h>

#include "rabema.h"

// ============================================================================
// Enums, Tags, Classes, Typedefs.
// ============================================================================

// Options specialization for the build gold standard subprogram.
template <>
class Options<BuildGoldStandard>
{
public:
    // True iff help is to be shown.
    bool showHelp;

    // Whether hits should be verified with MyersUkonnen search.
    bool verify;

    // Whether N should match all characters without penalty.
    bool matchN;

    // Whether or not to run in "oracle Sam mode", i.e. the Sam file was
    // generated by a read simulator.  If this is the case, the maximal
    // error rate for each read is determined by the alignment given in
    // the Sam file.  The generated error rate for each alignment is output
    // to be 0.
    bool oracleSamMode;
    
    // Maximum errors in percent, relative to the read length.
    double maxError;

    // Distance function to use, also see validDistanceFunction.
    String<char> distanceFunction;

    // Path to output file.
    CharString outFileName;

    // Name of reference contigs file name.
    String<char> referenceSeqFilename;

    // Name of Sam file with golden reads.
    String<char> perfectMapFilename;
    
    // Verbosity level.
    int verbosity;

    Options()
            : showHelp(false),
              verify(false),
              matchN(false),
              oracleSamMode(false),
              maxError(0),
              distanceFunction("edit"),
              outFileName("-"),
              verbosity(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Check whether the distanceFunction member variable of the build
// gold standard options specialization is valid.
bool validDistanceFunction(Options<BuildGoldStandard> const & options)
{
    if (options.distanceFunction == "hamming")
        return true;
    if (options.distanceFunction == "edit")
        return true;
    return false;
}

// Set up the command line parser, adding options for the gold
// standard building subprogram.
void setUpCommandLineParser(CommandLineParser & parser, BuildGoldStandard const & /*tag*/)
{
    setUpCommandLineParser(parser);
    
    // Add usage lines.
    addUsageLine(parser, "build_standard [OPTIONS] <REFERENCE SEQ> <PERFECT MAP>");

    // Set options.
    addOption(parser, CommandLineOption("e", "max-error-rate", "the maximal error in percent of read length, default: 0", OptionType::Double));
    addOption(parser, CommandLineOption("x", "verify", "Verify Result.  Default: true.", OptionType::Boolean));
    addOption(parser, CommandLineOption("d", "distance-function", ("the distance function to use, default: hamming"), OptionType::String | OptionType::Label));
    addHelpLine(parser, "hamming          = Hamming distance");
    addHelpLine(parser, "edit             = Edit distance");
    addOption(parser, CommandLineOption("o", "out-file", ("Path to the output file.  Use \"-\" for stdout, default is \"-\""), OptionType::String));
    addOption(parser, CommandLineOption("mN" , "match-N", "If specified, N characters match all other characters.  The default is for them to mismatch with CGAT.", OptionType::Boolean));
    addOption(parser, CommandLineOption("os" , "oracle-sam", "If specified, wit_builder is run in oracle Sam mode.  Use this for Sam files generated by a read simulator.  Default: false", OptionType::Boolean));
    
    requiredArguments(parser, 3);
}

// Parse command line parameters and validate them for the build
// standard building subprogram.
int parseCommandLineAndCheck(Options<BuildGoldStandard> & options, CommandLineParser & parser, int const argc, char const * argv[])
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

    // Get arguments.
    if (isSetLong(parser, "verify"))
        options.verify = true;
    if (isSetLong(parser, "max-error-rate"))
        getOptionValueLong(parser, "max-error-rate", options.maxError);
    if (isSetLong(parser, "distance-function"))
        getOptionValueLong(parser, "distance-function", options.distanceFunction);
    if (isSetLong(parser, "match-N"))
        options.matchN = true;
    if (isSetLong(parser, "out-file"))
        getOptionValueLong(parser, "out-file", options.outFileName);
    if (isSetLong(parser, "oracle-sam"))
        options.oracleSamMode = true;

    // Validate values.
    if (options.maxError < 0) {
        std::cerr << "ERROR: Invalid maximum error value: " << options.maxError << std::endl;
        return 1;
    }
    if (!validDistanceFunction(options)) {
      std::cerr << "ERROR: Invalid distance function: " << options.distanceFunction << std::endl;
      return 1;
    }
    
    // Get positional arguments.
    options.referenceSeqFilename = getArgumentValue(parser, 1);
    options.perfectMapFilename = getArgumentValue(parser, 2);

    return 0;
}

#endif  // #ifndef APPS_RABEMA_BUILD_GOLD_STANDARD_OPTIONS_H_
