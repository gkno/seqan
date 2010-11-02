#ifndef APPS_RABEMA_EVALUATION_H_
#define APPS_RABEMA_EVALUATION_H_

#include "rabema.h"

#include <seqan/misc/misc_cmdparser.h>

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

template <>
struct Options<EvaluateResults>
{
    // True iff help is to be shown.
    bool showHelp;

    // Maximum number or errors per read length in percent.
    int maxError;

    // Print the missed intervals to stderr for debugging purposes.
    bool showMissedIntervals;

    // Print superflous intervals (intervals found in SAM file but have too bad score).
    bool showSuperflousIntervals;

    // Print additional intervals (intervals found in SAM with good score that are not in WIT file).
    bool showAdditionalIntervals;

    // Print the hit intervals to stderr for debugging purposes.
    bool showHitIntervals;

    // Print each end position that we try to match agains the interval.
    bool showTryHitIntervals;

    // If enabled, don't panic on additional hits in non-weighted mode.
    bool dontPanic;

    // If enabled, the distance of alignments is considered to be 0.  This is used when comparing against WIT files that were generated by the read simulator.  In this case, the intervals will have a length of 1 and be the single original position of this read.
    bool oracleWitMode;

    // If true, N matches as a wildcard.  Otherwise it matches none.
    bool matchN;

    // If true, use weighted distances instead of unit ones.
    bool weightedDistances;

    // The benchmark category, one of {"all", "any-best", "all-best"}.
    CharString benchmarkCategory;

    // Distance function to use, also see validDistanceFunction.
    CharString distanceFunction;

    // Output filename.
    CharString outFileName;

    // Name of reference sequence file.
    CharString seqFileName;

    // Name of WIT file.
    CharString witFileName;

    // Name of SAM file.
    CharString samFileName;

    Options() : showHelp(false) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Return true iff distanceFunction is a valid distance function.
// Valid distances are one of {"hamming", "edit"}.
bool validDistanceFunction(Options<EvaluateResults> const & options)
{
    if (options.distanceFunction == "hamming")
        return true;
    if (options.distanceFunction == "edit")
        return true;
    return false;
}

// Return true iff benchmarkCategory is a valid benchmark category,
// i.e. one of {"all", "any-best", "all-best"}.
bool validBenchmarkCategory(Options<EvaluateResults> const & options)
{
    if (options.benchmarkCategory == "all")
        return true;
    if (options.benchmarkCategory == "any-best")
        return true;
    if (options.benchmarkCategory == "all-best")
        return true;
    return false;
}

// Set up the command line parser, adding options for the evaluation subprogram.
void setUpCommandLineParser(CommandLineParser & parser, EvaluateResults const & /*tag*/)
{
    setUpCommandLineParser(parser);

    // Add usage lines.
    addUsageLine(parser, "compare [OPTIONS] <SEQUENCE FILE> <SAM FILE> <WIT FILE>");

    // Set options.
    addOption(parser, CommandLineOption("e", "max-error-rate", "Maximal error rate in percent [0].", OptionType::Int));
    addOption(parser, CommandLineOption("d", "distance-function", "Distance function to use [hamming]", OptionType::String | OptionType::Label));
    addHelpLine(parser, "hamming          = Hamming distance");
    addHelpLine(parser, "edit             = Edit distance");
    addOption(parser, CommandLineOption("o", "out-file", "Path to the output file.  (\"-\")", OptionType::String));
    addOption(parser, CommandLineOption("sm", "show-missed-intervals", "Missed intervals are printed to stderr for debugging if this option is set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("sh", "show-hit-intervals", "Hit intervals are printed to stderr for debugging if this option is set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("st", "show-try-hit-intervals", "Last positions tried to hit against an intervals are printd to stderr if set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("ss", "show-superflous-intervals", "Intervals that are in the SAM file but have a too bad score are printed to stderr if set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("sa", "show-additional-intervals", "Intervals that are in the SAM with good score but not in WIT file to stderr if set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("mN", "match-N", "If set, N matches all as a wildcard character, otherwise it never matches.", OptionType::Boolean));
    addOption(parser, CommandLineOption("wm", "weighted-distances", "If set, use weighted distances instead of unit ones.", OptionType::Boolean));
    addOption(parser, CommandLineOption("c", "benchmark-category", "The benchmark category to compare for.  One of {all, any-best, all-best}.  Default: all.", OptionType::String));
    addOption(parser, CommandLineOption("DP", "DONT-PANIC", "Don't panic on additional hits in non-weightedmode.  Default: False.", OptionType::Boolean));
    addOption(parser, CommandLineOption("ow", "oracle-wit-mode", "Enable 'oracle wit mode', alignment distances are ignored, assumed to be 0.  Default: False.", OptionType::Boolean));
    
    // We require 4 command line options.
    requiredArguments(parser, 4);
}

// Parse command line parameters and validate them for the build
// standard building subprogram.
int parseCommandLineAndCheck(Options<EvaluateResults> & options, CommandLineParser & parser, int const argc, char const * argv[])
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
    if (isSetLong(parser, "DONT-PANIC"))
        options.dontPanic = true;
    if (isSetLong(parser, "oracle-wit-mode"))
        options.oracleWitMode = true;
    if (isSetLong(parser, "match-N"))
        options.matchN = true;
    if (isSetLong(parser, "weighted-distances"))
        options.weightedDistances = true;
    if (isSetLong(parser, "show-missed-intervals"))
        options.showMissedIntervals = true;
    if (isSetLong(parser, "show-hit-intervals"))
        options.showHitIntervals = true;
    if (isSetLong(parser, "show-try-hit-intervals"))
        options.showTryHitIntervals = true;
    if (isSetLong(parser, "show-superflous-intervals"))
        options.showSuperflousIntervals = true;
    if (isSetLong(parser, "show-additional-intervals"))
        options.showAdditionalIntervals = true;
    if (isSetLong(parser, "max-error-rate"))
        getOptionValueLong(parser, "max-error-rate", options.maxError);
    if (isSetLong(parser, "distance-function"))
        getOptionValueLong(parser, "distance-function", options.distanceFunction);
    if (isSetLong(parser, "out-file"))
        getOptionValueLong(parser, "out-file", options.outFileName);
    if (isSetLong(parser, "benchmark-category"))
        getOptionValueLong(parser, "benchmark-category", options.benchmarkCategory);

    // Validate values.
    if (options.maxError < 0) {
        std::cerr << "ERROR: Invalid maximum error value: " << options.maxError << std::endl;
        return 1;
    }
    if (!validDistanceFunction(options)) {      
      std::cerr << "ERROR: Invalid distance function: " << options.distanceFunction << std::endl;
      return 1;
    }
    if (!validBenchmarkCategory(options)) {
        std::cerr << "ERROR: Invalid benchmark category: " << options.benchmarkCategory << std::endl;
        return 1;
    }
    
    // Get positional arguments.
    options.seqFileName = getArgumentValue(parser, 1);
    options.samFileName = getArgumentValue(parser, 2);
    options.witFileName = getArgumentValue(parser, 3);

    return 0;
}

// Entry point for the read mapper evaluation subprogram.
int evaluateReadMapperResult(Options<EvaluateResults> const & options)
{
    return 0;
}

#endif  // #ifndef APPS_RABEMA_EVALUATION_H_
