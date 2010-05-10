/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2007
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  FUBerMarkRM -- SAM against WIT checker

  Compares the SAM file yielded by a read mapper against a golden standard
  WIT file.
  ===========================================================================*/

// Uncomment to disable debugging.
// #define SEQAN_ENABLE_DEBUG 0

#include "compare_sam_wit.h"

#include <cassert>
#include <cstdlib>    // exit()

#include <algorithm>  // sort()
#include <map>
#include <set>

#include <seqan/basic.h>
#include <seqan/store.h>

#include "wit_store.h"
#include "return_codes.h"
#include "find_myers_ukkonen_reads.h"

using namespace seqan;  // Remove some syntatic noise.


// The revision string of the program.
const char * kRevision = "0.0alpha";


// Set up the CommandLineParser options with options.
void setUpCommandLineParser(CommandLineParser &parser) {
    // Add a version string.
    std::string versionString = "compare_sam_wit ";
    versionString += kRevision;
    addVersionLine(parser, versionString);

    // Add usage lines.
    addTitleLine(parser, "Compare SAM hits against WIT file.");
    addUsageLine(parser, "[OPTIONS] <SEQUENCE FILE> <SAM FILE> <WIT FILE>");

    // Set options.
    addOption(parser, CommandLineOption("e", "max-error-rate", "the maximal error rate in percent, default: 0", OptionType::Int));
    addOption(parser, CommandLineOption("d", "distance-function", "the distance function to use, default: hamming", OptionType::String | OptionType::Label));
    addHelpLine(parser, "hamming          = Hamming distance");
    addHelpLine(parser, "edit             = Edit distance");
    addOption(parser, CommandLineOption("o", "out-file", "Path to the output file.  Use \"-\" for stdout, default is \"-\"", OptionType::String));
    addOption(parser, CommandLineOption("mi", "show-missed-intervals", "the missed intervals are printed to stderr for debugging if this option is set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("hi", "show-hit-intervals", "the hit intervals are printed to stderr for debugging if this option is set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("st", "show-try-hit-intervals", "The last positions tried to hit against an intervals are printd to stderr if set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("mN", "match-N", "If set, N matches all as a wildcard character, otherwise it never matches.", OptionType::Boolean));
    addOption(parser, CommandLineOption("wm", "weighted-distances", "If set, use weighted distances instead of unit ones.", OptionType::Boolean));
    
    // We require 4 command line options.
    requiredArguments(parser, 3);
}


// Actually parse the arguments.  Return kRetOk or the return value
// for main.
int parseCommandLineAndCheck(Options &options,
                             CommandLineParser &parser,
                             CharString &outFile,
                             const int argc, const char *argv[]) {
    // Show short help on invalid arguments and long help if help
    // argument was given.
    if (not parse(parser, argc, argv)) {
        shortHelp(parser, std::cerr);
        return kRetArgsErr;
    } else if (isSetShort(parser, 'h')) {
        exit(kRetOk);
    }

    // Get arguments.
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
    if (isSetLong(parser, "max-error-rate"))
        getOptionValueLong(parser, "max-error-rate", options.maxError);
    if (isSetLong(parser, "distance-function"))
        getOptionValueLong(parser, "distance-function", options.distanceFunction);
    if (isSetLong(parser, "out-file"))
        getOptionValueLong(parser, "out-file", outFile);

    // Validate values.
    if (options.maxError < 0) {
        std::cerr << "ERROR: Invalid maximum error value: " << options.maxError << std::endl;
        return kRetArgsErr;
    }
    if (not options.validDistanceFunction())     {      
      std::cerr << "ERROR: Invalid distance function: " << options.distanceFunction << std::endl;
      return kRetArgsErr;
    }
    
    // Get positional arguments.
    options.seqFileName = getArgumentValue(parser, 0);
    options.samFileName = getArgumentValue(parser, 1);
    options.witFileName = getArgumentValue(parser, 2);

    return kRetOk;
}


int main(int argc, const char *argv[]) {
    // =================================================================
    // Parse Options.
    // =================================================================
    // Initialize Options object with defaults.
    Options options;
    options.maxError = 0;
    options.distanceFunction = "edit";
    options.showMissedIntervals = false;
    options.showHitIntervals = false;
    options.showTryHitIntervals = false;
    options.matchN = false;
    options.weightedDistances = false;
    CharString outFile = "-";

    // Setup the parser, parse command line and return if an error occured.
    CommandLineParser parser;
    setUpCommandLineParser(parser);
    int ret = parseCommandLineAndCheck(options, parser, outFile, argc, argv);
    if (ret != kRetOk)
        return ret;

    // =================================================================
    // Load FASTA Sequence And SAM File Into FragmentStore.
    // =================================================================
    typedef FragmentStore<> TFragmentStore;
    TFragmentStore fragments;

    // Load Contigs.
    double startTime = sysTime();
    std::cerr << "Reading FASTA contigs sequence file " << options.seqFileName << " ..." << std::endl;
    if (!loadContigs(fragments, options.seqFileName)) {
        std::cerr << "Could not read contigs." << std::endl;
        return kRetIoErr;
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // Load SAM File.
    std::cerr << "Reading SAM file file " << options.samFileName << " ..." << std::endl;
    startTime = sysTime();
    {
        std::fstream fstrm(toCString(options.samFileName),
                           std::ios_base::in | std::ios_base::binary);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open SAM file." << std::endl;
            return kRetIoErr;
        }
        read(fstrm, fragments, SAM());
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Load WIT file.
    // =================================================================
    std::cerr << "Loading intervals from " << options.witFileName << std::endl;
    startTime = sysTime();
    WitStore witStore;
    loadWitFile(witStore, fragments.readNameStore, fragments.contigNameStore, options.witFileName);
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Compare The SAM Hits Against WIT Intervals.
    // =================================================================
    std::cerr << "Compare reader hits from SAM file against WIT file." << std::endl;
    startTime = sysTime();
    typedef Position<WitStore::TIntervalStore>::Type TPos;
    ComparisonResult result;
    if (options.distanceFunction == "edit")
        // TODO(holtgrew): Switch to "reads version" (force alignment of last characters) after RazerS can do this, too.  Tag is MyersUkkonenReads().
        compareAlignedReadsToReference(options, fragments, witStore, result, Myers<FindInfix>());
    else  // options.distanceFunction == "hamming"
        compareAlignedReadsToReference(options, fragments, witStore, result, HammingSimple());
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;
    

    // =================================================================
    // Write Output.
    // =================================================================
    std::cerr << "Write output..." << std::endl;
    startTime = sysTime();
    // The output consists of one line that describes the total and
    // found intervals as a JSON record with the entries
    // "total_intervals", "found_itervals", "superflous_intervals",
    // "additional_intervals".
    if (outFile == "-") {
        // Print to stdout.
        std::cout << result << std::endl;
    } else {
        // Write output to file.
        std::fstream fstrm(toCString(outFile), std::ios_base::out);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open output JSON file." << std::endl;
            return kRetIoErr;
        }
        fstrm << result << std::endl;
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    return kRetOk;
}
