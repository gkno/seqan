/*==========================================================================
   SeqAn - The Library for Sequence Analysis
   http://www.seqan.de 
  ==========================================================================
   Copyright (C) 2010
  
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.
  
  ==========================================================================
   Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ==========================================================================
   Usage: wit_builder [options] <contigs.fasta> <golden.sam>

   Call "wit_builder --help" to see a complete list of options.

   This file creates a reference WIT file from a "golden" SAM result file
   by a read mapper with full sensitivity (such as RazerS).  This file
   only parses the command line arguments.  The real flesh of the program
   is in wit_builder.h.
  ==========================================================================*/

#include "wit_builder.h"

#include <seqan/misc/misc_cmdparser.h>

#include "return_codes.h"

const char * kRevision = "0.0alpha";

using namespace seqan;  // Remove some syntatic noise.

// Set up the CommandLineParser options with options.
void setUpCommandLineParser(CommandLineParser &parser) {
    // Add a version string.
    std::string versionString = "Golden standard generator version ";
    versionString += kRevision;
    addVersionLine(parser, versionString);

    // Add usage lines.
    addTitleLine(parser, "Golden standard generator.");
    addUsageLine(parser, "[OPTIONS] <REFERENCE SEQ> <PERFECT MAP>");

    // Set options.
    addOption(parser, CommandLineOption("e", "max-error-rate", "the maximal error in percent of read length, default: 0", OptionType::Double));
    addOption(parser, CommandLineOption("x", "verify", "Verify Result.  Default: true.", OptionType::Boolean));
    addOption(parser, CommandLineOption("d", "distance-function", ("the distance function to use, default: hamming"), OptionType::String | OptionType::Label));
    addHelpLine(parser, "hamming          = Hamming distance");
    addHelpLine(parser, "edit             = Edit distance");
    addOption(parser, CommandLineOption("o", "out-file", ("Path to the output file.  Use \"-\" for stdout, default is \"-\""), OptionType::String));
    addOption(parser, CommandLineOption("mN" , "match-N", "If specified, N characters match all other characters.  The default is for them to mismatch with CGAT.", OptionType::Boolean));
    
    // We require 4 command line options.
    requiredArguments(parser, 2);
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
    if (isSetLong(parser, "verify"))
        options.verify = true;
    if (isSetLong(parser, "max-error-rate"))
        getOptionValueLong(parser, "max-error-rate", options.maxError);
    if (isSetLong(parser, "distance-function"))
        getOptionValueLong(parser, "distance-function", options.distanceFunction);
    if (isSetLong(parser, "match-N"))
        options.matchN = true;
    if (isSetLong(parser, "out-file"))
        getOptionValueLong(parser, "out-file", outFile);

    // Validate values.
    if (options.maxError < 0) {
        std::cerr << "ERROR: Invalid maximum error value: " << options.maxError << std::endl;
        return kRetArgsErr;
    }
    if (not options.validDistanceFunction()) {      
      std::cerr << "ERROR: Invalid distance function: " << options.distanceFunction << std::endl;
      return kRetArgsErr;
    }
    
    // Get positional arguments.
    options.referenceSeqFilename = getArgumentValue(parser, 0);
    options.perfectMapFilename = getArgumentValue(parser, 1);

    return kRetOk;
}


// Load the FASTA and SAM file with the contigs and the mapped result
// into the given FragmentStore.
//
// Returns 0 on success and an error code otherwise.
template <typename TFragmentStore>
int loadContigsAndSamFile(TFragmentStore & fragmentStore,
                          CharString const & contigsFilename,
                          CharString const & samFilename) {
    // Load contigs.
    if (!loadContigs(fragmentStore, contigsFilename)) {
        std::cerr << "Could not read contigs." << std::endl;
        return kRetIoErr;
    }
    
    // Load SAM file.
    {
        std::fstream fstrm(toCString(samFilename),
                           std::ios_base::in | std::ios_base::binary);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open SAM file." << std::endl;
            return kRetIoErr;
        }
        read(fstrm, fragmentStore, SAM());
    }

    return kRetOk;
}


int main(int argc, const char *argv[]) {
    // =================================================================
    // Initialize Options object with defaults.
    // =================================================================
    Options options;
    options.verify = false;
    options.maxError = 0;
    options.distanceFunction = "edit";
    options.matchN = false;
    CharString outFile = "-";

    seqan::printDebugLevel(std::cerr);

    // =================================================================
    // Setup the parser, parse command line and return if an error occured.
    // =================================================================
    CommandLineParser parser;
    setUpCommandLineParser(parser);
    int ret = parseCommandLineAndCheck(options, parser, outFile, argc, argv);
    if (ret != kRetOk)
        return ret;

    // =================================================================
    // Load contigs and SAM file.
    // =================================================================
    typedef FragmentStore<> TFragmentStore;
    TFragmentStore fragments;
    ret = loadContigsAndSamFile(fragments, options.referenceSeqFilename, options.perfectMapFilename);
    if (ret != kRetOk)
        return ret;

    // =================================================================
    // Build point-wise error curve.
    // =================================================================
    TErrorCurves errorCurves;
    if (options.distanceFunction == "edit")
        // TODO(holtgrew): Using non-read version of MyersUkkonen for now to check whether this causes problems with RazerS.
        matchesToErrorFunction(fragments, errorCurves, options, Myers<FindInfix>());
    else // options.distanceFunction == "hamming"
        matchesToErrorFunction(fragments, errorCurves, options, HammingSimple());

    // =================================================================
    // Verify result when enabled.
    // =================================================================
    if (options.verify) {
        bool valid;
        // TODO(holtgrew): Using non-read version of MyersUkkonen for now to check whether this causes problems with RazerS.
        if (options.distanceFunction == "edit")
            valid = verifyMatchesToErrorFunctionResults(fragments, errorCurves, options, Myers<FindInfix>());
        else // options.distanceFunction == "hamming"
            valid = verifyMatchesToErrorFunctionResults(fragments, errorCurves, options, HammingSimple());

        if (!valid) {
            std::cerr << "ERROR: Result does not validate!" << std::endl;
            return kFatalErr;
        }
        // else...
        std::cerr << "Result DID validate!" << std::endl;
    }

    // =================================================================
    // Convert points in error curves to intervals and write them to
    // stdout or a file.
    // =================================================================
    String<WitRecord> witRecords;
    typedef Iterator<String<WitRecord>, Standard>::Type TWitRecordIterator;
    intervalizeErrorCurves(witRecords, errorCurves, fragments, options);
    // The two alternatives are equivalent after opening the file.
    if (outFile == "-") {
        writeWitHeader(std::cout);
        writeWitComment(std::cout, WIT_COLUMN_NAMES);
        for (TWitRecordIterator it = begin(witRecords, Standard());
             it != end(witRecords, Standard()); ++it)
            std::cout << *it << std::endl;
    } else {
        std::fstream fstrm(toCString(outFile), std::ios_base::out);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open out file \"" << outFile << "\""
                      << std::endl;
            return kRetIoErr;
        }
        writeWitHeader(fstrm);
        writeWitComment(fstrm, WIT_COLUMN_NAMES);
        for (TWitRecordIterator it = begin(witRecords, Standard());
             it != end(witRecords, Standard()); ++it)
            fstrm << *it << std::endl;
    }

    return kRetOk;
}
