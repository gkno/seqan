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
   FUBerMarkRM -- WIT builder
  ==========================================================================*/

#include <algorithm>
#include <cstdlib>    // exit()

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/misc/misc_cmdparser.h>

#include "intervals.h"
#include "verification.h"
#include "wit_builder.h"
//#include "witio.h"


const char * kRevision = "0.0alpha";


// Define some return codes.
const int kRetOk = 0;       // OK, no errors.
const int kRetArgsErr = 1;  // Errors in arguments.
const int kRetIoErr = 2;    // I/O error, problem reading files.
const int kFatalErr = 3;    // Some other sort of fatal error.


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
    addHelpLine(parser, "hamming-weighted = Hamming distance, weighted by base call score");
    addHelpLine(parser, "edit             = Edit distance");
    addHelpLine(parser, "edit-weighted    = Edit distance, weighted by base call score");
    addOption(parser, CommandLineOption("o", "out-file", ("Path to the output file.  Use \"-\" for stdout, default is \"-\""), OptionType::String));
//     addOption(parser, CommandLineOption("rr", "reverse-complement-reads", "If set, the reads are reverse-complemented instead of the genome for finding hits on the reverse strand.  Default: not set.", OptionType::Boolean));
    
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
        getOptionValueLong(parser, "distance-function",
                           options.distanceFunction);
    if (isSetLong(parser, "out-file"))
        getOptionValueLong(parser, "out-file", outFile);
//     if (isSetLong(parser, "reverse-complement-reads"))
//         options.reverseComplementGenome = false;

    // Validate values.
    if (options.maxError < 0) {
        std::cerr << "ERROR: Invalid maximum error value: " <<
            options.maxError << std::endl;
        return kRetArgsErr;
    }
    if (not options.validDistanceFunction()) {      
      std::cerr << "ERROR: Invalid distance function: " <<
        options.distanceFunction << std::endl;
      return kRetArgsErr;
    }
    
    // Get positional arguments.
    options.referenceSeqFilename = getArgumentValue(parser, 0);
    options.perfectMapFilename = getArgumentValue(parser, 1);

    return kRetOk;
}



// Build intervals from the error curves.
template <typename TFragmentStore>
void intervalizeErrorCurves(String<WitRecord> & result,
                            TErrorCurves const & errorCurves,
                            TFragmentStore const & fragments,
                            Options const & options) {
    typedef typename TErrorCurves::const_iterator TErrorCurvesIter;
    for (TErrorCurvesIter it = errorCurves.begin(); it != errorCurves.end(); ++it) {
        size_t readId = it->first;
        TWeightedMatches const & matches = it->second;
        
        // Sort the matches.  Matches with high scores (negative score and
        // low absolute value) come first.
        TWeightedMatches sortedMatches(matches);
        std::sort(begin(sortedMatches, Standard()), end(sortedMatches, Standard()));
    
        // intervals[e] holds the intervals for error e of the current read.
        String<String<ContigInterval> > intervals;
        resize(intervals, options.maxError + 1);

        // Join the intervals stored in sortedMatches.
        //
        // The position of the previous match, so we can consider only the
        // ones with the smallest error.
        //
        // The following two vars should be != first pos and contigId.
        size_t previousPos = supremumValue<size_t>();
        size_t previousContigId = supremumValue<size_t>();
        typedef Iterator<TWeightedMatches>::Type TWeightedMatchesIter;
        for (TWeightedMatchesIter it = begin(sortedMatches);
             it != end(sortedMatches); ++it) {
            // Skip it if (it - 1) pointed to same pos (and must point to
            // one with smaller absolute distance.
            if (it->pos == previousPos and it->contigId == previousContigId)
                continue;
            // Consider all currently open intervals with a greater error
            // than the error in *it and extend them or create a new one.
            int error = abs(it->distance);
            SEQAN_ASSERT_LEQ(error, options.maxError);
            for (int e = error; e <= options.maxError; ++e) {
                // Handle base case of no open interval:  Create new one.
                if (length(intervals[e]) == 0) {
                    appendValue(intervals[e], ContigInterval(it->contigId, it->isForward, it->pos, it->pos));
                    continue;
                }
                ContigInterval &interval = back(intervals[e]);
                // Either extend the interval or create a new one.
                SEQAN_ASSERT(interval.last <= pos);
                if (interval.last + 1 == it->pos)
                    back(intervals[e]).last += 1;
                else
                    appendValue(intervals[e], ContigInterval(it->contigId, it->isForward, it->pos, it->pos));
            }
            // Book-keeping.
            previousPos = it->pos;
            previousContigId = it->contigId;
        }

        // Print the resulting intervals.
        typedef Iterator<String<String<ContigInterval> > >::Type TIntervalContainerIter;
        int distance = 0;
        for (TIntervalContainerIter it = begin(intervals);
             it != end(intervals); ++it, ++distance) {
            typedef Iterator<String<ContigInterval> >::Type TIntervalIter;
            for (TIntervalIter it2 = begin(*it); it2 != end(*it); ++it2)
                appendValue(result,
                            WitRecord(fragments.readNameStore[readId],
                                      distance, fragments.contigNameStore[it2->contigId],
                                      it2->isForward, it2->first, it2->last));
        }
    }
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
//     options.reverseComplementGenome = true;
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
    // TODO(holtgrew): CONTINUE to interpret options for weighted edit and hamming distance.
    if (options.distanceFunction == "edit")
        // TODO(holtgrew): Using non-read version of MyersUkkonen for now to check whether this causes problems with RazerS.
        matchesToErrorFunction(fragments, errorCurves, options, Myers<FindInfix>());
    else if (options.distanceFunction == "edit-weighted")
        matchesToErrorFunction(fragments, errorCurves, options, QualityDpSearch<FindInfix>());
    else if (options.distanceFunction == "hamming")
        matchesToErrorFunction(fragments, errorCurves, options, HammingSimple());
    else  // options.distanceFunction == "hamming-weighted"
        matchesToErrorFunction(fragments, errorCurves, options, HammingSimpleQuality());

    // =================================================================
    // Verify result when enabled.
    // =================================================================
    if (options.verify) {
        bool valid;
        // TODO(holtgrew): Using non-read version of MyersUkkonen for now to check whether this causes problems with RazerS.
        if (options.distanceFunction == "edit")
            valid = verifyMatchesToErrorFunctionResults(fragments, errorCurves, options, Myers<FindInfix>());
        else if (options.distanceFunction == "edit-weighted")
            valid = verifyMatchesToErrorFunctionResults(fragments, errorCurves, options, QualityDpSearch<FindInfix>());
        else if (options.distanceFunction == "hamming")
            valid = verifyMatchesToErrorFunctionResults(fragments, errorCurves, options, HammingSimple());
        else  // options.distanceFunction == "hamming-weighted"
            valid = verifyMatchesToErrorFunctionResults(fragments, errorCurves, options, HammingSimpleQuality());

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
