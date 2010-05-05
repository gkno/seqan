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

#include "intervals.h"
#include "witio.h"
#include "find_myers_ukkonen_reads.h"

// TODO(holtgrew): Hamming distance not implemented.
// TODO(holtgrew): We do not need to know the distance type here?!

using namespace seqan;  // Remove some syntatic noise.


// The revision string of the program.
const char * kRevision = "0.0alpha";

// Define some return codes.
const int kRetOk = 0;       // OK, no errors.
const int kRetArgsErr = 1;  // Errors in arguments.
const int kRetIoErr = 2;    // I/O error, problem reading files.
const int kFatalErr = 3;    // Some other sort of fatal error.


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
    addHelpLine(parser, "hamming-weighted = Hamming distance, weighted by quality value");
    addHelpLine(parser, "edit             = Edit distance");
    addHelpLine(parser, "edit-weighted    = Edit distance, weighted by quality value");
    addOption(parser, CommandLineOption("o", "out-file", "Path to the output file.  Use \"-\" for stdout, default is \"-\"", OptionType::String));
    addOption(parser, CommandLineOption("mi", "show-missed-intervals", "the missed intervals are printed to stderr for debugging if this option is set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("hi", "show-hit-intervals", "the hit intervals are printed to stderr for debugging if this option is set.", OptionType::Boolean));
    addOption(parser, CommandLineOption("st", "show-try-hit-intervals", "The last positions tried to hit against an intervals are printd to stderr if set.", OptionType::Boolean));
    
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
    if (isSetLong(parser, "show-missed-intervals"))
        options.showMissedIntervals = true;
    if (isSetLong(parser, "show-hit-intervals"))
        options.showHitIntervals = true;
    if (isSetLong(parser, "show-try-hit-intervals"))
        options.showTryHitIntervals = true;
    if (isSetLong(parser, "max-error-rate"))
        getOptionValueLong(parser, "max-error-rate", options.maxError);
    if (isSetLong(parser, "distance-function"))
        getOptionValueLong(parser, "distance-function",
                           options.distanceFunction);
    if (isSetLong(parser, "out-file"))
        getOptionValueLong(parser, "out-file", outFile);

    // Validate values.
    if (options.maxError < 0) {
        std::cerr << "ERROR: Invalid maximum error value: " <<
            options.maxError << std::endl;
        return kRetArgsErr;
    }
    if (not options.validDistanceFunction())     {      
      std::cerr << "ERROR: Invalid distance function: " <<
        options.distanceFunction << std::endl;
      return kRetArgsErr;
    }
    
    // Get positional arguments.
    options.seqFileName = getArgumentValue(parser, 0);
    options.samFileName = getArgumentValue(parser, 1);
    options.witFileName = getArgumentValue(parser, 2);

    return kRetOk;
}


// Read the WIT file at the given filename.
//
// Return return code for the program.
//
// filename -- Path to file to read.
// fragments -- The FragmentStore<> to use for read/contig identification
//              from string name.
// weightedMatches  -- The resulting weighted matches.
// ignoredReadCount -- Output parameter, number of reads in the WIT file
//                     that could not be found in the fragment store.
template <typename TFragmentStore>
int readWitFile(const CharString &filename, const TFragmentStore &fragments,
                String<WitRecord> &witRecords, size_t &ignoredReadCount) {
    // Initialize output parameters.
    clear(witRecords);
    ignoredReadCount = 0;

    // Open the file.
    std::fstream fstrm(toCString(filename), std::ios_base::in);
    if (not fstrm.is_open()) {
        std::cerr << "Could not open WIT file." << std::endl;
        return kRetIoErr;
    }

    // Build map from read name to read id.
    typedef typename TFragmentStore::TReadStore TReadStore;
    typedef typename Position<TReadStore>::Type TReadStorePos;
    std::map<CharString, TReadStorePos> readNameToId;
    for (TReadStorePos readId = 0; readId < length(fragments.readNameStore); ++readId) {
        readNameToId[fragments.readNameStore[readId]] = readId;
        // std::cout << readId << ": " << fragments.readNameStore[readId] << std::endl;
    }

    // Build map from contig name to contig id.
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Position<TContigStore>::Type TContigStorePos;
    std::map<CharString, TContigStorePos> contigNameToId;
    for (TContigStorePos contigId = 0; contigId < length(fragments.contigNameStore); ++contigId) {
        contigNameToId[fragments.contigNameStore[contigId]] = contigId;
    }

    // Read header.
    char c;
    readWitHeader(fstrm, c);
    // Read all records.
    bool wasRecord;
    std::set<CharString> ignoredReadNames;
    while (true) {  // Will break, depending on wasRecord.
        WitRecord witRecord;
        wasRecord = readWitRecord(fstrm, witRecord, c);
        if (not wasRecord)
            break;  // out of loop.

        // If we have a read with an id in the WIT file that is not in
        // the SAM file then we have ignore the WIT entry.  We record
        // all the ignored read names and return a count of them in an
        // output parameter so they can be tallied as "not hit" later.
        if (readNameToId.find(witRecord.readName) == readNameToId.end()) {
            ignoredReadNames.insert(witRecord.readName);
            continue;  // Next record.
        }

        // Resolve contig and read name to id.
        SEQAN_ASSERT_TRUE(contigNameToId.find(witRecord.contigName) != contigNameToId.end());
        witRecord.contigId = contigNameToId[witRecord.contigName];
        SEQAN_ASSERT_TRUE(contigNameToId.find(witRecord.readName) != readNameToId.end());
        witRecord.readId = readNameToId[witRecord.readName];

        // Append record to the output parameter witRecords.
        appendValue(witRecords, witRecord);
//         std::cerr << "append(witRecords, " << witRecord << ")" << std::endl;
    }

    // Write number of ignored reads out.
    ignoredReadCount = ignoredReadNames.size();

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
    std::cerr << "Reading WIT file " << options.witFileName << " ..." << std::endl;
    startTime = sysTime();
    typedef String<WitRecord> TWitRecords;
    TWitRecords witRecords;
    size_t ignoredReadCount;
    ret = readWitFile(options.witFileName, fragments, witRecords,
                      ignoredReadCount);
    if (ret != kRetOk)
        return ret;
    // Sort WIT records by read id.
    std::sort(begin(witRecords), end(witRecords));
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;

    // =================================================================
    // Compare The SAM Hits Against WIT Intervals.
    // =================================================================
    std::cerr << "Compare reader hits from SAM file against WIT file." << std::endl;
    startTime = sysTime();
    typedef Position<TWitRecords>::Type TPos;
    size_t relevantIntervalCount = 0;
    size_t foundIntervalCount = 0;
    if (options.distanceFunction == "edit")
        compareAlignedReadsToReference(options, fragments, witRecords, foundIntervalCount, relevantIntervalCount, Myers<FindInfix>());
        // TODO(holtgrew): Using non-read version of MyersUkkonen for now to check whether this causes problems with RazerS.
//                                              MyersUkkonenReads());
    else if (options.distanceFunction == "edit-weighted")
        compareAlignedReadsToReference(options, fragments, witRecords, foundIntervalCount, relevantIntervalCount, QualityDpSearch<FindInfix>());
    else if (options.distanceFunction == "hamming")
        compareAlignedReadsToReference(options, fragments, witRecords, foundIntervalCount, relevantIntervalCount, HammingSimple());
    else // options.distanceFunction == "hamming-weighted"
        compareAlignedReadsToReference(options, fragments, witRecords, foundIntervalCount, relevantIntervalCount, HammingSimpleQuality());
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;
    

    // =================================================================
    // Write Output.
    // =================================================================
    std::cerr << "Write output..." << std::endl;
    startTime = sysTime();
    // The output consists of one line that describes the total and
    // found intervals as a JSON record with the entries
    // "total_intervals" and "found_itervals".
    if (outFile == "-") {
        // Print to stdout.
        std::cout << "{\"total_intervals\": " << relevantIntervalCount
                  << ", \"found_intervals\": " << foundIntervalCount
                  << "}" << std::endl;
    } else {
        // Write output to file.
        std::fstream fstrm(toCString(outFile), std::ios_base::out);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open output JSON file." << std::endl;
            return kRetIoErr;
        }
        fstrm << "{\"total_intervals\": " << relevantIntervalCount + ignoredReadCount
              << ", \"found_intervals\": " << foundIntervalCount
              << "}" << std::endl;
    }
    std::cerr << "Took " << sysTime() - startTime << " s" << std::endl;
    return kRetOk;
}
