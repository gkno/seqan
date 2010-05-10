#include <iostream>
#include <fstream>
#include <cstdlib>
#include <map>
#include <cmath>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/score.h>
#include <seqan/misc/misc_cmdparser.h>

#include "return_codes.h"
#include "wit_store.h"
#include "reweight_wit.h"

const char * kRevision = "0.0alpha";

using namespace seqan;


// Set up the CommandLineParser options with options.
void setUpCommandLineParser(CommandLineParser &parser) {
    // Add a version string.
    std::string versionString = "WIT Reweighting Program  ";
    versionString += kRevision;
    addVersionLine(parser, versionString);

    // Add usage lines.
    addTitleLine(parser, "WIT Reweighting Program.");
    addUsageLine(parser, "[OPTIONS] <REFERENCE SEQ> <READS FILE> <INPUT WIT FILE>");

    // Set options.
    addOption(parser, CommandLineOption("d", "distance-function", ("the distance function to use, default: hamming"), OptionType::String | OptionType::Label));
    addHelpLine(parser, "hamming          = Hamming distance");
    addHelpLine(parser, "edit             = Edit distance");
    addOption(parser, CommandLineOption("o", "out-file", ("Path to the output file.  Use \"-\" for stdout, default is \"-\""), OptionType::String));
    addOption(parser, CommandLineOption("e", "max-weighted-error", "Maximal weighted error.  Default is 0.", OptionType::Int));

    // We require 3 command line options.
    requiredArguments(parser, 3);
}


// Actually parse the arguments.  Return kRetOk or the return value
// for main.
int parseCommandLineAndCheck(Options &options,
                             CommandLineParser &parser,
                             int argc, const char *argv[]) {
    // Show short help on invalid arguments and long help if help
    // argument was given.
    if (not parse(parser, argc, argv)) {
        shortHelp(parser, std::cerr);
        return kRetArgsErr;
    } else if (isSetShort(parser, 'h')) {
        exit(kRetOk);
    }

    // Get arguments.
    if (isSetLong(parser, "out-file"))
        getOptionValueLong(parser, "out-file", options.outWitFilename);
    if (isSetLong(parser, "distance-function"))
        getOptionValueLong(parser, "distance-function", options.distanceFunction);
    if (isSetLong(parser, "max-weighted-error"))
        getOptionValueLong(parser, "max-weighted-error", options.maxWeightedError);

    if (not options.validDistanceFunction()) {      
      std::cerr << "ERROR: Invalid distance function: " <<
        options.distanceFunction << std::endl;
      return kRetArgsErr;
    }
    
    // Get positional arguments.
    options.genomeFilename = getArgumentValue(parser, 0);
    options.readsFilename = getArgumentValue(parser, 1);
    options.inputWitFilename = getArgumentValue(parser, 2);

    return kRetOk;
}


int main(int argc, const char *argv[]) {
    // =================================================================
    // Initialize Options object with defaults.
    // =================================================================
    Options options;
    options.maxWeightedError = 0;
    options.outWitFilename = "-";
    options.distanceFunction = "edit";

    seqan::printDebugLevel(std::cerr);

    // =================================================================
    // Setup the parser, parse command line and return if an error occured.
    // =================================================================
    CommandLineParser parser;
    setUpCommandLineParser(parser);
    int ret = parseCommandLineAndCheck(options, parser, argc, argv);
    if (ret != kRetOk)
        return ret;

    // =================================================================
    // Read In Reference Sequence, Reads and WIT file.
    // =================================================================
    std::cerr << "Loading contigs from " << options.genomeFilename << std::endl;
    StringSet<String<Dna5> > contigs;
    StringSet<CharString> contigNames;
    loadFileIntoStringSet(contigs, contigNames, options.genomeFilename);
    std::cerr << "Loading reads from " << options.readsFilename << std::endl;
    StringSet<String<Dna5Q> > reads;
    StringSet<CharString> readNames;
    loadFileIntoStringSet(reads, readNames, options.readsFilename);
    std::cerr << "Loading intervals from " << options.inputWitFilename << std::endl;
    WitStore witStore;
    loadWitFile(witStore, readNames, contigNames, options.inputWitFilename);

    // =================================================================
    // Reweight Weighted Intervals.
    // =================================================================
    if (options.distanceFunction == "edit")
        reweightWitStore(witStore, contigs, reads, options, Myers<FindInfix>());
    else  // options.distanceFunction == "hamming"
        reweightWitStore(witStore, contigs, reads, options, HammingSimple());

    // =================================================================
    // Write Output
    // =================================================================
    if (options.outWitFilename == "-") {
        writeWitFile(std::cout, witStore);
    } else {
        std::fstream fstrm(toCString(options.outWitFilename), std::ios_base::out | std::ios_base::binary);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open WIT file " << options.outWitFilename << std::endl;
            exit(1);
        }
        writeWitFile(fstrm, witStore);
    }

    return kRetOk;
}
