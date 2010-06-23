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
  Author: Anne-Katrin Emde <anne-katrin.emde@fu-berlin.de>
  ===========================================================================
  A simple read simulator.

  Usage: read_simulator --help
         read_simulator illumina [options] source file
             Simulation of Illumina reads.
  ===========================================================================
  This file only contains code to parse command line parameters and delegates
  the actual simulation to the simulate_*.h headers.  The Illumina simulation
  code is based on "readsim" by Anne-Katrin Emde, which is itself based on
  the read simulator by Tobias Rausch.
  ===========================================================================*/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/misc/misc_random.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/store.h>

#include "read_simulator.h"
#include "simulate_illumina.h"

using namespace seqan;

void printHelpGlobal() {
    std::cerr << "SeqAn read simulator" << std::endl
              << std::endl
              << "Usage: read_simulator illumina [OPTIONS] SEQUENCE.fasta" << std::endl
              << std::endl
              << "Call with 'read_simulator READS-TYPE --help' to get detailed help." << std::endl;
}


int parseOptions(GlobalOptions & options, const int argc, const char * argv[]) {
    if (argc == 1 ||
        (argc >= 2 && CharString(argv[1]) == "--help") ||
        (argc >= 2 && CharString(argv[1]) == "-h") ||
        (CharString(argv[1]) != "illumina")) {
        printHelpGlobal();
        return 1;
    }

    if (CharString(argv[1]) == "illumina") {
        options.readsType = GlobalOptions::READS_TYPE_ILLUMINA;
    } else {
        printHelpGlobal();
        return 1;
    }

    return 0;
}


void setUpCommandLineParser(CommandLineParser & parser,
                            IlluminaReads const &) {
    addVersionLine(parser, "SeqAn read simulator 0.1");

    addTitleLine(parser, "SeqAn read simulator");
    addUsageLine(parser, "illumina [OPTIONS] SEQUENCE");
    addLine(parser, "");
    addLine(parser, "Use 'random' for the SEQUENCE file name to generate it randomly.");

    addSection(parser, "Main Options");
    
    addOption(parser, CommandLineOption("s",  "seed", "The seed for RNG.  Default: 0.", OptionType::Integer | OptionType::Label));
    addOption(parser, CommandLineOption("n",  "length", "The length of the reads to simulate.  Default: 36.", OptionType::Integer | OptionType::Label));
    addHelpLine(parser, "All resulting reads will have the same length.");
    addOption(parser, CommandLineOption("N",  "num-reads", "Number of reads (or mate pairs) to simulate.  Default: 1000.", OptionType::Integer));
    addOption(parser, CommandLineOption("sn", "source-length", "Length of random source sequence.  Default: 1,000,000.", OptionType::Integer));
    addOption(parser, CommandLineOption("e",  "errors", "Maximal number of errors per read.  Default: 4.", OptionType::Integer));
    addOption(parser, CommandLineOption("f",  "forward-only", "Simulate from forward strand only.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("r",  "reverse-only", "Simulate from reverse strand only.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("o",  "output-file", "Write results to PARAM.fasta file instead of SEQUENCE.reads.fasta.  Default: \"\".", OptionType::String));
    addOption(parser, CommandLineOption("sq", "simulate-qualities", "Simulate qualities, generate FASTQ instead of FASTA.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("v", "verbose", "Verbosity mode.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("vv", "very-verbose", "High verbosity mode, implies verbosity mode.  Default: false.", OptionType::Bool));

    addSection(parser, "Mate-Pair Options");

    addOption(parser, CommandLineOption("ll", "library-length", "Mate-pair library length.  Default: 1000.", OptionType::Integer));
    addOption(parser, CommandLineOption("le", "library-error", "Mate-pair library tolerance.  Default: 100.", OptionType::Integer));
    addOption(parser, CommandLineOption("mp", "mate-pairs", "Enable mate pair simulation.  Default: false.", OptionType::Bool));

    addSection(parser, "Probability Options");

    addOption(parser, CommandLineOption("d",  "error-distribution", "File containing mismatch qualities.  If left blank, defaults are used.  Defaults are available for n = 36, 50, 100.  Default: \"\".", OptionType::String));
    addOption(parser, CommandLineOption("pi", "prob-insert", "Probability of an insertion.  Default: 0.01.", OptionType::Double));
    addOption(parser, CommandLineOption("pd", "prob-delete", "Probability of a deletion.  Default: 0.01.", OptionType::Double));

    addSection(parser, "Haplotype Options");

    addOption(parser, CommandLineOption("hn", "num-haplotypes", "Number of haplotypes to simulate.  Default: 1.", OptionType::Integer));
    addOption(parser, CommandLineOption("hs", "haplotype-snp-rate", "Haplotype SNP rate.  Default: 0.01.", OptionType::Double));
    addOption(parser, CommandLineOption("hi", "haplotype-indel-rate", "Haplotype indel rate.  Default: 0.01.", OptionType::Double));
    addOption(parser, CommandLineOption("hm", "haplotype-indel-range-min", "Haplotype indel size min.  Default: 4.", OptionType::Integer));
    addOption(parser, CommandLineOption("hM", "haplotype-indel-range-max", "Haplotype indel size max.  Default: 6.", OptionType::Integer));

    // Need "illumina" and SEQUENCE.fasta.
    requiredArguments(parser, 2);
}


int parseCommandLineAndCheck(IlluminaOptions & options,
                             CharString & referenceFilename,
                             CommandLineParser & parser,
                             const int argc,
                             const char * argv[],
                             IlluminaReads const &) {
    if (!parse(parser, argc, argv)) {
        if (!isSetShort(parser, "h"))
            shortHelp(parser, std::cerr);
        return 1;
    }
    if (isSetShort(parser, "h")) {
        options.showHelp = true;
        return 0;
    }

    if (isSetLong(parser, "length"))
        getOptionValueLong(parser, "length", options.readLength);
    if (isSetLong(parser, "seed"))
        getOptionValueLong(parser, "seed", options.seed);
    if (isSetLong(parser, "num-reads"))
        getOptionValueLong(parser, "num-reads", options.numReads);
    if (isSetLong(parser, "source-length"))
        getOptionValueLong(parser, "source-length", options.randomSourceLength);
    if (isSetLong(parser, "errors"))
        getOptionValueLong(parser, "errors", options.maxErrorsPerRead);
    if (isSetLong(parser, "forward-only"))
        options.onlyForward = true;
    if (isSetLong(parser, "reverse-only"))
        options.onlyReverse = true;
    if (isSetLong(parser, "output-file"))
        getOptionValueLong(parser, "output-file", options.outputFile);
    if (isSetLong(parser, "simulate-qualities"))
        options.simulateQualities = true;
    if (isSetLong(parser, "verbose"))
        options.verbose = true;
    if (isSetLong(parser, "very-verbose")) {
        options.verbose = true;
        options.veryVerbose = true;
    }

    if (isSetLong(parser, "library-length"))
        getOptionValueLong(parser, "library-length", options.libraryLength);
    if (isSetLong(parser, "library-error"))
        getOptionValueLong(parser, "library-error", options.libraryLengthError);
    if (isSetLong(parser, "mate-pairs"))
        options.generateMatePairs = true;

    if (isSetLong(parser, "error-distribution"))
        getOptionValueLong(parser, "error-distribution", options.errorDistributionFile);
    if (isSetLong(parser, "prob-insert"))
        getOptionValueLong(parser, "prob-insert", options.probabilityInsert);
    if (isSetLong(parser, "prob-delete"))
        getOptionValueLong(parser, "prob-delete", options.probabilityDelete);

    if (isSetLong(parser, "num-haplotypes"))
        getOptionValueLong(parser, "num-haplotypes", options.numHaplotypes);
    if (isSetLong(parser, "haplotype-snp-rate"))
        getOptionValueLong(parser, "haplotype-snp-rate", options.haplotypeSnpRate);
    if (isSetLong(parser, "haplotype-indel-rate"))
        getOptionValueLong(parser, "haplotype-indel-rate", options.haplotypeIndelRate);
    if (isSetLong(parser, "haplotype-indel-range-min"))
        getOptionValueLong(parser, "haplotype-indel-range-min", options.haplotypeIndelRangeMin);
    if (isSetLong(parser, "haplotype-indel-range-max"))
        getOptionValueLong(parser, "haplotype-indel-range-max", options.haplotypeIndelRangeMax);


    // First argument is "illumina", second one name of reference file.
    referenceFilename = getArgumentValue(parser, 1);

    if (referenceFilename == "random")
        options.useRandomSequence = true;

    return 0;
}


int main(const int argc, const char * argv[]) {
    // Switch command type (which reads to simulate) or show global help.
    GlobalOptions globalOptions;
    int ret = parseOptions(globalOptions, argc, argv);
    if (ret != 0)
        return ret;

    // Kick off read simulation, depending on the chosen command.
    if (globalOptions.readsType == GlobalOptions::READS_TYPE_ILLUMINA) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, IlluminaReads());
        IlluminaOptions options;
        CharString referenceFilename;
        int ret = parseCommandLineAndCheck(options, referenceFilename, parser, argc, argv, IlluminaReads());
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return simulateReads(options, referenceFilename, IlluminaReads());
    } else {
        SEQAN_ASSERT_FAIL("Invalid reads type!");
    }

    return 0;
}
