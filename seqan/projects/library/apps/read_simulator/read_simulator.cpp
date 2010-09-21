/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010
  
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
  A simple read simulator.

  Usage: read_simulator --help
         read_simulator illumina [options] source file
             Simulation of Illumina reads.
         read_simulator 454 [options] source file
             Simulation of 454 reads.
         read_simulator sanger [options] source file
             Simulation of Sanger reads.
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
#include "simulate_454.h"
#include "simulate_sanger.h"

using namespace seqan;

void printHelpGlobal() {
    std::cerr << "SeqAn read simulator" << std::endl
              << std::endl
              << "Usage: read_simulator illumina [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << "       read_simulator 454 [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << "       read_simulator sanger [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << std::endl
              << "Call with 'read_simulator READS-TYPE --help' to get detailed help." << std::endl;
}


int parseOptions(Options<Global> & options, const int argc, const char * argv[]) {
    if (argc == 1 ||
        (argc >= 2 && CharString(argv[1]) == "--help") ||
        (argc >= 2 && CharString(argv[1]) == "-h") ||
        (CharString(argv[1]) != "illumina" && CharString(argv[1]) != "454" &&
         CharString(argv[1]) != "sanger")) {
        printHelpGlobal();
        return 1;
    }

    if (CharString(argv[1]) == "illumina") {
        options.readsType = READS_TYPE_ILLUMINA;
    } else if (CharString(argv[1]) == "454") {
        options.readsType = READS_TYPE_454;
    } else if (CharString(argv[1]) == "sanger") {
        options.readsType = READS_TYPE_SANGER;
    } else {
        printHelpGlobal();
        return 1;
    }

    return 0;
}

int main(const int argc, const char * argv[]) {
    printDebugLevel(std::cerr);
    
    // Switch command type (which reads to simulate) or show global help.
    Options<Global> globalOptions;
    int ret = parseOptions(globalOptions, argc, argv);
    if (ret != 0)
        return ret;

    // Kick off read simulation, depending on the chosen command.
    if (globalOptions.readsType == READS_TYPE_ILLUMINA) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, IlluminaReads());
        Options<IlluminaReads> options;
        CharString referenceFilename;
        int ret = parseCommandLineAndCheck(options, referenceFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return simulateReads(options, referenceFilename, IlluminaReads());
    } else if (globalOptions.readsType == READS_TYPE_454) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, LS454Reads());
        Options<LS454Reads> options;
        CharString referenceFilename;
        int ret = parseCommandLineAndCheck(options, referenceFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return simulateReads(options, referenceFilename, LS454Reads());
    } else if (globalOptions.readsType == READS_TYPE_SANGER) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, SangerReads());
        Options<SangerReads> options;
        CharString referenceFilename;
        int ret = parseCommandLineAndCheck(options, referenceFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return simulateReads(options, referenceFilename, SangerReads());
    } else {
        SEQAN_ASSERT_FAIL("Invalid reads type!");
    }

    return 0;
}
