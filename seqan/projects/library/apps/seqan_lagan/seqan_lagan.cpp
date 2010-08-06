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
  An implementation of LAGAN in SeqAn using the new seeds module.
 ===========================================================================*/

#include <seqan/basic.h>
#include <seqan/misc/misc_cmdparser.h>

#include "seqan_lagan.h"
#include "seqan_lagan_lagan.h"

void printHelpGlobal() {
    std::cerr << "SeqAn::LAGAN" << std::endl
              << std::endl
              << "Computes a pairwise sequence alignment." << std::endl
              << std::endl
              << "Usage: seqan_lagan [-h, --help]" << std::endl
              << "       seqan_lagan dp LEFT.fasta RIGHT.fasta" << std::endl
              << "       seqan_lagan lagan LEFT.fasta RIGHT.fasta" << std::endl;
}


int parseOptions(Options<Global> & options, const int argc, const char * argv[])
{
    // No argument, --help, -h or invalid command.
    if (argc == 1 ||
        (argc >= 2 && CharString(argv[1]) == "--help") ||
        (argc >= 2 && CharString(argv[1]) == "-h") ||
        (CharString(argv[1]) != "lagan" && CharString(argv[1]) != "dp")) {
        printHelpGlobal();
        return 1;
    }

    // Parse command.
    if (CharString(argv[1]) == "lagan") {
        options.command = COMMAND_LAGAN;
    } else if (CharString(argv[1]) == "dp") {
        options.command = COMMAND_CLASSIC_DP;
    }

    return 0;
}


int main(int argc, char const * argv[])
{
    Options<Global> globalOptions;
    int ret = parseOptions(globalOptions, argc, argv);
    if (ret != 0)
        return ret;

    if (globalOptions.command == COMMAND_LAGAN) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, Lagan());
        Options<Lagan> options;
        CharString leftFilename;
        CharString rightFilename;
        int ret = parseCommandLineAndCheck(options, leftFilename, rightFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return executeParwiseAlignmentCommand(globalOptions, options, leftFilename, rightFilename, Lagan());
    } else if (globalOptions.command == COMMAND_CLASSIC_DP) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, ClassicDP());
        Options<ClassicDP> options;
        CharString leftFilename;
        CharString rightFilename;
        int ret = parseCommandLineAndCheck(options, leftFilename, rightFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return executeParwiseAlignmentCommand(globalOptions, options, leftFilename, rightFilename, ClassicDP());
    } else {
        SEQAN_ASSERT_FAIL("Invalid command!");
    }
    
    return 0;
}
