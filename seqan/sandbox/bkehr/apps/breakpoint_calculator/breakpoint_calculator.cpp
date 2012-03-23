// ==========================================================================
//                           breakpoint_calculator
// ==========================================================================
// Copyright (C) 2012 by Birte Kehr
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
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#include <lemon/smart_graph.h>
#include <lemon/matching.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>

#include "breakpoint_calculator.h"

using namespace seqan;

// Program entry point
int main(int argc, char const ** argv)
{
    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    int ret = parseCommandLineAndCheck(options, parser, argc, argv);
    if (ret != 0)
        return ret;
    if (options.showHelp || options.showVersion)
        return 0;
    
    // Finally, launch the program.
    ret = mainWithOptions(options);
    return ret;
}
