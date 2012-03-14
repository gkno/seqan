// ==========================================================================
// %(TITLE)s
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: %(AUTHOR)s
// ==========================================================================

#ifndef %(HEADER_GUARD)s
#define %(HEADER_GUARD)s

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Options
{
    bool showHelp;
    bool showVersion;
    CharString inputFileName;
    CharString outputFileName;
    String<CharString> texts;
    
    Options()
    {
        // Set defaults.
        showHelp = false;
        showVersion = false;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

void
setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
    addVersionLine(parser, "0.1");
    
    addTitleLine(parser, "**********************");
    addTitleLine(parser, "* %(NAME)s *");
    addTitleLine(parser, "**********************");
    addTitleLine(parser, "");
    addTitleLine(parser, "(c) %(YEAR)d by %(AUTHOR)s");

    addUsageLine(parser, "[OPTIONS] TEXT+");
    
	addSection(parser, "Main Options");
	addOption(parser, CommandLineOption("if",  "input-file",  "Input file.", OptionType::String | OptionType::Label, options.inputFileName));
	addOption(parser, CommandLineOption("of",  "output-file", "Output file.", OptionType::String | OptionType::Label, options.outputFileName));
    
    requiredArguments(parser, 1);
}

int parseCommandLineAndCheck(Options & options,
                             CommandLineParser & parser,
                             int argc,
                             char const ** argv)
{
    bool stop = !parse(parser, argc, argv);
    if (stop)
        return 1;
    if (isSetLong(parser, "help"))
    {
        options.showHelp = true;
        return 0;
    }
    if (isSetLong(parser, "version"))
    {
        options.showVersion = true;
        return 0;
    }

    getOptionValueLong(parser, "input-file", options.inputFileName);
    getOptionValueLong(parser, "output-file", options.outputFileName);
    
    options.texts = getArgumentValues(parser);

	return 0;
}

int mainWithOptions(Options & options)
{
    typedef Iterator<String<CharString> >::Type TIterator;
    std::cout << "Option Arguments:" << std::endl;
    std::cout << "  input file:  \"" << options.inputFileName << "\"" << std::endl;
    std::cout << "  output file: \"" << options.outputFileName << "\"" << std::endl;
    std::cout << "Non-option Arguments:" << std::endl;
    for (TIterator it = begin(options.texts); it != end(options.texts); ++it)
    {
        std::cout << "  " << *it << std::endl;
    }
    
    return 0;
}

#endif  // #ifndef %(HEADER_GUARD)s
