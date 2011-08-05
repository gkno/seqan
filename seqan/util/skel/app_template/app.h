// ==========================================================================
// %(TITLE)s
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
    int i;
    String<CharString> texts;
    
    Options()
    {
        // Set defaults.
        showHelp = false;
        showVersion = false;
        i = 0;
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
	addOption(parser, CommandLineOption("i",  "integer",  "set an integer option", OptionType::Integer | OptionType::Label, options.i));
    
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
    if (isSetLong(parser, "help")) {
        options.showHelp = true;
        return 0;
    }
    if (isSetLong(parser, "version")) {
        options.showVersion = true;
        return 0;
    }
    
    options.texts = getArgumentValues(parser);

	return 0;
}

int mainWithOptions(Options & options)
{
    typedef Iterator<String<CharString> >::Type TIterator;
    std::cout << "Non-option Arguments:" << std::endl;
    for (TIterator it = begin(options.texts); it != end(options.texts); ++it) {
        std::cout << "  " << *it << std::endl;
    }
    
    return 0;
}

#endif  // #ifndef %(HEADER_GUARD)s
