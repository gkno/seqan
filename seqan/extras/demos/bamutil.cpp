// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/bam_io.h>
#include <seqan/ext_lh3.h>

#if SEQAN_HAS_ZLIB

using namespace seqan;

enum Format
{
    FORMAT_SAM,
    FORMAT_BAM
};

struct Options
{
    bool showHelp;
    bool showVersion;
    CharString inFile;
    CharString outFile;
    Format inFormat;
    Format outFormat;
    unsigned verbosity;

    Options()
    {
        showHelp = false;
        showVersion = false;
        inFile = "-";
        outFile = "-";
        inFormat = FORMAT_SAM;
        outFormat = FORMAT_SAM;
        verbosity = 2;
    }
};

void
setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
    addVersionLine(parser, "1.0");
    
    addTitleLine(parser, "***********");
    addTitleLine(parser, "* BAMUTIL *");
    addTitleLine(parser, "***********");
    addTitleLine(parser, "");
    addTitleLine(parser, "SeqAn SAM/BAM Utility.");
    addTitleLine(parser, "");
    addTitleLine(parser, "Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>");

    addUsageLine(parser, "bamutil [OPTIONS] <IN >OUT");
    
	addSection(parser, "General Options");
    addOption(parser, CommandLineOption("v", "verbose", "Enable verbose mode.", OptionType::Bool));
    addOption(parser, CommandLineOption("vv", "very-verbose", "Enable very verbose mode.", OptionType::Bool));

	addSection(parser, "Input Specification");
    addOption(parser, CommandLineOption("i", "input-file", "Path to input, '-' for stdin.", OptionType::String, options.inFile));
    addOption(parser, CommandLineOption("S", "input-sam", "Input file is SAM.", OptionType::Bool, options.inFormat == FORMAT_SAM));
    addOption(parser, CommandLineOption("B", "input-bam", "Input file is BAM.", OptionType::Bool, options.inFormat == FORMAT_BAM));

	addSection(parser, "Output Specification");
    addOption(parser, CommandLineOption("o", "output-file", "Path to output, '-' for stdout.", OptionType::String, options.outFile));
    addOption(parser, CommandLineOption("s", "output-sam", "Output file is SAM.", OptionType::Bool, options.outFormat == FORMAT_SAM));
    addOption(parser, CommandLineOption("b", "output-bam", "Output file is BAM.", OptionType::Bool, options.outFormat == FORMAT_BAM));
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

    if (isSetLong(parser, "verbose"))
        options.verbosity = 1;
    if (isSetLong(parser, "very-verbose"))
        options.verbosity = 2;

    getOptionValueLong(parser, "input-file", options.inFile);
    if (isSetLong(parser, "input-sam"))
        options.inFormat = FORMAT_SAM;
    if (isSetLong(parser, "input-bam"))
        options.inFormat = FORMAT_BAM;

    getOptionValueLong(parser, "output-file", options.outFile);
    if (isSetLong(parser, "output-sam"))
        options.outFormat = FORMAT_SAM;
    if (isSetLong(parser, "output-bam"))
        options.outFormat = FORMAT_BAM;

	return 0;
}

template <typename TInStreamOrRecordReader, typename TOutStream, typename TInTag, typename TOutTag>
int performConversion(TInStreamOrRecordReader & in, TOutStream & out, TInTag const & inTag, TOutTag const & outTag)
{
    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;

    TNameStore nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    BamIOContext<TNameStore> context(nameStore, nameStoreCache);

    BamHeader header;
    if (readRecord(header, context, in, inTag) != 0)
    {
        std::cerr << "Could not read header!" << std::endl;
        return 1;
    }
    if (write2(out, header, context, outTag) != 0)
    {
        std::cerr << "Could not write header!" << std::endl;
        return 1;
    }

    BamAlignmentRecord record;
    while (!atEnd(in))
    {
        if (readRecord(record, context, in, inTag) != 0)
        {
            std::cerr << "Could not read alignment record!" << std::endl;
            return 1;
        }
        if (write2(out, record, context, outTag) != 0)
        {
            std::cerr << "Could not write alignment record!" << std::endl;
            return 1;
        }
    }
    return 0;
}
#endif  // #if SEQAN_HAS_ZLIB

int main(int argc, char const * argv[])
{
#if SEQAN_HAS_ZLIB
    using namespace seqan;

    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    int res = parseCommandLineAndCheck(options, parser, argc, argv);
    if (res != 0)
        return 1;
    if (options.showHelp || options.showVersion)
        return 0;

    std::istream * inS = &std::cin;
    std::ostream * outS = &std::cout;
    int inF = STDIN_FILENO;
    int outF = STDOUT_FILENO;

    bool ioGood = true;
    if (options.inFormat == FORMAT_SAM)
    {
        if (options.inFile != "-")
        {
            inS = new std::fstream(toCString(options.inFile), std::ios::binary | std::ios::in);
            if (!inS->good())
            {
                std::cerr << "Could not open file " << options.inFile << std::endl;
                ioGood = false;
            }
        }
    }
    else
    {
        if (options.inFile != "-")
        {
            inF = open(toCString(options.inFile), O_RDONLY);
            if (inF == -1)
            {
                std::cerr << "Could not open file " << options.inFile << std::endl;
                ioGood = false;
            }
        }
    }
    if (options.outFormat == FORMAT_SAM)
    {
        if (options.outFile != "-")
        {
            outS = new std::fstream(toCString(options.outFile), std::ios::binary | std::ios::out);
            if (!outS->good())
            {
                std::cerr << "Could not open file " << options.inFile << std::endl;
                ioGood = false;
            }
        }
    }
    else
    {
        if (options.outFile != "-")
        {
            outF = open(toCString(options.outFile), O_CREAT | O_WRONLY | O_DIRECT, 00666);
            if (outF == -1)
            {
                std::cerr << "Could not open file " << options.outFile << std::endl;
                ioGood = false;
            }
        }
    }

    if (!ioGood)
        goto main_end;
    
    if (options.inFormat == FORMAT_SAM)
    {
        RecordReader<std::istream, SinglePass<> > reader(*inS);
        if (options.outFormat == FORMAT_SAM)
        {
            res = performConversion(reader, *outS, Sam(), Sam());
        }
        else
        {
            BGZF * outBgzf = bgzf_fdopen(outF, "w");
            Stream<Bgzf> bamOutStream(outBgzf);
            res = performConversion(reader, bamOutStream, Sam(), Bam());
            bgzf_close(outBgzf);
        }
    }
    else
    {
        if (options.outFormat == FORMAT_SAM)
        {
            BGZF * inBgzf = bgzf_fdopen(inF, "r");
            if (inBgzf == 0)
            {
                std::cerr << "Could not open BGZF file!" << std::endl;
                res = 1;
            }
            else
            {
                Stream<Bgzf> bamInStream(inBgzf);
                res = performConversion(bamInStream, *outS, Bam(), Sam());
                bgzf_close(inBgzf);
            }
        }
        else
        {
            BGZF * inBgzf = bgzf_fdopen(inF, "r");
            Stream<Bgzf> bamInStream(inBgzf);
            BGZF * outBgzf = bgzf_fdopen(outF, "w");
            Stream<Bgzf> bamOutStream(outBgzf);
            if (inBgzf == 0)
            {
                res = 1;
                std::cerr << "Could not open BGZF file!" << std::endl;
            }
            else if (outBgzf == 0)
            {
                res = 1;
                std::cerr << "Could not open BGZF file!" << std::endl;
            }
            else
            {
                res = performConversion(bamInStream, bamOutStream, Bam(), Bam());
            }
            bgzf_close(outBgzf);
            bgzf_close(inBgzf);
        }
    }
    if (res != 0)
    {
        std::cerr << "Error during conversion!" << std::endl;
        return 1;
    }

main_end:

    if (options.outFile != "-" && options.outFormat == FORMAT_SAM)
        delete outS;
    if (options.outFile != "-" && options.outFormat != FORMAT_SAM)
        close(outF);
    if (options.inFile != "-" && options.inFormat == FORMAT_SAM)
        delete inS;
    if (options.inFile != "-" && options.inFormat != FORMAT_SAM)
        close(inF);
    
    return ioGood ? 0 : 1;
#else  // #if SEQAN_HAS_ZLIB
    return 0;
#endif  // #if SEQAN_HAS_ZLIB
}
