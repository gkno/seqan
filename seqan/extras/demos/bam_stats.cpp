// ==========================================================================
//                                 bam_stats
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
// Read a BAM file and collect statistics on the alignments therein.
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.

#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <seqan/misc/misc_cmdparser.h>

#if SEQAN_HAS_ZLIB

using namespace seqan;

enum Format
{
    FORMAT_AUTO,
    FORMAT_SAM,
    FORMAT_BAM
};

struct Options
{
    bool showHelp;
    bool showVersion;
    unsigned verbosity;
    CharString refFile;
    CharString inFile;
    Format inFormat;

    Options()
    {
        showHelp = false;
        showVersion = false;
        verbosity = 1;
        inFormat = FORMAT_AUTO;
    }
};

void
setupCommandLineParser(CommandLineParser & parser, Options const & options)
{
    addVersionLine(parser, "1.0");
    
    addTitleLine(parser, "*************");
    addTitleLine(parser, "* bam_stats *");
    addTitleLine(parser, "*************");
    addTitleLine(parser, "");
    addTitleLine(parser, "BAM Statistics.");
    addTitleLine(parser, "");
    addTitleLine(parser, "Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>");

    addUsageLine(parser, "bam_stats [OPTIONS] REF.fasta ALIGN.bam");
    
	addSection(parser, "General Options");
    addOption(parser, CommandLineOption("v", "verbose", "Enable verbose mode (show steps).", OptionType::Bool));
    addOption(parser, CommandLineOption("vv", "very-verbose", "Enable very verbose mode (show SAM lines and actual aligments).", OptionType::Bool));

	addSection(parser, "Input Specification");
    addOption(parser, CommandLineOption("i", "input-file", "Path to input, '-' for stdin.", OptionType::String, options.inFile));
    addOption(parser, CommandLineOption("S", "input-sam", "Input file is SAM (default: auto).", OptionType::Bool, options.inFormat == FORMAT_SAM));
    addOption(parser, CommandLineOption("B", "input-bam", "Input file is BAM (default: auto).", OptionType::Bool, options.inFormat == FORMAT_BAM));

    requiredArguments(parser, 2);
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
        options.verbosity = 2;
    if (isSetLong(parser, "very-verbose"))
        options.verbosity = 3;

    getOptionValueLong(parser, "input-file", options.inFile);
    if (isSetLong(parser, "input-sam"))
        options.inFormat = FORMAT_SAM;
    if (isSetLong(parser, "input-bam"))
        options.inFormat = FORMAT_BAM;

    options.refFile = getArgumentValue(parser, 0);
    options.inFile = getArgumentValue(parser, 1);

	return 0;
}

struct Stats
{
    __uint64 numRecords;
    __uint64 alignedRecords;
    String<unsigned> editDistanceHisto;
    String<unsigned> mismatchHisto;
    String<unsigned> insertHisto;
    String<unsigned> deletionHisto;
    String<unsigned> avrgQuality;

    Stats() : numRecords(0), alignedRecords(0)
    {}
};


template <typename TStreamOrReader, typename TSeqString, typename TSpec, typename TFormat>
int doWork(TStreamOrReader & reader, StringSet<TSeqString, TSpec> & seqs, Options const & options, TFormat const & tag)
{
    StringSet<CharString> refNames;
    NameStoreCache<StringSet<CharString> > refNamesCache(refNames);
    BamIOContext<StringSet<CharString> > context(refNames, refNamesCache);
	String<__uint64> qualSum;

    // Read header.
    BamHeader header;
    if (options.verbosity >= 2)
        std::cerr << "Reading header" << std::endl;
    if (readRecord(header, context, reader, tag) != 0)
    {
        std::cerr << "Could not read header!" << std::endl;
        return 1;
    }

    Stats stats;

    // Read alignments.
    BamAlignmentRecord record;
    if (options.verbosity >= 2)
        std::cerr << "Reading alignments" << std::endl;
    Align<Dna5String> align;
    __int64 reads = 0;
    while (!atEnd(reader))
    {
        // Read alignment record.
        if (readRecord(record, context, reader, tag) != 0)
        {
            std::cerr << "Could not read alignment!" << std::endl;
            return 1;
        }
        if (options.verbosity >= 3)
            write2(std::cerr, record, context, Sam());

        // Now, do something with it ;)
        stats.numRecords += 1;  // One more record.
        stats.alignedRecords += !hasFlagUnmapped(record);
        // Compute alignment.
        if (record.rId != -1 && record.rId < static_cast<int>(length(seqs)))
        {
            bamRecordToAlignment(align, seqs[record.rId], record);
            if (options.verbosity >= 3)
                std::cerr << align << std::endl;
            typedef Align<Dna5String>             TAlign;
            typedef typename Row<TAlign>::Type    TRow;
            typedef typename Iterator<TRow>::Type TRowIter;
            unsigned editDistance = 0;
            unsigned posReadFwd = 0;
            for (TRowIter it0 = begin(row(align, 0)), it1 = begin(row(align, 1)); !atEnd(it0); goNext(it0), goNext(it1))
            {
                unsigned posRead = posReadFwd;
                // is read aligned to reverse strand?
                if ((record.flag & 0x10) != 0)
                    posRead = (length(record.seq) - 1) - posReadFwd;
                
                if (isGap(it0) && isGap(it1))
                    continue;
                if (isGap(it0) || isGap(it1))
                {
                    if (isGap(it0))
                    {
                        unsigned len = length(stats.insertHisto);
                        resize(stats.insertHisto, std::max(len, posRead + 1), 0);
                        stats.insertHisto[posRead] += 1;
                        posReadFwd += 1;
                    }
                    else
                    {
                        unsigned len = length(stats.deletionHisto);
                        resize(stats.deletionHisto, std::max(len, posRead + 1), 0);
                        stats.deletionHisto[posRead] += 1;
                    }
                    editDistance += 1;
                    continue;
                }
                if (*it0 != *it1)
                {
                    unsigned len = length(stats.mismatchHisto);
                    resize(stats.mismatchHisto, std::max(len, posRead + 1), 0);
                    stats.mismatchHisto[posRead] += 1;
                    editDistance += 1;
                }
                posReadFwd += 1;
            }
            
            resize(qualSum, std::max(length(qualSum), length(record.qual)), 0);
            if ((record.flag & 0x10) == 0)
            {
                // read aligns to forward strand
                for (unsigned i = 0; i < length(record.qual); ++i)
                    qualSum[i] += record.qual[i] - '!';
            } else
            {
                // read aligns to reverse strand
                for (unsigned i = 0; i < length(record.qual); ++i)
                    qualSum[(length(record.qual) - 1) - i] += record.qual[i] - '!';
            }
            ++reads;
            
            if (options.verbosity >= 3)
                std::cerr << "edit distance: " << editDistance << std::endl;
            unsigned len = length(stats.editDistanceHisto);
            resize(stats.editDistanceHisto, std::max(len, editDistance + 1), 0);
            stats.editDistanceHisto[editDistance] += 1;
        }
    }
    
    // compute error probabilities
    resize(stats.avrgQuality, length(qualSum));
    for (unsigned i = 0; i < length(qualSum); ++i)
        stats.avrgQuality[i] = (double)qualSum[i] / (double)reads;

    // Print results.
    std::cout << "RESULTS\n\n";
    std::cout << "num records     \t" << stats.numRecords << std::endl;
    std::cout << "aligned records \t" << stats.alignedRecords << std::endl;
    std::cout << "aligned record %\t" << 100.0 * stats.alignedRecords / stats.numRecords << std::endl;
    std::cout << "distance histogram" << std::endl;
    std::cout << "#distance\tnumber of reads" << std::endl;
    for (unsigned i = 0; i < length(stats.editDistanceHisto); ++i)
        std::cout << i << '\t' << stats.editDistanceHisto[i] << std::endl;
    std::cout << "read position histogram" << std::endl;
    std::cout << "#position\tmismatches\tinsertions\tdeletions\terror prob\tquality-based error prob\tquality" << std::endl;
    for (unsigned i = 0; i < length(stats.mismatchHisto); ++i)
    {
        std::cout << i << '\t';
        std::cout << stats.mismatchHisto[i] << '\t';
        std::cout << stats.insertHisto[i] << '\t';
        std::cout << stats.deletionHisto[i] << '\t';
        std::cout << (stats.mismatchHisto[i] + stats.insertHisto[i] + stats.deletionHisto[i]) / (double)reads << '\t';
		double e = stats.avrgQuality[i] * log(10.0) / -10.0;
        std::cout << exp(e) << '\t';
        std::cout << stats.avrgQuality[i] << std::endl;
    }
    return 0;
}

int main(int argc, char const ** argv)
{
    using namespace seqan;

    // -----------------------------------------------------------------------
    // Handle Command Line
    // -----------------------------------------------------------------------

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

    // -----------------------------------------------------------------------
    // Guess format.
    // -----------------------------------------------------------------------

    if (options.inFormat == FORMAT_AUTO)
    {
        std::ifstream guessIn(toCString(options.inFile), std::ios::binary | std::ios::in);
        if (!guessIn.good())
        {
            std::cerr << "Could not open " << options.inFile << std::endl;
            return 1;
        }
        CharString magic;
        resize(magic, 2);
        guessIn.read(&magic[0], 2);
        if (magic != "\x1f\x8b")  // Is not gzip-compressed.
            options.inFormat = FORMAT_SAM;
        else
            options.inFormat = FORMAT_BAM;
    }

    // -----------------------------------------------------------------------
    // Do Work.
    // -----------------------------------------------------------------------

    // Read reference.
    if (options.verbosity >= 2)
        std::cerr << "Reading references from " << options.refFile << std::endl;
    StringSet<CharString> seqIds;
    StringSet<Dna5String> seqs;
    String<char, MMap<> > seqMMapString;
    if (!open(seqMMapString, toCString(options.refFile), OPEN_RDONLY))
    {
        std::cerr << "Could not open " << options.refFile << std::endl;
        return 1;
    }
    RecordReader<String<char, MMap<> >, DoublePass<Mapped> > refReader(seqMMapString);
    if (read2(seqIds, seqs, refReader, Fasta()) != 0)
    {
        std::cerr << "Could not read reference from " << options.refFile << std::endl;
        return 1;
    }

    // Open SAM/BAM file and do work.
    if (options.inFormat == FORMAT_SAM)
    {
        if (options.verbosity >= 2)
            std::cerr << "Opening SAM file " << options.inFile << std::endl;
        String<char, MMap<> > samMMapString;
        if (!open(samMMapString, toCString(options.inFile), OPEN_RDONLY))
        {
            std::cerr << "Could not open " << options.inFile << std::endl;
            return 1;
        }
        RecordReader<String<char, MMap<> >, SinglePass<Mapped> > samReader(samMMapString);
        return doWork(samReader, seqs, options, Sam());
    }
    else  // options.inFormat == FORMAT_BAM
    {
        if (options.verbosity >= 2)
            std::cerr << "Opening BAM file " << options.inFile << std::endl;
        Stream<Bgzf> bamStream;
        if (!open(bamStream, toCString(options.inFile), "r"))
        {
            std::cerr << "Could not open " << options.inFile << std::endl;
            return 1;
        }
        return doWork(bamStream, seqs, options, Bam());
    }

    return 0;
}

#else

int main(int, char const **)
{
    std::cerr << "bam_stats can only be compiled correctly with zlib." << std::endl;
    return 0;
}

#endif  // #if SEQAN_HAS_ZLIB
