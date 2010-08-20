 /*==========================================================================
                     SwiftLocal - Fast Local Alignment

 ============================================================================
  Copyright (C) 2010 by Birte Kehr

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#define SEQAN_PROFILE

#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>
//#include "swift_local_types.h"
#include "swift_local.h"
#include "swift_local_output.h"

using namespace seqan;

template<typename TSequence, typename TId, typename TPattern, typename TFile>
inline bool
_swiftLocalOnOne(TSequence & database,
				TId & databaseID,
				TPattern & pattern_swift,
				StringSet<TId> & queryIDs,
				bool databaseStrand,
				SwiftLocalOptions & options,
				TFile & file) {
	std::cout << "  " << databaseID << std::endl;
    int numSwiftHits;

	// finder
    typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
	TFinder finder_swift(database, 1000, 1);

    // container for eps-matches
	//StringSet<String<Align<TSequence> > > matches;
	StringSet<String<SwiftLocalMatch<TSequence, TId> > > matches;

	// local swift
	if (options.fastOption == CharString("exact"))
		numSwiftHits = localSwift(finder_swift, pattern_swift, options.epsilon, options.minLength, options.xDrop, 
								  options.compactThresh, options.numMatches, queryIDs, matches, AllLocal());
	else if (options.fastOption == "bestLocal")
		numSwiftHits = localSwift(finder_swift, pattern_swift, options.epsilon, options.minLength, options.xDrop, 
								  options.compactThresh, options.numMatches, queryIDs, matches, BestLocal());
	else if (options.fastOption == "bandedGlobal")
		numSwiftHits = localSwift(finder_swift, pattern_swift, options.epsilon, options.minLength, options.xDrop, 
								  options.compactThresh, options.numMatches, queryIDs, matches, BandedGlobal());
	else if (options.fastOption == "bandedGlobalExtend")
		numSwiftHits = localSwift(finder_swift, pattern_swift, options.epsilon, options.minLength, options.xDrop, 
								  options.compactThresh, options.numMatches, queryIDs, matches, BandedGlobalExtend());
	else {
		std::cerr << "Unknown verification strategy: " << options.fastOption << std::endl;
		return false;
	}

	// file output
	_outputMatches(matches, numSwiftHits, databaseID, databaseStrand, queryIDs, file);

	return true;
}

template<typename TSequence, typename TId, typename TFile>
inline bool
_swiftLocalOnAll(StringSet<TSequence> & databases,
				StringSet<TId> & databaseIDs,
				StringSet<TSequence> & queries,
				StringSet<TId> & queryIDs,
				SwiftLocalOptions & options,
				TFile & file) {
    // pattern
    typedef Index<StringSet<TSequence, Dependent<> >, Index_QGram<SimpleShape, OpenAddressing> > TQGramIndex;
    TQGramIndex index_qgram(queries);
    resize(indexShape(index_qgram), options.qGram);
	Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_qgram);

	// positive database strand
	std::cout << "Aligning all query sequences to database sequence";
	std::cout << ((length(databases)>1)?"s...":"...") << std::endl;

    for(unsigned i = 0; i < length(databases); ++i) {
		if(!_swiftLocalOnOne(databases[i], databaseIDs[i], pattern_swift, queryIDs, true, options, file))
			return false;
    }
	std::cout << std::endl;

	// negative (reverse complemented) database strand
    if (options.reverse) {
        std::cout << "Aligning all query sequences to reverse complement of database sequence";
		std::cout << ((length(databases)>1)?"s...":"...") << std::endl;

        reverseComplementInPlace(databases);

        for(unsigned i = 0; i < length(databases); ++i) {
			if(!_swiftLocalOnOne(databases[i], databaseIDs[i], pattern_swift, queryIDs, false, options, file))
				return false;
        }
		std::cout << std::endl;
    }

    return true;
}

template<typename TSequence, typename TId>
inline bool
_importSequences(CharString const & fileName,
				 CharString const & name,
                 StringSet<TSequence> & seqs,
                 StringSet<TId> & ids) {
    MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY)) {
		std::cerr << "Failed to open " << name << " file." << std::endl;
        return false;
	}

    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);
    
    unsigned seqCount = length(multiSeqFile);
    reserve(seqs, seqCount, Exact());
    reserve(ids, seqCount, Exact());

    TSequence seq;
    TId id;
    for(unsigned i = 0; i < seqCount; ++i) {
        assignSeq(seq, multiSeqFile[i], format);
        assignSeqId(id, multiSeqFile[i], format);
        appendValue(seqs, seq, Generous());
        appendValue(ids, id, Generous());
    }

	std::cout << "Loaded " << seqCount << " " << name << " sequence" << ((seqCount>1)?"s.":".") << std::endl;
    return true;
}

void _writeCalculatedParams(SwiftLocalOptions & options) {
	// Calculate parameters
	int n = (int) ceil((floor(options.epsilon * options.minLength) + 1) / options.epsilon);
	int threshold = (int) _max(1, (int) _min(
		(n + 1) - options.qGram * (floor(options.epsilon * n) + 1),
		(options.minLength + 1) - options.qGram * (floor(options.epsilon * options.minLength) + 1)));
	int overlap = (int) floor((2 * threshold + options.qGram - 3) / (1 / options.epsilon - options.qGram));
	int distanceCut = (threshold - 1) + options.qGram * overlap + options.qGram;
	int logDelta = _max(4, (int) ceil(log((double)overlap + 1) / log(2.0)));
	int delta = 1 << logDelta;

	// Output calculated parameters
	std::cout << "Calculated parameters:" << std::endl;
	std::cout << "  threshold   : " << threshold << std::endl;
	std::cout << "  distance cut: " << distanceCut << std::endl;
	std::cout << "  delta       : " << delta << std::endl;
	std::cout << "  overlap     : " << overlap << std::endl;
	std::cout << std::endl;
}

template<typename TOptions>
void
_writeSpecifiedParams(TOptions & options) {
	// Output user specified parameters
	std::cout << "User specified parameters:" << std::endl;
	std::cout << "  minimal match length        : " << options.minLength << std::endl;
	std::cout << "  maximal error rate (epsilon): " << options.epsilon << std::endl;
	std::cout << "  maximal x-drop              : " << options.xDrop << std::endl;
	std::cout << "  q-gram length               : " << options.qGram << std::endl;
	std::cout << "  reverse complement database : " << ((options.reverse)?"yes":"no") << std::endl << std::endl;

	std::cout << "  verification strategy       : " << options.fastOption << std::endl;
	std::cout << "  maximal number of matches   : " << options.numMatches << std::endl;
	std::cout << "  duplicate removal every     : " << options.compactThresh << std::endl;
	std::cout << std::endl;
}

template<typename TOptions>
void
_writeFileNames(TOptions & options) {
	std::cout << "Database file: " << options.databaseFile << std::endl;
	std::cout << "Query file   : " << options.queryFile << std::endl;
	std::cout << "Output file  : " << options.outputFile << std::endl;
	std::cout << std::endl;
}

template<typename TParser, typename TOptions>
bool
_parseOptions(TParser & parser, TOptions & options) {
    // i/o options
	getOptionValueShort(parser, 'd', options.databaseFile);
    getOptionValueShort(parser, 'q', options.queryFile);
    if (isSetShort(parser, 'o')) getOptionValueShort(parser, 'o', options.outputFile);

	// main options
	if (isSetLong(parser, "kmer")) getOptionValueLong(parser, "kmer", options.qGram);
    if (isSetLong(parser, "minLength")) getOptionValueLong(parser, "minLength", options.minLength);
	if (isSetShort(parser, 'e')) getOptionValueShort(parser, 'e', options.epsilon);
    if (options.epsilon > 0.25) {
        std::cerr << "Please choose a smaller error rate." << std::endl;
        return 0;
    }
    if (isSetShort(parser, 'x')) getOptionValueShort(parser, 'x', options.xDrop);

	if (isSetShort(parser, 'r')) getOptionValueShort(parser, 'r', options.reverse);
	if (isSetShort(parser, 'v')) getOptionValueShort(parser, 'v', options.fastOption);
	if (isSetShort(parser, 'n')) getOptionValueShort(parser, 'n', options.numMatches);
	if (isSetShort(parser, 's')) getOptionValueShort(parser, 's', options.compactThresh);
	return 1;
}

template<typename TParser>
void
_setParser(TParser & parser) {
    addTitleLine(parser, "******************************************");
	addTitleLine(parser, "* Local alignment using the SWIFT filter *");
	addTitleLine(parser, "* (c) Copyright 2010 by Birte Kehr       *");
	addTitleLine(parser, "******************************************");

	addUsageLine(parser, "-d <FASTA sequence file> -q <FASTA sequence file> [Options]");

	addLine(parser, "");
    addLine(parser, "An implementation of the SWIFT filter algorithm (Rasmussen et al., 2006)");
    addLine(parser, "and subsequent verification of the SWIFT hits using local alignment,");
    addLine(parser, "gapped X-drop extension, and extraction of the longest epsilon-match.");

	addSection(parser, "Non-optional Arguments:");
    addOption(parser, CommandLineOption('d', "database", "fasta file containing the database sequences",
              (OptionType::String | OptionType::Mandatory)));
    addOption(parser, CommandLineOption('q', "query", "fasta file containing the query sequences", 
              (OptionType::String | OptionType::Mandatory)));
    
	addSection(parser, "Filtering Options:");
    addOption(parser, CommandLineOption('k', "kmer", "length of the q-grams", OptionType::Int, 10));
    addOption(parser, CommandLineOption('l', "minLength", "minimal length of epsilon-matches", OptionType::Int, 100));
    addOption(parser, CommandLineOption('e', "epsilon", "maximal error rate", OptionType::Double, "0.05"));
    addOption(parser, CommandLineOption('r', "reverseComplement", "search also in reverse complement of database",
              OptionType::Boolean, false));

	addSection(parser, "Verification Options:");
    addOption(parser, CommandLineOption('x', "xDrop", "maximal x-drop for extension", OptionType::Double, 5));
	addOption(parser, CommandLineOption('v', "verification", "verification strategy", OptionType::String, "exact"));
	addHelpLine(parser, "exact        = compute and extend all local alignments in SWIFT hits");
	addHelpLine(parser, "bestLocal    = compute and extend only best local alignment in SWIFT hits");
	addHelpLine(parser, "bandedGlobal = banded global alignment on SWIFT hits");
	addOption(parser, CommandLineOption('n', "numMatches", "maximal number of kept matches per query and database", OptionType::Int, 50));
	addOption(parser, CommandLineOption('s', "sort", "number of matches for removal of duplicates", OptionType::Int, 500));

	addSection(parser, "Output Options:");
    addOption(parser, CommandLineOption('o', "out", "name of output file", OptionType::String, "swift_local.gff"));
	addOption(parser, CommandLineOption('f', "outFormat", "output format", OptionType::String, "gff"));
	addHelpLine(parser, "possible formats: gff");
}

int main(int argc, const char *argv[]) {

//-d "Z:\GenomeData\NC_001405_short.fa" -q "Z:\GenomeData\NC_001460_short.fa" -k 5 -l 30 -e 0.1 -x 10 -r
//-d "Z:\GenomeData\adenoviruses\NC_001405.fa" -q "Z:\GenomeData\adenoviruses\NC_001460.fa" -k 5 -l 30 -e 0.1 -x 5 -r

    // command line parsing
    CommandLineParser parser("swift_local");

    _setParser(parser);
    if (!parse(parser, argc, argv)) {
		if (isSetShort(parser, 'h')) return 0; 
        shortHelp(parser, std::cerr);
        return 1;
    }

	SwiftLocalOptions options = SwiftLocalOptions();
	if (!_parseOptions(parser, options)) {
        shortHelp(parser, std::cerr);
		return 1;
	}

	typedef String<Dna5> TSequence;

	// output header
	_title(parser, std::cout);
	std::cout << std::endl;

	// output file names
	_writeFileNames(options);

	// output parameters
	_writeSpecifiedParams(options);
	_writeCalculatedParams(options);

    // import query sequences
    StringSet<TSequence > queries;
    StringSet<CharString> queryIDs;
	if (!_importSequences(options.queryFile, "query", queries, queryIDs)) return 1;

    // import database sequence
    StringSet<TSequence > databases;
    StringSet<CharString> databaseIDs;
	if (!_importSequences(options.databaseFile, "database", databases, databaseIDs)) return 1;

    // open output file
    std::ofstream file;
    file.open(toCString(options.outputFile));
	if (!file.is_open()) {
		std::cerr << "Could not open output file." << std::endl;
		return 1;
	}
    std::cout << std::endl;


	// local swift on all databases and queries writing results to file
    SEQAN_PROTIMESTART(timeLocalSwift);
	if(!_swiftLocalOnAll(databases, databaseIDs, queries, queryIDs, options, file))
		return 1;
    std::cout << "Running time: " << SEQAN_PROTIMEDIFF(timeLocalSwift) << "s" << std::endl;

    file.close();
	
	return 0;
}
