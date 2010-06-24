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

#include <iostream>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>
#include "swift_local.h"

using namespace seqan;

struct SwiftLocalOptions {
	// i/o options
	CharString databaseFile;	// name of database file
	CharString queryFile;		// name of query file
	CharString outputFile;		// name of result file
	unsigned outputFormat;		// 1..gff
								// 2..??

	// main options
	unsigned qGram;				// length of the q-grams
	double epsilon;				// maximal error rate
	int minLength;				// minimal length of an epsilon-match
	double xDrop;				// maximal x-drop

	// more options
	bool reverse;				// compute also matches to reverse complemented database
	CharString fastOption;		// verification strategy: exact, bestLocal, bandedGlobal


	SwiftLocalOptions() {
		outputFile = "swift_local.gff";
		outputFormat = 1;

		qGram = 10;
		epsilon = 0.05;
		minLength = 100;
		xDrop = 5;

		reverse = false;
		fastOption = "exact";		// exact verification
	}
};

template<typename TAlign, typename TString>
void
_getCigarLine(TAlign const & align, TString & cigar, TString & mutations) { 
    typedef typename Size<typename Row<TAlign>::Type >::Type TSize;

    TSize dbPos = beginPosition(row(align, 0));
    TSize queryPos = beginPosition(row(align, 1));

    TSize dbEndPos = endPosition(row(align, 0));
    TSize queryEndPos = endPosition(row(align, 1));

    bool first = true;
    TSize readBasePos = queryPos;
    TSize readPos = 0;
	while (dbPos != dbEndPos && queryPos != queryEndPos) {
		int matched = 0;
		int inserted = 0;
		int deleted = 0;
		while (dbPos != dbEndPos && queryPos != queryEndPos &&
               !isGap(row(align, 0), dbPos) && !isGap(row(align, 1), queryPos)) {
            ++readPos;
			if (value(row(align, 0), dbPos) != value(row(align, 1), queryPos)) {
				if (first) first = false;
				else mutations << ",";
				mutations << readPos << value(source(row(align, 1)), readBasePos);
			}
			++readBasePos;
			++dbPos;
			++queryPos;
			++matched;
		}
		if (matched > 0) cigar << matched << "M" ;
		while (queryPos != queryEndPos && isGap(row(align, 1), queryPos)) {
			++dbPos;
			++queryPos;
			++deleted;
		}
		if (deleted > 0) cigar << deleted << "D";
		while (dbPos != dbEndPos && isGap(row(align, 0), dbPos)) {
			++dbPos;
			++queryPos;
			++readPos;
			if (first) first = false;
			else mutations << ",";
			mutations << readPos << value(source(row(align, 1)), readBasePos);
			++readBasePos;
			++inserted;
		}
		if (inserted > 0) cigar << inserted << "I";
	}
}

template<typename TAlign>
double
_calculateIdentity(TAlign const & align) {
    typedef typename Size<typename Row<TAlign>::Type >::Type TSize;
    TSize matches = 0;
    TSize len = _max(length(row(align, 0)), length(row(align, 1)));

    TSize pos0 = beginPosition(row(align, 0));
    TSize pos1 = beginPosition(row(align, 1));

    TSize end0 = endPosition(row(align, 0));
    TSize end1 = endPosition(row(align, 1));

    while ((pos0 < end0) && (pos1 < end1)) {
        if (!isGap(row(align, 0), pos0) && !isGap(row(align, 1), pos1)) {
            if (value(row(align, 0), pos0) == value(row(align , 1), pos1)) {
                ++matches;
            }
        }
        ++pos0;
        ++pos1;
    }

    return ((double)matches/(double)len)*100.0;
}

template<typename TId, typename TAlign, typename TFile>
void
_writeGffLine(TId const & databaseID,
              TId const & patternID,
              bool const databaseStrand,
              TAlign const & match,
              TFile & file) {
	typedef typename Row<TAlign>::Type TRow;
	TRow row0 = row(match, 0);
	TRow row1 = row(match, 1);
    
    for (typename Position<TId>::Type i = 0; i < length(databaseID) && value(databaseID, i) > 32; ++i) {
        file << value(databaseID, i);
    }

    file << "\tSwiftLocal";
    file << "\teps-matches";

    if (databaseStrand) {
        file << "\t" << 
			toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0)) + 1;
        file << "\t" << 
			toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
    } else {
        file << "\t" << length(source(row0)) - 
            (toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0))) + 1;
        file << "\t" << length(source(row0)) - 
            (toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0)));
    }

    file << "\t" << _calculateIdentity(match);

    file << "\t" << (databaseStrand ? '+' : '-');

    file << "\t.\t";
    for (typename Position<TId>::Type i = 0; i < length(patternID) && value(patternID, i) > 32; ++i) {
        file << value(patternID, i);
    }

	//file << ";seq2Length=" << length(source(row1));

    file << ";seq2Range=" << 
		toSourcePosition(row1, beginPosition(row1)) + beginPosition(source(row1)) + 1;
    file << "," << 
		toSourcePosition(row1, endPosition(row1)) + beginPosition(source(row1));

    std::stringstream cigar, mutations;
    _getCigarLine(match, cigar, mutations);
    file << ";cigar=" << cigar.str();
    file << ";mutations=" << mutations.str();
    file << "\n";
}

template<typename TAlign, typename TStrand, typename TFile>
void
_writeMatch(TAlign & match,
			TStrand databaseStrand,
			TFile & aliFile) {
	typedef typename Row<TAlign>::Type TRow;
	TRow row0 = row(match, 0);
	TRow row1 = row(match, 1);

	// write database positions
	if (databaseStrand) {
		aliFile << "< " <<
			toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0));
		aliFile << " , " << 
			toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
	} else {
		aliFile << "< " << length(source(row0)) - 
			toSourcePosition(row0, beginPosition(row0)) + beginPosition(source(row0));
		aliFile << " , " << length(source(row0)) -
			toSourcePosition(row0, endPosition(row0)) + beginPosition(source(row0));
	}
	// write query positions
	aliFile << " >< " << 
		toSourcePosition(row1, beginPosition(row1)) + beginPosition(source(row1));
	aliFile << " , " << 
		toSourcePosition(row1, endPosition(row1)) + beginPosition(source(row1)) << " >\n";

	// write match
	aliFile << match;
}

template<typename TInfix, typename TNumber, typename TId, typename TIds, typename TFile>
void
_outputMatches(StringSet<String<Align<TInfix> > > const & matches, 
              TNumber const numSwiftHits,
              TId const & databaseID,
              bool const databaseStrand,
              TIds const & ids,
              TFile & file) {
    typedef typename Size<Align<TInfix> >::Type TSize;

    //std::ofstream aliFile;
    //aliFile.open("swift_local.align");

    TSize maxLength = 0;
    TSize totalLength = 0;
    TSize numMatches = 0;

    //aliFile << "Database sequence: " << databaseID;
    //if (!databaseStrand) aliFile << " complement\n";
    //else aliFile << "\n";

    for (unsigned i = 0; i < length(matches); i++) {
        if (length(value(matches, i)) == 0) continue;
        //std::cout << "Pattern sequence: " << ids[i] << "\n";
        //aliFile << "Pattern sequence: " << ids[i] << "\n\n";
        for (TSize j = 0; j < length(value(matches, i)); j++) {
            Align<TInfix> m = value(value(matches, i), j);

            TSize len = _max(length(row(m, 0)), length(row(m, 1)));
            totalLength += len;
            if(len > maxLength) maxLength = len;

            _writeGffLine(databaseID, ids[i], databaseStrand, m, file);
			//_writeMatch(m, databaseStrand, aliFile);
        }
        numMatches += length(value(matches, i));
        //std::cout << "  # Eps-matches: " << length(value(matches, i)) << std::endl;
    }

	if (numMatches > 0) {
		std::cout << "    Longest eps-match: " << maxLength << std::endl;
        std::cout << "    Avg match length : " << totalLength / numMatches << std::endl;
	}
    std::cout << "    # SWIFT hits     : " << numSwiftHits << std::endl;
    std::cout << "    # Eps-matches    : " << numMatches << std::endl;

    //aliFile.close();
}

template<typename TSequence, typename TId, typename TPattern, typename TFile>
inline bool
_swiftLocalOnOne(TSequence & database,
				TId & databaseID,
				TPattern & pattern_swift,
				StringSet<TId> & queryIDs,
				SwiftLocalOptions & options,
				TFile & file) {
	std::cout << "  " << databaseID << std::endl;
    int numSwiftHits;

	// finder
    typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
	TFinder finder_swift(database, 1000, 1);

    // container for eps-matches
	StringSet<String<Align<TSequence> > > matches;

	// local swift
	if (options.fastOption == CharString("exact"))
		numSwiftHits = localSwift(finder_swift, pattern_swift, options.epsilon,
								   options.minLength, options.xDrop, matches, AllLocal());
	else if (options.fastOption == "bestLocal")
		numSwiftHits = localSwift(finder_swift, pattern_swift, options.epsilon,
								   options.minLength, options.xDrop, matches, BestLocal());
	else if (options.fastOption == "bandedGlobal")
		numSwiftHits = localSwift(finder_swift, pattern_swift, options.epsilon,
								   options.minLength, options.xDrop, matches, BandedGlobal());
	else if (options.fastOption == "bandedGlobalExtend")
		numSwiftHits = localSwift(finder_swift, pattern_swift, options.epsilon,
								   options.minLength, options.xDrop, matches, BandedGlobalExtend());
	else {
		std::cerr << "Unknown verification strategy: " << options.fastOption << std::endl;
		return false;
	}

	// file output
	_outputMatches(matches, numSwiftHits, databaseID, true, queryIDs, file);

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
		if(!_swiftLocalOnOne(databases[i], databaseIDs[i], pattern_swift, queryIDs, options, file))
			return false;
    }
	std::cout << std::endl;

	// negative (reverse complemented) database strand
    if (options.reverse) {
        std::cout << "Aligning all query sequences to reverse complement of database sequence";
		std::cout << ((length(databases)>1)?"s...":"...") << std::endl;

        reverseComplementInPlace(databases);

        for(unsigned i = 0; i < length(databases); ++i) {
			if(!_swiftLocalOnOne(databases[i], databaseIDs[i], pattern_swift, queryIDs, options, file))
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
	std::cout << "  reverse complement database : " << ((options.reverse)?"yes":"no") << std::endl;
	std::cout << "  verification strategy       : " << options.fastOption << std::endl;
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
    addOption(parser, CommandLineOption('e', "epsilon", "maximal error rate", OptionType::Double, 0.05));
    addOption(parser, CommandLineOption('r', "reverseComplement", "search also in reverse complement of database",
              OptionType::Boolean, false));

	addSection(parser, "Verification Options:");
    addOption(parser, CommandLineOption('x', "xDrop", "maximal x-drop for extension", OptionType::Double, 5));
	addOption(parser, CommandLineOption('v', "verification", "verification strategy", OptionType::String, "exact"));
	addHelpLine(parser, "exact        = compute and extend all local alignments in SWIFT hits");
	addHelpLine(parser, "bestLocal    = compute and extend only best local alignment in SWIFT hits");
	addHelpLine(parser, "bandedGlobal = banded global alignment on SWIFT hits");

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
