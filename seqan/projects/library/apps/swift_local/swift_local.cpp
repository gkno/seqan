#define SEQAN_PROFILE

#include <iostream>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>
#include "swift_local.h"

using namespace seqan;


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
		while (dbPos != dbEndPos && queryPos != queryEndPos && !isGap(row(align, 0), dbPos) && !isGap(row(align, 1), queryPos)) {
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
    
    for (typename Position<TId>::Type i = 0; i < length(databaseID) && value(databaseID, i) > 32; ++i) {
        file << value(databaseID, i);
    }

    file << "\tSwiftLocal";
    file << "\teps-matches";

    if (databaseStrand) {
        file << "\t" << toSourcePosition(row(match, 0), beginPosition(row(match, 0))) + beginPosition(source(row(match, 0))) + 1;
        file << "\t" << toSourcePosition(row(match, 0), endPosition(row(match, 0))) + beginPosition(source(row(match, 0)));
    } else {
        file << "\t" << length(host(source(row(match, 0)))) - 
            (toSourcePosition(row(match, 0), endPosition(row(match, 0))) + beginPosition(source(row(match, 0)))) + 1;
        file << "\t" << length(host(source(row(match, 0)))) - 
            (toSourcePosition(row(match, 0), beginPosition(row(match, 0))) + beginPosition(source(row(match, 0))));
    }

    file << "\t" << _calculateIdentity(match);

    file << "\t" << (databaseStrand ? '+' : '-');

    file << "\t.\t";
    for (typename Position<TId>::Type i = 0; i < length(patternID) && value(patternID, i) > 32; ++i) {
        file << value(patternID, i);
    }

    file << ";seq2Range=" << toSourcePosition(row(match, 1), beginPosition(row(match, 1))) + beginPosition(source(row(match, 1))) + 1;
    file << "," << toSourcePosition(row(match, 1), endPosition(row(match, 1))) + beginPosition(source(row(match, 1)));

    std::stringstream cigar, mutations;
    _getCigarLine(match, cigar, mutations);
    file << ";cigar=" << cigar.str();
    file << ";mutations=" << mutations.str();
    file << "\n";
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

    std::ofstream aliFile;
    aliFile.open("swift_local.align");

    TSize maxLength = 0;
    TSize totalLength = 0;
    TSize numMatches = 0;

    aliFile << "Database sequence: " << databaseID;
    if (!databaseStrand) aliFile << " complement\n";
    else aliFile << "\n";

    for (unsigned i = 0; i < length(matches); i++) {
        if (length(value(matches, i)) == 0) continue;
        std::cout << "Pattern sequence: " << ids[i] << "\n";
        aliFile << "Pattern sequence: " << ids[i] << "\n\n";
        for (TSize j = 0; j < length(value(matches, i)); j++) {
            Align<TInfix> m = value(value(matches, i), j);

            aliFile << "< " << toSourcePosition(row(m, 0), beginPosition(row(m, 0))) + beginPosition(source(row(m, 0)));
            aliFile << " , " << toSourcePosition(row(m, 0), endPosition(row(m, 0))) + beginPosition(source(row(m, 0)));
            aliFile << " >< " << toSourcePosition(row(m, 1), beginPosition(row(m, 1))) + beginPosition(source(row(m, 1)));
            aliFile << " , " << toSourcePosition(row(m, 1), endPosition(row(m, 1))) + beginPosition(source(row(m, 1))) << " >\n";
            aliFile << m;

            TSize len = _max(length(row(m, 0)), length(row(m, 1)));
            totalLength += len;
            if(len > maxLength) maxLength = len;

            _writeGffLine(databaseID, ids[i], databaseStrand, m, file);
        }
        numMatches += length(value(matches, i));
        std::cout << "  # Eps-matches: " << length(value(matches, i)) << std::endl;
    }

    std::cout << "# SWIFT hits: " << numSwiftHits << std::endl;
    std::cout << "# Eps-matches: " << numMatches << std::endl;
    std::cout << "Longest eps-match: " << maxLength << " cols" << std::endl;
    if (numMatches > 1)
        std::cout << "Avg match length: " << totalLength / numMatches << " cols" << std::endl << std::endl;

    aliFile.close();
}

template<typename TStringSet, typename TIdSet>
inline int
_importSequences(CharString const & fileName,
                TStringSet & seqs,
                TIdSet & ids) {
    MultiSeqFile multiSeqFile;
    if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY))
        return 1;

    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);
    
    unsigned seqCount = length(multiSeqFile);
    reserve(seqs, seqCount, Exact());
    reserve(ids, seqCount, Exact());

    String<Dna5> seq;
    String<char> id;
    for(unsigned i = 0; i < seqCount; ++i) {
        assignSeq(seq, multiSeqFile[i], format);
        assignSeqId(id, multiSeqFile[i], format);
        appendValue(seqs, seq, Generous());
        appendValue(ids, id, Generous());
    }
    return seqCount;
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
    addLine(parser, "An implementation of the SWIFT filter algorithm (Rasmussen et al., 2006).");
    addLine(parser, "SWIFT hits are verified using local alignment, gapped X-drop extension");
    addLine(parser, "and extraction of the longest epsilon-match.");

	addSection(parser, "Non-optional Arguments:");
    addOption(parser, CommandLineOption('d', "database", "fasta file containing the database sequence", (OptionType::String | OptionType::Mandatory)));
    addOption(parser, CommandLineOption('q', "query", "file containing the query sequences", (OptionType::String | OptionType::Mandatory)));
    
	addSection(parser, "Main Options:");
    addOption(parser, CommandLineOption('o', "out", "output file", OptionType::String, "swift_local.gff"));
    addOption(parser, CommandLineOption('r', "reverseComplement", "search also in reverse complement of database", OptionType::Boolean, false));
    addOption(parser, CommandLineOption('k', "kmer", "length of the q-grams", OptionType::Int, 10));
    addOption(parser, CommandLineOption('l', "minLength", "minimal length of epsilon-matches", OptionType::Int, 100));
    addOption(parser, CommandLineOption('e', "epsilon", "maximal error rate", OptionType::Double, 0.05));
    addOption(parser, CommandLineOption('x', "x-drop", "maximal x-drop for extension", OptionType::Int, 5));
}

int main(int argc, const char *argv[]) {

//-d "Z:\GenomeData\NC_001405_short.fa" -q "Z:\GenomeData\NC_001460_short.fa" -k 5 -l 30 -e 0.1 -x 10 -r
//-d "Z:\GenomeData\adenoviruses\NC_001405.fa" -q "Z:\GenomeData\adenoviruses\NC_001460.fa" -k 5 -l 30 -e 0.1 -x 5 -r

    // Get rid of warnings for unused variables.
    (void)argc;
    (void)argv;

    // command line parsing
    CommandLineParser parser("swift_local");
    _setParser(parser);
    if (!parse(parser, argc, argv)) {
        shortHelp(parser, std::cerr);
        return 1;
    }
    else if (isSetShort(parser, 'h')) {
        return 0;
    }

    // import database sequence
    CharString databaseFile;
    getOptionValueShort(parser, 'd', databaseFile);
    StringSet<String<Dna5> > databases;
    StringSet<CharString> databaseIDs;
    std::cout << "Loaded " << _importSequences(databaseFile, databases, databaseIDs) << " database sequences." << std::endl;

    // import query sequences
    CharString queryFile;
    getOptionValueShort(parser, 'q', queryFile);
    StringSet<String<Dna5> > queries;
    StringSet<CharString> queryIDs;
    std::cout << "Loaded " << _importSequences(queryFile, queries, queryIDs) << " query sequences." << std::endl;

    // input parameters
    int q = 10;
    if (isSetLong(parser, "kmer")) getOptionValueLong(parser, "kmer", q);
    int minLength = 100;
    if (isSetLong(parser, "minLength")) getOptionValueLong(parser, "minLength", minLength);
    double eps = 0.05;
    if (isSetShort(parser, 'e')) getOptionValueShort(parser, 'e', eps);
    if (eps > 0.25) {
        std::cerr << "Please choose a smaller error rate." << std::endl;
        return 1;
    }
    int xDrop = 10;
    if (isSetShort(parser, 'x')) getOptionValueShort(parser, 'x', xDrop);
    CharString outFile = "swift_local.gff";
    if (isSetShort(parser, 'o')) getOptionValueShort(parser, 'o', outFile);
    std::ofstream file;
    file.open(toCString(outFile));
    bool complement = false;
    if (isSetShort(parser, 'r')) getOptionValueShort(parser, 'r', complement);

    std::cout << std::endl;

    SEQAN_PROTIMESTART(timeLocalSwift);

    // pattern
    typedef Index<StringSet<String<Dna5> >, Index_QGram<SimpleShape> > TQGramIndex;
    TQGramIndex index_qgram(queries);
    resize(indexShape(index_qgram), q);

    typedef Finder<String<Dna5>, Swift<SwiftLocal> > TFinder;
    
    // container for eps-matches
    typedef Infix<GetSequenceByNo<TQGramIndex const >::Type >::Type TInfix;
    StringSet<String<Align<TInfix> > > matches;

    int numSwiftHits;
    for(unsigned i = 0; i < length(databases); ++i) {
        numSwiftHits = 0;
        clear(matches);
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << databaseIDs[i] << std::endl;
	    //file << databaseIDs[i] << "\n";
        //pattern
        Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_qgram);
        // finder
        TFinder finder_swift(databases[i], 1000, 1);
        // local swift
        numSwiftHits += localSwift(finder_swift, pattern_swift, eps, minLength, xDrop, matches);
        // file output
        _outputMatches(matches, numSwiftHits, databaseIDs[i], true, queryIDs, file);
        std::cout << std::endl;
    }

    if (complement) {
        // local swift on reverse complement of database
        std::cout << "================================================" << std::endl;
        std::cout << "Database sequence(s) reverse complemented: " << std::endl;
        //file << "\n\nDatabase sequence(s) reverse complemented: \n";
        reverseComplementInPlace(databases);
        for(unsigned i = 0; i < length(databases); ++i) {
            clear(matches);
            numSwiftHits = 0;
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << databaseIDs[i] << std::endl;
            //file << databaseIDs[i] << "\n";
            // pattern
            Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_qgram);
            // finder
            TFinder finder_swift_compl(databases[i], 1000, 1);
            // local swift
            numSwiftHits += localSwift(finder_swift_compl, pattern_swift, eps, minLength, xDrop, matches);
            // file output
            _outputMatches(matches, numSwiftHits, databaseIDs[i], false, queryIDs, file);
            std::cout << std::endl;
        }
    }

    std::cout << "Running time: " << SEQAN_PROTIMEDIFF(timeLocalSwift) << "s" << std::endl;

    file.close();

    // read alignment file
    //std::fstream fstrm;
    //fstrm.open(/* some file*/ );
    //Align<DnaString> ali;
    //read(fstrm, ali, FastaAlign());
    //fstrm.close();

	return 0;
}
