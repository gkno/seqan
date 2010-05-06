#define SEQAN_PROFILE

#include <iostream>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>
#include <demos/swift_local.h>

using namespace seqan;

template<typename TInfix, typename TNumber, typename TIds, typename TFile>
void
outputMatches(StringSet<String<Align<TInfix> > > const & matches, 
              TNumber numSwiftHits,
              TIds ids,
              TFile & file) {
    typedef typename Size<Align<TInfix> >::Type TSize;

    TSize maxLength = 0;
    TSize totalLength = 0;
    TSize numMatches = 0;

    for (unsigned i = 0; i < length(matches); i++) {
        std::cout << "Pattern sequence: " << ids[i] << "\n";
        file << "Pattern sequence: " << ids[i] << "\n\n";
        for (TSize j = 0; j < length(value(matches, i)); j++) {
            Align<TInfix> m = value(value(matches, i), j);

            file << "< " << toSourcePosition(row(m, 0), beginPosition(row(m, 0))) + beginPosition(source(row(m, 0)));
            file << " , " << toSourcePosition(row(m, 0), endPosition(row(m, 0))) + beginPosition(source(row(m, 0)));
            file << " >< " << toSourcePosition(row(m, 1), beginPosition(row(m, 1))) + beginPosition(source(row(m, 1)));
            file << " , " << toSourcePosition(row(m, 1), endPosition(row(m, 1))) + beginPosition(source(row(m, 1))) << " >\n";
            file << m;

            totalLength += length(row(m, 0));
            if(length(row(m, 0)) > maxLength) maxLength = length(row(m, 0));
        }
        numMatches += length(value(matches, i));
        std::cout << "  # Eps-matches: " << length(value(matches, i)) << std::endl;
    }

    std::cout << "# SWIFT hits: " << numSwiftHits << std::endl;
    std::cout << "# Eps-matches: " << numMatches << std::endl;
    std::cout << "Longest eps-match: " << maxLength << " cols" << std::endl;
    if (numMatches > 1)
        std::cout << "Avg match length: " << totalLength / numMatches << " cols" << std::endl << std::endl;
}

template<typename TStringSet, typename TIdSet>
inline int
importSequences(CharString const & fileName,
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
setParser(TParser & parser) {
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
    addOption(parser, CommandLineOption('o', "out", "output file", OptionType::String));
    addOption(parser, CommandLineOption('r', "reverseComplement", "search also in reverse complement of database", OptionType::Boolean));
    addOption(parser, CommandLineOption('k', "kmer", "length of the q-grams", OptionType::Int));
    addOption(parser, CommandLineOption('l', "minLength", "minimal length of epsilon-matches", OptionType::Int));
    addOption(parser, CommandLineOption('e', "epsilon", "maximal error rate", OptionType::Double));
    addOption(parser, CommandLineOption('x', "x-drop", "maximal x-drop for extension", OptionType::Int));
}

int main(int argc, const char *argv[]) {
    // Get rid of warnings for unused variables.
    (void)argc;
    (void)argv;

    // command line parsing
    CommandLineParser parser("swift_local");
    setParser(parser);
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
    std::cout << "Loaded " << importSequences(databaseFile, databases, databaseIDs) << " database sequences." << std::endl;

    // import query sequences
    CharString queryFile;
    getOptionValueShort(parser, 'q', queryFile);
    StringSet<String<Dna5> > queries;
    StringSet<CharString> queryIDs;
    std::cout << "Loaded " << importSequences(queryFile, queries, queryIDs) << " query sequences." << std::endl;

    // input parameters
    int q = 10;
    if (isSetLong(parser, "kmer")) getOptionValueLong(parser, "kmer", q);
    int minLength = 100;
    if (isSetLong(parser, "minLength")) getOptionValueLong(parser, "minLength", minLength);
    double eps = 0.05;
    if (isSetShort(parser, 'e')) getOptionValueShort(parser, 'e', eps);
    int xDrop = 10;
    if (isSetShort(parser, 'x')) getOptionValueShort(parser, 'x', xDrop);
    CharString outFile = "swift_local.out";
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
    Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_qgram);

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
	    file << databaseIDs[i] << "\n";
        // finder
        TFinder finder_swift(databases[i], 1000, 1);
        // local swift
        numSwiftHits += localSwift(finder_swift, pattern_swift, eps, minLength, xDrop, matches);
        // file output
        outputMatches(matches, numSwiftHits, queryIDs, file);
        std::cout << std::endl;
    }

    if (complement) {
        // local swift on reverse complement of database
        std::cout << "================================================" << std::endl;
        std::cout << "Database sequence(s) reverse complemented: " << std::endl;
        file << "\n\nDatabase sequence(s) reverse complemented: \n";
        reverseComplementInPlace(databases);
        for(unsigned i = 0; i < length(databases); ++i) {
            clear(matches);
            numSwiftHits = 0;
            std::cout << databaseIDs[i] << std::endl;
            file << databaseIDs[i] << "\n";
            //finder
            TFinder finder_swift_compl(databases[0], 1000, 1);
            // local swift
            numSwiftHits += localSwift(finder_swift_compl, pattern_swift, eps, minLength, xDrop, matches);
            // file output
            outputMatches(matches, numSwiftHits, queryIDs, file);
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
