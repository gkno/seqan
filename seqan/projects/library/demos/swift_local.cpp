#define SEQAN_PROFILE

#include <iostream>
#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>
#include <demos/swift_local.h>

using namespace seqan;

template<typename TInfix, typename TNumber, typename TFile>
void
outputMatches(StringSet<String<Align<TInfix> > > matches, 
              TNumber numSwiftHits, 
              TFile outFile) {
    typedef typename Size<Align<TInfix> >::Type TSize;

    TSize maxLength = 0;
    TSize totalLength = 0;
    TSize numMatches = 0;

    std::ofstream file;
    file.open(outFile);
    for (unsigned i = 0; i < length(matches); i++) {
        std::cout << "Pattern sequence " << i << ":\n";
        file << "Pattern sequence " << i << ":\n\n";
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
    file.close();

    std::cout << "# SWIFT hits: " << numSwiftHits << std::endl;
    std::cout << "# Eps-matches: " << numMatches << std::endl;
    std::cout << "Longest eps-match: " << maxLength << " cols" << std::endl;
    if (numMatches > 1)
        std::cout << "Avg match length: " << totalLength / numMatches << " cols" << std::endl << std::endl;
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
    addOption(parser, CommandLineOption('q', "query", "fasta file containing the query sequences", (OptionType::String | OptionType::Mandatory)));
    
	addSection(parser, "Main Options:");
    addOption(parser, CommandLineOption('o', "out", "output file", OptionType::String));
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

    // input parameters

    CharString databaseFile, queryFile;
    getOptionValueShort(parser, 'd', databaseFile);
    getOptionValueShort(parser, 'q', queryFile);
    String<Dna5> database = String<Dna5, FileReader<Fasta> >(databaseFile);
    StringSet<String<Dna5> > queries;
    appendValue(queries, String<Dna5, FileReader<Fasta> >(queryFile));

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


    SEQAN_PROTIMESTART(timeLocalSwift);

    // finder
    typedef Finder<String<Dna5>, Swift<SwiftLocal> > TFinder;
  	TFinder finder_swift(database, 1000, 1);

    // pattern
    typedef Index<StringSet<String<Dna5> >, Index_QGram<SimpleShape> > TQGramIndex;
    TQGramIndex index_qgram(queries);
    resize(indexShape(index_qgram), q);
    Pattern<TQGramIndex, Swift<SwiftLocal> > pattern_swift(index_qgram);

    // container for eps-matches
    typedef Infix<GetSequenceByNo<TQGramIndex const >::Type >::Type TInfix;
    StringSet<String<Align<TInfix> > > matches;

    // local swift
    int numSwiftHits = localSwift(finder_swift, pattern_swift, eps, minLength, xDrop, matches);

    std::cout << "Running time: " << SEQAN_PROTIMEDIFF(timeLocalSwift) << "s" << std::endl;

    // file output
    outputMatches(matches, numSwiftHits, toCString(outFile));

    // read alignment file
    //std::fstream fstrm;
    //fstrm.open(/* some file*/ );
    //Align<DnaString> ali;
    //read(fstrm, ali, FastaAlign());
    //fstrm.close();

	return 0;
}

