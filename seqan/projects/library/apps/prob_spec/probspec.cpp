 /*==========================================================================
		ProbSpec - Computing special local matches for probe design

 ============================================================================
  Copyright (C) 2011 by Knut Reinert
  
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
#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>
#include "probspec.h"


using namespace seqan;


///////////////////////////////////////////////////////////////////////////////
// Imports sequences from a file, 
//  stores them in the StringSet seqs and their identifiers in the StringSet ids
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

///////////////////////////////////////////////////////////////////////////////
// Writes user specified parameters from options object to std::cout
template<typename TOptions>
void
_writeSpecifiedParams(TOptions & options) {
//IOREV _notio_
	// Output user specified parameters
	std::cout << "User specified parameters:" << std::endl;
	std::cout << "  minimal exact match length       : " << options.lengthExact << std::endl;
	std::cout << "  minimal match length             : " << options.minLength << std::endl;
	std::cout << "  maximal error rate (epsilon)     : " << options.epsilon << std::endl;
	std::cout << "  maximal x-drop                   : " << options.xDrop << std::endl;
	if (options.qGram != (unsigned)-1)
		std::cout << "  k-mer (q-gram) length            : " << options.qGram << std::endl;
	std::cout << "  search forward strand            : " << ((options.forward)?"yes":"no") << std::endl;
	std::cout << "  search reverse complement        : " << ((options.reverse)?"yes":"no") << std::endl;
	std::cout << std::endl;

//	std::cout << "  verification strategy            : " << options.fastOption << std::endl;
//	if (options.disableThresh != (unsigned)-1) {
//		std::cout << "  disable queries with more than   : " << options.disableThresh << " matches" << std::endl;
//	}
//	std::cout << "  maximal number of matches        : " << options.numMatches << std::endl;
//	std::cout << "  duplicate removal every          : " << options.compactThresh << std::endl;
//	if (options.maxRepeatPeriod != 1 || options.minRepeatLength != 1000) {
//		std::cout << "  max low complexity repeat period : " << options.maxRepeatPeriod << std::endl;
//		std::cout << "  min low complexity repeat length : " << options.minRepeatLength << std::endl;
//	}
//	if (options.qgramAbundanceCut != 1) {
//		std::cout << "  q-gram abundance cut ratio       : " << options.qgramAbundanceCut << std::endl;
//	}
//	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes file name from options object to std::cout
template<typename TOptions>
void
_writeFileNames(TOptions & options) {
//IOREV _notio_
	std::cout << "Database file   : " << options.databaseFile << std::endl;
	std::cout << "Query file      : " << options.queryFile << std::endl;
	std::cout << "Output file     : " << options.outputFile << std::endl;

	std::cout << std::endl;
}

inline void
_addVersion(CommandLineParser& parser) {
	::std::string rev = "$Revision: 8300 $";
	addVersionLine(parser, "Version 1.1 (February 5th 2011) SeqAn Revision: " + rev.substr(11, 4) + "");
}

///////////////////////////////////////////////////////////////////////////////
// Parses options from command line parser and writes them into options object
template<typename TParser, typename TOptions>
bool
_parseOptions(TParser & parser, TOptions & options) {
//IOREV _notio_
    // i/o options
	getOptionValueShort(parser, 'd', options.databaseFile);
    getOptionValueShort(parser, 'q', options.queryFile);
    if (isSetShort(parser, 'o')) getOptionValueShort(parser, 'o', options.outputFile);
 	if (isSetLong(parser, "lengthExact")) getOptionValueLong(parser, "lengthExact", options.lengthExact);
	if (isSetShort(parser, "le")) getOptionValueShort(parser, "le", options.lengthExact);

    if (isSetLong(parser, "minLength")) getOptionValueLong(parser, "minLength", options.minLength);
	if (isSetShort(parser, 'e')) getOptionValueShort(parser, 'e', options.epsilon);
    if (isSetShort(parser, 'x')) getOptionValueShort(parser, 'x', options.xDrop);
	if (isSetShort(parser, 'v')) options.verbose = 1;
	if (isSetShort(parser, 'f')) if (!isSetShort(parser, 'r')) options.reverse = false;
	if (isSetShort(parser, 'r')) if (!isSetShort(parser, 'f')) options.forward = false;

	if (isSetShort(parser, "le") && options.lengthExact > 32) {
		std::cerr << "Invalid parameter value: Please choose a smaller k-mer length." << std::endl;
		return 0;
	}

	if (options.epsilon < 0.0) {
		std::cerr << "Invalid parameter value: Please choose a greater error rate." << std::endl;
		return 0;
	}

	if (options.epsilon > 0.25) {
		std::cerr << "Invalid parameter value: Please choose a smaller error rate." << std::endl;
		return 0;
	}
	if (isSetShort(parser, 'k') && options.qGram >= 1/options.epsilon) {
		std::cerr << "Invalid parameter value: Please choose q-gram length lower than 1/epsilon." << std::endl; 
		return 0;
	}
	return 1;
}






///////////////////////////////////////////////////////////////////////////////
// Set-Up of Command Line Parser
template<typename TParser>
void
_setParser(TParser & parser) {
	_addVersion(parser);

    addTitleLine(parser, "***************************************************");
	addTitleLine(parser, "*  prob_spec - Computing special local alignments *");
	addTitleLine(parser, "*              to assess probe specificity        *");
	addTitleLine(parser, "* (c) Copyright 2011 by Knut Reinert              *");
	addTitleLine(parser, "***************************************************");

	addUsageLine(parser, "-d <FASTA sequence file> -q <FASTA sequence file> [Options]");

	addLine(parser, "");
    addLine(parser, "An implementation a two pass filter based on a q-gram index and STELLAR");
 
	addSection(parser, "Non-optional Arguments:");
    addOption(parser, CommandLineOption('d', "database", "Fasta file containing the database sequences",
              (OptionType::String | OptionType::Mandatory)));
    addOption(parser, CommandLineOption('q', "query", "Fasta file containing the query sequences", 
              (OptionType::String | OptionType::Mandatory)));

	addSection(parser, "Main Options:");
    addOption(parser, CommandLineOption('e', "epsilon", "Maximal error rate (max 0.25)", OptionType::Double, "0.05"));
    addOption(parser, CommandLineOption('l', "minLength", "Minimal length of epsilon-matches", OptionType::Int, 100));
    addOption(parser, CommandLineOption('v', "verbose", "Verbosity mode.", OptionType::Bool, "false"));
	addOption(parser, CommandLineOption('f', "forward", "Search only in forward strand of database",
										OptionType::Boolean, "both"));
	addOption(parser, CommandLineOption('r', "reverse", "Search only in reverse complement of database",
										OptionType::Boolean, "both"));
	
	addSection(parser, "Filtering Options:");

  	addOption(parser, CommandLineOption("le", "lengthExact",
										"length of required exact submatch", OptionType::Int, "12"));


	addSection(parser, "Verification Options:");
    addOption(parser, CommandLineOption('x', "xDrop", "Maximal x-drop for extension", OptionType::Double, 5));

	addSection(parser, "Output Options:");
    addOption(parser, CommandLineOption('o', "out", "Name of output file", OptionType::String, "stellar.gff"));

}




///////////////////////////////////////////////////////////////////////////////
// Parses and outputs parameters, calls _stellarOnAll
int main(int argc, const char *argv[]) {

    // command line parsing
    CommandLineParser parser("prob_spec");

    _setParser(parser);
    if (!parse(parser, argc, argv)) {
		if (isSetShort(parser, 'h') || isSetShort(parser, 'v')) return 0; 
        shortHelp(parser, std::cerr);
        return 1;
    }

	ProbSpecOptions options = ProbSpecOptions();
	if (!_parseOptions(parser, options)) {
		return 1;
	}

	typedef String<Dna5> TSequence;
	typedef Iterator<TSequence>::Type TSequenceIt;
	
	// output header
	_title(parser, std::cout);
	std::cout << std::endl;

	// output file names
	_writeFileNames(options);

	// output parameters
	_writeSpecifiedParams(options);


    // import query sequences
    StringSet<TSequence> queries;
    StringSet<CharString> queryIDs;
	if (!_importSequences(options.queryFile, "query", queries, queryIDs)) return 1;

    // import database sequence
    StringSet<TSequence > databases;
    StringSet<CharString> databaseIDs;
	if (!_importSequences(options.databaseFile, "database", databases, databaseIDs)) return 1;
	
	// allocate result Stringset. it contains the alignments for each query
	StringSet< String<Align<TSequence> > > resultf, resultr;
	resize(resultf,length(queries));
	resize(resultr,length(queries));

	if (options.forward){	
		std::cout << "Forward strand" << ::std::endl;
		if( ! _extendExactMatches(databases[0],queries,resultf,options)) return 1;
	}
	
	if (options.reverse){ 
		std::cout << "Reverse strand" << ::std::endl;
		reverseComplement(databases[0]);
		if( ! _extendExactMatches(databases[0],queries,resultr,options)) return 1;
		reverseComplement(databases[0]);
	}
	
	std::ofstream file(toCString(options.outputFile));	

	if (options.forward){		
		file << ::std::endl << "Reverse strand" << ::std::endl;
		for (unsigned i=0; i<length(resultf); ++i)
			for (unsigned j=0; j<length(resultf[i]); j++) {	
				file << "hits for query " << i << ::std::endl;
				unsigned db1 = clippedBeginPosition(row(resultf[i][j],0));
				unsigned db2 = clippedEndPosition(row(resultf[i][j],0))-1;
				unsigned q1  = clippedBeginPosition(row(resultf[i][j],1));
				unsigned q2  = clippedEndPosition(row(resultf[i][j],1))-1;
				file << "Aligns database [" << db1 << ":" << db2<< "]" << " and query " << i << "[" << q1 << ":" <<  q2 << "]" << ::std::endl; 
				file << resultf[i][j] << ::std::endl;
			}
	
		for (unsigned i=0; i<length(resultf); ++i)	
			file << "query " << i << " has " << length(resultf[i]) << " hits in database " <<  ::std::endl;
	}
	
	if (options.reverse){ 
		file << ::std::endl << "Reverse strand" << ::std::endl;
		reverseComplement(databases[0]);
		for (unsigned i=0; i<length(resultr); ++i)
			for (unsigned j=0; j<length(resultr[i]); j++) {
				file << "hits for query " << i << ::std::endl;
				unsigned db1 = clippedBeginPosition(row(resultr[i][j],0));
				unsigned db2 = clippedEndPosition(row(resultr[i][j],0))-1;
				unsigned q1  = clippedBeginPosition(row(resultr[i][j],1));
				unsigned q2  = clippedEndPosition(row(resultr[i][j],1))-1;
				file << "Aligns database [" << db1 << ":" << db2<< "]" << " and query " << i << "[" << q1 << ":" <<  q2 << "]" << ::std::endl; 
				file << resultr[i][j] << ::std::endl;
			}
		for (unsigned i=0; i<length(resultr); ++i)	
			file << "query " << i << " has " << length(resultr[i]) << " hits in database " <<  ::std::endl;
		reverseComplement(databases[0]);
	}
	file.close();
	
	return 0;
}
