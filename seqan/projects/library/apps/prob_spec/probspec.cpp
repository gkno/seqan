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

#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>
#include "stellar.h"


using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
// Initializes a Finder object for a database sequence,
//  calls stellar, and writes matches to file

/*
 
 template<typename TSequence, typename TId, typename TPattern, typename TMatches>
inline bool
_stellarOnOne(TSequence & database,
			  TId & databaseID,
			  TPattern & swiftPattern,
			  bool databaseStrand,
			  TMatches & matches,
			  ProbSpecOptions & options) {
	std::cout << "  " << databaseID;
	if (!databaseStrand) std::cout << ", complement";

	// finder
    typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
	TFinder swiftFinder(database, options.minRepeatLength, options.maxRepeatPeriod);

	// stellar
	if (options.fastOption == CharString("exact"))
		stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
							   options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
							   databaseID, databaseStrand, matches, AllLocal());
	else if (options.fastOption == "bestLocal")
		stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
							   options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
							   databaseID, databaseStrand, matches, BestLocal());
	else if (options.fastOption == "bandedGlobal")
		stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
							   options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
							   databaseID, databaseStrand, matches, BandedGlobal());
	else if (options.fastOption == "bandedGlobalExtend")
		stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
							   options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
							   databaseID, databaseStrand, matches, BandedGlobalExtend());
	else {
		std::cerr << "\nUnknown verification strategy: " << options.fastOption << std::endl;
		return false;
	}

	std::cout << std::endl;
	return true;
}
*/
//////////////////////////////////////////////////////////////////////////////
namespace SEQAN_NAMESPACE_MAIN
{

	/*
template <typename TStringSet, typename TShape, typename TSpec>
struct Cargo<Index<TStringSet, IndexQGram<TShape, TSpec> > > {
	typedef struct {
		double		abundanceCut;
	} Type;
};


//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TStringSet, typename TShape, typename TSpec>
inline bool _qgramDisableBuckets(Index<TStringSet, IndexQGram<TShape, TSpec> > &index) 
{
	typedef Index<TStringSet, IndexQGram<TShape, TSpec> >	TReadIndex;
	typedef typename Fibre<TStringSet, QGramDir>::Type		TDir;
	typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
	typedef typename Value<TDir>::Type						TSize;

	TDir &dir    = indexDir(index);
	bool result  = false;
	unsigned counter = 0;
	TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
	if (thresh < 100) thresh = 100;

	TDirIterator it = begin(dir, Standard());
	TDirIterator itEnd = end(dir, Standard());
	for (; it != itEnd; ++it)
		if (*it > thresh) 
		{
			*it = (TSize)-1;
			result = true;
			++counter;
		}

	if (counter > 0)
		std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

	return result;
}
*/
}
	 
	 /*
///////////////////////////////////////////////////////////////////////////////
// Initializes a Pattern object with the query sequences, 
//  and calls _stellarOnOne for each database sequence
template<typename TSequence, typename TId>
inline bool
_stellarOnAll(StringSet<TSequence> & databases,
			  StringSet<TId> & databaseIDs,
			  StringSet<TSequence> & queries,
			  StringSet<TId> & queryIDs,
			  StringSet<QueryMatches<StellarMatch<TSequence, TId> > >& matches,
			   ProbSpecOptions & options) {
    // pattern
    typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
    TQGramIndex qgramIndex(queries);
 
	resize(indexShape(qgramIndex), options.qGram);

	cargo(qgramIndex).abundanceCut = options.qgramAbundanceCut;
	Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);
	
	if (options.verbose) swiftPattern.params.printDots = true;

	// Construct index
	std::cout << "Constructing index..." << std::endl;
	indexRequire(qgramIndex, QGramSADir());
	std::cout << std::endl;

	std::cout << "Aligning query sequence to database sequences..." << std::endl;
	for(unsigned i = 0; i < length(databases); ++i) {
		// positive database strand
		if (options.forward) {
			if (!_stellarOnOne(databases[i], databaseIDs[i], swiftPattern, true, matches, options))
				return 1;
		}
		// negative (reverse complemented) database strand
		if (options.reverse) {
			reverseComplement(databases[i]);
			if (!_stellarOnOne(databases[i], databaseIDs[i], swiftPattern, false, matches, options))
				return 1;
			reverseComplement(databases[i]);
		}
	}
	std::cout << std::endl;
	
	
	return 0;
}
*/


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
// Calculates parameters from parameters in options object and from sequences and writes them to std::cout

/*
template<typename TStringSet>
void _writeMoreCalculatedParams(ProbSpecOptions & options, TStringSet & databases, TStringSet & queries) {
//IOREV _notio_
	typedef typename Size<TStringSet>::Type TSize;

	if (options.qgramAbundanceCut != 1) {
		std::cout << "Calculated parameters:" << std::endl;
	}

	TSize queryLength = length(concat(queries));
	if (options.qgramAbundanceCut != 1) {
		std::cout << "  q-gram expected abundance : ";
		std::cout << queryLength/(double)((long)1<<(options.qGram<<1)) << std::endl;
		std::cout << "  q-gram abundance threshold: ";
		std::cout << _max(100,(int)(queryLength*options.qgramAbundanceCut)) << std::endl;
		std::cout << std::endl;
	}

	// Computation of maximal E-value for this search

	TSize maxLengthQueries = 0;
	TSize maxLengthDatabases = 0;

	typename Iterator<TStringSet>::Type dbIt = begin(databases);
	typename Iterator<TStringSet>::Type dbEnd = end(databases);
	while (dbIt != dbEnd) {
		if (length(*dbIt) > maxLengthDatabases) {
			maxLengthDatabases = length(*dbIt);
		}
		++dbIt;
	}

	typename Iterator<TStringSet>::Type queriesIt = begin(queries);
	typename Iterator<TStringSet>::Type queriesEnd = end(queries);
	while (queriesIt != queriesEnd) {
		if (length(*queriesIt) > maxLengthQueries) {
			maxLengthQueries = length(*queriesIt);
		}
		++queriesIt;
	}

	TSize errors = static_cast<TSize>(options.minLength * options.epsilon);
	TSize minScore = options.minLength - 3*errors; // #matches - 2*#errors // #matches = minLenght - errors, 

	std::cout << "All matches matches resulting from your search have an E-value of: " << std::endl;
	std::cout << "        " << _computeEValue(minScore, maxLengthQueries, maxLengthDatabases) << " or smaller";
	std::cout << "  (match score = 1, error penalty = -2)" << std::endl;

	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and writes them to std::cout
void _writeCalculatedParams(ProbSpecOptions & options) {
//IOREV _notio_
	int errMinLen = (int) floor(options.epsilon * options.minLength);
	int n = (int) ceil((errMinLen + 1) / options.epsilon);
	int errN = (int) floor(options.epsilon * n);
	unsigned smin = (unsigned) _min(ceil((double)(options.minLength-errMinLen)/(errMinLen+1)),
		                            ceil((double)(n-errN)/(errN+1)));

	std::cout << "Calculated parameters:" << std::endl;
	if (options.qGram == (unsigned)-1) {
		options.qGram = (unsigned)_min(smin, 32u);
		std::cout << "  k-mer length: " << options.qGram << std::endl;
	}

	int threshold = (int) _max(1, (int) _min((n + 1) - options.qGram * (errN + 1),
											 (options.minLength + 1) - options.qGram * (errMinLen + 1)));
	int overlap = (int) floor((2 * threshold + options.qGram - 3) / (1 / options.epsilon - options.qGram));
	int distanceCut = (threshold - 1) + options.qGram * overlap + options.qGram;
	int logDelta = _max(4, (int) ceil(log((double)overlap + 1) / log(2.0)));
	int delta = 1 << logDelta;

	std::cout << "  s^min       : " << smin << std::endl;
	std::cout << "  threshold   : " << threshold << std::endl;
	std::cout << "  distance cut: " << distanceCut << std::endl;
	std::cout << "  delta       : " << delta << std::endl;
	std::cout << "  overlap     : " << overlap << std::endl;
	std::cout << std::endl;
}
*/

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
//	std::cout << "  search forward strand            : " << ((options.forward)?"yes":"no") << std::endl;
//	std::cout << "  search reverse complement        : " << ((options.reverse)?"yes":"no") << std::endl;
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
//	std::cout << "Output format   : " << options.outputFormat << std::endl;
//	if (options.disableThresh != (unsigned)-1) {
//		std::cout << "Disabled queries: " << options.disabledQueriesFile << std::endl;
//	}
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
  //  if (isSetShort(parser, "od")) getOptionValueShort(parser, "od", options.disabledQueriesFile);
  //	if (isSetShort(parser, "of")) getOptionValueShort(parser, "of", options.outputFormat);

	// main options
	if (isSetLong(parser, "kmer")) getOptionValueLong(parser, "kmer", options.qGram);
	if (isSetLong(parser, "lengthExact")) getOptionValueLong(parser, "lengthExact", options.lengthExact);
	if (isSetShort(parser, "le")) getOptionValueShort(parser, "le", options.lengthExact);

    if (isSetLong(parser, "minLength")) getOptionValueLong(parser, "minLength", options.minLength);
	if (isSetShort(parser, 'e')) getOptionValueShort(parser, 'e', options.epsilon);
    if (isSetShort(parser, 'x')) getOptionValueShort(parser, 'x', options.xDrop);

//	if (isSetShort(parser, 'f')) if (!isSetShort(parser, 'r')) options.reverse = false;
// if (isSetShort(parser, 'r')) if (!isSetShort(parser, 'f')) options.forward = false;

//	if (isSetShort(parser, "vs")) getOptionValueShort(parser, "vs", options.fastOption);
//	if (isSetShort(parser, "dt")) getOptionValueShort(parser, "dt", options.disableThresh);
//	if (isSetShort(parser, 'n')) getOptionValueShort(parser, 'n', options.numMatches);
//	if (isSetShort(parser, 's')) getOptionValueShort(parser, 's', options.compactThresh);
//	if (isSetShort(parser, "rp")) getOptionValueShort(parser, "rp", options.maxRepeatPeriod);
//	if (isSetShort(parser, "rl")) getOptionValueShort(parser, "rl", options.minRepeatLength);
//	if (isSetShort(parser, 'a')) getOptionValueShort(parser, 'a', options.qgramAbundanceCut);

	if (isSetShort(parser, 'v')) options.verbose = 1;

//	if (options.outputFormat != "gff" && options.outputFormat != "text") {
//		std::cerr << "Invalid parameter value: Unknown output format." << std::endl;
//		return 0;
//	}
/*
	if (options.fastOption != "exact" && options.fastOption != "bestLocal"
		 && options.fastOption != "bandedGlobal") {
		std::cerr << "Invalid parameter value: Unknown verification strategy." << std::endl;
		return 0;
	}
*/	if (isSetShort(parser, 'k') && options.qGram < 1) {
		std::cerr << "Invalid parameter value: Please choose a greater k-mer length." << std::endl;
		return 0;
	}

	if (isSetShort(parser, 'k') && options.qGram > 32) {
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
/*
	if (options.qgramAbundanceCut > 1 || options.qgramAbundanceCut < 0) {
		std::cerr << "Invalid parameter value: Please choose a k-mer overabundance cut ration between 0 and 1.\n";
		return 0;
	}

	if (options.numMatches > options.compactThresh) {
		std::cerr << "Invalid parameter values: Please choose numMatches <= sortThresh." << std::endl;
		return 0;
	}
 */
	return 1;
}





template<typename TInfix, typename TSize, typename TEps, typename TAlign>
bool
_extendKmer(TInfix & a, TInfix & b, TSize minLength, TEps eps, TAlign & align) {
	typedef Seed<Simple> TSeed;
	typedef int TScore;
	
	TSeed seed(beginPosition(a), beginPosition(b), endPosition(a), endPosition(b));
	TSeed seedOld(seed);
	
	TScore penalty = static_cast<TScore>(ceil(-1/eps) + 1);
	Score<TScore> scoreMatrix(1, penalty, penalty);
	TScore scoreDropOff = -penalty * static_cast<TScore>(minLength * eps);

	
	ExtensionDirection direction = EXTEND_BOTH;
	extendSeed(seed, host(a), host(b), direction, scoreMatrix, scoreDropOff, GappedXDrop());
	
	return _bestExtension(a, b, seed, seedOld, unsigned(length(a)), 0u, scoreMatrix, direction, minLength, eps, align);
}




template<typename TSequence>
inline bool
_extendExactMatches(TSequence & database,
			  StringSet<TSequence> & queries,
			  StringSet<String<Align<TSequence > > >& result,
			ProbSpecOptions & options) {
	
	StringSet<String<Triple<unsigned,unsigned,unsigned> > > ematches;
	resize(ematches, length(queries));
	
	// built qGram index on database
	typedef Index<TSequence, IndexQGram<SimpleShape> > TQGramIndex;
    TQGramIndex qgramIndex(database);
	resize(indexShape(qgramIndex), options.lengthExact);
	
	indexRequire(qgramIndex, QGramSADir());
	indexRequire(qgramIndex, QGramDir());
 
	typedef typename Iterator<TSequence>::Type TSequenceIt;

	for (unsigned int s=0; s < length(queries); s++) {
	//::std::cout << "Inspecting query " << s << " " << queries[s] << std::endl;
		
		TSequenceIt sit  = begin(queries[s]);
		TSequenceIt send = end(queries[s])-options.lengthExact;
		
		unsigned hv = hash(indexShape(qgramIndex),sit);
		typedef typename Infix<typename Fibre<TQGramIndex, QGramSA>::Type const>::Type TInfix;
		TInfix occs;
		
		int qp=0;
		for (; sit != send; ) { 
			occs = getOccurrences(qgramIndex, indexShape(qgramIndex));
			for (unsigned i = 0; i < length(occs); i++){
		//		::std::cout << s << " " << occs[i] << ::std::endl;
				Triple<unsigned, unsigned, unsigned> p(occs[i],s,qp);
		//		::std::cout << p << ::std::endl;
				appendValue(ematches[s], p);
			}
			sit++;
			qp++;
			hv = hashNext(indexShape(qgramIndex),sit);
/*			
			TSequence result;
			unhash(result, hv, options.lengthExact);	
			::std::cout << result <<  ::std::endl;
*/			
		}
	}
	
	// remove adjacent q-gram hits
	StringSet< String<Triple<unsigned,unsigned,unsigned> > > cematches;
	resize(cematches,length(ematches));
	
	
	for (unsigned s=0; s<length(ematches); s++) {
		unsigned cind = 0;		
		if( length(ematches[s]) > 0 ){
			std::sort(begin(ematches[s]), end(ematches[s]));
			resize(cematches[s],1);
			cematches[s][cind++] = ematches[s][0];
		//	::std::cout << "CMATCH "<< ematches[s][0] << ::std::endl;
			
		}
		if( length(ematches[s]) > 1 ){
			for (unsigned i=1; i<length(ematches[s]); ++i) {
				if( ematches[s][i].i2 != ematches[s][i-1].i2 || ematches[s][i].i1 != ematches[s][i-1].i1+1 ){
					appendValue(cematches[s],ematches[s][i]);
					::std::cout << "CMATCH "<< ematches[s][i] << ::std::endl;
				}
			}
		}
		//		::std::cout << "size of cematches " << cind <<  " " << length(cematches[s]) << std::endl;
	}
	
	
	typedef typename Infix<TSequence>::Type TInfix;	
	TSequenceIt dbbegin = begin(database);

	Align<TSequence>  oldalign;
	resize(rows(oldalign), 2);

	
	for (unsigned s=0; s<length(cematches); s++) {
		TSequenceIt qbegin = begin(queries[s]);

		// there are extensions to be made
		for (unsigned i=0; i<length(cematches[s]); i++) {
			unsigned dbpos = cematches[s][i].i1;
			unsigned qpos = cematches[s][i].i3;

			TInfix dbinfix = infix(database,dbbegin + dbpos,dbbegin + dbpos + options.lengthExact);
			TInfix qinfix = infix(queries[s],qbegin + qpos,qbegin + qpos + options.lengthExact);
	
			::std::cout << "positions " << dbpos << " " << qpos << ::std::endl;	
			::std::cout << dbinfix << ::std::endl;
			::std::cout << qinfix << ::std::endl;
	
			Align<TSequence> align;
			resize(rows(align), 2);
			setSource(row(align, 0), host(dbinfix));
			setSource(row(align, 1), host(qinfix));
				
	
			_extendKmer(dbinfix, qinfix, options.lengthExact, options.epsilon, align);
			
		
			if(length(row(align,0))  > options.minLength){
				unsigned db1 = clippedBeginPosition(row(align,0));
				unsigned db2 = clippedEndPosition(row(align,0))-1;
				unsigned q1  = clippedBeginPosition(row(align,1));
				unsigned q2  = clippedEndPosition(row(align,1))-1;
				unsigned odb1 = clippedBeginPosition(row(oldalign,0));
				unsigned odb2 = clippedEndPosition(row(oldalign,0))-1;
				unsigned oq1  = clippedBeginPosition(row(oldalign,1));
				unsigned oq2  = clippedEndPosition(row(oldalign,1))-1;
		
		//		::std::cout << "Old Aligns Seq1[" << odb1 << ":" << odb2<< "]" << " and Seq2[" << oq1 << ":" <<  oq2 << "]";
		//			::std::cout << ::std::endl; 
				
			//	::std::cout << oldalign << ::std::endl;
				
		//		::std::cout << "Aligns Seq1[" << db1 << ":" << db2<< "]" << " and Seq2[" << q1 << ":" <<  q2 << "]";
		//		::std::cout << ::std::endl; 
		
		//		::std::cout << align << ::std::endl;
	
				if (odb1 != db1 || odb2 != db2 || oq1 != q1 || oq2 != q2 ) {
					appendValue(result[s],align);
					oldalign = align;
				}
		//		else	
		//			::std::cout << "double match " << ::std::endl;
				
			}
		}
	}

			
	
	return true;
}


///////////////////////////////////////////////////////////////////////////////
// Set-Up of Command Line Parser
template<typename TParser>
void
_setParser(TParser & parser) {
	_addVersion(parser);

    addTitleLine(parser, "***************************************************");
	addTitleLine(parser, "*  prob_spec - Computing special local alignments *");
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
//	addOption(parser, CommandLineOption('f', "forward", "Search only in forward strand of database",
//		OptionType::Boolean, "both"));
//	addOption(parser, CommandLineOption('r', "reverse", "Search only in reverse complement of database",
//		OptionType::Boolean, "both"));
    addOption(parser, CommandLineOption('v', "verbose", "Verbosity mode.", OptionType::Bool, "false"));
    
	addSection(parser, "Filtering Options:");
    addOption(parser, CommandLineOption('k', "kmer", "Length of the q-grams (max 32)", OptionType::Int, 10));
  //  addOption(parser, CommandLineOption("rp", "repeatPeriod",
  //	"Maximal period of low complexity reapeats to be filtered", OptionType::Int, 1));
  //  addOption(parser, CommandLineOption("rl", "repeatLength",
  //		"Minimal length of low complexity reapeats to be filtered", OptionType::Int, 1000));
  //  addOption(parser, CommandLineOption('a', "abundanceCut",
  //		"k-mer overabundance cut ratio", OptionType::Double, "1"));
  	addOption(parser, CommandLineOption("le", "lengthExact",
										"length of required exact submatch", OptionType::Int, "12"));


	addSection(parser, "Verification Options:");
    addOption(parser, CommandLineOption('x', "xDrop", "Maximal x-drop for extension", OptionType::Double, 5));
//	addOption(parser, CommandLineOption("vs", "verification", "Verification strategy", OptionType::String, "exact"));
//	addHelpLine(parser, "exact        = compute and extend all local alignments in SWIFT hits");
//	addHelpLine(parser, "bestLocal    = compute and extend only best local alignment in SWIFT hits");
//	addHelpLine(parser, "bandedGlobal = banded global alignment on SWIFT hits");
//	addOption(parser, CommandLineOption("dt", "disableThresh",
//		"Maximal number of verified matches before disabling verification", OptionType::Int));
//	addHelpLine(parser, "for one query sequence (default infinity)");
//	addOption(parser, CommandLineOption('n', "numMatches",
//		"Maximal number of kept matches per query and database", OptionType::Int, 50));
//	addHelpLine(parser, "If there are more matches, only the longest ones are kept.");
//	addOption(parser, CommandLineOption('s', "sortThresh",
//		"Number of matches triggering removal of duplicates", OptionType::Int, 500));
//	addHelpLine(parser, "Choose a smaller value for saving space.");

	addSection(parser, "Output Options:");
    addOption(parser, CommandLineOption('o', "out", "Name of output file", OptionType::String, "stellar.gff"));
//	addOption(parser, CommandLineOption("of", "outFormat", "Output format", OptionType::String, "gff"));
//	addHelpLine(parser, "Possible formats: gff, text");
//	addOption(parser, CommandLineOption("od", "outDisabled",
//		"Name of output file containing disabled query sequences", OptionType::String));
//	addHelpLine(parser, "(default stellar.disabled.fasta)");
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
//	_writeCalculatedParams(options);

    // import query sequences
    StringSet<TSequence> queries;
    StringSet<CharString> queryIDs;
	if (!_importSequences(options.queryFile, "query", queries, queryIDs)) return 1;

    // import database sequence
    StringSet<TSequence > databases;
    StringSet<CharString> databaseIDs;
	if (!_importSequences(options.databaseFile, "database", databases, databaseIDs)) return 1;
	
	// allocate result Stringset. it contains the alignments for each query
	StringSet< String<Align<TSequence> > > result;
	resize(result,length(queries));

	if( ! _extendExactMatches(databases[0],queries,result,options)) return 1;
	
	for (unsigned i=0; i<length(result); ++i)
		for (unsigned j=0; j<length(result[i]); j++) {
			::std::cout << "hits for query " << i << ::std::endl;
			unsigned db1 = clippedBeginPosition(row(result[i][j],0));
			unsigned db2 = clippedEndPosition(row(result[i][j],0))-1;
			unsigned q1  = clippedBeginPosition(row(result[i][j],1));
			unsigned q2  = clippedEndPosition(row(result[i][j],1))-1;
			::std::cout << "Aligns database [" << db1 << ":" << db2<< "]" << " and query " << i << "[" << q1 << ":" <<  q2 << "]" << ::std::endl; 
			::std::cout << result[i][j] << ::std::endl;
		}
		for (unsigned i=0; i<length(result); ++i)	
			::std::cout << "query " << i << " has " << length(result[i]) << " hits in database " <<  ::std::endl;
	
	return 0;
}
