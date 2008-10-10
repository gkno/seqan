 /*==========================================================================
                     RazerS - Fast Mapping of Short Reads
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

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

#ifndef SEQAN_HEADER_RAZERS_H
#define SEQAN_HEADER_RAZERS_H

#include <iostream>
#include <fstream>
#include <sstream>

#include <seqan/find.h>
#include <seqan/align.h>
#include <seqan/index.h>
#include <seqan/find/find_swift.h>
#include "mmap_fasta.h"

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Default options

	template < bool _DONT_VERIFY = false, bool _DONT_DUMP_RESULTS = false >
	struct RazerSSpec 
	{
		enum { DONT_VERIFY = _DONT_VERIFY };				// omit verifying potential matches
		enum { DONT_DUMP_RESULTS = _DONT_DUMP_RESULTS };	// omit dumping results
	};

	template < typename TSpec = RazerSSpec<> >
	struct RazerSOptions
	{
	// main options
		TSpec		spec;
		bool		forward;			// compute forward oriented read matches
		bool		reverse;			// compute reverse oriented read matches
		double		errorRate;			// Criteria 1 threshold
		int			maxHits;			// hit count threshold
		const char	*output;			// name of result file
		int			_debugLevel;		// level of verbosity
		bool		printVersion;		// print version number
		bool		hammingOnly;		// no indels

	// output format options
		unsigned	outputFormat;		// 0..Razer format
										// 1..enhanced Fasta
		bool		dumpAlignment;		// compute and dump the match alignments in the result files
		unsigned	genomeNaming;		// 0..use Fasta id
										// 1..enumerate reads beginning with 1
		unsigned	readNaming;			// 0..use Fasta id
										// 1..enumerate reads beginning with 1
										// 2..use the read sequence (only for short reads!)
		unsigned	sortOrder;			// 0..sort keys: 1. read number, 2. genome position
										// 1..           1. genome pos50ition, 2. read number
		unsigned	positionFormat;		// 0..gap space
										// 1..position space

	// filtration parameters
		::std::string shape;			// shape (e.g. 11111111111)
		int			threshold;			// threshold
		int			tabooLength;		// taboo length
		int			repeatLength;		// repeat length threshold
		double		abundanceCut;		// abundance threshold

	// verification parameters
		bool		matchN;				// false..N is always a mismatch, true..N matches with all
		unsigned char compMask[5];

	// statistics
		__int64		FP;					// false positives (threshold reached, no match)
		__int64		TP;					// true positives (threshold reached, match)
		double		timeLoadFiles;		// time for loading input files
		double		timeMapReads;		// time for mapping reads
		double		timeDumpResults;	// time for dumping the results

		RazerSOptions() 
		{
			forward = true;
			reverse = true;
			errorRate = 0.08;
			maxHits = 100;
			output = "";
			_debugLevel = 0;
			printVersion = false;
			hammingOnly = false;

			outputFormat = 0;
			dumpAlignment = false;
			genomeNaming = 0;
			readNaming = 0;
			sortOrder = 0;
			positionFormat = 0;

			matchN = false;

			shape = "11111111111";
			threshold = 1;
			tabooLength = 1;
			repeatLength = 1000;
			abundanceCut = 1;
			
			for (unsigned i = 0; i < 4; ++i)
				compMask[i] = 1 << i;
			compMask[4] = 0;
		}
	};




//////////////////////////////////////////////////////////////////////////////
// Typedefs

	// definition of a Read match
	template <typename _TGPos>
	struct ReadMatch 
	{
		typedef typename _MakeSigned<_TGPos>::Type TGPos;

		unsigned		gseqNo;			// genome seqNo
		unsigned		rseqNo;			// read seqNo
		TGPos			gBegin;			// begin position of the match in the genome
		TGPos			gEnd;			// end position of the match in the genome
		short			editDist;		// Levenshtein distance
		char			orientation;	// 'F'..forward strand, 'R'..reverse comp. strand
	};

	// less-operators (to sort matches and remove duplicates with equal gBegin)
	template <typename TReadMatch>
	struct LessGPosRNo : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};

	template <typename TReadMatch>
	struct LessRNoGPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};

	// less-operators (to sort matches and remove duplicates with equal gEnd)
	template <typename TReadMatch>
	struct LessRNoGEndPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gEnd   < b.gEnd) return true;
			if (a.gEnd   > b.gEnd) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};

//////////////////////////////////////////////////////////////////////////////
// Definitions

	typedef Dna5String									TGenome;
	typedef StringSet<TGenome>							TGenomeSet;
	typedef Dna5String									TRead;
#ifdef RAZERS_CONCATREADS
	typedef StringSet<TRead, Owner<ConcatDirect<> > >	TReadSet;
#else
	typedef StringSet<TRead>							TReadSet;
#endif

	typedef ReadMatch<Difference<TGenome>::Type>		TMatch;		// a single match
	typedef String<TMatch/*, MMap<>*/ >					TMatches;	// array of matches


	template <typename TShape>
	struct Cargo< Index<TReadSet, TShape> > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};

//////////////////////////////////////////////////////////////////////////////
// Memory tuning

#ifdef RAZERS_MEMOPT

	template <typename TShape>
	struct SAValue< Index<TReadSet, TShape> > {
		typedef Pair<
			unsigned,				
			unsigned,
			BitCompressed<25, 7>	// max. 32M reads of length < 128
		> Type;
	};
	
	template <>
	struct Size<Dna5String>
	{
		typedef unsigned Type;
	};

	template <typename TShape>
	struct Size< Index<TReadSet, Index_QGram<TShape> > >
	{
		typedef unsigned Type;
	};
	
#else

	template <typename TShape>
	struct SAValue< Index<TReadSet, TShape> > {
		typedef Pair<
			unsigned,			// many reads
			unsigned,			// of arbitrary length
			Compressed
		> Type;
	};
	
#endif

#ifdef RAZERS_PRUNE_QGRAM_INDEX

	//////////////////////////////////////////////////////////////////////////////
	// Repeat masker
	template <typename TShape>
	inline bool _qgramDisableBuckets(Index<TReadSet, Index_QGram<TShape> > &index) 
	{
		typedef Index<TReadSet, Index_QGram<TShape>	>		TReadIndex;
		typedef typename Fibre<TReadIndex, QGram_Dir>::Type	TDir;
		typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
		typedef typename Value<TDir>::Type					TSize;

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

		if (counter > 0 && cargo(index)._debugLevel >= 1)
			::std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

		return result;
	}

#endif


//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TReadSet, typename TNameSet, typename TRazerSOptions>
bool loadReads(TReadSet &reads, TNameSet &fastaIDs, const char *fileName, TRazerSOptions &options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);

	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
	split(multiFasta, Fasta());

	unsigned seqCount = length(multiFasta);
#ifndef RAZERS_CONCATREADS
	resize(reads, seqCount);
#endif
	if (options.readNaming == 0)
		resize(fastaIDs, seqCount);
	
	Dna5String seq;
	unsigned kickoutcount = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0)
			assignSeqId(fastaIDs[i], multiFasta[i], Fasta());	// read Fasta id
		assignSeq(seq, multiFasta[i], Fasta());					// read Read sequence
		
		if (countN)
		{
			int count = 0;
			int cutoffCount = (int)(options.errorRate * length(seq));
			for (unsigned j = 0; j < length(seq); ++j)
				if (getValue(seq, j) == 'N')
					if (++count > cutoffCount)
					{
						clear(seq);
						++kickoutcount;
						break;
					}
		}
#ifdef RAZERS_CONCATREADS
		appendValue(reads, seq, Generous());
#else
		assign(reads[i], seq, Exact());
#endif
	}
#ifdef RAZERS_CONCATREADS
	reserve(reads.concat, length(reads.concat), Exact());
#endif

	if (options._debugLevel > 1 && kickoutcount > 0) 
		::std::cerr << "Ignoring " << kickoutcount << " low quality reads.\n";
	return (seqCount > 0);
}


//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TGenomeSet>
bool loadGenomes(TGenomeSet &genomes, const char *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
	split(multiFasta, Fasta());

	unsigned seqCount = length(multiFasta);
	resize(genomes, seqCount);
	for(unsigned i = 0; i < seqCount; ++i)
		assignSeq(genomes[i], multiFasta[i], Fasta());		// read Genome sequence

	return (seqCount > 0);
}


//////////////////////////////////////////////////////////////////////////////
// Hamming verification
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadIndex, 
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadIndex &readIndex,					// q-gram index
	TMyersPatterns const &,					// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,	// RazerS options
	SwiftSemiGlobalHamming)					// Hamming only
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

#ifdef RAZERS_DEBUG
	cout<<"Verify: "<<::std::endl;
	cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	cout<<"Read:   "<<host(myersPattern)<<::std::endl;
#endif

	unsigned ndlLength = sequenceLength(rseqNo, readIndex);
	if (length(inf) < ndlLength) return false;

	// verify
	TRead &read				= indexText(readIndex)[rseqNo];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	TGenomeIterator git		= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);

	int maxErrors = (int)(ndlLength * options.errorRate);
	int minErrors = maxErrors + 1;

	for (; git < gitEnd; ++git)
	{
		int errors = 0;
		TGenomeIterator g = git;
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
				if (++errors > maxErrors)
					break;
		if (minErrors > errors)
		{
			minErrors = errors;
			m.gBegin = git - begin(host(inf), Standard());
		}
	}

	if (minErrors > maxErrors) return false;

	m.gEnd = m.gBegin + ndlLength;
	m.editDist = minErrors;
	return true;
}	


//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadIndex, 
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadIndex &readIndex,					// q-gram index
	TMyersPatterns &forwardPatterns,		// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,	// RazerS options
	SwiftSemiGlobal)						// Swift specialization
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Value<TReadSet>::Type					TRead;

	// find read match end
	typedef Finder<TGenomeInfix>							TMyersFinder;
	typedef typename Value<TMyersPatterns>::Type			TMyersPattern;

	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>		TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>				TReadRev;
	typedef Finder<TGenomeInfixRev>							TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>			TMyersPatternRev;

	TMyersFinder myersFinder(inf);
	TMyersPattern &myersPattern = forwardPatterns[rseqNo];

#ifdef RAZERS_DEBUG
	cout<<"Verify: "<<::std::endl;
	cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	cout<<"Read:   "<<host(myersPattern)<<::std::endl;
#endif

	unsigned ndlLength = sequenceLength(rseqNo, readIndex);
	int maxScore = InfimumValue<int>::VALUE;
	int minScore = -(int)(ndlLength * options.errorRate);
	TMyersFinder maxPos;

	// find end of best semi-global alignment
	while (find(myersFinder, myersPattern, minScore))
		if (maxScore < getScore(myersPattern)) 
		{
			maxScore = getScore(myersPattern);
			maxPos = myersFinder;
		}
	
	if (maxScore < minScore) return false;
	m.editDist	= -maxScore;

	setEndPosition(inf, beginPosition(inf) + position(maxPos) + 1);
	
	TGenomeInfixRev		infRev(inf);
	TReadRev			readRev(indexText(readIndex)[rseqNo]);
	TMyersFinderRev		myersFinderRev(infRev);
	TMyersPatternRev	myersPatternRev(readRev);

	_patternMatchNOfPattern(myersPatternRev, options.matchN);
	_patternMatchNOfFinder(myersPatternRev, options.matchN);

	// find beginning of best semi-global alignment
	if (find(myersFinderRev, myersPatternRev, maxScore))
	{
		m.gEnd = endPosition(inf);
		m.gBegin = m.gEnd - (position(myersFinderRev) + 1);
	} else {
		// this case should never occur
		::std::cerr << "1GENOME: " << host(myersFinder) << ::std::endl;
		::std::cerr << "1READ:   " << indexText(readIndex)[rseqNo] << ::std::endl;
		::std::cerr << "2GENOME: " << inf << '\t' << m.gBegin << ',' << m.gEnd << ::std::endl;
		::std::cerr << "2READ:   " << indexText(readIndex)[rseqNo] << ::std::endl;
		::std::cerr << "3GENOME: " << infRev << ::std::endl;
		::std::cerr << "3READ:   " << readRev << ::std::endl;
		::std::cerr << "HUH?" << ::std::endl;
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename THitCount, 
	typename TSpec >
void findReads(
	TMatches &matches,			// resulting matches
	TGenome &genome,			// genome ...
	unsigned gseqNo,			// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPattern,
	TVerifier &forwardPatterns,
	char orientation,			// q-gram index of reads
#ifdef RAZERS_MAXHITS
	THitCount &hitCount,		// maximum number of hits for each read
#else
	THitCount &,
#endif
	RazerSOptions<TSpec> &options)
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Value<TReadSet>::Type					TRead;
	typedef typename Value<TMatches>::Type					TMatch;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinder;

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
	if (orientation == 'F')
		::std::cerr << "[fwd]";
	else
		::std::cerr << "[rev]";

	TReadIndex &readIndex = host(swiftPattern);
	TSwiftFinder swiftFinder(genome, options.repeatLength, 1);

	TMatch m;
	m.gBegin = 0;
	
	// iterate all verification regions returned by SWIFT
	while (find(swiftFinder, swiftPattern, options.errorRate, options._debugLevel)) 
	{
		unsigned rseqNo = (*swiftFinder.curHit).ndlSeqNo;
#ifdef RAZERS_MAXHITS
		if (hitCount[rseqNo] == -1) continue;
#endif			
		if (!options.spec.DONT_VERIFY && 
			matchVerify(m, range(swiftFinder, genome), rseqNo, readIndex, forwardPatterns, options, TSwiftSpec()))
		{
			// transform coordinates to the forward strand
			if (orientation == 'R') 
			{
				TSize gLength = length(genome);
				TSize temp = m.gBegin;
				m.gBegin = gLength - m.gEnd;
				m.gEnd = gLength - temp;
			}
			m.gseqNo = gseqNo;
			m.rseqNo = rseqNo;
			m.orientation = orientation;

			if (!options.spec.DONT_DUMP_RESULTS)
				appendValue(matches, m);
#ifdef RAZERS_MAXHITS				
			--hitCount[rseqNo];
#endif
			++options.TP;
/*			::std::cerr << "\"" << range(swiftFinder, genomeInf) << "\"  ";
			::std::cerr << hstkPos << " + ";
			::std::cerr << ::std::endl;
*/		} else {
			++options.FP;
/*			::std::cerr << "\"" << range(swiftFinder, genomeInf) << "\"   \"" << range(swiftPattern) << "\"  ";
			::std::cerr << rseqNo << " : ";
			::std::cerr << hstkPos << " + ";
			::std::cerr << bucketWidth << "  " << TP << ::std::endl;
*/		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences
template <
	typename TMatches, 
	typename TReadSet, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
bool mapReads(
	TMatches &				matches,
	const char *			genomeFileName,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	TReadSet const &		readSet,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type				TRead;
	typedef Index<TReadSet, Index_QGram<TShape> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

	// open genome file	
	::std::ifstream file;
	file.open(genomeFileName, ::std::ios_base::in | ::std::ios_base::binary);
	if (!file.is_open()) {
		::std::cerr << "Failed to load genomes" << ::std::endl;
		return false;
	}

	// configure q-gram index
	TIndex swiftIndex(readSet, shape);
	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPattern(swiftIndex);
	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;

	// init edit distance verifiers
	String<TMyersPattern> forwardPatterns;
	if (options.hammingOnly)
		options.compMask[4] = (options.matchN)? 15: 0;
	else
	{
		unsigned readCount = countSequences(swiftIndex);
		resize(forwardPatterns, readCount);
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	String<int>	hitCount;
#ifdef RAZERS_MAXHITS	
	// fill max number of hits per read
	// at most twice as much potential matches are allowed (overlapping parallelograms)
	fill(hitCount, length(readSet), 2*options.maxHits);	
#endif

	CharString	id;
	Dna5String	genome;

	// iterate over genome sequences
	SEQAN_PROTIMESTART(find_time);
	for(unsigned gseqNo = 0; !_streamEOF(file); ++gseqNo)
	{
		if (options.genomeNaming == 0)
		{
			readID(file, id, Fasta());			// read Fasta id
			appendValue(genomeNames, id, Generous());
		}
		read(file, genome, Fasta());			// read Fasta sequence

		if (options.forward)
			findReads(matches, genome, gseqNo, swiftPattern, forwardPatterns, 'F', hitCount, options);

		if (options.reverse)
		{
			reverseComplementInPlace(genome);
			findReads(matches, genome, gseqNo, swiftPattern, forwardPatterns, 'R', hitCount, options);
		}

	}
	options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
	file.close();

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Swift FP: " << options.FP << ::std::endl;
		::std::cerr << "Swift TP: " << options.TP << ::std::endl;
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// Wrapper for different template specializations
template <typename TMatches, typename TReadSet, typename TSpec>
bool mapReads(
	TMatches &				matches,
	const char *			genomeFileName,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	TReadSet const &		readSet, 
	RazerSOptions<TSpec> &	options)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GappedShape>		gapped;

	// 2x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting shape
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeFileName, genomeNames, readSet, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeFileName, genomeNames, readSet, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeFileName, genomeNames, readSet, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeFileName, genomeNames, readSet, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeFileName, genomeNames, readSet, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeFileName, genomeNames, readSet, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return false;
}

//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int mapReads(
	const char *genomeFileName,
	const char *readFileName,
	const char *errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	MultiFasta				genomeSet;
	TReadSet				readSet;
	StringSet<CharString>	genomeNames;	// genome names, taken from the Fasta file
	StringSet<CharString>	readNames;		// read names, taken from the Fasta file
	TMatches				matches;		// resulting forward/reverse matches

	// dump configuration in verbose mode
	if (options._debugLevel >= 1) 
	{
//		CharString bitmap;
//		shapeToString(bitmap, shape);
		::std::cerr << "___SETTINGS____________" << ::std::endl;
		::std::cerr << "Genome file:                     \t" << genomeFileName << ::std::endl;
		::std::cerr << "Read file:                       \t" << readFileName << ::std::endl;
		::std::cerr << "Compute forward matches:         \t";
		if (options.forward)	::std::cerr << "YES" << ::std::endl;
		else				::std::cerr << "NO" << ::std::endl;
		::std::cerr << "Compute reverse matches:         \t";
		if (options.reverse)		::std::cerr << "YES" << ::std::endl;
		else				::std::cerr << "NO" << ::std::endl;
		::std::cerr << "Error rate:                      \t" << options.errorRate << ::std::endl;
		::std::cerr << "Minimal threshold:               \t" << options.threshold << ::std::endl;
		::std::cerr << "Shape:                           \t" << options.shape << ::std::endl;
		::std::cerr << "Repeat threshold:                \t" << options.repeatLength << ::std::endl;
		::std::cerr << "Overabundance threshold:         \t" << options.abundanceCut << ::std::endl;
		::std::cerr << "Taboo length:                    \t" << options.tabooLength << ::std::endl;
		::std::cerr << ::std::endl;
	}
	
	// circumvent numerical obstacles
	options.errorRate += 0.0000001;

	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Load fasta files
	SEQAN_PROTIMESTART(load_time);

	if (!loadReads(readSet, readNames, readFileName, options)) {
		::std::cerr << "Failed to load reads" << ::std::endl;
		return 1;
	}
	if (options._debugLevel >= 1) ::std::cerr << lengthSum(readSet) << " bps of " << length(readSet) << " reads loaded." << ::std::endl;
	options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);

	if (options._debugLevel >= 1)
		::std::cerr << "Loading files took               \t" << options.timeLoadFiles << " seconds" << ::std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Find matches using SWIFT
	if (!mapReads(matches, genomeFileName, genomeNames, readSet, options))
		return 2;

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Remove duplicates and output matches
	if (!options.spec.DONT_DUMP_RESULTS)
		dumpMatches(matches, genomeNames, genomeFileName, readSet, readNames, readFileName, errorPrbFileName, options);

	return 0;
}	


//////////////////////////////////////////////////////////////////////////////
// Dump an alignment
template <typename TFile, typename TSource, typename TSpec>
inline void
dumpAlignment(TFile & target, Align<TSource, TSpec> const & source)
{
	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;

	TRowsPosition row_count = length(rows(source));
	TPosition begin_ = beginPosition(cols(source));
	TPosition end_ = endPosition(cols(source));
	
	// Print sequences
	for(TRowsPosition i=0;i<row_count;++i) {
		if (i == 0)
			_streamWrite(target, "#Read:   ");
		else
			_streamWrite(target, "#Genome: ");
		TRow& row_ = row(source, i);
		typedef typename Iterator<typename Row<TAlign>::Type const>::Type TIter;
		TIter begin1_ = iter(row_, begin_);
		TIter end1_ = iter(row_, end_);
		for (; begin1_ != end1_; ++begin1_) {
			if (isGap(begin1_)) _streamPut(target, gapValue<char>());
			else _streamPut(target, *begin1_);
		}
		_streamPut(target, '\n');
	}
}


//////////////////////////////////////////////////////////////////////////////
// Output matches
template <
	typename TMatchSet,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames,
	typename TSpec
>
void dumpMatches(
	TMatchSet &matches,					// forward/reverse matches
	TGenomeNames const &genomeIDs,		// Read names (read from Fasta file, currently unused)
	::std::string const &genomeFName,	// genome name (e.g. "hs_ref_chr1.fa")
	TReads const &reads,				// Read sequences
	TReadNames const &readIDs,			// Read names (read from Fasta file, currently unused)
	::std::string const &readFName,		// read name (e.g. "reads.fa")
	::std::string const &errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatchSet>::Type		TMatches;
	typedef typename Value<TMatches>::Type		TMatch;
	typedef typename Value<TReads>::Type		TRead;
	typedef typename TMatch::TGPos				TGPos;

	unsigned maxReadLength = 0;
	for (unsigned i = 0; i < length(reads); ++i)
		if (maxReadLength < length(reads[i]))
			maxReadLength = length(reads[i]);

	String<int> posError;
	fill(posError, maxReadLength, 0);

	SEQAN_PROTIMESTART(dump_time);

	// load Genome sequences for alignment dumps
	TGenomeSet genomes;
	if (options.dumpAlignment || (!errorPrbFileName.empty() && options.hammingOnly))

		if (!loadGenomes(genomes, toCString(genomeFName))) {
			::std::cerr << "Failed to load genomes" << ::std::endl;
			options.dumpAlignment = false;
		}

	// how many 0's should be padded?
	int pzeros = 0;
	for (unsigned l = length(reads); l > 9; l = l / 10)
		++pzeros;

	int gzeros = 0;
	for (unsigned l = length(genomes); l > 9; l = l / 10)
		++gzeros;


	// remove the directory prefix of genomeFName and readFName
	size_t lastPos = genomeFName.find_last_of('/') + 1;
	if (lastPos == genomeFName.npos) lastPos = genomeFName.find_last_of('\\') + 1;
	if (lastPos == genomeFName.npos) lastPos = 0;
	::std::string genomeName = genomeFName.substr(lastPos);

	lastPos = readFName.find_last_of('/') + 1;
	if (lastPos == readFName.npos) lastPos = readFName.find_last_of('\\') + 1;
	if (lastPos == readFName.npos) lastPos = 0;
	::std::string readName = readFName.substr(lastPos);
	

	Align<String<Dna5>, ArrayGaps> align;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)

	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);

	::std::ofstream file;
	::std::ostringstream fileName;
	if (*options.output != 0)
		fileName << options.output;
	else
		fileName << readFName << ".result";

	file.open(fileName.str().c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
	if (!file.is_open()) {
		::std::cerr << "Failed to open output file" << ::std::endl;
		return;
	}

	String<int> hitCount;
	fill(hitCount, length(reads), 0);

	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal ends

	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoGEndPos<TMatch>());

	typename	TMatch::TGPos gBegin = -1;
	typename	TMatch::TGPos gEnd = -1;
	unsigned	gseqNo = -1;
	unsigned	readNo = -1;
	char		orientation = '-';

	typename Iterator<TMatches, Standard>::Type it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type itEnd = end(matches, Standard());

	for(; it != itEnd; ++it) 
	{
		if (gEnd == (*it).gEnd && orientation == (*it).orientation &&
			gseqNo == (*it).gseqNo && readNo == (*it).rseqNo) 
		{
			(*it).gseqNo = (unsigned)-1;
			(*it).rseqNo = (unsigned)-1;
			(*it).orientation = '-';
			continue;
		}
		readNo = (*it).rseqNo;
		gseqNo = (*it).gseqNo;
		gEnd = (*it).gEnd;
		orientation = (*it).orientation;
	}

	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal begins

	switch (options.sortOrder) {
		case 0:
			::std::sort(
				begin(matches, Standard()),
				end(matches, Standard()), 
				LessRNoGPos<TMatch>());
			break;

		case 1:
			::std::sort(
				begin(matches, Standard()),
				end(matches, Standard()), 
				LessGPosRNo<TMatch>());
			break;
	}

	gBegin = -1;
	gEnd = -1;
	gseqNo = -1;
	readNo = -1;
	orientation = '-';

	it = begin(matches, Standard());
	itEnd = end(matches, Standard());

	for(; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (gBegin == (*it).gBegin && readNo == (*it).rseqNo &&
			gseqNo == (*it).gseqNo && orientation == (*it).orientation) 
		{
			(*it).gseqNo = (unsigned)-1;
			(*it).rseqNo = (unsigned)-1;
			(*it).orientation = '-';
			continue;
		}
		readNo = (*it).rseqNo;
		gseqNo = (*it).gseqNo;
		gBegin = (*it).gBegin;
		orientation = (*it).orientation;
		++hitCount[readNo];
	}

	it = begin(matches, Standard());
	itEnd = end(matches, Standard());

	switch (options.outputFormat) 
	{
		case 0:	// Razer Format
			for(; it != itEnd; ++it) 
			{
				if ((*it).orientation == '-') continue;

				readNo = (*it).rseqNo;
#ifdef RAZERS_MAXHITS				
				if (hitCount[readNo] > options.maxHits) 
				{
					(*it).orientation = '-';
					continue;
				}
#endif
				gseqNo = (*it).gseqNo;
				gBegin = (*it).gBegin;
				gEnd = (*it).gEnd;
				orientation = (*it).orientation;

				unsigned	readLen = length(reads[readNo]);
				double		percId;
		
				percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);

				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << readIDs[readNo];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << readName << '#' << ::std::setw(pzeros) << readNo + 1;
						break;

					// 2..filename is the read sequence itself
					case 2:
						file << reads[readNo];
				}

				file << ',' << options.positionFormat << ',' << readLen << ',' << orientation << ',';

				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[gseqNo];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << genomeName << '#' << ::std::setw(gzeros) << gseqNo + 1;
				}

				file << ',' << (gBegin + options.positionFormat) << ',' << gEnd << ',' << ::std::setprecision(5) << percId << ::std::endl;

				if (options.dumpAlignment) {
					assignSource(row(align, 0), reads[readNo]);
					assignSource(row(align, 1), infix(genomes[gseqNo], gBegin, gEnd));
					if (orientation == 'R')
						reverseComplementInPlace(source(row(align, 1)));
					globalAlignment(align, scoreType);
					dumpAlignment(file, align);
				}
			}
			break;


		case 1:	// Enhanced Fasta Format
			for(; it != itEnd; ++it) 
			{
				if ((*it).orientation == '-') continue;

				readNo = (*it).rseqNo;
#ifdef RAZERS_MAXHITS				
				if (hitCount[readNo] > options.maxHits)
				{
					(*it).orientation = '-';
					continue;
				}
#endif
				gseqNo = (*it).gseqNo;
				gBegin = (*it).gBegin;
				gEnd = (*it).gEnd;
				orientation = (*it).orientation;

				unsigned	readLen = length(reads[readNo]);
				double		percId;

				::std::string fastaID;
				assign(fastaID, readIDs[readNo]);

				int id = readNo;
				int fragId = readNo;

				size_t left = fastaID.find_first_of('[');
				size_t right = fastaID.find_last_of(']');
				if (left != fastaID.npos && right != fastaID.npos && left < right) 
				{
					fastaID.erase(right);
					fastaID.erase(0, left + 1);
					replace(fastaID.begin(), fastaID.end(), ',', ' ');
					size_t pos = fastaID.find("id=");
					if (pos != fastaID.npos) {
						::std::istringstream iss(fastaID.substr(pos + 3));
						iss >> id;
					}
					pos = fastaID.find("fragId=");
					if (pos != fastaID.npos) {
						::std::istringstream iss(fastaID.substr(pos + 7));
						iss >> fragId;
					}
				}

				percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);

				if (orientation == 'F')
					// forward strand
					file << '>' << (gBegin + options.positionFormat) << ',' << gEnd;
				else
					// reverse strand (switch begin and end)
					file << '>' << gEnd << ',' << (gBegin + options.positionFormat);
				file << "[id=" << id << ",fragId=" << fragId;
				file << ",errors=" << (*it).editDist << ",percId=" << ::std::setprecision(5) << percId;
				file << ",ambiguity=" << hitCount[readNo] << ']' << ::std::endl;

				file << reads[readNo] << ::std::endl;
			}
			break;
	}

	file.close();

	// get empirical error distribution
	if (!errorPrbFileName.empty() && options.hammingOnly)
	{
		it = begin(matches, Standard());
		itEnd = end(matches, Standard());

		unsigned unique = 0;
		for (; it != itEnd; ++it) 
		{
			if ((*it).orientation == '-') continue;

			Dna5String const &read = reads[(*it).rseqNo];
			Dna5String genome = infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd);
			if ((*it).orientation == 'R')
				reverseComplementInPlace(genome);

			for (unsigned i = 0; i < length(read); ++i)
				if (genome[i] != read[i])
				{
					if (!(genome[i] == 'N' && read[i] != 'N'))
						++posError[i];
				}
			++unique;
		}
		file.open(errorPrbFileName.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
		if (!file.is_open()) {
			::std::cerr << "Failed to open error distribution file" << ::std::endl;
			return;
		}
		for (unsigned i = 0; i < length(posError); ++i)
			file << (double)posError[i] / unique << ::std::endl;
		file.close();
	}


	options.timeDumpResults = SEQAN_PROTIMEDIFF(dump_time);

	if (options._debugLevel >= 1)
		::std::cerr << "Dumping results took             \t" << options.timeDumpResults << " seconds" << ::std::endl;
}

}

#endif
