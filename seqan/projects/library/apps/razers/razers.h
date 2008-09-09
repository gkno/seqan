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
			errorRate = 0.2;
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

	typedef Dna5String			TGenome;
	typedef StringSet<TGenome>	TGenomeSet;
	typedef Dna5String			TRead;
	typedef StringSet<TRead>	TReadSet;


	template <typename TShape>
	struct SAValue< Index<TReadSet, TShape> > {
		typedef Pair<
			unsigned,		// many
			unsigned,		// short reads
			Compressed
		> Type;
	};

	template <typename TShape>
	struct Cargo< Index<TReadSet, TShape> > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};


	typedef ReadMatch<Difference<TGenome>::Type>	TMatch;			// a single match
	typedef String<TMatch/*, MMap<>*/ >				TMatches;		// array of matches

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
template <typename TReadSet, typename TNameSet>
bool loadFasta(TReadSet &reads, TNameSet &fastaIDs, char const *fileName)
{
	// count sequences
	unsigned seqCount = 0;

	::std::ifstream file;
	file.open(fileName, ::std::ios_base::in | ::std::ios_base::binary);
	if (!file.is_open()) return false;
	while (!_streamEOF(file)) {
		goNext(file, Fasta());
		++seqCount;
	}

	// import sequences
	file.clear();
	file.seekg(0, ::std::ios_base::beg);
	resize(fastaIDs, seqCount);
	resize(reads, seqCount);
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 
	{
		readID(file, fastaIDs[i], Fasta());		// read Fasta id
		read(file, reads[i], Fasta());			// read Read sequence
	}
	file.close();
	return (seqCount > 0);
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches (hamming only)
template <
	typename TMatches, 
	typename TGenomeSet, 
	typename TReadIndex, 
	typename THitCount, 
	typename TSpec >
void findReads(
	TMatches &matches,			// resulting matches
	TGenomeSet &genomes,		// Genome
	TReadIndex &readIndex,
	char orientation,			// q-gram index of reads
#ifdef RAZERS_MAXHITS
	THitCount &hitCount,		// maximum number of hits for each read
#else
	THitCount &,
#endif
	RazerSOptions<TSpec> &options,
	Swift<SwiftSemiGlobalHamming> )
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Value<TGenomeSet>::Type				TGenome;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Infix<TGenome>::Type					TGenomeInfix;
	typedef typename Value<TReadSet>::Type					TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;
	typedef typename Value<TMatches>::Type					TMatch;

	// FILTRATION
	
	typedef Finder<TGenome, Swift<SwiftSemiGlobalHamming> >		TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<SwiftSemiGlobalHamming> >	TSwiftPattern;

	TSwiftPattern swiftPattern(readIndex);

	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;
	swiftPattern.params.minLog2Delta = 0;

	__int64 TP = 0;
	__int64 FP = 0;
	SEQAN_PROTIMESTART(find_time);

	// iterate all genomic sequences
	for(unsigned hstkSeqNo = 0; hstkSeqNo < length(genomes); ++hstkSeqNo) 
	{
		if (options._debugLevel >= 1)
			::std::cerr << ::std::endl << "Process genome seq #" << hstkSeqNo;
		TGenome &genome = genomes[hstkSeqNo];
		TSwiftFinder swiftFinder(genome, options.repeatLength, 1);

		// iterate all verification regions returned by SWIFT
		while (find(swiftFinder, swiftPattern, options.errorRate, options._debugLevel)) 
		{
			unsigned ndlSeqNo = (*swiftFinder.curHit).ndlSeqNo;
#ifdef RAZERS_MAXHITS
			if (hitCount[ndlSeqNo] == -1) continue;
#endif			
			unsigned ndlLength = sequenceLength(ndlSeqNo, readIndex);

			TGenomeInfix inf(range(swiftFinder, genome));
#ifdef RAZERS_DEBUG
			cout<<"Verify: "<<::std::endl;
			cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
			cout<<"Read:   "<<host(myersPattern)<<::std::endl;
#endif
			if (options.spec.DONT_VERIFY) 
			{
				++FP;
				continue;
			}

			// verify
			TMatch m = {
				hstkSeqNo,
				ndlSeqNo,
				beginPosition(swiftFinder),
				endPosition(swiftFinder),
				0,
				orientation
			};

			int errors = ndlLength - (m.gEnd - m.gBegin);
			int maxErrors = (int)(ndlLength * options.errorRate);
			TReadIterator rit = begin(indexText(readIndex)[ndlSeqNo], Standard()) + errors;
			TGenomeIterator git = begin(genome, Standard()) + m.gBegin;
			TGenomeIterator gitEnd;

			if ((TSize)m.gEnd > length(genome))
			{
				gitEnd = end(genome, Standard());
				errors += m.gEnd - length(genome);
				m.gEnd = length(genome);
			} else
				gitEnd = begin(genome, Standard()) + m.gEnd;

			if (errors <= maxErrors)
			{
				if (options.matchN)
				{
					for(; git != gitEnd; ++git, ++rit)
						if (*git != *rit && *git != 'N' && *rit != 'N')
							if (++errors > maxErrors)
								break;
				} else {
					for(; git != gitEnd; ++git, ++rit)
						if (*git != *rit)
							if (++errors > maxErrors)
								break;
				}
			}

			if (errors <= maxErrors)
			{
				m.editDist = errors;
				// transform coordinates to the forward strand
				if (orientation == 'R') 
				{
					TSize gLength = length(genome);
					TSize temp = m.gBegin;
					m.gBegin = gLength - m.gEnd;
					m.gEnd = gLength - temp;
				}
				if (!options.spec.DONT_DUMP_RESULTS)
					appendValue(matches, m);
#ifdef RAZERS_MAXHITS				
				--hitCount[ndlSeqNo];
#endif
				++TP;
	/*			::std::cerr << "\"" << inf << "\"  ";
				::std::cerr << hstkPos << " + ";
				::std::cerr << ::std::endl;
	*/		} else {
				++FP;
	/*			::std::cerr << "\"" << inf << "\"   \"" << range(swiftPattern) << "\"  ";
				::std::cerr << 'R' << ndlSeqNo << ", G";
				::std::cerr << beginPosition(swiftFinder) << ':' << endPosition(swiftFinder) << ::std::endl;
	*/		}
		}
	}

	options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;

	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Swift FP: " << FP << ::std::endl;
		::std::cerr << "Swift TP: " << TP << ::std::endl;
	}
	options.FP += FP;
	options.TP += TP;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches
template <
	typename TMatches, 
	typename TGenomeSet, 
	typename TReadIndex, 
	typename THitCount, 
	typename TSpec,
	typename TSwiftSpec >
void findReads(
	TMatches &matches,			// resulting matches
	TGenomeSet &genomes,		// Genome
	TReadIndex &readIndex,
	char orientation,			// q-gram index of reads
#ifdef RAZERS_MAXHITS
	THitCount &hitCount,		// maximum number of hits for each read
#else
	THitCount &,
#endif
	RazerSOptions<TSpec> &options,
	Swift<TSwiftSpec> )
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Value<TGenomeSet>::Type				TGenome;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Infix<TGenome>::Type					TGenomeInfix;
	typedef typename Value<TReadSet>::Type					TRead;
	typedef typename Value<TMatches>::Type					TMatch;

	// FILTRATION
	
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// VERIFICATION

	// find read match end
	typedef Finder<TGenomeInfix>							TMyersFinder;
#ifdef RAZERS_HAMMINGVERIFY
	typedef Score<int>										TScore;
	typedef Pattern<TRead, DPSearch<TScore> >				TMyersPattern;
#else
	typedef Pattern<TRead, MyersUkkonen>					TMyersPattern;
#endif

	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>	TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>			TReadRev;
	typedef Finder<TGenomeInfixRev>						TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>		TMyersPatternRev;

	TSwiftPattern swiftPattern(readIndex);
	String<TMyersPattern> forwardPatterns;

	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;

#ifdef RAZERS_HAMMINGVERIFY
	TScore	scoreType = Score<int>(0, -1, -1001, -1000); // levenshtein-score (match, mismatch, gapOpen, gapExtend)
#endif

	__int64 TP = 0;
	__int64 FP = 0;
	SEQAN_PROTIMESTART(find_time);

	// init forward verifiers
	unsigned readCount = countSequences(readIndex);
	resize(forwardPatterns, readCount);
	for(unsigned i = 0; i < readCount; ++i)
	{
#ifdef RAZERS_HAMMINGVERIFY
		setHost(forwardPatterns[i], indexText(readIndex)[i]);
		setScoringScheme(forwardPatterns[i], scoreType);
#else
		setHost(forwardPatterns[i], indexText(readIndex)[i]);
		if (options.matchN)
		{
			_patternMatchNOfPattern(forwardPatterns[i]);
			_patternMatchNOfFinder(forwardPatterns[i]);
		}
#endif
	}

	// iterate all genomic sequences
	for(unsigned hstkSeqNo = 0; hstkSeqNo < length(genomes); ++hstkSeqNo) 
	{
		if (options._debugLevel >= 1)
			::std::cerr << ::std::endl << "Process genome seq #" << hstkSeqNo;
		TGenome &genome = genomes[hstkSeqNo];
		TSwiftFinder swiftFinder(genome, options.repeatLength, 1);

		// iterate all verification regions returned by SWIFT
		while (find(swiftFinder, swiftPattern, options.errorRate, options._debugLevel)) 
		{
			unsigned ndlSeqNo = (*swiftFinder.curHit).ndlSeqNo;
#ifdef RAZERS_MAXHITS
			if (hitCount[ndlSeqNo] == -1) continue;
#endif			
			unsigned ndlLength = sequenceLength(ndlSeqNo, readIndex);

			TGenomeInfix inf(range(swiftFinder, genome));
			TMyersFinder myersFinder(inf);
			TMyersPattern &myersPattern = forwardPatterns[ndlSeqNo];
#ifdef RAZERS_DEBUG
			cout<<"Verify: "<<::std::endl;
			cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
			cout<<"Read:   "<<host(myersPattern)<<::std::endl;
#endif
			// find end of best semi-global alignment
			int maxScore = InfimumValue<int>::VALUE;
			int minScore = -(int)(ndlLength * options.errorRate);
			TMyersFinder maxPos;
			while (!options.spec.DONT_VERIFY && find(myersFinder, myersPattern, minScore))
				if (maxScore < getScore(myersPattern)) 
				{
					maxScore = getScore(myersPattern);
					maxPos = myersFinder;
				}
			
			if (maxScore >= minScore) 
			{
				TMatch m = {
					hstkSeqNo,
					ndlSeqNo,
					beginPosition(swiftFinder),
					beginPosition(swiftFinder) + position(maxPos) + 1,
					-maxScore,
					orientation
				};
				if (m.gBegin < 0) m.gBegin = 0;
				
#ifdef RAZERS_HAMMINGVERIFY
				m.gBegin = m.gEnd - ndlLength;
#else
				TGenomeInfixRev		infRev(infix(genome, m.gBegin, m.gEnd));
				TReadRev			readRev(indexText(readIndex)[ndlSeqNo]);
				TMyersFinderRev		myersFinderRev(infRev);
				TMyersPatternRev	myersPatternRev(readRev);

				if (options.matchN)
				{
					_patternMatchNOfPattern(myersPatternRev);
					_patternMatchNOfFinder(myersPatternRev);
				}
				// find beginning of best semi-global alignment
				if (find(myersFinderRev, myersPatternRev, maxScore))
					m.gBegin = m.gEnd - (position(myersFinderRev) + 1);
				else {
					// this case should never occur
					::std::cerr << "1GENOME: " << host(myersFinder) << ::std::endl;
					::std::cerr << "1READ:   " << indexText(readIndex)[ndlSeqNo] << ::std::endl;
					::std::cerr << "2GENOME: " << infix(genome, m.gBegin, m.gEnd) << '\t' << m.gBegin << ',' << m.gEnd << ::std::endl;
					::std::cerr << "2READ:   " << indexText(readIndex)[ndlSeqNo] << ::std::endl;
					::std::cerr << "3GENOME: " << infRev << ::std::endl;
					::std::cerr << "3READ:   " << readRev << ::std::endl;
					::std::cerr << "HUH?" << ::std::endl;
				}
#endif

				// transform coordinates to the forward strand
				if (orientation == 'R') 
				{
					TSize gLength = length(genome);
					TSize temp = m.gBegin;
					m.gBegin = gLength - m.gEnd;
					m.gEnd = gLength - temp;
				}
				if (!options.spec.DONT_DUMP_RESULTS)
					appendValue(matches, m);
#ifdef RAZERS_MAXHITS				
				--hitCount[ndlSeqNo];
#endif
				++TP;
	/*			::std::cerr << "\"" << range(swiftFinder, genomeInf) << "\"  ";
				::std::cerr << hstkPos << " + ";
				::std::cerr << ::std::endl;
	*/		} else {
				++FP;
	/*			::std::cerr << "\"" << range(swiftFinder, genomeInf) << "\"   \"" << range(swiftPattern) << "\"  ";
				::std::cerr << ndlSeqNo << " : ";
				::std::cerr << hstkPos << " + ";
				::std::cerr << bucketWidth << "  " << TP << ::std::endl;
	*/		}
		}
	}

	options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;

	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Swift FP: " << FP << ::std::endl;
		::std::cerr << "Swift TP: " << TP << ::std::endl;
	}
	options.FP += FP;
	options.TP += TP;
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
	typename TGenome,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames,
	typename TSpec
>
void dumpMatches(
	TMatchSet &matches,					// forward/reverse matches
	TGenome const &genomes,				// Genome sequences
	TGenomeNames const &genomeIDs,		// Read names (read from Fasta file, currently unused)
	TReads const &reads,				// Read sequences
	TReadNames const &readIDs,			// Read names (read from Fasta file, currently unused)
	::std::string const &genomeFName,	// genome name (e.g. "hs_ref_chr1.fa")
	::std::string const &readFName,		// read name (e.g. "reads.fa")
	RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatchSet>::Type		TMatches;
	typedef typename Value<TMatches>::Type		TMatch;
	typedef typename Value<TReads>::Type		TRead;
	typedef typename TMatch::TGPos				TGPos;

	SEQAN_PROTIMESTART(dump_time);

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
	

	Align<TRead, ArrayGaps> align;
#ifdef RAZERS_HAMMINGVERIFY
	Score<int> scoreType = Score<int>(0, -1, -1001, -1000);		// levenshtein-score (match, mismatch, gapOpen, gapExtend)
#else
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
#endif
	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);

	bool multipleGenomes = countSequences(genomes) > 1;
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
				if (hitCount[readNo] > options.maxHits) continue;
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
						if (multipleGenomes) {
							file.fill('0');
							file << genomeName << '#' << ::std::setw(gzeros) << gseqNo + 1;
						} else
							file << genomeName;
				}

				file << ',' << gBegin + options.positionFormat << ',' << gEnd << ',' << ::std::setprecision(5) << percId << ::std::endl;

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
				if (hitCount[readNo] > options.maxHits) continue;
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
					file << '>' << gBegin + options.positionFormat << ',' << gEnd;
				else
					// reverse strand (switch begin and end)
					file << '>' << gEnd << ',' << gBegin + options.positionFormat;
				file << "[id=" << id << ",fragId=" << fragId;
				file << ",errors=" << (*it).editDist << ",percId=" << ::std::setprecision(5) << percId;
				file << ",ambiguity=" << hitCount[readNo] << ']' << ::std::endl;

				file << reads[readNo] << ::std::endl;
			}
			break;
	}

	file.close();

	options.timeDumpResults = SEQAN_PROTIMEDIFF(dump_time);

	if (options._debugLevel >= 1)
		::std::cerr << "Dumping results took             \t" << options.timeDumpResults << " seconds" << ::std::endl;
}

//////////////////////////////////////////////////////////////////////////////
// Main Part

template <
	typename TMatches, 
	typename TGenomeSet, 
	typename TReadSet, 
	typename TSpec, 
	typename TShape >
void mapReads(
	TMatches &matches, 
	TGenomeSet &genomeSet, 
	TReadSet const &readSet, 
	RazerSOptions<TSpec> &options,
	TShape const &shape)
{
	Index<TReadSet, Index_QGram<TShape> > swiftIndex(readSet, shape);
	String<int> hitCount;

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeLoadFiles = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

#ifdef RAZERS_MAXHITS	
	// fill max number of hits per read
	// at most twice as much potential matches are allowed (overlapping parallelograms)
	fill(hitCount, length(readSet), 2*options.maxHits);	
#endif

	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;

	if (options.forward)
	{
		if (options._debugLevel >= 1)
			::std::cerr << ::std::endl << "___FORWARD_STRAND______";
		if (options.hammingOnly)
			findReads(matches, genomeSet, swiftIndex, 'F', hitCount, options, Swift<SwiftSemiGlobalHamming>());
		else
			findReads(matches, genomeSet, swiftIndex, 'F', hitCount, options, Swift<SwiftSemiGlobal>());
	}

	if (options.reverse) 
	{
		if (options._debugLevel >= 1)
			::std::cerr << ::std::endl << "___BACKWARD_STRAND_____";
		reverseComplementInPlace(genomeSet);			// build reverse-compl of genome
		if (options.hammingOnly)
			findReads(matches, genomeSet, swiftIndex, 'R', hitCount, options, Swift<SwiftSemiGlobalHamming>());
		else
			findReads(matches, genomeSet, swiftIndex, 'R', hitCount, options, Swift<SwiftSemiGlobal>());
		reverseComplementInPlace(genomeSet);			// restore original genome seqs
	}
}

template <typename TMatches, typename TGenomeSet, typename TReadSet, typename TSpec>
bool mapReads(
	TMatches &matches,
	TGenomeSet &genomeSet, 
	TReadSet const &readSet, 
	RazerSOptions<TSpec> &options)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GappedShape>		gapped;

	// select best-fitting shape

	if (stringToShape(ungapped, options.shape)) {
		mapReads(matches, genomeSet, readSet, options, ungapped);
		return true;
	}
	
	if (stringToShape(onegapped, options.shape)) {
		mapReads(matches, genomeSet, readSet, options, onegapped);
		return true;
	}

	if (stringToShape(gapped, options.shape)) {
		mapReads(matches, genomeSet, readSet, options, gapped);
		return true;
	}

	return false;
}

template <typename TSpec>
int mapReads(
	const char *genomeFileName,
	const char *readFileName,
	RazerSOptions<TSpec> &options)
{
	TGenomeSet				genomeSet;
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
	if (!loadFasta(genomeSet, genomeNames, genomeFileName)) {
		::std::cerr << "Failed to load genomes" << ::std::endl;
		return 1;
	}
	if (options._debugLevel >= 1) ::std::cerr << lengthSum(genomeSet) << " bps of " << length(genomeSet) << " genomes loaded." << ::std::endl;

	if (!loadFasta(readSet, readNames, readFileName)) {
		::std::cerr << "Failed to load reads" << ::std::endl;
		return 1;
	}
	if (options._debugLevel >= 1) ::std::cerr << lengthSum(readSet) << " bps of " << length(readSet) << " reads loaded." << ::std::endl;
	options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);

	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Find matches using SWIFT
	if (!mapReads(matches, genomeSet, readSet, options))
		return 2;

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Remove duplicates and output matches
	if (!options.spec.DONT_DUMP_RESULTS)
		dumpMatches(matches, genomeSet, genomeNames, readSet, readNames, genomeFileName, readFileName, options);

	return 0;
}	

}

#endif
