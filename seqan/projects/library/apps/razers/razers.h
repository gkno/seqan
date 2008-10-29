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

#include <seqan/find.h>
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
		unsigned	maxHits;			// output at most maxHits many matches
		unsigned	distanceRange;		// output only the best, second best, ..., distanceRange best matches
		bool		purgeAmbiguous;		// true..remove reads with more than maxHits best matches, false..keep them
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

	// misc
		unsigned	compactThresh;		// compact match array if larger than compactThresh
#ifdef RAZERS_MASK_READS
		String<unsigned long> readMask;	// bit-string of bool (1..verify read, 0..ignore read)
		enum { WORD_SIZE = BitsPerValue<unsigned long>::VALUE };
#endif

		RazerSOptions() 
		{
			forward = true;
			reverse = true;
			errorRate = 0.08;
			maxHits = 100;
			distanceRange = 0;	// disabled
			purgeAmbiguous = false;
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

			compactThresh = 1024;
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
		unsigned short	editDist;		// Levenshtein distance
		char			orientation;	// 'F'..forward strand, 'R'..reverse comp. strand
	};
	
	enum RAZERS_ERROR {
		RAZERS_READS_FAILED = 1,
		RAZERS_GENOME_FAILED = 2,
		RAZERS_INVALID_SHAPE = 3
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

	// ... to sort matches and remove duplicates with equal gEnd
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

	template <typename TReadMatch>
	struct LessErrors : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};
	
//////////////////////////////////////////////////////////////////////////////
// Remove duplicate matches and leave at most maxHits many distanceRange
// best matches per read
template < typename TMatches >
void maskDuplicates(TMatches &matches)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
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

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());

	for (; it != itEnd; ++it) 
	{
		if (gEnd == (*it).gEnd && orientation == (*it).orientation &&
			gseqNo == (*it).gseqNo && readNo == (*it).rseqNo) 
		{
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

	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoGPos<TMatch>());

	orientation = '-';

	it = begin(matches, Standard());
	itEnd = end(matches, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (gBegin == (*it).gBegin && readNo == (*it).rseqNo &&
			gseqNo == (*it).gseqNo && orientation == (*it).orientation) 
		{
			(*it).orientation = '-';
			continue;
		}
		readNo = (*it).rseqNo;
		gseqNo = (*it).gseqNo;
		gBegin = (*it).gBegin;
		orientation = (*it).orientation;
	}
	
	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessErrors<TMatch>());
}

//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template < typename TMatches, typename TCounts >
void countMatches(TMatches &matches, TCounts &cnt)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	typedef typename Value<TCounts>::Type					TRow;
	typedef typename Value<TRow>::Type						TValue;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	
	unsigned readNo = -1;
	short editDist = -1;
	__int64 count = 0;
	__int64 maxVal = SupremumValue<TValue>::VALUE;

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo && editDist == (*it).editDist)
			++count;
		else
		{
			if (readNo != (unsigned)-1 && (unsigned)editDist < length(cnt))
				cnt[editDist][readNo] = (maxVal < count)? maxVal : count;
			readNo = (*it).rseqNo;
			editDist = (*it).editDist;
			count = 1;
		}
	}
	if (readNo != (unsigned)-1 && (unsigned)editDist < length(cnt))
		cnt[editDist][readNo] = count;
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TMatches, typename TSpec >
void compactMatches(TMatches &matches, RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int editDistCutOff = SupremumValue<int>::VALUE;

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo)
		{ 
			if ((*it).editDist >= editDistCutOff) continue;
			if (++hitCount >= hitCountCutOff)
			{
				if (hitCount == hitCountCutOff)
				{
#ifdef RAZERS_MASK_READS
					if (options.purgeAmbiguous)
					{
						dit = ditBeg;
						if (options._debugLevel >= 2)
							::std::cerr << "(read #" << readNo << " disabled)";
						options.readMask[readNo / options.WORD_SIZE] &= ~(1ul << (readNo % options.WORD_SIZE));
					} else
						if ((*it).editDist == 0)
						{
							if (options._debugLevel >= 2)
								::std::cerr << "(read #" << readNo << " disabled)";
							options.readMask[readNo / options.WORD_SIZE] &= ~(1ul << (readNo % options.WORD_SIZE));
						}
#endif
				}
				continue;
			}
		} else
		{
			readNo = (*it).rseqNo;
			hitCount = 0;
			if (options.distanceRange > 0)
				editDistCutOff = (*it).editDist + options.distanceRange;
			ditBeg = dit;
		}
		*dit = *it;
		++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
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

	unsigned maxErrors = (unsigned)(ndlLength * options.errorRate);
	unsigned minErrors = maxErrors + 1;

	for (; git < gitEnd; ++git)
	{
		unsigned errors = 0;
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
		if (maxScore <= getScore(myersPattern)) 
		{
			maxScore = getScore(myersPattern);
			maxPos = myersFinder;
		}
	
	if (maxScore < minScore) return false;
	m.editDist	= (unsigned)-maxScore;
	setEndPosition(inf, m.gEnd = (beginPosition(inf) + position(maxPos) + 1));

	// limit the beginning to needle length plus errors (== -maxScore)
	if (length(inf) > ndlLength - maxScore)
		setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
	
	// find beginning of best semi-global alignment
	TGenomeInfixRev		infRev(inf);
	TReadRev			readRev(indexText(readIndex)[rseqNo]);
	TMyersFinderRev		myersFinderRev(infRev);
	TMyersPatternRev	myersPatternRev(readRev);

/*	cout<<"Verify: "<<::std::endl;
	cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	cout<<"Read:   "<<indexText(readIndex)[rseqNo]<<::std::endl;

	cout<<"Genome: "<<infRev << ::std::endl;
	cout<<"Read:   "<<readRev<<::std::endl;
*/
	_patternMatchNOfPattern(myersPatternRev, options.matchN);
	_patternMatchNOfFinder(myersPatternRev, options.matchN);
	while (find(myersFinderRev, myersPatternRev, maxScore))
		m.gBegin = m.gEnd - (position(myersFinderRev) + 1);

	return true;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename TSpec >
void findReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPattern,
	TVerifier &forwardPatterns,
	char orientation,				// q-gram index of reads
	RazerSOptions<TSpec> &options)
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Value<TMatches>::Type					TMatch;

	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinder;

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

	TReadIndex &readIndex = host(swiftPattern);
	TSwiftFinder swiftFinder(genome, options.repeatLength, 1);

	TMatch m = { 0, 0, 0, 0, 0, 0 };	// to supress uninitialized warnings
	
	// iterate all verification regions returned by SWIFT
	while (find(swiftFinder, swiftPattern, options.errorRate, options._debugLevel)) 
	{
		unsigned rseqNo = (*swiftFinder.curHit).ndlSeqNo;
		if (!options.spec.DONT_VERIFY && 
#ifdef RAZERS_MASK_READS
			((options.readMask[rseqNo / options.WORD_SIZE] & (1ul << (rseqNo % options.WORD_SIZE))) != 0) &&
#endif
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
			{
				appendValue(matches, m);
				if (length(matches) > options.compactThresh)
				{
					typename Size<TMatches>::Type oldSize = length(matches);
					maskDuplicates(matches);
					compactMatches(matches, options);
					options.compactThresh += (options.compactThresh >> 1);
					if (options._debugLevel >= 2)
						::std::cerr << '(' << oldSize - length(matches) << " matches removed)";
				}
			}

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
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TMatches, 
	typename TReadSet, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapReads(
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
	if (!file.is_open())
		return RAZERS_GENOME_FAILED;

	// configure q-gram index
	TIndex swiftIndex(readSet, shape);
	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;

	// configure Swift
	TSwiftPattern swiftPattern(swiftIndex);
	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;

	// init edit distance verifiers
	unsigned readCount = countSequences(swiftIndex);
	String<TMyersPattern> forwardPatterns;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
	{
		resize(forwardPatterns, readCount);
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}

#ifdef RAZERS_MASK_READS
	// init read mask
	clear(options.readMask);
	fill(options.readMask, (readCount + options.WORD_SIZE - 1) / options.WORD_SIZE, (unsigned long)-1);
#endif

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

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
			findReads(matches, genome, gseqNo, swiftPattern, forwardPatterns, 'F', options);

		if (options.reverse)
		{
			reverseComplementInPlace(genome);
			findReads(matches, genome, gseqNo, swiftPattern, forwardPatterns, 'R', options);
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
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (given as StringSet)
template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type				TRead;
	typedef Index<TReadSet, Index_QGram<TShape> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

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
	options.compMask[4] = (options.matchN)? 15: 0;
	if (!options.hammingOnly)
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

	CharString	id;
	
	// iterate over genome sequences
	SEQAN_PROTIMESTART(find_time);
	for(unsigned gseqNo = 0; gseqNo < length(genomeSet); ++gseqNo)
	{
		if (options.forward)
			findReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, 'F', options);

		if (options.reverse)
		{
			reverseComplementInPlace(genomeSet[gseqNo]);
			findReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, 'R', options);
			reverseComplementInPlace(genomeSet[gseqNo]);
		}

	}
	options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);

	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Swift FP: " << options.FP << ::std::endl;
		::std::cerr << "Swift TP: " << options.TP << ::std::endl;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Wrapper for different template specializations
template <typename TMatches, typename TReadSet, typename TSpec>
int mapReads(
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

	return RAZERS_INVALID_SHAPE;
}

template <typename TMatches, typename TGenomeSet, typename TReadSet, typename TSpec>
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
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
			return mapReads(matches, genomeSet, readSet, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

}

#endif
