 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
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

#ifdef RAZERS_PARALLEL
#include "tbb/spin_mutex.h"
#endif

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
		int			trimLength;			// if >0, cut reads to #trimLength characters
		
	// output format options
		unsigned	outputFormat;		// 0..Razer format
										// 1..enhanced Fasta
										// 2..ELAND format
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
		const char	*runID;				// runID needed for gff output

	// filtration parameters
		::std::string shape;			// shape (e.g. 11111111111)
		int			threshold;			// threshold
		int			tabooLength;		// taboo length
		int			repeatLength;		// repeat length threshold
		double		abundanceCut;		// abundance threshold

	// mate-pair parameters
		int			libraryLength;		// offset between two mates
		int			libraryError;		// offset tolerance
		unsigned	nextMatePairId;		// use this id for the next mate-pair

	// verification parameters
		bool		matchN;				// false..N is always a mismatch, true..N matches with all
		unsigned char compMask[5];

	// statistics
		__int64		FP;					// false positives (threshold reached, no match)
		__int64		TP;					// true positives (threshold reached, match)
		double		timeLoadFiles;		// time for loading input files
		double		timeMapReads;		// time for mapping reads
		double		timeDumpResults;	// time for dumping the results
		
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		bool		maqMapping;
		int		maxMismatchQualSum;
		int		mutationRateQual;
		int		absMaxQualSumErrors;
		unsigned	artSeedLength;
		bool		noBelowIdentity;
#endif

		bool		lowMemory;		// set maximum shape weight to 13 to limit size of q-gram index
		bool		fastaIdQual;		// hidden option for special fasta+quality format we use


	// misc
		unsigned	compactThresh;		// compact match array if larger than compactThresh
#ifdef RAZERS_MASK_READS
		String<unsigned long> readMask;	// bit-string of bool (1..verify read, 0..ignore read)
		enum { WORD_SIZE = BitsPerValue<unsigned long>::VALUE };
#endif

	// multi-threading

#ifdef RAZERS_PARALLEL
		typedef ::tbb::spin_mutex	TMutex;

		TMutex		*patternMutex;
		TMutex		optionsMutex;
		TMutex		matchMutex;
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
			trimLength = 0;
			
			outputFormat = 0;
			dumpAlignment = false;
			genomeNaming = 0;
			readNaming = 0;
			sortOrder = 0;
			positionFormat = 0;
			runID = "s"; 	//

			matchN = false;

			shape = "11111111111";
			threshold = 1;
			tabooLength = 1;
			repeatLength = 1000;
			abundanceCut = 1;

			libraryLength = 2000;
			libraryError = 200;
			nextMatePairId = 1;
			
			for (unsigned i = 0; i < 4; ++i)
				compMask[i] = 1 << i;
			compMask[4] = 0;

			compactThresh = 1024;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			maqMapping = false;
			maxMismatchQualSum = 70;
			mutationRateQual = 30;
			artSeedLength = 28;	// the "artificial" seed length that is used for mapping quality assignment 
						// (28bp is maq default)
			absMaxQualSumErrors = 100; // maximum for sum of mism qualities in total readlength
			noBelowIdentity = false;
#endif
			lowMemory = false;		// set maximum shape weight to 13 to limit size of q-gram index
			fastaIdQual = false;

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
#ifdef RAZERS_MATEPAIRS
		unsigned		pairId;			// unique id for the two mate-pair matches (0 if unpaired)
		int				mateDelta:24;	// outer coordinate delta to the other mate 
		int				pairScore:8;	// combined score of both mates
#endif
		unsigned short	editDist;		// Levenshtein distance
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		short	 		mScore;
		short			seedEditDist;
#endif
		char			orientation;	// 'F'..forward strand, 'R'..reverse comp. strand
	};
	
	enum RAZERS_ERROR {
		RAZERS_READS_FAILED = -1,
		RAZERS_GENOME_FAILED = -2,
		RAZERS_INVALID_SHAPE = -3
	};

//////////////////////////////////////////////////////////////////////////////
// Definitions

	typedef Dna5String									TGenome;
	typedef StringSet<TGenome>							TGenomeSet;
//	typedef Dna5String									TRead;
	typedef String<Dna5Q>								TRead;
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
			BitCompressed<24, 8>	// max. 16M reads of length < 256
		> Type;
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
// Load multi-Fasta sequences with or w/o quality values
template <typename TReadSet, typename TNameSet, typename TRazerSOptions>
bool loadReads(
	TReadSet &reads, 
	TNameSet &fastaIDs, 
	const char *fileName, 
	TRazerSOptions &options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);

	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);

	unsigned seqCount = length(multiFasta);
#ifndef RAZERS_CONCATREADS
	resize(reads, seqCount, Exact());
#endif
	if (options.readNaming == 0)
		resize(fastaIDs, seqCount, Exact());
	
	String<Dna5Q> seq;
	CharString qual;
	
	unsigned kickoutcount = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			|| options.fastaIdQual
#endif
			)
			assignSeqId(fastaIDs[i], multiFasta[i], format);	// read Fasta id
		assignSeq(seq, multiFasta[i], format);					// read Read sequence
		assignQual(qual, multiFasta[i], format);				// read ascii quality values  
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		if(options.fastaIdQual)
		{
			qual = suffix(fastaIDs[i],length(fastaIDs[i])-length(seq));
			if(options.readNaming == 0)
				fastaIDs[i] = prefix(fastaIDs[i],length(fastaIDs[i])-length(seq));
			else clear(fastaIDs[i]);
		}
#endif
		if (countN)
		{
			int count = 0;
			int cutoffCount = (int)(options.errorRate * length(seq));
			for (unsigned j = 0; j < length(seq); ++j)
				if (getValue(seq, j) == 'N')
					if (++count > cutoffCount)
					{
						clear(seq);
						clear(qual);
						++kickoutcount;
						break;
					}
// low qual. reads are empty to output them and their id later as LQ reads
//			if (count > cutoffCount) continue;
		}

		// store dna and quality together
		for (unsigned j = 0; j < length(qual) && j < length(seq); ++j)
			assignQualityValue(seq[j], (int)(ordValue(qual[j]) - 33));
/*			if(ordValue(seq[j])>3)
				setValue(hybridSeq[j],(unsigned int)164); //N=164 such that &3 == 0
			else
				setValue(hybridSeq[j],(unsigned int)((ordValue(qual[j]) - 33) * 4 + ordValue(seq[j])));
*/
		if (options.trimLength > 0 && length(seq) > (unsigned)options.trimLength)
			resize(seq, options.trimLength);

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
// Read the first sequence of a multi-sequence file
// and return its length
inline int estimateReadLength(char const *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY))	// open the whole file
		return RAZERS_READS_FAILED;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);					// guess file format
	split(multiFasta, format);								// divide into single sequences

	if (length(multiFasta) == 0)
		return 0;

	Dna5String firstRead;
	assignSeq(firstRead, multiFasta[0], format);			// read the first sequence
	return length(firstRead);
}


/*
//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TGenomeSet>
bool loadGenomes(TGenomeSet &genomes, const char *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
	split(multiFasta, Fasta());

	unsigned seqCount = length(multiFasta);
	resize(genomes, seqCount, Exact());
	for(unsigned i = 0; i < seqCount; ++i)
		assignSeq(genomes[i], multiFasta[i], Fasta());		// read Genome sequence

	return (seqCount > 0);
}*/

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences from multiple files
template <typename TGenomeSet>
bool loadGenomes(TGenomeSet &genomes, StringSet<CharString> &fileNameList)
{
	unsigned gSeqNo = 0;
	unsigned filecount = 0;
	while(filecount < length(fileNameList))
	{
		MultiFasta multiFasta;
		if (!open(multiFasta.concat, toCString(fileNameList[filecount]), OPEN_RDONLY)) return false;
		split(multiFasta, Fasta());

		unsigned seqCount = length(multiFasta);
		if(length(genomes) < gSeqNo+seqCount) 
			resize(genomes,gSeqNo+seqCount);
		for(unsigned i = 0; i < seqCount; ++i)
			assignSeq(genomes[gSeqNo+i], multiFasta[i], Fasta());		// read Genome sequence
		gSeqNo += seqCount;
		++filecount;
	}
	resize(genomes,gSeqNo);
	return (gSeqNo > 0);
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
#ifdef RAZERS_MATEPAIRS
			return a.pairScore > b.pairScore;
#else
			return a.editDist < b.editDist;
#endif
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
#ifdef RAZERS_MATEPAIRS
			return a.pairScore > b.pairScore;
#else
			return a.editDist < b.editDist;
#endif
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
#ifdef RAZERS_MATEPAIRS
			return a.pairScore > b.pairScore;
#else
			return a.editDist < b.editDist;
#endif
		}
	};
	
#ifdef RAZERS_DIRECT_MAQ_MAPPING

	struct QualityBasedScoring{};

	template <typename TReadMatch>
	struct LessRNoMQ : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;
			
			// quality
			if (a.mScore < b.mScore) return true; // sum of quality values of mismatches (the smaller the better)
			if (a.mScore > b.mScore) return false;
			
			return (a.editDist < b.editDist); // seedEditDist?
			// genome position and orientation
	/*		if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;
	*/		
		}
	};
#endif

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
#ifdef RAZERS_MATEPAIRS
		if ((*it).pairId != 0) continue;
#endif
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
		if ((*it).orientation == '-'
#ifdef RAZERS_MATEPAIRS
			|| ((*it).pairId != 0)
#endif
			) continue;
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

template < typename TReadNo, typename TMaxErrors >
inline void 
setMaxErrors(Nothing &, TReadNo, TMaxErrors)
{
}

template < typename TSwift, typename TReadNo, typename TMaxErrors >
inline void 
setMaxErrors(TSwift &swift, TReadNo readNo, TMaxErrors maxErrors)
{
	int minT = _qgramLemma(swift, readNo, maxErrors);
	if (minT > 1)
	{
//		::std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
		setMinThreshold(swift, readNo, (unsigned)minT);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TMatches, typename TCounts, typename TSpec, typename TSwift >
void compactMatches(TMatches &matches, TCounts & 
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		cnts
#endif
	, RazerSOptions<TSpec> &options, TSwift &swift)
{
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping) compactMatches(matches, cnts,options,swift,true);
#endif
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
		if (readNo == (*it).rseqNo
#ifdef RAZERS_MATEPAIRS
			&& (*it).pairId == 0
#endif
			)
		{ 
			if ((*it).editDist >= editDistCutOff) continue;
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					int maxErrors = (*it).editDist - 1;
					if (options.purgeAmbiguous && (options.distanceRange == 0 || (*it).editDist < options.distanceRange))
						maxErrors = -1;

					setMaxErrors(swift, readNo, maxErrors);

					if (maxErrors == -1 && options._debugLevel >= 2)
						::std::cerr << "(read #" << readNo << " disabled)";

					if (options.purgeAmbiguous && (options.distanceRange == 0 || (*it).editDist < options.distanceRange))
						dit = ditBeg;
				}
#endif
				continue;
			}
		}
		else
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


#ifdef RAZERS_DIRECT_MAQ_MAPPING
//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template < typename TMatches, typename TCounts, typename TSpec, typename TSwift >
void compactMatches(TMatches &matches, TCounts &cnts, RazerSOptions<TSpec> &, TSwift &, bool dontCountFirstTwo)
{
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoMQ<TMatch>());
	
	unsigned readNo = -1;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;

	//number of errors may not exceed 31!
	bool second = true;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo)
		{
			//second best match
			if (second)
			{
				second = false;
				if((cnts[1][(*it).rseqNo] & 31)  > (*it).editDist)
				{
					//this second best match is better than any second best match before
					cnts[1][(*it).rseqNo] = (*it).editDist; // second best dist is this editDist
										// count is 0 (will be updated if countFirstTwo)
				}
				if(!dontCountFirstTwo) 
					if((cnts[1][(*it).rseqNo]>>5) != 2047) cnts[1][(*it).rseqNo] += 32;
			}
			else
			{
				if ((*it).editDist <= (cnts[0][(*it).rseqNo] & 31) )
					if(cnts[0][(*it).rseqNo]>>5 != 2047)
						cnts[0][(*it).rseqNo] +=32;
				if ((*it).editDist <= (cnts[1][(*it).rseqNo] & 31) )
					if((cnts[1][(*it).rseqNo]>>5) != 2047)
						cnts[1][(*it).rseqNo] +=32;
				continue;
			}
		} else
		{	//best match
			second = true;
			readNo = (*it).rseqNo;
			//cnts has 16bits, 11:5 for count:dist
			if((cnts[0][(*it).rseqNo] & 31)  > (*it).editDist)
			{
				//this match is better than any match before
				cnts[1][(*it).rseqNo] = cnts[0][(*it).rseqNo]; // best before is now second best 
									       // (count will be updated when match is removed)
				cnts[0][(*it).rseqNo] = (*it).editDist; // best dist is this editDist
									// count is 0 (will be updated if countFirstTwo)
			}
			if(!dontCountFirstTwo) 
				if((cnts[0][(*it).rseqNo]>>5) != 2047) cnts[0][(*it).rseqNo] += 32;	// shift 5 to the right, add 1, shift 5 to the left, and keep count
		}
		*dit = *it;
		++dit;
	}

	resize(matches, dit - begin(matches, Standard()));
}
#endif

//////////////////////////////////////////////////////////////////////////////
// Hamming verification
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadSet &readSet,						// reads
	TMyersPatterns const &,					// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,	// RazerS options
	SwiftSemiGlobalHamming)					// Hamming only
{
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping) 
		return matchVerify(m,inf,rseqNo,readSet,readSet /*pseudoparamter*/,options,SwiftSemiGlobalHamming(),QualityBasedScoring());
#endif

	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[rseqNo] << ::std::endl;
#endif

	unsigned ndlLength = sequenceLength(rseqNo, readSet);
	if (length(inf) < ndlLength) return false;

	// verify
	TRead &read				= readSet[rseqNo];
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
	typename TReadSet, 
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadSet &readSet,	    				// reads
	TMyersPatterns &forwardPatterns,		// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,	// RazerS options
	SwiftSemiGlobal)						// Swift specialization
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
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
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[rseqNo]<<::std::endl;
#endif

    unsigned ndlLength = sequenceLength(rseqNo, readSet);
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
	TReadRev			readRev(readSet[rseqNo]);
	TMyersFinderRev		myersFinderRev(infRev);
	TMyersPatternRev	myersPatternRev(readRev);

	_patternMatchNOfPattern(myersPatternRev, options.matchN);
	_patternMatchNOfFinder(myersPatternRev, options.matchN);
	while (find(myersFinderRev, myersPatternRev, maxScore))
		m.gBegin = m.gEnd - (position(myersFinderRev) + 1);

	return true;
}


#ifdef RAZERS_DIRECT_MAQ_MAPPING
//////////////////////////////////////////////////////////////////////////////
// Hamming verification recording sum of mismatching qualities in m.mScore
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet, 
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,					// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,				// read number
	TReadSet &readSet,				// reads
	TMyersPatterns const &,				// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,		// RazerS options
	SwiftSemiGlobalHamming const &,				// Hamming only
	QualityBasedScoring)					// MaqMapping
{
	
	typedef Segment<TGenome, InfixSegment>				TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

// #ifdef RAZERS_DEBUG
// 	cout<<"Verify: "<<::std::endl;
// 	cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
// 	cout<<"Read:   "<<host(myersPattern)<<::std::endl;
// #endif

//	bool derBesgte = false;
	//if(rseqNo == 2) derBesgte = true;
//	if(derBesgte) ::std::cout << "der besagte\n";
	unsigned ndlLength = sequenceLength(rseqNo, readSet);
	if (length(inf) < ndlLength) return false;

	// verify
	TRead &read		= readSet[rseqNo];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	TGenomeIterator git	= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);
	TGenomeIterator bestIt	= begin(inf, Standard());

	// this is max number of errors the 28bp 'seed' should have
	//assuming that maxErrors-1 matches can be found with 100% SN 
	unsigned maxErrorsSeed = (unsigned)(options.artSeedLength * options.errorRate) + 1;	
	unsigned maxErrorsTotal = (unsigned)(ndlLength * 0.25); //options.maxErrorRate);
	unsigned minErrors = maxErrorsTotal + 1;
	int minQualSumErrors = options.absMaxQualSumErrors + 10;
	unsigned minSeedErrors = maxErrorsSeed + 1;

	for (; git < gitEnd; ++git)
	{
		bool hit = true;
		unsigned errors = 0;
		unsigned count = 0;
		unsigned seedErrors = 0;
		int qualSumErrors = 0;
		TGenomeIterator g = git;	//maq would count errors in the first 28bp only (at least initially. not for output)
		for(TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
		{
			if ((options.compMask[ordValue(*g)] & options.compMask[ordValue(*r)]) == 0)
			{
			//	::std::cout << count << "<-";
				if (++errors > maxErrorsTotal)
				{
					hit = false;
					break;
				}
				//int qualityValue = (int)((unsigned char)*r >> 3);
				//qualSumErrors += (qualityValue < options.mutationRateQual) ? qualityValue : options.mutationRateQual;
				qualSumErrors += (getQualityValue(*r) < options.mutationRateQual) ? getQualityValue(*r) : options.mutationRateQual;
				if(qualSumErrors > options.absMaxQualSumErrors || qualSumErrors > minQualSumErrors)
				{
					hit = false;
					break;
				}
				if (count < options.artSeedLength)		// the maq (28bp-)seed
				{
					if(++seedErrors > maxErrorsSeed)
					{
						hit = false;
						break;
					}
					if(qualSumErrors > options.maxMismatchQualSum)
					{
						hit = false;
						break;							
					}// discard match, if 'seed' is bad (later calculations are done with the quality sum over the whole read)
				}
			}
			++count;
		}
		if (hit && (qualSumErrors < minQualSumErrors /*&& seedErrors <=maxErrorsSeed*/) ) //oder (seedErrors < minSeedErrors)
		{
			minErrors = errors;
			minSeedErrors = seedErrors;
			minQualSumErrors = qualSumErrors;
			m.gBegin = git - begin(host(inf), Standard());
			bestIt = git;
		}
	}
	if (minQualSumErrors > options.absMaxQualSumErrors || minSeedErrors > maxErrorsSeed || minErrors > maxErrorsTotal) return false;
	
	m.gEnd = m.gBegin + ndlLength;
	m.editDist = minErrors;			// errors in seed or total number of errors?
	m.mScore = minQualSumErrors;
	m.seedEditDist = minSeedErrors;
	return true;
}	


//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
	typename TMatch, 
	typename TGenome, 
	typename TReadSet,
	typename TMyersPatterns,
	typename TSpec >
inline bool
matchVerify(
	TMatch &m,								// resulting match
	Segment<TGenome, InfixSegment> inf,		// potential match genome region
	unsigned rseqNo,						// read number
	TReadSet &readSet,	    				// reads
	TMyersPatterns const & pat,				// MyersBitVector preprocessing data
	RazerSOptions<TSpec> const &options,		// RazerS options
	SwiftSemiGlobal const &swiftsemi,				// Hamming only
	QualityBasedScoring)						// Swift specialization
{
	//if(!options.maqMapping) 
		return matchVerify(m,inf,rseqNo,readSet,pat,options,swiftsemi);
	//else
	//	return matchVerify(m,inf,rseqNo,readSet,pat,options,swiftsemi); // todo!
}
#endif




#ifndef RAZERS_PARALLEL
//////////////////////////////////////////////////////////////////////////////
// Find read matches in one genome sequence
template <
	typename TMatches, 
	typename TGenome,
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TVerifier,
	typename TCounts,
	typename TSpec >
void mapSingleReads(
	TMatches &matches,				// resulting matches
	TGenome &genome,				// genome ...
	unsigned gseqNo,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> > &swiftPattern,
	TVerifier &forwardPatterns,
	TCounts & cnts,
	char orientation,				// q-gram index of reads
	RazerSOptions<TSpec> &options)
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Value<TMatches>::Type					TMatch;

	
	// FILTRATION
	typedef Finder<TGenome, Swift<TSwiftSpec> >				TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;

	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << gseqNo;
		if (orientation == 'F')
			::std::cerr << "[fwd]";
		else
			::std::cerr << "[rev]";
	}

	TReadSet &readSet = host(host(swiftPattern));
	TSwiftFinder swiftFinder(genome, options.repeatLength, 1);
	
	TMatch m = {	// to supress uninitialized warnings
		0, 0, 0, 0,
#ifdef RAZERS_MATEPAIRS
		0, 0, 0,
#endif
		0,
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		0, 0,
#endif
		0
	};

	// iterate all verification regions returned by SWIFT
	while (find(swiftFinder, swiftPattern, options.errorRate, options._debugLevel)) 
	{
		unsigned rseqNo = (*swiftFinder.curHit).ndlSeqNo;
		if (!options.spec.DONT_VERIFY && 
			matchVerify(m, range(swiftFinder, genome), rseqNo, readSet, forwardPatterns, options, TSwiftSpec()))
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
#ifdef RAZERS_MATEPAIRS
			m.pairId = 0;
			m.pairScore = 0 - m.editDist;
#endif

			if (!options.spec.DONT_DUMP_RESULTS)
			{
				appendValue(matches, m, Generous());
				if (length(matches) > options.compactThresh)
				{
					typename Size<TMatches>::Type oldSize = length(matches);
					maskDuplicates(matches);	// overlapping parallelograms cause duplicates
#ifdef RAZERS_DIRECT_MAQ_MAPPING
					if(options.maqMapping)
						compactMatches(matches, cnts, options, swiftPattern, true);
					else	
#endif
						compactMatches(matches, cnts, options, swiftPattern);
					options.compactThresh += (options.compactThresh >> 1);
					if (options._debugLevel >= 2)
						::std::cerr << '(' << oldSize - length(matches) << " matches removed)";
				}
			}

			++options.TP;
//			::std::cerr << "\"" << range(swiftFinder, genomeInf) << "\"  ";
//			::std::cerr << hstkPos << " + ";
//			::std::cerr << ::std::endl;
		} else {
			++options.FP;
//			::std::cerr << "\"" << range(swiftFinder, genomeInf) << "\"   \"" << range(swiftPattern) << "\"  ";
//			::std::cerr << rseqNo << " : ";
//			::std::cerr << hstkPos << " + ";
//			::std::cerr << bucketWidth << "  " << TP << ::std::endl;
		}
	}
}
#endif

//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences (import from Fasta)
template <
	typename TMatches, 
	typename TReadSet, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSingleReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet const &		readSet,
	TCounts & cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type				TRead;
	typedef Index<TReadSet, Index_QGram<TShape> >		TIndex;			// q-gram index
	typedef Pattern<TIndex, Swift<TSwiftSpec> >			TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier

/*	// try opening each genome file once before running the whole mapping procedure
	int filecount = 0;
	int numFiles = length(genomeFileNameList);
	while(filecount < numFiles)
	{
		::std::ifstream file;
		file.open(toCString(genomeFileNameList[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return RAZERS_GENOME_FAILED;
		file.close();
		++filecount;
	}
	*/

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
		resize(forwardPatterns, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatterns[i], indexText(swiftIndex)[i]);
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}

#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping)
	{
		resize(cnts, 2);
		for (unsigned i = 0; i < length(cnts); ++i)
			fill(cnts[i], readCount, 31); //initialize with maxeditDist, 11:5 for count:dist
	}
#endif

	// clear stats
	options.FP = 0;
	options.TP = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;

	unsigned filecount = 0;
	unsigned numFiles = length(genomeFileNameList);
	unsigned gseqNo = 0;

	// open genome files, one by one	
	while (filecount < numFiles)
	{
		// open genome file	
		::std::ifstream file;
		file.open(toCString(genomeFileNameList[filecount]), ::std::ios_base::in | ::std::ios_base::binary);
		if (!file.is_open())
			return RAZERS_GENOME_FAILED;

		// remove the directory prefix of current genome file
		::std::string genomeFile(toCString(genomeFileNameList[filecount]));
		size_t lastPos = genomeFile.find_last_of('/') + 1;
		if (lastPos == genomeFile.npos) lastPos = genomeFile.find_last_of('\\') + 1;
		if (lastPos == genomeFile.npos) lastPos = 0;
		::std::string genomeName = genomeFile.substr(lastPos);
		

		CharString	id;
		Dna5String	genome;
		unsigned gseqNoWithinFile = 0;
		// iterate over genome sequences
		SEQAN_PROTIMESTART(find_time);
		for(; !_streamEOF(file); ++gseqNo)
		{
			if (options.genomeNaming == 0)
			{
				//readID(file, id, Fasta());			// read Fasta id
				readShortID(file, id, Fasta());			// read Fasta id up to first whitespace
				appendValue(genomeNames, id, Generous());
			}
			read(file, genome, Fasta());			// read Fasta sequence
			
			gnoToFileMap.insert(::std::make_pair<unsigned,::std::pair< ::std::string,unsigned> >(gseqNo,::std::make_pair< ::std::string,unsigned>(genomeName,gseqNoWithinFile)));
			
			if (options.forward)
				mapSingleReads(matches, genome, gseqNo, swiftPattern, forwardPatterns, cnts, 'F', options);

			if (options.reverse)
			{
				reverseComplementInPlace(genome);
				mapSingleReads(matches, genome, gseqNo, swiftPattern, forwardPatterns, cnts, 'R', options);
			}
			++gseqNoWithinFile;

		}
		options.timeMapReads += SEQAN_PROTIMEDIFF(find_time);
		file.close();
		++filecount;
	}

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
	typename TCounts, 
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapSingleReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
	typedef typename Value<TReadSet>::Type				TRead;
	typedef Index<TReadSet, Index_QGram<TShape> >			TIndex;			// q-gram index
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
		resize(forwardPatterns, readCount, Exact());
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
#ifdef RAZERS_DIRECT_MAQ_MAPPING
 	if(options.maqMapping)
 	{
		resize(cnts, 2);
		for (unsigned i = 0; i < length(cnts); ++i)
			fill(cnts[i], length(readSet), 31);
	}
#endif
	
	
	
	// iterate over genome sequences
	SEQAN_PROTIMESTART(find_time);
	for(unsigned gseqNo = 0; gseqNo < length(genomeSet); ++gseqNo)
	{
		if (options.forward)
			mapSingleReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, cnts, 'F', options);

		if (options.reverse)
		{
			reverseComplementInPlace(genomeSet[gseqNo]);
			mapSingleReads(matches, genomeSet[gseqNo], gseqNo, swiftPattern, forwardPatterns, cnts, 'R', options);
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
// Wrapper for single/mate-pair mapping
template <
	typename TMatches, 
	typename TReadSet,
	typename TCounts,
	typename TSpec,
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet const &		readSet, 
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return mapMatePairReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, shape, Swift<TSwiftSpec>());
	else
#endif
		return mapSingleReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, shape, Swift<TSwiftSpec>());
}

template <
	typename TMatches, 
	typename TGenomeSet,
	typename TReadSet, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TSwiftSpec >
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet,
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options,
	TShape const &			shape,
	Swift<TSwiftSpec> const)
{
#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return mapMatePairReads(matches, genomeSet, readSet, cnts, options, shape, Swift<TSwiftSpec>());
	else
#endif
		return mapSingleReads(matches, genomeSet, readSet, cnts, options, shape, Swift<TSwiftSpec>());
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different template specializations
template <typename TMatches, typename TReadSet, typename TCounts, typename TSpec>
int mapReads(
	TMatches &				matches,
	StringSet<CharString> &	genomeFileNameList,
	StringSet<CharString> &	genomeNames,	// genome names, taken from the Fasta file
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap,
	TReadSet const &		readSet, 
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GenericShape>	gapped;

	// 2x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting shape
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, cnts, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

template <typename TMatches, typename TGenomeSet, typename TReadSet, typename TCounts, typename TSpec>
int mapReads(
	TMatches &				matches,
	TGenomeSet &			genomeSet,
	TReadSet const &		readSet, 
	TCounts &				cnts,
	RazerSOptions<TSpec> &	options)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GenericShape>	gapped;

	// 2x3 SPECIALIZATION

	if (options.hammingOnly)
	{
		// select best-fitting shape
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobalHamming>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobalHamming>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, gapped, Swift<SwiftSemiGlobalHamming>());
	} 
	else 
	{
		if (stringToShape(ungapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, ungapped, Swift<SwiftSemiGlobal>());
		
		if (stringToShape(onegapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, onegapped, Swift<SwiftSemiGlobal>());

		if (stringToShape(gapped, options.shape))
			return mapReads(matches, genomeSet, readSet, cnts, options, gapped, Swift<SwiftSemiGlobal>());
	}

	return RAZERS_INVALID_SHAPE;
}

}

#endif
