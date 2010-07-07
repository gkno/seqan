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
#include <seqan/store.h>

#ifdef RAZERS_PARALLEL_READS
	#ifdef _OPENMP
		#include <omp.h>
	#endif
#endif

#ifdef RAZERS_PARALLEL
#include "tbb/spin_mutex.h"
#endif

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// RazerS modes

	// Alignment mode
	struct RazerSLocal;
	struct RazerSGlobal;
	struct RazerSPrefix;

	// Gap mode
	struct RazerSGapped;
	struct RazerSUngapped;
	
	// Score mode
	struct RazerSErrors;
	struct RazerSScore;
	struct RazerSMAQ;
	
	template <typename TSpec = Default>
	struct RazerSQuality;

	template <typename _TAlignMode, typename _TGapMode, typename _TScoreMode>
	struct RazerSMode
	{
		typedef _TAlignMode	TAlignMode;
		typedef _TGapMode	TGapMode;
		typedef _TScoreMode	TScoreMode;
	};
	
	enum AlignMode			{ RAZERS_LOCAL, RAZERS_PREFIX, RAZERS_GLOBAL };
	enum GapMode			{ RAZERS_GAPPED, RAZERS_UNGAPPED };
	enum ScoreMode			{ RAZERS_ERRORS, RAZERS_SCORE, RAZERS_QUALITY };
	enum CompactMatchesMode	{ COMPACT, COMPACT_FINAL };

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
	// major options
		 AlignMode	alignMode;
		 GapMode	gapMode;
		 ScoreMode	scoreMode;
	
	// main options
		TSpec		spec;
		bool		forward;			// compute forward oriented read matches
		bool		reverse;			// compute reverse oriented read matches
		double		errorRate;			// Criteria 1 threshold
		unsigned	maxHits;			// output at most maxHits many matches
		unsigned	scoreDistanceRange;	// output only the best, second best, ..., scoreDistanceRange best matches
		unsigned	errorDistanceRange;	// output only matches with errors in [e..e+errorDistanceRange) according 
										// to a best match with e errors
		bool		purgeAmbiguous;		// true..remove reads with more than maxHits best matches, false..keep them
		CharString	output;				// name of result file
		int			_debugLevel;		// level of verbosity
		bool		printVersion;		// print version number
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
		unsigned	nextPairMatchId;	// use this id for the next mate-pair

	// verification parameters
		unsigned	prefixSeedLength;	// length of the prefix seed
		bool		matchN;				// false..N is always a mismatch, true..N matches with all
		unsigned char compMask[5];
		Score<int, Simple> scoringScheme;
		int			minScore;			// minimal alignment score

	// statistics
		__int64		countFiltration;	// matches returned by the filter
		__int64		countVerification;	// matches returned by the verifier
		double		timeLoadFiles;		// time for loading input files
		double		timeMapReads;		// time for mapping reads
		double		timeDumpResults;	// time for dumping the results
		
		bool		maqMapping;
		int			absMaxQualSumErrors;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		int			mutationRateQual;
		int			artSeedLength;
#endif

#ifdef RAZERS_MICRO_RNA
		bool		microRNA;
		bool 		exactSeed;
#endif			

		bool		lowMemory;		// set maximum shape weight to 13 to limit size of q-gram index
		bool		fastaIdQual;		// hidden option for special fasta+quality format we use


	// misc
		unsigned	compactThresh;		// compact match array if larger than compactThresh

	// multi-threading

#ifdef RAZERS_PARALLEL_READS
        unsigned	windowSize;
		unsigned	numberOfCores;
        double		blocksPerCore;
		unsigned	numberOfBlocks;
		unsigned	blockSize;
		unsigned	accuracy;
		unsigned	collect;
		unsigned	splitThreshold;
#endif
#ifdef RAZERS_OPENADDRESSING
		double		loadFactor;
#endif

#ifdef RAZERS_PARALLEL
		typedef ::tbb::spin_mutex	TMutex;

		TMutex		*patternMutex;
		TMutex		optionsMutex;
		TMutex		matchMutex;
#endif

		RazerSOptions() 
		{
			alignMode = RAZERS_GLOBAL;
			gapMode = RAZERS_GAPPED;
			scoreMode = RAZERS_ERRORS;
		
			forward = true;
			reverse = true;
			errorRate = 0.08;
			maxHits = 100;
			scoreDistanceRange = 0;	// disabled
			errorDistanceRange = 0; // disabled
			purgeAmbiguous = false;
			output = "";
			_debugLevel = 0;
			printVersion = false;
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

			libraryLength = 220;
			libraryError = 50;
			nextPairMatchId = 0;
			
			prefixSeedLength = 28;	// the "artificial" seed length that is used for mapping quality assignment 
			for (unsigned i = 0; i < 4; ++i)
				compMask[i] = 1 << i;
			compMask[4] = 0;

//			compactThresh = 1024;
			compactThresh = 40;

			absMaxQualSumErrors = 100;	// maximum for sum of mism qualities in total readlength
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			maqMapping = false;
			mutationRateQual = 30;		// (28bp is maq default)
			maqSeedLength = 28;
#endif

#ifdef RAZERS_MICRO_RNA
			microRNA = false;
			exactSeed = true;
#endif			

			lowMemory = false;		// set maximum shape weight to 13 to limit size of q-gram index
			fastaIdQual = false;
            
#ifdef RAZERS_PARALLEL_READS
            windowSize = 10000;
#ifdef _OPENMP
			numberOfCores = omp_get_num_procs();
#endif
			blocksPerCore = 1.0;
			numberOfBlocks = numberOfCores * blocksPerCore;
			blockSize = 0;
			accuracy = 200;
			collect = 1000000;
			splitThreshold = 100; // should be high enough to justify the overhead it is causing (scheduling, sorting)
#endif
#ifdef RAZERS_OPENADDRESSING
            loadFactor = 1.6;
#endif
		}
	};

struct MicroRNA{};	

#ifdef RAZERS_MICRO_RNA
#define RAZERS_EXTENDED_MATCH
#endif

#ifdef RAZERS_DIRECT_MAQ_MAPPING 
#define RAZERS_EXTENDED_MATCH
#endif
	
	
//////////////////////////////////////////////////////////////////////////////
// Typedefs
/*
	// definition of a Read match
	template <typename _TGPos>
	struct ReadMatch 
	{
		typedef typename _MakeSigned<_TGPos>::Type TGPos;

		unsigned		gseqNo;			// genome seqNo
		unsigned		rseqNo;			// read seqNo
		TGPos			beginPos;			// begin position of the match in the genome
		TGPos			gEnd;			// end position of the match in the genome
#ifdef RAZERS_MATEPAIRS
		unsigned		pairId;			// unique id for the two mate-pair matches (0 if unpaired)
		int				mateDelta:24;	// outer coordinate delta to the other mate 
		int				pairScore:8;	// combined score of both mates
#endif
		unsigned short	editDist;		// Levenshtein distance
#ifdef RAZERS_EXTENDED_MATCH
		short	 		mScore;
		short			seedEditDist;
#endif
		char			orientation;	// 'F'..forward strand, 'R'..reverse comp. strand
	};
*/	
	enum RAZERS_ERROR 
	{
		RAZERS_INVALID_OPTIONS = -1,
		RAZERS_READS_FAILED    = -2,
		RAZERS_GENOME_FAILED   = -3,
		RAZERS_INVALID_SHAPE   = -4
	};

//////////////////////////////////////////////////////////////////////////////
// Definitions

	typedef Dna5String									TGenome;
	typedef StringSet<TGenome>							TGenomeSet;
//	typedef Dna5String									TRead;
	typedef String<Dna5Q>								TRead;
/*#ifdef RAZERS_CONCATREADS
	typedef StringSet<TRead, Owner<ConcatDirect<> > >	TReadSet;
#else
	typedef StringSet<TRead>							TReadSet;
#endif
*/
/*	typedef ReadMatch<Difference<TGenome>::Type>		TMatch;		// a single match
	typedef String<TMatch>								TMatches;	// array of matches
*/

	template <typename TReadSet, typename TShape, typename TSpec>
	struct Cargo< Index<TReadSet, Index_QGram<TShape, TSpec> > > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};

//////////////////////////////////////////////////////////////////////////////
// Memory tuning

#ifdef RAZERS_MEMOPT

	template <typename TReadSet, typename TShape, typename TSpec>
	struct SAValue< Index<TReadSet, Index_QGram<TShape, TSpec> > > 
	{
		typedef Pair<
			unsigned,				
			unsigned,
			BitCompressed<24, 8>	// max. 16M reads of length < 256
		> Type;
	};
	
#else

	template <typename TReadSet, typename TShape, typename TSpec>
	struct SAValue< Index<TReadSet, Index_QGram<TShape, TSpec> > > 
	{
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

	template <typename TReadSet, typename TShape>
	struct Size< Index<TReadSet, Index_QGram<TShape> > >
	{
		typedef unsigned Type;
	};

	template <typename TReadSet, typename TShape>
	struct Size< Index<TReadSet, Index_QGram<TShape, OpenAddressing> > >
	{
		typedef unsigned Type;
	};
	

#ifdef RAZERS_PRUNE_QGRAM_INDEX

	//////////////////////////////////////////////////////////////////////////////
	// Repeat masker
	template <typename TReadSet, typename TShape, typename TSpec>
	inline bool _qgramDisableBuckets(Index<TReadSet, Index_QGram<TShape, TSpec> > &index) 
	{
		typedef Index<TReadSet, Index_QGram<TShape, TSpec> >	TReadIndex;
		typedef typename Fibre<TReadIndex, QGram_Dir>::Type		TDir;
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

		if (counter > 0 && cargo(index)._debugLevel >= 1)
			::std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

		return result;
	}

#endif


	template <
		typename _TFragmentStore, 
		typename _TRazerSOptions,
		typename _TRazerSMode,
		typename _TPreprocessing,
		typename _TSwiftPattern,
		typename _TCounts
	>
	struct MatchVerifier
	{
		typedef _TFragmentStore									TFragmentStore;
		typedef _TRazerSOptions									TOptions;
		typedef _TRazerSMode									TRazerSMode;
		typedef _TPreprocessing									TPreprocessing;
		typedef _TSwiftPattern									TSwiftPattern;
		typedef _TCounts										TCounts;
		
		typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
		typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
		typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
		typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
		typedef typename Size<TGenome>::Type					TSize;
		
		typedef Pattern<TRead, Myers<FindInfix, False, void>  >  TMyersPattern;
		typedef typename PatternState<TMyersPattern>::Type      TPatternState;

		TFragmentStore	*store;
		TOptions		*options;			// RazerS options
		TPreprocessing	*preprocessing;
		TSwiftPattern	*swiftPattern;
		TCounts			*cnts;

		TAlignedRead	m;
		TAlignQuality	q;
		bool			onReverseComplement;
		TSize			genomeLength;
		bool			oneMatchPerBucket;
		// TPatternState	patternState;
		
		MatchVerifier() {}
               
		MatchVerifier(_TFragmentStore &_store, TOptions &_options, TPreprocessing &_preprocessing, TSwiftPattern &_swiftPattern, TCounts &_cnts):
			store(&_store),
			options(&_options),
			preprocessing(&_preprocessing),
			swiftPattern(&_swiftPattern),
			cnts(&_cnts)
			// , patternState()
		{
			onReverseComplement = false;
			genomeLength = 0;
			oneMatchPerBucket = false;
		}

		inline void push()
		{			
			if (onReverseComplement) 
			{
				// transform coordinates to the forward strand
				m.beginPos = genomeLength - m.beginPos;
				m.endPos = genomeLength - m.endPos;
			}
						
//#pragma omp critical
// begin of critical section
			{
				if (!options->spec.DONT_DUMP_RESULTS)
				{
					m.id = length(store->alignedReadStore);
					appendValue(store->alignedReadStore, m, Generous());
					appendValue(store->alignQualityStore, q, Generous());
					if (length(store->alignedReadStore) > options->compactThresh)
					{
						typename Size<TAlignedReadStore>::Type oldSize = length(store->alignedReadStore);

						if (TYPECMP<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
							maskDuplicates(*store, TRazerSMode());	// overlapping parallelograms cause duplicates
		
						compactMatches(*store, *cnts, *options, TRazerSMode(), *swiftPattern, COMPACT);
						
						if (length(store->alignedReadStore) * 4 > oldSize)			// the threshold should not be raised
							options->compactThresh += (options->compactThresh >> 1);	// if too many matches were removed
						
//						if (options._debugLevel >= 2)
//							::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
					}
				}
				++options->countVerification;
			}
// end of critical section
		}
	};
 


//////////////////////////////////////////////////////////////////////////////
// Read a list of genome file names
template<typename TSpec>
int getGenomeFileNameList(CharString filename, StringSet<CharString> & genomeFileNames, RazerSOptions<TSpec> &options)
{
	std::fstream file;
	file.open(toCString(filename), std::ios_base::in | std::ios_base::binary);
	if(!file.is_open())
		return RAZERS_GENOME_FAILED;

	clear(genomeFileNames);
	char c = _streamGet(file);
	if (c != '>' && c != '@')	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		if (options._debugLevel >=1)
			std::cout << std::endl << "Reading multiple genome files:" << std::endl;
		
/*		//locations of genome files are relative to list file's location
		string tempGenomeFile(filename);
		size_t lastPos = tempGenomeFile.find_last_of('/') + 1;
		if (lastPos == tempGenomeFile.npos) lastPos = tempGenomeFile.find_last_of('\\') + 1;
		if (lastPos == tempGenomeFile.npos) lastPos = 0;
		string filePrefix = tempGenomeFile.substr(0,lastPos);*/

		unsigned i = 1;
		while(!_streamEOF(file))
		{ 
			_parse_skipWhitespace(file, c);
			appendValue(genomeFileNames,_parse_readFilepath(file,c));
			//CharString currentGenomeFile(filePrefix);
			//append(currentGenomeFile,_parse_readFilepath(file,c));
			//appendValue(genomeFileNames,currentGenomeFile);
			if(options._debugLevel >=2)
				std::cout <<"Genome file #"<< i <<": " << genomeFileNames[length(genomeFileNames)-1] << std::endl;
			++i;
			_parse_skipWhitespace(file, c);
		}
		if(options._debugLevel >=1)
			std::cout << i-1 << " genome files total." << std::endl;
	}
	else		//if file starts with a fasta header --> regular one-genome-file input
		appendValue(genomeFileNames, filename, Exact());
	file.close();
	return 0;

}

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences with or w/o quality values
template <typename TFSSpec, typename TFSConfig, typename TRazerSOptions>
bool loadReads(
	FragmentStore<TFSSpec, TFSConfig> &store,
	const char *fileName, 
	TRazerSOptions &options)
{
	bool countN = !(options.matchN || options.outputFormat == 1);
#ifdef RAZERS_MICRO_RNA
	if (options.microRNA) countN = false;
#endif

	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);

	unsigned seqCount = length(multiFasta);

	String<Dna5Q>	seq;
	CharString		qual;
	CharString		id;
	
	unsigned kickoutcount = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			|| options.fastaIdQual
#endif
			)
			assignSeqId(id, multiFasta[i], format);	// read Fasta id
		assignSeq(seq, multiFasta[i], format);					// read Read sequence
		assignQual(qual, multiFasta[i], format);				// read ascii quality values  
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		if(options.fastaIdQual)
		{
			qual = suffix(id, length(id) - length(seq));
			if (options.readNaming == 0)
				id = prefix(id,length(id) - length(seq));
			else 
				clear(id);
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
						clear(id);
						++kickoutcount;
						break;
					}
// low qual. reads are empty to output them and their id later as LQ reads
//			if (count > cutoffCount) continue;
		}

		// store dna and quality together
		for (unsigned j = 0; j < length(qual) && j < length(seq); ++j)
			assignQualityValue(seq[j], (int)(ordValue(qual[j]) - 33));
		if (options.trimLength > 0 && length(seq) > (unsigned)options.trimLength)
			resize(seq, options.trimLength);

		appendRead(store, seq, id);
	}
	// memory optimization
	reserve(store.readSeqStore.concat, length(store.readSeqStore.concat), Exact());
//	reserve(store.readNameStore.concat, length(store.readNameStore.concat), Exact());

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


	template <typename TAlignedReadStore, typename TLessScore>
	struct LessRNoGPos : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TLessScore lessScore;
		
		LessRNoGPos(TLessScore const &_lessScore):
			lessScore(_lessScore) {}
		
		inline bool operator() (TAlignedRead const &a, TAlignedRead const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// contig number
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// beginning position
			typename TAlignedRead::TPos ba = _min(a.beginPos, a.endPos);
			typename TAlignedRead::TPos bb = _min(b.beginPos, b.endPos);
			if (ba < bb) return true;
			if (ba > bb) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			if (oa != ob) return oa;
			
			int result = lessScore.compare(a, b);
			if (result == 0)
			{
				// prefer reads that support more of the reference
				return _max(a.beginPos, a.endPos) > _max(b.beginPos, b.endPos);
			}
			return result == -1;
		}
	};

	// ... to sort matches and remove duplicates with equal gEnd
	template <typename TAlignedReadStore, typename TLessScore>
	struct LessRNoGEndPos : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		typedef typename Value<TAlignedReadStore>::Type TAlignedRead;		
		TLessScore lessScore;
		
		LessRNoGEndPos(TLessScore const &_lessScore):
			lessScore(_lessScore) {}
		
		inline bool operator() (
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// contig number
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;

			// end position
			typename TAlignedRead::TPos ea = _max(a.beginPos, a.endPos);
			typename TAlignedRead::TPos eb = _max(b.beginPos, b.endPos);
			if (ea < eb) return true;
			if (ea > eb) return false;

			// orientation
			bool oa = a.beginPos < a.endPos;
			bool ob = b.beginPos < b.endPos;
			if (oa != ob) return oa;

			int result = lessScore.compare(a, b);
			if (result == 0)
			{
				// prefer reads that support more of the reference
				return _min(a.beginPos, a.endPos) < _min(b.beginPos, b.endPos);
			}
			return result == -1;
		}
	};

	template <typename TAlignedReadStore, typename TAlignedReadQualityStore, typename TRazerSMode>
	struct LessScore : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		TAlignedReadQualityStore &qualStore;
		
		LessScore(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline int compare(
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			typedef typename Value<TAlignedReadStore>::Type TAlignedRead;

			// read number
			if (a.readId < b.readId) return -1;
			if (a.readId > b.readId) return 1;

			// quality
			if (a.id == TAlignedRead::INVALID_ID) return 1;
			if (b.id == TAlignedRead::INVALID_ID) return -1;

			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.pairScore > qb.pairScore) return -1;
			if (qa.pairScore < qb.pairScore) return 1;
			if (qa.score > qb.score) return -1;
			if (qb.score > qa.score) return 1;
			return 0;
		}
		
		inline bool operator() (
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			return compare(a, b) == -1;
		}
	};
	
	// longest prefix mapping
	template <typename TAlignedReadStore, typename TAlignedReadQualityStore, typename TGapMode, typename TScoreMode>
	struct LessScore<TAlignedReadStore, TAlignedReadQualityStore, RazerSMode<RazerSPrefix, TGapMode, TScoreMode> > : 
		public ::std::binary_function < typename Value<TAlignedReadStore>::Type, typename Value<TAlignedReadStore>::Type, bool >
	{
		TAlignedReadQualityStore &qualStore;
		
		LessScore(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		inline int compare(
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			typedef typename Value<TAlignedReadStore>::Type TAlignedRead;

			// read number
			if (a.readId < b.readId) return -1;
			if (a.readId > b.readId) return 1;

			// quality
			if (a.id == TAlignedRead::INVALID_ID) return 1;
			if (b.id == TAlignedRead::INVALID_ID) return -1;

			typename GetValue<TAlignedReadQualityStore>::Type qa = getValue(qualStore, a.id);
			typename GetValue<TAlignedReadQualityStore>::Type qb = getValue(qualStore, b.id);
			if (qa.errors < qb.errors) return -1;
			if (qa.errors > qb.errors) return 1;
			if (qa.score > qb.score) return -1;
			if (qb.score > qa.score) return 1;
			return 0;
		}

		inline bool operator() (
			typename Value<TAlignedReadStore>::Type const &a, 
			typename Value<TAlignedReadStore>::Type const &b) const 
		{
			return compare(a, b) == -1;
		}
	};

//////////////////////////////////////////////////////////////////////////////

	template <typename TAlignedReadQualityStore, typename TRazerSMode>
	struct BinFunctorDefault
	{
		TAlignedReadQualityStore &qualStore;
		
		BinFunctorDefault(TAlignedReadQualityStore &_qualStore):
			qualStore(_qualStore) {}
		
		template <typename TAlignedRead>
		inline int operator() (TAlignedRead &alignedRead) const
		{
			return qualStore[alignedRead.id].errors;
		}
	};


//////////////////////////////////////////////////////////////////////////////
// Mark duplicate matches for deletion
template <typename TFragmentStore, typename TRazerSMode>
void maskDuplicates(TFragmentStore &store, TRazerSMode)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename TFragmentStore::TContigPos						TContigPos;
	
	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal ends

	typedef LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>	TLessScore;
	typedef LessRNoGPos<TAlignedReadStore, TLessScore>						TLessBeginPos;
	typedef LessRNoGEndPos<TAlignedReadStore, TLessScore>					TLessEndPos;
	
	sortAlignedReads(store.alignedReadStore, TLessEndPos(TLessScore(store.alignQualityStore)));

	TContigPos	beginPos = -1;
	TContigPos	endPos = -1;
	unsigned	contigId = TAlignedRead::INVALID_ID;
	unsigned	readId = TAlignedRead::INVALID_ID;
	bool		orientation = false;

	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).pairMatchId != TAlignedRead::INVALID_ID) continue;
		TContigPos itEndPos = _max((*it).beginPos, (*it).endPos);
		if (endPos == itEndPos && orientation == ((*it).beginPos < (*it).endPos) &&
			contigId == (*it).contigId && readId == (*it).readId) 
		{
			(*it).id = TAlignedRead::INVALID_ID;	// mark this alignment for deletion
			continue;
		}
		readId = (*it).readId;
		contigId = (*it).contigId;
		endPos = itEndPos;
		orientation = (*it).beginPos < (*it).endPos;
	}

	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal begins

	sortAlignedReads(store.alignedReadStore, TLessBeginPos(TLessScore(store.alignQualityStore)));

	contigId = TAlignedRead::INVALID_ID;
	it = begin(store.alignedReadStore, Standard());
	itEnd = end(store.alignedReadStore, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).id == TAlignedRead::INVALID_ID || (*it).pairMatchId != TAlignedRead::INVALID_ID) continue;

		TContigPos itBeginPos = _min((*it).beginPos, (*it).endPos);
		if (beginPos == itBeginPos && readId == (*it).readId &&
			contigId == (*it).contigId && orientation == ((*it).beginPos < (*it).endPos))
		{
			(*it).id = TAlignedRead::INVALID_ID;	// mark this alignment for deletion
			continue;
		}
		readId = (*it).readId;
		contigId = (*it).contigId;
		beginPos = itBeginPos;
		orientation = (*it).beginPos < (*it).endPos;
	}
}
/*
//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template <typename TFragmentStore, typename TCounts, typename TBinFunctor, typename TRazerSMode>
void countMatches(TFragmentStore &store, TCounts &cnt, TBinFunctor &binF, TRazerSMode)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
	typedef typename Value<TCounts>::Type							TRow;
	typedef typename Value<TRow>::Type								TValue;
	
	sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));

	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	
	unsigned readId = TAlignedRead::INVALID_ID;
	int lastBin = -1;
	__int64 count = 0;
	
	String<TValue> row, empty;
	for (; it != itEnd; ++it) 
	{
		if ((*it).id == TAlignedRead::INVALID_ID) continue;
		int bin = binF((*it).id);
		
		if (readId == (*it).readId)
		{
			if (lastBin == bin)
				++count;
			else
			{
				appendValue(row, TValue(bin, count), Generous());
				lastBin = bin;
				count = 1;
			}
		}
		else
		{
			while (length(cnt) < readId)
				appendValue(cnt, empty, Generous());
			appendValue(cnt, row, Generous());
			clear(row);
			readId = (*it).readId;
			lastBin = bin;
			count = 1;
		}
	}
	while (length(cnt) < readId)
		appendValue(cnt, empty, Generous());
	appendValue(cnt, row, Generous());
}
*/
//////////////////////////////////////////////////////////////////////////////
// Count matches for each number of errors
template <typename TFragmentStore, typename TCounts, typename TRazerSMode>
void countMatches(TFragmentStore &store, TCounts &cnt, TRazerSMode)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
	typedef typename Value<TCounts>::Type							TRow;
	typedef typename Value<TRow>::Type								TValue;
	
	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	
	unsigned readId = TAlignedRead::INVALID_ID;
	short errors = -1;
	__int64 count = 0;
	__int64 maxVal = SupremumValue<TValue>::VALUE;

	sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));

	for (; it != itEnd; ++it) 
	{
		if ((*it).id == TAlignedRead::INVALID_ID) continue;
		if (readId == (*it).readId && errors == store.alignQualityStore[(*it).id].errors)
			++count;
		else
		{
			if (readId != TAlignedRead::INVALID_ID && (unsigned)errors < length(cnt))
				cnt[errors][readId] = (maxVal < count)? (TValue)maxVal : (TValue)count;
			readId = (*it).readId;
			errors = store.alignQualityStore[(*it).id].errors;
			count = 1;
		}
	}
	if (readId != TAlignedRead::INVALID_ID && (unsigned)errors < length(cnt))
		cnt[errors][readId] = (TValue)count;
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
	if (readNo==643)
		std::cout<<"dman"<<std::endl;
	int minT = _qgramLemma(swift, readNo, maxErrors);
	if (minT > 1)
	{
//		::std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
		if (maxErrors < 0) minT = SupremumValue<int>::VALUE;
		setMinThreshold(swift, readNo, (unsigned)minT);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template <
	typename TFragmentStore,
	typename TCounts,
	typename TSpec,
	typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
	typename TSwift 
>
void compactMatches(
	TFragmentStore &store,
	TCounts &, 
	RazerSOptions<TSpec> &options, 
	RazerSMode<TAlignMode, TGapMode, TScoreMode> const,
	TSwift & swift, 
	CompactMatchesMode compactMode)
{
	typedef typename TFragmentStore::TAlignedReadStore				TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore				TAlignQualityStore;
	typedef typename Value<TAlignedReadStore>::Type					TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type	TIterator;
	typedef RazerSMode<TAlignMode, TGapMode, TScoreMode>			TRazerSMode;

	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreCutOff = InfimumValue<int>::VALUE;
	int errorCutOff = SupremumValue<int>::VALUE;
	int errorRangeBest = options.errorDistanceRange;// (TYPECMP<TScoreMode, RazerSErrors>::VALUE)? options.scoreDistanceRange: 0;
	int scoreRangeBest = (TYPECMP<TAlignMode, RazerSGlobal>::VALUE && !TYPECMP<TScoreMode, RazerSScore>::VALUE)? -(int)options.scoreDistanceRange : SupremumValue<int>::VALUE;

	sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));

	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;
	
	for (; it != itEnd; ++it) 
	{
		if ((*it).id == TAlignedRead::INVALID_ID) continue;
		int score = store.alignQualityStore[(*it).id].score;
		int errors = store.alignQualityStore[(*it).id].errors;
										if (readNo == 643) std::cerr <<"["<<score<<","<<errors<<"] "<<::std::flush;
		if (readNo == (*it).readId && (*it).pairMatchId == TAlignedRead::INVALID_ID)
		{ 
			if (score <= scoreCutOff || errors >= errorCutOff) continue;
			if (++hitCount >= hitCountCutOff)
			{
#ifdef RAZERS_MASK_READS
				if (hitCount == hitCountCutOff)
				{
					// we have enough, now look for better matches
					if (options.purgeAmbiguous && (options.scoreDistanceRange == 0 || errors < errorRangeBest || score > scoreRangeBest))
					{
						setMaxErrors(swift, readNo, -1);
						if (options._debugLevel >= 2)
							::std::cerr << "(read #" << readNo << " disabled)";
					}
					else
						// we only need better matches
						if (TYPECMP<TScoreMode, RazerSErrors>::VALUE)
						{
							setMaxErrors(swift, readNo, errors - 1);
							if (errors == 0 && options._debugLevel >= 2)
								::std::cerr << "(read #" << readNo << " disabled)";
						}

					if (options.purgeAmbiguous)
					{
						if (options.scoreDistanceRange == 0 || errors < errorRangeBest || score > scoreRangeBest || compactMode == COMPACT_FINAL)
							dit = ditBeg;
						else {
							*dit = *it;
							++dit;
						}
					}
				}
#endif
				continue;
			}
		}
		else
		{
			readNo = (*it).readId;
			hitCount = 0;
			if (options.scoreDistanceRange > 0)	scoreCutOff = score - options.scoreDistanceRange;
			if (options.errorDistanceRange > 0)	errorCutOff = errors + options.errorDistanceRange;
			ditBeg = dit;
		}
		*dit = *it;
		++dit;
	}
	resize(store.alignedReadStore, dit - begin(store.alignedReadStore, Standard()));
	compactAlignedReads(store);
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template <
	typename TFragmentStore,
	typename TCounts,
	typename TSpec,
	typename TGapMode,
	typename TSwift 
>
void compactMatches(
	TFragmentStore &store,
	TCounts & cnts, 
	RazerSOptions<TSpec> &options, 
	RazerSMode<RazerSGlobal, TGapMode, RazerSQuality<RazerSMAQ> > const,
	TSwift & swift, 
	CompactMatchesMode compactMode)
{
	typedef typename TFragmentStore::TAlignedReadStore						TAlignedReadStore;
	typedef typename TFragmentStore::TAlignQualityStore						TAlignQualityStore;
	typedef typename Value<TAlignedReadStore>::Type							TAlignedRead;
	typedef typename Iterator<TAlignedReadStore, Standard>::Type			TIterator;
	typedef RazerSMode<RazerSGlobal, TGapMode, RazerSQuality<RazerSMAQ> >	TRazerSMode;
	
	unsigned readNo = -1;

	sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));
		
	TIterator it = begin(store.alignedReadStore, Standard());
	TIterator itEnd = end(store.alignedReadStore, Standard());
	TIterator dit = it;

	//number of errors may not exceed 31!
	bool second = true;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).readId)
		{
			//second best match
			if (second)
			{
				second = false;
				if((cnts[1][(*it).readId] & 31)  > (*it).editDist)
				{
					//this second best match is better than any second best match before
					cnts[1][(*it).readId] = (*it).editDist; // second best dist is this editDist
										// count is 0 (will be updated if countFirstTwo)
				}
				if (compactMode == COMPACT_FINAL) 
					if((cnts[1][(*it).readId]>>5) != 2047) cnts[1][(*it).readId] += 32;
			}
			else
			{
				if ((*it).editDist <= (cnts[0][(*it).readId] & 31) )
					if(cnts[0][(*it).readId]>>5 != 2047)
						cnts[0][(*it).readId] +=32;
				if ((*it).editDist <= (cnts[1][(*it).readId] & 31) )
					if((cnts[1][(*it).readId]>>5) != 2047)
						cnts[1][(*it).readId] +=32;
				continue;
			}
		} else
		{	//best match
			second = true;
			readNo = (*it).readId;
			//cnts has 16bits, 11:5 for count:dist
			if((cnts[0][(*it).readId] & 31)  > (*it).editDist)
			{
				//this match is better than any match before
				cnts[1][(*it).readId] = cnts[0][(*it).readId]; // best before is now second best 
									       // (count will be updated when match is removed)
				cnts[0][(*it).readId] = (*it).editDist; // best dist is this editDist
									// count is 0 (will be updated if countFirstTwo)
			}
			if (compactMode == COMPACT_FINAL) 
				if((cnts[0][(*it).readId]>>5) != 2047) cnts[0][(*it).readId] += 32;	// shift 5 to the right, add 1, shift 5 to the left, and keep count
		}
		*dit = *it;
		++dit;
	}

	resize(store.alignedReadStore, dit - begin(store.alignedReadStore, Standard()));
	resize(store.alignQualityStore, length(store.alignedReadStore));
}

/* // fallback
template <
	typename TFragmentStore,
	typename TCounts,
	typename TSpec,
	typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
	typename TSwift 
>
void compactMatches(
	TFragmentStore &,
	TCounts &, 
	RazerSOptions<TSpec> &, 
	RazerSMode<TAlignMode, TGapMode, TScoreMode> const,
	TSwift &, 
	CompactMatchesMode)
{
}
*/

//////////////////////////////////////////////////////////////////////////////
// Best Hamming prefix verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,									// potential match genome region
	unsigned readId,													// read number
	TReadSet &readSet,													// reads
	RazerSMode<RazerSPrefix, RazerSUngapped, RazerSErrors> const)		// Hamming only
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

	// verify
	TRead &read				= readSet[readId];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	unsigned ndlLength		= ritEnd - ritBeg;

	if (length(inf) < ndlLength) return false;
	TGenomeIterator git		= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);
	
	unsigned maxErrors = (unsigned)(verifier.options->prefixSeedLength * verifier.options->errorRate);
	unsigned minErrors = maxErrors + 2;
	unsigned errorThresh = (verifier.oneMatchPerBucket)? SupremumValue<unsigned>::VALUE: maxErrors;
	int bestHitLength = 0;

	for (; git < gitEnd; ++git)
	{
		unsigned errors = 0;
		TReadIterator r = ritBeg;
		TGenomeIterator g = git;
		for (; r != ritEnd; ++r, ++g)
			if ((verifier.options->compMask[ordValue(*g)] & verifier.options->compMask[ordValue(*r)]) == 0)
			{
				if (r - ritBeg < (int)verifier.options->prefixSeedLength)	// seed
				{
					if (++errors > maxErrors)				// doesn't work for islands with errorThresh > maxErrors
						break;
				}
				else
					break;
			}
			
		if (errors < minErrors)
		{
			minErrors = errors;
			bestHitLength = r - ritBeg;
			verifier.m.beginPos = git - begin(host(inf), Standard());
		} else if (errors == minErrors && bestHitLength < r - ritBeg)
		{
			bestHitLength = r - ritBeg;
			verifier.m.beginPos = git - begin(host(inf), Standard());
		} else if (errorThresh < errors)
		{
			if (minErrors <= maxErrors)
			{
				verifier.m.endPos = verifier.m.beginPos + bestHitLength;
				verifier.q.pairScore = verifier.q.score = bestHitLength;
				verifier.q.errors = minErrors;
				verifier.push();
				minErrors = maxErrors + 2;
			}
		}
	}

	if (minErrors <= maxErrors)
	{
		verifier.m.endPos = verifier.m.beginPos + bestHitLength;
		verifier.q.pairScore = verifier.q.score = bestHitLength;
		verifier.q.errors = minErrors;
		if (!verifier.oneMatchPerBucket)
			verifier.push();
		return true;
	}
	return false;
}	


template <typename TRazerSMode>
struct __UseQualityValues { enum { VALUE = false }; };
template <typename TAlignMode, typename TGapMode, typename TSpec>
struct __UseQualityValues<RazerSMode<TAlignMode, TGapMode, RazerSQuality<TSpec> > > {  enum { VALUE = true }; };

//////////////////////////////////////////////////////////////////////////////
// Hamming verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet,
	typename TScoreMode >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,								// potential match genome region
	unsigned readId,												// read number
	TReadSet &readSet,												// reads
	RazerSMode<RazerSGlobal, RazerSUngapped, TScoreMode> const)		// Semi-global, no gaps
{
	typedef Segment<TGenome, InfixSegment>							TGenomeInfix;
	typedef typename Value<TReadSet>::Type const					TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type			TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type				TReadIterator;
	typedef RazerSMode<RazerSGlobal, RazerSUngapped, TScoreMode>	TRazerSMode;

#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[readId] << ::std::endl;
#endif

	// verify
	TRead &read				= readSet[readId];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	unsigned ndlLength		= ritEnd - ritBeg;

	if (length(inf) < ndlLength) return false;
	TGenomeIterator git		= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);
	
	int mismatchDelta, scoreInit;
	int minScore;
	if (TYPECMP<TScoreMode, RazerSErrors>::VALUE)
		minScore = -(int)(ndlLength * verifier.options->errorRate);
	else if (__UseQualityValues<TRazerSMode>::VALUE)
		minScore = -verifier.options->absMaxQualSumErrors;
	else if (TYPECMP<TScoreMode, RazerSScore>::VALUE)
	{
		minScore = verifier.options->minScore;
		mismatchDelta = scoreMatch(verifier.options->scoringScheme) - scoreMismatch(verifier.options->scoringScheme);
		scoreInit = scoreMatch(verifier.options->scoringScheme) * ndlLength;
	}
	
	int maxScore = minScore - 1;
	int scoreThresh = (verifier.oneMatchPerBucket)? SupremumValue<int>::VALUE: minScore;
	int score, errors;
	
	for (; git < gitEnd; ++git)
	{
		if (!TYPECMP<TScoreMode, RazerSScore>::VALUE)
			score = 0;
		else
			score = scoreInit;
		
		if (!TYPECMP<TScoreMode, RazerSErrors>::VALUE)
			errors = 0;
		
		TGenomeIterator g = git;
		for (TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
			if ((verifier.options->compMask[ordValue(*g)] & verifier.options->compMask[ordValue(*r)]) == 0)
			{
				if (TYPECMP<TScoreMode, RazerSErrors>::VALUE)
				{
					// A. Count mismatches only
					--score;
				} else
				{
					++errors;
					if (__UseQualityValues<TRazerSMode>::VALUE)
						// B. Count mismatches and mismatch qualities
						score -= getQualityValue(*g);
					else if (TYPECMP<TScoreMode, RazerSScore>::VALUE)
						// C. Count mismatches and alignment score
						score -= mismatchDelta;
					else
						std::cerr << "Unsupported score mode" << std::endl;
				}
				if (score < minScore)	// doesn't work for islands with errorThresh > maxErrors
					break;
			}
		
		if (score > maxScore)
		{
			maxScore = score;
			if (TYPECMP<TScoreMode, RazerSErrors>::VALUE)
				verifier.q.errors = -score;
			else
				verifier.q.errors = errors;
			verifier.m.beginPos = git - begin(host(inf), Standard());
		} else if (scoreThresh > score)
		{
			if (maxScore >= minScore)
			{
				// for RazerSErrors bestErrors == -maxScore
				verifier.m.endPos = verifier.m.beginPos + ndlLength;
				verifier.q.pairScore = verifier.q.score = maxScore;
				verifier.push();
				maxScore = minScore - 1;
			}
		}
	}

	if (maxScore >= minScore)
	{
		verifier.m.endPos = verifier.m.beginPos + ndlLength;
		verifier.q.pairScore = verifier.q.score = maxScore;
		if (!verifier.oneMatchPerBucket)
			verifier.push();
		return true;
	}
	return false;
}	


//////////////////////////////////////////////////////////////////////////////
// Edit distance verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,									// potential match genome region
	unsigned readId,														// read number
	TReadSet &readSet,													// reads
	RazerSMode<RazerSGlobal, RazerSGapped, RazerSErrors> const)		// Mismatches and Indels
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type					TRead;
	typedef typename Position<TGenomeInfix>::Type			TPosition;

	// find read match end
	typedef Finder<TGenomeInfix>								TMyersFinder;
	typedef typename TMatchVerifier::TPreprocessing		TPreprocessing;
	typedef typename Value<TPreprocessing>::Type			TMyersPattern;
	typedef typename PatternState<TMyersPattern>::Type		TPatternState;

	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>		TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>				TReadRev;
	typedef Finder<TGenomeInfixRev>							TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>			TMyersPatternRev;

	TMyersFinder myersFinder(inf);
	TMyersPattern &myersPattern = (*verifier.preprocessing)[readId];
	// TPatternState & state = verifier.patternState;
	// TPatternState state;
	
#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[readId]<<::std::endl;
#endif

    unsigned ndlLength = sequenceLength(readId, readSet);
	int maxScore = InfimumValue<int>::VALUE;
	int minScore = -(int)(ndlLength * verifier.options->errorRate);
	TPosition maxPos = 0;
	TPosition lastPos = length(inf);
	unsigned minDistance = (verifier.oneMatchPerBucket)? lastPos: 1;

	// find end of best semi-global alignment
	// while (find(myersFinder, myersPattern, state, minScore))
	while (find(myersFinder, myersPattern, minScore))
	{
		TPosition pos = position(hostIterator(myersFinder));
		if (lastPos + minDistance < pos)
		{
			if (minScore <= maxScore)
			{
				verifier.m.endPos = beginPosition(inf) + maxPos + 1;
				verifier.q.errors = -maxScore;

				if ((verifier.q.pairScore = verifier.q.score = maxScore) == 0)
					verifier.m.beginPos = verifier.m.endPos - ndlLength;
				else
				{
					TPosition infBeginPos = beginPosition(inf);
					TPosition infEndPos = endPosition(inf);

					// find beginning of best semi-global alignment
					TGenomeInfixRev infRev(inf);
					setEndPosition(inf, verifier.m.endPos = (beginPosition(inf) + maxPos + 1));

					// limit the beginning to needle length plus errors (== -maxScore)
					if (length(inf) > ndlLength - maxScore)
						setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
					
					TReadRev			readRev(readSet[readId]);
					TMyersFinderRev		myersFinderRev(infRev);
					TMyersPatternRev	myersPatternRev(readRev);

					_patternMatchNOfPattern(myersPatternRev, verifier.options->matchN);
					_patternMatchNOfFinder(myersPatternRev, verifier.options->matchN);
					while (find(myersFinderRev, myersPatternRev, maxScore))
						verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);

					setBeginPosition(inf, infBeginPos);
					setEndPosition(inf, infEndPos);
				}
				verifier.push();
				maxScore = minScore - 1;
			}
		}
		if (getScore(myersPattern) >= maxScore)
		// if (getScore(state) >= maxScore)
		{
			// maxScore = getScore(state);
			maxScore = getScore(myersPattern);
			maxPos = pos;
		}
		lastPos = pos;
	}
	
	if (minScore <= maxScore)
	{
		verifier.m.endPos = beginPosition(inf) + maxPos + 1;
		verifier.q.errors = -maxScore;

		if ((verifier.q.pairScore = verifier.q.score = maxScore) == 0)
			verifier.m.beginPos = verifier.m.endPos - ndlLength;
		else
		{
			setEndPosition(inf, verifier.m.endPos);

			// limit the beginning to needle length plus errors (== -maxScore)
			if (length(inf) > ndlLength - maxScore)
				setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
			
			// find beginning of best semi-global alignment
			TGenomeInfixRev		infRev(inf);
			TReadRev			readRev(readSet[readId]);
			TMyersFinderRev		myersFinderRev(infRev);
			TMyersPatternRev	myersPatternRev(readRev);

			_patternMatchNOfPattern(myersPatternRev, verifier.options->matchN);
			_patternMatchNOfFinder(myersPatternRev, verifier.options->matchN);
			while (find(myersFinderRev, myersPatternRev, maxScore))
				verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);
		}

		if (!verifier.oneMatchPerBucket)
			verifier.push();
        
#ifdef RAZERS_DEBUG
        std::cout << "OK" << std::endl; 
#endif
		return true;
	}
#ifdef RAZERS_DEBUG
    std::cout << "FAILED" << std::endl; 
#endif
	return false;
}

#ifdef RAZERS_DIRECT_MAQ_MAPPING

//////////////////////////////////////////////////////////////////////////////
// Best Hamming prefix verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,												// potential match genome region
	unsigned readId,																	// read number
	TReadSet &readSet,																// reads
	RazerSMode<RazerSGlobal, RazerSUngapped, RazerSQuality<RazerSMAQ> > const)	// Hamming only
{
	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type	TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type		TReadIterator;

	// verify
	TRead &read				= readSet[readId];
	TReadIterator ritBeg	= begin(read, Standard());
	TReadIterator ritEnd	= end(read, Standard());
	unsigned ndlLength		= ritEnd - ritBeg;

	if (length(inf) < ndlLength) return false;
	TGenomeIterator git		= begin(inf, Standard());
	TGenomeIterator gitEnd	= end(inf, Standard()) - (ndlLength - 1);
	
	unsigned maxErrors = (unsigned)(verifier.options->prefixSeedLength * verifier.options->errorRate);
	unsigned minErrors = 0;
	unsigned minQualSum = verifier.options->absMaxQualSumErrors + 1;
	unsigned qualSumThresh = (verifier.oneMatchPerBucket)? SupremumValue<unsigned>::VALUE: verifier.options->absMaxQualSumErrors;
	unsigned bestHitLength = 0;

	for (; git < gitEnd; ++git)
	{
		unsigned seedErrors = 0;
		unsigned errors = 0;
		unsigned qualSum = 0;
		TReadIterator r = ritBeg;
		TGenomeIterator g = git;
		for (; r != ritEnd; ++r, ++g)
			if ((verifier.options->compMask[ordValue(*g)] & verifier.options->compMask[ordValue(*r)]) == 0)
			{
				if (r - ritBeg < verifier.options->prefixSeedLength)		// seed
				{
					if (++seedErrors > maxErrors)
					{
						qualSum = verifier.options->absMaxQualSumErrors + 1;
						break;
					}
				}
				qualSum += (getQualityValue(*r) < verifier.options->mutationRateQual) ? getQualityValue(*r) : verifier.options->mutationRateQual;
				if (qualSum > verifier.options->absMaxQualSumErrors)
					break;
				++errors;
			}
			
		if (verifier.options->prefixSeedLength != 0) errors = seedErrors;
		if (qualSum < minQualSum)
		{
			minQualSum = qualSum;
			minErrors = errors;
			verifier.m.beginPos = git - begin(host(inf), Standard());
		} else if (qualSum == minQualSum && errors < minErrors)
		{
			minErrors = errors;
			verifier.m.beginPos = git - begin(host(inf), Standard());
		} else if (qualSumThresh < qualSum)
		{
			if (minQualSum <= verifier.options->absMaxQualSumErrors)
			{
				verifier.m.endPos = verifier.m.beginPos + ndlLength;
				verifier.q.pairScore = verifier.q.score = -(int)minQualSum;
				verifier.q.errors = minErrors;
				verifier.push();
				minQualSum = verifier.options->absMaxQualSumErrors + 1;
			}
		}
	}

	if (minErrors <= maxErrors)
	{
		verifier.m.endPos = verifier.m.beginPos + ndlLength;
		verifier.q.pairScore = verifier.q.score = -(int)minQualSum;
		verifier.q.errors = minErrors;
		if (!verifier.oneMatchPerBucket)
			verifier.push();
		return true;
	}
	return false;
}

#endif

template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet,
	typename TAlignMode,
	typename TGapMode,
	typename TScoreMode >
inline bool
matchVerify(
	TMatchVerifier &,
	Segment<TGenome, InfixSegment>,								// potential match genome region
	unsigned,													// read number
	TReadSet &,													// reads
	RazerSMode<TAlignMode, TGapMode, TScoreMode> const)
{
	std::cerr << "Verification not implemenented!" << std::endl;
	return false;
}


#ifndef RAZERS_PARALLEL
//////////////////////////////////////////////////////////////////////////////
// Find read matches in a single genome sequence
template <
	typename TFragmentStore, 
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TPreprocessing,
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode >
void _mapSingleReadsToContig(
	TFragmentStore							& store,
	unsigned								  contigId,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> >	& swiftPattern,
	TPreprocessing							& preprocessing,
	TCounts									& cnts,
	char									  orientation,				// q-gram index of reads
	TRazerSOptions							& options,
	TRazerSMode						  const & mode)
{
	// FILTRATION
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;
	
	// VERIFICATION
	typedef MatchVerifier <
		TFragmentStore, 
		TRazerSOptions, 
		TRazerSMode,
		TPreprocessing, 
		TSwiftPattern,
		TCounts >											TVerifier;
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	
	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F') ::std::cerr << "[fwd]";
		else                    ::std::cerr << "[rev]";
	}
	lockContig(store, contigId);
	TContigSeq &contigSeq = store.contigStore[contigId].seq;
	if (orientation == 'R')	reverseComplementInPlace(contigSeq);

	TReadSet		&readSet = host(host(swiftPattern));
	TSwiftFinder	swiftFinder(contigSeq, options.repeatLength, 1);
	TVerifier		verifier(store, options, preprocessing, swiftPattern, cnts);
    
	// initialize verifier
	verifier.onReverseComplement = (orientation == 'R');
	verifier.genomeLength = length(contigSeq);
	verifier.m.contigId = contigId;
		
	// iterate all verification regions returned by SWIFT
	while (find(swiftFinder, swiftPattern, options.errorRate)) 
	{
		verifier.m.readId = (*swiftFinder.curHit).ndlSeqNo;
		if (!options.spec.DONT_VERIFY)
			matchVerify(verifier, infix(swiftFinder), verifier.m.readId, readSet, mode);
		++options.countFiltration;
	}
	if (!unlockAndFreeContig(store, contigId))							// if the contig is still used
		if (orientation == 'R')	reverseComplementInPlace(contigSeq);	// we have to restore original orientation
}
#endif


//////////////////////////////////////////////////////////////////////////////
// Find read matches in many genome sequences
template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TCounts,
	typename TSpec, 
	typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
	typename TReadIndex>
int _mapSingleReads(
	FragmentStore<TFSSpec, TFSConfig>					& store,
	TCounts												& cnts,
	RazerSOptions<TSpec>								& options,
	RazerSMode<TAlignMode, TGapMode, TScoreMode>  const & mode,
	TReadIndex											& readIndex)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename IF<
				TYPECMP<TGapMode,RazerSGapped>::VALUE,
				SwiftSemiGlobal,
				SwiftSemiGlobalHamming>::Type			TSwiftSpec;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >		TSwiftPattern;	// filter
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier
	// typedef Pattern<TRead, Myers<FindInfix, False, void> >	TMyersPattern;	// verifier
	
	// configure Swift pattern
	TSwiftPattern swiftPattern(readIndex);
	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;
	swiftPattern.params.printDots = options._debugLevel > 0; 

	// init edit distance verifiers
	unsigned readCount = countSequences(readIndex);
	String<TMyersPattern> forwardPatterns;
	options.compMask[4] = (options.matchN)? 15: 0;
	if (options.gapMode == RAZERS_GAPPED)
	{
		resize(forwardPatterns, readCount, Exact());
		for(unsigned i = 0; i < readCount; ++i)
		{
			setHost(forwardPatterns[i], indexText(readIndex)[i]);
			_patternMatchNOfPattern(forwardPatterns[i], options.matchN);
			_patternMatchNOfFinder(forwardPatterns[i], options.matchN);
		}
	}
	
	if (options.maqMapping)
	{
		resize(cnts, 2);
		for (unsigned i = 0; i < length(cnts); ++i)
			fill(cnts[i], readCount, 31); //initialize with maxeditDist, 11:5 for count:dist
	}

	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;
	SEQAN_PROTIMESTART(find_time);
	
	// iterate over genome sequences
    #ifdef RAZERS_PARALLEL_CONTIGS
    #pragma omp parallel for private(swiftPattern)
    #endif
    for (int contigId = 0; contigId < (int)length(store.contigStore); ++contigId)
	{
		// lock to prevent releasing and loading the same contig twice
		// (once per _mapSingleReadsToContig call)
		lockContig(store, contigId);
#ifndef RAZERS_WINDOW
		if (options.forward)
			_mapSingleReadsToContig(store, contigId, swiftPattern, forwardPatterns, cnts, 'F', options, mode);
		if (options.reverse)
			_mapSingleReadsToContig(store, contigId, swiftPattern, forwardPatterns, cnts, 'R', options, mode);
#else
		if (options.forward)
			_mapSingleReadsToContigWindow(store, contigId, swiftPattern, forwardPatterns, cnts, 'F', options, mode);
		if (options.reverse)
			_mapSingleReadsToContigWindow(store, contigId, swiftPattern, forwardPatterns, cnts, 'R', options, mode);
#endif
		unlockAndFreeContig(store, contigId);
	}
	
	options.timeMapReads = SEQAN_PROTIMEDIFF(find_time);
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
	if (options._debugLevel >= 2) {
		::std::cerr << ::std::endl;
		::std::cerr << "___FILTRATION_STATS____" << ::std::endl;
		::std::cerr << "Filtration counter:      " << options.countFiltration << ::std::endl;
		::std::cerr << "Successful verfications: " << options.countVerification << ::std::endl;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Wrapper for SWIFT (default)
template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TAlignMode,
	typename TGapMode,
	typename TScoreMode >
int _mapSingleReads(
	FragmentStore<TFSSpec, TFSConfig>					& store,
	TCounts												& cnts,
	RazerSOptions<TSpec>								& options,
	TShape const										& shape,
	RazerSMode<TAlignMode, TGapMode, TScoreMode>  const & mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
#ifndef RAZERS_OPENADDRESSING
	typedef Index<TReadSeqStore, Index_QGram<TShape> >	TIndex;			// q-gram index
#else
	typedef Index<TReadSeqStore, Index_QGram<TShape, OpenAddressing> >	TIndex;
#endif
	
	// configure q-gram index
	TIndex swiftIndex(store.readSeqStore, shape);
#ifdef RAZERS_OPENADDRESSING
	swiftIndex.alpha = options.loadFactor;
#endif
	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;

	return _mapSingleReads(store, cnts, options, mode, swiftIndex);
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for SWIFT with Micro RNA
template <
	typename TFSSpec, 
	typename TFSConfig, 
	typename TCounts,
	typename TSpec, 
	typename TShape,
	typename TGapMode,
	typename TScoreMode >
int _mapSingleReads(
	FragmentStore<TFSSpec, TFSConfig>						& store,
	TCounts													& cnts,
	RazerSOptions<TSpec>									& options,
	TShape const											& shape,
	RazerSMode<RazerSPrefix, TGapMode, TScoreMode>    const & mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>				TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore			TReadSeqStore;
	
	typedef typename Value<TReadSeqStore>::Type				TRead;
	typedef typename Infix<TRead>::Type						TReadInfix;
	typedef StringSet<TReadInfix>							TReadSet;
	typedef Index<TReadSet, Index_QGram<TShape> >			TIndex;			// q-gram index

	TReadSet readSet;
	unsigned readCount = length(store.readSeqStore);
	resize(readSet, readCount, Exact());

	for (unsigned i = 0; i < readCount; ++i)
		assign(readSet[i], prefix(store.readSeqStore[i], _min(length(store.readSeqStore[i]), options.prefixSeedLength)));

	// configure q-gram index
	TIndex swiftIndex(readSet, shape);
	cargo(swiftIndex).abundanceCut = options.abundanceCut;
	cargo(swiftIndex)._debugLevel = options._debugLevel;
	
	return _mapSingleReads(store, cnts, options, mode, swiftIndex);
}


//////////////////////////////////////////////////////////////////////////////
// Wrapper for single/mate-pair mapping
template <
	typename TFSSpec, 
	typename TFSConfig,
	typename TCounts,
	typename TSpec,
	typename TShape,
	typename TRazerSMode >
int _mapReads(
	FragmentStore<TFSSpec, TFSConfig>		& store,
	TCounts									& cnts,
	RazerSOptions<TSpec>					& options,
	TShape const							& shape,
	TRazerSMode						  const & mode)
{
#ifdef RAZERS_MATEPAIRS
	if (options.libraryLength >= 0)
		return _mapMatePairReads(store, cnts, options, shape, mode);
	else
#endif
#ifndef RAZERS_PARALLEL_READS
        return _mapSingleReads(store, cnts, options, shape, mode);
#else
        return _mapSingleReadsParallel(store, cnts, options, shape, mode);
#endif
}


//////////////////////////////////////////////////////////////////////////////
// Wrapper for different shapes
template <typename TFSSpec, typename TFSConfig, typename TCounts, typename TSpec, typename TRazersMode>
int _mapReads(
	FragmentStore<TFSSpec, TFSConfig>		& store,
	TCounts									& cnts,
	RazerSOptions<TSpec>					& options,
	TRazersMode						  const & mode)
{
	Shape<Dna, SimpleShape>		ungapped;
	Shape<Dna, OneGappedShape>	onegapped;
	Shape<Dna, GenericShape>	gapped;

	// 2x3 SPECIALIZATION

	// select best-fitting shape
	if (stringToShape(ungapped, options.shape))
		return _mapReads(store, cnts, options, ungapped, mode);
	if (stringToShape(onegapped, options.shape))
		return _mapReads(store, cnts, options, onegapped, mode);
	if (stringToShape(gapped, options.shape))
		return _mapReads(store, cnts, options, gapped, mode);
	return RAZERS_INVALID_SHAPE;
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different score modes
template <typename TFSSpec, typename TFSConfig, typename TCounts, typename TSpec, typename TAlignMode, typename TGapMode>
int _mapReads(
	FragmentStore<TFSSpec, TFSConfig>		& store,
	TCounts									& cnts,
	RazerSOptions<TSpec>					& options,
	RazerSMode<TAlignMode, TGapMode, Nothing> const)
{
	if (options.scoreMode == RAZERS_ERRORS)
		return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSErrors>());
	if (options.scoreMode == RAZERS_SCORE)
		return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSScore>());
	if (options.scoreMode == RAZERS_QUALITY)
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		if (options.maqMapping)
			return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSQuality<RazerSMAQ> >());
		else
#endif
			return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSQuality<> >());
	return RAZERS_INVALID_OPTIONS;
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper for different gap and align modes
template <typename TFSSpec, typename TFSConfig, typename TCounts, typename TSpec>
int _mapReads(
	FragmentStore<TFSSpec, TFSConfig>		& store,
	TCounts									& cnts,
	RazerSOptions<TSpec>					& options)
{
	if (options.gapMode == RAZERS_GAPPED)
	{
		if (options.alignMode == RAZERS_LOCAL)
			return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSGapped, Nothing>());
		if (options.alignMode == RAZERS_PREFIX)
			return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSGapped, Nothing>());
		if (options.alignMode == RAZERS_GLOBAL)
			return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSGapped, Nothing>());
	} else 
	{
		if (options.alignMode == RAZERS_LOCAL)
			return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSUngapped, Nothing>());
		if (options.alignMode == RAZERS_PREFIX)
			return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSUngapped, Nothing>());
		if (options.alignMode == RAZERS_GLOBAL)
			return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSUngapped, Nothing>());
	}
	return RAZERS_INVALID_OPTIONS;
}

}

#endif
