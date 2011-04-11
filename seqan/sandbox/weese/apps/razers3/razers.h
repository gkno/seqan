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

// TODO(holtgrew): I probably broke maq mapping.
// TODO(holtgrew): What about longest prefix mapping stuff?

#include <iostream>
#include <fstream>

#include <omp.h>

#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/pipe.h>

#ifdef RAZERS_PROFILE
#include "profile_timeline.h"
#endif  // #ifdef RAZERS_PROFILE

// No parallelism for less than MIN_PARALLEL_WORK reads.
const unsigned MIN_PARALLEL_WORK = 1;//100/*0*/; // TODO(holtgrew): Set to some useful value after development.

namespace SEQAN_NAMESPACE_MAIN
{

// Compact representation of a match.
template <typename TContigPos_>
struct MatchRecord 
{
    typedef typename MakeSigned_<TContigPos_>::Type TContigPos;

    unsigned        contigId;       // genome seqNo
    unsigned        readId;         // read seqNo
    TContigPos      beginPos;       // begin position of the match in the genome
    TContigPos      endPos;         // end position of the match in the genome
    char            orientation;    // 'F', 'R', '-'
    short int       score;          // Levenshtein distance / score.
    unsigned		pairMatchId;			// unique id for the two mate-pair matches (0 if unpaired)
#ifndef RAZERS_DEFER_COMPACTION
    int				mateDelta:24;	// outer coordinate delta to the other mate 
    int				pairScore:8;	// combined score of both mates
#else
    bool            isRegistered:1; // registered in masking process.
    int				mateDelta:23;	// outer coordinate delta to the other mate 
    int				pairScore:8;	// combined score of both mates
#endif  // #ifndef RAZERS_DEFER_COMPACTION

	static const unsigned INVALID_ID;

    MatchRecord()
            : contigId(MaxValue<unsigned>::VALUE), readId(MaxValue<unsigned>::VALUE),
              beginPos(0), endPos(0), orientation('-'), score(0),
              pairMatchId(MaxValue<unsigned>::VALUE),
#ifdef RAZERS_DEFER_COMPACTION
              isRegistered(false),
#endif  // #ifndef RAZERS_DEFER_COMPACTION
              mateDelta(0), pairScore(0)
    {}
};

template <typename TStream, typename TPos>
TStream & 
operator<<(TStream & stream, MatchRecord<TPos> & record)
{
    stream << "(contigId=" << record.contigId << ", readId=" << record.readId << ", beginPos=" << record.beginPos << ", endPos = " << record.endPos << ", orientation=" << record.orientation << ", score=" << record.score << ", pairMatchId=" << record.pairMatchId << ", mateDelta=" << record.mateDelta << ", pairScore=" << record.pairScore << ")";
    return stream;
}

template <typename TGPos_>
const unsigned MatchRecord<TGPos_>::INVALID_ID = MaxValue<unsigned>::VALUE;

#ifdef RAZERS_PROFILE
enum {
    TASK_WAIT,
    TASK_ON_CONTIG,
    TASK_INIT,
    TASK_REVCOMP,
    TASK_FILTER,
    TASK_VERIFY,
    TASK_WRITEBACK,
    TASK_COMPACT,
    TASK_DUMP_MATCHES,
    TASK_LOAD,
    TASK_SORT
};
#endif  // #ifdef RAZERS_PROFILE

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

    template <typename TAlignMode_, typename TGapMode_, typename TScoreMode_, typename TMatchNPolicy_>
	struct RazerSMode
	{
		typedef TAlignMode_	TAlignMode;
		typedef TGapMode_	TGapMode;
		typedef TScoreMode_	TScoreMode;
		typedef TMatchNPolicy_	TMatchNPolicy;
	};
	
	enum AlignMode			{ RAZERS_LOCAL, RAZERS_PREFIX, RAZERS_GLOBAL };
	enum GapMode			{ RAZERS_GAPPED, RAZERS_UNGAPPED };
	enum ScoreMode			{ RAZERS_ERRORS, RAZERS_SCORE, RAZERS_QUALITY };
    enum CompactMatchesMode	{ COMPACT, COMPACT_FINAL, COMPACT_FINAL_EXTERNAL };

//////////////////////////////////////////////////////////////////////////////
// Default options

	template < bool DONT_VERIFY_ = false, bool DONT_DUMP_RESULTS_ = false >
	struct RazerSSpec 
	{
		enum { DONT_VERIFY = DONT_VERIFY_ };				// omit verifying potential matches
		enum { DONT_DUMP_RESULTS = DONT_DUMP_RESULTS_ };	// omit dumping results
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
        // TODO(holtgrew): SAM export should imply --read-naming 3
		unsigned	readNaming;			// 0..use Fasta id
										// 1..enumerate reads beginning with 1
										// 2..use the read sequence (only for short reads!)
										// 3..use Fasta id, do not append /L and /R for mate pairs.
		unsigned	sortOrder;			// 0..sort keys: 1. read number, 2. genome position
										// 1..           1. genome pos50ition, 2. read number
		unsigned	positionFormat;		// 0..gap space
										// 1..position space
		const char	*runID;				// runID needed for gff output
        bool        computeGlobal;      // compute global alignment in SAM output

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
        double      timeBuildQGramIndex;  // time for q-gram index building.
        double      timeCompactMatches;     // time for compacting reads
        double      timeMaskDuplicates; // time spent masking duplicates
        double      timeFsCopy; // time spent copying alignments back into the fragment store

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
	    double      noCompactFrac;      // If in last noCompactFrac of genome, don't compact.
        double      compactMult;        // Multiplicator for compaction threshold.
		unsigned	compactThresh;		// compact match array if larger than compactThresh

	// multi-threading

        unsigned    threadCount;  // Number of threads to use in the parallel version.
        unsigned    windowSize;  // Collect SWIFT hits in windows of this length.
        unsigned    verificationPackageSize;  // This number of SWIFT hits per verification.
        unsigned    maxVerificationPackageCount;  // Maximum number of verification packages to create.
        __int64     availableMatchesMemorySize;  // Memory available for matches.  Used for switching to external memory algorithms. -1 for always external, 0 for never.

#ifdef RAZERS_OPENADDRESSING
		double		loadFactor;
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
            computeGlobal = false;

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

            noCompactFrac = 0.05;
            compactMult = 2.2;
			compactThresh = 1024;
			// compactThresh = 40;

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

#ifdef _OPENMP
            threadCount = omp_get_max_threads();
#else  // #ifdef _OPENMP
            threadCount = 1;
#endif  // #ifdef _OPENMP
            // TODO(holtgrew): Tune this!
            windowSize = 500000;
            verificationPackageSize = 100;
            maxVerificationPackageCount = 100;
            availableMatchesMemorySize = 0;

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

	enum RAZERS_ERROR 
	{
		RAZERS_INVALID_OPTIONS = -1,
		RAZERS_READS_FAILED    = -2,
		RAZERS_GENOME_FAILED   = -3,
		RAZERS_INVALID_SHAPE   = -4
	};

//////////////////////////////////////////////////////////////////////////////
// Definitions

	template <typename TReadSet, typename TShape, typename TSpec>
	struct Cargo< Index<TReadSet, IndexQGram<TShape, TSpec> > > {
		typedef struct {
			double		abundanceCut;
			int			_debugLevel;
		} Type;
	};

//////////////////////////////////////////////////////////////////////////////
// Memory tuning

#ifdef RAZERS_MEMOPT

	template <typename TReadSet, typename TShape, typename TSpec>
	struct SAValue< Index<TReadSet, IndexQGram<TShape, TSpec> > > 
	{
		typedef Pair<
			unsigned,				
			unsigned,
			BitCompressed<24, 8>	// max. 16M reads of length < 256
		> Type;
	};
	
#else

	template <typename TReadSet, typename TShape, typename TSpec>
	struct SAValue< Index<TReadSet, IndexQGram<TShape, TSpec> > > 
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
	struct Size< Index<TReadSet, IndexQGram<TShape> > >
	{
		typedef unsigned Type;
	};

	template <typename TReadSet, typename TShape>
	struct Size< Index<TReadSet, IndexQGram<TShape, OpenAddressing> > >
	{
		typedef unsigned Type;
	};
	

#ifdef RAZERS_PRUNE_QGRAM_INDEX

	//////////////////////////////////////////////////////////////////////////////
	// Repeat masker
	template <typename TReadSet, typename TShape, typename TSpec>
	inline bool _qgramDisableBuckets(Index<TReadSet, IndexQGram<TShape, TSpec> > &index) 
	{
		typedef Index<TReadSet, IndexQGram<TShape, TSpec> >	TReadIndex;
		typedef typename Fibre<TReadIndex, QGramDir>::Type		TDir;
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
        typename TFragmentStore_,
		typename TMatches_, 
		typename TRazerSOptions_,
		typename TRazerSMode_,
		typename TSwiftPattern_,
		typename TCounts_,
        typename _TPreprocessing	
	>
	struct MatchVerifier
	{
		typedef TFragmentStore_									TFragmentStore;
		typedef TMatches_									    TMatches;
		typedef TRazerSOptions_									TOptions;
		typedef TRazerSMode_									TRazerSMode;
		typedef TSwiftPattern_									TSwiftPattern;
		typedef TCounts_										TCounts;
        typedef _TPreprocessing                                 TPreprocessing;
		
        typedef typename TRazerSMode::TMatchNPolicy             TMatchNPolicy;
		
		typedef typename TFragmentStore::TReadSeqStore			TReadSeqStore;
		typedef typename Value<TReadSeqStore>::Type const		TRead;
		typedef typename TFragmentStore::TAlignedReadStore		TAlignedReadStore;
		typedef typename TFragmentStore::TAlignQualityStore		TAlignQualityStore;
		typedef typename TFragmentStore::TContigSeq             TContigSeq;
		typedef typename Value<TAlignedReadStore>::Type			TAlignedRead;
		typedef typename Value<TAlignQualityStore>::Type		TAlignQuality;
		typedef typename Size<TContigSeq>::Type					TSize;
		typedef ModifiedString<TRead, ModReverse>				TRevRead;

        typedef typename Value<TMatches>::Type TMatchRecord;
		
#ifdef RAZERS_BANDED_MYERS
		typedef PatternState_<TRead,	Myers<AlignTextBanded<FindInfix, TMatchNPolicy, TMatchNPolicy>, True, void> > TPatternState;
		typedef PatternState_<TRevRead, Myers<AlignTextBanded<FindPrefix, TMatchNPolicy, TMatchNPolicy>, True, void> > TRPatternState;
#else  // #ifdef RAZERS_BANDED_MYERS
		typedef Pattern<TRead, Myers<FindInfix, False, void> >		TMyersPattern; 
		typedef Pattern<TRevRead, Myers<FindInfix, False, void> >	TRevMyersPattern; 
		typedef typename PatternState<TMyersPattern>::Type			TPatternState;
		typedef typename PatternState<TRevMyersPattern>::Type		TRPatternState;
#endif  // #ifdef RAZERS_BANDED_MYERS

		TMatches	    *matches;
		TOptions		*options;			// RazerS options
		TSwiftPattern	*swiftPattern;
		TCounts			*cnts;
		TPreprocessing  *preprocessing;

		TMatchRecord	m;
		bool			onReverseComplement;
		TSize			genomeLength;
		bool			oneMatchPerBucket;
		TPatternState	patternState;
		TRPatternState  revPatternState;

        double compactionTime;
		
		MatchVerifier() : onReverseComplement(false), genomeLength(0), oneMatchPerBucket(false), compactionTime(0) {}
               
		MatchVerifier(TMatches_ &_matches, TOptions &_options, TSwiftPattern &_swiftPattern, TCounts &_cnts):
			matches(&_matches),
			options(&_options),
			swiftPattern(&_swiftPattern),
			cnts(&_cnts),
            compactionTime(0)
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
                std::swap(m.beginPos, m.endPos);
                m.orientation = 'R';
			} else {
                m.orientation = 'F';
            }
						
//#pragma omp critical
// begin of critical section
			{
				if (!options->spec.DONT_DUMP_RESULTS)
				{
					appendValue(*matches, m, Generous());

					if (length(*matches) > options->compactThresh)
					{
                        double beginTime = sysTime();
						typename Size<TMatches>::Type oldSize = length(*matches);

						if (IsSameType<typename TRazerSMode::TGapMode, RazerSGapped>::VALUE)
							maskDuplicates(*matches, *options, TRazerSMode());	// overlapping parallelograms cause duplicates
                        // SEQAN_ASSERT_MSG((back(*matches).endPos - back(*matches).beginPos == 100), "len == %d", int(m.endPos - m.beginPos));
		
						compactMatches(*matches, *cnts, *options, TRazerSMode(), *swiftPattern, COMPACT);
                        // SEQAN_ASSERT_MSG((back(*matches).endPos - back(*matches).beginPos == 100), "len == %d", int(m.endPos - m.beginPos));
						
						if (length(*matches) * 4 > oldSize) {			// the threshold should not be raised
                            // fprintf(stderr, "[raising threshold]");
							// options->compactThresh += (options->compactThresh >> 1);	// if too many matches were removed
							options->compactThresh *= options->compactMult;
                        }
						
//						if (options._debugLevel >= 2)
//							::std::cerr << '(' << oldSize - length(store.alignedReadStore) << " matches removed)";
                        double endTime = sysTime();
                        compactionTime += (endTime - beginTime);
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
			_parseSkipWhitespace(file, c);
			appendValue(genomeFileNames,_parseReadFilepath(file,c));
			//CharString currentGenomeFile(filePrefix);
			//append(currentGenomeFile,_parseReadFilepath(file,c));
			//appendValue(genomeFileNames,currentGenomeFile);
			if(options._debugLevel >=2)
				std::cout <<"Genome file #"<< i <<": " << genomeFileNames[length(genomeFileNames)-1] << std::endl;
			++i;
			_parseSkipWhitespace(file, c);
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
	unsigned maxReadLength = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		if (options.readNaming == 0 || options.readNaming == 3
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
			if (options.readNaming == 0 || options.readNaming == 3)
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
		assignQualities(seq, qual); 
		if (options.trimLength > 0 && length(seq) > (unsigned)options.trimLength)
			resize(seq, options.trimLength);

		appendRead(store, seq, id);
		if (maxReadLength < length(seq))
			maxReadLength = length(seq);
	}
	// memory optimization
	reserve(store.readSeqStore.concat, length(store.readSeqStore.concat), Exact());
//	reserve(store.readNameStore.concat, length(store.readNameStore.concat), Exact());

	typedef Shape<Dna, SimpleShape> TShape;
	typedef typename SAValue< Index<StringSet<Dna5String>, IndexQGram<TShape, OpenAddressing> > >::Type TSAValue;
	TSAValue sa(0, 0);
	sa.i1 = ~sa.i1;
	sa.i2 = ~sa.i2;
	
	if ((unsigned)sa.i1 < length(store.readSeqStore) - 1)
	{
		::std::cerr << "Maximal read number of " << (unsigned)sa.i1 + 1 << " exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}
	if ((unsigned)sa.i2 < maxReadLength - 1)
	{
		::std::cerr << "Maximal read length of " << (unsigned)sa.i2 + 1 << " bps exceeded. Please remove \"#define RAZERS_MEMOPT\" in razers.cpp and recompile." << ::std::endl;
		seqCount = 0;
	}

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

// Comparators for RazerS1-style matches.

// TODO(holtgrew): Slightly different comparators than in previous RazerS 3 version, add back the additional checks?

	template <typename TReadMatch>
	struct LessRNoBeginPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// genome position and orientation
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;
			if (a.beginPos < b.beginPos) return true;
			if (a.beginPos > b.beginPos) return false;
            if (a.orientation == '-') return false;
            if (b.orientation == '-') return true;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// quality
            if (a.pairScore > b.pairScore) return true;
            if (a.pairScore < b.pairScore) return false;
            if (a.score > b.score) return true;
            if (b.score > a.score) return false;

            if (a.endPos > b.endPos) return true;
            return false;
		}
	};

	// ... to sort matches and remove duplicates with equal gEnd
	template <typename TReadMatch>
	struct LessRNoEndPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

			// genome position and orientation
			if (a.contigId < b.contigId) return true;
			if (a.contigId > b.contigId) return false;
			if (a.endPos   < b.endPos) return true;
			if (a.endPos   > b.endPos) return false;
            if (a.orientation == '-') return false;
            if (b.orientation == '-') return true;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// quality
            if (a.pairScore > b.pairScore) return true;
            if (a.pairScore < b.pairScore) return false;
            if (a.score > b.score) return true;
            if (b.score > a.score) return false;

            if (a.beginPos < b.beginPos) return true;
            return false;
		}
	};

	template <typename TReadMatch>
	struct LessScoreBackport : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;

            // quality
            if (a.orientation == '-') return false;
            if (b.orientation == '-') return true;

            if (a.pairScore > b.pairScore) return true;
            if (a.pairScore < b.pairScore) return false;
            if (a.score > b.score) return true;
            if (b.score > a.score) return false;

            // Sort by leftmost begin pos, longest end pos on ties.
            if (a.contigId < b.contigId) return true;
            if (a.contigId > b.contigId) return false;
            if (a.orientation < b.orientation) return true;
            if (a.orientation > b.orientation) return false;
            if (a.beginPos < b.beginPos) return true;
            if (a.beginPos > b.beginPos) return false;
            if (a.endPos < b.endPos) return false;
            if (a.endPos > b.endPos) return true;
            
            return false;
		}
	};

#ifdef RAZERS_DIRECT_MAQ_MAPPING

	template <typename TReadMatch>
	struct LessRNoMQ : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.readId < b.readId) return true;
			if (a.readId > b.readId) return false;
			
			// quality
			if (a.score < b.score) return true; // sum of quality values of mismatches (the smaller the better)
			if (a.score > b.score) return false;
			
			return (a.pairScore < b.pairScore); // seedEditDist?
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

// Comparators for Fragment Store

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
    template <typename TAlignedReadStore, typename TAlignedReadQualityStore, typename TGapMode, typename TScoreMode, typename TMatchNPolicy>
	struct LessScore<TAlignedReadStore, TAlignedReadQualityStore, RazerSMode<RazerSPrefix, TGapMode, TScoreMode, TMatchNPolicy> > : 
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
template <typename TMatches, typename TIterator, typename TOptions, typename TRazerSMode>
void maskDuplicates(TMatches &, TIterator const itBegin, TIterator const itEnd, TOptions & options, TRazerSMode)
{
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename TMatch::TContigPos			TContigPos;
	
	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal ends

    double beginTime = sysTime();
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
	::std::sort(itBegin, itEnd, LessRNoEndPos<TMatch>());
	// sortAlignedReads(matches, TLessEndPos(TLessScore(store.alignQualityStore)));
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

	TContigPos	beginPos = -1;
	TContigPos	endPos = -1;
	unsigned	contigId = TMatch::INVALID_ID;
	unsigned	readId = TMatch::INVALID_ID;
	char        orientation = '-';
    unsigned    masked = 0;

	TIterator it = itBegin;

	for (; it != itEnd; ++it) 
	{
		if ((*it).pairMatchId != TMatch::INVALID_ID) continue;
		TContigPos itEndPos = _max((*it).beginPos, (*it).endPos);
		if (endPos == itEndPos && orientation == (*it).orientation &&
			contigId == (*it).contigId && readId == (*it).readId) 
		{
			(*it).orientation = '-';
            masked += 1;
			continue;
		}
		readId = (*it).readId;
		contigId = (*it).contigId;
		endPos = itEndPos;
		orientation = (*it).orientation;
	}

	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal begins

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
    ::std::sort(itBegin, itEnd, LessRNoBeginPos<TMatch>());
	// sortAlignedReads(store.alignedReadStore, TLessBeginPos(TLessScore(store.alignQualityStore)));
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

    orientation = '-';
	contigId = TMatch::INVALID_ID;
	it = itBegin;

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-' || (*it).pairMatchId != TMatch::INVALID_ID) continue;

		TContigPos itBeginPos = _min((*it).beginPos, (*it).endPos);
		if (beginPos == itBeginPos && readId == (*it).readId &&
			contigId == (*it).contigId && orientation == ((*it).beginPos < (*it).endPos))
		{
			(*it).orientation = '-';
            masked += 1;
			continue;
		}
		readId = (*it).readId;
		contigId = (*it).contigId;
		beginPos = itBeginPos;
		orientation = (*it).beginPos < (*it).endPos;
	}
    options.timeMaskDuplicates += sysTime() - beginTime;
    if (options._debugLevel >= 2)
        fprintf(stderr, " [%u matches masked]", masked);
}

template <typename TMatches, typename TOptions, typename TRazerSMode>
void maskDuplicates(TMatches &matches, TOptions & options, TRazerSMode const & mode)
{
    maskDuplicates(matches, begin(matches, Standard()), end(matches, Standard()), options, mode);
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
void countMatches(TFragmentStore &store, TCounts &cnt, TRazerSMode const &)
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
	__int64 maxVal = MaxValue<TValue>::VALUE;

#ifdef RAZERS_PROFILE
  timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
	sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));
#ifdef RAZERS_PROFILE
  timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

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
	// if (readNo==643)
	// 	std::cout<<"dman"<<std::endl;
	int minT = _qgramLemma(swift, readNo, maxErrors);
	if (minT > 1)
	{
//		::std::cout<<" read:"<<readNo<<" newThresh:"<<minT;
		if (maxErrors < 0) minT = MaxValue<int>::VALUE;
		setMinThreshold(swift, readNo, (unsigned)minT);
	}
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template <
	typename TMatches,
	typename TCounts,
	typename TSpec,
	typename TAlignMode,
	typename TGapMode,
	typename TScoreMode,
	typename TSwift,
    typename TMatchNPolicy
>
void compactMatches(
	TMatches & matches,
	TCounts &, 
	RazerSOptions<TSpec> &options, 
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const &,
	TSwift & swift, 
	CompactMatchesMode compactMode)
{
    // fprintf(stderr, "[compact]");
    double beginTime = sysTime();
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type	TIterator;
	typedef RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> TRazerSMode;

	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int scoreCutOff = MinValue<int>::VALUE;
	int errorCutOff = MaxValue<int>::VALUE;
	int errorRangeBest = options.errorDistanceRange;// (IsSameType<TScoreMode, RazerSErrors>::VALUE)? options.scoreDistanceRange: 0;
	int scoreRangeBest = (IsSameType<TAlignMode, RazerSGlobal>::VALUE && !IsSameType<TScoreMode, RazerSScore>::VALUE)? -(int)options.scoreDistanceRange : MaxValue<int>::VALUE;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
#ifdef RAZERS_EXTERNAL_MATCHES
    if (compactMode == COMPACT_FINAL_EXTERNAL) {
        typedef Pipe<TMatches, Source<> > TSource;
        typedef LessScoreBackport<TMatch> TLess;
        typedef Pool<TMatch, SorterSpec<SorterConfigSize<TLess, typename Size<TSource>::Type> > > TSorterPool;

        TSource source(matches);
        TSorterPool sorter;
        sorter << matches;
        beginRead(sorter);
        SEQAN_ASSERT_EQ(length(sorter), length(matches));
        TIterator it = begin(matches, Standard());
        TIterator itEnd = end(matches, Standard());
        (void) itEnd;
        // bool first = true;
        for (__int64 leftToRead = length(sorter); leftToRead > 0; --leftToRead, ++sorter, ++it) {
            *it = *sorter;
            // if (!first)
            //     SEQAN_ASSERT(!TLess()(*it, *(it - 1)));
            // first = false;
        }
        SEQAN_ASSERT(it == itEnd);
        endRead(sorter);
    } else {
#endif  // #ifdef RAZERS_EXTERNAL_MATCHES
        ::std::sort(begin(matches, Standard()), end(matches, Standard()), LessScoreBackport<TMatch>());
        // sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));
#ifdef RAZERS_EXTERNAL_MATCHES
    }
#endif  // #ifdef RAZERS_EXTERNAL_MATCHES
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;
	TIterator ditBeg = it;
    // fprintf(stderr, "[%u matches to compact]", unsigned(itEnd - it));
    unsigned disabled = 0;
	
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		int score = (*it).score;
		int errors = -(*it).score;
        //if (readNo == 643) std::cerr <<"["<<score<<","<<errors<<"] "<<::std::flush;
		if (readNo == (*it).readId && (*it).pairMatchId == TMatch::INVALID_ID)  // Only compact unpaired matches.
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
                        // std::cerr << "PURGED " << readNo << std::endl;
						setMaxErrors(swift, readNo, -1);
                        disabled += 1;
						// if (options._debugLevel >= 2)
						// 	::std::cerr << "(read #" << readNo << " disabled)";
					}
					else
						// we only need better matches
						if (IsSameType<TScoreMode, RazerSErrors>::VALUE)
						{
                            // std::cerr << "LIMITED " << readNo << std::endl;
							setMaxErrors(swift, readNo, errors - 1);
                            disabled += 1;
							// if (errors == 0 && options._debugLevel >= 2)
							// 	::std::cerr << "(read #" << readNo << " disabled)";
						}

					if (options.purgeAmbiguous && (compactMode == COMPACT_FINAL || compactMode == COMPACT_FINAL_EXTERNAL))
					{
						if (options.scoreDistanceRange == 0 || errors < errorRangeBest || score > scoreRangeBest || compactMode == COMPACT_FINAL || compactMode == COMPACT_FINAL_EXTERNAL){
							dit = ditBeg;
						}
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
	unsigned origSize = length(matches);
	resize(matches, dit - begin(matches, Standard()));
	// compactAlignedReads(store);
    options.timeCompactMatches += sysTime() - beginTime;
    // fprintf(stderr, "[compacted in %f s]", endTime - beginTime);
	unsigned newSize = length(matches);
    if (options._debugLevel >= 2) {
        fprintf(stderr, " [%u matches removed]", unsigned(origSize - newSize));
        fprintf(stderr, " [%u reads disabled]", disabled);
    }
}

//////////////////////////////////////////////////////////////////////////////
// Remove low quality matches
template <
	typename TMatches,
	typename TCounts,
	typename TSpec,
	typename TGapMode,
	typename TSwift,
    typename TMatchNPolicy
>
void compactMatches(
	TMatches & matches,
	TCounts & cnts, 
	RazerSOptions<TSpec> &options, 
	RazerSMode<RazerSGlobal, TGapMode, RazerSQuality<RazerSMAQ>, TMatchNPolicy> const &,
	TSwift & swift, 
	CompactMatchesMode compactMode)
{
	typedef typename Value<TMatches>::Type						TMatch;
	typedef typename Iterator<TMatches, Standard>::Type			TIterator;
	typedef RazerSMode<RazerSGlobal, TGapMode, RazerSQuality<RazerSMAQ>, TMatchNPolicy> TRazerSMode;
	
	unsigned readNo = -1;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessScoreBackport<TMatch>());
	// sortAlignedReads(store.alignedReadStore, LessScore<TAlignedReadStore, TAlignQualityStore, TRazerSMode>(store.alignQualityStore));
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_SORT);
#endif  // #ifdef RAZERS_PROFILE
		
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
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

	resize(matches, dit - begin(matches, Standard()));
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
	typename TReadSet,
    typename TMatchNPolicy >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,									// potential match genome region
	unsigned readId,													// read number
	TReadSet &readSet,													// reads
	RazerSMode<RazerSPrefix, RazerSUngapped, RazerSErrors, TMatchNPolicy> const &)	// Hamming only
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
	unsigned errorThresh = (verifier.oneMatchPerBucket)? MaxValue<unsigned>::VALUE: maxErrors;
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
struct UseQualityValues__ { enum { VALUE = false }; };
template <typename TAlignMode, typename TGapMode, typename TSpec, typename TMatchNPolicy>
struct UseQualityValues__<RazerSMode<TAlignMode, TGapMode, RazerSQuality<TSpec>, TMatchNPolicy> > {  enum { VALUE = true }; };

//////////////////////////////////////////////////////////////////////////////
// Hamming verification
template <
	typename TMatchVerifier,
	typename TGenome, 
	typename TReadSet,
	typename TScoreMode,
    typename TMatchNPolicy >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,								// potential match genome region
	unsigned readId,												// read number
	TReadSet &readSet,												// reads
	RazerSMode<RazerSGlobal, RazerSUngapped, TScoreMode, TMatchNPolicy> const &) // Semi-global, no gaps
{
	typedef Segment<TGenome, InfixSegment>							TGenomeInfix;
	typedef typename Value<TReadSet>::Type const					TRead;
	typedef typename Iterator<TGenomeInfix, Standard>::Type			TGenomeIterator;
	typedef typename Iterator<TRead, Standard>::Type				TReadIterator;
	typedef RazerSMode<RazerSGlobal, RazerSUngapped, TScoreMode, TMatchNPolicy> TRazerSMode;

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
	if (IsSameType<TScoreMode, RazerSErrors>::VALUE)
		minScore = -(int)(ndlLength * verifier.options->errorRate);
	else if (UseQualityValues__<TRazerSMode>::VALUE)
		minScore = -verifier.options->absMaxQualSumErrors;
	else if (IsSameType<TScoreMode, RazerSScore>::VALUE)
	{
		minScore = verifier.options->minScore;
		mismatchDelta = scoreMatch(verifier.options->scoringScheme) - scoreMismatch(verifier.options->scoringScheme);
		scoreInit = scoreMatch(verifier.options->scoringScheme) * ndlLength;
	}
	
	int maxScore = minScore - 1;
	int scoreThresh = (verifier.oneMatchPerBucket)? MaxValue<int>::VALUE: minScore;
	int score, errors;
	
	for (; git < gitEnd; ++git)
	{
		if (!IsSameType<TScoreMode, RazerSScore>::VALUE)
			score = 0;
		else
			score = scoreInit;
		
		if (!IsSameType<TScoreMode, RazerSErrors>::VALUE)
			errors = 0;
		
		TGenomeIterator g = git;
		for (TReadIterator r = ritBeg; r != ritEnd; ++r, ++g)
			if ((verifier.options->compMask[ordValue(*g)] & verifier.options->compMask[ordValue(*r)]) == 0)
			{
				if (IsSameType<TScoreMode, RazerSErrors>::VALUE)
				{
					// A. Count mismatches only
					--score;
				} else
				{
					++errors;
					if (UseQualityValues__<TRazerSMode>::VALUE)
						// B. Count mismatches and mismatch qualities
						score -= getQualityValue(*g);
					else if (IsSameType<TScoreMode, RazerSScore>::VALUE)
						// C. Count mismatches and alignment score
						score -= mismatchDelta;
					else
						SEQAN_FAIL("Unsupported score mode!");
				}
				if (score < minScore)	// doesn't work for islands with errorThresh > maxErrors
					break;
			}
		
		if (score > maxScore)
		{
			maxScore = score;
			if (IsSameType<TScoreMode, RazerSErrors>::VALUE)
				verifier.m.score = score;
			else
				verifier.m.score = errors;
			verifier.m.beginPos = git - begin(host(inf), Standard());
		}
#ifdef RAZERS_ISLAND_CRITERION
        else if (scoreThresh > score)
		{
			if (maxScore >= minScore)
			{
				// for RazerSErrors bestErrors == -maxScore
				verifier.m.endPos = verifier.m.beginPos + ndlLength;
				verifier.m.pairScore = verifier.m.score = maxScore;
				if (!verifier.oneMatchPerBucket)
					verifier.push();
				maxScore = minScore - 1;
			}
		}
#else
        (void)scoreThresh;
#endif
	}

	if (maxScore >= minScore)
	{
		verifier.m.endPos = verifier.m.beginPos + ndlLength;
		verifier.m.pairScore = verifier.m.score = maxScore;
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
	typename TReadSet,
    typename TMatchNPolicy >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,									// potential match genome region
	unsigned readId,														// read number
	TReadSet &readSet,													// reads
	RazerSMode<RazerSGlobal, RazerSGapped, RazerSErrors, TMatchNPolicy> const &) // Mismatches and Indels
{
    if (length(inf) == 0u)
      return false;

	typedef Segment<TGenome, InfixSegment>					TGenomeInfix;
	typedef typename Value<TReadSet>::Type const			TRead;
	typedef typename Position<TGenomeInfix>::Type			TPosition;

	// find read match end
	typedef Finder<TGenomeInfix>							TMyersFinder;
	typedef typename TMatchVerifier::TPatternState			TPatternState;
    typedef typename TMatchVerifier::TPreprocessing         TPreprocessing; 
	typedef typename Value<TPreprocessing>::Type			TMyersPattern;

	// find read match begin
    // TODO(holtgrew): Use reverse-search here, as well!
	typedef ModifiedString<TGenomeInfix, ModReverse>		TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>				TReadRev;
	typedef Finder<TGenomeInfixRev>							TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>			TMyersPatternRev;

	TMyersFinder myersFinder(inf);
#ifndef RAZERS_BANDED_MYERS
    TMyersPattern &myersPattern = (*verifier.preprocessing)[readId]; 
#endif  // #ifdef RAZERS_BANDED_MYERS
	TPatternState & state = verifier.patternState;
	
#ifdef RAZERS_DEBUG
	::std::cout<<"Verify: "<<::std::endl;
	::std::cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << ::std::endl;
	::std::cout<<"Read:   "<<readSet[readId] << "(id: " << readId << ")" <<::std::endl;
#endif

    unsigned ndlLength = sequenceLength(readId, readSet);
	int maxScore = MinValue<int>::VALUE;
	int minScore = -(int)(ndlLength * verifier.options->errorRate);
	TPosition maxPos = 0;
	TPosition lastPos = length(inf);
#ifdef RAZERS_ISLAND_CRITERION
	unsigned minDistance = (verifier.oneMatchPerBucket)? lastPos: 1;
#else  // #ifdef RAZERS_ISLAND_CRITERION
	unsigned minDistance = lastPos;
    (void)minDistance;
#endif  // #ifdef RAZERS_ISLAND_CRITERION

	// find end of best semi-global alignment
#ifdef RAZERS_BANDED_MYERS
    TRead read(readSet[readId]);  // here only infixes (no sequence) is copied
	while (find(myersFinder, read, state, minScore))
#else  // #ifdef RAZERS_BANDED_MYERS
    while (find(myersFinder, myersPattern, state, minScore)) 
#endif  // #ifdef RAZERS_BANDED_MYERS
	{
		TPosition const pos = position(hostIterator(myersFinder));
#ifdef RAZERS_ISLAND_CRITERION
		if (lastPos + minDistance < pos)
		{
			if (minScore <= maxScore)
			{
				verifier.m.endPos = beginPosition(inf) + maxPos + 1;
				verifier.m.score = -maxScore;

				verifier.m.pairScore = verifier.m.score =  maxScore;
				if (maxScore == 0)
					verifier.m.beginPos = verifier.m.endPos - ndlLength;
				else
				{
					TPosition infBeginPos = beginPosition(inf);
					TPosition infEndPos = endPosition(inf);

					// find beginning of best semi-global alignment
					TGenomeInfixRev infRev(inf);
					
					TPosition newInfEndPos = infBeginPos + maxPos + 1;
#ifdef RAZERS_BANDED_MYERS
					verifier.revPatternState.leftClip = infEndPos - newInfEndPos;
#endif
					setEndPosition(inf, verifier.m.endPos = newInfEndPos);

					// limit the beginning to needle length plus errors (== -maxScore)
					if (length(inf) > ndlLength - maxScore)
						setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
					
					TReadRev			readRev(readSet[readId]);
					TMyersFinderRev		myersFinderRev(infRev);

#ifdef RAZERS_BANDED_MYERS
					verifier.m.beginPos = verifier.m.endPos;
					while (find(myersFinderRev, readRev, verifier.revPatternState, maxScore)) {
#else
					TMyersPatternRev	myersPatternRev(readRev);

					_patternMatchNOfPattern(myersPatternRev, verifier.options->matchN);
					_patternMatchNOfFinder(myersPatternRev, verifier.options->matchN);
					while (find(myersFinderRev, myersPatternRev, maxScore)) {
#endif
						verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);
                    }

					setBeginPosition(inf, infBeginPos);
					setEndPosition(inf, infEndPos);
				}
				// minDistance implicitly forbids to get here with verifier.oneMatchPerBucket == true
				SEQAN_ASSERT_NOT(verifier.oneMatchPerBucket);
#ifdef RAZERS_BANDED_MYERS
				if (verifier.m.beginPos != verifier.m.endPos)
					verifier.push();
#else  // RAZERS_BANDED_MYERS
				SEQAN_ASSERT_LT(verifier.m.beginPos, verifier.m.endPos);
				verifier.push();
#endif // RAZERS_BANDED_MYERS

				maxScore = minScore - 1;
			}
		}
#endif  // #ifdef RAZERS_ISLAND_CRITERION
		// if (getScore(myersPattern) >= maxScore)
		if (getScore(state) >= maxScore)
		{
			maxScore = getScore(state);
			// maxScore = getScore(myersPattern);
			maxPos = pos;
		}
		lastPos = pos;
	}
	
	if (minScore <= maxScore)
	{
		verifier.m.endPos = beginPosition(inf) + maxPos + 1;
		verifier.m.score = maxScore;

		verifier.m.pairScore = verifier.m.score = maxScore;
		if (maxScore == 0)
			verifier.m.beginPos = verifier.m.endPos - ndlLength;
		else
		{
            TPosition newInfEndPos = beginPosition(inf) + maxPos + 1;
#ifdef RAZERS_BANDED_MYERS
            verifier.revPatternState.leftClip = endPosition(inf) - newInfEndPos;
#endif
            setEndPosition(inf, verifier.m.endPos = newInfEndPos);

			// limit the beginning to needle length plus errors (== -maxScore)
			if (length(inf) > ndlLength - maxScore)
				setBeginPosition(inf, endPosition(inf) - ndlLength + maxScore);
			
			// find beginning of best semi-global alignment
			TGenomeInfixRev		infRev(inf);
			TReadRev			readRev(readSet[readId]);
			TMyersFinderRev		myersFinderRev(infRev);
			TMyersPatternRev	myersPatternRev(readRev);

#ifdef RAZERS_BANDED_MYERS
			verifier.m.beginPos = verifier.m.endPos;
            while (find(myersFinderRev, readRev, verifier.revPatternState, maxScore)) {
#else
            _patternMatchNOfPattern(myersPatternRev, verifier.options->matchN);
            _patternMatchNOfFinder(myersPatternRev, verifier.options->matchN);
            while (find(myersFinderRev, myersPatternRev, maxScore)) {
#endif
				verifier.m.beginPos = verifier.m.endPos - (position(myersFinderRev) + 1);
            }
#ifdef RAZERS_BANDED_MYERS
			if (verifier.m.beginPos == verifier.m.endPos)
			{
#ifdef RAZERS_DEBUG
			    std::cout << "FAILED" << std::endl; 
#endif
				return false;
			}
#else  // RAZERS_BANDED_MYERS
			SEQAN_ASSERT_LT(verifier.m.beginPos, verifier.m.endPos);
#endif // RAZERS_BANDED_MYERS
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
	typename TReadSet,
    typename TMatchNPolicy >
inline bool
matchVerify(
	TMatchVerifier &verifier,
	Segment<TGenome, InfixSegment> inf,												// potential match genome region
	unsigned readId,																	// read number
	TReadSet &readSet,																// reads
	RazerSMode<RazerSGlobal, RazerSUngapped, RazerSQuality<RazerSMAQ>, TMatchNPolicy> const &)	// Hamming only
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
	unsigned qualSumThresh = (verifier.oneMatchPerBucket)? MaxValue<unsigned>::VALUE: verifier.options->absMaxQualSumErrors;
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
				verifier.m.pairScore = verifier.m.score = - (int)minQualSum;
				verifier.m.errors = minErrors;
				verifier.push();
				minQualSum = verifier.options->absMaxQualSumErrors + 1;
			}
		}
	}

	if (minErrors <= maxErrors)
	{
		verifier.m.endPos = verifier.m.beginPos + ndlLength;
		verifier.m.pairScore = verifier.m.score = -(int)minQualSum;
		verifier.m.errors = minErrors;
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
	typename TScoreMode,
    typename TMatchNPolicy >
inline bool
matchVerify(
	TMatchVerifier &,
	Segment<TGenome, InfixSegment>,								// potential match genome region
	unsigned,													// read number
	TReadSet &,													// reads
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const &)
{
    SEQAN_FAIL("Verification not implemented!");
	return false;
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches in a single genome sequence
template <
    typename TMatches,
	typename TFragmentStore, 
	typename TReadIndex, 
	typename TSwiftSpec, 
	typename TCounts,
	typename TRazerSOptions,
	typename TRazerSMode,
	typename TPreprocessing>
void _mapSingleReadsToContig(
    TMatches                                & matches,
	TFragmentStore							& store,
	unsigned								  contigId,				// ... and its sequence number
	Pattern<TReadIndex, Swift<TSwiftSpec> >	& swiftPattern,
	TCounts									& cnts,
	char									  orientation,				// q-gram index of reads
	TRazerSOptions							& options,
	TRazerSMode						  const & mode,
#ifdef RAZERS_BANDED_MYERS
	 TPreprocessing							&)
#else  // #ifdef RAZERS_BANDED_MYERS
	 TPreprocessing							& preprocessing)
#endif  // #ifdef RAZERS_BANDED_MYERS
{
	// FILTRATION
	typedef typename TFragmentStore::TContigSeq				TContigSeq;
	typedef Finder<TContigSeq, Swift<TSwiftSpec> >			TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >			TSwiftPattern;
	
	// VERIFICATION
	typedef MatchVerifier <
        TFragmentStore,
		TMatches, 
		TRazerSOptions, 
		TRazerSMode,
		TSwiftPattern,
		TCounts,
		TPreprocessing>											TVerifier;
	typedef typename Fibre<TReadIndex, FibreText>::Type	TReadSet;
	
	// iterate all genomic sequences
	if (options._debugLevel >= 1)
	{
		::std::cerr << ::std::endl << "Process genome seq #" << contigId;
		if (orientation == 'F') ::std::cerr << "[fwd]";
		else                    ::std::cerr << "[rev]";
	}
	lockContig(store, contigId);
	TContigSeq &contigSeq = store.contigStore[contigId].seq;
	if (orientation == 'R')	reverseComplement(contigSeq);

	TReadSet		&readSet = host(host(swiftPattern));
	TSwiftFinder	swiftFinder(contigSeq, options.repeatLength, 1);
	TVerifier		verifier(matches, options, swiftPattern, cnts);

#ifndef RAZERS_BANDED_MYERS
	verifier.preprocessing = &preprocessing;
#endif  // #ifdef RAZERS_BANDED_MYERS

	// initialize verifier
	verifier.onReverseComplement = (orientation == 'R');
	verifier.genomeLength = length(contigSeq);
	verifier.m.contigId = contigId;

    double beginTime = sysTime();
    // Build q-gram index separately, so we can better compute the time for it.
    indexRequire(host(swiftPattern), QGramSADir());
    options.timeBuildQGramIndex += sysTime() - beginTime;

	// iterate all verification regions returned by SWIFT
	while (find(swiftFinder, swiftPattern, options.errorRate))
	{
        // std::cout << "read id = " << (*swiftFinder.curHit).ndlSeqNo << ", " << beginPosition(swiftFinder) << std::endl;
        
//        if (length(infix(swiftFinder)) < length(readSet[(*swiftFinder.curHit).ndlSeqNo]))
//            continue;  // Skip if hit length < read length.  TODO(holtgrew): David has to fix something in banded myers to make this work.
#ifdef RAZERS_BANDED_MYERS
		verifier.patternState.leftClip = (beginPosition(swiftFinder) >= 0)? 0: -beginPosition(swiftFinder);	// left clip if match begins left of the genome
#endif
		verifier.m.readId = (*swiftFinder.curHit).ndlSeqNo;
		if (!options.spec.DONT_VERIFY)
			matchVerify(verifier, infix(swiftFinder), verifier.m.readId, readSet, mode);
		++options.countFiltration;
	}
	if (!unlockAndFreeContig(store, contigId))							// if the contig is still used
		if (orientation == 'R')	reverseComplement(contigSeq);	// we have to restore original orientation
}


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
	typename TReadIndex,
    typename TMatchNPolicy >
int _mapSingleReads(
	FragmentStore<TFSSpec, TFSConfig>					& store,
	TCounts												& cnts,
	RazerSOptions<TSpec>								& options,
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy>  const & mode,
	TReadIndex											& readIndex)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
	typedef typename If<
				IsSameType<TGapMode,RazerSGapped>::VALUE,
				SwiftSemiGlobal,
				SwiftSemiGlobalHamming>::Type			TSwiftSpec;
	typedef Pattern<TReadIndex, Swift<TSwiftSpec> >		TSwiftPattern;	// filter
	
	typedef typename Value<TReadSeqStore>::Type	const	TRead;
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;	// verifier
	// typedef Pattern<TRead, Myers<FindInfix, False, void> >	TMyersPattern;	// verifier

    typedef typename TFragmentStore::TContigSeq TContigSeq;
    typedef typename Position<TContigSeq>::Type TContigPos;
    typedef MatchRecord<TContigPos> TMatchRecord;
	
	// configure Swift pattern
	TSwiftPattern swiftPattern(readIndex);
	swiftPattern.params.minThreshold = options.threshold;
	swiftPattern.params.tabooLength = options.tabooLength;
	swiftPattern.params.printDots = options._debugLevel > 0; 

	// init edit distance verifiers
	options.compMask[4] = (options.matchN)? 15: 0;
#ifdef RAZERS_BANDED_MYERS
	Nothing preprocessing;
#else  // #ifdef RAZERS_BANDED_MYERS
    String<TMyersPattern> preprocessing;
	if (options.gapMode == RAZERS_GAPPED) 
	{
		unsigned readCount = countSequences(readIndex);
		resize(preprocessing, readCount, Exact()); 
		for(unsigned i = 0; i < readCount; ++i) 
		{ 
			setHost(preprocessing[i], indexText(readIndex)[i]); 
			_patternMatchNOfPattern(preprocessing[i], options.matchN); 
			_patternMatchNOfFinder(preprocessing[i], options.matchN); 
		} 
	}
#endif  // #ifdef RAZERS_BANDED_MYERS
	
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	unsigned readCount = countSequences(readIndex);
	if (options.maqMapping)
	{
		resize(cnts, 2);
		for (unsigned i = 0; i < length(cnts); ++i)
			resize(cnts[i], readCount, 31); //initialize with maxeditDist, 11:5 for count:dist
	}
#endif  // #ifdef RAZERS_DIRECT_MAQ_MAPPING

	// clear stats
	options.countFiltration = 0;
	options.countVerification = 0;
	options.timeMapReads = 0;
	options.timeDumpResults = 0;
    options.timeBuildQGramIndex = 0;
    options.timeCompactMatches = 0;
    options.timeMaskDuplicates = 0;
    options.timeFsCopy = 0;
	SEQAN_PROTIMESTART(find_time);

    // We collect the matches in a more compact data structure than the
    // AlignedReadStoreElement from FragmentStore.
    String<TMatchRecord> matches;
	
	// iterate over genome sequences
    for (unsigned contigId = 0; contigId < length(store.contigStore); ++contigId)
	{
		// lock to prevent releasing and loading the same contig twice
		// (once per _mapSingleReadsToContig call)
		lockContig(store, contigId);
#ifndef RAZERS_WINDOW
		if (options.forward)
			_mapSingleReadsToContig(matches, store, contigId, swiftPattern, cnts, 'F', options, mode, preprocessing);
		if (options.reverse)
			_mapSingleReadsToContig(matches, store, contigId, swiftPattern, cnts, 'R', options, mode, preprocessing);
#else
		if (options.forward)
			_mapSingleReadsToContigWindow(store, contigId, swiftPattern, cnts, 'F', options, mode, preprocessing);
		if (options.reverse)
			_mapSingleReadsToContigWindow(store, contigId, swiftPattern, cnts, 'R', options, mode, preprocessing);
#endif
		unlockAndFreeContig(store, contigId);
	}

	options.timeMapReads = SEQAN_PROTIMEDIFF(find_time);
	if (options._debugLevel >= 1)
		::std::cerr << ::std::endl << "Finding reads took               \t" << options.timeMapReads << " seconds" << ::std::endl;
    
    double beginCopyTime = sysTime();
    // Final mask duplicates and compact matches.
    Nothing nothing;
    if (IsSameType<TGapMode, RazerSGapped>::VALUE)
        maskDuplicates(matches, options, mode);
    compactMatches(matches, cnts, options, mode, nothing, COMPACT_FINAL);
    // Write back to store.
    reserve(store.alignedReadStore, length(matches));
    reserve(store.alignQualityStore, length(matches));
    typedef typename Iterator<String<TMatchRecord>, Standard>::Type TIterator;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedReadStoreElem;
	typedef typename Value<typename TFragmentStore::TAlignQualityStore>::Type TAlignedQualStoreElem;
    for (TIterator it = begin(matches), itEnd = end(matches); it != itEnd; ++it) {
        SEQAN_ASSERT_NEQ(it->orientation, '-');
        SEQAN_ASSERT_LEQ(it->beginPos, it->endPos);
        if (it->orientation == 'R')
            ::std::swap(it->beginPos, it->endPos);
        appendValue(store.alignedReadStore, TAlignedReadStoreElem(length(store.alignQualityStore), it->readId, it->contigId, it->beginPos, it->endPos));
        appendValue(store.alignQualityStore, TAlignedQualStoreElem(it->pairScore, it->score, -it->score));
    }
    options.timeFsCopy = sysTime() - beginCopyTime;

	if (options._debugLevel >= 2) {
		::std::cerr << "Masking duplicates took          \t" << options.timeMaskDuplicates << " seconds" << ::std::endl;
		::std::cerr << "Compacting matches took          \t" << options.timeCompactMatches << " seconds" << ::std::endl;
		::std::cerr << "Building q-gram index took       \t" << options.timeBuildQGramIndex << " seconds" << ::std::endl;
		::std::cerr << "Time for copying back            \t" << options.timeFsCopy << " seconds" << ::std::endl;
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
	typename TScoreMode,
    typename TMatchNPolicy >
int _mapSingleReads(
	FragmentStore<TFSSpec, TFSConfig>					& store,
	TCounts												& cnts,
	RazerSOptions<TSpec>								& options,
	TShape const										& shape,
	RazerSMode<TAlignMode, TGapMode, TScoreMode, TMatchNPolicy> const & mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>			TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore		TReadSeqStore;
#ifndef RAZERS_OPENADDRESSING
	typedef Index<TReadSeqStore, IndexQGram<TShape> >	TIndex;			// q-gram index
#else
	typedef Index<TReadSeqStore, IndexQGram<TShape, OpenAddressing> >	TIndex;
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
	typename TScoreMode,
    typename TMatchNPolicy >
int _mapSingleReads(
	FragmentStore<TFSSpec, TFSConfig>						& store,
	TCounts													& cnts,
	RazerSOptions<TSpec>									& options,
	TShape const											& shape,
	RazerSMode<RazerSPrefix, TGapMode, TScoreMode, TMatchNPolicy> const & mode)
{
	typedef FragmentStore<TFSSpec, TFSConfig>				TFragmentStore;
	typedef typename TFragmentStore::TReadSeqStore			TReadSeqStore;
	
	typedef typename Value<TReadSeqStore>::Type				TRead;
	typedef typename Infix<TRead>::Type						TReadInfix;
	typedef StringSet<TReadInfix>							TReadSet;
	typedef Index<TReadSet, IndexQGram<TShape> >			TIndex;			// q-gram index

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
    if (options.threadCount == 0 || length(store.readNameStore) < MIN_PARALLEL_WORK) {
        // Sequential RazerS
        #ifdef RAZERS_MATEPAIRS
        if (options.libraryLength >= 0)
            return _mapMatePairReads(store, cnts, options, shape, mode);
        else
        #endif  // #ifndef RAZERS_MATEPAIRS
            return _mapSingleReads(store, cnts, options, shape, mode);
    } else {
        // Parallel RazerS
        #ifdef RAZERS_MATEPAIRS
        if (options.libraryLength >= 0)
            return _mapMatePairReadsParallel(store, cnts, options, shape, mode);
        else
        #endif  // #ifndef RAZERS_MATEPAIRS
            return _mapSingleReadsParallel(store, cnts, options, shape, mode);
    }
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
template <typename TFSSpec, typename TFSConfig, typename TCounts, typename TSpec, typename TAlignMode, typename TGapMode, typename TMatchNPolicy>
int _mapReads(
	FragmentStore<TFSSpec, TFSConfig>		& store,
	TCounts									& cnts,
	RazerSOptions<TSpec>					& options,
	RazerSMode<TAlignMode, TGapMode, Nothing, TMatchNPolicy> const)
{
	if (options.scoreMode == RAZERS_ERRORS)
		return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSErrors, TMatchNPolicy>());
	if (options.scoreMode == RAZERS_SCORE)
		return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSScore, TMatchNPolicy>());
	if (options.scoreMode == RAZERS_QUALITY)
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		if (options.maqMapping)
			return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSQuality<RazerSMAQ>, TMatchNPolicy>());
		else
#endif
			return _mapReads(store, cnts, options, RazerSMode<TAlignMode, TGapMode, RazerSQuality<>, TMatchNPolicy>());
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
//    if (options.matchN) {
//        if (options.gapMode == RAZERS_GAPPED)
//         {
//             if (options.alignMode == RAZERS_LOCAL)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSGapped, Nothing, NMatchesAll_>());
//             if (options.alignMode == RAZERS_PREFIX)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSGapped, Nothing, NMatchesAll_>());
//             if (options.alignMode == RAZERS_GLOBAL)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSGapped, Nothing, NMatchesAll_>());
//         } else {
//             if (options.alignMode == RAZERS_LOCAL)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSUngapped, Nothing, NMatchesAll_>());
//             if (options.alignMode == RAZERS_PREFIX)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSUngapped, Nothing, NMatchesAll_>());
//             if (options.alignMode == RAZERS_GLOBAL)
//                 return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSUngapped, Nothing, NMatchesAll_>());
//         }
//    } else {
        if (options.gapMode == RAZERS_GAPPED)
        {
        //     if (options.alignMode == RAZERS_LOCAL)
        //         return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSGapped, Nothing, NMatchesNone_>());
        //     if (options.alignMode == RAZERS_PREFIX)
        //         return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSGapped, Nothing, NMatchesNone_>());
           if (options.alignMode == RAZERS_GLOBAL)
               return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSGapped, Nothing, NMatchesNone_>());
        } else {
             // if (options.alignMode == RAZERS_LOCAL)
             //     return _mapReads(store, cnts, options, RazerSMode<RazerSLocal, RazerSUngapped, Nothing, NMatchesNone_>());
             // if (options.alignMode == RAZERS_PREFIX)
             //     return _mapReads(store, cnts, options, RazerSMode<RazerSPrefix, RazerSUngapped, Nothing, NMatchesNone_>());
             if (options.alignMode == RAZERS_GLOBAL)
                 return _mapReads(store, cnts, options, RazerSMode<RazerSGlobal, RazerSUngapped, Nothing, NMatchesNone_>());
        }
//    }
	return RAZERS_INVALID_OPTIONS;
}

}

#endif
