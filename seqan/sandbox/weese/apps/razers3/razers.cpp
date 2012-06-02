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

#define SEQAN_PROFILE					// enable time measuring
//#define SEQAN_DEBUG_SWIFT				// test SWIFT correctness and print bucket parameters
//#define RAZERS_DEBUG					// print verification regions
#define RAZERS_PRUNE_QGRAM_INDEX		// ignore highly abundant q-grams
//#define RAZERS_MEMOPT					// optimize memory consumption
#define RAZERS_MASK_READS				// remove matches with max-hits optimal hits on-the-fly
//#define NO_PARAM_CHOOSER				// disable loss-rate parameter choosing
#define RAZERS_ISLAND_CRITERION         // island match criterion
#define RAZERS_NOOUTERREADGAPS			// enforce the alignment of the first and last base (determines the lakes)

#define RAZERS_OPENADDRESSING			// enables open addressing for the q-gram index as well as the possibility to set the load factor (-lf)
#define RAZERS_BANDED_MYERS				// uses a banded version of Myers bitvector algorithm (analogous to H. Hyyr\"o, 2001)
//#define SEQAN_OPENADDRESSING_COMPACT	// saves some memory for the openaddressing index / faster hash table access (if undefined)s
//#define RAZERS_DEBUG_MATEPAIRS
#define RAZERS_DEFER_COMPACTION         // mask duplicates on the fly and defer compaction
#define RAZERS_EXTERNAL_MATCHES         // use external memory algorithms for managing matches

#ifdef _OPENMP
#include <omp.h>
//#define _GLIBCXX_PARALLEL               // parallel STL if available
#define RAZERS_OPENADDRESSING
#endif
//#define RAZERS_PROFILE                // Extensive profiling information.
//#define RAZERS_TIMER					// output information on how fast filtration and verification as well as waiting times
//#define RAZERS_WINDOW					// use the findWindownext function on the "normal" index

#define RAZERS_MATEPAIRS				// enable paired-end matching
//#define SEQAN_USE_SSE2_WORDS			// use SSE2 128-bit integers for MyersBitVector

#include <iostream>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

#include "config.h"

#include "razers.h"
#include "outputFormat.h"
#include "paramChooser.h"

#ifdef RAZERS_WINDOW
#include "razers_window.h"
#endif

#ifdef RAZERS_MATEPAIRS
#include "razers_matepairs.h"
#include "razers_matepairs_parallel.h"
#endif

#include "razers_parallel.h"

#include "profile_timeline.h"

#ifdef RAZERS_PROFILE
#include "profile_timeline.h"
#endif  // #ifdef RAZERS_PROFILE

using namespace std;
using namespace seqan;

struct MyFragStoreConfig 
{
	typedef String<Dna5Q>	TReadSeq;
	typedef String<Dna5>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void					TReadStoreElementSpec;
	typedef Owner<ConcatDirect<> >	TReadSeqStoreSpec;
	typedef void					TMatePairStoreElementSpec;
	typedef void					TLibraryStoreElementSpec;
	typedef void					TContigStoreElementSpec;
	typedef void					TContigFileSpec;
	typedef void					TAlignedReadStoreElementSpec;
	typedef Owner<ConcatDirect<> >	TAlignedReadTagStoreSpec;
	typedef void					TAnnotationStoreElementSpec;
};

//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int mapReads(
	StringSet<CharString> & genomeFileNames,
	StringSet<CharString> & readFileNames,	// NULL terminated list of one/two read files (single/mate-pairs)
	CharString & errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	FragmentStore<MyFragStoreConfig>    store;			// stores all of the tables
	MultiFasta                          genomeSet;
	String<String<unsigned short> > 	stats;		// needed for mapping quality calculation 

	// dump configuration in verbose mode
	if (options._debugLevel >= 1) 
	{
		CharString bitmap;
		Shape<Dna, GenericShape> shape;
		stringToShape(shape, options.shape);
		shapeToString(bitmap, shape);
		
		cerr << "___SETTINGS____________" << endl;
		cerr << "Genome file:                     \t" << genomeFileNames[0] << endl;
		if (length(readFileNames) < 2)
        {
			cerr << "Read file:                       \t" << readFileNames[0] << endl;
        } 
        else
		{
			cerr << "Read files:                      \t" << readFileNames[0] << endl;
			for (unsigned i = 1; i < length(readFileNames); ++i)
				cerr << "                                 \t" << readFileNames[i] << endl;
		}
		cerr << "Compute forward matches:         \t";
		if (options.forward)	cerr << "YES" << endl;
		else				cerr << "NO" << endl;
		cerr << "Compute reverse matches:         \t";
		if (options.reverse)		cerr << "YES" << endl;
		else				cerr << "NO" << endl;
		cerr << "Allow Indels:                    \t";
		if (options.gapMode == RAZERS_GAPPED)
		    cerr << "YES" << endl;
    else
        cerr << "NO" << endl;
		cerr << "Error rate:                      \t" << options.errorRate << endl;
        if (options.threshold > 0)
            cerr << "Minimal threshold:               \t" << options.threshold << endl;
        else
            cerr << "Pigeonhole mode with overlap:    \t" << options.overlap << endl;
		cerr << "Shape:                           \t" << bitmap << endl;
		cerr << "Repeat threshold:                \t" << options.repeatLength << endl;
		cerr << "Overabundance threshold:         \t" << options.abundanceCut << endl;
        if (options.threshold > 0)
            cerr << "Taboo length:                    \t" << options.tabooLength << endl;
        if (options._debugLevel >= 1)
        {
#ifdef PLATFORM_WINDOWS
            int pid = _getpid();
#else // #ifdef PLATFORM_WINDOWS
            int pid = getpid();
#endif // #ifdef PLATFORM_WINDOWS
            cerr << "Program PID:                     \t" << pid << endl;
        }
		cerr << endl;
	}
	
	// circumvent numerical obstacles
	options.errorRate += 0.0000001;

#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_LOAD);
#endif  // #ifdef RAZERS_PROFILE
	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Load reads
	SEQAN_PROTIMESTART(load_time);

#ifdef RAZERS_MATEPAIRS
	if (length(readFileNames) == 2)
	{
		if (!loadReads(store, toCString(readFileNames[0]), toCString(readFileNames[1]), options)) {
		//if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
			cerr << "Failed to load reads" << endl;
			return RAZERS_READS_FAILED;
		}
	}
	else
#endif
	{
		if (!loadReads(store, toCString(readFileNames[0]), options)) {
			cerr << "Failed to load reads" << endl;
			return RAZERS_READS_FAILED;
		}
	} 

	if (options._debugLevel >= 1) cerr << lengthSum(store.readSeqStore) << " bps of " << length(store.readSeqStore) << " reads loaded." << endl;
	options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);

	if (options._debugLevel >= 1)
		cerr << "Loading reads took               \t" << options.timeLoadFiles << " seconds" << endl;
		
	#ifdef RAZERS_MEMOPT
		if (length(store.readSeqStore) > 16777216) {
			cerr << "more than 2^24 reads. Switch of RAZERS_MEMOPT in razers.cpp or use less." << std::endl;
			return 1;
		}
	#endif

	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Load genomes
	if (length(genomeFileNames) == 1)
	{
		int result = getGenomeFileNameList(genomeFileNames[0], genomeFileNames, options);
		if (result == RAZERS_GENOME_FAILED)
		{
			cerr << "Failed to open genome file " << genomeFileNames[0] << endl;
			return result;
		}
	}
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_LOAD);
#endif  // #ifdef RAZERS_PROFILE

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Find matches using SWIFT
	loadContigs(store, genomeFileNames, false); // add filenames to the contig store (they are loaded on-demand)
	int error = _mapReads(store, stats, options);
	if (error != 0)
	{
		switch (error)
		{
			case RAZERS_GENOME_FAILED:
				cerr << "Failed to load genomes" << endl;
				break;
			
			case RAZERS_INVALID_SHAPE:
				cerr << "Invalid Shape" << endl;
				break;
		}
		return error;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Step 4: Remove duplicates and output matches
#ifdef RAZERS_PROFILE
    timelineBeginTask(TASK_DUMP_MATCHES);
#endif  // #ifdef RAZERS_PROFILE
	if (!options.spec.DONT_DUMP_RESULTS)
		dumpMatches(store, stats, readFileNames[0], errorPrbFileName, options);
#ifdef RAZERS_PROFILE
    timelineEndTask(TASK_DUMP_MATCHES);
#endif  // #ifdef RAZERS_PROFILE

	return 0;
}	

inline void whichMacros(){
#ifdef RAZERS_OPENADDRESSING
	std::cerr << "Index:    Open addressing" << std::endl;
#else
	std::cerr << "Index:    Normal" << std::endl;
#endif
	
#ifdef RAZERS_TIMER
	std::cerr << "Timer:    ON" << std::endl;
#else
	std::cerr << "Timer:    OFF" << std::endl;
#endif

#ifdef _OPENMP
	std::cerr << "OpenMP:   ON" << std::endl;
#else
	std::cerr << "OpenMP:   OFF" << std::endl;
#endif

#ifdef RAZERS_BANDED_MYERS
	std::cerr << "Myers:    Banded" << std::endl;
#else
	std::cerr << "Myers:    Unbanded" << std::endl;
#endif

#ifdef RAZERS_PROFILE
	std::cerr << "Timeline: ON" << std::endl;
#else
	std::cerr << "Timeline: OFF" << std::endl;
#endif
	
	std::cerr << std::endl;
}


//////////////////////////////////////////////////////////////////////////////
// Command line parsing and parameter choosing
int main(int argc, const char *argv[]) 
{
	//whichMacros();

#ifdef RAZERS_PROFILE
    initTimeline();
    unsigned x = timelineAddTaskType("ON_CONTIG", "Work on contig.");
    (void)x;  // Disable warning if assertions are disable.
    SEQAN_ASSERT_EQ(x, 1u);  // The following will be OK, too.
    timelineAddTaskType("INIT", "Initialization.");
    timelineAddTaskType("REVCOMP", "Reverse-complementing contig.");
    timelineAddTaskType("FILTER", "Filtration using SWIFT.");
    timelineAddTaskType("VERIFY", "Verification of SWIFT hits.");
    timelineAddTaskType("WRITEBACK", "Write back to block-local store.");
    timelineAddTaskType("COMPACT", "Compaction");
    timelineAddTaskType("DUMP_MATCHES", "Dump matches.");
    timelineAddTaskType("LOAD", "Load input.");
    timelineAddTaskType("SORT", "Sorting.");
    timelineAddTaskType("COPY_FINDER", "Copy SWIFT Finder.");
#endif  // #ifndef RAZERS_PROFILE
	
	RazerSOptions<>			options;
	ParamChooserOptions		pm_options;

#ifdef RAZERS_MATEPAIRS
	const unsigned			maxFiles = 3;
#else
	const unsigned			maxFiles = 2;
#endif
	StringSet<CharString>	genomeFileNames;
	StringSet<CharString>	readFileNames;
	CharString				errorPrbFileName;
	

	// Change defaults
	options.forward = false;
	options.reverse = false;
	
	CommandLineParser parser;
	string rev = "$Revision$";
	addVersionLine(parser, "RazerS version 2.1 20120101 [" + rev.substr(11, 5) + "]");
	addVersionLine(parser, "Build date: " __DATE__ ", " __TIME__);

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "***********************************************************");
	addTitleLine(parser, "*** RazerS - Fast Read Mapping with Sensitivity Control ***");
	addTitleLine(parser, "***          (c) Copyright 2009 by David Weese          ***");
	addTitleLine(parser, "***********************************************************");
	addUsageLine(parser, "[OPTION]... <GENOME FILE> <READS FILE>");
#ifdef RAZERS_MATEPAIRS
	addUsageLine(parser, "[OPTION]... <GENOME FILE> <MP-READS FILE1> <MP-READS FILE2>");
#endif
	addSection(parser, "Main Options:");
	addOption(parser, CommandLineOption("f",  "forward",           "only compute forward matches", OptionType::Boolean));
	addOption(parser, CommandLineOption("r",  "reverse",           "only compute reverse complement matches", OptionType::Boolean));
	addOption(parser, CommandLineOption("i",  "percent-identity",  "set the percent identity threshold", OptionType::Double | OptionType::Label, 100 - (100.0 * options.errorRate)));
#ifndef NO_PARAM_CHOOSER
	addOption(parser, CommandLineOption("rr", "recognition-rate",  "set the percent recognition rate", OptionType::Double | OptionType::Label, 100 - (100.0 * pm_options.optionLossRate)));
	addOption(parser, CommandLineOption("mr", "mutation-rate",     "set the percent mutation rate", OptionType::Double | OptionType::Label, 100.0 * options.mutationRate));
	addOption(parser, addArgumentText(CommandLineOption("pd", "param-dir",         "folder containing user-computed parameter files (optional)", OptionType::String | OptionType::Label), "DIR"));
#endif
	addOption(parser, CommandLineOption("id", "indels",            "allow indels (default: mismatches only)", OptionType::Boolean));
#ifdef RAZERS_MATEPAIRS
	addOption(parser, CommandLineOption("ll", "library-length",    "mate-pair library length", OptionType::Int | OptionType::Label, options.libraryLength));
	addOption(parser, CommandLineOption("le", "library-error",     "mate-pair library length tolerance", OptionType::Int | OptionType::Label, options.libraryError));
#endif
	addOption(parser, CommandLineOption("m",  "max-hits",          "output only NUM of the best hits", OptionType::Int | OptionType::Label, options.maxHits));
	addOption(parser, CommandLineOption("",   "unique",            "output only unique best matches (-m 1 -dr 0 -pa)", OptionType::Boolean));
	addOption(parser, CommandLineOption("tr", "trim-reads",        "trim reads to given length (default off)", OptionType::Int | OptionType::Label));
	addOption(parser, addArgumentText(CommandLineOption("o",  "output",            "change output filename (default <READS FILE>.result)", OptionType::String), "FILE"));
	addOption(parser, CommandLineOption("v",  "verbose",           "verbose mode", OptionType::Boolean));
	addOption(parser, CommandLineOption("vv", "vverbose",          "very verbose mode", OptionType::Boolean));
	addSection(parser, "Output Format Options:");
	addOption(parser, CommandLineOption("a",  "alignment",         "dump the alignment for each match", OptionType::Boolean));
	addOption(parser, CommandLineOption("pa", "purge-ambiguous",   "purge reads with more than max-hits best matches", OptionType::Boolean));
	addOption(parser, CommandLineOption("dr", "distance-range",    "only consider matches with at most NUM more errors compared to the best (default output all)", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("of", "output-format",     "set output format", OptionType::Int | OptionType::Label, options.outputFormat));
	addHelpLine(parser, "0 = Razer format");
	addHelpLine(parser, "1 = enhanced Fasta format");
	addHelpLine(parser, "2 = Eland format");
	addHelpLine(parser, "3 = Gff format");
	addHelpLine(parser, "4 = Sam format");
	addHelpLine(parser, "5 = Amos AFG format");
	addOption(parser, CommandLineOption("gn", "genome-naming",     "select how genomes are named", OptionType::Int | OptionType::Label, options.genomeNaming));
	addHelpLine(parser, "0 = use Fasta id");
	addHelpLine(parser, "1 = enumerate beginning with 1");
	addOption(parser, CommandLineOption("rn", "read-naming",       "select how reads are named", OptionType::Int | OptionType::Label, options.readNaming));
	addHelpLine(parser, "0 = use Fasta id");
	addHelpLine(parser, "1 = enumerate beginning with 1");
	addHelpLine(parser, "2 = use the read sequence (only for short reads!)");
	addHelpLine(parser, "3 = use the Fasta id, do NOT append '/L' or '/R' for mate pairs");
	addOption(parser, CommandLineOption("",   "full-readid",       "use the whole read id (don't clip after whitespace)", OptionType::Boolean));
	addOption(parser, CommandLineOption("so", "sort-order",        "select how matches are sorted", OptionType::Int | OptionType::Label, options.sortOrder));
	addHelpLine(parser, "0 = 1. read number, 2. genome position");
	addHelpLine(parser, "1 = 1. genome position, 2. read number");
	addOption(parser, CommandLineOption("pf", "position-format",   "select begin/end position numbering", OptionType::Int | OptionType::Label, options.sortOrder));
	addHelpLine(parser, "0 = gap space");
	addHelpLine(parser, "1 = position space");
	addOption(parser, CommandLineOption("ga", "global-alignment",   "compute global alignment (in SAM output)", OptionType::Bool, options.computeGlobal));
	addSection(parser, "Misc Options:");
	addOption(parser, CommandLineOption("cm", "compact-mult", "multiply compaction treshold by this value after reaching and compacting", OptionType::Double | OptionType::Label, options.compactMult));
	addOption(parser, CommandLineOption("ncf", "no-compact-frac", "don't compact if in this last fraction of genome", OptionType::Double | OptionType::Label, options.noCompactFrac));
	addSection(parser, "Filtration Options:");
	addOption(parser, addArgumentText(CommandLineOption("s",  "shape",             "set k-mer shape", OptionType::String | OptionType::Label, options.shape), "BITSTRING"));
	addOption(parser, CommandLineOption("t",  "threshold",         "set minimum k-mer threshold (0=pigeonhole principle)", OptionType::Int | OptionType::Label, options.threshold));
	addOption(parser, CommandLineOption("ol", "overlap-length",    "set the overlap length of adjacent q-grams (pigeonhole mode)", OptionType::Int | OptionType::Label, options.overlap));
	addOption(parser, CommandLineOption("oc", "overabundance-cut", "set k-mer overabundance cut ratio", OptionType::Int | OptionType::Label, options.abundanceCut));
	addOption(parser, CommandLineOption("rl", "repeat-length",     "set simple-repeat length threshold", OptionType::Int | OptionType::Label, options.repeatLength));
	addOption(parser, CommandLineOption("tl", "taboo-length",      "set taboo length", OptionType::Int | OptionType::Label, options.tabooLength));
#ifdef RAZERS_OPENADDRESSING
	addOption(parser, CommandLineOption("lf", "load-factor", "set the load factor for the open addressing q-gram index", OptionType::Double | OptionType::Label, options.loadFactor));
#endif
	addSection(parser, "Verification Options:");
	addOption(parser, CommandLineOption("mN", "match-N",           "\'N\' matches with all other characters", OptionType::Boolean));
	addOption(parser, addArgumentText(CommandLineOption("ed", "error-distr",       "write error distribution to FILE", OptionType::String), "FILE"));
	addOption(parser, addArgumentText(CommandLineOption("mf", "mismatch-file",     "write mismatch patterns to FILE", OptionType::String), "FILE"));
	addSection(parser, "Parallelism Options:");
	addOption(parser, CommandLineOption("tc", "thread-count",   "Set the number of threads to use (0 to force sequential mode).", OptionType::Int | OptionType::Label, options.threadCount));
	addOption(parser, CommandLineOption("pws", "parallel-window-size",   "Collect SWIFT hits in windows of this length.", OptionType::Int | OptionType::Label, options.windowSize));
	addOption(parser, CommandLineOption("pvs", "parallel-verification-size",   "Verify SWIFT hits in packages of this size.", OptionType::Int | OptionType::Label, options.verificationPackageSize));
	addOption(parser, CommandLineOption("pvmpc", "parallel-verification-max-package-count",   "Largest number of packages to create for verification per thread-1, go over package size if this limit is reached..", OptionType::Int | OptionType::Label, options.maxVerificationPackageCount));
	addOption(parser, CommandLineOption("amms", "available-matches-memory-size",   "Bytes of main memory available for storing matches.  Used to switch to external sorting.  -1 for always external, 0 for never, other value as threshold.", OptionType::Int | OptionType::Label, options.availableMatchesMemorySize));
	addOption(parser, CommandLineOption("mhst", "match-histo-start-threshold",   "When to start histogram.", OptionType::Int | OptionType::Label, options.matchHistoStartThreshold));
	bool stop = !parse(parser, argc, argv, cerr);
  std::cerr << "COMMAND LINE\t";
  for (int i = 0; i < argc; ++i)
    std::cerr << " " << argv[i];
  std::cerr << "\n";
	version(parser);
  std::cerr << "DEBUG\t";
  printDebugLevel(std::cerr);
	
	//////////////////////////////////////////////////////////////////////////////
	// Extract options
	getOptionValueLong(parser, "forward", options.forward);
	getOptionValueLong(parser, "reverse", options.reverse);
	getOptionValueLong(parser, "percent-identity", options.errorRate);
	getOptionValueLong(parser, "mutation-rate", options.mutationRate);
#ifndef NO_PARAM_CHOOSER
	getOptionValueLong(parser, "recognition-rate", pm_options.optionLossRate);
	getOptionValueLong(parser, "param-dir", pm_options.paramFolder);
#endif
	bool _indels = false;
	getOptionValueLong(parser, "indels", _indels);
	options.gapMode = (_indels)? RAZERS_GAPPED: RAZERS_UNGAPPED;
#ifdef RAZERS_MATEPAIRS
	getOptionValueLong(parser, "library-length", options.libraryLength);
	getOptionValueLong(parser, "library-error", options.libraryError);
#endif
	getOptionValueLong(parser, "max-hits", options.maxHits);
	getOptionValueLong(parser, "purge-ambiguous", options.purgeAmbiguous);
	getOptionValueLong(parser, "distance-range", options.scoreDistanceRange);
	if (isSetLong(parser, "distance-range")) options.scoreDistanceRange++;
	getOptionValueLong(parser, "alignment", options.dumpAlignment);
	getOptionValueLong(parser, "output", options.output);
	getOptionValueLong(parser, "output-format", options.outputFormat);
	getOptionValueLong(parser, "sort-order", options.sortOrder);
	getOptionValueLong(parser, "genome-naming", options.genomeNaming);
	getOptionValueLong(parser, "read-naming", options.readNaming);
	getOptionValueLong(parser, "full-readid", options.fullFastaId);    
	getOptionValueLong(parser, "position-format", options.positionFormat);
	getOptionValueLong(parser, "compact-mult", options.compactMult);
	getOptionValueLong(parser, "no-compact-frac", options.noCompactFrac);
	if (isSetLong(parser, "global-alignment")) options.computeGlobal = true;
	getOptionValueLong(parser, "shape", options.shape);
	getOptionValueLong(parser, "threshold", options.threshold);
	getOptionValueLong(parser, "overlap-length", options.overlap);
	getOptionValueLong(parser, "overabundance-cut", options.abundanceCut);
	getOptionValueLong(parser, "repeat-length", options.repeatLength);
    getOptionValueLong(parser, "thread-count", options.threadCount);
    getOptionValueLong(parser, "parallel-window-size", options.windowSize);
    getOptionValueLong(parser, "parallel-verification-size", options.verificationPackageSize);
    getOptionValueLong(parser, "parallel-verification-max-package-count", options.maxVerificationPackageCount);
	getOptionValueLong(parser, "available-matches-memory-size", options.availableMatchesMemorySize);
    getOptionValueLong(parser, "match-histo-start-threshold", options.matchHistoStartThreshold);
#ifdef RAZERS_OPENADDRESSING
	getOptionValueLong(parser, "load-factor", options.loadFactor);
#endif 
	getOptionValueLong(parser, "trim-reads", options.trimLength);
	getOptionValueLong(parser, "taboo-length", options.tabooLength);
	getOptionValueLong(parser, "match-N", options.matchN);
	getOptionValueLong(parser, "error-distr", errorPrbFileName);
	getOptionValueLong(parser, "mismatch-file", options.mismatchFilename);
	if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit
	if (isSetLong(parser, "verbose")) options._debugLevel = max(options._debugLevel, 1);
	if (isSetLong(parser, "vverbose")) options._debugLevel = max(options._debugLevel, 3);
	if (isSetLong(parser, "unique"))
	{
		options.maxHits = 1;
		options.scoreDistanceRange = 1;
		options.purgeAmbiguous = true;
	}
	if (!options.forward && !options.reverse)  // enable both per default
	{
		options.forward = true;
		options.reverse = true;
	}
    // don't append /L/R in SAM mode
    if (!isSetLong(parser, "read-naming") && options.outputFormat == 4)
        options.readNaming = 3;
	appendValue(genomeFileNames, getArgumentValue(parser, 0));
	for (unsigned i = 1; i < argumentCount(parser) && i < maxFiles; ++i)
		appendValue(readFileNames, getArgumentValue(parser, i), Generous());

	//////////////////////////////////////////////////////////////////////////////
	// Check options
	if ((options.errorRate < 50 || options.errorRate > 100) && (stop = true))
		cerr << "Percent identity threshold must be a value between 50 and 100" << endl;
	if ((pm_options.optionLossRate < 80 || pm_options.optionLossRate > 100) && (stop = true))
		cerr << "Recognition rate must be a value between 80 and 100" << endl;
	if ((options.mutationRate < 0 || options.mutationRate > 20) && (stop = true))
		cerr << "Mutation rate must be a value between 0 and 20" << endl;
#ifdef RAZERS_MATEPAIRS
	if ((options.libraryLength <= 0) && (stop = true))
		cerr << "Library length must be a value greater 0" << endl;
	if ((options.libraryError <= 0) && (stop = true))
		cerr << "Library error must be a value greater or equal 0" << endl;
#endif
	if ((options.maxHits < 1) && (stop = true))
		cerr << "Maximum hits threshold must be greater than 0" << endl;
	if ((options.outputFormat > 5) && (stop = true))
		cerr << "Invalid output format option." << endl;
	if ((options.sortOrder > 1) && (stop = true))
		cerr << "Invalid sort order options." << endl;
	if ((options.genomeNaming > 1) && (stop = true))
		cerr << "Invalid genome naming options." << endl;
	if ((options.readNaming > 3) && (stop = true))
		cerr << "Invalid read naming options." << endl;
	if ((options.positionFormat > 1) && (stop = true))
		cerr << "Invalid position format options." << endl;
	if ((options.threshold < 0) && (stop = true))
		cerr << "Threshold must be a value greater or equal 0" << endl;
	if ((options.overlap < 0) && (stop = true))
		cerr << "Overlap length must be a value greater or equal 0" << endl;
	if (isSetLong(parser, "threshold") && (options.threshold > 0) && isSetLong(parser, "overlap-length") && (stop = true))
		cerr << "Overlap length can only be set in pigeonhole mode (threshold == 0)" << endl;
	if (isSetLong(parser, "shape"))
	{
		unsigned ones = 0;
		unsigned zeros = 0;
		for(unsigned i = 0; i < length(options.shape); ++i)
			switch (options.shape[i])
			{
				case '0':
					++zeros;
					break;
				case '1':
					++ones;
					break;
				default:
					cerr << "Shape must be a binary string" << endl;
					stop = true;
					i = length(options.shape);
			}
		if ((ones == 0 || ones > 31) && !stop) 
		{
			cerr << "Invalid Shape" << endl;
			stop = true;
		}

        unsigned maxOnes = 14;
#ifdef RAZERS_OPENADDRESSING
        maxOnes = 31;
#endif
        if ((ones < 7 || ones > maxOnes) && !stop)
			cerr << "Warning: Shape should contain at least 7 and at most " << maxOnes << " '1's" << endl;
        options.delta = ones + zeros;
	}
	if ((options.abundanceCut <= 0 || options.abundanceCut > 1) && (stop = true))
		cerr << "Overabundance cut ratio must be a value >0 and <=1. Set to 1 to disable." << endl;
	if ((options.repeatLength <= 1) && (stop = true))
		cerr << "Repeat length must be a value greater than 1" << endl;
	if ((options.trimLength != 0 && options.trimLength < 14) && (stop = true))
		cerr << "Minimum read length is 14" << endl;
	if ((options.tabooLength < 1) && (stop = true))
		cerr << "Taboo length must be a value greater than 0" << endl;
	if (argumentCount(parser) == 2)
		options.libraryLength = -1;		// only 1 readset -> disable mate-pair mapping
	if ((argumentCount(parser) > maxFiles) && (stop = true))
		cerr << "More than " << maxFiles << " input files specified." << endl;
	if ((argumentCount(parser) < 2) && (stop = true))
	{
		if (argc > 1)
			cerr << "Less than 2 input files specified." << endl;
		else
		{
			shortHelp(parser, cerr);	// print short help and exit
			return 0;
		}
	}

	options.errorRate = (100.0 - options.errorRate) / 100.0;
	options.mutationRate = options.mutationRate / 100.0;
	options.lossRate = pm_options.optionLossRate = (100.0 - pm_options.optionLossRate) / 100.0;
	if (stop)
	{
		cerr << "Exiting ..." << endl;
		return RAZERS_INVALID_OPTIONS;
	}

#ifdef _OPENMP
    // Set maximal number of threads.
    int oldMaxThreads = omp_get_max_threads();
    omp_set_num_threads(options.threadCount == 0 ? 1 : options.threadCount);
#endif  // #ifdef _OPENMP

	//////////////////////////////////////////////////////////////////////////////
	// get read length
	int readLength = estimateReadLength(toCString(readFileNames[0]));
	if (readLength == RAZERS_READS_FAILED)
	{
		cerr << "Failed to open reads file " << readFileNames[0] << endl;
		cerr << "Exiting ..." << endl;
		return RAZERS_READS_FAILED;
	}
	if (readLength == 0) {
		cerr << "Failed to read the first read sequence." << endl;
		cerr << "Exiting ..." << endl;
		return RAZERS_READS_FAILED;
	}

	if (options.trimLength > readLength)
		options.trimLength = readLength;
	
#ifndef NO_PARAM_CHOOSER
	if (!(isSetLong(parser, "shape") || isSetLong(parser, "threshold")))
	{
		if (options.lowMemory) pm_options.maxWeight = 13;
		pm_options.verbose = (options._debugLevel >= 1);
		pm_options.optionErrorRate = options.errorRate;
		if (options.gapMode == RAZERS_UNGAPPED)
		{
			pm_options.optionProbINSERT = (ParamChooserOptions::TFloat)0.0;
			pm_options.optionProbDELETE = (ParamChooserOptions::TFloat)0.0;
		}
		else
		{
			pm_options.optionProbINSERT = (ParamChooserOptions::TFloat)0.01;	//this number is basically without meaning, any value > 0 will lead to
			pm_options.optionProbDELETE = (ParamChooserOptions::TFloat)0.01;	//edit distance parameter choosing
		}

		if (empty(pm_options.paramFolder)) 
		{
			string razersFolder = argv[0];
			size_t lastPos = razersFolder.find_last_of('/') + 1;
			if (lastPos == razersFolder.npos + 1) lastPos = razersFolder.find_last_of('\\') + 1;
			if (lastPos == razersFolder.npos + 1) lastPos = 0;
			razersFolder.erase(lastPos); 
			pm_options.paramFolderPath = razersFolder;
		}
		if (options.trimLength > 0) readLength = options.trimLength;
		if (readLength > 0)
		{
/*			if(options.maqMapping && readLength != options.artSeedLength)
				pm_options.totalN = options.artSeedLength;
			else*/
				pm_options.totalN = readLength;
			if (options._debugLevel >= 1)
				cerr << "___PARAMETER_CHOOSING__" << endl;
			if (!chooseParams(options,pm_options))
			{
				if (pm_options.verbose) 
					cerr << "Couldn't find preprocessed parameter files. Please configure manually (options --shape and --threshold)." << endl;
                if (options._debugLevel >= 1)
                    cerr << "Using default configurations (shape = " << options.shape << " and q-gram lemma)." << endl;
			}
			if (options._debugLevel >= 1) cerr << endl;
		} else
		{
			cerr << "Failed to load reads" << endl;
			cerr << "Exiting ..." << endl;
			return RAZERS_READS_FAILED;
		}
	}
#endif	

	int result = mapReads(genomeFileNames, readFileNames, errorPrbFileName, options);
	if (result != 0)
		cerr << "Exiting ..." << endl;

#ifdef _OPENMP
    // Restoring number of threads for side-effect freeness.
    omp_set_num_threads(oldMaxThreads);
#endif  // #ifdef _OPENMP

#ifdef RAZERS_PROFILE
    dumpTimeline("razers.profile.txt", true);
#endif  // #ifndef RAZERS_PROFILE

	return result;
}
