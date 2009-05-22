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

#define SEQAN_PROFILE			// enable time measuring
//#define SEQAN_DEBUG_SWIFT		// test SWIFT correctness and print bucket parameters
//#define RAZERS_DEBUG			// print verification regions
#define RAZERS_PRUNE_QGRAM_INDEX
#define RAZERS_CONCATREADS		// use <ConcatDirect> StringSet to store reads
#define RAZERS_MEMOPT			// optimize memory consumption
#define RAZERS_MASK_READS		// remove matches with max-hits optimal hits on-the-fly
//#define NO_PARAM_CHOOSER
//#define RAZERS_PARALLEL			// parallelize using Intel's Threading Building Blocks
#define RAZERS_MATEPAIRS
//#define RAZERS_DIRECT_MAQ_MAPPING
//#define SEQAN_USE_SSE2_WORDS	// use SSE2 128-bit integers for MyersBitVector

#include "seqan/platform.h"
#ifdef PLATFORM_WINDOWS
	#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
	#define SEQAN_DEFAULT_TMPDIR "./"
#endif

#include "razers.h"
#include "outputFormat.h"
#include "paramChooser.h"

#ifdef RAZERS_PARALLEL
#include "razers_parallel.h"
#endif

#ifdef RAZERS_MATEPAIRS
#include "razers_matepairs.h"
#endif


#include <iostream>
#include <sstream>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int mapReads(
	const char *genomeFileName,
	const char *readFileNames[],	// NULL terminated list of one/two read files (single/mate-pairs)
	const char *errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	MultiFasta				genomeSet;
	TReadSet				readSet;
	StringSet<CharString>	genomeNames;	// genome names, taken from the Fasta file
	StringSet<CharString>	readNames;		// read names, taken from the Fasta file
	TMatches				matches;		// resulting forward/reverse matches
	String<String<unsigned short> > 	stats;		// needed for mapping quality calculation 

	// dump configuration in verbose mode
	if (options._debugLevel >= 1) 
	{
//		CharString bitmap;
//		shapeToString(bitmap, shape);
		::std::cerr << "___SETTINGS____________" << ::std::endl;
		::std::cerr << "Genome file:                     \t" << genomeFileName << ::std::endl;
		if (empty(readFileNames[1]))
			::std::cerr << "Read file:                       \t" << readFileNames[0] << ::std::endl;
		else
		{
			::std::cerr << "Read files:                      \t" << readFileNames[0] << ::std::endl;
			for (const char **i = readFileNames + 1; !empty(*i); ++i)
				::std::cerr << "                                 \t" << *i << ::std::endl;
		}
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
	// Step 1: Load fasta files and determine genome file type
	SEQAN_PROTIMESTART(load_time);

#ifdef RAZERS_MATEPAIRS
	if (!empty(readFileNames[1]))
	{
		if (!loadReads(readSet, readNames, readFileNames[0], readFileNames[1], options)) {
		//if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
			::std::cerr << "Failed to load reads" << ::std::endl;
			return RAZERS_READS_FAILED;
		}
	}
	else
#endif
	{
		if (!loadReads(readSet, readNames, readFileNames[0], options)) {
		//if (!loadReads(readSet, readQualities, readNames, readFileNames[0], readFileNames[1], options)) {
			::std::cerr << "Failed to load reads" << ::std::endl;
			return RAZERS_READS_FAILED;
		}
	} 

	if (options._debugLevel >= 1) ::std::cerr << lengthSum(readSet) << " bps of " << length(readSet) << " reads loaded." << ::std::endl;
	options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);

	if (options._debugLevel >= 1)
		::std::cerr << "Loading reads took               \t" << options.timeLoadFiles << " seconds" << ::std::endl;

	StringSet<CharString> genomeFileNameList;
	int result = getGenomeFileNameList(genomeFileName, genomeFileNameList, options);
	if (result == RAZERS_GENOME_FAILED)
	{
		::std::cerr << "Failed to open genome file " << genomeFileName << ::std::endl;
		return result;
	}


	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Find matches using SWIFT
#ifdef RAZERS_PARALLEL
    typedef typename RazerSOptions<TSpec>::TMutex TMutex;
    options.patternMutex = new TMutex[length(readSet)];
#endif

	::std::map<unsigned,::std::pair< ::std::string,unsigned> > gnoToFileMap; //map to retrieve genome filename and sequence number within that file
	int error = mapReads(matches, genomeFileNameList, genomeNames, gnoToFileMap, readSet, stats, options);
	if (error != 0)
	{
		switch (error)
		{
			case RAZERS_GENOME_FAILED:
				::std::cerr << "Failed to load genomes" << ::std::endl;
				break;
			
			case RAZERS_INVALID_SHAPE:
				::std::cerr << "Invalid Shape" << endl << endl;
				break;
		}
		return error;
	}

#ifdef RAZERS_PARALLEL
    delete[] options.patternMutex;
#endif

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Remove duplicates and output matches
	if (!options.spec.DONT_DUMP_RESULTS)
		dumpMatches(matches, genomeNames, genomeFileNameList, gnoToFileMap, readSet, stats, readNames, readFileNames[0], errorPrbFileName, options);



	return 0;
}	


//////////////////////////////////////////////////////////////////////////////
// Print usage
void printVersion() 
{
	string rev = "$Revision$";
	cerr << "RazerS version 0.3 20090223 [" << rev.substr(11, 4) << "] (prerelease)" << endl;
}

template <typename TSpec>
void printHelp(int, const char *[], RazerSOptions<TSpec> &options, ParamChooserOptions &pm_options, bool longHelp = false) 
{
	if (options.printVersion)
	{
		printVersion();
		return;
	}

	cerr << "********************************************" << endl;
	cerr << "*** RazerS - Fast Mapping of Short Reads ***" << endl;
	cerr << "*** written by David Weese (c) Aug 2008  ***" << endl;
	cerr << "********************************************" << endl << endl;
	cerr << "Usage: razers [OPTION]... <GENOME FILE> <READS FILE>" << endl;
#ifdef RAZERS_MATEPAIRS
	cerr << "       razers [OPTION]... <GENOME FILE> <MP-READS FILE1> <MP-READS FILE2>" << endl;
#endif
	if (longHelp) {
		cerr << endl << "Main Options:" << endl;
		cerr << "  -f,  --forward               \t" << "only compute forward matches" << endl;
		cerr << "  -r,  --reverse               \t" << "only compute reverse complement matches" << endl;
		cerr << "  -i,  --percent-identity NUM  \t" << "set the percent identity threshold (default " << 100 - (100.0 * options.errorRate) << ')' << endl;
#ifndef NO_PARAM_CHOOSER
		cerr << "  -rr, --recognition-rate NUM  \t" << "set the percent recognition rate (default " << 100 - (100.0 * pm_options.optionLossRate) << ')' << endl;
#endif
		cerr << "  -id, --indels                \t" << "allow indels (default: mismatches only)" << endl;
#ifdef RAZERS_MATEPAIRS
		cerr << "  -ll, --library-length NUM    \t" << "mate-pair library length (default " << options.libraryLength << ')' << endl;
		cerr << "  -le, --library-error NUM     \t" << "mate-pair library length tolerance (default " << options.libraryError << ')' << endl;
#endif
		cerr << "  -m,  --max-hits NUM          \t" << "output only NUM of the best hits (default " << options.maxHits << ')' << endl;
		cerr << "  -tr, --trim-reads NUM        \t" << "trim reads to length NUM (default off)" << endl;
		cerr << "  -o,  --output FILE           \t" << "change output filename (default <READS FILE>.result)" << endl;
		cerr << "  -v,  --verbose               \t" << "verbose mode" << endl;
		cerr << "  -vv, --vverbose              \t" << "very verbose mode" << endl;
		cerr << "  -V,  --version               \t" << "print version number" << endl;
		cerr << "  -h,  --help                  \t" << "print this help" << endl;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
		cerr << "  -lm, --low-memory            \t" << "decrease memory usage at the expense of runtime" << endl;
		cerr << endl << "Mapping Quality Options:" << endl;
		cerr << "  -mq, --mapping-quality       \t" << "switch on mapping quality mode" << endl;
		cerr << "  -nbi,--no-below-id           \t" << "do not report matches with seed identity < percent id (default off)" << endl;
		cerr << "  -qsl,--mq-seed-length NUM    \t" << "seed length used for mapping quality assignment ("<<options.artSeedLength<<')' << endl;
#endif
		cerr << endl << "Output Format Options:" << endl;
		cerr << "  -a,  --alignment             \t" << "dump the alignment for each match" << endl;
#ifdef RAZERS_MASK_READS
		cerr << "  -pa, --purge-ambiguous       \t" << "purge reads with more than max-hits best matches" << endl;
#endif
		cerr << "  -dr, --distance-range NUM    \t" << "only consider matches with at most NUM more errors compared to the best (default output all)" << endl;
		cerr << "  -of, --output-format NUM     \t" << "set output format" << endl;
		cerr << "                               \t" << "0 = Razer format (default, see README)" << endl;
		cerr << "                               \t" << "1 = enhanced Fasta format" << endl;
		cerr << "                               \t" << "2 = Eland format" << endl;
		cerr << "                               \t" << "3 = GFF format" << endl;
		cerr << "  -gn, --genome-naming NUM     \t" << "select how genomes are named" << endl;
		cerr << "                               \t" << "0 = use Fasta id (default)" << endl;
		cerr << "                               \t" << "1 = enumerate beginning with 1" << endl;
		cerr << "  -rn, --read-naming NUM       \t" << "select how reads are named" << endl;
		cerr << "                               \t" << "0 = use Fasta id (default)" << endl;
		cerr << "                               \t" << "1 = enumerate beginning with 1" << endl;
		cerr << "                               \t" << "2 = use the read sequence (only for short reads!)" << endl;
		cerr << "  -so, --sort-order NUM        \t" << "select how matches are sorted" << endl;
		cerr << "                               \t" << "0 = 1. read number, 2. genome position (default)" << endl;
		cerr << "                               \t" << "1 = 1. genome position, 2. read number" << endl;
		cerr << "  -pf, --position-format       \t" << "0 = gap space (default)" << endl;
		cerr << "                               \t" << "1 = position space" << endl;
		cerr << endl << "Filtration Options:" << endl;
		cerr << "  -s,  --shape BITSTRING       \t" << "set k-mer shape (binary string, default " << options.shape << ')' << endl;
		cerr << "  -t,  --threshold NUM         \t" << "set minimum k-mer threshold (default " << options.threshold << ')' << endl;
		cerr << "  -oc, --overabundance-cut NUM \t" << "set k-mer overabundance cut ratio (default " << options.abundanceCut << ')' << endl;
		cerr << "  -rl, --repeat-length NUM     \t" << "set simple-repeat length threshold (default " << options.repeatLength << ')' << endl;
		cerr << "  -tl, --taboo-length NUM      \t" << "set taboo length (default " << options.tabooLength << ')' << endl;
		cerr << "  -pd, --param-dir STR         \t" << "folder containing user-computed parameter files (optional)" << endl;
		cerr << endl << "Verification Options:" << endl;
		cerr << "  -mN, --match-N               \t" << "\'N\' matches with all other characters" << endl;
		cerr << "  -ed, --error-distr FILE      \t" << "write error distribution to FILE" << endl;
	} else {
		cerr << "Try 'razers --help' for more information." << endl;
	}
}


template<typename TSpec>
int getGenomeFileNameList(char const * filename, StringSet<CharString> & genomeFileNames, RazerSOptions<TSpec> &options)
{
	::std::ifstream file;
	file.open(filename,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return RAZERS_GENOME_FAILED;

	char c = _streamGet(file);
	if (c != '>' && c != '@')	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		if(options._debugLevel >=1)
			::std::cout << ::std::endl << "Reading multiple genome files:" <<::std::endl;
		
/*		//locations of genome files are relative to list file's location
		::std::string tempGenomeFile(filename);
		size_t lastPos = tempGenomeFile.find_last_of('/') + 1;
		if (lastPos == tempGenomeFile.npos) lastPos = tempGenomeFile.find_last_of('\\') + 1;
		if (lastPos == tempGenomeFile.npos) lastPos = 0;
		::std::string filePrefix = tempGenomeFile.substr(0,lastPos);*/
		unsigned i = 1;
		while(!_streamEOF(file))
		{ 
			_parse_skipWhitespace(file, c);
			appendValue(genomeFileNames,_parse_readFilepath(file,c));
			//CharString currentGenomeFile(filePrefix);
			//append(currentGenomeFile,_parse_readFilepath(file,c));
			//appendValue(genomeFileNames,currentGenomeFile);
			if(options._debugLevel >=2)
				::std::cout <<"Genome file #"<< i <<": " << genomeFileNames[length(genomeFileNames)-1] << ::std::endl;
			++i;
			_parse_skipWhitespace(file, c);
		}
		if(options._debugLevel >=1)
			::std::cout << i-1 << " genome files total." <<::std::endl;
	}
	else		//if file starts with a fasta header --> regular one-genome-file input
		appendValue(genomeFileNames,filename);
	file.close();
	return 0;

}


int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line

	RazerSOptions<>		options;
	ParamChooserOptions	pm_options;

	bool				paramChoosing = true;
	unsigned			fnameCount = 0;
	string				errorPrbFileName;
#ifdef RAZERS_MATEPAIRS
	const unsigned		maxFiles = 3;
#else
	const unsigned		maxFiles = 2;
#endif
	const char			*fname[maxFiles + 1] = { NULL, };

	options.forward = false;
	options.reverse = false;
	options.hammingOnly = true;


	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse options

			if (strcmp(argv[arg], "-f") == 0 || strcmp(argv[arg], "--forward") == 0) {
				options.forward = true;
				continue;
			}
			if (strcmp(argv[arg], "-r") == 0 || strcmp(argv[arg], "--reverse") == 0) {
				options.reverse = true;
				continue;
			}
			if (strcmp(argv[arg], "-i") == 0 || strcmp(argv[arg], "--percent-identity") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.errorRate;
					options.errorRate = (100.0 - options.errorRate) / 100.0;
					if (!istr.fail()) 
					{
						if (options.errorRate < 0 || options.errorRate > 0.5)
							cerr << "Percent identity threshold must be a value between 50 and 100" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-rr") == 0 || strcmp(argv[arg], "--recognition-rate") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionLossRate;
					if (!istr.fail())
					{
						if (pm_options.optionLossRate < 80 || pm_options.optionLossRate > 100)
							cerr << "Recognition rate must be a value between 80 and 100" << endl << endl;
						else
						{
							pm_options.optionLossRate = 100.0-pm_options.optionLossRate;
							pm_options.optionLossRate /= 100.0;
							continue;
						}
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}

#ifdef RAZERS_MATEPAIRS
			if (strcmp(argv[arg], "-ll") == 0 || strcmp(argv[arg], "--library-length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.libraryLength;
					if (!istr.fail())
					{
						if (options.libraryLength <= 0)
							cerr << "Library length must be a value greater 0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-le") == 0 || strcmp(argv[arg], "--library-error") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.libraryError;
					if (!istr.fail())
					{
						if (options.libraryError <= 0)
							cerr << "Library error must be a value greater or equal 0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
#endif

			if (strcmp(argv[arg], "-pa") == 0 || strcmp(argv[arg], "--purge-ambiguous") == 0) {
				options.purgeAmbiguous = true;
				continue;
			}
			if (strcmp(argv[arg], "-m") == 0 || strcmp(argv[arg], "--max-hits") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.maxHits;
					if (!istr.fail()) 
					{
						if (options.maxHits < 1)
							cerr << "Maximum hits threshold must be greater than 0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-dr") == 0 || strcmp(argv[arg], "--distance-range") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.distanceRange;
					if (!istr.fail()) 
					{
						options.distanceRange++;
						continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-pd") == 0 || strcmp(argv[arg], "--param-dir") == 0) { //prefix for previously computed param files
				if (arg + 1 < argc) {
					++arg;
					pm_options.paramFolder = argv[arg];
					continue;
//					cout << "Session id prefix specified\n";
				}
				else 
				{
					printHelp(argc, argv, options, pm_options);
					return 0;
				}
			}

			if (strcmp(argv[arg], "-id") == 0 || strcmp(argv[arg], "--indels") == 0) {
				options.hammingOnly = false;
				continue;
			}
			if (strcmp(argv[arg], "-a") == 0 || strcmp(argv[arg], "--alignment") == 0) {
				options.dumpAlignment = true;
				continue;
			}
			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, pm_options);
					return 0;
				}
				++arg;
				options.output = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-of") == 0 || strcmp(argv[arg], "--output-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.outputFormat;
					if (!istr.fail())
					{
						if (options.outputFormat > 3)
							cerr << "Invalid output format options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-so") == 0 || strcmp(argv[arg], "--sort-order") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.sortOrder;
					if (!istr.fail())
					{
						if (options.sortOrder > 1)
							cerr << "Invalid sort order options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-gn") == 0 || strcmp(argv[arg], "--genome-naming") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.genomeNaming;
					if (!istr.fail())
					{
						if (options.genomeNaming > 1)
							cerr << "Invalid genome naming options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-rn") == 0 || strcmp(argv[arg], "--read-naming") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.readNaming;
					if (!istr.fail())
					{
						if (options.readNaming > 2)
							cerr << "Invalid read naming options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-pf") == 0 || strcmp(argv[arg], "--position-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.positionFormat;
					if (!istr.fail())
					{
						if (options.positionFormat > 1)
							cerr << "Invalid position format options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-s") == 0 || strcmp(argv[arg], "--shape") == 0){
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					paramChoosing = false;
					istr >> options.shape;
					if (istr.fail())
						cerr << "Could not parse shape" << endl << endl;
					else
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
									cerr <<	"Shape must be a binary string" << endl << endl;
									printHelp(argc, argv, options, pm_options);
									return 0;
							}
						if (ones == 0 || ones > 20) 
						{
							cerr <<	"Invalid Shape" << endl << endl;
							printHelp(argc, argv, options, pm_options);
							return 0;
						}

						if (ones < 7 || ones > 14)
							cerr <<	"Warning: Shape should contain at least 7 and at most 14 '1's" << endl << endl;

						continue;				
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-oc") == 0 || strcmp(argv[arg], "--overabundance-cut") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.abundanceCut;
					if (!istr.fail())
					{
						if (options.abundanceCut <= 0 || options.abundanceCut > 1)
							cerr << "Overabundance cut ratio must be a value >0 and <=1. Set to 1 to disable." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-rl") == 0 || strcmp(argv[arg], "--repeat-length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.repeatLength;
					if (!istr.fail())
					{
						if (options.repeatLength <= 1)
							cerr << "Repeat length must be a value greater than 1" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-t") == 0 || strcmp(argv[arg], "--threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					paramChoosing = false;
					istringstream istr(argv[arg]);
					istr >> options.threshold;
					if (!istr.fail())
					{
						if (options.threshold < 1)
							cerr << "Threshold must be a value greater than 0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
#ifdef RAZERS_DIRECT_MAQ_MAPPING
			if (strcmp(argv[arg], "--quality-in-header") == 0) {
				options.fastaIdQual = true;
				continue;
			}
			if (strcmp(argv[arg], "-mq") == 0 || strcmp(argv[arg], "--mapping-quality") == 0) {
				options.maqMapping = true;
				continue;
			}
			if (strcmp(argv[arg], "-nbi") == 0 || strcmp(argv[arg], "--no-below-id") == 0) {
				options.noBelowIdentity = true;
				continue;
			}
			if (strcmp(argv[arg], "-qsl") == 0 || strcmp(argv[arg], "--mq-seed-length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.artSeedLength;
					if (!istr.fail()) 
					{
						if (options.artSeedLength < 24)
							cerr << "Minimum seed length is 24" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-smq") == 0 || strcmp(argv[arg], "--seed-mism-quality") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.maxMismatchQualSum;
					if (!istr.fail()) 
					{
						if (options.maxMismatchQualSum < 0  )
							cerr << "Max seed mismatch quality needs to be 0 or a positive value" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-tmq") == 0 || strcmp(argv[arg], "--total-mism-quality") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.absMaxQualSumErrors;
					if (!istr.fail()) 
					{
						if (options.absMaxQualSumErrors < 0  )
							cerr << "Max total mismatch quality needs to be 0 or a positive value" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
#endif
			if (strcmp(argv[arg], "-lm") == 0 || strcmp(argv[arg], "--low-memory") == 0) {
				options.lowMemory = true;
				continue;
			}
	
			if (strcmp(argv[arg], "-tr") == 0 || strcmp(argv[arg], "--trim-reads") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.trimLength;
					if (!istr.fail()) 
					{
						if (options.trimLength < 14)
							cerr << "Minimum read length is 14" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}

			if (strcmp(argv[arg], "-tl") == 0 || strcmp(argv[arg], "--taboo-length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.tabooLength;
					if (!istr.fail())
					{
						if (options.tabooLength < 1)
							cerr << "Taboo length must be a value greater than 0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, options, pm_options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 1);
				continue;
			}
			if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--vverbose") == 0) {
				options._debugLevel = 3;
				continue;
			}
			if (strcmp(argv[arg], "-V") == 0 || strcmp(argv[arg], "--version") == 0) {
				options.printVersion = true;
				continue;
			}
			if (strcmp(argv[arg], "-mN") == 0 || strcmp(argv[arg], "--match-N") == 0) {
				options.matchN = true;
				continue;
			}
			if (strcmp(argv[arg], "-ed") == 0 || strcmp(argv[arg], "--error-distr") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, pm_options);
					return 0;
				}
				++arg;
				errorPrbFileName = argv[arg];
				continue;
			}
			cerr << "Unknown option: " << argv[arg] << endl << endl;
			printHelp(argc, argv, options, pm_options);
			return 0;
		} else {
			// parse file name
			if (fnameCount == maxFiles) {
				cerr << "More than " << maxFiles << " input files specified." << endl << endl;
				printHelp(argc, argv, options, pm_options);
				return 0;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 2) {
		if (argc > 1 && !options.printVersion)
			cerr << "Less than 2 input files specified." << endl << endl;
		printHelp(argc, argv, options, pm_options);
		return 0;
	}
	if (!options.forward && !options.reverse) { // enable both per default
		options.forward = true;
		options.reverse = true;
	}

	if (fnameCount == 2)
		options.libraryLength = -1;	// only 1 readset -> disable mate-pair mapping
	
	if (options.printVersion)
		printVersion();
		
	// get read length
	int readLength = estimateReadLength(fname[1]);
	if (readLength == RAZERS_READS_FAILED)
	{
		::std::cerr << "Failed to open reads file " << fname[1] << ::std::endl;
		return 0;
	}
	if (readLength == 0) {
		::std::cerr << "Failed to read the first read sequence.";
		return 0;
	}

	if (options.trimLength > readLength)
		options.trimLength = readLength;
	
#ifndef NO_PARAM_CHOOSER
	if (paramChoosing)
	{
		if(options.lowMemory) pm_options.maxWeight = 13;
		pm_options.verbose = (options._debugLevel >= 1);
		pm_options.optionErrorRate = options.errorRate;
		if (options.hammingOnly)
		{
			pm_options.optionProbINSERT = 0.0;
			pm_options.optionProbDELETE = 0.0;
		}
		else
		{
			pm_options.optionProbINSERT = 0.01;	//this number is basically without meaning, any value > 0 will lead to
			pm_options.optionProbDELETE = 0.01;	//edit distance parameter choosing
		}

		if(length(pm_options.paramFolder) == 0) 
		{
			pm_options.paramFolderPath = argv[0];
			size_t lastPos = pm_options.paramFolderPath.find_last_of('/') + 1;
			if (lastPos == pm_options.paramFolderPath.npos + 1) lastPos = pm_options.paramFolderPath.find_last_of('\\') + 1;
			if (lastPos == pm_options.paramFolderPath.npos + 1) lastPos = 0;
			pm_options.paramFolderPath.erase(lastPos); 
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
				cerr << "Using default configurations (shape = " << options.shape << " and q-gram lemma)." << endl;
			}
			cerr << endl;
		} else
		{
			cerr << "Failed to load reads" << endl;
			return 0;
		}
	}
#endif	

#ifdef RAZERS_PARALLEL
	tbb::task_scheduler_init scheduler;
#endif

	int result = mapReads(fname[0], fname + 1, errorPrbFileName.c_str(), options);
	if (result == RAZERS_INVALID_SHAPE) 
	{
		printHelp(argc, argv, options, pm_options);
		return 0;
	}
	return result;
}
