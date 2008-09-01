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

#define SEQAN_PROFILE			// enable time measuring
//#define SEQAN_DEBUG_SWIFT		// test SWIFT correctness and print bucket parameters
//#define RAZERS_PROFILE		// omit dumping results
//#define RAZERS_DEBUG			// print verification regions
#define RAZERS_PRUNE_QGRAM_INDEX
//#define RAZERS_HAMMINGVERIFY	// allow only mismatches, no indels
#define RAZERS_MAXHITS			// drop reads with too many matches

#include "seqan/platform.h"
#ifdef PLATFORM_WINDOWS
	#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
	#define SEQAN_DEFAULT_TMPDIR "./"
#endif

#include "razers.h"

#include <iostream>
#include <sstream>

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
// Print usage
template <typename TSpec>
void printHelp(int, const char *[], RazerSOptions<TSpec> &options, bool longHelp = false) 
{
	cerr << "********************************************" << endl;
	cerr << "*** RazerS - Fast Mapping of Short Reads ***" << endl;
	cerr << "*** written by David Weese (c) Aug 2008  ***" << endl;
	cerr << "********************************************" << endl << endl;
	cerr << "Usage: razers [OPTION]... <GENOME FILE> <READS FILE>" << endl;
	if (longHelp) {
		cerr << endl << "Main Options:" << endl;
		cerr << "  -f,  --forward               \t" << "only compute forward matches" << endl;
		cerr << "  -r,  --reverse               \t" << "only compute reverse complement matches" << endl;
		cerr << "  -i,  --percent-identity NUM  \t" << "set the percent identity threshold" << endl;
		cerr << "                               \t" << "default value is 80" << endl;
		cerr << "  -m,  --max-hits              \t" << "ignore reads with more than a maximum number of hits" << endl;
		cerr << "                               \t" << "default value is " << options.maxHits << endl;
		cerr << "  -o,  --output FILE           \t" << "change output filename (default <READS FILE>.result)" << endl;
		cerr << "  -v,  --verbose               \t" << "verbose mode" << endl;
		cerr << "  -vv, --vverbose              \t" << "very verbose mode" << endl;
		cerr << "  -V,  --version               \t" << "print version number" << endl;
		cerr << "  -h,  --help                  \t" << "print this help" << endl;
		cerr << endl << "Output Format Options:" << endl;
		cerr << "  -a,  --alignment             \t" << "dump the alignment for each match" << endl;
		cerr << "  -of, --output-format NUM     \t" << "set output format" << endl;
		cerr << "                               \t" << "0 = Razer format (default, see README)" << endl;
		cerr << "                               \t" << "1 = enhanced Fasta format" << endl;
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
		cerr << "  -s,  --shape                 \t" << "set k-mer shape (binary string, default " << options.shape << ')' << endl;
		cerr << "  -t,  --threshold NUM         \t" << "set minimum k-mer threshold (default " << options.threshold << ')' << endl;
		cerr << "  -oc, --overabundance-cut NUM \t" << "set k-mer overabundance cut ratio (default " << options.abundanceCut << ')' << endl;
		cerr << "  -rl, --repeat-length NUM     \t" << "set simple-repeat length threshold (default " << options.repeatLength << ')' << endl;
		cerr << "  -tl, --taboo-length NUM      \t" << "set taboo length (default " << options.tabooLength << ')' << endl;
	} else {
		cerr << "Try 'razers --help' for more information." << endl;
	}
}

void printVersion() 
{
	cerr << "RazerS version 0.3 20080828 (prerelease)" << endl;
}

int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line

	RazerSOptions<>		options;
	unsigned			fnameCount = 0;
	const char			*fname[2] = { "", "" };

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
				printHelp(argc, argv, options);
				return 0;
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
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-a") == 0 || strcmp(argv[arg], "--alignment") == 0) {
				options.dumpAlignment = true;
				continue;
			}
			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
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
						if (options.outputFormat > 1)
							cerr << "Invalid output format options." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options);
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
				printHelp(argc, argv, options);
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
				printHelp(argc, argv, options);
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
				printHelp(argc, argv, options);
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
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-s") == 0 || strcmp(argv[arg], "--shape") == 0){
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
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
									printHelp(argc, argv, options);
									return 0;
							}
						if (ones == 0 || ones > 20) 
						{
							cerr <<	"Invalid Shape" << endl << endl;
							printHelp(argc, argv, options);
							return 0;
						}

						if (ones < 7 || ones > 14)
							cerr <<	"Warning: Shape should contain at least 7 and at most 14 '1's" << endl << endl;

						continue;				
					}
				}
				printHelp(argc, argv, options);
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
				printHelp(argc, argv, options);
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
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-t") == 0 || strcmp(argv[arg], "--threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
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
				printHelp(argc, argv, options);
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
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, options, true);
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
		} else {
			// parse file name
			if (fnameCount == 2) {
				printHelp(argc, argv, options);
				return 0;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 2) {
		if (options.printVersion)
			printVersion();
		else
			printHelp(argc, argv, options);
		return 0;
	}
	if (!options.forward && !options.reverse) { // enable both per default
		options.forward = true;
		options.reverse = true;
	}

	if (options.printVersion)
		printVersion();

	int result = mapReads(fname[0], fname[1], options);
	if (result == -1) 
	{
		cerr <<	"Invalid Shape" << endl << endl;
		printHelp(argc, argv, options);
		return 0;
	}
	return result;
}
