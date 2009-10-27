#include <fstream>
#include <iostream>
#include <sstream>


#include <seqan/basic/basic_debug.h> 	
#define SEQAN_PROFILE
#ifndef RELEASE	
			
//#define SEQAN_DEBUG			
//#define SEQAN_TEST	
	
					
#endif

#include <string>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/map.h>
#include <seqan/refinement.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/basic/basic_profile.h>

#include "base.h"
#include <seqan/store/store_align_intervals.h>
#include "read_gff.h"
#include <seqan/store/store_intervaltree.h>
#include "create_gff.h"
#include "overlap_module_o.h"

using namespace seqan;
using namespace std;


///////////////////////////////////////////////////////////////////////////////
inline void
test_createCombinations()
{
	StringSet<String<unsigned> > annoIds;
	StringSet<String<unsigned> > tuple;
	clear(annoIds);
	
	String<unsigned> str1;
	String<unsigned> str2;
	String<unsigned> str3;
	
	clear(str1);
	clear(str2);
	clear(str3);
	
	appendValue(str1, 1);
	appendValue(str1, 2);
	appendValue(str1, 3);
	appendValue(str2, 4);
	appendValue(str2, 5);
	appendValue(str3, 6);
	appendValue(str3, 7);
	appendValue(str3, 8);
	
	appendValue(annoIds, str1);
	appendValue(annoIds, str2);
	appendValue(annoIds, str3);
	
	createCombinations(tuple, annoIds);
	
	std::cout << "Tuple: " << std::endl;
	for (int i = 0; i < length(tuple); ++i)
	{
		for (int j = 0; j < length(getValue(tuple, i)); ++j)
		{
			std::cout << getValue(getValue(tuple, i), j);
		}
		std::cout << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

inline void
test_interSec()
{
	
	String<int> list1;
	String<int> list2;
	
	clear(list1);
	clear(list2);
	
	appendValue(list1, 3);
	appendValue(list1, 4);
	appendValue(list1, 4);
	appendValue(list1, 8);
	appendValue(list1, 9);
	appendValue(list1, 9);
	appendValue(list1, 12);
	appendValue(list1, 13);
	
	appendValue(list2, 4);
	appendValue(list2, 6);
	appendValue(list2, 9);
	appendValue(list2, 9);
	appendValue(list2, 10);
	appendValue(list2, 12);
	appendValue(list2, 12);
	appendValue(list2, 12);
	appendValue(list2, 13);
	appendValue(list2, 19);

	String<int> result;
	interSec(result, list1, list2);
	std::cout << "interSec: " << std::endl;
	for (int i = 0; i < length(result); ++i)
	{
		std::cout << getValue(result, i) << ", ";
	}
	std::cout << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
////// main
///////////////////////////////////////////////////////////////////////////////

int main( int argc, const char *argv[] ) 
{
	CharString nameSAM;
	CharString nameFASTA;
	CharString nameGFF;
	CharString outputPath;
	unsigned nTuple;
	unsigned offsetInterval;
	unsigned thresholdGaps;
	unsigned thresholdCount;
	double thresholdRPKM;
	bool maxTuple = 0;
	bool exact_nTuple = 0;
	bool unknownO = 0;
	
	CommandLineParser parser;
	
	addSection(parser, "Main Options:");
	addOption(parser, addArgumentText(CommandLineOption("s", "sam", "SAM-file with aligned reads", OptionType::String), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("f", "fa", "FASTA-file with contig sequence (only names, sequences should be empty)", OptionType::String), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("g", "gff", "GFF_file with annotations", OptionType::String), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("p", "outputPath", "path for output-files", (int)OptionType::String, ""), "<String>"));
	addOption(parser, addArgumentText(CommandLineOption("n", "nTuple", "nTuple", (int)OptionType::Int, 2), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("o", "offsetInterval", "offset to short alignment-intervals for search", (int)OptionType::Int, 5), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("t", "thresholdGaps", "threshold for allowed gaps in alignment (not introns)", (int)OptionType::Int, 5), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("c", "thresholdCount", "threshold for min. count of tuple for output", (int)OptionType::Int, 1), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("r", "thresholdRPKM", "threshold for min. RPKM of tuple for output", (double)OptionType::Double, 0.0), "<Double>"));
	addOption(parser, CommandLineOption("m", "maxTuple", "create only maxTuple", OptionType::Boolean));
	addOption(parser, CommandLineOption("e", "exact_nTuple", "create only Tuple of exact length n", OptionType::Boolean));
	addOption(parser, CommandLineOption("u", "unknown_orientation", "orientation of reads is unknown", OptionType::Boolean));
	
	if (!parse(parser, argc, argv, ::std::cerr)) return 1;
	
	getOptionValueLong(parser, "sam", nameSAM);
	getOptionValueLong(parser, "fa", nameFASTA);
	getOptionValueLong(parser, "gff", nameGFF);
	getOptionValueLong(parser, "outputPath", outputPath);
	getOptionValueLong(parser, "nTuple", nTuple);
	getOptionValueLong(parser, "offsetInterval", offsetInterval);
	getOptionValueLong(parser, "thresholdGaps", thresholdGaps);
	getOptionValueLong(parser, "thresholdCount", thresholdCount);
	getOptionValueLong(parser, "thresholdRPKM", thresholdRPKM);

	if (isSetLong(parser, "maxTuple")) maxTuple = 1;
	if (isSetLong(parser, "exact_nTuple")) exact_nTuple = 1;
	if (isSetLong(parser, "unknown_orientation")) unknownO = 1;
	
	if (maxTuple) 
	{
		nTuple = 0;		// sign for maxTuple
		exact_nTuple = 0;	// not both possible: maxTuple is prefered over exact_nTuple and n
	}
	
	ngsOverlapper(nameSAM, nameFASTA, nameGFF, outputPath, nTuple, exact_nTuple, thresholdGaps, offsetInterval, thresholdCount, thresholdRPKM, unknownO);
	

	return 0;
}



