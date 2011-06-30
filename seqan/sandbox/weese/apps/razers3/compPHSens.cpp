#include <fstream>
#include <iostream>
#include <string>

#include <seqan/misc/misc_cmdparser.h>
#include "razers.h"

using namespace seqan;

int main(int argc, const char *argv[])
{

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	typedef RazerSOptions<> TOptions;
	CommandLineParser parser;
	TOptions options;

	addUsageLine(parser, "[OPTION]... <error_dist file>");

	unsigned maxOverlap = 6;
	unsigned maxErrors = 12;
	unsigned numReads = 1000000;

	addOption(parser, CommandLineOption("n",   "num-reads",    "number of reads", OptionType::Int | OptionType::Label, numReads));
	addOption(parser, CommandLineOption("mo",  "max-overlap",  "estimate for overlaps 0,1,...,max-overlap", OptionType::Int | OptionType::Label, maxOverlap));
	addOption(parser, CommandLineOption("me",  "max-errors",   "estimate for errors 0,1,...,max-errors", OptionType::Int | OptionType::Label, maxErrors));
	requiredArguments(parser, 1);

	bool stop = !parse(parser, argc, argv, std::cerr);
	if (stop) return 0;

	//////////////////////////////////////////////////////////////////////////////
	// Extract and check options
    getOptionValueLong(parser, "num-reads", numReads);
    getOptionValueLong(parser, "max-overlap", maxOverlap);
    getOptionValueLong(parser, "max-errors", maxErrors);

	std::ifstream errorDistFile(toCString(getArgumentValue(parser, 0)), std::ios_base::in | std::ios_base::binary);
	double errorProb;
	
	while (true)
	{
		errorDistFile >> errorProb;
		if (!errorDistFile.good()) break;
		appendValue(options.errorProb, errorProb);
	}
	resize(options.readLengths, length(options.errorProb), 0);
	appendValue(options.readLengths, numReads);
	options.lossRate = 1.0;
	
//	for(unsigned k=0;k<length(options.errorProb);++k)
//		std::cout<<"  "<<options.errorProb[k];
		
	//////////////////////////////////////////////////////////////////////////////
	// Output losses to file
	
	unsigned errorLimit = 10;
	std::cout << "error_dist_file\terrors\toverlap\tq";
	for (unsigned i = 0; i <= errorLimit; ++i)
		std::cout << '\t' << i << "-error_loss";
	std::cout << std::endl;
	
	//////////////////////////////////////////////////////////////////////////////
	// Compute losses for different errors, overlaps
	
	Shape<Dna, GenericShape> dummyShape;
	unsigned maxLength = length(options.readLengths) - 1;
	for (unsigned errors = 0; errors <= maxErrors; ++errors)
	{
		options.errorRate = (double)errors / (double)maxLength + 0.00001;
		
		typedef TOptions::TProb TFloat;
		String<TFloat> estLosses;
		estimatePigeonholeLosses(estLosses, _pigeonholeMaxShapeWeight(dummyShape), options);

		for (unsigned idx = 3 + errors; idx < length(estLosses); idx += 3 + errors)
		{
			unsigned ol = estLosses[idx];
			unsigned q  = estLosses[idx + 1];

			std::cout << getArgumentValue(parser, 0) << '\t' << errors << '\t' << ol << '\t' << q;
				
			for (unsigned e = 0; e <= errorLimit; ++e)
				if (e <= errors)
					std::cout << '\t' << estLosses[idx + 2 + e];
			
			std::cout << std::endl;
		}
	}

	return 0;
}
