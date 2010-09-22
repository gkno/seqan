///A simple annotation converter. Convert a GFF to GTF or vice versa.
#include <fstream>
#include <iostream>
#include <string>

#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

int main(int argc, const char *argv[])
{
	typedef FragmentStore<> TFragStore;
	
	//////////////////////////////////////////////////////////////////////////////
	// Define options
	CommandLineParser parser;
	addUsageLine(parser, "[OPTION]... <infile> <outfile>");
	
	addOption(parser, CommandLineOption("gff",  "",    "write annotation in GFF format", OptionType::Bool));
	addOption(parser, CommandLineOption("gtf",  "",    "write annotation in GTF format (default)", OptionType::Bool));
	requiredArguments(parser, 2);

	bool stop = !parse(parser, argc, argv, std::cerr);
	if (stop) return 0;

	//////////////////////////////////////////////////////////////////////////////
	// Extract and check options	
	TFragStore store;
	std::ifstream inFile(toCString(getArgumentValue(parser, 0)), std::ios_base::in);
	std::ofstream outFile(toCString(getArgumentValue(parser, 1)), std::ios_base::out);

	if (!read(inFile, store, GFF()) && (stop = true))
		std::cerr << "Failed to load annotation." << std::endl;

	if (!stop)
		if (isSetLong(parser, "gff"))
		{
			if (!write(outFile, store, GFF()) && (stop = true))
				std::cerr << "Failed to write annotation." << std::endl;
		} else {
			if (!write(outFile, store, GTF()) && (stop = true))
				std::cerr << "Failed to write annotation." << std::endl;
		}
	
	if (stop)
	{
		cerr << "Exiting ..." << endl;
		return 1;
	}

	return 0;
}
