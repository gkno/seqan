///A tutorial about importing/exporting annotations
#include <fstream>
#include <iostream>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace std;
using namespace seqan;

int main(int argc, const char *argv[])
{
	CommandLineParser parser;

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "****************************");
	addTitleLine(parser, "*** Annotation Converter ***");
	addTitleLine(parser, "****************************");
	addHelpLine(parser, "");	
	addOption(parser, CommandLineOption("rg",  "read-gff",      "read annotation in GFF/GTF format", OptionType::String | OptionType::Label));
	addOption(parser, CommandLineOption("ru",  "read-ucsc",  2, "read annotation in UCSC format (knownGenes.txt knownIsoforms.txt)", OptionType::String | OptionType::Label));
	addHelpLine(parser, "");
	addOption(parser, CommandLineOption("wg",  "write-gff",     "write annotation in GFF format", OptionType::String | OptionType::Label));
	addOption(parser, CommandLineOption("wt",  "write-gtf",     "write annotation in GTF format", OptionType::String | OptionType::Label));
	addOption(parser, CommandLineOption("wu",  "write-ucsc",    "write annotation in UCSC format (knownGenes.txt knownIsoforms.txt)", OptionType::String | OptionType::Label));

	bool stop = !parse(parser, argc, argv, cerr);
	
	if (argc == 1 || stop)
	{
		if (argc == 1)
			shortHelp(parser, cerr);	// print short help and exit
		return 0;
	}

	FragmentStore<> store;
	String<char, CStyle> fileName;

	// IMPORT
	if (isSetLong(parser, "read-gff"))
	{
		getOptionValueLong(parser, "read-gff", fileName);
		ifstream file(toCString(fileName), ios_base::in | ios_base::binary);
		read(file, store, GFF());
	}
	if (isSetLong(parser, "read-ucsc"))
	{
		ifstream file1(toCString(getOptionValuesLong(parser, "read-ucsc")[0]), ios_base::in | ios_base::binary);
		ifstream file2(toCString(getOptionValuesLong(parser, "read-ucsc")[1]), ios_base::in | ios_base::binary);
		read(file1, store, UCSC());
		read(file2, store, UCSC());
	}

	// EXPORT
	if (isSetLong(parser, "write-gff"))
	{
		getOptionValueLong(parser, "write-gff", fileName);
		ofstream file(toCString(fileName), ios_base::out | ios_base::binary);
		write(file, store, GFF());
	}
	if (isSetLong(parser, "write-gtf"))
	{
		getOptionValueLong(parser, "write-gtf", fileName);
		ofstream file(toCString(fileName), ios_base::out | ios_base::binary);
		write(file, store, GTF());
	}
	if (isSetLong(parser, "write-ucsc"))
	{
		ofstream file1(toCString(getOptionValuesLong(parser, "write-ucsc")[0]), ios_base::in | ios_base::binary);
		ofstream file2(toCString(getOptionValuesLong(parser, "write-ucsc")[1]), ios_base::in | ios_base::binary);
		write(file1, store, UCSC());
		write(file2, store, UCSC_ISOFORMS());
	}
	
	return 0;
}
