#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace seqan;

struct ReadEntry
{
	unsigned readId;
    string chromosome;
    char strand;
    String<__int64> positions;
};

void daveSuessSauer(FragmentStore<> &store, std::ofstream &fasta_file_stream, const ReadEntry & read)
{
	unsigned contigId;
	CharString contigName = read.chromosome;
	if (!getIdByName(store.contigNameStore, contigName, contigId))
	{
		cerr << "contig " << read.chromosome << " not found." << endl;
		return;
	}
	Dna5String readSeq;
	unsigned half = length(read.positions) / 2;
	for (unsigned i = 0; i < half; ++i)
		append(readSeq, infix(store.contigStore[contigId].seq, read.positions[i], read.positions[i+half]));
	if (read.strand == '-')
		reverseComplementInPlace(readSeq);
	stringstream sstr;
	sstr << read.readId;
	CharString id = sstr.str();
	write(fasta_file_stream, readSeq, id, Fasta());
}

int parseChinaReads(FragmentStore<> &store, string infile_name, string outfile_name)
{
	ReadEntry read;
	reserve(read.positions, 10);
	char tmp;
	long unsigned pos;
	string line;
	ifstream infile(infile_name.c_str());
	ofstream outfile(outfile_name.c_str());
	
	if(!infile.is_open())
	{
		cerr<<"cannot open input file: "<<infile_name<<endl;
		return 1;
	}
	
	if(!outfile.is_open())
	{
		cerr<<"cannot open output file: "<<outfile_name<<endl;
		return 1;
	}

	while(!infile.eof())
	{
		clear(read.positions);
		getline(infile, line);
		stringstream line_s(line);
		line_s.seekg(3,ios::beg);
		line_s>>read.chromosome>>read.strand;
		while(line_s>>pos)
		{
			appendValue(read.positions, pos);
			line_s>>tmp;
		}        
		daveSuessSauer(store, outfile, read);
	}
	infile.close();
	outfile.close();
	return 0;
}

int main (int argc, char const * argv[])
{
	CommandLineParser	parser;
	FragmentStore<>		store;			// stores all of the tables

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "*****************************************");
	addTitleLine(parser, "***     Chop Suey Read Generator      ***");
	addTitleLine(parser, "*** (c) Copyright 2010 by David Weese ***");
	addTitleLine(parser, "*****************************************");
	addUsageLine(parser, "[OPTION]... <genome file> <china file> <output file>");
	
	addHelpLine(parser, "");
	
	if (argc == 1 || !parse(parser, argc, argv, cerr))
	{
		shortHelp(parser, cerr);	// print short help and exit
		return 0;
	}
	
	unsigned argCount = argumentCount(parser);
	if (argCount < 2 || !loadContigs(store, getArgumentValues(parser)[0]))
		return 1;
	
	parseChinaReads(store, toCString(getArgumentValues(parser)[1]), toCString(getArgumentValues(parser)[2]));
	
}
