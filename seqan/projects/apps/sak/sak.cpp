#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include "../library/apps/razers/mmap_fasta.h"

using namespace std;
using namespace seqan;



//____________________________________________________________________________
// Global Parameters

	int			optionSeqNo = 0;
	int			optionInfStart = -1;
	int			optionInfEnd = -1;
	bool		optionRevComp = false;
	const char	*optionOutput = NULL;					

	typedef Dna5 TAlphabet;
	typedef String<TAlphabet> TSeqString;

//____________________________________________________________________________


template<typename TChar>
inline bool
parse_isDigit(TChar const c)
{
	//return (((int ) c >=  48) && ((int) c <=  57));
	return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
		    (c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
}


template<typename TFile, typename TChar>
inline long double
parse_readDouble(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!parse_isDigit(c) && (c != '.')) break;
		append(str, c);
	}
 	return atof(toCString(str));
}

template<typename TFile, typename TChar>
inline void 
parse_skipWhitespace(TFile& file, TChar& c)
{
	if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) break;
	}
}



//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TSeqSet, typename TIDs>
bool loadSeqs(TSeqSet &seqs, TIDs &ids, const char *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
	split(multiFasta, Fasta());

	unsigned seqCount = length(multiFasta);
	resize(seqs, seqCount, Exact());
	resize(ids, seqCount, Exact());
	for(unsigned i = 0; i < seqCount; ++i)
	{
		assignSeq(seqs[i], multiFasta[i], Fasta());		// read Genome sequence
		assignSeqId(ids[i], multiFasta[i], Fasta());	// read Genome ids
	}

	return (seqCount > 0);
}

template < 
	typename TReadSet, 
	typename TReadIDs >
void saveFasta(
	TReadSet const &readSet,		// generated read sequences
	TReadIDs const &readIDs,		// corresponding Fasta ids
	string const &readFName)
{
	ostringstream fileName;
	if (*optionOutput != 0)
		fileName << optionOutput;
	else
		fileName << readFName << ".reads";

	ofstream file;
	
	file.open(fileName.str().c_str(), ios_base::out | ios_base::trunc);
	if (!file.is_open()) {
		cerr << "\nFailed to open output file" << endl;
		return;
	}
	else
		cout << "\nWriting reads to " << fileName.str() << "\n";

	unsigned reads = length(readSet);
	unsigned ids = length(readIDs);
	for(unsigned i = 0; i < reads; ++i)
	{
		if (i < ids)
			file << '>' << readIDs[i] << std::endl;
		else
			file << '>' << std::endl;
		file << readSet[i] << std::endl;
	}
	file.close();
}
/*
//temporary fake fastq output each base has quality 25
template < 
	typename TReadSet, 
	typename TReadIDs >
void saveFastq(
	TReadSet const &readSet,		// generated read sequences
	TReadIDs const &readIDs,		// corresponding Fasta ids
	string const &readFName)
{
	ostringstream fileName;
	if (*optionOutputFastq != 0)
		fileName << optionOutputFastq;
	else
		fileName << readFName << ".reads.fq";

	ofstream file;
	
	file.open(fileName.str().c_str(), ios_base::out | ios_base::trunc);
	if (!file.is_open()) {
		cerr << "\nFailed to open output file" << endl;
		return;
	}
	else
		cout << "\nWriting reads to " << fileName.str() << "\n";

	CharString toAscii;
	resize(toAscii,200);
	for (int i = -100; i < 100; i++) 
	{
	    	////convert to phred
        	////int q = (int) (10.0 * log(1.0 + pow(10.0,(double) (i / 10.0))) / log(10.0));
		int q = i;
		//get ascii character
	        toAscii[i+100] = (char)((q <= 93? q : 93) + 33);
	}

	unsigned reads = length(readSet);
	unsigned ids = length(readIDs);
	for(unsigned i = 0; i < reads; ++i)
	{
		if (i < ids)
			file << '@' << readIDs[i] << std::endl;
		else
			file << '@' << std::endl;
		file << readSet[i] << std::endl;
		file << "+" << std::endl;
		for(unsigned j=0; j<length(readSet[i]); ++j) file << toAscii[25+100];
		file << std::endl;
	}
	file.close();
}

*/

//////////////////////////////////////////////////////////////////////////////
// Print usage
void printHelp(int, const char *[], bool longHelp = false) 
{
	cerr << "************************" << endl;
	cerr << "*** Swiss Army Knife ***" << endl;
	cerr << "************************" << endl << endl;
	cerr << "Usage: sak [OPTION]... <SOURCE SEQUENCE FILE>" << endl;
	cerr << "\n";
	if (longHelp) {
		cerr << endl << "Main Options:" << endl;
		cerr << "  -o,  --output FILE            \t" << "set output filename (default: use stdout)" << endl;
		cerr << "  -h,  --help                   \t" << "print this help" << endl;
		cerr << endl << "Extract Options:" << endl;
		cerr << "  -sn, --seqno NUM              \t" << "select a sequence (default: select all)" << endl;
		cerr << "  -i,  --infix START END        \t" << "extract infix" << endl;
		cerr << "  -rc, --revcomp                \t" << "reverse complement" << endl;
	} else {
		cerr << "Try 'sak --help' for more information." << endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// Main part
int main(int argc, const char *argv[])
{
	unsigned fnameCount = 0;
	const char *fname[2] = { "" , "" };
	
	// Command line parsing
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option

			if (strcmp(argv[arg], "-sn") == 0 || strcmp(argv[arg], "--seqno") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionSeqNo;
					if (!istr.fail())
					{
						if (optionSeqNo < 0)
							cerr << "sequence number must be a value >=0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-i") == 0 || strcmp(argv[arg], "--infix") == 0) {
				if (arg + 2 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionInfStart;
					++arg;
					if (!istr.fail())
					{
						if (optionInfStart < 0)
							cerr << "infix start a value >=0" << endl << endl;
						else
						{
							istringstream istr(argv[arg]);
							if (!istr.fail())
							{
								istr >> optionInfEnd;
								if (optionInfEnd < 0)
									cerr << "infix end a value >=0" << endl << endl;
								else
									continue;
							}
						}
					}
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-rc") == 0 || strcmp(argv[arg], "--revcomp") == 0) {
				optionRevComp = true;
				continue;
			}

			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv);
					return 0;
				}
				++arg;
				optionOutput = argv[arg];
				continue;
			}

			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, true);
				return 0;
			}
		}
		else {
			// parse file name
			if (fnameCount == 1) {
				printHelp(argc, argv);
				return 0;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 1) {
		printHelp(argc, argv);
		return 0;
	}
	

	StringSet<TSeqString>	seqsIn, seqsOut;
	StringSet<CharString>	seqNamesIn, seqNamesOut;	// genome names, taken from the Fasta file

	if (!loadSeqs(seqsIn, seqNamesIn, fname[0]))
	{
		cerr << "Failed to open file" << fname[0] << endl;
		return 0;
	}

	cout << "\n"<<lengthSum(seqsIn) << " bps of " << length(seqsIn) << " source sequence loaded." << endl;

//____________________________________________________________________________

	if (optionSeqNo >= 0)
	{
		TSeqString tempSeq = seqsIn[optionSeqNo];
		CharString tempSeqName = seqNamesIn[optionSeqNo];
		clear(seqsIn);
		clear(seqNamesIn);
		appendValue(seqsIn, tempSeq);
		appendValue(seqNamesIn, tempSeqName);
	}

	resize(seqsOut, length(seqsIn));
	seqNamesOut = seqNamesIn;

	for (unsigned i = 0; i < length(seqsOut); ++i)
	{
		seqsOut[i] = infix(seqsIn[i], optionInfStart, optionInfEnd);
		if (optionRevComp)
			reverseComplementInPlace(seqsOut[i]);
	}

//____________________________________________________________________________

	if (optionOutput == NULL)
	{
		for (unsigned i = 0; i < length(seqsOut); ++i)
			cout << seqsOut[i] << endl;
	}
	else
		saveFasta(seqsOut, seqNamesOut, optionOutput);
	
//____________________________________________________________________________

	return 0;
}
