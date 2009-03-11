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

	int			optionSeqStart = 0;
	int			optionSeqEnd = SupremumValue<int>::VALUE;
	int			optionInfStart = 0;
	int			optionInfEnd = -1;
	bool		optionRevComp = false;
	const char	*optionOutput = NULL;
	bool		optionFastQ = false;

	typedef Dna5 TAlphabet;
	typedef String<TAlphabet> TSeqString;

//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TSeqSet, typename TQuals, typename TIDs>
bool loadSeqs(TSeqSet &seqs, TQuals &quals, TIDs &ids, const char *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;

	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);

	int seqCount = length(multiFasta);

	if (optionSeqEnd > seqCount) optionSeqEnd = seqCount;
	if (optionSeqStart < 0) optionSeqStart = 0;
		
	if (optionSeqStart < optionSeqEnd)
	{
		resize(seqs, optionSeqEnd - optionSeqStart, Exact());
		resize(ids, optionSeqEnd - optionSeqStart, Exact());
		resize(quals, optionSeqEnd - optionSeqStart, Exact());
		for(int i = 0, j = optionSeqStart; j < optionSeqEnd; ++i, ++j)
		{
			assignSeq(seqs[i], multiFasta[j], format);		// read Genome sequence
			assignSeqId(ids[i], multiFasta[j], format);		// read Genome ids
			assignQual(quals[i], multiFasta[j], format);	// read qualities
			fill(quals[i], length(seqs[i]), 33 + 80);
		}
	}

	return (seqCount > 0);
}

template < 
	typename TStream, 
	typename TSeq >
void dumpFastaSeq(
	TStream &out,
	TSeq &seq)
{
	unsigned size = length(seq);
	unsigned i;
	for (i = 0; i + 60 < size; i += 60)
		out << infix(seq, i, i + 60) << endl;
	out << infix(seq, i, size) << endl;
}


template < 
	typename TReadSet, 
	typename TReadIDs >
void saveFasta(
	TReadSet const &readSet,		// generated read sequences
	TReadIDs const &readIDs)		// corresponding Fasta ids
{
	ostream *out = &cout;
	ofstream file;
	
	if (optionOutput != NULL)
	{
		file.open(optionOutput, ios_base::out | ios_base::trunc);
		if (!file.is_open()) {
			cerr << "Failed to open output file" << endl;
			return;
		}
		else
			cout << "Writing reads to " << optionOutput << "\n";
		out = &file;
	}

	unsigned reads = length(readSet);
	for(unsigned i = 0; i < reads; ++i)
	{
		(*out) << '>' << readIDs[i] << std::endl;
		dumpFastaSeq(*out, readSet[i]);
	}
	
	file.close();
}

template < 
	typename TReadSet,
	typename TQualSet,
	typename TReadIDs >
void saveFastq(
	TReadSet const &readSet,		// generated read sequences
	TQualSet const &qualSet,		// qualities
	TReadIDs const &readIDs)		// corresponding Fasta ids
{
	ostream *out = &cout;
	ofstream file;
	
	if (optionOutput != NULL)
	{
		file.open(optionOutput, ios_base::out | ios_base::trunc);
		if (!file.is_open()) {
			cerr << "Failed to open output file" << endl;
			return;
		}
		else
			cout << "Writing reads to " << optionOutput << "\n";
		out = &file;
	}

	unsigned reads = length(readSet);
	for(unsigned i = 0; i < reads; ++i)
	{
		(*out) << '@' << readIDs[i] << std::endl;
		dumpFastaSeq(*out, readSet[i]);
		(*out) << '+' << readIDs[i] << std::endl;
		dumpFastaSeq(*out, qualSet[i]);
	}
	
	file.close();
}

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
		cerr << "  -q,  --qual                   \t" << "enable Fastq output (default: Fasta)" << endl;
		cerr << "  -h,  --help                   \t" << "print this help" << endl;
		cerr << endl << "Extract Options:" << endl;
		cerr << "  -s,  --sequence NUM           \t" << "select a single sequence" << endl;
		cerr << "  -ss, --sequences START END    \t" << "select sequences (default: select all)" << endl;
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

			if (strcmp(argv[arg], "-s") == 0 || strcmp(argv[arg], "--sequence") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionSeqStart;
					if (!istr.fail())
					{
						if (optionSeqStart < 0)
							cerr << "sequence number must be a value >=0" << endl << endl;
						else
						{
							optionSeqEnd = optionSeqStart + 1;
							continue;
						}
					}
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-ss") == 0 || strcmp(argv[arg], "--sequences") == 0) {
				if (arg + 2 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionSeqStart;
					++arg;
					if (!istr.fail())
					{
						if (optionSeqStart < 0)
							cerr << "first sequence number must be a value >=0" << endl << endl;
						else
						{
							istringstream istr(argv[arg]);
							if (!istr.fail())
							{
								istr >> optionSeqEnd;
								if (optionSeqEnd < 0)
									cerr << "last sequence number must be a value >=0" << endl << endl;
								else
								{
									++optionSeqEnd;
									continue;
								}
							}
						}
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

			if (strcmp(argv[arg], "-q") == 0 || strcmp(argv[arg], "--qual") == 0) {
				optionFastQ = true;
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
	
//____________________________________________________________________________
// input

	StringSet<TSeqString>	seqsIn, seqsOut;
	StringSet<CharString>	qualsIn, qualsOut;
	StringSet<CharString>	seqNamesIn, seqNamesOut;	// genome names, taken from the Fasta file

	if (!loadSeqs(seqsIn, qualsIn, seqNamesIn, fname[0]))
	{
		cerr << "Failed to open file" << fname[0] << endl;
		return 0;
	}
//	cout << lengthSum(seqsIn) << " bps of " << length(seqsIn) << " source sequence loaded." << endl;

//____________________________________________________________________________
// data processing

	resize(seqsOut, length(seqsIn));
	for (unsigned i = 0; i < length(seqsIn); ++i)
	{
		unsigned end = optionInfEnd;
		if (end > length(seqsIn[i]))
			end = length(seqsIn[i]);

		if (optionInfStart < (int)length(seqsIn[i]))
		{
			seqsOut[i] = infix(seqsIn[i], optionInfStart, end);
			if (optionRevComp)
				reverseComplementInPlace(seqsOut[i]);
		}
	}
	seqNamesOut = seqNamesIn;
	qualsOut = qualsIn;

//____________________________________________________________________________
// output

	if (optionFastQ)
		saveFastq(seqsOut, qualsOut, seqNamesOut);
	else
		saveFasta(seqsOut, seqNamesOut);
	
//____________________________________________________________________________

	return 0;
}
