 /*==========================================================================
                     DFI - The Deferred Frequency Index
                   http://www.seqan.de/projects/dfi.html

 ============================================================================
  This is an application of the DFI algorithm in
  "Efficient string mining under constraints via the deferred frequency index"

 ============================================================================
  Copyright (C) 2008 by David Weese and Marcel H. Schulz

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#include <seqan/index.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////
// Custom frequency predicates

	static const double DFI_EPSILON = 0.0000001;

	// minfreq predicate for D0
	struct PredMinFreq 
	{	
		unsigned minFreq;

		template <typename TDataSet>
		PredMinFreq(unsigned _minFreq, TDataSet const &):
			minFreq(_minFreq) {}
			
		template <typename TDataSet>
		PredMinFreq(double _minSupp, TDataSet const &ds)
		{
			// emerging substring mode
			if (_minSupp * ds[1] < 1.0) {
				cerr << "Support must be at least 1/|db_1|... exit!" << endl;
				exit(1);
			}
			// adapt parameters from support to frequency
			minFreq = (unsigned) ceil((double) _minSupp * (ds[2] - ds[1]) - DFI_EPSILON);
		}
			
		inline bool operator()(_DFIEntry const &entry) const {
			return entry.freq[0] >= minFreq;
		}
	};

	// minsupp predicate for at least one dataset
	struct PredMinAllSupp
	{	
		double minSupp;
		String<int> dsLen;

		template <typename TDataSet>
		PredMinAllSupp(double _minSupp, TDataSet const &ds):
			minSupp(_minSupp) 
		{
			resize(dsLen, length(ds) - 1, Exact());
			for (unsigned i = 1; i < length(ds); ++i)
				dsLen[i - 1] = ds[i] - ds[i - 1];
		}
			
		inline bool operator()(_DFIEntry const &entry) const {
			for (unsigned i = 0; i < length(entry.freq); ++i)
				if (entry.freq[i] >= dsLen[i] * minSupp)
					return true;
			return false;
		}
	};

	// predicate for the Frequent Pattern Mining Problem
	struct PredFrequent 
	{	
		unsigned maxFreq;

		template <typename TDataSet>
		PredFrequent(unsigned _maxFreq, TDataSet const &):
			maxFreq(_maxFreq) {}
			
		inline bool operator()(_DFIEntry const &entry) const {
			return entry.freq[1] <= maxFreq;
		}
	};

	// predicate for the Emerging Substring Mining Problem
	struct PredEmerging
	{
		double growthRate;

		template <typename TDataSet>
		PredEmerging(double _growthRate, TDataSet const &ds) {
			growthRate = _growthRate * ((double) ds[1] / (double) (ds[2] - ds[1]) - DFI_EPSILON);
		}
			
		// HINT: here growthRate is frequency-related, not support-related
		inline bool operator()(_DFIEntry const &entry) const {
			return (double)entry.freq[0] >= (double)entry.freq[1] * growthRate;
		}
	};

	// predicate for the Maximum Entropy Problem
	struct PredEntropy
	{
		double maxEntropy;

		template <typename TDataSet>
		PredEntropy(double _maxEntropy, TDataSet const &):
			maxEntropy(_maxEntropy) {}
			
		inline bool operator()(_DFIEntry const &entry) const 
		{
			int sum = 0;
			double H = 0;

			for (unsigned i = 0; i < length(entry.freq); ++i)
				sum += entry.freq[i];
			
			double lSum = log((double)sum);					// sum cannot be zero
				
			for (unsigned i = 0; i < length(entry.freq); ++i)
				if (entry.freq[i])
				{
					double freq = entry.freq[i];
					H += freq * (log(freq) - lSum);
				}
			H /= -sum * log((double)length(entry.freq));	// normalize by datasets (divide by log m)
		
			return H <= maxEntropy;
		}
	};



//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta files
//
// seq ........ StringSet containing all sequences
// fileName ... array of database file names
// fileCount .. number of databases
// ds ......... seq[ds[i]],seq[ds[i]+1],...,seq[ds[i+1]-1] are the seqs. of database i
//
template <typename TSequences, typename TFileNames, typename TDatasets>
bool loadDatasets(
	TSequences		&seq, 
	TFileNames		const &fileNames, 
	TDatasets		&ds)
{
	// count sequences
	resize(ds, length(fileNames) + 1);
	unsigned seqCount = 0;
	for(unsigned s = 0; s < length(fileNames); ++s)
	{
		ds[s] = seqCount;
		ifstream file;
		file.open(fileNames[s], ios_base::in | ios_base::binary);
		if (!file.is_open()) return false;
		while (!_streamEOF(file)) {
			goNext(file, Fasta());
			++seqCount;
		}
	}
	ds[length(fileNames)] = seqCount;

	// import sequences
	resize(seq, seqCount);
	for(unsigned i = 0, s = 0; s < length(fileNames); ++s)	// for each database
	{
		ifstream file;
		file.open(fileNames[s], ios_base::in | ios_base::binary);
		if (!file.is_open()) return false;

		for(; (i < seqCount) && !_streamEOF(file); ++i)		// and each sequence
			read(file, seq[i], Fasta());					// read sequence
		file.close();
	}
	return (seqCount > 0);
}



//////////////////////////////////////////////////////////////////////////////
// Create DFI and output substrings within constraints band
//
// TPred .......... frequency predicate
// TPredHull ...... monotonic hull of the predicate above
// TAlphabet ...... sequence alphabet (Dna, Dna5, AminoAcid, char, ...)
//
// fileName ....... array of database file names
// paramPred ... .. parameter for the frequency predicate
// paramPredHull .. parameter for the monotonic hull
//
template <
	typename TPredHull, 
	typename TPred, 
	typename TAlphabet,
	typename TParamPredHull,
	typename TParamPred,
	typename TFileNames
>
int runDFI(
	TFileNames		const &fileNames,
	TParamPredHull	paramPredHull,
	TParamPred		paramPred)
{
	typedef String<TAlphabet, Alloc<> >						TString;
	typedef StringSet<TString>								TStringSet;
	typedef Index<TStringSet, Index_Wotd<
		WotdDFI<TPredHull, TPred> > >						TIndex;
	typedef Iter<TIndex, VSTree<TopDown<ParentLinks<> > > >	TIter;

	String<unsigned>	ds;
	TStringSet			mySet;

	if (!loadDatasets(mySet, fileNames, ds)) {
		cerr << "Database read error... exit!" << endl;
		return 1;
	}
	
/*	unsigned i=0;
	for(unsigned j=1;j<length(ds);++j)
	{
		cout<<"dataset "<<j<<':'<<ds[j]-ds[j-1]<<endl;
		for(;i<ds[j];++i)
			cout<<mySet[i]<<endl;
	}
*/
	TPred			pred(paramPred, ds);
	TPredHull		predHull(paramPredHull, ds);
	TIndex			index(mySet, predHull, pred);
	Pair<unsigned>	lPos;

	// set index partition of sequences into datasets
	index.ds = ds;
	TIter it(index);

	goBegin(it);
	while (!atEnd(it)) 
	{
        posLocalize(lPos, getOccurrence(it), stringSetLimits(container(it)));
		unsigned len = repLength(it);
//		for(unsigned l = parentRepLength(it) + 1; l <= len; ++l)
		unsigned l = len; 
		{
			cout << infix(
				mySet[getSeqNo(lPos)], 
				getSeqOffset(lPos),
				getSeqOffset(lPos) + l) << endl;
		}
		goNext(it);
	}

	return 0;
}

void printHelp(int, const char *[], bool longHelp = false)
{
	if (longHelp) {
		cerr << "**************************************************************" << endl;
		cerr << "***           DFI - The Deferred Frequency Index           ***" << endl;
		cerr << "*** (c) Copyright 2009 by David Weese and Marcel H. Schulz ***" << endl;
		cerr << "**************************************************************" << endl;
		cerr << endl;
		cerr << "Usage: dfi [OPTIONS] ... --minmax  <min_1> <max_2> <database 1> <database 2>" << endl;
		cerr << "       dfi [OPTIONS] ... --growth  <rho_s> <rho_g> <database 1> <database 2>" << endl;
		cerr << "       dfi [OPTIONS] ... --entropy <rho_s> <alpha> <database 1> <database 2> ... <database m>" << endl;
		cerr << endl;
		cerr << "Options:" << endl;
		cerr << "-f,  --minmax  \tsolve Frequent Pattern Mining Problem" << endl;
		cerr << "-g,  --growth  \tsolve Emerging Substring Mining Problem" << endl;
		cerr << "-e,  --entropy \tsolve Emerging Substring Mining Problem" << endl;
		cerr << endl;
		cerr << "-p,  --protein \tuse AminoAcid alphabet (for proteomes)" << endl;
		cerr << "-d,  --dna     \tuse DNA alphabet (for genomes)" << endl;
		cerr << "               \tThe default is byte alphabet" << endl;
	} else {
		cerr << "Try 'dfi --help' for more information." << endl;
	}
}

int main(int argc, const char *argv[])
{
	int optionAlphabet = 0;   // 0..char, 1..protein, 2..dna
	int optionPredicate = -1; // 0..minmax, 1..growth, 2..entropy
	unsigned optionMinFreq = 0;
	unsigned optionMaxFreq = 0;
	double optionMinSupp = 0;
	double optionGrowthRate = 0;
	double optionEntropy = 0;
	
	basic_string<const char*> fileNames;

	int arg = 1;
	for(; arg < argc; ++arg) 
	{
		if (argv[arg][0] == '-') 
		{
			// parse options

			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, true);
				return 0;
			}

			if (strcmp(argv[arg], "-f") == 0 || strcmp(argv[arg], "--minmax") == 0) 
			{
				optionPredicate = 0;
				if (arg + 2 < argc) 
				{
					istringstream istr1(argv[++arg]);
					istringstream istr2(argv[++arg]);
					istr1 >> optionMinFreq;
					istr2 >> optionMaxFreq;
					if (!istr1.fail() && !istr2.fail())
					{
						if (optionMinFreq < 1)
							cerr << "Minimum frequency threshold must be at least 1. Exit." << endl;
						else
						{
							if (optionMaxFreq < 1)
								cerr << "Maximum frequency threshold must be greater than 0. Exit." << endl;
							else
								continue;
						}
					}
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-g") == 0 || strcmp(argv[arg], "--growth") == 0) 
			{
				optionPredicate = 1;
				if (arg + 2 < argc) 
				{
					istringstream istr1(argv[++arg]);
					istringstream istr2(argv[++arg]);
					istr1 >> optionMinSupp;
					istr2 >> optionGrowthRate;
					if (!istr1.fail() && !istr2.fail())
					{
						if (optionMinSupp <= 0.0 || optionMinSupp > 1.0)
							cerr << "Support threshold must be greater than 0 and less than or equal to 1. Exit." << endl;
						else
						{
							if (optionGrowthRate < 1.0)
								cerr << "Growth rate must not be less than 1. Exit." << endl;
							else
								continue;
						}
					}
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-e") == 0 || strcmp(argv[arg], "--entropy") == 0) 
			{
				optionPredicate = 2;
				if (arg + 2 < argc) 
				{
					istringstream istr1(argv[++arg]);
					istringstream istr2(argv[++arg]);
					istr1 >> optionMinSupp;
					istr2 >> optionEntropy;
					if (!istr1.fail() && !istr2.fail())
					{
						if (optionMinSupp <= 0.0 || optionMinSupp > 1.0)
							cerr << "Support threshold must be greater than 0 and less than or equal to 1. Exit." << endl;
						else
						{
							if (optionEntropy <= 0.0 || optionEntropy > 1.0)
								cerr << "Entropy must not be grater than 0 and less or equal to 1. Exit." << endl;
							else
								continue;
						}
					}
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-p") == 0 || strcmp(argv[arg], "--protein") == 0) 
			{
				optionAlphabet = 1;
				continue;
			}

			if (strcmp(argv[arg], "-d") == 0 || strcmp(argv[arg], "--dna") == 0) 
			{
				optionAlphabet = 2;
				continue;
			}
		} else
			fileNames += argv[arg];
	}
	
	switch (optionPredicate)
	{
		case 0:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinFreq, PredFrequent, unsigned char> (fileNames, optionMinFreq, optionMaxFreq);
				case 1: return runDFI<PredMinFreq, PredFrequent, AminoAcid> (fileNames, optionMinFreq, optionMaxFreq);
				case 2: return runDFI<PredMinFreq, PredFrequent, Dna> (fileNames, optionMinFreq, optionMaxFreq);	
			}
		case 1:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinFreq, PredEmerging, unsigned char> (fileNames, optionMinSupp, optionGrowthRate);
				case 1: return runDFI<PredMinFreq, PredEmerging, AminoAcid> (fileNames, optionMinSupp, optionGrowthRate);
				case 2: return runDFI<PredMinFreq, PredEmerging, Dna> (fileNames, optionMinSupp, optionGrowthRate);
			}
		case 2:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinAllSupp, PredEntropy, unsigned char> (fileNames, optionMinSupp, optionEntropy);
				case 1: return runDFI<PredMinAllSupp, PredEntropy, AminoAcid> (fileNames, optionMinSupp, optionEntropy);
				case 2: return runDFI<PredMinAllSupp, PredEntropy, Dna> (fileNames, optionMinSupp, optionEntropy);
			}
	}
	printHelp(argc, argv);
	return 0;
}
