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

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/index.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
using namespace seqan;

#define DEBUG_ENTROPY

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
			
		static inline double
		getEntropy(_DFIEntry const &entry)
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
			return H;
		}
			
		inline bool operator()(_DFIEntry const &entry) const 
		{
			return getEntropy(entry) <= maxEntropy;
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
		file.open(toCString(fileNames[s]), ios_base::in | ios_base::binary);
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
		file.open(toCString(fileNames[s]), ios_base::in | ios_base::binary);
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

#ifdef DEBUG_ENTROPY	
	// database lookup table
	typedef typename Infix< typename Fibre<TIndex, Fibre_SA>::Type const >::Type TOccs;
	typedef typename Iterator<TOccs, Standard>::Type TOccIter;
	String<unsigned>	dbLookup;
	String<bool>		seen;
	_DFIEntry			entry;
	PredEntropy			entrp(0, mySet);
	
	resize(dbLookup, length(mySet));
	resize(seen, length(mySet));
	resize(entry.freq, length(ds) - 1);
	for (unsigned d = 0, i = 0; i < length(mySet); ++i)
	{
		while (ds[d + 1] == i) ++d;
		dbLookup[i] = d;
	}
#endif

	TIter it(index);
	goBegin(it);
	while (!atEnd(it)) 
	{
        posLocalize(lPos, getOccurrence(it), stringSetLimits(container(it)));
		unsigned len = repLength(it);
//		for(unsigned l = parentRepLength(it) + 1; l <= len; ++l)
		unsigned l = len; 
		{
#ifdef DEBUG_ENTROPY
			// count frequencies (debug)
			TOccs occs = getOccurrences(it);
			TOccIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());
			arrayFill(begin(seen, Standard()), end(seen, Standard()), false);
			arrayFill(begin(entry.freq, Standard()), end(entry.freq, Standard()), 0);
			for (; oc != ocEnd; ++oc)
			{
				unsigned seqNo = getSeqNo(*oc, stringSetLimits(index));
				if (!seen[seqNo])
				{
					seen[seqNo] = true;
					++entry.freq[dbLookup[seqNo]];
				}
			}
				
			double H = entrp.getEntropy(entry);
			if (H <= 0.0) H = 0.0;
			cout << left << setw(14) << H << "[";
			for (unsigned i = 0; i < length(entry.freq); ++i)
				cout << right << setw(6) << entry.freq[i];
			cout << "]      ";
#endif
			cout << infix(
				mySet[getSeqNo(lPos)], 
				getSeqOffset(lPos),
				getSeqOffset(lPos) + l) << endl;
		}
		goNext(it);
	}

	return 0;
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
		
	CommandLineParser parser;
	string rev = "$Revision 0001 $";
	addVersionLine(parser, "DFI version 2.0 20090715 [" + rev.substr(11, 4) + "]");

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "**************************************************************");
	addTitleLine(parser, "***           DFI - The Deferred Frequency Index           ***");
	addTitleLine(parser, "*** (c) Copyright 2009 by David Weese and Marcel H. Schulz ***");
	addTitleLine(parser, "**************************************************************");
	addUsageLine(parser, "[OPTION]... --minmax  <min_1> <max_2> <database 1> <database 2>");
	addUsageLine(parser, "[OPTION]... --growth  <rho_s> <rho_g> <database 1> <database 2>");
	addUsageLine(parser, "[OPTION]... --entropy <rho_s> <alpha> <database 1> <database 2> ... <database m>");
	
	addOption(parser, CommandLineOption("f",  "minmax",  2, "solve Frequent Pattern Mining Problem", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("g",  "growth",  2, "solve Emerging Substring Mining Problem", OptionType::Double | OptionType::Label));
	addOption(parser, CommandLineOption("e",  "entropy", 2, "solve Entropy Mining Problem", OptionType::Double | OptionType::Label));
	addHelpLine(parser, "");
	addOption(parser, CommandLineOption("p",  "protein",    "use AminoAcid alphabet (for proteomes)", OptionType::Boolean));
	addOption(parser, CommandLineOption("d",  "dna",        "use DNA alphabet (for genomes)", OptionType::Boolean));
	addHelpLine(parser, "The default is byte alphabet");

	if (argc == 1)
	{
		shortHelp(parser, cerr);	// print short help and exit
		return 0;
	}

	bool stop = !parse(parser, argc, argv, cerr);

	//////////////////////////////////////////////////////////////////////////////
	// Extract options
	if (isSetLong(parser, "protein")) optionAlphabet = 1;
	if (isSetLong(parser, "dna")) optionAlphabet = 2;
	if (isSetLong(parser, "help")) return 0;	// print help and exit	
	
	//////////////////////////////////////////////////////////////////////////////
	// Check options
	if (isSetLong(parser, "minmax"))
	{
		optionPredicate = 0;
		getOptionValueLong(parser, "minmax", 0, optionMinFreq);
		getOptionValueLong(parser, "minmax", 1, optionMaxFreq);
		if ((optionMinFreq < 1) && (stop = true))
			cerr << "Minimum frequency threshold must be at least 1." << endl;
		if ((optionMaxFreq < 1) && (stop = true))
			cerr << "Maximum frequency threshold must be greater than 0." << endl;
		if ((argumentCount(parser) != 2) && (stop = true))
			cerr << "Please specify 2 databases." << endl;
	}
	if (isSetLong(parser, "growth"))
	{
		optionPredicate = 1;
		getOptionValueLong(parser, "growth", 0, optionMinSupp);
		getOptionValueLong(parser, "growth", 1, optionGrowthRate);
		if ((optionMinSupp <= 0.0 || optionMinSupp > 1.0) && (stop = true))
			cerr << "Support threshold must be greater than 0 and less than or equal to 1." << endl;
		if ((optionGrowthRate < 1.0) && (stop = true))
			cerr << "Growth rate must not be less than 1." << endl;
		if ((argumentCount(parser) != 2) && (stop = true))
			cerr << "Please specify 2 databases." << endl;
	}
	if (isSetLong(parser, "entropy"))
	{
		optionPredicate = 2;
		getOptionValueLong(parser, "entropy", 0, optionMinSupp);
		getOptionValueLong(parser, "entropy", 1, optionEntropy);
		if ((optionMinSupp <= 0.0 || optionMinSupp > 1.0) && (stop = true))
			cerr << "Support threshold must be greater than 0 and less than or equal to 1." << endl;
		if ((optionEntropy <= 0.0 || optionEntropy > 1.0) && (stop = true))
			cerr << "Entropy must not be grater than 0 and less or equal to 1." << endl;
		if ((argumentCount(parser) < 2) && (stop = true))
			cerr << "Please specify at least 2 databases." << endl;
	}

	if (stop)
	{
		cerr << "Exiting ..." << endl;
		return -1;
	}

	switch (optionPredicate)
	{
		case 0:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinFreq, PredFrequent, unsigned char> (getArgumentValues(parser), optionMinFreq, optionMaxFreq);
				case 1: return runDFI<PredMinFreq, PredFrequent, AminoAcid> (getArgumentValues(parser), optionMinFreq, optionMaxFreq);
				case 2: return runDFI<PredMinFreq, PredFrequent, Dna> (getArgumentValues(parser), optionMinFreq, optionMaxFreq);	
			}
		case 1:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinFreq, PredEmerging, unsigned char> (getArgumentValues(parser), optionMinSupp, optionGrowthRate);
				case 1: return runDFI<PredMinFreq, PredEmerging, AminoAcid> (getArgumentValues(parser), optionMinSupp, optionGrowthRate);
				case 2: return runDFI<PredMinFreq, PredEmerging, Dna> (getArgumentValues(parser), optionMinSupp, optionGrowthRate);
			}
		case 2:
			switch (optionAlphabet)
			{
				case 0: return runDFI<PredMinAllSupp, PredEntropy, unsigned char> (getArgumentValues(parser), optionMinSupp, optionEntropy);
				case 1: return runDFI<PredMinAllSupp, PredEntropy, AminoAcid> (getArgumentValues(parser), optionMinSupp, optionEntropy);
				case 2: return runDFI<PredMinAllSupp, PredEntropy, Dna> (getArgumentValues(parser), optionMinSupp, optionEntropy);
			}
	}
	cerr << "Please choose a mining problem." << endl;
	cerr << "Exiting ..." << endl;
	return -1;
}
