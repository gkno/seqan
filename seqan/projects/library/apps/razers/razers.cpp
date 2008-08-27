 /*==========================================================================
                     RazerS - Fast Mapping of Short Reads
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

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

#define SEQAN_PROFILE			// enable time measuring
//#define SEQAN_DEBUG_SWIFT		// test SWIFT correctness and print bucket parameters
//#define RAZERS_PROFILE		// omit dumping results
//#define RAZERS_DEBUG			// print verification regions
#define RAZERS_PRUNE_QGRAM_INDEX
//#define RAZERS_HAMMINGVERIFY	// allow only mismatches, no indels

#include "seqan/platform.h"
#ifdef PLATFORM_WINDOWS
	#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
	#define SEQAN_DEFAULT_TMPDIR "./"
#endif

#include <iostream>
#include <fstream>
#include <seqan/pipe.h>
#include <seqan/find.h>
#include <seqan/align.h>
#include <seqan/index.h>
#include <seqan/find/find_swift.h>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////
// Fixed parameters

/*
	// Text setup (replace the 1-hull enumerator with the standard swift)
	typedef	unsigned char	TAlphabet;
	const int				QGRAM_LEN = 3;
	#define OMIT_REVERSECOMPLEMENT
*/
	// DNA setup
	typedef	Dna				TAlphabet;
#ifndef QGRAM_LEN
	const int				QGRAM_LEN = 11;
#endif


//////////////////////////////////////////////////////////////////////////////
// Default options

	static bool			optionForward = false;				// compute forward oriented read matches
	static bool			optionRev = false;					// compute reverse oriented read matches
	static double		optionErrorRate = 0.2;				// Criteria 1 threshold
	static bool			optionDumpAlignment = false;		// compute and dump the match alignments in the result files
	static const char	*optionOutput = "";					// prefix of output files (e.g. "results/run02.", default: "")
	static unsigned		optionOutputFormat = 0;				// 0..Razer format
															// 1..enhanced Fasta
	static unsigned		optionGenomeNaming = 0;				// 0..use Fasta id
															// 1..enumerate reads beginning with 1
	static unsigned		optionReadNaming = 0;				// 0..use Fasta id
															// 1..enumerate reads beginning with 1
															// 2..use the read sequence (only for short reads!)
	static unsigned		optionSortOrder = 0;				// 0..sort keys: 1. read number, 2. genome position
															// 1..           1. genome position, 2. read number
	static unsigned		optionPositionFormat = 0;			// 0..gap space
															// 1..position space
	static string		optionShape = "11111111111";		// shape (e.g. 11111111111)
	static int			optionThreshold = 1;				// threshold
	static int			optionTabooLength = 1;				// taboo length
	static int			optionRepeatLength = 1000;			// repeat length threshold
	static double		optionAbundanceCut = 1;				// abundance threshold
	static int			_debugLevel = 0;					// level of verbosity
	static bool			optionPrintVersion = false;			// print version number


//////////////////////////////////////////////////////////////////////////////
// Typedefs

	// definition of a Read match
	template <typename _TGPos>
	struct ReadMatch {
		typedef _TGPos TGPos;
		unsigned gseqNo;		// genome seqNo
		unsigned rseqNo;		// read seqNo
		TGPos	gBegin;			// begin position of the match in the genome
		TGPos	gEnd;			// end position of the match in the genome
		short	editDist;		// Levenshtein distance
		char	orientation;	// 'F'..forward strand, 'R'..reverse comp. strand
	};

	// less-operators (to sort matches and remove duplicates)
	template <typename TReadMatch>
	struct LessGPosRNo : public binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.rseqNo < b.rseqNo) return true;
			return a.editDist < b.editDist;
		}
	};

	template <typename TReadMatch>
	struct LessRNoGPos : public binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			return a.editDist < b.editDist;
		}
	};


//////////////////////////////////////////////////////////////////////////////
// Definitions

	typedef String<Dna5>		TGenome;
	typedef StringSet<TGenome>	TGenomeSet;
	typedef String<TAlphabet>	TRead;
	typedef StringSet<TRead>	TReadSet;

namespace seqan 
{
	template <typename TShape>
	struct SAValue< Index<TReadSet, TShape> > {
		typedef Pair<
			unsigned,		// many
			unsigned,		// short reads
			Compressed
		> Type;
	};
}

	typedef ReadMatch<Difference<TGenome>::Type>	TMatch;			// a single match
	typedef String<TMatch/*, Block<>*/ >			TMatches;		// array of matches

#ifdef RAZERS_PRUNE_QGRAM_INDEX
namespace seqan 
{
	//////////////////////////////////////////////////////////////////////////////
	// Repeat masker
	template <typename TShape>
	inline bool _qgramDisableBuckets(Index<TReadSet, Index_QGram<TShape> > &index) 
	{
		typedef Index<TReadSet, Index_QGram<TShape>	>		TReadIndex;
		typedef typename Fibre<TReadIndex, QGram_Dir>::Type	TDir;
		typedef typename Iterator<TDir, Standard>::Type		TDirIterator;
		typedef typename Value<TDir>::Type					TSize;

		TDir &dir    = indexDir(index);
		bool result  = false;
		unsigned counter = 0;
		TSize thresh = (TSize)(length(index) * optionAbundanceCut);
		if (thresh < 100) thresh = 100;

		TDirIterator it = begin(dir, Standard());
		TDirIterator itEnd = end(dir, Standard());
		for (; it != itEnd; ++it)
			if (*it > thresh) {
//				String<TAlphabet> qgram;
//				unhash(qgram, it - begin(dir, Standard()), QGRAM_LEN);
//				cerr << qgram << " " << *it << endl;
				*it = (TSize)-1;
				result = true;
				++counter;
			}

		if (counter > 0 && _debugLevel >= 1)
			cerr << "Removed " << counter << " k-mers" << endl;


		return result;
	}
}
#endif

//////////////////////////////////////////////////////////////////////////////
// Load multi-Fasta sequences
template <typename TReadSet, typename TNameSet>
bool loadFasta(TReadSet &reads, TNameSet &fastaIDs, char const *fileName)
{
	// count sequences
	unsigned seqCount = 0;

	ifstream file;
	file.open(fileName, ios_base::in | ios_base::binary);
	if (!file.is_open()) return false;
	while (!_streamEOF(file)) {
		goNext(file, Fasta());
		++seqCount;
	}

	// import sequences
	file.clear();
	file.seekg(0, ios_base::beg);
	resize(fastaIDs, seqCount);
	resize(reads, seqCount);
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 
	{
		readID(file, fastaIDs[i], Fasta());		// read Fasta id
		read(file, reads[i], Fasta());			// read Read sequence
	}
	file.close();
	return (seqCount > 0);
}


//////////////////////////////////////////////////////////////////////////////
// Find read matches
template <typename TMatches, typename TGenomeSet, typename TReadIndex>
void findReads(
	TMatches &matches,			// resulting matches
	TGenomeSet &genomes,		// Genome
	TReadIndex &readIndex,
	char orientation)			// q-gram index of Reads
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Value<TGenomeSet>::Type				TGenome;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Infix<TGenome>::Type					TGenomeInfix;
	typedef typename Value<TReadSet>::Type					TRead;
	typedef typename Value<TMatches>::Type					TMatch;

	// FILTRATION
	
	typedef Finder<TGenome, Swift<> >						TSwiftFinder;
	typedef Pattern<TReadIndex, Swift<> >					TSwiftPattern;

	// VERIFICATION

	// find read match end
	typedef Finder<TGenomeInfix>							TMyersFinder;
#ifdef RAZERS_HAMMINGVERIFY
	typedef Score<int>										TScore;
	typedef Pattern<TRead, DPSearch<TScore> >				TMyersPattern;
#else
	typedef Pattern<TRead, MyersUkkonen>					TMyersPattern;
#endif

	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>	TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>			TReadRev;
	typedef Finder<TGenomeInfixRev>						TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>		TMyersPatternRev;

	TSwiftPattern swiftPattern(readIndex);
	String<TMyersPattern> forwardPatterns;

	swiftPattern.params.minThreshold = optionThreshold;
	swiftPattern.params.tabooLength = optionTabooLength;

#ifdef RAZERS_HAMMINGVERIFY
	TScore	scoreType = Score<int>(0, -1, -1001, -1000); // levenshtein-score (match, mismatch, gapOpen, gapExtend)
#endif

	__int64 TP = 0;
	__int64 FP = 0;
	SEQAN_PROTIMESTART(find_time);

	// init forward verifiers
	unsigned readCount = countSequences(readIndex);
	resize(forwardPatterns, readCount);
	for(unsigned i = 0; i < readCount; ++i)
	{
#ifdef RAZERS_HAMMINGVERIFY
		setHost(forwardPatterns[i], indexText(readIndex)[i]);
		setScoringScheme(forwardPatterns[i], scoreType);
#else
		setHost(forwardPatterns[i], indexText(readIndex)[i]);
#endif
	}

	// iterate all genomic sequences
	for(unsigned hstkSeqNo = 0; hstkSeqNo < length(genomes); ++hstkSeqNo) 
	{
		if (_debugLevel >= 1)
			cerr << endl << "Process genome seq #" << hstkSeqNo;
		TGenome &genome = genomes[hstkSeqNo];
		TSwiftFinder swiftFinder(genome, optionRepeatLength, 1);

		// iterate all verification regions returned by SWIFT
		while (find(swiftFinder, swiftPattern, optionErrorRate, _debugLevel)) 
		{
			unsigned ndlSeqNo = (*swiftFinder.curHit).ndlSeqNo;
			unsigned ndlLength = sequenceLength(ndlSeqNo, readIndex);

			TGenomeInfix inf(range(swiftFinder, genome));
			TMyersFinder myersFinder(inf);
			TMyersPattern &myersPattern = forwardPatterns[ndlSeqNo];
#ifdef RAZERS_DEBUG
			cout<<"Verify: "<<endl;
			cout<<"Genome: "<<inf<<"\t" << beginPosition(inf) << "," << endPosition(inf) << endl;
			cout<<"Read:   "<<host(myersPattern)<<endl;
#endif
			// find end of best semi-global alignment
			int maxScore = InfimumValue<int>::VALUE;
			int minScore = -(int)(ndlLength * optionErrorRate);
			TMyersFinder maxPos;
			while (find(myersFinder, myersPattern, minScore))
				if (maxScore < getScore(myersPattern)) 
				{
					maxScore = getScore(myersPattern);
					maxPos = myersFinder;
				}
			
			if (maxScore >= minScore) 
			{
				TMatch m = {
					hstkSeqNo,
					ndlSeqNo,
					beginPosition(swiftFinder),
					beginPosition(swiftFinder) + position(maxPos) + 1,
					-maxScore,
					orientation
				};
				if (m.gBegin < 0) m.gBegin = 0;
				
#ifdef RAZERS_HAMMINGVERIFY
				m.gBegin = m.gEnd - ndlLength;
#else
				TGenomeInfixRev		infRev(infix(genome, m.gBegin, m.gEnd));
				TReadRev			readRev(indexText(readIndex)[ndlSeqNo]);
				TMyersFinderRev		myersFinderRev(infRev);
				TMyersPatternRev	myersPatternRev(readRev);

				// find beginning of best semi-global alignment
				if (find(myersFinderRev, myersPatternRev, maxScore))
					m.gBegin = m.gEnd - (position(myersFinderRev) + 1);
				else {
					// this case should never occur
					cerr << "1GENOME: " << host(myersFinder) << endl;
					cerr << "1READ:   " << indexText(readIndex)[ndlSeqNo] << endl;
					cerr << "2GENOME: " << infix(genome, m.gBegin, m.gEnd) << '\t' << m.gBegin << ',' << m.gEnd << endl;
					cerr << "2READ:   " << indexText(readIndex)[ndlSeqNo] << endl;
					cerr << "3GENOME: " << infRev << endl;
					cerr << "3READ:   " << readRev << endl;
					cerr << "HUH?" << endl;
				}
#endif

#ifndef RAZERS_PROFILE
				// transform coordinates to the forward strand
				if (orientation == 'R') 
				{
					TSize gLength = length(genome);
					TSize temp = m.gBegin;
					m.gBegin = gLength - m.gEnd;
					m.gEnd = gLength - temp;
				}
				appendValue(matches, m);
#endif

				++TP;
	/*			cerr << "\"" << range(swiftFinder, genomeInf) << "\"  ";
				cerr << hstkPos << " + ";
				cerr << endl;
	*/		} else {
				++FP;
	/*			cerr << "\"" << range(swiftFinder, genomeInf) << "\"   \"" << range(swiftPattern) << "\"  ";
				cerr << ndlSeqNo << " : ";
				cerr << hstkPos << " + ";
				cerr << bucketWidth << "  " << TP << endl;
	*/		}
		}
	}
	if (_debugLevel >= 1)
		cerr << endl << "Finding reads took               \t" << SEQAN_PROTIMEDIFF(find_time) << " seconds" << endl;

	if (_debugLevel >= 2) {
		cerr << endl;
		cerr << "___FILTRATION_STATS____" << endl;
		cerr << "Swift FP: " << FP << endl;
		cerr << "Swift TP: " << TP << endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// Dump an alignment
template <typename TFile, typename TSource, typename TSpec>
inline void
dumpAlignment(TFile & target, Align<TSource, TSpec> const & source)
{
SEQAN_CHECKPOINT
	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;

	TRowsPosition row_count = length(rows(source));
	TPosition begin_ = beginPosition(cols(source));
	TPosition end_ = endPosition(cols(source));
	
	// Print sequences
	for(TRowsPosition i=0;i<row_count;++i) {
		if (i == 0)
			_streamWrite(target, "#Read:   ");
		else
			_streamWrite(target, "#Genome: ");
		TRow& row_ = row(source, i);
		typedef typename Iterator<typename Row<TAlign>::Type const>::Type TIter;
		TIter begin1_ = iter(row_, begin_);
		TIter end1_ = iter(row_, end_);
		for (; begin1_ != end1_; ++begin1_) {
			if (isGap(begin1_)) _streamPut(target, gapValue<char>());
			else _streamPut(target, *begin1_);
		}
		_streamPut(target, '\n');
	}
}


//////////////////////////////////////////////////////////////////////////////
// Output matches
template <
	typename TMatchSet,
	typename TGenome,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames
>
void dumpMatches(
	TMatchSet &matches,				// forward/reverse matches
	TGenome const &genomes,			// Genome sequences
	TGenomeNames const &genomeIDs,	// Read names (read from Fasta file, currently unused)
	TReads const &reads,			// Read sequences
	TReadNames const &readIDs,		// Read names (read from Fasta file, currently unused)
	string const &genomeFName,		// genome name (e.g. "hs_ref_chr1.fa")
	string const &readFName)		// read name (e.g. "reads.fa")
{
	typedef typename Value<TMatchSet>::Type		TMatches;
	typedef typename Value<TMatches>::Type		TMatch;
	typedef typename Value<TReads>::Type		TRead;
	typedef typename TMatch::TGPos				TGPos;

	SEQAN_PROTIMESTART(dump_time);

	// how many 0's should be padded?
	int pzeros = 0;
	for (unsigned l = length(reads); l > 9; l = l / 10)
		++pzeros;

	int gzeros = 0;
	for (unsigned l = length(genomes); l > 9; l = l / 10)
		++gzeros;


	// remove the directory prefix of genomeFName and readFName
	size_t lastPos = genomeFName.find_last_of('/') + 1;
	if (lastPos == genomeFName.npos) lastPos = genomeFName.find_last_of('\\') + 1;
	if (lastPos == genomeFName.npos) lastPos = 0;
	string genomeName = genomeFName.substr(lastPos);

	lastPos = readFName.find_last_of('/') + 1;
	if (lastPos == readFName.npos) lastPos = readFName.find_last_of('\\') + 1;
	if (lastPos == readFName.npos) lastPos = 0;
	string readName = readFName.substr(lastPos);
	

	Align<TRead, ArrayGaps> align;
#ifdef RAZERS_HAMMINGVERIFY
	Score<int> scoreType = Score<int>(0, -1, -1001, -1000);		// levenshtein-score (match, mismatch, gapOpen, gapExtend)
#else
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
#endif
	resize(rows(align), 2);

	bool multipleGenomes = countSequences(genomes) > 1;
	ofstream file;

	ostringstream fileName;
	if (*optionOutput != 0)
		fileName << optionOutput;
	else
		fileName << readFName << ".result";

	file.open(fileName.str().c_str(), ios_base::out | ios_base::trunc);
	if (!file.is_open()) {
		cerr << "Failed to open output file" << endl;
		return;
	}

	switch (optionSortOrder) {
		case 0:
			sort(
				begin(matches, Standard()),
				end(matches, Standard()), 
				LessRNoGPos<TMatch>());
			break;

		case 1:
			sort(
				begin(matches, Standard()),
				end(matches, Standard()), 
				LessGPosRNo<TMatch>());
			break;
	}

	typename Iterator<TMatches, Standard>::Type it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type itEnd = end(matches, Standard());

	typename TMatch::TGPos _gBegin = -1;
	unsigned gseqNo = -1;
	unsigned readNo = -1;
	char orientation;
	typename TMatch::TGPos gBegin, gEnd;


	switch (optionOutputFormat) 
	{
		case 0:	// Razer Format
			for(; it != itEnd; ++it) 
			{
				if (readNo == (*it).rseqNo && _gBegin == (*it).gBegin && gseqNo == (*it).gseqNo)
					continue;

				readNo = (*it).rseqNo;
				gseqNo = (*it).gseqNo;
				_gBegin = (*it).gBegin;
				gEnd = (*it).gEnd;
				orientation = (*it).orientation;

				unsigned	readLen = length(reads[readNo]);
				double		percId;
		
				gBegin = _gBegin;

				if (optionPositionFormat == 1)
					++gBegin;

				percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);

				switch (optionReadNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << readIDs[readNo];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << readName << '#' << setw(pzeros) << readNo + 1;
						break;

					// 2..filename is the read sequence itself
					case 2:
						file << reads[readNo];
				}

				file << ',' << optionPositionFormat << ',' << readLen << ',' << orientation << ',';

				switch (optionGenomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[gseqNo];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						if (multipleGenomes) {
							file.fill('0');
							file << genomeName << '#' << setw(gzeros) << gseqNo + 1;
						} else
							file << genomeName;
				}

				file << ',' << gBegin << ',' << gEnd << ',' << setprecision(5) << percId << endl;

				if (optionDumpAlignment) {
					assignSource(row(align, 0), reads[readNo]);
					assignSource(row(align, 1), infix(genomes[gseqNo], gBegin, gEnd));
#ifndef OMIT_REVERSECOMPLEMENT
					if (orientation == 'R')
						reverseComplementInPlace(source(row(align, 1)));
#endif
					globalAlignment(align, scoreType);
					dumpAlignment(file, align);
				}
			}
			break;


		case 1:	// Enhanced Fasta Format
			for(; it != itEnd; ++it) 
			{
				if (readNo == (*it).rseqNo && _gBegin == (*it).gBegin && gseqNo == (*it).gseqNo)
					continue;

				readNo = (*it).rseqNo;
				gseqNo = (*it).gseqNo;
				_gBegin = (*it).gBegin;
				gEnd = (*it).gEnd;

				unsigned	readLen = length(reads[readNo]);
				double		percId;

				gBegin = _gBegin;

				if (optionPositionFormat == 1)
					++gBegin;

				string fastaID;
				assign(fastaID, readIDs[readNo]);

				int id = readNo;
				int fragId = readNo;

				size_t left = fastaID.find_first_of('[');
				size_t right = fastaID.find_last_of(']');
				if (left != fastaID.npos && right != fastaID.npos && left < right) 
				{
					fastaID.erase(right);
					fastaID.erase(0, left + 1);
					replace(fastaID.begin(), fastaID.end(), ',', ' ');
					size_t pos = fastaID.find("id=");
					if (pos != fastaID.npos) {
						istringstream iss(fastaID.substr(pos + 3));
						iss >> id;
					}
					pos = fastaID.find("fragId=");
					if (pos != fastaID.npos) {
						istringstream iss(fastaID.substr(pos + 7));
						iss >> fragId;
					}
				}

				percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);

				if ((*it).orientation == 'F')		
					file << '>' << gBegin << ',' << gEnd;	// forward strand
				else
					file << '>' << gEnd << ',' << gBegin;	// reverse strand (switch begin and end)
				file << "[id=" << id << ",fragId=" << fragId;
				file << ",errors=" << (*it).editDist << ",percId=" << setprecision(5) << percId << ']' << endl;

				file << reads[readNo] << endl;
			}
			break;
	}

	file.close();

	if (_debugLevel >= 1)
		cerr << "Dumping results took             \t" << SEQAN_PROTIMEDIFF(dump_time) << " seconds" << endl;
}

//////////////////////////////////////////////////////////////////////////////
// Main Part

template <typename TShape>
int mapReads(const char *genomeFileName, const char *readFileName, TShape &shape)
{
	typedef Index<TReadSet, Index_QGram<TShape> >	TReadIndex;

	TGenomeSet				genomeSet;
	StringSet<CharString>	genomeNames;	// genome names, taken from the Fasta file
	TReadSet				readSet;
	StringSet<CharString>	readNames;		// read names, taken from the Fasta file
	
	TReadIndex				swiftIndex(readSet, shape);
	TMatches				matches;		// resulting forward/reverse matches

	static const double epsilon = 0.0000001;

	// dump configuration in verbose mode
	if (_debugLevel >= 1) 
	{
		CharString bitmap;
		shapeToString(bitmap, indexShape(swiftIndex));
		cerr << "___SETTINGS____________" << endl;
		cerr << "Compute forward matches:         \t";
		if (optionForward)	cerr << "YES" << endl;
		else				cerr << "NO" << endl;
		cerr << "Compute reverse matches:         \t";
		if (optionRev)		cerr << "YES" << endl;
		else				cerr << "NO" << endl;
		cerr << "Error rate:                      \t" << optionErrorRate << endl;
		cerr << "Minimal threshold:               \t" << optionThreshold << endl;
		cerr << "Shape:                           \t" << bitmap << endl;
		cerr << "Repeat threshold:                \t" << optionRepeatLength << endl;
		cerr << "Overabundance threshold:         \t" << optionAbundanceCut << endl;
		cerr << "Taboo length:                    \t" << optionTabooLength << endl;
		cerr << endl;
	}

	// circumvent numerical obstacles
	optionErrorRate += epsilon;


	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Load fasta files
	if (!loadFasta(genomeSet, genomeNames, genomeFileName)) {
		cerr << "Failed to load genomes" << endl;
		return 1;
	}
	if (_debugLevel >= 1) cerr << lengthSum(genomeSet) << " bps of " << length(genomeSet) << " genomes loaded." << endl;

	if (!loadFasta(readSet, readNames, readFileName)) {
		cerr << "Failed to load reads" << endl;
		return 1;
	}
	if (_debugLevel >= 1) cerr << lengthSum(readSet) << " bps of " << length(readSet) << " reads loaded." << endl;


	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Find matches using SWIFT
	if (optionForward)
	{
		if (_debugLevel >= 1)
			cerr << endl << "___FORWARD_STRAND______";
		findReads(matches, genomeSet, swiftIndex, 'F');
	}

#ifndef OMIT_REVERSECOMPLEMENT
	if (optionRev) 
	{
		if (_debugLevel >= 1)
			cerr << endl << "___BACKWARD_STRAND_____";
		reverseComplementInPlace(genomeSet);			// build reverse-compl of genome
		findReads(matches, genomeSet, swiftIndex, 'R');
		reverseComplementInPlace(genomeSet);			// restore original genome seqs
	}
#endif

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Remove duplicates and output matches
	dumpMatches(matches, genomeSet, genomeNames, readSet, readNames, genomeFileName, readFileName);
	return 0;
}	


//////////////////////////////////////////////////////////////////////////////
// Print usage
void printHelp(int, const char *[], bool longHelp = false) 
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
		cerr << "  -o,  --output FILE           \t" << "change output filename (default: <READS FILE>.result)" << endl;
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
		cerr << "  -s,  --shape                 \t" << "set k-mer shape (binary string, default " << optionShape << ')' << endl;
		cerr << "  -t,  --threshold NUM         \t" << "set minimum k-mer threshold (default " << optionThreshold << ")" << endl;
		cerr << "  -oc, --overabundance-cut NUM \t" << "set k-mer overabundance cut ratio (default " << optionAbundanceCut << ")" << endl;
		cerr << "  -rl, --repeat-length NUM     \t" << "set simple-repeat length threshold (default " << optionRepeatLength << ")" << endl;
		cerr << "  -tl, --taboo-length NUM      \t" << "set taboo length (default " << optionTabooLength << ")" << endl;
	} else {
		cerr << "Try 'razers --help' for more information." << endl;
	}
}

void printVersion() 
{
	cerr << "RazerS version 0.3 20080821 (prerelease)" << endl;
}

int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line

	unsigned			fnameCount = 0;
	const char			*fname[2] = { "", "" };

	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option

			if (strcmp(argv[arg], "-f") == 0 || strcmp(argv[arg], "--forward") == 0) {
				optionForward = true;
				continue;
			}
			if (strcmp(argv[arg], "-r") == 0 || strcmp(argv[arg], "--reverse") == 0) {
				optionRev = true;
				continue;
			}
			if (strcmp(argv[arg], "-i") == 0 || strcmp(argv[arg], "--percent-identity") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionErrorRate;
					optionErrorRate = (100.0 - optionErrorRate) / 100.0;
					if (!istr.fail()) 
					{
						if (optionErrorRate < 0 || optionErrorRate > 0.5)
							cerr << "Percent identity threshold must be a value between 50 and 100" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-a") == 0 || strcmp(argv[arg], "--alignment") == 0) {
				optionDumpAlignment = true;
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
			if (strcmp(argv[arg], "-of") == 0 || strcmp(argv[arg], "--output-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionOutputFormat;
					if (!istr.fail())
					{
						if (optionOutputFormat > 1)
							cerr << "Invalid output format option" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-so") == 0 || strcmp(argv[arg], "--sort-order") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionSortOrder;
					if (!istr.fail())
					{
						if (optionSortOrder > 1)
							cerr << "Invalid sort order option" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-gn") == 0 || strcmp(argv[arg], "--genome-naming") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionGenomeNaming;
					if (!istr.fail())
					{
						if (optionGenomeNaming > 1)
							cerr << "Invalid genome naming option" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-rn") == 0 || strcmp(argv[arg], "--read-naming") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionReadNaming;
					if (!istr.fail())
					{
						if (optionReadNaming > 2)
							cerr << "Invalid read naming option" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-pf") == 0 || strcmp(argv[arg], "--position-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionPositionFormat;
					if (!istr.fail())
					{
						if (optionPositionFormat > 1)
							cerr << "Invalid position format option" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-s") == 0 || strcmp(argv[arg], "--shape") == 0){
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionShape;
					if (istr.fail())
						cerr << "Could not parse shape" << endl << endl;
					else
					{
						unsigned ones = 0;
						unsigned zeros = 0;
						for(unsigned i = 0; i < length(optionShape); ++i)
							switch (optionShape[i])
							{
								case '0':
									++zeros;
									break;
								case '1':
									++ones;
									break;
								default:
									cerr <<	"Shape must be a binary string" << endl << endl;
									printHelp(argc, argv);
									return 0;
							}
						if (ones == 0 || ones > 20) 
						{
							cerr <<	"Invalid Shape" << endl << endl;
							printHelp(argc, argv);
							return 0;
						}

						if (ones < 7 || ones > 14)
							cerr <<	"Warning: Shape should contain at least 7 and at most 14 '1's" << endl << endl;

						continue;				
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-oc") == 0 || strcmp(argv[arg], "--overabundance-cut") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionAbundanceCut;
					if (!istr.fail())
					{
						if (optionAbundanceCut <= 0 || optionAbundanceCut > 1)
							cerr << "Overabundance cut ratio must be a value >0 and <=1. Set to 1 to disable." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-rl") == 0 || strcmp(argv[arg], "--repeat-length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionRepeatLength;
					if (!istr.fail())
					{
						if (optionRepeatLength <= 1)
							cerr << "Repeat length must be a value greater than 1" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-t") == 0 || strcmp(argv[arg], "--threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionThreshold;
					if (!istr.fail())
					{
						if (optionThreshold < 1)
							cerr << "Threshold must be a value greater than 0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-tl") == 0 || strcmp(argv[arg], "--taboo-length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionTabooLength;
					if (!istr.fail())
					{
						if (optionTabooLength < 1)
							cerr << "Taboo length must be a value greater than 0" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, true);
				return 0;
			}
			if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
				_debugLevel = max(_debugLevel, 1);
				continue;
			}
			if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--vverbose") == 0) {
				_debugLevel = 3;
				continue;
			}
			if (strcmp(argv[arg], "-V") == 0 || strcmp(argv[arg], "--version") == 0) {
				optionPrintVersion = true;
				continue;
			}
		} else {
			// parse file name
			if (fnameCount == 2) {
				printHelp(argc, argv);
				return 0;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 2) {
		if (optionPrintVersion)
			printVersion();
		else
			printHelp(argc, argv);
		return 0;
	}
	if (!optionForward && !optionRev) { // enable both per default
		optionForward = true;
		optionRev = true;
	}

	if (optionPrintVersion)
		printVersion();

	Shape<TAlphabet, SimpleShape>		ungapped;
	Shape<TAlphabet, OneGappedShape>	onegapped;
	Shape<TAlphabet, GappedShape>		gapped;

	if (stringToShape(ungapped, optionShape))
		return mapReads(fname[0], fname[1], ungapped);
	
	if (stringToShape(onegapped, optionShape))
		return mapReads(fname[0], fname[1], onegapped);

	if (stringToShape(gapped, optionShape))
		return mapReads(fname[0], fname[1], gapped);

	cerr <<	"Invalid Shape" << endl << endl;
	printHelp(argc, argv);
	return 0;
}
