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
//#define RAZERS_HAMMINGVERIFY	// 
//#define WITH_1HULL			// allow 1-error per qgram

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
	static int			optionThreshold = 1;				// threshold
	static int			optionTabooLength = 1;				// taboo length
	static int			_debugLevel = 0;					// level of verbosity
	static bool			optionPrintVersion = false;			// print version number


//////////////////////////////////////////////////////////////////////////////
// Typedefs

	// definition of a Read match
	template <typename _TGPos>
	struct ReadMatch {
		typedef _TGPos TGPos;
		unsigned gseqNo;	// genome seqNo
		unsigned rseqNo;	// read seqNo
		TGPos	gBegin;		// begin position of the match in the genome
		TGPos	gEnd;		// end position of the match in the genome
		short	editDist;	// Levenshtein distance
	};

	// (to sort and remove duplicates)
	template <typename TReadMatch>
	struct ReadMatchLess : public binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const {
			return (a.gBegin < b.gBegin) || (a.gBegin == b.gBegin && a.gEnd < b.gEnd);
		}
	};

	template <typename TReadMatch>
	struct ReadMatchGreater : public binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const {
			return (a.gBegin > b.gBegin) || (a.gBegin == b.gBegin && a.gEnd > b.gEnd);
		}
	};

	typedef Pair<short,	short, Compressed> PrimerSAValue;


//////////////////////////////////////////////////////////////////////////////
// Definitions

	typedef String<Dna5>		TGenome;
	typedef StringSet<TGenome>	TGenomeSet;
	typedef String<TAlphabet>	TRead;
	typedef StringSet<TRead>	TReadSet;

	typedef Index<TReadSet, Index_QGram<FixedShape<QGRAM_LEN> > >	TReadIndex;

namespace seqan 
{
	template <>
	struct SAValue<TReadIndex> {
		typedef Pair<
			unsigned,		// many
			unsigned,		// short reads
			Compressed
		> Type;
	};
}

	typedef SAValue<TGenomeSet>::Type				TSAValue;
	typedef Repeat<TSAValue, unsigned>				TRepeat;

	typedef ReadMatch<Difference<TGenome>::Type>	TMatch;			// a single match
	typedef String<TMatch/*, Block<>*/ >			TMatches[2];	// array[orientation][seqNo] of matches

#ifdef RAZERS_PRUNE_QGRAM_INDEX
namespace seqan 
{
	//////////////////////////////////////////////////////////////////////////////
	// Repeat masker
	template <>
	inline bool _qgramDisableBuckets(TReadIndex &index) 
	{
		typedef Fibre<TReadIndex, QGram_Dir>::Type	TDir;
		typedef Iterator<TDir, Standard>::Type		TDirIterator;
		typedef Value<TDir>::Type					TSize;

		TDir &dir    = indexDir(index);
		bool result  = false;
		TSize thresh = length(index) / 5000;
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
			}

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
template <typename TMatches, typename TGenomeSet, typename TReadIndex, typename TRepeats>
void findReads(
	TMatches &matches,			// resulting matches
	TGenomeSet &genomes,		// Genome
	TReadIndex &readIndex,		// q-gram index of Reads
	TRepeats &repeats)
{
	typedef typename Fibre<TReadIndex, Fibre_Text>::Type	TReadSet;
	typedef typename Value<TGenomeSet>::Type				TGenome;
	typedef typename Size<TGenome>::Type					TSize;
	typedef typename Infix<TGenome>::Type					TGenomeInfix;
	typedef typename Value<TReadSet>::Type					TRead;
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Value<TRepeats>::Type					TRepeat;
	typedef typename Iterator<TRepeats, Standard>::Type		TRepeatIterator;

	// *** SPECIALIZATION ***

	typedef Pipe< TGenomeInfix, Source<> >							TSource;
	typedef Pipe< TSource, Caster<Dna> >							TCaster;
	typedef Pipe< TCaster, Tupler<QGRAM_LEN, true, Compressed> >	TTupler;

#ifdef WITH_1HULL
	// with 1-hull enumeration
	typedef Pipe< TTupler, EditEnvironment<LevenshteinDistance> >	TEnumerator;

	typedef Finder<TEnumerator, Swift<> >				TSwiftFinder;
#else
	// without 1-hull enumeration
	typedef Finder<TTupler, Swift<> >					TSwiftFinder;
#endif

	typedef Pattern<TReadIndex, Swift<> >				TSwiftPattern;


	// find read match end
	typedef Finder<TGenomeInfix>						TMyersFinder;
#ifdef RAZERS_HAMMINGVERIFY
	typedef Score<int>									TScore;
	typedef Pattern<TRead, DPSearch<TScore> >			TMyersPattern;
#else
	typedef Pattern<TRead, MyersUkkonen>				TMyersPattern;
#endif

	// find read match begin
	typedef ModifiedString<TGenomeInfix, ModReverse>	TGenomeInfixRev;
	typedef ModifiedString<TRead, ModReverse>			TReadRev;
	typedef Finder<TGenomeInfixRev>						TMyersFinderRev;
	typedef Pattern<TReadRev, MyersUkkonenGlobal>		TMyersPatternRev;

	TSwiftPattern swiftPattern(readIndex);
	swiftPattern.params.minThreshold = optionThreshold;
	swiftPattern.params.tabooLength = optionTabooLength;

	// VERIFICATION
/*	String<TMyersPattern> myersPattern;
	resize(myersPattern, countSequences(readIndex));
	for(unsigned i = 0; i < countSequences(readIndex); ++i)
		setHost(myersPattern[i], indexText(readIndex)[i]);
*//*
	String<char> text2 = "xxgenerationsxx";
	Finder<String<char> > find2(text2);

	while (find(find2, myersPattern[1], -(int)3)) {
		cerr << getScore(myersPattern[1]) << ": ";
		cerr << prefix(text2, position(find2));
		cerr << endl;
	}
*/

#ifdef RAZERS_HAMMINGVERIFY
	TScore	scoreType = Score<int>(0, -1, -1001, -1000); // levenshtein-score (match, mismatch, gapOpen, gapExtend)
#endif

	__int64 TP = 0;
	__int64 FP = 0;
	SEQAN_PROTIMESTART(find_time);

	TRepeatIterator repeatIt = begin(repeats, Standard());
	TRepeatIterator repeatItEnd = end(repeats, Standard());
	String<Pair<TSize> > borders;

	// iterate all genomic sequences
	for(unsigned hstkSeqNo = 0; hstkSeqNo < length(genomes); ++hstkSeqNo) 
	{
		TGenome &genome = genomes[hstkSeqNo];

		// get intervals of non-repeat regions
		resize(borders, 0);
		TSize last = 0;
		while (repeatIt != repeatItEnd && (*repeatIt).beginPosition.i1 == hstkSeqNo)
		{
			appendValue(borders, Pair<TSize>(last, (*repeatIt).beginPosition.i2));
			last = (*repeatIt).endPosition.i2;
			++repeatIt;
		}
		appendValue(borders, Pair<TSize>(last, length(genome)));

		// iterate all non-repeat regions
		for(unsigned segment = 0; segment < length(borders); ++segment)
		{
			if (_debugLevel >= 1)
				::std::cerr << ::std::endl << "Process genome seq:" << hstkSeqNo << " range:" << borders[segment] << ::std::endl;
			// 1-HULL
			TGenomeInfix	genomeInf(genome, borders[segment].i1, borders[segment].i2);
			TSource			src(genomeInf);
			TCaster			caster(src);
			TTupler			tupler(caster);
#ifdef WITH_1HULL
			TEnumerator		enumerator(tupler);
			// SWIFT
			TSwiftFinder	swiftFinder(enumerator);
#else
			TSwiftFinder	swiftFinder(tupler);
#endif
			// iterate all verification regions returned bz SWIFT
			while (find(swiftFinder, swiftPattern, optionErrorRate, (_debugLevel >= 1))) 
			{
				unsigned ndlSeqNo = (*swiftFinder.curHit).ndlSeqNo;
				unsigned ndlLength = sequenceLength(ndlSeqNo, readIndex);
		//		unsigned hstkPos = (*swiftFinder.curHit).hstkPos;
		//		unsigned bucketWidth = (*swiftFinder.curHit).bucketWidth;

		//		TGenomeInfix inf(genomeInf, hstkPos, hstkPos + bucketWidth);
				TGenomeInfix inf(range(swiftFinder, genomeInf));
				TMyersFinder myersFinder(inf);
#ifdef RAZERS_HAMMINGVERIFY
				TMyersPattern myersPattern(indexText(readIndex)[ndlSeqNo], scoreType);
#else
				TMyersPattern myersPattern(indexText(readIndex)[ndlSeqNo]);
#endif


#ifdef RAZERS_DEBUG
				cout<<"Verify: "<<endl;
				cout<<"Genome: "<<inf<<endl;
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
						beginPosition(genomeInf) + beginPosition(swiftFinder),
						beginPosition(genomeInf) + beginPosition(swiftFinder) + position(maxPos) + 1,
						-maxScore
					};
					if (m.gBegin < 0) m.gBegin = 0;


					TGenomeInfixRev		infRev(infix(genome, m.gBegin, m.gEnd));
					TReadRev			readRev(indexText(readIndex)[ndlSeqNo]);
					TMyersFinderRev		myersFinderRev(infRev);
					TMyersPatternRev	myersPatternRev(readRev);

					// find beginning of best semi-global alignment
					if (find(myersFinderRev, myersPatternRev, maxScore))
						m.gBegin = m.gEnd - (position(myersFinderRev) + 1);
					else {
						// this case should never occur
						::std::cerr << "1GENOME: " << host(myersFinder) << ::std::endl;
						::std::cerr << "1READ:   " << indexText(readIndex)[ndlSeqNo] << ::std::endl;
						::std::cerr << "2GENOME: " << infix(genome, m.gBegin, m.gEnd) << ::std::endl;
						::std::cerr << "2READ:   " << indexText(readIndex)[ndlSeqNo] << ::std::endl;
						::std::cerr << "3GENOME: " << infRev << ::std::endl;
						::std::cerr << "3READ:   " << readRev << ::std::endl;
						::std::cerr << "HUH?" << ::std::endl;
					}

#ifndef RAZERS_PROFILE
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
	}
	if (_debugLevel >= 1)
		cerr << endl << "Finding Reads took               \t" << SEQAN_PROTIMEDIFF(find_time) << " seconds" << endl;

	if (_debugLevel >= 2) {
		cerr << endl;
		cerr << "___Filtration_Stats____" << endl;
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
	Score<int> scoreType = Score<int>(0, -1000, -1001, -1000); // levenshtein-score (match, mismatch, gapOpen, gapExtend)
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

	sort(
		begin(matches[0], Standard()),
		end(matches[0], Standard()), 
		ReadMatchLess<TMatch>());

	sort(
		begin(matches[1], Standard()),
		end(matches[1], Standard()), 
		ReadMatchGreater<TMatch>());

	for(unsigned orientation = 0; orientation < 2; ++orientation)
	{
		typename Iterator<TMatches, Standard>::Type it = begin(matches[orientation], Standard());
		typename Iterator<TMatches, Standard>::Type itEnd = end(matches[orientation], Standard());

		typename TMatch::TGPos lastBegin = -1;
		typename TMatch::TGPos gBegin, gEnd;


		switch (optionOutputFormat) 
		{
			case 0:	// Razer Format
				for(; it != itEnd; ++it) 
				{
					if (lastBegin != (*it).gBegin) 
					{
						unsigned	gseqNo = (*it).gseqNo;
						unsigned	readNo = (*it).rseqNo;
						TGPos		gLength = length(genomes[gseqNo]);
						unsigned	readLen = length(reads[readNo]);
						double		percId;

						lastBegin = (*it).gBegin;

						if (orientation == 0) {
							gBegin = lastBegin;
							gEnd = (*it).gEnd;
						} else {
							gBegin = gLength - (*it).gEnd;
							gEnd = gLength - lastBegin;
						}

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
								file << readName << "#" << setw(pzeros) << readNo + 1;
								break;

							// 2..filename is the read sequence itself
							case 2:
								file << reads[readNo];
						}

						if (orientation == 0)
							file << ",0," << readLen << ",F,";
						else
							file << ",0," << readLen << ",R,";

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
									file << genomeName << "#" << setw(gzeros) << gseqNo + 1;
								} else
									file << genomeName;
						}

						file << "," << gBegin << "," << gEnd << "," << setprecision(5) << percId << endl;

						if (optionDumpAlignment) {
							assignSource(row(align, 0), reads[readNo]);
							assignSource(row(align, 1), infix(genomes[gseqNo], gBegin, gEnd));
#ifndef OMIT_REVERSECOMPLEMENT
							if (orientation == 1)
								reverseComplementInPlace(source(row(align, 1)));
#endif
							globalAlignment(align, scoreType);
							dumpAlignment(file, align);
						}
					}
				}
				break;


			case 1:	// Enhanced Fasta Format
				for(; it != itEnd; ++it) 
				{
					if (lastBegin != (*it).gBegin) 
					{
						unsigned	gseqNo = (*it).gseqNo;
						unsigned	readNo = (*it).rseqNo;
						TGPos		gLength = length(genomes[gseqNo]);
						unsigned	readLen = length(reads[readNo]);
						double		percId;

						lastBegin = (*it).gBegin;

						if (orientation == 0) {
							gBegin = lastBegin + 1;
							gEnd = (*it).gEnd;
						} else {
							gBegin = gLength - lastBegin;
							gEnd = gLength - (*it).gEnd + 1;
						}

						percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);

						file << ">" << gBegin << "," << gEnd;
						file << "[id=" << readNo << ",fragId=" << readNo;
						file << ",errors=" << (*it).editDist << ",percId=" << setprecision(5) << percId << "]" << endl;

						file << reads[readNo] << endl;
					}
				}
				break;
		}
	}

	file.close();

	if (_debugLevel >= 1)
		cerr << "Dumping results took             \t" << SEQAN_PROTIMEDIFF(dump_time) << " seconds" << endl;
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
		cerr << endl << "Options:" << endl;
		cerr << "  -f,  --forward               \t" << "only compute forward matches" << endl;
		cerr << "  -r,  --reverse               \t" << "only compute reverse complement matches" << endl;
		cerr << "  -i,  --percent-identity NUM  \t" << "set the percent identity threshold" << endl;
		cerr << "                               \t" << "default value is 80" << endl;
		cerr << "  -a,  --alignment             \t" << "also dump the alignment for each match" << endl;
		cerr << "  -o,  --output FILE           \t" << "change output filename (default: <READS FILE>.result)" << endl;
		cerr << "  -of, --output-format NUM     \t" << "change output format" << endl;
		cerr << "                               \t" << "0 = Razer format (default, see README)" << endl;
		cerr << "                               \t" << "1 = enhanced Fasta format" << endl;
		cerr << "  -GN, --genome-naming NUM     \t" << "select how genomes are named" << endl;
		cerr << "                               \t" << "0 = use Fasta id (default)" << endl;
		cerr << "                               \t" << "1 = enumerate beginning with 1" << endl;
		cerr << "  -RN, --read-naming NUM       \t" << "select how reads are named" << endl;
		cerr << "                               \t" << "0 = use Fasta id (default)" << endl;
		cerr << "                               \t" << "1 = enumerate beginning with 1" << endl;
		cerr << "                               \t" << "2 = use the read sequence (only for short reads!)" << endl;
		cerr << "  -t,  --threshold NUM         \t" << "set minimum threshold (default " << optionThreshold << ")" << endl;
		cerr << "  -tl, --taboo-length NUM      \t" << "set taboo length (default " << optionTabooLength << ")" << endl;
		cerr << "  -v,  --verbose               \t" << "verbose mode" << endl;
		cerr << "  -vv, --vverbose              \t" << "very verbose mode" << endl;
		cerr << "  -V,  --version               \t" << "print version number" << endl;
		cerr << "  -h,  --help                  \t" << "print this help" << endl;
	} else {
		cerr << "Try 'razers --help' for more information." << endl;
	}
}

void printVersion() 
{
	cerr << "RazerS version 1.0 20080805 (prerelease)" << endl;
}

int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line

	unsigned			fnameCount = 0;
	const char			*fname[2] = { "", "" };
	static const double epsilon = 0.0000001;

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
						if (optionErrorRate < 0 || optionErrorRate > 0.5)
							cerr << "Percent identity threshold must be a value between 50 and 100" << endl << endl;
						else
							continue;
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
						if (optionOutputFormat > 1)
							cerr << "Invalid output format option" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-GN") == 0 || strcmp(argv[arg], "--genome-naming") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionGenomeNaming;
					if (!istr.fail())
						if (optionGenomeNaming > 1)
							cerr << "Invalid genome naming option" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-RN") == 0 || strcmp(argv[arg], "--read-naming") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionReadNaming;
					if (!istr.fail())
						if (optionReadNaming > 2)
							cerr << "Invalid read naming option" << endl << endl;
						else
							continue;
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
						if (optionThreshold < 1)
							cerr << "Threshold must be a value greater than 0" << endl << endl;
						else
							continue;
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
						if (optionTabooLength < 1)
							cerr << "TabooLength must be a value greater than 0" << endl << endl;
						else
							continue;
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

	// dump configuration in verbose mode
	if (_debugLevel >= 1) {
		cerr << "Compute forward matches:         \t";
		if (optionForward)	cerr << "YES" << endl;
		else				cerr << "NO" << endl;
		cerr << "Compute reverse matches:         \t";
		if (optionRev)		cerr << "YES" << endl;
		else				cerr << "NO" << endl;
		cerr << "Error rate:                      \t" << optionErrorRate << endl;
		cerr << "Minimal threshold:               \t" << optionThreshold << endl;
		cerr << "Taboo length:                    \t" << optionTabooLength << endl;
		cerr << endl;
	}

	// circumvent numerical obstacles
	optionErrorRate += epsilon;


	TGenomeSet				genomeSet;
	StringSet<CharString>	genomeNames;	// genome names, taken from the Fasta file
	TReadSet				readSet;
	StringSet<CharString>	readNames;		// read names, taken from the Fasta file
	
	TReadIndex				swiftIndex(readSet);
	TMatches				matches;		// resulting forward/reverse matches


	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Load fasta files
	if (!loadFasta(genomeSet, genomeNames, fname[0])) {
		cerr << "Failed to load Genome" << endl;
		return 1;
	}
	if (_debugLevel >= 1) cerr << lengthSum(genomeSet) << " bps in " << length(genomeSet) << " Genomes loaded." << endl;

	String<TRepeat> repeats;
	findRepeats(repeats, genomeSet, 1000, 1);

	if (!loadFasta(readSet, readNames, fname[1])) {
		cerr << "Failed to load Primers" << endl;
		return 1;
	}
	if (_debugLevel >= 1) cerr << lengthSum(readSet) << " bps in " << length(readSet) << " Reads loaded." << endl;


	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Find matches using SWIFT
	if (optionForward)
		findReads(matches[0], genomeSet, swiftIndex, repeats);

#ifndef OMIT_REVERSECOMPLEMENT
	if (optionRev) {
		reverseComplementInPlace(genomeSet);			// build reverse-compl of genome
		findRepeats(repeats, genomeSet, 1000, 1);
		findReads(matches[1], genomeSet, swiftIndex, repeats);
		reverseComplementInPlace(genomeSet);			// restore original genome seqs
	}
#endif

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Remove duplicates and output matches
	dumpMatches(matches, genomeSet, genomeNames, readSet, readNames, fname[0], fname[1]);
	
	return 0;
}
