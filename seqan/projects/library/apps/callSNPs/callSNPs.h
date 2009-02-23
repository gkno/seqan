 /*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#ifndef SEQAN_HEADER_CALLSNPS_H
#define SEQAN_HEADER_CALLSNPS_H

#include <iostream>
#include <fstream>

//#include "../razers/razers_utils.h"


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Default options

	template < bool _HAMMING_ONLY = true >
	struct SNPCallingSpec 
	{
		enum { HAMMING_ONLY = _HAMMING_ONLY };				// omit verifying potential matches
	};

	template < typename TSpec = SNPCallingSpec<> >
	struct SNPCallingOptions
	{
	// main options
		TSpec		spec;
		bool		forward;			// compute forward oriented read matches
		bool		reverse;			// compute reverse oriented read matches
		double		errorRate;			// Criteria 1 threshold
		unsigned	maxHits;			// output at most maxHits many matches
		unsigned	distanceRange;		// output only the best, second best, ..., distanceRange best matches
		bool		purgeAmbiguous;		// true..remove reads with more than maxHits best matches, false..keep them
		const char	*output;			// name of result file
		int			_debugLevel;		// level of verbosity
		bool		printVersion;		// print version number
		bool		hammingOnly;		// no indels

	// output format options
		unsigned	outputFormat;		// 0..Razer format
		unsigned	inputFormat;		// 0 = razers, 1=eland, 2 = maq
										// 1..enhanced Fasta
		bool		dumpAlignment;		// compute and dump the match alignments in the result files
		unsigned	genomeNaming;		// 0..use Fasta id
										// 1..enumerate reads beginning with 1
		unsigned	readNaming;			// 0..use Fasta id
										// 1..enumerate reads beginning with 1
										// 2..use the read sequence (only for short reads!)
		unsigned	sortOrder;			// 0..sort keys: 1. read number, 2. genome position
										// 1..           1. genome pos50ition, 2. read number
		unsigned	positionFormat;		// 0..gap space
										// 1..position space

	// filtration parameters
		::std::string shape;			// shape (e.g. 11111111111)
		int			threshold;			// threshold
		int			tabooLength;		// taboo length
		int			repeatLength;		// repeat length threshold
		double		abundanceCut;		// abundance threshold

	// verification parameters
		bool		matchN;				// false..N is always a mismatch, true..N matches with all
		unsigned char compMask[5];

	// statistics
		__int64		FP;					// false positives (threshold reached, no match)
		__int64		TP;					// true positives (threshold reached, match)
		double		timeLoadFiles;		// time for loading input files
		double		timeMapReads;		// time for mapping reads
		double		timeDumpResults;	// time for dumping the results

		unsigned	readLength;		// specify read length, for direct trimming
		const char	*runID;			// runID needed for gff output
		bool		fastQ;			// reads are given in fastq format (-->use quality values)

		const char 	*positionFile;
		const char 	*outputPositionAnalysis;

		unsigned method;
		unsigned	maxPile;			// 
		double testLevel;
		const char	*outputSNP;			// name of result file
		const char	*errorPrbFileName;			// name of result file
	
		bool		qualityFile;			// name of .qual input file
		const char	*tabFile;			// name of .tab.txt result file
		bool 		useBaseQuality;		// use base quality instead of min{base,mapping}
		float		avgQualT;
		float		percentageT;
		unsigned	minMutT;
#ifdef RAZERS_MAQ_MAPPING
		bool maqMapping;
#endif
	
	// misc
		unsigned	compactThresh;		// compact match array if larger than compactThresh

		unsigned	minCoverage;
		unsigned	extractMinCov;
		unsigned	extractMaxCov;
		bool 		schnellUndDreckig;
		bool 		orientationAware;

		unsigned 	rLen; //temporary		
		SNPCallingOptions() 
		{
			forward = true;
			reverse = true;
			errorRate = 0.08;
			maxHits = 100;
			distanceRange = 0;	// disabled
			purgeAmbiguous = false;
			output = "";
			qualityFile = false;
			tabFile = "";
			
			_debugLevel = 0;
			printVersion = false;
			hammingOnly = true;
			
			outputFormat = 0;
			dumpAlignment = false;
			genomeNaming = 0;
			readNaming = 0;
			sortOrder = 0;
			positionFormat = 1;
			
			matchN = false;
			
			shape = "11111111111";
			threshold = 1;
			tabooLength = 1;
			repeatLength = 1000;
			abundanceCut = 1;
			
			for (unsigned i = 0; i < 4; ++i)
				compMask[i] = 1 << i;
			compMask[4] = 0;
			
			compactThresh = 1024;
			
			readLength = 0; //-> use entire read length
			runID = "run"; //
			fastQ = true;
			
			rLen = 32;
			method = 0;
			testLevel = 0.01;
			outputSNP = "";
			positionFile = "";
			outputPositionAnalysis = "";
			errorPrbFileName = "";
			inputFormat = 0;
			
			maxPile = 4;
			avgQualT = 10;
			percentageT = 0.3;
			minMutT = 5;
			useBaseQuality = false;
			
			minCoverage = 5;
			extractMinCov = 0;
			extractMaxCov = 1024*4;	//sollte reichen
			schnellUndDreckig = true;
			orientationAware = false;
		
			
#ifdef RAZERS_MAQ_MAPPING
			maqMapping = false;
#endif

		}
	};

//////////////////////////////////////////////////////////////////////////////
// Typedefs

	// definition of a Read match
	template <typename _TGPos>
	struct MappedReadMatch 
	{
		typedef typename _MakeSigned<_TGPos>::Type TGPos;
		
		unsigned		gseqNo;			// genome seqNo
		unsigned		rseqNo;			// read seqNo
		TGPos			gBegin;			// begin position of the match in the genome
		TGPos			gEnd;			// end position of the match in the genome
		unsigned short		editDist;		// Levenshtein distance
		char			orientation;		// 'F'..forward strand, 'R'..reverse comp. strand

		short			mScore;			// mapping quality
		short			avgQuality;		// avg read quality
	};
	
	enum CALLSNPS_ERROR {
		CALLSNPS_GFF_FAILED = 1,
		CALLSNPS_GENOME_FAILED = 2,
		CALLSNPS_POS_FAILED = 3,
		CALLSNPS_QUALITY_FAILED = 4
	};

//////////////////////////////////////////////////////////////////////////////
// Definitions
	// ... to sort matches and remove duplicates with equal gBegin
	template <typename TReadMatch>
	struct LessGPosMQ : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;

			if (a.mScore > b.mScore) return true;
			if (a.mScore < b.mScore) return false;

			if (a.avgQuality > b.avgQuality) return true;
			if (a.avgQuality < b.avgQuality) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};

	template <typename TReadMatch>
	struct LessGPosOaMQ : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
 			if (a.orientation < b.orientation) return true;
 			if (a.orientation > b.orientation) return false;

			if (a.mScore > b.mScore) return true;
			if (a.mScore < b.mScore) return false;

			if (a.avgQuality > b.avgQuality) return true;
			if (a.avgQuality < b.avgQuality) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};

template <typename TGenomeSet>
bool loadGenomes(TGenomeSet &genomes, StringSet<CharString> &fileNameList, ::std::map<CharString,unsigned> &gIdStringToIdNumMap)
{
	unsigned gSeqNo = 0;
	unsigned filecount = 0;
	CharString temp;
	while(filecount < length(fileNameList))
	{
		clear(temp);
		MultiFasta multiFasta;
		if (!open(multiFasta.concat, toCString(fileNameList[filecount]), OPEN_RDONLY)) return false;
		split(multiFasta, Fasta());
		
		unsigned seqCount = length(multiFasta);
		if(length(genomes) < gSeqNo+seqCount) 
			resize(genomes,gSeqNo+seqCount);
		for(unsigned i = 0; i < seqCount; ++i)
		{
			assignSeq(genomes[gSeqNo+i], multiFasta[i], Fasta());		// read Genome sequence
			assignSeqId(temp, multiFasta[i], Fasta());
			for (unsigned pos = 0; pos < length(temp); ++pos)
			{
				if(temp[pos]=='\t' || temp[pos]=='\b' || temp[pos]==' ')
				{
					resize(temp,pos);
					break;
				}
			}
			gIdStringToIdNumMap.insert(::std::make_pair<CharString,unsigned>(temp,gSeqNo+i)); //TODO shortID
		}
		gSeqNo += seqCount;
		++filecount;
	}
	resize(genomes,gSeqNo);
	return (gSeqNo > 0);
}



double cumulated_normal(const double x)
{
	const double b1 =  0.319381530;
	const double b2 = -0.356563782;
	const double b3 =  1.781477937;
	const double b4 = -1.821255978;
	const double b5 =  1.330274429;
	const double p  =  0.2316419;
	const double c  =  0.39894228;

	if(x >= 0.0) 
	{
		double t = 1.0 / ( 1.0 + p * x );
		return (1.0 - c * exp( -x * x / 2.0 ) * t * ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
	}
	else 
	{
		double t = 1.0 / ( 1.0 - p * x );
		return ( c * exp( -x * x / 2.0 ) * t * ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
	}
}

double cumulated_normal(const double x, const double mu, const double sd)
{
	return cumulated_normal(x/sd -  mu/sd);
}

int factorial(const int n)
{
	if(n < 2) return 1;
	else
	{
		int f = 2;
		for(int s = 3; s <= n;++s) f *= s;
		return f;
	}
}

double poisson(const int k, const double lambda)
{
	return (1.0/static_cast<double>(factorial(k))) * pow(lambda,static_cast<double>(k)) * exp((-1.0) * lambda);    
}

double p_value(const double x,const double mu, const double sd)
{
	return (1.0 - cumulated_normal(x-1,mu,sd));
}

double p_value(const unsigned n_u, const double lambda)
{
	double p=0.0;
	for(unsigned i = 0; i < n_u;++i) p += poisson(i,lambda); 
	return (1 - p);
}

double over(const int n, const int k)
{
	return ((double)factorial(n)/(factorial(n-k)*factorial(k)));
}

double over2(int n, int k)
{
    double cnk = 1.0;
    
    if (k * 2 > n) 
        k = n - k;
    
    for (int i = 1; i <= k; n--, i++)
    {
        cnk /= i;
        cnk *= n;
    }
    return cnk;
}


//////////////////////////////////////////////////////////////////////////////
// Output SNPs
// Simple SNP calling, under construction
template <
	typename TMatches,
	typename TGenomeSet,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames,
	typename TReadQualities,
	typename TOptions
>
void dumpSNPs(
	TMatches &matches,							// forward/reverse matches
	TGenomeSet & genomes,
	TGenomeNames const &genomeIDs,				// Read names (read from Fasta file, currently unused)
	StringSet<CharString> &genomeFileNameList,	// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & , //map to retrieve genome filename and sequence number within that file
	TReads const &reads,						// Read sequences
	TReadNames const &,					// Read names (read from Fasta file, currently unused)
	TReadQualities & readQualities,
	::std::string readFName,					// read name (e.g. "reads.fa")
	TOptions &options)
{
	typedef typename Value<TMatches>::Type		TMatch;
	typedef typename Value<TReads>::Type		TRead;
	typedef typename Value<TGenomeSet>::Type	TGenome;
	typedef typename TMatch::TGPos				TGPos;

	if(!options.hammingOnly)
	{
		::std::cout << "SNP calling only implemented for Hamming distance mapping." << ::std::endl;
		return;
	}
	
	bool useCharQual = false;
	if(length(readQualities)>0) useCharQual = true;

	// matches need to be ordered accordign to genome position
	if(options.orientationAware) 
		::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPosOaMQ<TMatch>());		
	else
		::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPosMQ<TMatch>());		

	//	options.maxHits = 1;		// only take into account best unique matches

	SEQAN_PROTIMESTART(dump_time);

	// load Genome sequences for alignment dumps
	if (length(genomes)==0)
		if(!loadGenomes(genomes, genomeFileNameList)) 
		{
			::std::cerr << "Failed to load genomes" << ::std::endl;
			return;
		}
	int gzeros = 0;
	for (unsigned l = length(genomes); l > 9; l = l / 10)
		++gzeros;

	// remove the directory prefix of readFName
	size_t lastPos = readFName.find_last_of('/') + 1;
	if (lastPos == readFName.npos) lastPos = readFName.find_last_of('\\') + 1;
	if (lastPos == readFName.npos) lastPos = 0;
	::std::string readName = readFName.substr(lastPos);

	Align<String<Dna5>, ArrayGaps> align;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
	scoreType.data_mismatch = -1;
	resize(rows(align), 2);

	::std::ofstream file;
	::std::ostringstream fileName;
	if (*options.outputSNP != 0)
		fileName << options.outputSNP;
	else
		fileName << readFName << ".snp";

	file.open(fileName.str().c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
	if (!file.is_open()) 
	{
		::std::cerr << "Failed to open output file" << ::std::endl;
		return;
	}
	
	unsigned rLen = options.rLen;
//	if(length(reads)>0) rLen = length(reads[0]);
	
	
//	maskDuplicates(matches);
//	compactMatches(matches, options);

	typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;

	TMatchIterator matchIt = begin(matches, Standard());
	TMatchIterator matchItEnd = end(matches, Standard());	
	
	double priorHet = 0.001;
	double priorHomo = (1.0 - priorHet)/2.0;
	Dna5 candidateBase;
	int quality;
	
	if(options.outputFormat == 0)
	{
		if(options.orientationAware)
			file << "chr\tpos\tref\t[A+]\t[C+]\t[G+]\t[T+]\t[A-]\t[C-]\t[G-]\t[T-]\tcov\tcall\n";
		else
			file << "chr\tpos\tref\t[A]\t[C]\t[G]\t[T]\tcov\tcall\n";
	}

	if(options.orientationAware) options.maxPile/=2;
	//collect candidate SNP positions  
	::std::set<unsigned> candidates;
	while(matchIt != matchItEnd) 
	{
		candidates.clear();	
//		posToCovMapForward.clear();	
//		posToCovMapReverse.clear();	
		unsigned currSeqNo = (*matchIt).gseqNo;
		if(options._debugLevel > 0) ::std::cout << "Scanning genome #" << currSeqNo << " for SNPs..." << ::std::endl;
		TMatchIterator currSeqMatchItBegin = matchIt;
		while(matchIt != matchItEnd)
		{
			if ((*matchIt).gseqNo != currSeqNo) break;
			if ((*matchIt).editDist > 0) 
			{
				Dna5String gInf = infix(genomes[(*matchIt).gseqNo], (*matchIt).gBegin, (*matchIt).gEnd);
//				::std::cout << gInf << "g\n";
//				::std::cout << (String<Dna5>)reads[(*matchIt).rseqNo] << "r\n";
				if ((*matchIt).orientation == 'R')
					reverseComplementInPlace(gInf);
				//	::std::cout << (Dna5String)reads[(*matchIt).rseqNo] << "r\n";	
				for (unsigned i = 0; i < length(gInf); ++i)
				{
					if ((options.compMask[ordValue((Dna5)reads[(*matchIt).rseqNo][i])] & options.compMask[ordValue(gInf[i])]) == 0)
					{
						if ((*matchIt).orientation == 'R')
							candidates.insert((*matchIt).gEnd-i-1);
						else
							candidates.insert((*matchIt).gBegin+i);
					}
				}
				//::std::cout << candidates.size() << "c\n";
			}
			++matchIt;
		}
		TMatchIterator currSeqMatchItEnd = matchIt;

		if(options._debugLevel > 0) ::std::cout << "Verify candidates...\n";
		if(options._debugLevel > 0) ::std::cout << matchItEnd-currSeqMatchItBegin << " matches left.\n";
		
		matchIt = currSeqMatchItBegin;
		typename ::std::set<unsigned>::iterator candidateIt = candidates.begin();
		typename ::std::set<unsigned>::iterator candidateItEnd = candidates.end();

		String<int> columnQualityF;
		resize(columnQualityF,5);
		String<unsigned> countF;
		resize(countF,5);
		String<CharString> qualityStringF;
		resize(qualityStringF,5);

		String<int> columnQualityR;
		resize(columnQualityR,5);
		String<unsigned> countR;
		resize(countR,5);
		String<CharString> qualityStringR;
		resize(qualityStringR,5);

		String<unsigned> count;
		resize(count,5);
		String<unsigned> columnQuality;
		resize(columnQuality,5);
		
		//bool extraV = false;
		if(options._debugLevel>0)::std::cout << candidates.size() << " candidates\n";
		for(; candidateIt != candidateItEnd; ++candidateIt)
		{

			//candidate SNP position
//			if(extraV) ::std::cout << "100\n";
//			extraV = false;
			unsigned candidatePos = *candidateIt;
			//std::cout << candidatePos << "candiadte\n";
			//if(candidatePos == 48255422) extraV = true;
//			if(extraV) ::std::cout << "3\n";

			for(unsigned t=0;t<5;++t) 
			{
				countF[t] = 0;
				columnQualityF[t] = 0;
				qualityStringF[t] = "";
				countR[t] = 0;
				columnQualityR[t] = 0;
				qualityStringR[t] = "";
			}

//			::std::cout << candidatePos<< "\n";
//			if(candidatePos > length(genomes[currSeqNo])) ::std::cerr<< "ohoh\n";
			//find range of relevant read matches
			while(matchIt != currSeqMatchItEnd && (*matchIt).gEnd <= candidatePos)
				++matchIt;
			TMatchIterator matchRangeBegin = matchIt;
			while(matchIt != currSeqMatchItEnd && (*matchIt).gBegin <= candidatePos)
				++matchIt;
			TMatchIterator matchRangeEnd = matchIt;

			int coverage = matchRangeEnd-matchRangeBegin;
			if(coverage<(int)options.minCoverage)
			{
				matchIt = matchRangeBegin;
				continue;
			
			}
//			if(extraV) ::std::cout << "5\n";
			//reference sequence nucleotide
//			::std::cout << "Match range:" << matchRangeEnd - matchRangeBegin << ::std::endl;
			matchIt = matchRangeBegin;
			
			Dna5 refBase = genomes[(*matchRangeBegin).gseqNo][candidatePos];
			if(refBase=='N') continue;
//			if(refBase > 4) ::std::cerr <<"ohohohi\n";
			unsigned lastBegin,currentBegin;
			char lastOrientation,currentOrientation;
//			lastBegin  =  length(genomes[currSeqNo]) + 1;
//			lastOrientation  = '-';//(*matchIt).orientation;
			
			while(matchIt != matchRangeEnd)
			{
				currentBegin = (*matchIt).gBegin;
				currentOrientation = (*matchIt).orientation;
				unsigned currPile = 0;
				while(matchIt != matchRangeEnd && (*matchIt).gBegin == currentBegin && ((*matchIt).orientation == currentOrientation || !options.orientationAware) && (currPile < options.maxPile || options.maxPile == 0))
				{
					++currPile;
					if ((*matchIt).orientation == 'R')
					{
						FunctorComplement<Dna5> f;
						candidateBase = f((Dna5)reads[(*matchIt).rseqNo][(*matchIt).gEnd - candidatePos - 1]);
						if(useCharQual)
							quality = ordValue(readQualities[(*matchIt).rseqNo][(*matchIt).gEnd - candidatePos - 1]) - 33;
						else
							quality = ((int)reads[(*matchIt).rseqNo][(*matchIt).gEnd - candidatePos - 1] >> 3);
						if(!options.useBaseQuality && quality > (*matchIt).mScore) quality = (*matchIt).mScore;
				//		if(ordValue(candidateBase)>4)::std::cerr << "bledsinn\n";
						columnQualityF[ordValue(candidateBase)] += quality;
						++countF[ordValue(candidateBase)];
						appendValue(qualityStringF[ordValue(candidateBase)],(char)(quality+33));
					}
					else
					{
					//	std::cout << reads[(*matchIt).rseqNo] <<"<-read"<<candidatePos - (*matchIt).gBegin<< "<-pos\n";
					//	std::cout <<(Dna5)reads[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin]<<"\n";
						candidateBase = (Dna5)reads[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin];
						if(useCharQual)
							quality = ordValue(readQualities[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin]) - 33;
						else
							quality = ((int)reads[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin] >> 3);
					//	std::cout << ordValue((Dna5)candidateBase) <<"<-ordvalue\n";
						if(!options.useBaseQuality && quality > (*matchIt).mScore) quality = (*matchIt).mScore;
				//		if(ordValue(candidateBase)>4)
				//		{
				//			::std::cerr << "cand="<<candidatePos <<" gBegin=" << (*matchIt).gBegin;
				//			::std::cerr << "gend="<< (*matchIt).gEnd << " rSeqNo=" << (*matchIt).rseqNo << "\n";
				//			::std::cerr << "bledsini2n\n";
				//			::std::cerr << "read="<<reads[(*matchIt).rseqNo]<<"\n";
				//		}
						++countR[ordValue(candidateBase)];
						columnQualityR[ordValue(candidateBase)] += quality;
						appendValue(qualityStringR[ordValue(candidateBase)],(char)(quality+33));
					}
					++matchIt;
				}
				//if(matchRangeEnd > matchItEnd) ::std::cerr <<"neeeeeeee\n";
				while(matchIt != matchRangeEnd && (*matchIt).gBegin == currentBegin && ((*matchIt).orientation == currentOrientation || !options.orientationAware))
					++matchIt;

//				::std::cout << "curB=" << currentBegin << " curO="<< currentOrientation << " curP=" << currPile << std::endl;
//				::std::cout << "lasB=" << lastBegin << " curO="<< lastOrientation << std::endl;
//				if(doTabFile)
//				{	
//					if(/*options.orientationAware &&*/ lastBegin != currentBegin)
//					{
//				//		::std::cout <<"insert 0!\n";
//						if(currentOrientation=='R')
//							posToCovMapForward.insert(std::make_pair<unsigned,short>(currentBegin,0));
//						else
//							posToCovMapReverse.insert(std::make_pair<unsigned,short>(currentBegin,0));
//					}
//					if(currentOrientation == 'F')
//					{
//						typename std::map<unsigned,short>::iterator it = posToCovMapForward.find(currentBegin);
//						if(it == posToCovMapForward.end()) 
//							posToCovMapForward.insert(std::make_pair<unsigned,short>(currentBegin,currPile));
//						else
//							it->second = currPile;
//					}
//					else
//					{
//						typename std::map<unsigned,short>::iterator it = posToCovMapReverse.find(currentBegin);
//						if(it == posToCovMapReverse.end()) 
//							posToCovMapReverse.insert(std::make_pair<unsigned,short>(currentBegin,currPile));
//						else
//							it->second = currPile;
//					}
//				}
				lastBegin  = currentBegin;
				lastOrientation  = currentOrientation;
			}
			unsigned realCoverageF = countF[0] + countF[1] +countF[2] +countF[3] +countF[4];
			unsigned realCoverageR = countR[0] + countR[1] +countR[2] +countR[3] +countR[4];
			unsigned realCoverage = realCoverageF + realCoverageR;

			for(unsigned k=0;k<5;++k)
			{
				columnQuality[k] = columnQualityF[k]+columnQualityR[k];
				count[k] = countF[k]+countR[k];
			}

/*		if(realCoverage<options.minCoverage)
			{
				::std::cout << "Coverage " << realCoverage << " after applying maximum of 2 reads per genome pos" << ::std::endl;
				matchIt = matchRangeBegin;
				continue;
			}*/
			//do statistics and output snp
			int allele1 = -1;
			int allele2 = -1;
			int errorAllele = -1;
			
//			::std::cout << "refAllele\n";
			int refAllele = ordValue(genomes[(*matchIt).gseqNo][candidatePos]);
			unsigned maxCount=0;
			for(int k=0; k < 5; ++k)
			{
				if(count[k] > maxCount)
				{
					maxCount = count[k];
					allele1 = k;
				}
			}
			maxCount = 0;
			for(int k=0; k < 5; ++k)
			{
				if(k != allele1 && count[k] > maxCount)
				{
					maxCount = count[k];
					allele2 = k;
				}
			}
			int n = count[allele1];
			if(allele1 < 0) // reads all vote for one nucleotide (non-ref)
				allele1 = refAllele;
			if(allele2 < 0) // reads all vote for one nucleotide (non-ref)
				allele2 = refAllele;
			else n += count[allele2]; // ignore all but the two most common bases
			if(allele1==refAllele && allele2==refAllele)
			{
				if(options._debugLevel > 1)::std::cout << "No evidence of non-ref bases left after applying max pile filter.\n" << ::std::endl;
				matchIt = matchRangeBegin;
				//				::std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1Happens.!!!!!!!!!!!!!!!!!!!!!!!!1\n";
				continue;
			}
			if(allele1==refAllele)
				errorAllele=allele2;
			else
				errorAllele=allele1;
		
//			if(errorAllele <0 || allele1 <0 || allele2 <0 ||refAllele<0)
//				::std::cerr << "keckig\n";
			
//			if(errorAllele > 4 || allele1 > 4 || allele2 > 4 ||refAllele>4)
//				::std::cerr << "keckig\n";
			double avgErrProb = 0.08;
			// TODO: minimum coverage on the two most common alleles??
			/*if(n<10)
			{
				matchIt = matchRangeBegin;
				continue;
			
			}*/
			
			// 
			// argmax P(g|D)=P(D|g)*P(g)/P(D)
			//    g
			// 
			//todo: alphas so wie in maq
			double pHomoAllele1 = 0;//alpha_n_k;
			double pHomoAllele2 = 0;//alpha_n_nk;
			double pHet = 0;//...
			bool hetSNP = false;
			bool homoSNP = false;
			
		//	std::cout << "calling snp?\n";
			if(options.method==0)
			{
				int mutAllele=allele1;
				if(allele1==refAllele) mutAllele=allele2;
			//	if(mutAllele > 4 || mutAllele < 0) ::std::cerr <<"geht nicht\n";
				homoSNP=false; 
				hetSNP= false;
				if(count[mutAllele] > options.minMutT && (float)count[mutAllele]/realCoverage > options.percentageT && (float)(columnQualityF[mutAllele]+columnQualityR[mutAllele])/count[mutAllele] > options.avgQualT) homoSNP=true;
			}
			else
			{
				__int64 possib = (1 << n); //vorsicht!
				if(n>30) possib = (1 << 30); //temporary solution to circumvent overflows
			
				pHet = over2(n,count[allele2])/possib; 
				//double pHomoAllele1 = alpha_n_k;
				//double pHomoAllele1 = binomial(n,count[allele2],avgErrProb);//uniform and independent
				pHomoAllele1 = over2(n,count[allele2]) * pow(avgErrProb,(double)count[allele2]) * pow(1.0-avgErrProb,(double)count[allele1]) ;//uniform and independent
				//double pHomoAllele2 = alpha_n_nk;
				//double pHomoAllele2 = binomial(n,count[allele1],avgErrProb);//uniform and independent
				pHomoAllele2 = over2(n,count[allele1]) * pow(avgErrProb,(double)count[allele1]) * pow(1.0-avgErrProb,(double)count[allele2]);//uniform and independent

				pHet *= priorHet; 
				pHomoAllele1 *= priorHomo;
				pHomoAllele2 *= priorHomo;
			
				if(pHomoAllele1 > pHomoAllele2)
				{
					if(pHet > pHomoAllele1)
						hetSNP = true;
					else
						if(allele1 != refAllele)
							homoSNP = true;
				}
				else
				{
					if(pHet > pHomoAllele2)
						hetSNP = true;
					else
						if(allele2 != refAllele)
							homoSNP = true;
					if(homoSNP && count[allele2] != count[allele1]) ::std::cout << "Vorsicht! HomoSnp auf allele2!\n";
				}
	//			if(allele1 != refAllele && allele2 != refAllele)
	//				::std::cout <<"+";
				//double lambda = coverage * 0.05;//=avg_err_prob;
				//double pHet = p_value(count[allele2],lambda);
				//if(coverage > 4 && pValHet <= options.testLevel) // heterozygote SNP
				//	hetSNP = true;
				//if(coverage > 4 && pValHet > options.testLevel && (unsigned)allele1 != ordValue(genomes[(*matchIt).gseqNo][candidatePos])) // homozygote SNP
				//	homoSNP = true;
				
				//temporär!!!!
				//hetSNP = false; homoSNP = true;
				//temporär!!!!
				if(hetSNP && (int)(count[allele1]/count[allele2])>=3) hetSNP = false;
				if(allele1 == refAllele) homoSNP = false;
			}
//			::std::cout << "outputtung\n";
			if(options.outputFormat == 0) //
			{
				//chromosome
				file << genomeIDs[currSeqNo] << '\t';
				file << candidatePos + options.positionFormat<< '\t';
				file << (Dna5)genomes[currSeqNo][candidatePos] <<'\t';
				if(options.orientationAware)
				{
					file << "["<<qualityStringF[0] <<"]\t";
					file << "["<<qualityStringF[1] <<"]\t";
					file << "["<<qualityStringF[2] <<"]\t";
					file << "["<<qualityStringF[3] <<"]\t";
					file << "["<<qualityStringR[0] <<"]\t";
					file << "["<<qualityStringR[1] <<"]\t";
					file << "["<<qualityStringR[2] <<"]\t";
					file << "["<<qualityStringR[3] <<"]\t";
				}
				file << "["<<qualityStringF[0]<<qualityStringR[0] <<"]\t";
				file << "["<<qualityStringF[1]<<qualityStringR[1] <<"]\t";
				file << "["<<qualityStringF[2]<<qualityStringR[2] <<"]\t";
				file << "["<<qualityStringF[3]<<qualityStringR[3] <<"]\t";
				file << realCoverage;
				int mutAllele=allele1;
				if(allele1==refAllele) mutAllele=allele2;
				if(homoSNP) file  << '\t' << (Dna)mutAllele;
				file << std::endl;
			}
			else if(hetSNP || homoSNP)
			{
				//
				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[currSeqNo] << " ";
						break;
		//			// 1..filename is the read filename + seqNo
		//			case 1:
		//				file.fill('0');
		//				file << gnoToFileMap[currSeqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[currSeqNo].second + 1 << " ";
				}
				file << candidatePos + options.positionFormat<< " " << genomes[(*matchIt).gseqNo][candidatePos] <<" ";
				file << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << " " << count[4];
				file <<  " "<< pHomoAllele1 << " " << pHet << " " << coverage <<" " << realCoverage << " " << count[errorAllele];
				file << " " << columnQuality[errorAllele] << " " << columnQuality[0] <<" " << columnQuality[1] <<" " << columnQuality[2] <<" " << columnQuality[3] << ::std::endl;
//				file << candidatePos << " " << genomes[(*matchIt).gseqNo][candidatePos] <<" ";
//				file << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << " " << count[4];
//				if(homoSNP) file <<  " <" << Dna(allele1) << "," << Dna(allele1) << "> "  << pHomoAllele1; //vorsicht!
//				else file <<  " <" << Dna(allele1) << "," << Dna(allele2) << "> "<< pHet;
//				file << " " << columnQuality << ::std::endl;
//				
				if(options._debugLevel > 0) 
				{
					::std::cout << "CandidateSNP: " << candidatePos + options.positionFormat << " (covered by " << coverage << " matches):" << ::std::endl;
					::std::cout << "Genome is " << genomes[(*matchIt).gseqNo][candidatePos];
					::std::cout << "\treads are "<<count[0] << "," << count[1] << "," << count[2] << "," << count[3] << "," << count[4];
					if(homoSNP) ::std::cout <<  " <" << Dna(allele1) << "," << Dna(allele1) << "> " << pHomoAllele1; //vorsicht!
					else ::std::cout <<  " <" << Dna(allele1) << "," << Dna(allele2) << "> " << pHet;
					::std::cout << " " << columnQuality[0] << ::std::endl;
				}
			}
			matchIt = matchRangeBegin;
//			if(extraV) ::std::cout << "9\n";
		}

//		if (doTabFile)
//		{
//			//::std::cout << "tabfile\n";
//			typename std::map<unsigned,short>::iterator posToCovMapFIt = posToCovMapForward.begin();
//			typename std::map<unsigned,short>::iterator posToCovMapRIt = posToCovMapReverse.begin();
//			std::cout << posToCovMapReverse.size() << "Rev"<< posToCovMapForward.size()<<"For" << "NEIN!!\n";
//			if(posToCovMapReverse.size() != posToCovMapForward.size()) std::cout << "NEIN!!\n"<<std::flush;
//			else
//			{
//				for(; posToCovMapRIt != posToCovMapReverse.end() && posToCovMapFIt != posToCovMapForward.end(); ++posToCovMapRIt, ++posToCovMapFIt)
//				{
//					tabfile <<  genomeIDs[currSeqNo] << '\t';
//					tabfile <<  posToCovMapFIt->first + options.positionFormat<< '\t';
//					tabfile <<  posToCovMapFIt->first + rLen << '\t';
//					if(options.orientationAware)
//					{
//						tabfile <<  posToCovMapFIt->second << '\t';
//						tabfile <<  posToCovMapRIt->second << std::endl;
//					}
//					else tabfile <<  posToCovMapRIt->second + posToCovMapFIt->second << std::endl;
//				}
//			}
//		}


//		if(extraV) ::std::cout << "1000\n";
		matchIt = currSeqMatchItEnd;
		if(options._debugLevel>0) std::cout <<"Finished scanning chromosome.\n"<<std::flush;
	}
	file.close();

	if (*options.tabFile != 0)
	{
		::std::ofstream tabfile;
		std::map<unsigned,short> posToCovMapForward;
		std::map<unsigned,short> posToCovMapReverse;
		tabfile.open(options.tabFile, ::std::ios_base::out | ::std::ios_base::trunc);
		if (!tabfile.is_open()) 
		{
			::std::cerr << "Failed to open output tab file" << ::std::endl;
			return;
		}
		if(options.orientationAware)
			tabfile << "chr\tbegin\tend\tcount+\tcount-"<<std::endl;
		else
			tabfile << "chr\tbegin\tend\tcount"<<std::endl;
		matchIt = begin(matches,Standard());
		while(matchIt != matchItEnd) 
		{
			posToCovMapForward.clear();	
			posToCovMapReverse.clear();	
			unsigned currSeqNo = (*matchIt).gseqNo;
			if(options._debugLevel > 0) ::std::cout << "Scanning genome #" << currSeqNo << " for coverage..." << ::std::endl;
			TMatchIterator currSeqMatchItBegin = matchIt;
			while(matchIt != matchItEnd)
			{
				if ((*matchIt).gseqNo != currSeqNo) break;
				++matchIt;
			}
			TMatchIterator currSeqMatchItEnd = matchIt;
//			bool extraV = false;
		
			matchIt = currSeqMatchItBegin;
			unsigned lastBegin,currentBegin;
			char lastOrientation,currentOrientation;
			lastBegin  =  length(genomes[currSeqNo]) + 1;
			lastOrientation  = '-';//(*matchIt).orientation;
			
			while(matchIt != currSeqMatchItEnd)
			{
				currentBegin = (*matchIt).gBegin;
				currentOrientation = (*matchIt).orientation;
//				if(currentBegin==215293) extraV = true;
//				else extraV = false;
				unsigned currPile = 0;
				while(matchIt != currSeqMatchItEnd && (*matchIt).gBegin == currentBegin && ((*matchIt).orientation == currentOrientation || !options.orientationAware) && (currPile < options.maxPile || options.maxPile == 0))
				{
//					if(extraV) std::cout << "++currPile\n";
					++currPile;
					++matchIt;
				}
				while(matchIt != currSeqMatchItEnd && (*matchIt).gBegin == currentBegin && ((*matchIt).orientation == currentOrientation || !options.orientationAware))
				{
					++matchIt;
//					if(extraV) std::cout << "++currPile\n";
				}

				if(/*options.orientationAware &&*/ lastBegin != currentBegin)
				{
//					if(extraV) ::std::cout <<"insert 0!\n";
					if(currentOrientation=='R')
						posToCovMapForward.insert(std::make_pair<unsigned,short>(currentBegin,0));
					else
						posToCovMapReverse.insert(std::make_pair<unsigned,short>(currentBegin,0));
				}
				if(currentOrientation == 'F')
				{
					typename std::map<unsigned,short>::iterator it = posToCovMapForward.find(currentBegin);
					if(it == posToCovMapForward.end()) 
						posToCovMapForward.insert(std::make_pair<unsigned,short>(currentBegin,currPile));
					else
						it->second = currPile;
				}
				else
				{
					typename std::map<unsigned,short>::iterator it = posToCovMapReverse.find(currentBegin);
					if(it == posToCovMapReverse.end()) 
						posToCovMapReverse.insert(std::make_pair<unsigned,short>(currentBegin,currPile));
					else
						it->second = currPile;
				}
//				if(extraV) std::cout << "lastBeg=" << lastBegin << " curBeg=" << currentBegin << "\n";
				lastBegin  = currentBegin;
				lastOrientation  = currentOrientation;
			}
			typename std::map<unsigned,short>::iterator posToCovMapFIt = posToCovMapForward.begin();
			typename std::map<unsigned,short>::iterator posToCovMapRIt = posToCovMapReverse.begin();
//			std::cout << posToCovMapReverse.size() << "Rev"<< posToCovMapForward.size()<<"For" << "NEIN!!\n";
			if(posToCovMapReverse.size() != posToCovMapForward.size()) std::cout << "NEIN!!\n"<<std::flush;
			else
			{
				for(; posToCovMapRIt != posToCovMapReverse.end() && posToCovMapFIt != posToCovMapForward.end(); ++posToCovMapRIt, ++posToCovMapFIt)
				{
					tabfile <<  genomeIDs[currSeqNo] << '\t';
					tabfile <<  posToCovMapFIt->first + options.positionFormat<< '\t';
					tabfile <<  posToCovMapFIt->first + rLen << '\t';
					if(options.orientationAware)
					{
						tabfile <<  posToCovMapFIt->second << '\t';
						tabfile <<  posToCovMapRIt->second << std::endl;
					}
					else tabfile <<  posToCovMapRIt->second + posToCovMapFIt->second << std::endl;
				}
			}

			matchIt = currSeqMatchItEnd;
			if(options._debugLevel>1) std::cout <<"Finished writing tabfile for sequence "<<currSeqNo<<".\n"<<std::flush;
		}
		if(options._debugLevel>0) std::cout <<"Finished writing tabfile.\n"<<std::flush;
	
		tabfile.close();
	//	posToCovMapForward.clear();
	//	posToCovMapReverse.clear();
	}



	if(options._debugLevel >= 1) ::std::cout << "Detecting SNPs took " << SEQAN_PROTIMEDIFF(dump_time) << " seconds." << ::std::endl;
	return;

}







///////////////////////////////////////////////////////////////////////////////////////////////7
// simple position stats analysis
template <typename TPositions, typename TOptions>
bool loadPositions(TPositions & positions, char const * filename, TOptions & options)
{
	::std::ifstream file;
	file.open(filename,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return 1;	

	clear(positions);
	char c = _streamGet(file);
	while(!_streamEOF(file))
	{
		//currPos = _parse_readNumber(file);
 		appendValue(positions,_parse_readNumber(file,c)-options.positionFormat);
		_parse_skipLine(file,c);
	}
	file.close();
	return (length(positions) > 0);
}


//////////////////////////////////////////////////////////////////////////////
// check out snp positions
template <
	typename TMatches,
	typename TPositions,
	typename TGenomeSet,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames,
	typename TReadQualities,
	typename TOptions
>
void dumpPositionStats(
	TMatches &matches,						// forward/reverse matches
	TPositions & positions,
	TGenomeSet & genomes,
	TGenomeNames const &,						// Read names (read from Fasta file, currently unused)
	StringSet<CharString> &,					// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & ,	//map to retrieve genome filename and sequence number within that file
	TReads const &reads,						// Read sequences
	TReadNames const &,						// Read names (read from Fasta file, currently unused)
	TReadQualities & readQualities,
	::std::string readFName,					// read name (e.g. "reads.fa")
	TOptions &options)
{
	typedef typename Value<TMatches>::Type		TMatch;
	typedef typename Value<TReads>::Type		TRead;
	typedef typename Value<TGenomeSet>::Type	TGenome;
	typedef typename TMatch::TGPos				TGPos;

	// matches need to be ordered accordign to genome position
	::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPosMQ<TMatch>());		

	//	options.maxHits = 1;		// only take into account best unique matches
	bool useCharQual = false;
	if(length(readQualities)>0)
		useCharQual = true;
	SEQAN_PROTIMESTART(dump_time);
	
	
	//if (!loadGenomes(genomes, genomeFileNameList)) 
	//{
	//	::std::cerr << "Failed to load genomes" << ::std::endl;
	//	return;
	//}
	int gzeros = 0;
	for (unsigned l = length(genomes); l > 9; l = l / 10)
		++gzeros;

	// remove the directory prefix of readFName
	size_t lastPos = readFName.find_last_of('/') + 1;
	if (lastPos == readFName.npos) lastPos = readFName.find_last_of('\\') + 1;
	if (lastPos == readFName.npos) lastPos = 0;
	::std::string readName = readFName.substr(lastPos);

	::std::ofstream file;
	::std::ostringstream fileName;
	if (*options.output != 0)
		fileName << options.output;
	else
		fileName << readFName << ".stats";

	file.open(fileName.str().c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
	if (!file.is_open()) 
	{
		::std::cerr << "Failed to open output file" << ::std::endl;
		return;
	}
	
	
	typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;
	typedef typename Iterator<TPositions, Standard>::Type TPosIterator;
	
	TMatchIterator matchIt = begin(matches, Standard());
	TMatchIterator matchItEnd = end(matches, Standard());	
	
	TPosIterator candidateIt = begin(positions, Standard());
	TPosIterator candidateItEnd = end(positions, Standard());
	String<int> countForward;
	fill(countForward,5,0);
	String<int> countReverse;
	fill(countReverse,5,0);
	String<float> columnQualityForward;
	fill(columnQualityForward,5,0.0); //qualitySum for each base
	String<float> columnQualityReverse;
	fill(columnQualityReverse,5,0.0); //qualitySum for each reverse base
	Dna5 candidateBase;
	int quality;
		
	for(; candidateIt != candidateItEnd; ++candidateIt)
	{
		unsigned candidatePos = *candidateIt;
		if(options._debugLevel > 1) std::cout <<"candidatePos="<<candidatePos<<"\n";
		for(unsigned t=0; t < 5; ++t)
		{
			countForward[t]=0;
			countReverse[t]=0;
			columnQualityForward[t]=0.0;
			columnQualityReverse[t]=0.0;
		}
		while(matchIt != matchItEnd && (*matchIt).gEnd <= candidatePos)
			++matchIt;
		TMatchIterator matchRangeBegin = matchIt;
		while(matchIt != matchItEnd && (*matchIt).gBegin <= candidatePos)
			++matchIt;
		TMatchIterator matchRangeEnd = matchIt;
		
		//unsigned coverage = matchRangeEnd-matchRangeBegin;
		
//		::std::cout << "Match range:" << matchRangeEnd - matchRangeBegin << ::std::endl;
		
		matchIt = matchRangeBegin;
		//reference sequence nucleotide
		Dna5 refBase = genomes[(*matchRangeBegin).gseqNo][candidatePos];
		if(refBase=='N') continue;
		
		
		while(matchIt != matchRangeEnd)
		{
			unsigned currentBegin = (*matchIt).gBegin;
			unsigned currPile = 0;
			while(matchIt != matchRangeEnd && (*matchIt).gBegin == currentBegin && (currPile < options.maxPile || options.maxPile ==0))
			{
				++currPile;
				if ((*matchIt).orientation == 'R')
				{
					FunctorComplement<Dna5> f;
					candidateBase = f((Dna5)reads[(*matchIt).rseqNo][(*matchIt).gEnd - candidatePos - 1]);
					if(useCharQual)
						quality = ordValue(readQualities[(*matchIt).rseqNo][(*matchIt).gEnd - candidatePos - 1]) - 33;
					else
						quality = ((int)reads[(*matchIt).rseqNo][(*matchIt).gEnd - candidatePos - 1] >> 3);
					if(!options.useBaseQuality && quality > (*matchIt).mScore) quality = (*matchIt).mScore;
					columnQualityReverse[ordValue(candidateBase)] += quality;
					++countReverse[ordValue(candidateBase)];
				}
				else
				{
				//	std::cout << reads[(*matchIt).rseqNo] <<"<-read"<<candidatePos - (*matchIt).gBegin<< "<-pos\n";
				//	std::cout <<(Dna5)reads[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin]<<"\n";
					candidateBase = (Dna5)reads[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin];
					if(useCharQual)
						quality = ordValue(readQualities[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin]) - 33;
					else
						quality = ((int)reads[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin] >> 3);
				//	std::cout << ordValue((Dna5)candidateBase) <<"<-ordvalue\n";
					++countForward[ordValue((Dna5)candidateBase)];
					if(!options.useBaseQuality && quality > (*matchIt).mScore) quality = (*matchIt).mScore;
					columnQualityForward[ordValue((Dna5)candidateBase)] += quality;
				}
				++matchIt;
			}
		}
		matchIt = matchRangeBegin;
		
		unsigned realCoverageForward = countForward[0] + countForward[1] +countForward[2] +countForward[3] +countForward[4];
		unsigned realCoverageReverse = countReverse[0] + countReverse[1] +countReverse[2] +countReverse[3] +countReverse[4];
	
		//unsigned qualitySumForward = columnQualityForward[0] +columnQualityForward[1] +columnQualityForward[2]+columnQualityForward[3]+columnQualityForward[4];
		//unsigned qualitySumReverse = columnQualityReverse[0] +columnQualityReverse[1] +columnQualityReverse[2]+columnQualityReverse[3]+columnQualityReverse[4];
		unsigned realCov = realCoverageForward + realCoverageReverse;
		
		if(realCov < options.extractMinCov || realCov >= options.extractMaxCov)
			continue;
		
//		file << genomeIDs[currSeqNo] << " ";
		file << candidatePos + options.positionFormat<<  "\t" << genomes[(*matchIt).gseqNo][candidatePos] <<"\t";
	//	if((float)(countForward[ordValue(genomes[(*matchIt).gseqNo][candidatePos])] + countReverse[ordValue(genomes[(*matchIt).gseqNo][candidatePos])]) /(realCoverageForward + realCoverageReverse) < 0.5) std::cout << candidatePos << "\n";
		if(options.orientationAware)
		{
			if(realCoverageForward ==0) 
				file << "0\t0\t0\t0\t0\t0\t0\t0 ";
			else
			{
				file << countForward[0] << "\t";
				file << countForward[1] << "\t";
				file << countForward[2] << "\t";
				file << countForward[3] << "\t";
				file << columnQualityForward[0]/realCoverageForward << "\t";
				file << columnQualityForward[1]/realCoverageForward << "\t";
				file << columnQualityForward[2]/realCoverageForward << "\t";
				file << columnQualityForward[3]/realCoverageForward << "\t";
			}
			if(realCoverageReverse ==0) 
				file << "0\t0\t0\t0\t0\t0\t0\t0" << std::endl;
			else
			{
				file << countReverse[0] << "\t";
				file << countReverse[1] << "\t";
				file << countReverse[2] << "\t";
				file << countReverse[3] << "\t";
				file << columnQualityReverse[0]/realCoverageReverse << "\t";
				file << columnQualityReverse[1]/realCoverageReverse << "\t";
				file << columnQualityReverse[2]/realCoverageReverse << "\t";
				file << columnQualityReverse[3]/realCoverageReverse << std::endl;
			}
		}
		else {
			//int qualitySum = qualitySumForward + qualitySumReverse;
			if(realCov==0) 
				file << "0\t0\t0\t0\t0\t0\t0\t0" << std::endl;
			else
			{
				if(options.schnellUndDreckig)
				{
					unsigned snpAllele = 0;
					unsigned refAllele = ordValue(genomes[(*matchIt).gseqNo][candidatePos]);
					unsigned countRefAllele = countForward[refAllele]+countReverse[refAllele];
					unsigned countSnpAllele = 0;
					for(unsigned k = 0; k < 4; ++k)
					{
						if(k!= refAllele)
							if(countForward[k]+countReverse[k] > (int)countSnpAllele)
							{
								snpAllele = k;
								countSnpAllele = countForward[k]+countReverse[k];
							}
					}
					file << "\t"<<(float)(columnQualityForward[refAllele]+columnQualityReverse[refAllele])/(countRefAllele) <<"\t";
					file <<( Dna5)snpAllele<<"\t" << "\t"<<(float)(columnQualityForward[snpAllele]+columnQualityReverse[snpAllele])/(countSnpAllele) <<"\t" << ((float)countSnpAllele/realCov) <<std::endl;
					
				}
				else
				{
				file << countForward[0] + countReverse[0] << "\t";
				file << countForward[1] + countReverse[1] << "\t";
				file << countForward[2] + countReverse[2] << "\t";
				file << countForward[3] + countReverse[3]<< "\t";
				// per-obseved-base quality amount
/*				file << (float)(columnQualityForward[0]+columnQualityReverse[0])/realCov << "\t";
				file << (float)(columnQualityForward[1]+columnQualityReverse[1])/realCov << "\t";
				file << (float)(columnQualityForward[2]+columnQualityReverse[2])/realCov << "\t";
				file <<(float)(columnQualityForward[3]+columnQualityReverse[3])/realCov << std::endl;*/
				// fractional quality amount
/*				file << (int)(100*(columnQualityForward[0]+columnQualityReverse[0])/qualitySum) << "\t";
				file << (int)(100*(columnQualityForward[1]+columnQualityReverse[1])/qualitySum) << "\t";
				file << (int)(100*(columnQualityForward[2]+columnQualityReverse[2])/qualitySum) << "\t";
				file << (int)(100*(columnQualityForward[3]+columnQualityReverse[3])/qualitySum) << std::endl;*/
				// absolute quality sum
				file << (columnQualityForward[0]+columnQualityReverse[0]) << "\t";
				file << (columnQualityForward[1]+columnQualityReverse[1]) << "\t";
				file << (columnQualityForward[2]+columnQualityReverse[2]) << "\t";
				file << (columnQualityForward[3]+columnQualityReverse[3]) << std::endl;
				// average nucleotide quality
/*				if(countForward[0] + countReverse[0] > 0)
					file << (float)(columnQualityForward[0]+columnQualityReverse[0])/(countForward[0] + countReverse[0]) << "\t";
				else file <<"0\t";
				if(countForward[1] + countReverse[1] > 0)
					file << (float)(columnQualityForward[1]+columnQualityReverse[1])/(countForward[1] + countReverse[1]) << "\t";
				else file <<"0\t";
				if(countForward[2] + countReverse[2] > 0)
					file << (float)(columnQualityForward[2]+columnQualityReverse[2])/(countForward[2] + countReverse[2])<< "\t";
				else file <<"0\t";
				if(countForward[3] + countReverse[3] > 0)
					file <<(float)(columnQualityForward[3]+columnQualityReverse[3])/(countForward[3] + countReverse[3])<< std::endl;
				else file <<std::endl;
*/				}
			}
		}
	}
//	clear(countReverse);clear(countForward);clear(columnQualityForward);clear(columnQualityReverse);
	//clear(countForward);clear(countForward);
	file.close();
	if(options._debugLevel >= 1) ::std::cout << "Inspecting positions took " << SEQAN_PROTIMEDIFF(dump_time) << " seconds." << ::std::endl;
	return;

}




}

#endif
