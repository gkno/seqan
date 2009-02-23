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

//#define RAZERS_CONCATREADS		// use <ConcatDirect> StringSet to store reads
#define RAZERS_DUMP_SNPS







#include "seqan/platform.h"
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/graph_types/graph_utility_parsing.h>

#ifdef PLATFORM_WINDOWS
	#define SEQAN_DEFAULT_TMPDIR "C:\\TEMP\\"
#else
	#define SEQAN_DEFAULT_TMPDIR "./"
#endif

//#include "../razers/razers_utils.h"
#include "../razers/razers.h"
#include "../razers/outputFormat.h"
#include "callSNPs.h"


#include <iostream>
#include <sstream>
#include <map>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////
// Main read mapper function
template <typename TSpec>
int detectSNPs(
	const char *genomeFileName,
	String<CharString> & readFNames,
	String<CharString> & qualityFNames,
	SNPCallingOptions<TSpec> &options)
{
	typedef Dna5String					TGenome;
	typedef StringSet<TGenome>				TGenomeSet;
	typedef String<Dna5>					TRead;
#ifdef RAZERS_CONCATREADS
	typedef StringSet<TRead, Owner<ConcatDirect<> > >	TReadSet;
#else
	typedef StringSet<TRead>				TReadSet;
#endif
	typedef StringSet<CharString>				TReadQualities;
	typedef MappedReadMatch<Difference<TGenome>::Type>	TMatch;		// a single match
	typedef String<TMatch/*, MMap<>*/ >			TMatches;	// array of matches

	typedef ::std::map<CharString,unsigned> TGenomeMap;
	typedef typename TGenomeMap::iterator TMapIter;

	TReadSet			reads;
	StringSet<CharString>		readNames;		// read names, taken from the Fasta file
	TGenomeSet 			genomes;
	StringSet<CharString> 		genomeFileNameList;
	StringSet<CharString> 		genomeNames;
	TGenomeMap			gIdStringToIdNumMap;
	TMatches			matches;		// resulting forward/reverse matches
	TReadQualities			readQualities;
	// dump configuration in verbose mode
	if (options._debugLevel >= 1) 
	{
		::std::cerr << "___SETTINGS____________" << ::std::endl;
		::std::cerr << "Genome file:                             \t" << genomeFileName << ::std::endl;
		::std::cerr << "Read files:                              \t" << readFNames[0] << ::std::endl;
		for(unsigned i = 1; i < length(readFNames); ++i)
			::std::cerr << "                                         \t" << readFNames[i] << ::std::endl;
		::std::cerr << "MaxPile:                                 \t" << options.maxPile << ::std::endl;
		::std::cerr << "MinCoverage:                             \t" << options.minCoverage << ::std::endl;
		::std::cerr << "MinMutThreshold:                         \t" << options.minMutT << ::std::endl;
		::std::cerr << "MinPercentageThreshold:                  \t" << options.percentageT << ::std::endl;
		::std::cerr << "MinQualityThreshold:                     \t" << options.avgQualT << ::std::endl;
		::std::cerr << ::std::endl;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// Step 1: Determine genome file type and load genomes
	SEQAN_PROTIMESTART(load_time);

	int result = getGenomeFileNameList(genomeFileName, genomeFileNameList, options);
	if(result == CALLSNPS_GENOME_FAILED || !loadGenomes(genomes, genomeFileNameList,gIdStringToIdNumMap))
	{
		::std::cerr << "Failed to open genome file " << genomeFileName << ::std::endl;
		return result;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Step 2: Load reads and matches
	for(unsigned i = 0; i < length(readFNames); ++i)
	{
		if(options.inputFormat ==0)
			result = readMatchesFromGFF(matches,gIdStringToIdNumMap,genomes,reads,readNames,readQualities,toCString(readFNames[i]),options);
		else
		{
			if (options.inputFormat ==1)
			{
				String<bool> skip;
				result = readMatchesFromEland(matches,gIdStringToIdNumMap,genomes,reads,readNames,toCString(readFNames[i]), skip, options);
				options.useBaseQuality = true;
				// Step 2a: Load read qualities
				if(result != CALLSNPS_GFF_FAILED && options.qualityFile)
					result = readQualityValues(reads,readNames,readQualities,toCString(qualityFNames[i]),skip,options);
				if(result == CALLSNPS_QUALITY_FAILED)
				{
					::std::cerr << "Failed to open quality file " << qualityFNames[i] << ::std::endl;
					return result;
				}
			}
			else 
				result = readMatchesFromMaq(matches,gIdStringToIdNumMap,genomes,reads,readNames,readQualities,toCString(readFNames[i]),options);
		}
		if(result == CALLSNPS_GFF_FAILED)
		{
			::std::cerr << "Failed to open read file " << readFNames[i] << ::std::endl;
			return result;
		}
		
	}
	if (options._debugLevel >= 1) 
		::std::cerr << lengthSum(readQualities) << " chars of " << length(readQualities) << " read qualities loaded." << ::std::endl;
	//update matches with average read quality information
	addReadQualityToMatches(matches,reads,readQualities,options);
	
	if (options._debugLevel >= 1) 
		::std::cerr << lengthSum(reads) << " bps of " << length(reads) << " reads loaded." << ::std::endl;
	options.timeLoadFiles = SEQAN_PROTIMEDIFF(load_time);
	
	if (options._debugLevel >= 1)
		::std::cerr << "Loading reads took               \t" << options.timeLoadFiles << " seconds" << ::std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Step 3: Do SNP calling
	resize(genomeNames,gIdStringToIdNumMap.size()); //prepare genomeNames
	TMapIter gIt=gIdStringToIdNumMap.begin();
	for(unsigned i=0; i < length(genomeNames); ++i,++gIt)
		genomeNames[i] = gIt->first;
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > dummy; 
	dumpSNPs(matches, genomes, genomeNames, genomeFileNameList, dummy, reads, readNames, readQualities,toCString(readFNames[0]), options);
	
	if(*options.positionFile != 0)
	{
		String<unsigned> positions;
		result = loadPositions(positions,options.positionFile,options);
		if(result == CALLSNPS_POS_FAILED)
		{
			::std::cerr << "Failed to open position file " << options.positionFile << ::std::endl;
			return result;
		}
		dumpPositionStats(matches,positions, genomes, genomeNames, genomeFileNameList, dummy, reads, readNames, readQualities, toCString(readFNames[0]), options);

	}

	if(*options.errorPrbFileName != 0)
	{
		::std::ofstream file;
		file.open(options.errorPrbFileName, ::std::ios_base::out | ::std::ios_base::trunc);
		if (file.is_open())
		{
			String<int> posError;
			if(length(reads)> 0)fill(posError,length(reads[0]),0);
			
			int unique = getErrorDistribution(posError, matches, reads, genomes, options);
			
			file << (double)posError[0] / (double)unique;
			for (unsigned i = 1; i < length(posError); ++i)
				file << '\t' << (double)posError[i] / (double)unique;
			file << ::std::endl;
			file.close();
		} else
			::std::cerr << "Failed to open error distribution file" << ::std::endl;
		
	}


	return 0;
}	


template<typename TMatches,typename TReads,typename TReadQualities, typename TOptions>
void
addReadQualityToMatches(TMatches &matches, TReads &reads, TReadQualities &readQualities, TOptions &)
{
	typedef typename Value<TReads>::Type				TRead;
	typedef typename Value<TReadQualities>::Type			TReadQuality;
	typedef typename Value<TMatches>::Type				TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());

	int avgRQ;
	for (; it != itEnd; ++it) 
	{
		TRead &read = reads[(*it).rseqNo];
		TReadQuality &readQuality = readQualities[(*it).rseqNo];
		avgRQ = 0;
		for(unsigned i = 0; i < length(read); ++i)
			avgRQ += readQuality[i];
		(*it).avgQuality = avgRQ/length(read);
	}

}



template<typename TOptions>
int getGenomeFileNameList(char const * filename, StringSet<CharString> & genomeFileNames, TOptions &options)
{
	::std::ifstream file;
	file.open(filename,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GENOME_FAILED;
	//	return 1;
	//	return RAZERS_GENOME_FAILED;

	CharString nameStr;
	char c = _streamGet(file);
	if (c != '>' && c != '@')	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		if(options._debugLevel >=1)
			::std::cout << ::std::endl << "Reading multiple genome files:" <<::std::endl;
/*		//locations of genome files are relative to list file's location
		::std::string tempGenomeFile(filename);
		size_t lastPos = tempGenomeFile.find_last_of('/') + 1;
		if (lastPos == tempGenomeFile.npos) lastPos = tempGenomeFile.find_last_of('\\') + 1;
		if (lastPos == tempGenomeFile.npos) lastPos = 0;
		::std::string filePrefix = tempGenomeFile.substr(0,lastPos);*/
		unsigned i = 1;
		while(!_streamEOF(file))
		{
			clear(nameStr);
			_parse_skipWhitespace(file, c);
			while ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r'))
			{
				appendValue(nameStr,c);
				c = _streamGet(file);
			}
			appendValue(genomeFileNames,nameStr);
			//CharString currentGenomeFile(filePrefix);
			//append(currentGenomeFile,_parse_readFilepath(file,c));
			//appendValue(genomeFileNames,currentGenomeFile);
			if(options._debugLevel >=2)
				::std::cout <<"Genome file #"<< i <<": " << genomeFileNames[length(genomeFileNames)-1] << ::std::endl;
			++i;
			_parse_skipWhitespace(file, c);
		}
		if(options._debugLevel >=1)
			::std::cout << i-1 << " genome files total." <<::std::endl;
	}
	else		//if file starts with a fasta header --> regular one-genome-file input
		appendValue(genomeFileNames,filename);
	file.close();
	return 0;

}

template <typename TReads, typename TReadNames, typename TReadQualities,typename TSkip,typename TOptions>
int
readQualityValues(TReads & reads, TReadNames &, TReadQualities & readQualities, char const * qualityFilename, TSkip &skip, TOptions & )
{
	::std::ifstream file;
	file.open(qualityFilename,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_QUALITY_FAILED;
	//	return 1;

	unsigned rSeq = 0;
	unsigned matchCount = length(readQualities);
	resize(readQualities, length(reads));
//	std::cout << length(readQualities) << "lenreadqual\n";
	CharString temp_str;
	//typename Iterator<TReads,Rooted>::Type rIt = begin(reads,Rooted());
	char c = _streamGet(file);
	if(c == '@')
	{
		// fastq format
		//vorsicht! need to make sure that identifiers match!!! for gff
		cerr << ".qual file format expected.\n";
		return CALLSNPS_QUALITY_FAILED; //TODO!
	}
	else
	{ //qual format: AG_120_NA_4_1_165_753   14.8    IDI-*II&65'$I&$#+&10()+"*@,?$--)
		while(!_streamEOF(file))
		{
			temp_str = _parse_readIdentifier(file,c);
			_parse_skipWhitespace(file, c); 
			temp_str = _parse_readIdentifier(file,c);
			_parse_skipWhitespace(file, c); 
			temp_str = _parse_readFilepath(file,c);
			if(!skip[rSeq])
			{
				readQualities[matchCount] = prefix(temp_str,length(reads[matchCount]));
				++matchCount;
			}
			++rSeq;
			_parse_skipWhitespace(file, c); 
		}
	}
	file.close();
	if(matchCount == length(reads)) return 0;
	else return CALLSNPS_QUALITY_FAILED;
}


template<typename TFile, typename TChar, typename TWord>
void
_parse_readWord(TFile & file, TWord & str, TChar& c)
{
	// Read word
	clear(str);
	append(str,c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(str, c);
	}
}




/////////////////////////////////////////////////////////////
// read GFF input file containing mapped reads
template <
	typename TMatches,
	typename TGenomeIdMap,
	typename TGenomes,
	typename TReads,
	typename TReadNames,
	typename TReadQualities,
	typename TOptions
>
int readMatchesFromGFF(
	TMatches &matches,							// forward/reverse matches
	TGenomeIdMap const &gIdStringToIdNumMap,				// Read names (read from Fasta file, currently unused)
	TGenomes const &genomes,					// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	TReads &reads,						// Read sequences
	TReadNames &,					// Read names (read from Fasta file, currently unused)
	TReadQualities &readQualities,
	char const * readFName,					// read name (e.g. "reads.fa")
	TOptions &options)
{

	::std::ifstream file;
	file.open(readFName,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GFF_FAILED;
	// count lines and resize matches
	char c = _streamGet(file);
	unsigned count = 0;
	while (!_streamEOF(file))	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		c = _streamGet(file);
		_parse_skipLine(file, c);
		++count;
	}
	file.close();
	file.open(readFName,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GFF_FAILED;
	
	resize(matches,count);
	resize(reads,count);
	resize(readQualities,count);
	//resize(readNames,count);
	
	if(options._debugLevel > 1)::std::cout << count << " lines\n";
	typename Iterator<TMatches,Standard>::Type mIt = begin(matches,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
	typename Iterator<TReads,Standard>::Type rIt = begin(reads,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
	//typename Iterator<TReadNames,Standard>::Type rnIt = begin(readNames,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
	
	unsigned rSeq = 0, pos;
	bool qualityFound, readFound;
	typename TGenomeIdMap::const_iterator it;
	Dna5String gInf;


	CharString readName, temp_str;
	c = _streamGet(file);
	while (!_streamEOF(file))	//if file does not start with a fasta header --> list of multiple reference genome files
	{
//X       run_razers      read            100919085       100919120       2       +       .       ID=s_3_1_3;unique=1;mutations=34A;quality=I)IEIIII-7IA>IIIIII07,-%I>)&#029.2-.
		qualityFound = false;
		readFound = false;
		clear(temp_str);

		(*mIt).rseqNo = rSeq;
		_parse_skipWhitespace(file, c);
	
		clear(temp_str);
		_parse_readIdentifier(file,temp_str,c); //genomeID
		it = gIdStringToIdNumMap.find(temp_str);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		if(it != gIdStringToIdNumMap.end()) (*mIt).gseqNo = it->second;

		_parse_skipWhitespace(file, c);
		temp_str = _parse_readIdentifier(file,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";

		_parse_skipWhitespace(file, c);
		temp_str = _parse_readIdentifier(file,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";

		_parse_skipWhitespace(file, c);
		(*mIt).gBegin = _parse_readNumber(file,c) - options.positionFormat;
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).gBegin << "\t";
		_parse_skipWhitespace(file, c);
		(*mIt).gEnd  = _parse_readNumber(file,c);
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).gEnd << "\t";
		
		_parse_skipWhitespace(file, c);
		(*mIt).mScore = _parse_readNumber(file,c);
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).mScore << "\t";
		_parse_skipWhitespace(file, c);
		if (c=='+')
			(*mIt).orientation = 'F';
		else
			(*mIt).orientation = 'R';
		c = _streamGet(file);
		_parse_skipWhitespace(file, c);
		c = _streamGet(file);
		_parse_skipWhitespace(file, c);
	
                int rLen = (*mIt).gEnd - (*mIt).gBegin;	

		// find ID, check if read with same id has already been stored (only possible when suboptimal/non-unique matches
		// are considered) --> qualityFound = true; readFound = true;
		//temp_str = _parse_readWord(file,c);
		_parse_readWord(file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\n";
		if(temp_str!="ID") ::std::cout << "first feature field should be 'ID'"<<::std::endl;
		c = _streamGet(file);
		clear(temp_str);
		_parse_readIdentifier(file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << "myID = "<<temp_str << "\n";
		//if(rSeq > 0 && temp_str == readNames[rSeq-1])
		//	qualityFound = true;
		//else
		//{
	//		readNames[rSeq] = temp_str;
			gInf = infix(genomes[(*mIt).gseqNo], (*mIt).gBegin, (*mIt).gEnd);
			if(options._debugLevel > 1) cout << gInf << "\n";
			if ((*mIt).orientation == 'R')
				reverseComplementInPlace(gInf);
			reads[rSeq] = gInf;
			//			resize(reads[rSeq],rLen);
			//for(int i = 0; i < rLen; ++i) //TODO vorsicht! funktioniert so nicht fuer indels bzw rLen muss woanders herkommen
			      //     reads[rSeq][i] = (unsigned char)gInf[i];
			
			      //   cout << (Dna5String)reads[rSeq] << "\n",		
		//}
		(*mIt).editDist = 0;
		//::std::cout << "hier\n";

		_parse_skipWhitespace(file,c);
		while(!(c == '\n' || (c == '\r' && _streamPeek(file) != '\n'))) // while in same line
		{
		//::std::cout << "hier\n";
			while(c != ';')
			{
				if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
					break;
				c = _streamGet(file);
			}
			if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
				break;
			c = _streamGet(file);
			_parse_readWord(file,temp_str,c);
			//c = _streamGet(file);
			if(options._debugLevel > 1) 
				::std::cout << temp_str << " in schleife\n";
			if(!options.qualityFile && temp_str=="quality")
			{
				CharString qual;
				//::std::cout << rLen<<"rLen\n";
				qualityFound = true;
				for(int i = 0; i < rLen; ++i) //TODO vorsicht! funktioniert so nicht fuer indels bzw rLen muss woanders herkommen
				{
					c = _streamGet(file);
					//cout << c << "\n";
//					int q = (ordValue(c)>64) ? 31 : (ordValue(c)-33);
					//cout << rSeq<<"\n" << length(reads) << "\n" << ordValue((Dna5)reads[rSeq][i]) << " \n";
			//		reads[rSeq][i] = (unsigned char)(ordValue((Dna5)reads[rSeq][i]) |(q<<3));
				        //cout << q<<"\n";
					append(qual, c);
				}
		//		::std::cout << "hier2\n";
				readQualities[rSeq] = qual;
				qualityFound = true;
//				if(options._debugLevel > 1) 
//					::std::cout << "qualityString=" << qual << "\n";
//				if(length(readQualities) != rSeq) ::std::cout << "hier ist was faul.\n";
			}
			else
			{
				 if(temp_str=="mutations")
				{
					readFound = true;
					while(c==',' || c=='=')
					{
						c = _streamGet(file);
						pos = _parse_readNumber(file,c);
//						int q = (int)(reads[rSeq][pos-1]) >>3;
						reads[rSeq][pos-1]=c;
						//::std::cout<<readNames[rSeq]<<"=" << reads[rSeq] << " with edit="<<(*mIt).editDist<<" ("<<(Dna)c << " at pos " << pos<<")\n";
						c = _streamGet(file);
						++(*mIt).editDist;
					}
				}
			}
			if(qualityFound) {_parse_skipLine(file,c); break;}
			
		}
//		if ((rSeq%10000)==0)cout <<rSeq<<".."<<std::flush;
		//(*mIt)=m;
		if((*mIt).mScore > 0)
		{++mIt;
                ++rSeq;
		++rIt;
		//++rnIt;
		}
		if(options._debugLevel > 1) 
			::std::cout<<reads[rSeq-1]<<" with edit="<<(*mIt).editDist<<" at position "<< matches[rSeq-1].gBegin<<"\n";
		_parse_skipWhitespace(file, c);
//_parse_skipLine(file, c);
	}
	file.close();
	resize(matches,rSeq);
	resize(reads,rSeq);
	//resize(readNames,rSeq);
	resize(readQualities,rSeq);
	::std::cout << rSeq << " matches with mapping quality > 0.\n";
	if(options._debugLevel > 0)
		::std::cout << "Parsed "<<length(matches)<<" matches of "<<rSeq<<" reads." << ::std::endl;
	return 0;
}


/////////////////////////////////////////////////////////////
// read GFF input file containing mapped reads
template <
	typename TMatches,
	typename TGenomeIdMap,
	typename TGenomes,
	typename TReads,
	typename TReadNames,
	typename TReadQualities,
	typename TOptions
>
int readMatchesFromMaq(
	TMatches &matches,							// forward/reverse matches
	TGenomeIdMap const &gIdStringToIdNumMap,				// Read names (read from Fasta file, currently unused)
	TGenomes const &,					// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	TReads &reads,						// Read sequences
	TReadNames &,					// Read names (read from Fasta file, currently unused)
	TReadQualities & readQualities,
	char const * readFName,					// read name (e.g. "reads.fa")
	TOptions &options)
{


	::std::ifstream file;
	file.open(readFName,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GFF_FAILED;
	// count lines and resize matches
	char c = _streamGet(file);
	unsigned count = 0;
	while (!_streamEOF(file))	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		c = _streamGet(file);
		_parse_skipLine(file, c);
		++count;
	}
	file.close();
	file.open(readFName,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GFF_FAILED;
		
	resize(matches,count);
	resize(reads,count);
	resize(readQualities,count);
	//resize(readNames,count);
	
	if(options._debugLevel > 1)::std::cout << count << "lines\n";
	typename Iterator<TMatches,Standard>::Type mIt = begin(matches,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
	typename Iterator<TReads,Standard>::Type rIt = begin(reads,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
//	typename Iterator<TReadNames,Standard>::Type rnIt = begin(readNames,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
	
	unsigned rSeq = 0;
	typename TGenomeIdMap::const_iterator it;
	Dna5String gInf;
	CharString readName, temp_str;//, qual;


	c = _streamGet(file);
	while (!_streamEOF(file))	//if file does not start with a fasta header --> list of multiple reference genome files
	{
//s_4_100_67233   X       4       +       0       0       0       0       0       2       13      0       60      36      aCCcgAacCctaaccctaaccctaaccctaacccta    7<E"+I-5I#%6+3&2#,4&/3'+')2+#&)&"(#$
		clear(temp_str);

		(*mIt).rseqNo = rSeq;
		_parse_skipWhitespace(file, c);
	
		clear(temp_str);
		_parse_readIdentifier(file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << "myID = "<<temp_str << "\n";
	//	readNames[rSeq] = temp_str;
		
		_parse_skipWhitespace(file, c);
		clear(temp_str);
		_parse_readIdentifier(file,temp_str,c); //genomeID
		it = gIdStringToIdNumMap.find(temp_str);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		if(it != gIdStringToIdNumMap.end()) (*mIt).gseqNo = it->second;
		
		_parse_skipWhitespace(file, c);
		(*mIt).gBegin = _parse_readNumber(file,c) - options.positionFormat;
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).gBegin << "\t";

		
		_parse_skipWhitespace(file, c);
		if (c=='+')
			(*mIt).orientation = 'F';
		else
			(*mIt).orientation = 'R';
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).orientation << "\t";

		//if(rSeq>10) options._debugLevel = 0;

		c = _streamGet(file);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file,c);

		_parse_skipWhitespace(file, c);
		_parse_readNumber(file,c);

		
		_parse_skipWhitespace(file, c);
		(*mIt).mScore = _parse_readNumber(file,c);
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).mScore << "\t";


		_parse_skipWhitespace(file, c);
		_parse_readNumber(file,c);

		_parse_skipWhitespace(file, c);
		_parse_readNumber(file,c);


		_parse_skipWhitespace(file, c);
		(*mIt).editDist = _parse_readNumber(file,c);
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).editDist << "\t";


		_parse_skipWhitespace(file, c);
		_parse_readNumber(file,c);

		_parse_skipWhitespace(file, c);
		_parse_readNumber(file,c);

		_parse_skipWhitespace(file, c);
		_parse_readNumber(file,c);

		_parse_skipWhitespace(file, c);
		unsigned rLen = _parse_readNumber(file,c);
		(*mIt).gEnd = (*mIt).gBegin + rLen;	
		
		
		_parse_skipWhitespace(file, c);
		clear(temp_str);
		temp_str = _parse_readIdentifier(file,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		if((*mIt).orientation == 'R')
		{
		//	std::cout << "hierrlen="<<rLen<<"\n";
			String<Dna5> tempR = temp_str;
			reverseComplementInPlace(tempR);
			resize(reads[rSeq],rLen);
			for(unsigned k=0; k < rLen; ++k)
				reads[rSeq][k] = (Dna5/*Q*/)tempR[k];
		}
		else
			reads[rSeq] = (String<Dna5>) temp_str;
		
		_parse_skipWhitespace(file, c);
		clear(temp_str);
		temp_str = _parse_readFilepath(file,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\n";
		if((*mIt).orientation == 'R')
		{
			resize(readQualities[rSeq],rLen);
			for(unsigned k=0; k < rLen; ++k)
				readQualities[rSeq][k] = temp_str[rLen-k-1];
		}
		else
			readQualities[rSeq] = temp_str;

		if(options._debugLevel>1)::std::cout << (String<Dna>) reads[rSeq]<< "\t" << readQualities[rSeq] << "\n";
		
		if((*mIt).mScore > 0)
		{
			++mIt;
                	++rSeq;
			++rIt;
//			++rnIt;
		}
	//	if(options._debugLevel > 1) 
	//		::std::cout<<*rnIt<<"=" << (Dna5String)(*rIt) << " with edit="<<(*mIt).editDist<<"\n";
		_parse_skipWhitespace(file, c);
	}
	file.close();
	resize(matches,rSeq);
	resize(reads,rSeq);
	resize(readQualities,rSeq);
//	resize(readNames,rSeq);
	if(options._debugLevel > 0) ::std::cout << rSeq << " matches with mapping quality > 0.\n";
	if(options._debugLevel > 0)
		::std::cout << "Parsed "<<length(matches)<<" matches of "<<rSeq<<" reads." << ::std::endl;
	return 0;
}


/////////////////////////////////////////////////////////////
// read Eland input file containing mapped reads
template <
	typename TMatches,
	typename TGenomeIdMap,
	typename TGenomes,
	typename TReads,
	typename TReadNames,
	typename TSkip,
	typename TOptions
>
int readMatchesFromEland(
	TMatches &matches,							// forward/reverse matches
	TGenomeIdMap const &gIdStringToIdNumMap,				// Read names (read from Fasta file, currently unused)
	TGenomes const &,					// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	TReads &reads,						// Read sequences
	TReadNames &,					// Read names (read from Fasta file, currently unused)
	char const * readFName,					// read name (e.g. "reads.fa")
	TSkip &skip,
	TOptions &options)
{

	unsigned matchCount = length(matches);	

	::std::ifstream file;
	file.open(readFName,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GFF_FAILED;
	// count lines and resize reads
	char c = _streamGet(file);
	unsigned count = 0;
	while (!_streamEOF(file))	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		c = _streamGet(file);
		_parse_skipLine(file, c);
		++count;
	}
	file.close();
	file.open(readFName,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GFF_FAILED;
		
//	resize(matches,count);
	resize(reads,matchCount+count);
	//resize(readNames,count);
	reserve(matches,matchCount+count);
//	reserve(reads,count);
//	reserve(readNames,count);
	fill(skip,count,true);
	
	if(options._debugLevel > 1) ::std::cout << count << " lines\n";
	
	unsigned rSeq = 0;
	typename TGenomeIdMap::const_iterator it;
	Dna5String gInf;
	CharString readName, temp_str;//, qual;
	typename Value<TMatches>::Type m = { 0, 0, 0, 0, 0, 0, 0, 0 };	// extra value for sum of mismatching qualities
	
	bool first = true;
	unsigned rLen = 0;
	int idOffset = 0;
	c = _streamGet(file);
	if(c=='>') idOffset = 1;
	while (!_streamEOF(file))	//if file does not start with a fasta header --> list of multiple reference genome files
	{
//>AG_120_NA_4_1_1032_922 GTTCCTGCACATCTTTTTCATGCCTCCCATTC        U1      0       1       0       Homo_sapiens.NCBI36.49.dna.chromosome.X.fa   47154562 F       DD      24A	
		clear(temp_str);
	//	clear(readName);
		m.rseqNo = matchCount;
	
		/*temp_str = */ _parse_readWordUntilWhitespace(file,c);
		//temp_str = _parse_readIdentifier(file,c); //readID
		//readNames[matchCount] = suffix(temp_str,idOffset);
//		clear(temp_str);

		_parse_skipWhitespace(file, c);
		temp_str = _parse_readIdentifier(file,c);
		if(first){
			rLen = length(temp_str);
//			cout << "rLem="<<rLen<<"\n";
		first =false;}
//		reads[rSeq] = temp_str;
		_parse_skipWhitespace(file, c);
		//c = _streamGet(file);	//U
		unsigned tempRLen = length(temp_str);
		if(c=='U' && tempRLen==rLen)
		{
			//appendValue(readNames,readName);
			if(tempRLen>rLen)
			{
				CharString temptempstr = temp_str;	
				if(temp_str[0]=='.')
					temp_str = infix(temptempstr,1,1+rLen);
			//	 cout << temp_str <<"<-\n";
			}
			//appendValue(reads,temp_str);
			reads[matchCount] = temp_str;

			c = _streamGet(file);	//eg 2
			//std::cout << "c="<<(*mIt).editDist<<"\n";
			m.editDist = _parse_readNumber(file, c);
			//std::cout << "edit="<<(*mIt).editDist<<"\n";
			_parse_skipWhitespace(file, c);
			_parse_readNumber(file, c);
			_parse_skipWhitespace(file, c);
			_parse_readNumber(file, c);
			_parse_skipWhitespace(file, c);
			_parse_readNumber(file, c);
			_parse_skipWhitespace(file, c);
			temp_str = _parse_readIdentifier(file,c); //chrID
				
			// remove the file name of chromosome
			::std::string temp_str_c = toCString(temp_str);
			size_t lastPos1 = (temp_str_c).find_last_of('.');
			size_t lastPos2 = (temp_str_c).find_last_of('.',lastPos1-2);
			if (lastPos2 == temp_str_c.npos) lastPos2 = 0;
			//::std::strin = readFName.substr(lastPos2,lastPos);
			temp_str = infix(temp_str,lastPos2+1,lastPos1);
			it = gIdStringToIdNumMap.find(temp_str);
	//		if(options._debugLevel > 1) 
	//			::std::cout << temp_str <<"<-chr\n";
			if(it != gIdStringToIdNumMap.end())
			{
//			
//				if(options._debugLevel > 1)
//					std::cout << readNames[matchCount] << "<-readId\n";
				m.gseqNo = it->second;
				_parse_skipWhitespace(file, c);
				m.gBegin = _parse_readNumber(file, c) - options.positionFormat;
//				if(m.gBegin==215293) std::cout <<"daaaaaaaaaaaaaaaaaa\n";//extraV = true;
				m.gEnd = m.gBegin + rLen;
			//	if(rLen>32)
			//	{
			//		::std::cerr <<"eheeeeeeeee\n";
			//		::std::cerr << infix(genomes[0],m.gBegin,m.gEnd) << "<-g";
			//	}
				_parse_skipWhitespace(file, c);
				m.orientation = c;
				m.mScore = 0;  //no mapping quality given
				//matches[matchCount] = m;
				appendValue(matches,m);
				++matchCount;
				skip[rSeq] = false;
				if(options._debugLevel > 1)
				{
					//std::cout << readNames[matchCount] << '\t';
					std::cout << reads[matchCount] << '\t';
					std::cout << m.gBegin << '\t' << m.editDist <<'\t' << m.orientation <<std::endl;
					
				}
			}//
		}
		++rSeq;
		if(options._debugLevel==1 && (rSeq%1000000)==0) cout << rSeq<<".."<<std::flush;
		_parse_skipLine(file,c);
	}
	file.close();
	if(options._debugLevel==1) cout << endl;
//	reserve(matches, matchCount, Exact());
	resize(reads,matchCount);
//	resize(readNames,matchCount);
	if(options._debugLevel > 0)
		::std::cout << "Parsed "<<length(matches)<<" matches of "<<rSeq<<" reads." << ::std::endl;
	return 0;
}

/*
/////////////////////////////////////////////////////////////
// read GFF input file containing mapped reads
template <
	typename TMatches,
	typename TGenomeIdMap,
	typename TGenomes,
	typename TReads,
	typename TReadNames,
	typename TOptions
>
int readMatchesFromGFF(
	TMatches &matches,							// forward/reverse matches
	TGenomeIdMap const &gIdStringToIdNumMap,				// Read names (read from Fasta file, currently unused)
	TGenomes const &genomes,					// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	TReads &reads,						// Read sequences
	TReadNames &readNames,					// Read names (read from Fasta file, currently unused)
	char const * readFName,					// read name (e.g. "reads.fa")
	TOptions &options)
{

	::std::ifstream file;
	file.open(readFName,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GFF_FAILED;
	// count lines and resize matches
	char c = _streamGet(file);
	unsigned count = 0;
	while (!_streamEOF(file))	//if file does not start with a fasta header --> list of multiple reference genome files
	{
		c = _streamGet(file);
		_parse_skipLine(file, c);
		++count;
	}
	file.close();
	file.open(readFName,::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
		return CALLSNPS_GFF_FAILED;
		
	resize(matches,count);
	resize(reads,count);
	resize(readNames,count);
	
	::std::cout << count << "many matches\n";
	typename Iterator<TMatches,Standard>::Type mIt = begin(matches,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
	typename Iterator<TReads,Standard>::Type rIt = begin(reads,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
	typename Iterator<TReadNames,Standard>::Type rnIt = begin(readNames,Standard());// = { 0, 0, 0, 0, 0, 0, 0};
	
	unsigned rSeq = 0, pos;
	bool qualityFound, readFound;
	typename TGenomeIdMap::const_iterator it;
	Dna5String gInf;


	c = _streamGet(file);
	while (!_streamEOF(file))	//if file does not start with a fasta header --> list of multiple reference genome files
	{
//X       run_razers      read            100919085       100919120       2       +       .       ID=s_3_1_3;unique=1;mutations=34A;quality=I)IEIIII-7IA>IIIIII07,-%I>)&#029.2-.
	CharString readName, temp_str;//, qual;
		qualityFound = false;
		readFound = false;
		clear(temp_str);

		(*mIt).rseqNo = rSeq;
		_parse_skipWhitespace(file, c);
	
		clear(temp_str);
		_parse_readIdentifier(file,temp_str,c); //genomeID
		it = gIdStringToIdNumMap.find(temp_str);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		if(it != gIdStringToIdNumMap.end()) (*mIt).gseqNo = it->second;

		_parse_skipWhitespace(file, c);
		temp_str = _parse_readIdentifier(file,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";

		_parse_skipWhitespace(file, c);
		temp_str = _parse_readIdentifier(file,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";

		_parse_skipWhitespace(file, c);
		(*mIt).gBegin = _parse_readNumber(file,c) - options.positionFormat;
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).gBegin << "\t";
		_parse_skipWhitespace(file, c);
		(*mIt).gEnd  = _parse_readNumber(file,c);
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).gEnd << "\t";
		
		_parse_skipWhitespace(file, c);
		(*mIt).mScore = _parse_readNumber(file,c);
		if(options._debugLevel > 1) 
			::std::cout << (*mIt).mScore << "\t";
		_parse_skipWhitespace(file, c);
		if (c=='+')
			(*mIt).orientation = 'F';
		else
			(*mIt).orientation = 'R';
		c = _streamGet(file);
		_parse_skipWhitespace(file, c);
		c = _streamGet(file);
		_parse_skipWhitespace(file, c);
	
                int rLen = (*mIt).gEnd - (*mIt).gBegin;	

		// find ID, check if read with same id has already been stored (only possible when suboptimal/non-unique matches
		// are considered) --> qualityFound = true; readFound = true;
		//temp_str = _parse_readWord(file,c);
		_parse_readWord(file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\n";
		if(temp_str!="ID") ::std::cout << "first feature field should be 'ID'"<<::std::endl;
		c = _streamGet(file);
		clear(temp_str);
		_parse_readIdentifier(file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << "myID = "<<temp_str << "\n";
		//if(rSeq > 0 && temp_str == readNames[rSeq-1])
		//	qualityFound = true;
		//else
		//{
			value(rnIt) = temp_str;
			gInf = infix(genomes[(*mIt).gseqNo], (*mIt).gBegin, (*mIt).gEnd);
			if(options._debugLevel > 1) cout << gInf << "\n";
			if ((*mIt).orientation == 'R')
				reverseComplementInPlace(gInf);
			reads[rSeq] = gInf;
			//			resize(reads[rSeq],rLen);
			//for(int i = 0; i < rLen; ++i) //TODO vorsicht! funktioniert so nicht fuer indels bzw rLen muss woanders herkommen
			      //     reads[rSeq][i] = (unsigned char)gInf[i];
			
			      //   cout << (Dna5String)reads[rSeq] << "\n",		
		//}
		(*mIt).editDist = 0;
		//::std::cout << "hier\n";

		_parse_skipWhitespace(file,c);
		while(!(c == '\n' || (c == '\r' && _streamPeek(file) != '\n'))) // while in same line
		{
		//::std::cout << "hier\n";
			while(c != ';')
			{
				if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
					break;
				c = _streamGet(file);
			}
			if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
				break;
			c = _streamGet(file);
			_parse_readWord(file,temp_str,c);
			//c = _streamGet(file);
			if(options._debugLevel > 1) 
				::std::cout << temp_str << " in schleife\n";
			if(temp_str=="quality")
			{
//		//		clear(qual);
				//::std::cout << rLen<<"rLen\n";
				qualityFound = true;
				for(int i = 0; i < rLen; ++i) //TODO vorsicht! funktioniert so nicht fuer indels bzw rLen muss woanders herkommen
				{
					c = _streamGet(file);
					//cout << c << "\n";
					int q = (ordValue(c)>64) ? 31 : (ordValue(c)-33);
					//cout << rSeq<<"\n" << length(reads) << "\n" << ordValue((Dna5)reads[rSeq][i]) << " \n";
					reads[rSeq][i] = (unsigned char)(ordValue((Dna5)reads[rSeq][i]) |(q<<3));
				        //cout << q<<"\n";
//					append(qual, c);
				}
		//		::std::cout << "hier2\n";
//				appendValue(readQualities,qual);
				      qualityFound = true;
//				if(options._debugLevel > 1) 
//					::std::cout << "qualityString=" << qual << "\n";
//				if(length(readQualities) != rSeq) ::std::cout << "hier ist was faul.\n";
			}
			else
			{
				 if(temp_str=="mutations")
				{
					readFound = true;
					while(c==',' || c=='=')
					{
						c = _streamGet(file);
						pos = _parse_readNumber(file,c);
						int q = (int)(reads[rSeq][pos-1]) >>3;
						reads[rSeq][pos-1]=(unsigned char)(q<<3|ordValue((Dna)c));
						//::std::cout<<readNames[rSeq]<<"=" << reads[rSeq] << " with edit="<<m.editDist<<" ("<<(Dna)c << " at pos " << pos<<")\n";
						c = _streamGet(file);
						++(*mIt).editDist;
					}
				}
			}
			if(qualityFound) {_parse_skipLine(file,c); break;}
			
		}
		if (options._debugLevel > 0 && (rSeq%10000)==0)cout <<rSeq<<".."<<std::flush;
		//(*mIt)=m;
		if((*mIt).mScore> 0 )
		{++mIt;
                ++rSeq;
		++rIt;
		++rnIt;
		}if(options._debugLevel > 1) 
			::std::cout<<*rnIt<<"=" << (Dna5String)(*rIt) << " with edit="<<(*mIt).editDist<<"\n";
		_parse_skipWhitespace(file, c);
//_parse_skipLine(file, c);
	}
	file.close();
	if(options._debugLevel > 0)
		::std::cout << "Parsed "<<length(matches)<<" matches of "<<rSeq<<" reads." << ::std::endl;
	return 0;
}

*/


template <typename TOptions>
void printHelp(int, const char *[], TOptions &options, bool longHelp = false)
{

	cerr << "Usage: callSNPs [OPTION]... <GENOME FILE> <MAPPED READ FILE>" << endl;
	if (longHelp) {
		cerr << endl << "Options:" << endl;
		cerr << "  -if, --input-format NUM      \t" << "input format:" << endl;
		cerr << "                               \t" << "0 = razers gff (default)" << endl;
		cerr << "                               \t" << "1 = eland" << endl;
		cerr << "                               \t" << "2 = maq" << endl;
		cerr << "  -q,  --qual-file FILE        \t" << "file containing read qualities (obligatory if -i 1)" << endl;
		cerr << "  -o,  --output FILE           \t" << "change output filename (default <READ FILE>.snp)" << endl;
//		cerr << "  -of, --output-format NUM     \t" << "output format:" << endl;
//		cerr << "                               \t" << "0 = quality strings (default)" << endl;
//		cerr << "                               \t" << "1 = integer quality sums + counts" << endl;
		cerr << "  -t,  --tab-file FILE         \t" << "produce tab file" << endl;
/*		cerr << "  -m,  --method NUM            \t" << "set method used for SNP calling" << endl;
		cerr << "                               \t" << "0 = threshold method" << endl;
		cerr << "                               \t" << "1 = binomial" << endl;
		cerr << "                               \t" << "2 = maq (requires mapping qualities in GFF file)" << endl;
		cerr << "                               \t" << "(default = "<<options.method << ")" << endl;*/
		cerr << "  -mp, --max-pile NUM          \t" << "maximal number of reads allowed to pile up at the same genome position ("<<options.maxPile<<")" << endl;
		cerr << "  -mc, --min-coverage NUM      \t" << "minimal number of reads covering a candidate position ("<< options.minCoverage<<")" << endl;
		cerr << "  -oa, --orientation-aware     \t" << "distinguish between forward and reverse matches (off)" << endl;
		cerr << endl;
		cerr << "SNP calling options: " << endl;
		cerr << "  -mm, --min-mutations NUM     \t" << "minimal number of observed mutations for mutation to be called ("<<options.minMutT<<")" << endl;
		cerr << "  -pt, --perc-threshold NUM    \t" << "minimal percentage of mutational base for mutation to be called (" << options.percentageT << ")" << endl;
		cerr << "  -mq, --min-Quality NUM       \t" << "minimal average quality of mutational base for mutation to be called ("<<options.avgQualT <<")" << endl;
		cerr << endl;
		cerr << "Other options: " << endl;
		cerr << "  -p,  --position-file FILE    \t" << "file containing positions to inspect (assumes only 1 input chromosome)" << endl;
		cerr << "  -op, --output-position FILE  \t" << "change output filename for position inspection (default <READ FILE>.stats)" << endl;
		cerr << "  -ed, --error-distr FILE      \t" << "write error distribution to FILE" << endl;
		cerr << "  -v,  --verbose               \t" << "verbose mode" << endl;
		cerr << "  -vv, --very-verbose          \t" << "very verbose mode" << endl;
		cerr << "  -h,  --help                  \t" << "print this help" << endl;
	}
	else {
		cerr << "Try 'callSNPs --help' for more information." << endl <<endl;
	}
}



int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line

	SNPCallingOptions<>		options;
	options.genomeNaming = 0;
	options.readNaming = 0;

	unsigned			fnameCount = 0;
	const char			*genomeFName = "";
	String<CharString> 		readFNames;
	String<CharString> 		qualityFNames;


	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse options
			if (strcmp(argv[arg], "-m") == 0 || strcmp(argv[arg], "--method") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.method;
					if (!istr.fail())
					{
						if (options.method > 1)
							cerr << "Invalid method option." << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mq") == 0 || strcmp(argv[arg], "--min-quality") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.avgQualT;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mm") == 0 || strcmp(argv[arg], "--min-mutations") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.minMutT;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-pt") == 0 || strcmp(argv[arg], "--perc-threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.percentageT;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mp") == 0 || strcmp(argv[arg], "--max-pile") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.maxPile;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-mc") == 0 || strcmp(argv[arg], "--min-coverage") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.minCoverage;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-oa") == 0 || strcmp(argv[arg], "--orientation-aware") == 0) {
				options.orientationAware = true;
				continue;
			}
			if (strcmp(argv[arg], "-if") == 0 || strcmp(argv[arg], "--input-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.inputFormat;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-of") == 0 || strcmp(argv[arg], "--output-format") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.outputFormat;
					if (!istr.fail())
						continue;
				}
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-q") == 0 || strcmp(argv[arg], "--qual-file") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				options.qualityFile = true;
				++arg;
				String<char> tempStr = argv[arg];
				if(argv[arg][0]=='[')
				{
					appendValue(qualityFNames,suffix(tempStr,1));
					++arg;
					bool inList = true;
					if(arg < argc && (argv[arg][0] == '-' || tempStr[length(tempStr)-1]==']')) inList = false;
					while(arg < argc && inList)
					{
						tempStr = argv[arg];
						appendValue(qualityFNames,tempStr);
						++arg;
						if(arg < argc && (argv[arg][0] == '-' || tempStr[length(tempStr)-1]==']')) inList = false;
					}
					if(tempStr[length(tempStr)-1] != ']') cerr << "Something wrong with quality file list?" << endl;
					resize(qualityFNames[length(qualityFNames)-1],length(qualityFNames[length(qualityFNames)-1])-1);
					--arg;
				}
				else appendValue(qualityFNames,tempStr);

				continue;
			}
			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.outputSNP = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-p") == 0 || strcmp(argv[arg], "--position-file") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.positionFile = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-t") == 0 || strcmp(argv[arg], "--tab-file") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.tabFile = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-ed") == 0 || strcmp(argv[arg], "--error-distr") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.errorPrbFileName = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-op") == 0 || strcmp(argv[arg], "--output-position") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options, true);
					return 0;
				}
				++arg;
				options.outputPositionAnalysis = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, options, true);
				return 0;
			}
			if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 1);
				continue;
			}
			if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--very-verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 2);
				continue;
			}
			cerr << "Unknown option: " << argv[arg] << endl << endl;
			printHelp(argc, argv, options);
			return 0;
		} else {
			// parse file name
			if (fnameCount == 0)
				genomeFName = argv[arg];
			if (fnameCount == 1)
			{
				String<char> tempStr = argv[arg];
				if(argv[arg][0]=='[')
				{
					appendValue(readFNames,suffix(tempStr,1));
					++arg;
					while(arg < argc && argv[arg][0] != '-')
					{
						tempStr = argv[arg];
						appendValue(readFNames,tempStr);
						++arg;
					}
					if(readFNames[length(readFNames)-1][length(readFNames[length(readFNames)-1])-1] != ']') cerr << "Something wrong with read file list?" << endl;
					resize(readFNames[length(readFNames)-1],length(readFNames[length(readFNames)-1])-1);
					--arg;
				}
				else appendValue(readFNames,tempStr);
			}
			if (fnameCount == 2) {
				cerr << "More than 2 input files specified." <<endl;
				cerr << "If more than 2 mapped read files are to be parsed, use quotation marks directly before first file name and directly after last file name (e.g. \"lane1.gff lane2.gff\")." << endl << endl;
				printHelp(argc, argv, options);
				return 0;
			}
			++fnameCount;
		}
	}
	if (fnameCount != 2) {
		if (argc > 1 && !options.printVersion)
			cerr << "Exactly 2 input files need to be specified." << endl << endl;
		printHelp(argc, argv, options);
		return 0;
	}
	if(options.inputFormat == 1 && (!options.qualityFile || (length(qualityFNames)!=length(readFNames))))
	{
		cerr << "If mapped read file is in Eland format, a .qual or .fastq file containing read qualities needs to be specified." << endl << endl;
		return 0;
	}
	if(*(options.positionFile) == 0 && *(options.outputPositionAnalysis) != 0)
	{
		cerr << "Position analysis output specified, but no position file given." << endl << endl;
		return 0;
	}

//	if (options.printVersion)
//		printVersion();

	int result = detectSNPs(genomeFName, readFNames, qualityFNames, options);
	if (result > 0)
	{
		printHelp(argc, argv, options);
		return 0;
	}
	return result;
}
