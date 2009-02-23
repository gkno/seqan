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

#ifndef SEQAN_HEADER_OUTPUT_FORMAT_H
#define SEQAN_HEADER_OUTPUT_FORMAT_H

#include <iostream>
#include <fstream>
#include <sstream>

#include "razers.h"
#include <seqan/align.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Less-operators ...

	// ... to sort matches and remove duplicates with equal gBegin
	template <typename TReadMatch>
	struct LessGPosRNo : public ::std::binary_function < TReadMatch, TReadMatch, bool >
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

			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};


//////////////////////////////////////////////////////////////////////////////
// Determine error distribution
template <typename TErrDistr, typename TMatches, typename TReads, typename TGenomes, typename TOptions>
inline unsigned
getErrorDistribution(
	TErrDistr &posError, 
	TMatches &matches, 
	TReads &reads, 
	TGenomes &genomes, 
	TOptions &options)
{
	typename Iterator<TMatches, Standard>::Type	it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type	itEnd = end(matches, Standard());

	Dna5String genome;
	unsigned unique = 0;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;

		Dna5String const &read = reads[(*it).rseqNo];
		genome = infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd);
		if ((*it).orientation == 'R')
			reverseComplementInPlace(genome);
		for (unsigned i = 0; i < length(posError) && i < length(read); ++i)
			if ((options.compMask[ordValue(genome[i])] & options.compMask[ordValue(read[i])]) == 0)
				++posError[i]; 
		++unique;
	}
	return unique;
}

template <typename TErrDistr, typename TCount1, typename TCount2, typename TMatches, typename TReads, typename TGenomes, typename TSpec>
inline unsigned
getErrorDistribution(
	TErrDistr &posError,
	TCount1 &insertions,
	TCount2 &deletions,
	TMatches &matches, 
	TReads &reads, 
	TGenomes &genomes, 
	RazerSOptions<TSpec> &options)
{
	typedef Align<String<Dna5>, ArrayGaps> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TIter;

	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;

	typename Iterator<TMatches, Standard>::Type	it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type	itEnd = end(matches, Standard());

	Align<Dna5String, ArrayGaps> align;
	Score<int> scoreType = Score<int>(0, -999, -1001, -1000);	// levenshtein-score (match, mismatch, gapOpen, gapExtend)
	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);

	unsigned unique = 0;
	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;

		assignSource(row(align, 0), reads[(*it).rseqNo]);
		assignSource(row(align, 1), infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd));
		if ((*it).orientation == 'R')
			reverseComplementInPlace(source(row(align, 1)));
		globalAlignment(align, scoreType);
		
		TRow& row0 = row(align, 0);
		TRow& row1 = row(align, 1);
		
		TPosition begin = beginPosition(cols(align));
		TPosition end = endPosition(cols(align));
		
		TIter it0 = iter(row0, begin);
		TIter it1 = iter(row1, begin);
		TIter end0 = iter(row0, end);
		
		unsigned pos = 0;
		for (; it0 != end0 && pos < length(posError); ++it0, ++it1)
		{
			if (isGap(it0))
				++insertions;
			else
			{
				if (isGap(it1))
					++deletions;
				else
					if ((options.compMask[ordValue(getValue(it0))] & options.compMask[ordValue(getValue(it1))]) == 0)
						++posError[pos];
				++pos;
			}
		}
		++unique;
	}
	return unique;
}


//////////////////////////////////////////////////////////////////////////////
// Dump an alignment
template <typename TFile, typename TSource, typename TSpec>
inline void
dumpAlignment(TFile & target, Align<TSource, TSpec> const & source)
{
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

template<typename TMatches, typename TCounts, typename TOptions>
void
countCoocurrences(TMatches & matches, TCounts & cooc, TOptions & options)
{
	clear(cooc);
	int maxSeedErrors = (int)(options.errorRate * options.artSeedLength) + 1;
	fill(cooc,maxSeedErrors+1,0);
	for (int i = 0; i < maxSeedErrors+1; ++i)
		cooc[i] = 1;
	
	int count = 0;
	unsigned readNo = -1;
	int preEditDist = -1;
	typename Iterator<TMatches>::Type it = begin(matches,Standard());
	typename Iterator<TMatches>::Type itEnd = end(matches,Standard());
	
	for(; it != itEnd; ++it)
	{
		if ((*it).rseqNo == readNo)
		{
			if(preEditDist > 1) continue;// || dist > options.errorRate * maxReadLength + 1) continue;
			int dist = (*it).seedEditDist - preEditDist;
			if(dist > maxSeedErrors) continue;
			if(dist < 0) ++cooc[0];
			else ++cooc[dist];
		}
		else
		{
			readNo = (*it).rseqNo;
			preEditDist = (*it).seedEditDist;
			if(preEditDist <= 1) ++count;
		}
	}
/*	fprintf(stderr, "[abs_mapping_count] %d, %d, %d, %d\n", cooc[0], cooc[1], cooc[2], cooc[3]);
	printf("n=%i\n",count);*/
	for (unsigned i = 0; i < length(cooc); ++i)
	{
		cooc[i] = (int)(-4.343 * log((double)cooc[i]/count) );
		if (cooc[i] < 0) cooc[i] = 0;
	}
	if(options._debugLevel > 0)
	{
		::std::cerr << "[mapping_count] ";
		for(unsigned j = 0; j < length(cooc); ++j)
			::std::cerr << cooc[j] << " ";
		::std::cerr << ::std::endl;
	}

}


//////////////////////////////////////////////////////////////////////////////
// Output matches
template <
	typename TMatches,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames,
	typename TCounts,
	typename TSpec
>
void dumpMatches(
	TMatches &matches,							// forward/reverse matches
	TGenomeNames const &genomeIDs,				// Read names (read from Fasta file, currently unused)
	StringSet<CharString> &genomeFileNameList,	// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > &gnoToFileMap, //map to retrieve genome filename and sequence number within that file
	TReads const &reads,
	TCounts & stats,						// Match statistics (possibly empty)
	TReadNames const &readIDs,					// Read names (read from Fasta file, currently unused)
	::std::string readFName,					// read name (e.g. "reads.fa")
	::std::string errorPrbFileName,
	RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatches>::Type		TMatch;
	typedef typename Value<TReads>::Type		TRead;
	typedef typename Value<TGenomeSet>::Type	TGenome;
	typedef typename TMatch::TGPos				TGPos;

	if (options.outputFormat == 2)
	{
		options.maxHits = 1;		// Eland outputs at most one match
		options.sortOrder = 0;		// read numbers are increasing
		options.positionFormat = 1;	// bases in file are numbered starting at 1
		options.dumpAlignment = options.hammingOnly;
	}
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if (options.maqMapping) options.outputFormat = 3;
#endif
	if (options.outputFormat == 3)
	{
		options.sortOrder = 0;		// read numbers are increasing
		options.positionFormat = 1;	// bases in file are numbered starting at 1
		options.dumpAlignment = options.hammingOnly;
	}
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping) 
		options.dumpAlignment = false;
#endif

	// error profile
	unsigned maxReadLength = 0;
	for (unsigned i = 0; i < length(reads); ++i)
		if (maxReadLength < length(reads[i]))
			maxReadLength = length(reads[i]);

	SEQAN_PROTIMESTART(dump_time);

	// load Genome sequences for alignment dumps
	TGenomeSet genomes;
	if (options.dumpAlignment || !errorPrbFileName.empty())
		if (!loadGenomes(genomes, genomeFileNameList)) {
			::std::cerr << "Failed to load genomes" << ::std::endl;
			options.dumpAlignment = false;
		}

	// how many 0's should be padded?
	int pzeros = 0;
	for (unsigned l = length(reads); l > 9; l = l / 10)
		++pzeros;

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

	if (options.hammingOnly)
		scoreType.data_mismatch = -1;
	resize(rows(align), 2);

	::std::ofstream file;
	::std::ostringstream fileName;
	if (*options.output != 0)
		fileName << options.output;
	else
		fileName << readFName << ".result";

	file.open(fileName.str().c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
	if (!file.is_open()) {
		::std::cerr << "Failed to open output file" << ::std::endl;
		return;
	}

	
	maskDuplicates(matches);
	if (options.outputFormat > 0
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	 && !options.maqMapping
#endif
	)
	{
		// match statistics
		unsigned maxErrors = (int)(options.errorRate * maxReadLength);
		if (maxErrors > 10) maxErrors = 10;
		resize(stats, maxErrors + 1);
		for (unsigned i = 0; i <= maxErrors; ++i)
			fill(stats[i], length(reads), 0);
		countMatches(matches, stats);
	}
#ifdef RAZERS_DIRECT_MAQ_MAPPING
	if(options.maqMapping)
		compactMatches(matches, stats, options, false);
	else	 
#endif
	compactMatches(matches, stats, options);

	switch (options.sortOrder) {
		case 0:
			::std::sort(
				begin(matches, Standard()),
				end(matches, Standard()), 
				LessRNoGPos<TMatch>());
			break;

		case 1:
			::std::sort(
				begin(matches, Standard()),
				end(matches, Standard()), 
				LessGPosRNo<TMatch>());
			break;
			
	}
	
	typename Iterator<TMatches, Standard>::Type	it = begin(matches, Standard());
	typename Iterator<TMatches, Standard>::Type	itEnd = end(matches, Standard());

	

#ifdef RAZERS_DIRECT_MAQ_MAPPING
	typename Iterator<TMatches, Standard>::Type     bestIt;
	typename Iterator<TMatches, Standard>::Type     endIt;

	int maxSeedErrors = (int)(options.errorRate*options.artSeedLength)+1;
	String<int> cooc;
	if(options.maqMapping)
		countCoocurrences(matches,cooc,options);

#endif
	
	Dna5String gInf;
	

	switch (options.outputFormat) 
	{
		case 0:	// Razer Format
			for(; it != itEnd; ++it) 
			{
				unsigned	readLen = length(reads[(*it).rseqNo]);
				double		percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);

				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << readIDs[(*it).rseqNo];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << readName << '#' << ::std::setw(pzeros) << (*it).rseqNo + 1;
						break;

					// 2..filename is the read sequence itself
					case 2:
						file << reads[(*it).rseqNo];
				}

				file << ',' << options.positionFormat << ',' << readLen << ',' << (*it).orientation << ',';

				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[(*it).gseqNo];
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1;
				}

				file << ',' << ((*it).gBegin + options.positionFormat) << ',' << (*it).gEnd << ',' << ::std::setprecision(5) << percId << ::std::endl;

				if (options.dumpAlignment) {
					assignSource(row(align, 0), reads[(*it).rseqNo]);
					assignSource(row(align, 1), infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd));
					if ((*it).orientation == 'R')
						reverseComplementInPlace(source(row(align, 1)));
					globalAlignment(align, scoreType);
					dumpAlignment(file, align);
				}

			}
			break;


		case 1:	// Enhanced Fasta Format
			for(; it != itEnd; ++it) 
			{
				unsigned	readLen = length(reads[(*it).rseqNo]);
				double		percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);

				::std::string fastaID;
				assign(fastaID, readIDs[(*it).rseqNo]);

				int id = (*it).rseqNo;
				int fragId = id;

				size_t left = fastaID.find_first_of('[');
				size_t right = fastaID.find_last_of(']');
				if (left != fastaID.npos && right != fastaID.npos && left < right) 
				{
					fastaID.erase(right);
					fastaID.erase(0, left + 1);
					replace(fastaID.begin(), fastaID.end(), ',', ' ');
					size_t pos = fastaID.find("id=");
					if (pos != fastaID.npos) {
						::std::istringstream iss(fastaID.substr(pos + 3));
						iss >> id;
					}
					pos = fastaID.find("fragId=");
					if (pos != fastaID.npos) {
						::std::istringstream iss(fastaID.substr(pos + 7));
						iss >> fragId;
					}
				}

				if ((*it).orientation == 'F')
					// forward strand
					file << '>' << ((*it).gBegin + options.positionFormat) << ',' << (*it).gEnd;
				else
					// reverse strand (switch begin and end)
					file << '>' << (*it).gEnd << ',' << ((*it).gBegin + options.positionFormat);
					
				unsigned ambig = 0;
				for (unsigned i = 0; i <= (*it).editDist && i < length(stats); ++i)
					ambig += stats[i][(*it).rseqNo];
				
				file << "[id=" << id << ",fragId=" << fragId;
				file << ",errors=" << (*it).editDist << ",percId=" << ::std::setprecision(5) << percId;
				file << ",ambiguity=" << ambig << ']' << ::std::endl;

				file << reads[(*it).rseqNo] << ::std::endl;
			}
			break;


		case 2:	// Eland Format
			for(unsigned readNo = 0; readNo < length(reads); ++readNo)
			{
				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << '>' << readIDs[readNo] << '\t';
						break;

					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << readName << '#' << ::std::setw(pzeros) << readNo + 1  << '\t';
						break;
				}

				if (it == itEnd || readNo < (*it).rseqNo)
				{
					if (!empty(reads[readNo]))
						file << reads[readNo] << "\tNM\t0\t0\t0" << ::std::endl;
					else
					{
						for (unsigned i = 0; i < maxReadLength; ++i)
							file << '.';
						file << "\tQC\t0\t0\t0" << ::std::endl;
					}
				} else
				{
					file << reads[readNo] << '\t';
					unsigned bestMatches = 1;
					if ((unsigned)(*it).editDist < length(stats))
						bestMatches = stats[(*it).editDist][readNo];
					
					if (bestMatches == 0) file << '?';	// impossible
					if (bestMatches == 1) file << 'U';	// unique best match
					if (bestMatches >  1) file << 'R';	// non-unique best matches
					
					file << (*it).editDist << '\t' << stats[0][readNo] << '\t' << stats[1][readNo] << '\t' << stats[2][readNo];
					
					if (bestMatches == 1)
					{
						file << '\t';
						switch (options.genomeNaming)
						{
							// 0..filename is the read's Fasta id
							case 0:
								file << genomeIDs[(*it).gseqNo];
								break;

							// 1..filename is the read filename + seqNo
							case 1:
								file.fill('0');
								file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1;
						}
						
						if ((*it).orientation == 'F')
							file << '\t' << ((*it).gBegin + options.positionFormat) << "\tF\t..";
						else
							file << '\t' << (*it).gEnd << "\tR\t..";

						if ((*it).editDist > 0 && options.dumpAlignment && options.hammingOnly) 
						{
							gInf = infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd);
							if ((*it).orientation == 'R')
								reverseComplementInPlace(gInf);
							for (unsigned i = 0; i < length(gInf); ++i)
								if ((options.compMask[ordValue(reads[readNo][i])] & 
									options.compMask[ordValue(gInf[i])]) == 0)
									file << '\t' << i + 1 << gInf[i];
						}
					}
					file << ::std::endl;
					++it;
				}
			}
			break;
		case 3: // GFF:  printf "$chr $name_$format read $pos %ld . $dir . ID=$col[0]$unique$rest\n",$pos+$len-1;
			for(unsigned readNo = 0; readNo < length(reads); ++readNo)
			{
				if (it == itEnd || readNo < (*it).rseqNo)
				{
					continue; //no matches
				}
				int unique = 1;
#ifdef RAZERS_DIRECT_MAQ_MAPPING
				int mappingQuality = 0;
				bestIt = it;
				endIt = it;
				int bestQualSum = (*it).mScore;
				
				if(options.maqMapping)
				{
					bool mappingQualityFound = false;
					double qualTerm1 = 1000.0,qualTerm2 = 1000.0;
					
					
					int bestDist = (*it).editDist;
					int kPrime = (*it).seedEditDist;
					
					int secondBestDist = -1;
					int secondBestQualSum = 1000;
					unsigned secondBestMatches = 0;
					
					++it;
					if(it != itEnd && (*it).rseqNo == readNo)
					{
						secondBestQualSum = (*it).mScore;
						secondBestDist = (*it).editDist;
						secondBestMatches = stats[secondBestDist][readNo];
						if(secondBestDist<=bestDist) unique=0;
						endIt = it;
					}
					else --it;

					if((bestQualSum==secondBestQualSum) || (kPrime>maxSeedErrors))
						mappingQualityFound = true;   //mq=0
					else{
						if(secondBestDist == -1) qualTerm1 = 99;
						else
						{
							qualTerm1 = secondBestQualSum - bestQualSum - 4.343 * log(secondBestMatches);
							if (secondBestDist - bestDist <= 1 && qualTerm1 > options.mutationRateQual) qualTerm1 = options.mutationRateQual;
						}
						float avgSeedQual = 0.0;
						if(!mappingQualityFound)
						{
							//TODO!!! generalize and adapt to razers lossrates
							int kPrimeLoss = 4;
							qualTerm2 = kPrimeLoss + cooc[maxSeedErrors-kPrime];
							
							for(unsigned j = 0; j<options.artSeedLength; ++j)
							{
								int q = (int)((unsigned char)(reads[readNo][j])>>3);
								if(q>options.mutationRateQual) q = options.mutationRateQual;
								avgSeedQual+=q;
							}
							avgSeedQual/=options.artSeedLength;
							if((avgSeedQual-=13)>0) qualTerm2 += ((maxSeedErrors-kPrime)*(avgSeedQual));
						}
					}
					if (!mappingQualityFound)
						mappingQuality = (qualTerm1<qualTerm2) ? qualTerm1:qualTerm2;
					if (mappingQuality < 0) mappingQuality = 0;
				}
				it=bestIt;
#else
				unsigned bestMatches = 1;
				if ((unsigned)(*it).editDist < length(stats))
					bestMatches = stats[(*it).editDist][readNo];
				if (bestMatches !=  1)
				{
					++it;
					continue; // TODO: output non-unique matches
				}
				unsigned readLen = length(reads[(*it).rseqNo]);
				double percId = 100.0 * (1.0 - (double)(*it).editDist / (double)readLen);
#endif
				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[(*it).gseqNo] <<'\t';
						break;
					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << gnoToFileMap[(*it).gseqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[(*it).gseqNo].second + 1 << '\t';
						break;
				}
				file <<  options.runID << "_razers\tread\t";
				file << '\t' << ((*it).gBegin + options.positionFormat) << '\t' << (*it).gEnd << '\t';
//				if ((*it).orientation == 'F')
//					file << '\t' << ((*it).gBegin + options.positionFormat) << '\t' << (*it).gEnd <<'\t';
//				else
//					file << '\t' << (*it).gEnd << '\t'<<((*it).gBegin + options.positionFormat)<< '\t';
#ifdef RAZERS_DIRECT_MAQ_MAPPING
				file << mappingQuality << "\t";
#else
				file << percId << "\t";
#endif
				if ((*it).orientation == 'F')
					file << '+' << '\t' << '.' <<'\t';
				else
					file << '-' << '\t' << '.' <<'\t';

				switch (options.readNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << "ID=" <<readIDs[readNo];
						break;
					
					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << "ID=" << readName << '#' << ::std::setw(pzeros) << readNo + 1;
						break;
				}
				file << ";unique=" << unique; 
				if ((*it).editDist > 0) 
				{
					file << ";mutations=";
					if (options.dumpAlignment && options.hammingOnly)
					{
						gInf = infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd);
						if ((*it).orientation == 'R')
							reverseComplementInPlace(gInf);
						bool first = true;
						for (unsigned i = 0; i < length(gInf); ++i)
							if ((options.compMask[ordValue((Dna5)reads[readNo][i])] & 
								options.compMask[ordValue(gInf[i])]) == 0)
							{
		//						if(first){ file << i + 1 << gInf[i]; first = false;}
		//						else file <<','<< i + 1 << gInf[i];
								if(first){ file << i + 1 << (Dna5)reads[readNo][i]; first = false;}
								else file <<','<< i + 1 << (Dna5)reads[readNo][i];
							}
					}
					else
					{
						file << (*it).editDist;
					}
				}
#ifdef RAZERS_DIRECT_MAQ_MAPPING
				if(options.maqMapping)
				{
					file << ";mqs=" << bestQualSum;
					file << ";read=";
					for(unsigned j=0;j<length(reads[readNo]);++j)
					{
						file << (Dna5)reads[readNo][j];
					}
					file << ";quality=";
					for(unsigned j=0;j<length(reads[readNo]);++j)
					{
						file << (char) ((int)((unsigned char)reads[readNo][j] >> 3) + 33);
					}
				}	
				it = endIt;
#endif
				file << ::std::endl;
				++it;
			}
			break;
	}

	file.close();

	// get empirical error distribution
	if (!errorPrbFileName.empty() && maxReadLength > 0)
	{
		file.open(errorPrbFileName.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
		if (file.is_open())
		{
			String<long double> posError;
			unsigned unique = 0;
			unsigned insertions = 0;
			unsigned deletions = 0;
			fill(posError, maxReadLength, 0);
			
			if (options.hammingOnly)
				unique = getErrorDistribution(posError, matches, reads, genomes, options);
			else
			{
				unique = getErrorDistribution(posError, insertions, deletions, matches, reads, genomes, options);
				::std::cerr << "insertProb: " << (double)insertions / ((double)length(posError) * (double)unique) << ::std::endl;
				::std::cerr << "deleteProb: " << (double)deletions / ((double)length(posError) * (double)unique) << ::std::endl;
			}

			file << (double)posError[0] / (double)unique;
			for (unsigned i = 1; i < length(posError); ++i)
				file << '\t' << (double)posError[i] / (double)unique;
			file << ::std::endl;
			file.close();
		} else
			::std::cerr << "Failed to open error distribution file" << ::std::endl;
	}

	options.timeDumpResults = SEQAN_PROTIMEDIFF(dump_time);

	if (options._debugLevel >= 1)
		::std::cerr << "Dumping results took             \t" << options.timeDumpResults << " seconds" << ::std::endl;
}




///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef RAZERS_DUMP_SNPS
/*
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
	typename TGenomeNames,
	typename TReads,
	typename TReadQualities,
	typename TReadNames,
	typename TSpec
>
void dumpSNPs(
	TMatches &matches,							// forward/reverse matches
	TGenomeNames const &genomeIDs,				// Read names (read from Fasta file, currently unused)
	StringSet<CharString> &genomeFileNameList,	// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > & gnoToFileMap, //map to retrieve genome filename and sequence number within that file
	TReads const &reads,						// Read sequences
	TReadQualities const &readQualities,				// Read quality values sequences
	TReadNames const &,					// Read names (read from Fasta file, currently unused)
	::std::string readFName,					// read name (e.g. "reads.fa")
	RazerSOptions<TSpec> &options)
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

	// matches need to be ordered accordign to genome position
	if(options.sortOrder != 1) 
		::std::sort(begin(matches, Standard()),	end(matches, Standard()), LessGPosRNo<TMatch>());		

	//	options.maxHits = 1;		// only take into account best unique matches

	SEQAN_PROTIMESTART(dump_time);

	// load Genome sequences for alignment dumps
	TGenomeSet genomes;
	if (!loadGenomes(genomes, genomeFileNameList)) 
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

//	maskDuplicates(matches);
//	compactMatches(matches, options);

	typedef typename Iterator<TMatches, Standard>::Type TMatchIterator;

	TMatchIterator matchIt = begin(matches, Standard());
	TMatchIterator matchItEnd = end(matches, Standard());	
	
	double priorHet = 0.001;
	double priorHomo = (1.0 - priorHet)/2.0;

	//collect candidate SNP positions  
	::std::set<unsigned> candidates;
	while(matchIt != matchItEnd) 
	{
		unsigned currSeqNo = (*matchIt).gseqNo;
		if(options._debugLevel > 0) ::std::cout << "Scanning genome #" << currSeqNo << " for SNPs..." << ::std::endl;
		TMatchIterator currSeqMatchItBegin = matchIt;
		while(matchIt != matchItEnd)
		{
			if((*matchIt).gseqNo != currSeqNo) break;
			if ((*matchIt).editDist > 0) 
			{
				Dna5String gInf = infix(genomes[(*matchIt).gseqNo], (*matchIt).gBegin, (*matchIt).gEnd);
				if ((*matchIt).orientation == 'R')
					reverseComplementInPlace(gInf);
				
				for (unsigned i = 0; i < length(gInf); ++i)
				{
					if ((options.compMask[ordValue(reads[(*matchIt).rseqNo][i])] & options.compMask[ordValue(gInf[i])]) == 0)
						if ((*matchIt).orientation == 'R')
							candidates.insert((*matchIt).gEnd-i-1);
						else
							candidates.insert((*matchIt).gBegin+i);
				}
			}
			++matchIt;
		}
		TMatchIterator currSeqMatchItEnd = matchIt;

		if(options._debugLevel > 1) ::std::cout << "Verify candidates...\n";
		if(options._debugLevel > 1) ::std::cout << matchItEnd-currSeqMatchItBegin << " matches left.\n";
		
		matchIt = currSeqMatchItBegin;
		typename ::std::set<unsigned>::iterator candidateIt = candidates.begin();
		typename ::std::set<unsigned>::iterator candidateItEnd = candidates.end();
		
		for(; candidateIt != candidateItEnd; ++candidateIt)
		{
			//candidate SNP position
			unsigned candidatePos = *candidateIt;

			String<unsigned> count;
			fill(count,5,0);

			//find range of relevant read matches
			while(matchIt != currSeqMatchItEnd && (*matchIt).gEnd < candidatePos) ++matchIt;	
			TMatchIterator matchRangeBegin = matchIt;
			while(matchIt != currSeqMatchItEnd && (*matchIt).gBegin <= candidatePos) ++matchIt;
			TMatchIterator matchRangeEnd = matchIt;

			unsigned coverage = matchRangeEnd-matchRangeBegin;
			//reference sequence nucleotide
			Dna5 refBase = genomes[(*matchRangeBegin).gseqNo][candidatePos];
			for(matchIt = matchRangeBegin; matchIt != matchRangeEnd; ++matchIt)
			{
				Dna5 candidateBase;
				if ((*matchIt).orientation == 'R')
				{
					FunctorComplement<Dna5> f;
					candidateBase = f((Dna5)reads[(*matchIt).rseqNo][(*matchIt).gEnd - candidatePos - 1]);
				}
				else
					candidateBase = reads[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin];
				++count[ordValue(candidateBase)];
			}
			//do statistics and output snp
			int allele1 = -1;
			int allele2 = -1;
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
			if(allele2 < 0) // reads all vote for one nucleotide (non-ref)
				allele2 = refAllele;


			double avgErrProb = 0.05;
			int n = count[allele1] + count[allele2]; // ignore all but the two most common bases
			// TODO: minimum coverage??
			if(n<10)
			{
				matchIt = matchRangeBegin;
				continue;
			
			}
			// 
			// argmax P(g|D)=P(D|g)*P(g)/P(D)
			//    g
			// 
			//todo: alphas so wie in maq
	
			double pHet = over2(n,count[allele2])/(1 << n); //vorsicht!
			//double pHomoAllele1 = alpha_n_k;
			//double pHomoAllele1 = binomial(n,count[allele2],avgErrProb);//uniform and independent
			double pHomoAllele1 = over2(n,count[allele2]) * pow(avgErrProb,(double)count[allele2]) * pow(1.0-avgErrProb,(double)count[allele1]) ;//uniform and independent
			//double pHomoAllele2 = alpha_n_nk;
			//double pHomoAllele2 = binomial(n,count[allele1],avgErrProb);//uniform and independent
			double pHomoAllele2 = over2(n,count[allele1]) * pow(avgErrProb,(double)count[allele1]) * pow(1.0-avgErrProb,(double)count[allele2]);//uniform and independent

			pHet *= priorHet; 
			pHomoAllele1 *= priorHomo;
			pHomoAllele2 *= priorHomo;
	
			bool hetSNP = false;
			bool homoSNP = false;
			
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
				if(homoSNP) ::std::cout << "Vorsicht! HomoSnp auf allele2!\n";
			}
			if(allele1 != refAllele && allele2 != refAllele)
				::std::cout <<"+";
			//double lambda = coverage * 0.05;//=avg_err_prob;  
			//double pHet = p_value(count[allele2],lambda);
			//if(coverage > 4 && pValHet <= options.testLevel) // heterozygote SNP
			//	hetSNP = true;
			//if(coverage > 4 && pValHet > options.testLevel && (unsigned)allele1 != ordValue(genomes[(*matchIt).gseqNo][candidatePos])) // homozygote SNP
			//	homoSNP = true;

			if(hetSNP || homoSNP)
			{
				//
				switch (options.genomeNaming)
				{
					// 0..filename is the read's Fasta id
					case 0:
						file << genomeIDs[currSeqNo] << " ";
						break;
					// 1..filename is the read filename + seqNo
					case 1:
						file.fill('0');
						file << gnoToFileMap[currSeqNo].first << '#' << ::std::setw(gzeros) << gnoToFileMap[currSeqNo].second + 1 << " ";
				}
				file << candidatePos << " " << genomes[(*matchIt).gseqNo][candidatePos] <<" ";
				file << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << " " << count[4];
				if(homoSNP) file <<  " <" << Dna(allele1) << "," << Dna(allele1) << "> "  << pHomoAllele1<< ::std::endl; //vorsicht!
				else file <<  " <" << Dna(allele1) << "," << Dna(allele2) << "> "<< pHet << ::std::endl; 

				if(options._debugLevel > 0) 
				{
					::std::cout << "CandidateSNP: " << candidatePos << " (covered by " << coverage << " matches):" << ::std::endl;
					::std::cout << "Genome is " << genomes[(*matchIt).gseqNo][candidatePos];
					::std::cout << "\treads are "<<count[0] << "," << count[1] << "," << count[2] << "," << count[3] << "," << count[4];
					if(homoSNP) ::std::cout <<  " <" << Dna(allele1) << "," << Dna(allele1) << "> " << pHomoAllele1 << ::std::endl; //vorsicht!
					else ::std::cout <<  " <" << Dna(allele1) << "," << Dna(allele2) << "> " << pHet << ::std::endl;
				}
			}
			matchIt = matchRangeBegin;
		}
		matchIt = currSeqMatchItEnd;
	}
	file.close();
	if(options._debugLevel > 1) ::std::cout << "Detecting SNPs took " << SEQAN_PROTIMEDIFF(dump_time) << " seconds." << ::std::endl;
	return;

}
*/
#endif
}
#endif

