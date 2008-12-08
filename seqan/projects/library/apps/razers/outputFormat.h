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
template <typename TErrDistr, typename TMatches, typename TReads, typename TGenomes, typename TSpec>
inline unsigned
getErrorDistribution(
	TErrDistr &posError, 
	TMatches &matches, 
	TReads &reads, 
	TGenomes &genomes, 
	RazerSOptions<TSpec> &options)
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


//////////////////////////////////////////////////////////////////////////////
// Output matches
template <
	typename TMatches,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames,
	typename TSpec
>
void dumpMatches(
	TMatches &matches,							// forward/reverse matches
	TGenomeNames const &genomeIDs,				// Read names (read from Fasta file, currently unused)
	StringSet<CharString> &genomeFileNameList,	// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > &gnoToFileMap, //map to retrieve genome filename and sequence number within that file
	TReads const &reads,						// Read sequences
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

	String<String<unsigned> > stats;

	maskDuplicates(matches);
	if (options.outputFormat > 0)
	{
		// match statistics
		unsigned maxErrors = (int)(options.errorRate * maxReadLength);
		if (maxErrors > 10) maxErrors = 10;
		resize(stats, maxErrors + 1);
		for (unsigned i = 0; i <= maxErrors; ++i)
			fill(stats[i], length(reads), 0);
		countMatches(matches, stats);
	}
	compactMatches(matches, options);

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
			Dna5String gInf;
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





//////////////////////////////////////////////////////////////////////////////
// Output SNPs
// Simple SNP calling, under construction
template <
	typename TMatches,
	typename TGenomeNames,
	typename TReads,
	typename TReadNames,
	typename TSpec
>
void dumpSNPs(
	TMatches &matches,							// forward/reverse matches
	TGenomeNames const &genomeIDs,				// Read names (read from Fasta file, currently unused)
	StringSet<CharString> &genomeFileNameList,	// list of genome names (e.g. {"hs_ref_chr1.fa","hs_ref_chr2.fa"})
	::std::map<unsigned,::std::pair< ::std::string,unsigned> > &, //map to retrieve genome filename and sequence number within that file
	TReads const &reads,						// Read sequences
	TReadNames const &readIDs,					// Read names (read from Fasta file, currently unused)
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
	
	//collect candidate SNP positions   //TODO: consider multiple genomes!!!
	::std::set<unsigned> candidates;
	for(; matchIt != matchItEnd; ++matchIt) 
	{
		if ((*matchIt).editDist > 0) 
		{
			::std::cout << "dist > 0" << ::std::endl;
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
	}


	typename ::std::set<unsigned>::iterator candidateIt = candidates.begin();
	typename ::std::set<unsigned>::iterator candidateItEnd = candidates.end();
	matchIt = begin(matches, Standard());
	
	for(; candidateIt != candidateItEnd; ++candidateIt)
	{
		//candidate SNP position
		unsigned candidatePos = *candidateIt;
		::std::cout << "CandidateSNP: " << candidatePos << ::std::endl;

		String<unsigned> count;
		fill(count,5,0);

		//find range of relevant read matches
		while(matchIt != matchItEnd && (*matchIt).gEnd < candidatePos) ++matchIt;	
		TMatchIterator matchRangeBegin = matchIt;
		while(matchIt != matchItEnd && (*matchIt).gBegin <= candidatePos) ++matchIt;	
		TMatchIterator matchRangeEnd = matchIt;

		unsigned coverage = matchRangeEnd-matchRangeBegin;
		::std::cout << "SNP covered by " << coverage << " matches." << ::std::endl;
		//reference sequence nucleotide
		Dna5 refBase = genomes[(*matchRangeBegin).gseqNo][candidatePos];
		for(matchIt = matchRangeBegin; matchIt != matchRangeEnd; ++matchIt)
		{
			Dna5 candidateBase;
			if ((*matchIt).orientation == 'R')
			{
				candidateBase = reads[(*matchIt).rseqNo][(*matchIt).gEnd - candidatePos - 1];
				::std::cout << "Position on read " << readIDs[(*matchIt).rseqNo]<< ": " << (*matchIt).gEnd - candidatePos - 1 << ::std::endl;
				Dna5String temp;
				appendValue(temp,candidateBase);
				reverseComplementInPlace(temp);
				candidateBase = temp[0];
				//                              reverseComplementInPlace(candidateBase);
			}
			else
				candidateBase = reads[(*matchIt).rseqNo][candidatePos - (*matchIt).gBegin];
			++count[ordValue(candidateBase)];
		}
		//do statistics and output snp
		::std::cout << count[0] << ::std::endl;
		::std::cout << count[1] << ::std::endl;
		::std::cout << count[2] << ::std::endl;
		::std::cout << count[3] << ::std::endl;
		::std::cout << count[4] << ::std::endl;
		matchIt = matchRangeBegin;
		
	}
	::std::cout << "Detecting SNPs took " << SEQAN_PROTIMEDIFF(dump_time) << " seconds." << ::std::endl;
	return;

}



}
#endif

