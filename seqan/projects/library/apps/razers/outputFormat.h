 /*==========================================================================
                     RazerS - Fast Mapping of Short Reads
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

	template <typename TReadMatch>
	struct LessRNoGPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gBegin < b.gBegin) return true;
			if (a.gBegin > b.gBegin) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};

	// ... to sort matches and remove duplicates with equal gEnd
	template <typename TReadMatch>
	struct LessRNoGEndPos : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// genome position and orientation
			if (a.gseqNo < b.gseqNo) return true;
			if (a.gseqNo > b.gseqNo) return false;
			if (a.gEnd   < b.gEnd) return true;
			if (a.gEnd   > b.gEnd) return false;
			if (a.orientation < b.orientation) return true;
			if (a.orientation > b.orientation) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};

	template <typename TReadMatch>
	struct LessErrors : public ::std::binary_function < TReadMatch, TReadMatch, bool >
	{
		inline bool operator() (TReadMatch const &a, TReadMatch const &b) const 
		{
			// read number
			if (a.rseqNo < b.rseqNo) return true;
			if (a.rseqNo > b.rseqNo) return false;

			// quality
			return a.editDist < b.editDist;
		}
	};


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
// Remove duplicate matches and leave at most maxHits many distanceRange
// best matches per read
template < typename TMatches >
void maskDuplicates(TMatches &matches)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal ends

	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoGEndPos<TMatch>());

	typename	TMatch::TGPos gBegin = -1;
	typename	TMatch::TGPos gEnd = -1;
	unsigned	gseqNo = -1;
	unsigned	readNo = -1;
	char		orientation = '-';

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());

	for (; it != itEnd; ++it) 
	{
		if (gEnd == (*it).gEnd && orientation == (*it).orientation &&
			gseqNo == (*it).gseqNo && readNo == (*it).rseqNo) 
		{
			(*it).orientation = '-';
			continue;
		}
		readNo = (*it).rseqNo;
		gseqNo = (*it).gseqNo;
		gEnd = (*it).gEnd;
		orientation = (*it).orientation;
	}

	//////////////////////////////////////////////////////////////////////////////
	// remove matches with equal begins

	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessRNoGPos<TMatch>());

	orientation = '-';

	it = begin(matches, Standard());
	itEnd = end(matches, Standard());

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (gBegin == (*it).gBegin && readNo == (*it).rseqNo &&
			gseqNo == (*it).gseqNo && orientation == (*it).orientation) 
		{
			(*it).orientation = '-';
			continue;
		}
		readNo = (*it).rseqNo;
		gseqNo = (*it).gseqNo;
		gBegin = (*it).gBegin;
		orientation = (*it).orientation;
	}
	
	//////////////////////////////////////////////////////////////////////////////
	// remove low quality matches

	::std::sort(
		begin(matches, Standard()),
		end(matches, Standard()), 
		LessErrors<TMatch>());
}

template < typename TMatches, typename TCounts >
void countMatches(TMatches &matches, TCounts &cnt)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	typedef typename Value<TCounts>::Type					TRow;
	typedef typename Value<TRow>::Type						TValue;
	
	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	
	unsigned readNo = -1;
	short editDist = -1;
	__int64 count = 0;
	__int64 maxVal = SupremumValue<TValue>::VALUE;

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo && editDist == (*it).editDist)
			++count;
		else
		{
			if (readNo != (unsigned)-1 && (unsigned)editDist < length(cnt))
				cnt[editDist][readNo] = (maxVal < count)? maxVal : count;
			readNo = (*it).rseqNo;
			editDist = (*it).editDist;
			count = 1;
		}
	}
	if (readNo != (unsigned)-1)
		cnt[editDist][readNo] = count;
}	

template < typename TMatches, typename TSpec >
void compactMatches(TMatches &matches, RazerSOptions<TSpec> &options)
{
	typedef typename Value<TMatches>::Type					TMatch;
	typedef typename Iterator<TMatches, Standard>::Type		TIterator;
	
	unsigned readNo = -1;
	unsigned hitCount = 0;
	unsigned hitCountCutOff = options.maxHits;
	int editDistCutOff = SupremumValue<int>::VALUE;

	TIterator it = begin(matches, Standard());
	TIterator itEnd = end(matches, Standard());
	TIterator dit = it;

	for (; it != itEnd; ++it) 
	{
		if ((*it).orientation == '-') continue;
		if (readNo == (*it).rseqNo)
		{ 
			if ((*it).editDist >= editDistCutOff || ++hitCount >= hitCountCutOff)
				continue;
		} else
		{
			readNo = (*it).rseqNo;
			hitCount = 0;
			if (options.distanceRange > 0)
				editDistCutOff = (*it).editDist + options.distanceRange;
		}
		*dit = *it;
		++dit;
	}
	resize(matches, dit - begin(matches, Standard()));
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
	TMatches &matches,					// forward/reverse matches
	TGenomeNames const &genomeIDs,		// Read names (read from Fasta file, currently unused)
	::std::string genomeFName,			// genome name (e.g. "hs_ref_chr1.fa")
	TReads const &reads,				// Read sequences
	TReadNames const &readIDs,			// Read names (read from Fasta file, currently unused)
	::std::string readFName,			// read name (e.g. "reads.fa")
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
	String<int> posError;
	fill(posError, maxReadLength, 0);

	// match statistics
	unsigned maxErrors = (int)(options.errorRate * maxReadLength);
	if (maxErrors > 10) maxErrors = 10;
	String<String<unsigned> > stats;
	resize(stats, maxErrors + 1);
	for (unsigned i = 0; i <= maxErrors; ++i)
		fill(stats[i], length(reads), 0);

	SEQAN_PROTIMESTART(dump_time);

	// load Genome sequences for alignment dumps
	TGenomeSet genomes;
	if (options.dumpAlignment || (!errorPrbFileName.empty() && options.hammingOnly))
		if (!loadGenomes(genomes, genomeFName.c_str())) {
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


	// remove the directory prefix of genomeFName and readFName
	size_t lastPos = genomeFName.find_last_of('/') + 1;
	if (lastPos == genomeFName.npos) lastPos = genomeFName.find_last_of('\\') + 1;
	if (lastPos == genomeFName.npos) lastPos = 0;
	::std::string genomeName = genomeFName.substr(lastPos);

	lastPos = readFName.find_last_of('/') + 1;
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
	countMatches(matches, stats);
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
						file << genomeName << '#' << ::std::setw(gzeros) << (*it).gseqNo + 1;
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
				for (int i = 0; i <= (*it).editDist; ++i)
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
								file << genomeName << '#' << ::std::setw(gzeros) << (*it).gseqNo + 1;
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
	if (!errorPrbFileName.empty() && options.hammingOnly)
	{
		it = begin(matches, Standard());
		itEnd = end(matches, Standard());

		unsigned unique = 0;
		for (; it != itEnd; ++it) 
		{
			if ((*it).orientation == '-') continue;

			Dna5String const &read = reads[(*it).rseqNo];
			Dna5String genome = infix(genomes[(*it).gseqNo], (*it).gBegin, (*it).gEnd);
			if ((*it).orientation == 'R')
				reverseComplementInPlace(genome);

			for (unsigned i = 0; i < length(read); ++i)
				if (genome[i] != read[i])
				{
					if (!(genome[i] == 'N' && read[i] != 'N'))
						++posError[i];
				}
			++unique;
		}
		file.open(errorPrbFileName.c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
		if (!file.is_open()) {
			::std::cerr << "Failed to open error distribution file" << ::std::endl;
			return;
		}
		for (unsigned i = 0; i < length(posError); ++i)
			file << (double)posError[i] / unique << ::std::endl;
		file.close();
	}


	options.timeDumpResults = SEQAN_PROTIMEDIFF(dump_time);

	if (options._debugLevel >= 1)
		::std::cerr << "Dumping results took             \t" << options.timeDumpResults << " seconds" << ::std::endl;
}

}

#endif

