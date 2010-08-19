 /*==========================================================================
  SNP Calling routine of RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by Anne-Katrin Emde

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

#ifndef SEQAN_HEADER_COMPAREVAR_H
#define SEQAN_HEADER_COMPAREVAR_H


#include <iostream>
#include <fstream>
#include <cmath>

#define EPS_PREC 0.00000001

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Default options


//____________________________________________________________________________
// Global Parameters



struct IndelCheck{};

template<typename TSpec = IndelCheck>
struct IndelCompareOptions
{
	const char *output; 	    // output file for statistics and shared indels
	const char *inputReference;	// reference indels in gff format
	const char *inputPredicted;	// predicted indels in gff format

	int positionTolerance;		// max. deviation of predicted indel position from reference indel position
	double sizeTolerance;		    // max. deviation of predicted indel size from reference indel size (in percent)
	bool sequenceContext;
	
	int _debugLevel;
	
	IndelCompareOptions()
	{
		output = "";		
		inputReference = "";		
		inputPredicted = "";		
		
		positionTolerance  = 10;  
		sizeTolerance = 0.0;
		sequenceContext = false;
		
		_debugLevel = 0;
		
	}
	
};

struct IndelInfo{
	
	unsigned genomeId;
	unsigned originalPos;
	unsigned simPos;
	int indelSize;
	bool duplication;
	CharString idStr;
	
};

//____________________________________________________________________________

template <typename TIndel>
struct LessGPos : public ::std::binary_function < TIndel, TIndel, bool >
{
	inline bool operator() (TIndel const &a, TIndel const &b) const
	{
		if (a.genomeId < b.genomeId ) return true;
		if (a.genomeId > b.genomeId ) return false;

		if (a.originalPos < b.originalPos ) return true;
		if (a.originalPos > b.originalPos ) return false;

		return a.indelSize < b.indelSize;
	}
};

//____________________________________________________________________________

template <typename TIndel, typename TOptions>
bool compareIndelPair(TIndel &refIndel, TIndel &predIndel, TOptions &options)
{
	
	int sizeTol = abs(refIndel.indelSize) * options.sizeTolerance;
	
	if(refIndel.indelSize - sizeTol <= predIndel.indelSize && predIndel.indelSize <= refIndel.indelSize + sizeTol ) 
		return true;
	return false;
	
}

// write one gff line
template <typename TFile, typename TOptions>
bool write(TFile &file, IndelInfo &refIndel, CharString &genomeID, CharString &tagAppend, TOptions &)
{
	if(!file.is_open()) return false;
	file << genomeID << "\tvariantCmp\t";
	
	if(refIndel.indelSize < 0) file << "insertion\t" << refIndel.originalPos + 1 << '\t' << refIndel.originalPos + 1 ;
	else file << "deletion\t"<< refIndel.originalPos + 1 << '\t' << refIndel.originalPos + refIndel.indelSize ;
	
	file << "\t+\t.\t.\t";
	file << "ID=" << refIndel.idStr << ";size=" << refIndel.indelSize;
	if(refIndel.duplication) file << ";duplication=1";
	
	file << toCString(tagAppend) << std::endl;
	
	return true;
	
}

//////////////////////////////////////////////////////////////////////////////
// 1) build interval tree of reference indels
// 2) go through predicted indels one by one
template <
	typename TIndelSet,
	typename TGenome,
	typename TGenomeIDs,
	typename TOptions
>
int compareIndels(
	TIndelSet			&refIndels,		  // reference indels
	TIndelSet			&predIndels, // predicted indels
	TGenome				&genomes,
	TGenomeIDs			&genomeIDs,
	TOptions 			&options)	  	  // options
{
	
	typedef typename Value<TIndelSet>::Type    TIndel;
	typedef typename Iterator<TIndelSet>::Type TIndelIt;

	::std::ostringstream fileName;
	if (*options.output != 0)
		fileName << options.output;
	else
		fileName << options.inputPredicted << ".overlap.gff";
	
	::std::ofstream file;
	file.open(fileName.str().c_str(), ::std::ios_base::out | ::std::ios_base::trunc);
	if (!file.is_open()) {
		::std::cerr << "\nFailed to open output file" << ::std::endl;
		return 1;
	}


	// a set to keep track of which of the reference indels were matched
	::std::set<int> refFoundSet;
	
	::std::sort(begin(predIndels),end(predIndels),LessGPos<TIndel>());
	::std::sort(begin(refIndels),end(refIndels),LessGPos<TIndel>());

	TIndelIt refIndelIt = begin(refIndels);
	TIndelIt refIndelsEnd = end(refIndels);
	TIndelIt predIndelIt = begin(predIndels);
	TIndelIt predIndelsEnd = end(predIndels);
	TIndelIt currPredGenomeEnd,currPredGenomeBegin,currRefGenomeEnd,currRefGenomeBegin;
	
	int TP = 0;
	
	// foreach genomeID
	//   build interval tree from ref indels
	//   while predictedIndels on current genomeID
	//     check indel interval +- positionTolerance
	//     foreach relevant reference indel
	//       check size +- sizeTolerance
	//       if(match) output refIndel-predictedIndel-combi
	//                 and increase TP counter
	
	for(unsigned i = 0; i < length(genomeIDs); ++i)
	{
		//skip ahead if necessary
		while(refIndelIt != refIndelsEnd && (*refIndelIt).genomeId < i)
			++refIndelIt;
		while(predIndelIt != predIndelsEnd && (*predIndelIt).genomeId < i)
			++predIndelIt;
		
		// get range of ref variants that are on current chromosome
		TIndelIt currRefGenomeBegin = refIndelIt;
		while(refIndelIt != refIndelsEnd && (*refIndelIt).genomeId == i)
			++refIndelIt;
		TIndelIt currRefGenomeEnd = refIndelIt;
		if(currRefGenomeEnd - currRefGenomeBegin == 0)
			continue;
			
		// get range of predicted variants that are on current chromosome
		TIndelIt currPredGenomeBegin = predIndelIt;
		while(predIndelIt != predIndelsEnd && (*predIndelIt).genomeId == i)
			++predIndelIt;
		TIndelIt currPredGenomeEnd = predIndelIt;
		if(currPredGenomeEnd - currPredGenomeBegin == 0)
			continue;

			
		// build interval tree of all variant intervals
		typedef int TCargo;
		typedef int TValue;
		String<IntervalAndCargo<TValue,TCargo> > intervals;
		int ii = currRefGenomeBegin - begin(refIndels);
		for(refIndelIt = currRefGenomeBegin; refIndelIt != currRefGenomeEnd; ++refIndelIt)
		{
			if((*refIndelIt).indelSize > 0) 
				appendValue(intervals, IntervalAndCargo<TValue,TCargo>((*refIndelIt).originalPos,(*refIndelIt).originalPos + (*refIndelIt).indelSize, ii ));
			else
				appendValue(intervals, IntervalAndCargo<TValue,TCargo>((*refIndelIt).originalPos,(*refIndelIt).originalPos + 1, ii ));
			++ii;
		}
		IntervalTree<TValue,TCargo> itree(intervals);
		
		// for each predicted variant check if it overlaps with a reference variant
		for(predIndelIt = currPredGenomeBegin; predIndelIt != currPredGenomeEnd; ++predIndelIt)
		{
			TIndel &predIndel = *predIndelIt;
			
			// "core" interval
			TValue beginPoint = predIndel.originalPos;
			TValue endPoint = predIndel.originalPos;
			if(predIndel.indelSize > 0) endPoint += predIndel.indelSize;
			else endPoint += 1;
			
			// add tolerance in genomic position
			if(options.sequenceContext)
			{
				//TODO: do sth
			}
			else
			{
				beginPoint -= options.positionTolerance;
				endPoint += options.positionTolerance;
			}
			
			// retrieve overlapping intervals
			String<TCargo> intersectingIntervals;
			findIntervals(itree,beginPoint,endPoint,intersectingIntervals);
			
			// check whether one of the overlapping intervals is a variant match
			bool predFound = false;
			for(unsigned j = 0; j < length(intersectingIntervals); ++j)
			{
				TIndel &refIndel = refIndels[intersectingIntervals[j]];
				bool res = compareIndelPair(refIndel,predIndel,options);
				if(res)
				{
					CharString tagAppend = ";matchID=";
					append(tagAppend,predIndel.idStr);
					write(file,refIndel,genomeIDs[i],tagAppend,options);
					predFound = true;
					refFoundSet.insert(intersectingIntervals[j]);
				}
				
			}
			if(predFound) ++TP;
			
		
		}
		
	}

	// optionally append non-overlapped reference indels to output file
	file << "###################################################" << ::std::endl;
	file << "# Stats: " << ::std::endl;
	file << "# Number of reference indels: " << length(refIndels) << ::std::endl;
	file << "# Number of predicted indels: " << length(predIndels) << ::std::endl;
	file << "# Number of TP predictions  : " << TP << ::std::endl;
	file << "# Number of FP predictions  : " << length(predIndels)-TP << ::std::endl;
	file << "# Number of FN predictions  : " << length(refIndels)-TP << ::std::endl;

	
	file.close();
	// FN = length(referenceIndels) - TP
	// FP = length(predictedIndels) - TP
	// output statistics at end of file

	
	return 0;
	
}




}

#endif
