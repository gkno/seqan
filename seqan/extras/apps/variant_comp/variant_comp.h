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
#define MAX_DUPLICATION_SIZE 400000


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
	int annotateRepeats;
	
	const char *attachTag;          // tag to attach to overlap output
	int _debugLevel;
	String<Pair<int,int> > ranges;  // seperate statistics for each range
	
	IndelCompareOptions()
	{
		output = "";		
		inputReference = "";		
		inputPredicted = "";		
		
		positionTolerance  = 10;  
		sizeTolerance = 0.0;
		sequenceContext = false;
		
		annotateRepeats = 50;	// 0 -> dont add "percentRepeat"-tag in result file
								// otherwise take +- x flanking sequence
		attachTag = "";
		_debugLevel = 0;
		appendValue(ranges,Pair<int,int>(-50,-6));
		appendValue(ranges,Pair<int,int>(-6,0));
		appendValue(ranges,Pair<int,int>(1,7));
		appendValue(ranges,Pair<int,int>(7,51));
		appendValue(ranges,Pair<int,int>(51,501));
		appendValue(ranges,Pair<int,int>(501,5001));
		
	}
	
};

enum SVTypes {
		INSERTION = 1,
		DELETION = 2,
		INVERSION = 3,
        TRANSLOCATION = 4
	};

struct IndelInfo{
	
	unsigned genomeId;
	unsigned originalPos;
	unsigned secondGenomeId;
	unsigned secondOriginalPos;
	unsigned simPos;
	int indelSize;
    int type;
	bool duplication;
	Dna5String insertionSeq;
	CharString idStr;
	CharString ninethCol;
	CharString field2;
	CharString field3;
	
	// int percentRepeat;
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
    if(refIndel.type != predIndel.type) return false;
	
	int sizeTol = int((double)abs(refIndel.indelSize) * options.sizeTolerance);
	if(!(refIndel.indelSize - sizeTol <= predIndel.indelSize && predIndel.indelSize <= refIndel.indelSize + sizeTol))
		return false; // --> doesnt match
	
	// check if begin position matches
	if((int)refIndel.originalPos - options.positionTolerance <= (int) predIndel.originalPos
			&& predIndel.originalPos <= refIndel.originalPos + options.positionTolerance)
		return true; // matches
	
	// if it is a deletion or inversion
	if(refIndel.indelSize > 0) // check if end position matches
	{
		if (((int)refIndel.originalPos+refIndel.indelSize-options.positionTolerance
				<= (int)predIndel.originalPos+predIndel.indelSize)
			&& ((int) predIndel.originalPos+predIndel.indelSize
				<= (int)refIndel.originalPos+refIndel.indelSize+options.positionTolerance))
			return true;
	}
	else // if it is an insertion
	{
		if(refIndel.duplication)
		{
			if (abs((int)refIndel.originalPos - (int)predIndel.originalPos) < -refIndel.indelSize + options.positionTolerance)
				return true;
	
		}
	}
	return false;
	
}

template <typename TIndel, typename TPos, typename TOptions>
bool compareIndelPair(TIndel &refIndel, TIndel &predIndel, TPos beginPoint, TPos endPoint, TOptions &options)
{
	int sizeTol = int((double)abs(refIndel.indelSize) * options.sizeTolerance);
	if(!(refIndel.indelSize - sizeTol <= predIndel.indelSize && predIndel.indelSize <= refIndel.indelSize + sizeTol))
		return false; // --> doesnt match
	
	// check if begin position matches
	if(refIndel.originalPos >= (unsigned)beginPoint
		&& refIndel.originalPos <= (unsigned)endPoint)
		return true; // matches
	return false;
	
}


// write one gff line
template <typename TFile, typename TOptions>
bool write(TFile &file, IndelInfo &refIndel, CharString &genomeID, CharString &tagAppend, TOptions &)
{
	if(!file.is_open()) return false;
	file << genomeID << '\t';
	
	if(!empty(refIndel.field2)) file << refIndel.field2 << '\t';
	else file << "variantCmp\t";

	if(!empty(refIndel.field3)) file << refIndel.field3 << '\t';
	else
	{
		if(refIndel.indelSize < 0) file << "insertion\t";
		if(refIndel.indelSize == 0) file << "baseexchange\t";
		if(refIndel.indelSize > 0) file << "deletion\t";
	}

	if(refIndel.indelSize <= 0) file << refIndel.originalPos + 1 << '\t' << refIndel.originalPos + 1 ;
	else file << refIndel.originalPos + 1 << '\t' << refIndel.originalPos + refIndel.indelSize ;
	
	file << "\t+\t.\t.\t";
	if(!empty(refIndel.ninethCol)) file << refIndel.ninethCol;
	else
	{
		file << "ID=" << refIndel.idStr << ";size=" << refIndel.indelSize;
		if(refIndel.duplication) file << ";duplication=1";
	}
	
	file << toCString(tagAppend) << std::endl;
	
	return true;
	
}


template<typename TIndel, typename TGenome, typename TPos>
void 
computeEir(TIndel &indel, TGenome & genome, TPos & beginPoint,TPos & endPoint)
{

	SEQAN_ASSERT(indel.indelSize > 0 || (int)length(indel.insertionSeq) == -indel.indelSize );

	if(indel.indelSize > 0)	// deletion
	{
		unsigned iDel = indel.originalPos;
		unsigned iRef = indel.originalPos + indel.indelSize;
		while(iRef < length(genome) && genome[iRef]==genome[iDel])
		{
			++iRef; ++iDel;
		}
		endPoint = iRef;
		iDel = indel.originalPos + indel.indelSize - 1;
		iRef = _max(0,(int)indel.originalPos - 1);
		while(iRef > 0 && genome[iRef] == genome[iDel])
		{
			--iRef; --iDel;
		}
		beginPoint = iRef;
	}
	else
	{
		SEQAN_ASSERT_EQ((int)length(indel.insertionSeq),-indel.indelSize);
		unsigned iIns = 0; // relative to insertion sequence
		unsigned iRef = indel.originalPos;
//		std::cout << "right: iRef=" << iRef << " ->" << genome[iRef] <<"   iIns=" <<iIns << " ->" << indel.insertionSeq[iIns] << std::endl;
		while(iRef < length(genome) && genome[iRef]==indel.insertionSeq[iIns])
		{
//			std::cout << "right: iRef=" << iRef << " ->" << genome[iRef] <<"   iIns=" <<iIns << " ->" << indel.insertionSeq[iIns] << std::endl;
			++iRef; ++iIns;
			if(iIns >= length(indel.insertionSeq)) iIns = 0;
		}
		endPoint = iRef;
		iIns = -indel.indelSize - 1;
		iRef = _max(0,(int)indel.originalPos - 1);
//		std::cout << "left: iRef=" << iRef << " ->" << genome[iRef] <<"   iIns=" <<iIns << " ->" << indel.insertionSeq[iIns] << std::endl;
		while(iRef > 0 && genome[iRef] == indel.insertionSeq[iIns])
		{
//			std::cout << "left: iRef=" << iRef << " ->" << genome[iRef] <<"   iIns=" <<iIns << " ->" << indel.insertionSeq[iIns] << std::endl;
			if(iIns == 0) iIns = -indel.indelSize;
			--iRef; --iIns;
		}
		beginPoint = iRef;
	}
//	std::cout <<"begin=" << beginPoint << " end=" << endPoint << std::endl;
	

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
	typename ::std::set<int>::iterator foundSetIt;
	
	::std::sort(begin(predIndels),end(predIndels),LessGPos<TIndel>());
	::std::sort(begin(refIndels),end(refIndels),LessGPos<TIndel>());

	TIndelIt refIndelIt = begin(refIndels);
	TIndelIt refIndelsEnd = end(refIndels);
	TIndelIt predIndelIt = begin(predIndels);
	TIndelIt predIndelsEnd = end(predIndels);
	TIndelIt currPredGenomeEnd,currPredGenomeBegin,currRefGenomeEnd,currRefGenomeBegin;
	
	int TP = 0;

	// rememer stats for each range
	String<unsigned> rangeCountRef,rangeCountPred,rangeTP,rangeFound;
    unsigned predInversions = 0, refInversions = 0, inversionsTP = 0, inversionsFound = 0;
	resize(rangeCountRef,length(options.ranges),0);
	resize(rangeCountPred,length(options.ranges),0);
	resize(rangeTP,length(options.ranges),0);         // to calculate FP
	resize(rangeFound,length(options.ranges),0);  // to calculate FN
	for(unsigned j = 0; j < length(options.ranges); ++j)
	{
		for(TIndelIt it = begin(refIndels); it != end(refIndels); ++it)
        {
            if((*it).type == INVERSION) ++refInversions;
            else
    			if(options.ranges[j].i1 <= (*it).indelSize && (*it).indelSize < options.ranges[j].i2)
	    			++rangeCountRef[j];
        }
		for(TIndelIt it = begin(predIndels); it != end(predIndels); ++it)
        {
            if((*it).type == INVERSION) ++predInversions;
            else
	    		if(options.ranges[j].i1 <= (*it).indelSize && (*it).indelSize < options.ranges[j].i2)
    				++rangeCountPred[j];
        }
	}
	CharString tagAttach = ";";
	if(*options.attachTag != 0) append(tagAttach,options.attachTag);
	else append(tagAttach,"matchID");
	append(tagAttach,"=");

	// foreach genomeID
	//   build interval tree from ref indels
	//   while predictedIndels on current genomeID
	//     check indel interval +- positionTolerance
	//     foreach relevant reference indel
	//       check size +- sizeTolerance
	//       if(match) output refIndel-predictedIndel-combi
	//                 and increase TP counter
	
	if(options._debugLevel > 0) std::cout << "Starting to compare indels..." <<std::endl;
	for(unsigned i = 0; i < length(genomeIDs); ++i)
	{
		if(options._debugLevel > 1 ) std::cout << "." <<std::flush;

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
//		if(currRefGenomeEnd - currRefGenomeBegin == 0)
//			continue;
			
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
			{
				if((*refIndelIt).duplication)
					appendValue(intervals, IntervalAndCargo<TValue,TCargo>((*refIndelIt).originalPos + (*refIndelIt).indelSize,(*refIndelIt).originalPos + 1 - (*refIndelIt).indelSize, ii ));
				else
					appendValue(intervals, IntervalAndCargo<TValue,TCargo>((*refIndelIt).originalPos,(*refIndelIt).originalPos + 1, ii ));
			}
			++ii;
		}
		IntervalTree<TValue,TCargo> itree(intervals);
		if(options._debugLevel > 1 ) std::cout << "Interval tree done.\n" <<std::flush;
		
		// for each predicted variant check if it overlaps with a reference variant
		for(predIndelIt = currPredGenomeBegin; predIndelIt != currPredGenomeEnd; ++predIndelIt)
		{
			TIndel &predIndel = *predIndelIt;
			if(options._debugLevel > 1 ) std::cout << "." <<std::flush;

			// hack!!! needs to be fixed in snpStore, insertion position is one position to the left!
			if(predIndel.indelSize < 0) predIndel.originalPos += 1;
			
			// "core" interval
			TValue beginPoint = predIndel.originalPos;
			TValue endPoint = predIndel.originalPos;
			if(predIndel.indelSize > 0) endPoint += predIndel.indelSize;
			else endPoint += 1;
			
			// add tolerance in genomic position
			if(options.sequenceContext)
			{
				//compute eir, equivalent indel region, Krawitz et. al
				computeEir(predIndel,genomes[i],beginPoint,endPoint);
				beginPoint -= options.positionTolerance;
				endPoint += options.positionTolerance;
				
			}
			else
			{
				beginPoint -= options.positionTolerance;
				endPoint += options.positionTolerance;
			}
			
			// retrieve overlapping intervals
			String<TCargo> intersectingIntervals;
			findIntervals(itree,beginPoint,endPoint,intersectingIntervals);

			if(options._debugLevel > 1) std::cout <<"begin=" << beginPoint << " end=" << endPoint << std::endl;
			
			// check whether one of the overlapping intervals is a variant match
			bool predFound = false;
			CharString tagAppend = "";
			for(unsigned j = 0; j < length(intersectingIntervals); ++j)
			{
				if(options._debugLevel > 1 ) std::cout << "+" <<std::flush;
				TIndel &refIndel = refIndels[intersectingIntervals[j]];
				bool res = false;
				if(!options.sequenceContext) res = compareIndelPair(refIndel,predIndel,options);
				else res = compareIndelPair(refIndel,predIndel,beginPoint,endPoint,options);
				if(res)
				{
					if(options.annotateRepeats>0)
					{
						int countN = 0;
						for(TValue k = (TValue)_max((int)0,(int)refIndel.originalPos - (int)options.annotateRepeats); k < (TValue)refIndel.originalPos && k < (TValue)length(genomes[i]); k++) // count 'N's in upstream flanking sequence
							if(genomes[i][k] == 'N') ++countN;
						for(TValue k = refIndel.originalPos; k < _min((TValue)length(genomes[i]),(TValue)refIndel.originalPos+options.annotateRepeats); k++) // count 'N's in upstream flanking sequence
							if(genomes[i][k] == 'N') ++countN;
						int total = refIndel.originalPos - _max((int)0,(int)refIndel.originalPos - (int)options.annotateRepeats);
						total += _min((TValue)length(genomes[i]),(TValue)refIndel.originalPos+options.annotateRepeats) - refIndel.originalPos;
						std::stringstream tagAppendStr;
						tagAppendStr << ";percRep="  << int(100*countN/total);
						append(tagAppend,tagAppendStr.str());
							
					}
//					CharString tagAppend = ";matchID=";
					if(!predFound) append(tagAppend,tagAttach);
					else append(tagAppend,",");
					append(tagAppend,refIndel.idStr);
					predFound = true;
					refFoundSet.insert(intersectingIntervals[j]);
				}
				
			}
			if(predFound)
			{
				++TP;
                if(predIndel.type == INVERSION)
                    ++inversionsTP;
                else
    				for(unsigned j = 0; j < length(options.ranges); ++j)
	    				if(options.ranges[j].i1 <= predIndel.indelSize && predIndel.indelSize < options.ranges[j].i2)
		    				++rangeTP[j];
			}
			else
			{
				if(options._debugLevel>1)std::cout << "FP:" << predIndel.ninethCol << " begin=" << beginPoint << " end=" << endPoint << std::endl;
			}
			write(file,predIndel,genomeIDs[i],tagAppend,options);
			
		}
	
		if(options._debugLevel > 1 ) std::cout << std::endl;
	}
        // check number of false negatives

	for(foundSetIt = refFoundSet.begin(); foundSetIt != refFoundSet.end(); ++foundSetIt)
	{
        if(refIndels[*foundSetIt].type == INVERSION)
            ++inversionsFound;
        else   
            for(unsigned j = 0; j < length(options.ranges); ++j)
	    		if(options.ranges[j].i1 <= refIndels[*foundSetIt].indelSize && refIndels[*foundSetIt].indelSize < options.ranges[j].i2)
		    		++rangeFound[j];
                
	}

	// optionally append non-overlapped reference indels to output file
	file << "###################################################" << ::std::endl;
	file << "# Total Stats: " << ::std::endl;
	file << "# Number of reference indels: " << length(refIndels) << ::std::endl;
	file << "# Number of predicted indels: " << length(predIndels) << ::std::endl;
	file << "# Number of TP predictions  : " << TP << ::std::endl;
	file << "# Number of FP predictions  : " << length(predIndels)-TP << ::std::endl;
	file << "# Number of FN predictions  : " << length(refIndels)-refFoundSet.size() << ::std::endl;

	for(unsigned j = 0; j < length(options.ranges); ++j)
	{
		file << "###################################################" << ::std::endl;
		file << "# Bucket Stats [" << options.ranges[j].i1 << "," << options.ranges[j].i2 << "): " << ::std::endl;
		file << "# Number of reference indels: " << rangeCountRef[j] << ::std::endl;
		file << "# Number of predicted indels: " << rangeCountPred[j] << ::std::endl;
		file << "# Number of TP predictions  : " << rangeTP[j] << ::std::endl;
		file << "# Number of FP predictions  : " << rangeCountPred[j]-rangeTP[j] << ::std::endl;
		file << "# Number of FN predictions  : " << rangeCountRef[j]-rangeFound[j] << ::std::endl;
	}

    file << "###################################################" << ::std::endl;
	file << "# Inversion Stats: " << ::std::endl;
	file << "# Number of reference inversions: " << refInversions << ::std::endl;
	file << "# Number of predicted inversions: " << predInversions << ::std::endl;
	file << "# Number of TP predictions  : " << inversionsTP << ::std::endl;
	file << "# Number of FP predictions  : " << predInversions-inversionsTP << ::std::endl;
	file << "# Number of FN predictions  : " << refInversions-inversionsFound << ::std::endl;

	file.close();

	if(options._debugLevel > 0)
	{
		std::cout << "# Total Stats: " << ::std::endl;
		std::cout << "# Number of reference indels: " << length(refIndels) << ::std::endl;
		std::cout << "# Number of predicted indels: " << length(predIndels) << ::std::endl;
		std::cout << "# Number of TP predictions  : " << TP << ::std::endl;
		std::cout << "# Number of FP predictions  : " << length(predIndels)-TP << ::std::endl;
		std::cout << "# Number of FN predictions  : " << length(refIndels)-refFoundSet.size() << ::std::endl;
	}
	// FN = length(referenceIndels) - TP
	// FP = length(predictedIndels) - TP
	// output statistics at end of file

	
	return 0;
	
}




}

#endif
