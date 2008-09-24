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

#ifndef SEQAN_HEADER_PARAMCHOOSER_H
#define SEQAN_HEADER_PARAMCHOOSER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <seqan/sequence.h>
#include "razers.h"

namespace SEQAN_NAMESPACE_MAIN
{

	struct ParamChooserOptions
	{
            typedef float TFloat;
            unsigned minThreshold;					// minimum value for threshold parameter 
            unsigned maxWeight;                                           // maximum value of q
            bool optionChooseOneGappedOnly;      // choose onegapped (or ungapped) shape (discard all other gapped shapes)
            
            
            // global input parameters
            unsigned totalN;				// sequence length
            unsigned totalK;					// errors
            TFloat optionLossRate;		// in
            TFloat chosenLossRate;			// out
            TFloat optionErrorRate;		// 
            bool optionHammingOnly;
            bool doUngapped;
            bool doAllOneGapped;
            bool doSelectedGapped;
            
            TFloat optionProbINSERT;
            TFloat optionProbDELETE;
            
            // output parameters
//            unsigned chosenQ;			//
//            unsigned chosenThreshold;		//
//            CharString chosenShape;			//
            
            ::std::string       paramFolderPath;
            bool		fnameCount0;
            bool		fnameCount1;
            bool		prefixCount;
            const char	*fname[2];
            const char	*fprefix[1];
            ::std::string		fparams;
            ::std::string		fgparams;
            bool		verbose;

            bool		solexaQual;
    //        char 		best_shape_folder[200];
      //      bool 		best_shape_helpFolder = false;
            String<bool> firstTimeK;
            
            ParamChooserOptions()
            {
                minThreshold = 1;					// minimum value for threshold parameter 
                maxWeight = 14;                                           // maximum value of q
                optionChooseOneGappedOnly = false;      // choose onegapped (or ungapped) shape (discard all other gapped shapes)
                
                
                // global input parameters
                totalN = 32;				// sequence length
                totalK = 2;					// errors
                optionLossRate = 0.01;		// in
                chosenLossRate = 0.0;			// out
                optionErrorRate = 0.05;		// 
                optionHammingOnly = false;
#ifdef LOSSRATE_VALIDATION
                doUngapped = true;
                doAllOneGapped = false;
                doSelectedGapped = false;
#else
                doUngapped = false;
                doAllOneGapped = false;
                doSelectedGapped = true;
#endif
                
                optionProbINSERT = 0.0;
                optionProbDELETE = 0.0;
                
                // output parameters
               // unsigned chosenQ = 4;			//
               // unsigned chosenThreshold = 2;		//
               // CharString chosenShape;			//
                
                paramFolderPath = "";
                fnameCount0 = 0;
                fnameCount1 = 0;
                prefixCount = 0;
                fname[0] = "";
                fname[1] = "";
                fprefix[0] =  "" ;
                verbose = true;
                solexaQual = true;
            //    best_shape_folder[200];
              //  best_shape_helpFolder = false;
            }
            
	// main options
	};



template<typename TSStr>
void
getParamsFilename(TSStr & paramsfile, ParamChooserOptions & pm_options)
{

	paramsfile.str("");
	if(pm_options.doSelectedGapped || pm_options.doAllOneGapped)
	{
		paramsfile << pm_options.fgparams<< pm_options.fprefix[0]<<"_N" << pm_options.totalN << "_E" << pm_options.totalK;
		//if(prefixCount) paramsfile << fgparams<< fprefix[0]<<"_N" << totalN << "_E" << totalK;
		//else paramsfile << fgparams<<"userdef_N" << totalN << "_E" << totalK;
		if(pm_options.optionHammingOnly) paramsfile << "_H";
		else paramsfile << "_L";
		if(pm_options.doAllOneGapped) paramsfile << "_onegapped.dat";
		else paramsfile << ".dat";
	}
	else
	{
		paramsfile << pm_options.fparams<< pm_options.fprefix[0]<<"_QE0_N" << pm_options.totalN << "_E" << pm_options.totalK << ".dat";
		//if(prefixCount) paramsfile << fparams<< fprefix[0]<<"_QE0_N" << totalN << "_E" << totalK << ".dat";
		//else paramsfile << fparams<<"userdef_QE0_N" << totalN << "_E" << totalK << ".dat";
	}
}



//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline void 
_parse_skipWhitespace(TFile& file, TChar& c)
{
	if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) break;
	}
}

template<typename TFile, typename TChar>
inline void 
_parse_skipBlanks(TFile& file, TChar& c)
{
	if ((c != ' ') && (c != '\t')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c != ' ') && (c != '\t')) break;
	}
}

template<typename TFile, typename TChar>
inline void 
_parse_skipLine2(TFile& file, TChar& c)
{
	if (c != '\n' && c != '\r')
		while (!_streamEOF(file)) {
			c = _streamGet(file);
			if (c == '\n' || c == '\r') break;
		}
	if (!_streamEOF(file))
		c = _streamGet(file);
}



//////////////////////////////////////////////////////////////////////////////
template<typename TFile, typename TChar>
inline double
_parse_readEValue(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
        typedef double TFloat;
	// Read number
	String<char> str(c);
	bool e = false;
	TFloat val1 = 0;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if(!e && c == 'e'){
			e = true;
			val1 = atof(toCString(str));
			c = _streamGet(file);
			resize(str,0);
		}
		if (!_parse_isDigit(c) && c != '.' && c != '-' && c != '+') break;
		append(str, c);
	}
	if(e)
	{
		return val1 * pow((TFloat)10.0,(TFloat)atof(toCString(str)));
	}	
 	else 
		return (TFloat)atof(toCString(str));
}



//////////////////////////////////////////////////////////////////////////////
template<typename TFile, typename TChar>
inline double
_parse_readDouble(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isDigit(c) && (c != '.')) break;
		append(str, c);
	}
 	return atof(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////
template<typename TChar>
inline bool
_parse_isDigit(TChar const c)
{
	return (c >= '0') && (c <= '9');
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline int
_parse_readNumber(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isDigit(c)) break;
		append(str, c);
	}
 	return atoi(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline void 
_parse_skipLine(TFile& file, TChar& c)
{
	if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
		c = _streamGet(file);
		return;
	}
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
	}
	c = _streamGet(file);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline void
_parse_readShape(TFile & file, TChar& c, CharString & str)
{
	// Read word
	append(str, c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!(c == '1' || c == '0')) break;
		append(str, c);
	}
}

//////////////////////////////////////////////////////////////////////////////


template<typename TShape>
inline int
numGaps(TShape & currShape)
{
    int count = 0;
    unsigned j=0;
    bool ingap = false;
    while(j<length(currShape))
    {
        if (currShape[j]=='0')
        {
            if(ingap) ++j;
            else ++count;
            ingap = true;
        }
        else ingap = false;
        ++j;
    }

    return count;

}


//////////////////////////////////////////////////////////////////////////////
// Get parameters q and t optimal for given loss rate
template<typename TFile, typename TSpec>
bool
parseGappedParams(RazerSOptions<TSpec> & r_options,TFile & file, ParamChooserOptions & pm_options)
{
        typedef float TFloat;
	String<CharString> shapes;
	resize(shapes,14); //best shape for each possible value of q
	String<unsigned> thresholds;
	resize(thresholds,14); //corresponding t
	String<unsigned> measure;
	resize(measure,14); //potential matches (or mincov if doAllOneGapped==true)
	String<TFloat> lossrates;
	resize(lossrates,14); //lossrates
	
	char c = _streamGet(file);
	if(_streamEOF(file))
        { 
            if(pm_options.verbose) ::std::cerr << "Loss rate file is empty!" << ::std::endl;
            return false;
        }
	else
	{
		_parse_skipLine(file,c);
		_parse_skipLine(file,c);
	}

	bool atLeastOneFound = false;
	while(!_streamEOF(file))
	{
		CharString currShape;
		_parse_readShape(file, c, currShape);
                if(pm_options.optionChooseOneGappedOnly && numGaps(currShape)>1)
                {
                    _parse_skipLine(file,c); 
                    continue;
                }
                _parse_skipWhitespace(file,c);
                unsigned currThreshold = _parse_readNumber(file,c);
        	_parse_skipWhitespace(file,c);
                TFloat currLossrate = _parse_readEValue(file,c);
		_parse_skipWhitespace(file,c);
		unsigned currMeasure = _parse_readNumber(file,c); //minCov in the case of oneGapped

#ifdef RUN_RAZERS
		if(!pm_options.doAllOneGapped) 
		{
			_parse_skipWhitespace(file,c);
			currMeasure = _parse_readNumber(file,c); //PM in the case of selectedGapped
		}
#endif
		if(currThreshold >= pm_options.minThreshold && currLossrate <= pm_options.optionLossRate /*&& val > bestSoFar*/)
		{
		
			unsigned weight = 0;
			for(unsigned pos = 0; pos < length(currShape) ; ++pos)
				if(currShape[pos] == '1')
					++weight;
			if(length(shapes[weight-1]) > 0)  // if this is not the first shape with weight weight
			{				  // compare currShape to the best one found so far
#ifndef RUN_RAZERS
				if(currMeasure >= measure[weight-1]) //if neither pm nor runtime available -> use mincov (approximation)
#else
				if((pm_options.doAllOneGapped && currMeasure >= measure[weight-1]) || (!pm_options.doAllOneGapped && currMeasure <= measure[weight-1]))
#endif
				{
					if(currMeasure == measure[weight-1])
					{
						bool undecided = false;
						//next measure: threshold
						if(thresholds[weight-1] > currThreshold) 
						{
							_parse_skipLine(file,c); 
							continue;
						}
						else if(thresholds[weight-1] == currThreshold) undecided = true;

						//if still undecided: next measure: span
						if(undecided && length(shapes[weight-1]) > length(currShape))
						{
							_parse_skipLine(file,c); 
							continue;
						}
						else if(undecided && length(shapes[weight-1]) < length(currShape)) undecided = false;

						//if still undecided: next measure: lossrate
						if(undecided && lossrates[weight-1] < currLossrate)
						{
							_parse_skipLine(file,c); 
							continue;
						}
					}
					shapes[weight-1] = currShape;
					measure[weight-1] = currMeasure;
					thresholds[weight-1] = currThreshold;
					lossrates[weight-1] = currLossrate;
					atLeastOneFound = true;
				}
				
			}
			else
			{
				shapes[weight-1] = currShape;
				measure[weight-1] = currMeasure;
				thresholds[weight-1] = currThreshold;
				lossrates[weight-1] = currLossrate;
				atLeastOneFound = true;
			
			}
                }
		_parse_skipLine(file,c);

        }
	if(!atLeastOneFound)
	{
		if(pm_options.verbose) ::std::cout << "\n!!! Something wrong with file? !!!" << ::std::endl;
		return false;
	}
	int i;
	for(i = pm_options.maxWeight-1; i >= 0; --i )
		if(length(shapes[i]) > 0)  // if a shape of weight i+1 has been found
			break;
	pm_options.chosenLossRate = lossrates[i];
	assign(r_options.shape, shapes[i]);
	r_options.threshold = thresholds[i];
	// suggest a suitable combination of q and t

        return true;
//        return false;
}



template<typename TSpec>
bool
chooseParams(RazerSOptions<TSpec> & r_options, ParamChooserOptions & pm_options)
{
        typedef float TFloat;
	static const TFloat epsilon = 0.0000001;	
	pm_options.optionLossRate += epsilon;

#ifdef LOSSRATE_VALIDATION	
	pm_options.fparams = pm_options.paramFolderPath;
	//pm_options.fparams = 
//	strcat(fparams,"");
#else
	pm_options.fparams = pm_options.paramFolderPath + "params/";
#endif
	pm_options.fgparams = pm_options.paramFolderPath + "gapped_params/";

	
	if(pm_options.optionProbINSERT <= epsilon && pm_options.optionProbDELETE <= epsilon)
		pm_options.optionHammingOnly=true;

	fill(pm_options.firstTimeK,20,true);//set maximal number of errors considered in parameter computation to <10

// compute data specific loss rates
	if (pm_options.fnameCount0 || pm_options.fnameCount1) 
	{
		if(!pm_options.prefixCount)
		{
			pm_options.fprefix[0] = "userdef";
			::std::cout << "\nNo session id given, using prefix 'userdef'"<<::std::endl;
		}
		String<TFloat> errorDistribution;
		resize(errorDistribution,pm_options.totalN);
		//error distribution given --> read file containing error distr and compute loss rates
		if(pm_options.fnameCount1)
		{
			::std::fstream file;
			file.open(pm_options.fname[1],::std::ios_base::in | ::std::ios_base::binary);
			if(!file.is_open())
			{
				::std::cerr << "Couldn't open file "<<pm_options.fname[1]<<"\n";
				return false;
			}
			unsigned count = 0;
			char c = _streamGet(file);
			while(!_streamEOF(file) && count < pm_options.totalN)
			{
				_parse_skipWhitespace(file,c);
				errorDistribution[count] = _parse_readEValue(file,c);// + (TFloat) 1.0/maxN;
				++count;
			}
			file.close();
			if(count != pm_options.totalN)
			{
				::std::cerr << "Error distribution file must contain at least " << pm_options.totalN << " probability values (one value per line).\n";
				return false;
			}
		}
		else // read qualtiy files and compute position dependent avg error probabilites
		{
//			getAvgFromPrbDirectory(pm_options.fname[0],errorDistribution);
		}

		::std::fstream file;
		//if(prefixCount)
//		if(pm_options.doAllOneGapped) makeOneGappedStatsFile(errorDistribution);
//		if(pm_options.doSelectedGapped) makeSelectedStatsFile(errorDistribution);
//		if(pm_options.doUngapped) makeUngappedStatsFile(errorDistribution);
	}
	else if(!pm_options.prefixCount) pm_options.fprefix[0] = "results";

	pm_options.totalK = (int)(pm_options.optionErrorRate * pm_options.totalN);
	
	//prioritize
	if(pm_options.doSelectedGapped)
	{
		pm_options.doAllOneGapped = false;
		pm_options.doUngapped = false;
	}
	if(pm_options.doAllOneGapped) pm_options.doUngapped = false;
	
	// decide on which loss rate file to parse
	::std::stringstream paramsfile;
	getParamsFilename(paramsfile,pm_options);
        if(pm_options.verbose)
        {
	   ::std::cout << "\nRead length      = " << pm_options.totalN << "bp\n";
	   ::std::cout << "Max num errors   = " << pm_options.totalK << "\n";
	   ::std::cout << "Recognition rate = " <<  100.0*(1.0-pm_options.optionLossRate) << "%\n";
	}		
	// parse loss rate file and find appropriate filter criterium
	if(pm_options.verbose) ::std::cout << "\n--> Reading " <<  paramsfile.str()<<"\n";
	::std::fstream file;
	file.open(paramsfile.str().c_str(),::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
	{
		if(pm_options.verbose)::std::cout << "Couldn't open file "<<paramsfile.str()<<"\n";
		return false;
	}
	else
	{
		/*if(pm_options.doSelectedGapped || pm_options.doAllOneGapped)*/ parseGappedParams(r_options,file,pm_options);
	//	else parseParams(file);
		if(pm_options.verbose) ::std::cout << "\n Choose \nshape: " << r_options.shape << "\n and \nthreshold: " << r_options.threshold<< "\n to achieve optimal performance for expected recognition rate >= " << (100.0-100.0*pm_options.optionLossRate) << "% (expected recognition = " << (100.0-pm_options.chosenLossRate*100.0) <<"%)\n\n";
		file.close();
	}

	return true;
}
}

#endif
