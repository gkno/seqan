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
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include <seqan/sequence.h>
#include "razers.h"
#include "recognitionRateDP.h"
#include "readSimulator.h"











namespace SEQAN_NAMESPACE_MAIN
{
// ls in directory dir, store filenames in files
template<typename TPath, typename TFilenameString>
int 
getDir(TPath dir, TFilenameString &files)
{
    typedef typename Value<TFilenameString>::Type TFilename;
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(toCString(dir))) == NULL) {
        ::std::cout << "Error(" << errno << ") opening " << dir << ::std::endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
	TFilename name = (::std::string(dirp->d_name)).c_str();
        appendValue(files,name);
//	cout <<  files[length(files)-1] << " ?\n";
    }
    closedir(dp);
    return 0;
}

	struct ParamChooserOptions
	{
             typedef double TFloat;
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
            char 		best_shape_folder[200];
           bool 		best_shape_helpFolder;
            String<bool> firstTimeK;
            
            ParamChooserOptions()
            {
                minThreshold = 1;					// minimum value for threshold parameter 
                maxWeight = 13;                                           // maximum value of q
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
                best_shape_helpFolder = false;
            }
            
	// main options
	};




template<typename TValue>
inline TValue
_convertSolexaQual2ErrProb(TValue sq)
{
	return pow((TValue)10, sq / (TValue)-10) / ((TValue)1 + pow((TValue)10, sq / (TValue)-10));
}

template<typename TValue>
inline TValue
_convertPhredQual2ErrProb(TValue sq)
{
	return pow((TValue)10, sq / (TValue)-10);
}

template<typename TValue>
inline TValue
_convertSolexaQual2PhredQual(TValue sq)
{
	return (TValue)10 * log((TValue)1 + pow((TValue)10, sq / (TValue)10)) / log((TValue)10);
}



//////////////////////////////////////////////////////////////////////////////
// Get parameters q and t optimal for given loss rate
template<typename TFile, typename TSpec>
bool
parseParams(RazerSOptions<TSpec> & r_options, TFile & file, ParamChooserOptions & pm_options)
{
        typedef double TFloat;
        unsigned countQ = 0;
        unsigned countT = 0;
        TFloat bestSoFar = 0.0;
        unsigned bestQ = 1;
        unsigned bestT = 0;
	unsigned secondBestT = 0;
        //String<TFloat> loss;
        //reserve(loss, 20*readLen);

        char c = _streamGet(file);
	if(_streamEOF(file)) ::std::cout << "Loss rate file is empty!\n";

        while(!_streamEOF(file))
        {
                TFloat val = _parse_readEValue(file,c);
 //            ::std::cout << "\n"<<val;
//              appendValue(loss,val);  //t=0
                countT = 1;
                while(!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')))
                {
                        _parse_skipWhitespace(file,c);
                        val = _parse_readEValue(file,c);
  //                  ::std::cout << " " << val;
//                      appendValue(loss,val);
                        if( countT >= pm_options.minThreshold && val <= pm_options.optionLossRate /*&& val > bestSoFar*/)
                        {
                                bestSoFar=val;
                                bestQ=countQ+1;
                                bestT=countT;
				if(bestQ==13)
				{
					secondBestT=countT;
				}
                        }
                        ++countT;
                }
                ++countQ;
		if(countQ>=pm_options.maxWeight)
			break;
                _parse_skipWhitespace(file,c);
        }
	if(bestT<1)
        { ::std::cerr << "\n!!! Something wrong with file? !!!\n"; return false;}
	pm_options.chosenLossRate = bestSoFar;
        CharString chosenShape;
        fill(chosenShape, bestQ, '1');
        assign(r_options.shape,chosenShape);
        r_options.threshold = bestT;

        return true;
}
//compute average position dependent error distribution (assumes solexa qualtiy values in prb.txt format)
template<typename TFile, typename TDistribution>
void
qualityDistributionFromPrbFile(TFile & file, TDistribution & avg, ParamChooserOptions & pm_options)
{
        typedef typename Value<TDistribution>::Type TFloat;

	String<TFloat> qualitySum;
	String<int> count;
	fill(qualitySum,pm_options.totalN,0);
	fill(count,pm_options.totalN,0);

	if (_streamEOF(file)) return;

	char c = _streamGet(file);
	_parse_skipWhitespace(file, c);

	while (!_streamEOF(file))
	{
		for (unsigned pos = 0; (!_streamEOF(file)) && (pos < pm_options.totalN); ++pos)
		{
			_parse_skipBlanks(file,c);
			int qualA = (int) _parse_readDouble(file,c);
			_parse_skipBlanks(file,c);
			int qualC = (int) _parse_readDouble(file,c);
			_parse_skipBlanks(file,c);
			int qualG = (int) _parse_readDouble(file,c);
			_parse_skipBlanks(file,c);
			int qualT = (int) _parse_readDouble(file,c);
			int qual = ::std::max(::std::max(qualA, qualC), ::std::max(qualG, qualT));

//			::std::cout << qual << " ";
//			f = _convertSolexaQual2ErrProb(f);

			qualitySum[pos] += _convertSolexaQual2ErrProb((TFloat)qual);
			++count[pos];
		}
//		::std::cout << ::std::endl;
			
		_parse_skipLine2(file, c);
	}
	::std::cout << " Readcount = " << count[0] << "\n";

	fill(avg,pm_options.totalN,0.0);
	for(unsigned t = 0; t < pm_options.totalN; ++t)
	{
		TFloat f = (TFloat) qualitySum[t] / (TFloat)count[t];
//		f = _convertSolexaQual2ErrProb(f);
		avg[t] = f;
	}
}


template<typename TFile, typename TDistribution>
void
qualityDistributionFromFastQFile(TFile & file, TDistribution & avg, ParamChooserOptions & pm_options)
{
        typedef typename Value<TDistribution>::Type TFloat;

	String<int> qualitySum, count;
	fill(qualitySum,pm_options.totalN,0);
	fill(count,pm_options.totalN,0);

	if (_streamEOF(file)) return;

	signed char c = _streamGet(file);
	_parse_skipWhitespace(file, c);

	while (!_streamEOF(file))
	{
		_parse_skipLine2(file, c);
		if (_streamEOF(file) || c != '+') continue;

		_parse_skipLine2(file, c);

		unsigned i = 0;
		while (!(_streamEOF(file) || c == '\n' || c == '\r'))
		{
			qualitySum[i] += c - 33;
			c = _streamGet(file);
			++count[i];
			if (++i == pm_options.totalN) break;
		};
	}
	::std::cout << " Readcount = " << count[0] << "\n";

	fill(avg,pm_options.totalN,0.0);
	for(unsigned t = 0; t < pm_options.totalN; ++t)
	{
		TFloat f = (TFloat) qualitySum[t] / (TFloat)count[t];
		if (pm_options.solexaQual) f = _convertSolexaQual2PhredQual(f);
		f = _convertPhredQual2ErrProb(f);
		avg[t] = f;
	}
}

template<typename TFile, typename TDistribution>
void
qualityDistributionFromFastQIntFile(TFile & file, TDistribution & avg, ParamChooserOptions & pm_options)
{
        typedef typename Value<TDistribution>::Type TFloat;

	String<int> qualitySum, count;
	fill(qualitySum,pm_options.totalN,0);
	fill(count,pm_options.totalN,0);

	if (_streamEOF(file)) return;

	signed char c = _streamGet(file);
	_parse_skipWhitespace(file, c);

	while (!_streamEOF(file))
	{
		_parse_skipLine2(file, c);
		if (_streamEOF(file) || c != '+') continue;

		_parse_skipLine2(file, c);

		unsigned i = 0;
		while (!(_streamEOF(file) || c == '\n' || c == '\r'))
		{
			int num = _parse_readNumber(file, c);
			qualitySum[i] += num;
			++count[i];
			_parse_skipBlanks(file,c);
			if (++i == pm_options.totalN) break;
		};
	}
	::std::cout << " Readcount = " << count[0] << "\n";

	fill(avg,pm_options.totalN,0.0);
	for(unsigned t = 0; t < pm_options.totalN; ++t)
	{
		TFloat f = (TFloat) qualitySum[t] / (TFloat)count[t];
		if (pm_options.solexaQual)	f = _convertSolexaQual2PhredQual(f);
		f = _convertPhredQual2ErrProb(f);
		avg[t] = f;
	}
}


// find all *_prb.txt files in directory prbPath and compute average position dependent quality distribution
// compute average over all averages and store in errorDistribution
template<typename TPath, typename TError>
void
getAvgFromPrbDirectory(TPath prbPath, TError & errorDistribution, ParamChooserOptions & pm_options)
{
        typedef typename Value<TError>::Type TFloat;

	fill(errorDistribution,pm_options.totalN,0.0);
	
	String< ::std::string > files;
	getDir(prbPath,files);
	unsigned countPrbs = 0;
	for (unsigned int i = 0;i < length(files);i++) 
	{
		if(suffix(files[i],length(files[i])-6) == ".fastq")
		{
			::std::cout << "Processing "<< files[i] << "...\n";
			TError avg_act;
			resize(avg_act,pm_options.totalN);
			::std::fstream filestrm;
			::std::stringstream sstrm;
			sstrm << prbPath << files[i];
			filestrm.open(sstrm.str().c_str(),::std::ios_base::in);
			qualityDistributionFromFastQFile(filestrm,avg_act,pm_options);
			filestrm.close();
			for(unsigned j=0; j < pm_options.totalN; ++j)
			{
//				::std::cout << " " << avg_act[j];
				errorDistribution[j] += avg_act[j];
			}
			++countPrbs;
			continue;
		}
		if(suffix(files[i],length(files[i])-9) == ".fastqint")
		{
			::std::cout << "Processing "<< files[i] << "...\n";
			TError avg_act;
			resize(avg_act,pm_options.totalN);
			::std::fstream filestrm;
			::std::stringstream sstrm;
			sstrm << prbPath << files[i];
			filestrm.open(sstrm.str().c_str(),::std::ios_base::in);
			qualityDistributionFromFastQIntFile(filestrm,avg_act,pm_options);
			filestrm.close();
			for(unsigned j=0; j < pm_options.totalN; ++j)
			{
				::std::cout << " " << avg_act[j];
				errorDistribution[j] += avg_act[j];
			}
			++countPrbs;
			continue;
		}
		if(suffix(files[i],length(files[i])-8) == "_prb.txt")
		{
			::std::cout << "Processing "<< files[i] << "...\n";
			TError avg_act;
			resize(avg_act,pm_options.totalN);
			::std::fstream filestrm;
			::std::stringstream sstrm;
			sstrm << prbPath << files[i];
			filestrm.open(sstrm.str().c_str(),::std::ios_base::in);
			qualityDistributionFromPrbFile(filestrm,avg_act,pm_options);
			filestrm.close();
			for(unsigned j=0; j < pm_options.totalN; ++j)
			{
//				::std::cout << " " << avg_act[j];
				errorDistribution[j] += avg_act[j];
			}
			++countPrbs;
			continue;
		}
	}
	for(unsigned j=0; j < pm_options.totalN; ++j)
		errorDistribution[j] /= (TFloat)countPrbs;
	::std::cout << "Writing average error probabilities to " << pm_options.fprefix[0] << "_errorProb.dat\n";
	::std::fstream out;
	::std::stringstream avgOut;
	avgOut << pm_options.fprefix[0] << "_errorProb.dat";
	out.open(avgOut.str().c_str(),::std::ios_base::out);
	if(!out.is_open()) ::std::cout << "Couldn't write to file "<<avgOut.str()<<"\n";
	else
		for(unsigned j=0; j < pm_options.totalN; ++j)
			out << errorDistribution[j] << "\n";
	out.close();

	
}


template<typename TError>
void
makeUngappedStatsFile(TError & errorDistr, ParamChooserOptions & pm_options)
{
        typedef typename Value<TError>::Type TFloat;

	unsigned maxErrors = (unsigned) pm_options.totalN / 10;
	if(maxErrors<5)maxErrors=5;
	unsigned maxT = pm_options.totalN;
	
	// prepare log error distribution 
	String<TFloat> logErrorDistribution;
	resize(logErrorDistribution, 4*pm_options.totalN);

	// transformed probs for seeing 1s at positions 0...optionMaxN-1
	double remainingProb = 1.0 - pm_options.optionProbINSERT - pm_options.optionProbDELETE;
	for(unsigned j = 0; j < pm_options.totalN; ++j) 
	{
		logErrorDistribution[SEQAN_MISMATCH*pm_options.totalN+j] = _transform(errorDistr[j]);
		logErrorDistribution[SEQAN_INSERT*pm_options.totalN+j]   = _transform(pm_options.optionProbINSERT);
		logErrorDistribution[SEQAN_DELETE*pm_options.totalN+j]   = _transform(pm_options.optionProbDELETE);
		logErrorDistribution[SEQAN_MATCH*pm_options.totalN+j]    = _transform(remainingProb - errorDistr[j]);
	}
	
	CharString shape;
	for(int qLen = 8; qLen < 15; ++qLen)
	{
		clear(shape);
		fill(shape, qLen, '1');
		
		String<TFloat> found;
		resize(found,maxT*maxErrors);
		
		String< State<TFloat> > states;
		initPatterns(states, shape, maxErrors-1, logErrorDistribution, pm_options.optionHammingOnly);
		computeFilteringLoss(found, states, length(shape), maxT, maxErrors,  logErrorDistribution);
		
#ifdef LOSSRATE_VALIDATION	
	//for loss rate sanity check
		if(qLen > 7)
		{
			for(unsigned e=0; e < maxErrors; ++e)
			{
				::std::stringstream filename;
				filename << pm_options.fparams <<pm_options.fprefix[0]<<"_QE0_N" << pm_options.totalN<< "_E" << e<<".dat";
				::std::ofstream outfile;
				if(qLen==8){
					outfile.open(filename.str().c_str(),::std::ios_base::out);
					::std::cout << "Creating file "<<filename.str() << "\n";
				}
				else outfile.open(filename.str().c_str(),::std::ios_base::app);
				for(unsigned t = 20; t > 0; --t)
				{
					outfile.precision(10);
					//write(outfiles[e],1.0-exp(found[e*totalN+t]));
					//outfiles[e].write(1.0-exp(found[e*totalN+t]));
					outfile << (1.0 - _transformBack(found[e*maxT+t]));
					if(t>1)outfile << "\n";
				//      ::std::cout.precision(3);
				//      if (t > 0) ::std::cout << "\t";
				//      ::std::cout << (1.0-exp(found[e*totalN+t]));
				}
	//                     fprintf(outfiles[e],"%c",'\n');
	//      		               outfiles[e].write("\n");
				outfile << ::std::endl;
				outfile.close();
			}
		}
#endif
#ifndef LOSSRATE_VALIDATION	
		//regular output
		for(unsigned e=0; e < maxErrors; ++e)
		{
			::std::stringstream filename;
			filename << pm_options.fparams <<pm_options.fprefix[0]<<"_QE0_N" << pm_options.totalN<< "_E" << e<<".dat";
			::std::ofstream outfile;
			if(qLen==1){
				outfile.open(filename.str().c_str(),::std::ios_base::out);
				::std::cout << "Creating file "<<filename.str() << "\n";
			}
			else outfile.open(filename.str().c_str(),::std::ios_base::app);
			for(unsigned t = 0; t < maxT; ++t) 
			{
				if (t > 0) outfile << "\t";
				outfile.precision(8);
				//write(outfiles[e],1.0-_transformBack(found[e*maxT+t]));
				//outfiles[e].write(1.0-_transformBack(found[e*maxT+t]));
				outfile << (1.0-_transformBack(found[e*maxT+t]));
			//	::std::cout.precision(3);
			//	if (t > 0) ::std::cout << "\t";
			//	::std::cout << (1.0-_transformBack(found[e*maxT+t]));
			}
//			fprintf(outfiles[e],"%c",'\n');
//			outfiles[e].write("\n");
			outfile << ::std::endl;	
			outfile.close();
		}
//		::std::cout << ::std::endl;	
#endif
	}
// 	for(unsigned e=0; e < maxErrors; ++e)optionHammingOnly
// 	{
// //		outfiles[e].close();
// 		fclose(outfiles[e]);
// 	}

}

//////////////////////////////////////////////////////////////////////////////
// Returns the estimated minimum coverage of a shape with weight q, span s at threshold t
template<typename TValueQ, typename TValueS, typename TValueT>
inline TValueS getMinCov(TValueQ q, TValueS s, TValueT t)
{
	TValueS mincov;
	if(t > s - q + 1){
		mincov = q + 2 * (t - 1) - (t - (s - q + 1));
	}
	else mincov = q + 2 * (t - 1);

	return mincov;
}



template<typename TError>
void
makeSelectedStatsFile(TError & errorDistr, ParamChooserOptions & pm_options)
{

        typedef typename Value<TError>::Type TFloat;
	unsigned maxErrors = (unsigned) 1 + pm_options.totalN / 10;
	unsigned minErrors = 0;
        if(maxErrors<5 && pm_options.totalN > 30) maxErrors = 5;
	unsigned minQ = 7;
	unsigned maxT = pm_options.totalN-minQ+1;
	unsigned minT = 0;//totalN-minQ+1;

	typedef typename Value<TError>::Type TErrorValue;
	String<TErrorValue> logErrorDistribution;
	
	String<CharString> shapeStrings;

	if(pm_options.totalN < 32)
	{
	//q=6
	appendValue(shapeStrings,"111111");
	appendValue(shapeStrings,"1111100001");
	appendValue(shapeStrings,"11000000100100101");
	}
	
	if(pm_options.totalN < 36)
	{
	//q=7
	appendValue(shapeStrings,"1111111");
	appendValue(shapeStrings,"1111100011");
	appendValue(shapeStrings,"10110000001100101");
	}

	if(pm_options.totalN < 40)
	{
	//q=8
	appendValue(shapeStrings,"11111111");
	appendValue(shapeStrings,"11111100011");
	appendValue(shapeStrings,"101001111000101");  //median shape
	}
	
	if(pm_options.totalN < 50)
	{
	//q=9
	appendValue(shapeStrings,"111111111");
	appendValue(shapeStrings,"111111100011");
	appendValue(shapeStrings,"111001001010001011");
	}

	//q=10
	appendValue(shapeStrings,"1111111111");
	appendValue(shapeStrings,"1111111000111");
	appendValue(shapeStrings,"111001001010011101");

	//q=11
	appendValue(shapeStrings,"11111111111");
	appendValue(shapeStrings,"1111111001111");
	appendValue(shapeStrings,"11111101110101");  //median shape
	
	//q=12
	appendValue(shapeStrings,"111111111111");
	appendValue(shapeStrings,"11111111100111");
	appendValue(shapeStrings,"1110100111010011101");
	
	//q=13
 	appendValue(shapeStrings,"1111111111111");
 	appendValue(shapeStrings,"11111111110000111");
	appendValue(shapeStrings,"110101111001100010111");  //made this one up

	//q=14
 	appendValue(shapeStrings,"11111111111111");
	appendValue(shapeStrings,"1111111111100000111");
 	appendValue(shapeStrings,"11101110110001110110001");
 	appendValue(shapeStrings,"1111011010001110011011"); //all made up


	String<unsigned> weights;
	fill(weights,length(shapeStrings),0);
	for(unsigned i = 0; i < length(shapeStrings) ; ++i)
		for(unsigned pos = 0; pos < length(shapeStrings[i]) ; ++pos)
			if(shapeStrings[i][pos] == '1')
				++weights[i];
		
	// prepare log error distribution 
	resize(logErrorDistribution, 4*pm_options.totalN);
	// transformed probs for seeing 1s at positions 0...optionMaxN-1
	double remainingProb = 1.0 - pm_options.optionProbINSERT - pm_options.optionProbDELETE;
	for(unsigned j = 0; j < pm_options.totalN; ++j) 
	{
		logErrorDistribution[SEQAN_MISMATCH*pm_options.totalN+j] = _transform(errorDistr[j]);
		logErrorDistribution[SEQAN_INSERT*pm_options.totalN+j]   = _transform(pm_options.optionProbINSERT);
		logErrorDistribution[SEQAN_DELETE*pm_options.totalN+j]   = _transform(pm_options.optionProbDELETE);
		logErrorDistribution[SEQAN_MATCH*pm_options.totalN+j]    = _transform(remainingProb - errorDistr[j]);
	}

#ifdef RUN_RAZERS
	// generate genome and reads
	StringSet<Dna5String> testGenome;
	StringSet<Dna5String> testReads;
	StringSet<CharString> dummyIDs;
	resize(testGenome, 1);
	simulateGenome(testGenome[0], 1000000);					// generate 1Mbp genomic sequence
	simulateReads(
		testReads, dummyIDs, testGenome, 
		50000, maxErrors, logErrorDistribution, 0.5);	// generate 50K reads
#endif


	
	for(int i = length(shapeStrings)-1; i >= 0; --i)
	{
                if(length(shapeStrings[i])>pm_options.totalN) continue;
		String<TFloat> found;
		resize(found,maxT*maxErrors);
		
		String< State<TFloat> > states;
		if(pm_options.verbose)::std::cout << "do DP\n";
		initPatterns(states, shapeStrings[i], maxErrors-1, logErrorDistribution, pm_options.optionHammingOnly);
		computeFilteringLoss(found, states, length(shapeStrings[i]), maxT, maxErrors,  logErrorDistribution);
		
		for(unsigned e = minErrors; e < maxErrors; ++e) {
			bool highestOptimalFound = false;
			for(unsigned t = maxT-1; t > minT; --t) {
				TFloat lossrate = 1.0 - (TFloat) _transformBack(found[e*maxT+t]);
			//	typename ::std::map<TFloat, unsigned>::iterator it, itlow, itup, itend, itbegin;
				if(lossrate <= 0.0){
					if(highestOptimalFound) break;
					else highestOptimalFound = true;
				}
				if(lossrate > 0.2) continue;

				unsigned gminCov = getMinCov(weights[i], length(shapeStrings[i]), t);

				// create the whole file name
				::std::stringstream datName;
				if(pm_options.best_shape_helpFolder) datName << pm_options.best_shape_folder;
				else datName << "gapped_params";
				datName << "/"<<pm_options.fprefix[0]<<"_N" << pm_options.totalN << "_E" << e << "_";
				if(!pm_options.optionHammingOnly) datName << "L.dat";
				else datName <<"H.dat";
				
			
				// if datName-file doesnt exist, write the title on it
				if(pm_options.firstTimeK[e]==true){
					pm_options.firstTimeK[e] = false;
					::std::ofstream fout(datName.str().c_str(), ::std::ios::out);
					fout << "shape\t\tt\t\tlossrate\t\tminCoverage";
					fout << "\tPM\t\truntime";
					fout << ::std::endl << ::std::endl;
					fout.close();
				}
				
#ifdef RUN_RAZERS
				// count verifications
				String<ReadMatch<unsigned> > matches;
				RazerSOptions<RazerSSpec<false, true> > razersOptions;
				razersOptions.errorRate = (double)e / (double)pm_options.totalN;
				razersOptions.errorRate += 0.0000001;
				razersOptions.threshold = t;
				razersOptions._debugLevel = 2;
				razersOptions.hammingOnly = pm_options.optionHammingOnly;

				assign(razersOptions.shape, shapeStrings[i]);
				mapReads(matches, testGenome, testReads, razersOptions);
#endif

				// write best shape with its properties on the file
				::std::ofstream fout(datName.str().c_str(), ::std::ios::app | ::std::ios::out);
				fout << shapeStrings[i] << "\t\t";
				fout << t << "\t\t";
				fout << lossrate << "\t\t";
				fout << gminCov;
#ifdef RUN_RAZERS
				fout << "\t\t" << razersOptions.FP + razersOptions.TP;
				fout << "\t\t" << razersOptions.timeMapReads;
#else
				fout << "\t\t0\t\t0";
#endif
				fout << ::std::endl;
				fout.close();
				
			} // t-loop
		}

	
	}


}



template<typename TError>
void
makeOneGappedStatsFile(TError & errorDistr, ParamChooserOptions & pm_options)
{
        typedef typename Value<TError>::Type TFloat;

	unsigned maxE = (unsigned) pm_options.totalN / 10;
	if(maxE<5) maxE = 5;
	unsigned minE = 0;				
	unsigned maxQ = 14;				// weights are considered from minQ..maxQ
	unsigned minQ = 10;
	unsigned minGap = 0; 
	unsigned maxGap = 6;				// spans are considered from minQ..maxS
	unsigned maxT = pm_options.totalN-minQ+1;
	
	typedef typename Value<TError>::Type TErrorValue;
	String<TErrorValue> logErrorDistribution;


	resize(logErrorDistribution, 4*pm_options.totalN);

	// transformed probs for seeing 1s at positions 0...optionMaxN-1
	double remainingProb = 1.0 - pm_options.optionProbINSERT - pm_options.optionProbDELETE;
	for(unsigned j = 0; j < pm_options.totalN; ++j) 
	{
		logErrorDistribution[SEQAN_MISMATCH*pm_options.totalN+j] = _transform(errorDistr[j]);
		logErrorDistribution[SEQAN_INSERT*pm_options.totalN+j]   = _transform(pm_options.optionProbINSERT);
		logErrorDistribution[SEQAN_DELETE*pm_options.totalN+j]   = _transform(pm_options.optionProbDELETE);
		logErrorDistribution[SEQAN_MATCH*pm_options.totalN+j]    = _transform(remainingProb - errorDistr[j]);
	}
	
#ifdef RUN_RAZERS_ONEGAPPED
	// generate genome and reads
	StringSet<Dna5String> testGenome;
	StringSet<Dna5String> testReads;
	StringSet<CharString> dummyIDs;
	resize(testGenome, 1);
	simulateGenome(testGenome[0], 1000000);					// generate 1Mbp genomic sequence
	simulateReads(
		testReads, dummyIDs, testGenome, 
		50000, maxE, logErrorDistribution, 0.5);	// generate 10M reads
#endif


	String<TErrorValue> found;
	resize(found,maxT*maxE);

	// for each weight q...
	for(unsigned q = minQ; q <= maxQ; ++q){

		// j = span of shape
		for(unsigned j = q+minGap; j <= q+maxGap /*|| j < q+q-1 */; ++j){

                        if(j > pm_options.totalN) continue;
			// k = position of gap
//			for(unsigned k = (q/2); k < q; ++k){
			for(unsigned k = q-3; k < q-2; ++k){
				
				CharString shapeString;
				fill(shapeString,j,'0');
				for(unsigned pos = 0; pos < k; ++pos)
					shapeString[pos] = '1';
				for(unsigned pos = 0; pos < j-q; ++pos)
					shapeString[k+pos] = '0';
				for(unsigned pos = k+j-q; pos < j; ++pos)
					shapeString[pos] = '1';
				
				if(pm_options.verbose) ::std::cout << "doDP...\n";
				String< State<TFloat> > states;
				initPatterns(states, shapeString, maxE-1, logErrorDistribution, pm_options.optionHammingOnly, true);
				computeFilteringLoss(found, states, j, maxT, maxE, logErrorDistribution, false, true);

				// go through found and find loss rate
				for(unsigned e = minE; e < maxE; ++e) {
					bool highestOptimalFound = false;
					for(unsigned t = maxT-1; t > 0; --t) {
						TFloat lossrate = 1.0 - (TFloat) _transformBack(found[e*maxT+t]);
						
						if(lossrate <= 0.0){
							if(highestOptimalFound) break;
							else highestOptimalFound = true;
						}
						if(lossrate > 0.2)
							continue;
						
						unsigned gminCov = getMinCov(q, j, t);
	

						// create the whole file name
						::std::stringstream datName;
						if(pm_options.best_shape_helpFolder) datName << pm_options.best_shape_folder;
						else datName << "gapped_params";
						datName << "/"<<pm_options.fprefix[0]<<"_N" << pm_options.totalN << "_E" << e << "_";
						if(!pm_options.optionHammingOnly) datName << "L_";
						else datName <<"H_";
						datName << "onegapped.dat";
						//datName << q<<"_onegapped.dat";
					
						// if datName-file doesnt exist, write the title on it
						if(pm_options.firstTimeK[e]==true){
							pm_options.firstTimeK[e] = false;
							::std::ofstream fout(datName.str().c_str(), ::std::ios::out);
							fout << "shape\t\tt\t\tlossrate\t\tminCoverage";
							fout << "\t\tPM\t\truntime";
							fout << ::std::endl << ::std::endl;
							fout.close();
						}
						
				
#ifdef RUN_RAZERS_ONEGAPPED
						// count verifications
						String<ReadMatch<unsigned> > matches;
						RazerSOptions<RazerSSpec<false, true> > razersOptions;
						razersOptions.errorRate = (double)e / (double)pm_options.totalN;
						razersOptions.errorRate += 0.0000001;
						razersOptions.threshold = t;
						razersOptions._debugLevel = 2;
						razersOptions.hammingOnly = pm_options.optionHammingOnly;
		
						assign(razersOptions.shape, shapeString);
						mapReads(matches, testGenome, testReads, razersOptions);
#endif
						// write shape with its properties into file
						::std::ofstream fout(datName.str().c_str(), ::std::ios::app | ::std::ios::out);
						fout << shapeString << "\t\t";
						fout << t << "\t\t";
						fout << lossrate << "\t\t";
						fout << gminCov;
#ifdef RUN_RAZERS_ONEGAPPED				
						fout << "\t\t" << razersOptions.FP + razersOptions.TP;
						fout << "\t\t" << razersOptions.timeMapReads;
#else
						fout << "\t\t0\t\t0";
#endif
						fout << ::std::endl;
						fout.close();
									//}
						
					} // t-loop
				} //e-loop
					
			}// k-loop

		} // j-loop
				
	}// q-loop

}



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
        typedef double TFloat;
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
                if(!r_options.hammingOnly && numGaps(currShape)>0)
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
		if(pm_options.verbose) ::std::cerr << "\n!!! Something wrong with file? !!!" << ::std::endl;
		return false;
	}
	int i;
	for(i = pm_options.maxWeight-1; i >= 0; --i )
		if(length(shapes[i]) > 0)  // if a shape of weight i+1 has been found
			break;
	if(i<0)
	{
		if(pm_options.verbose) ::std::cerr << "\n!!! Something wrong with file? !!!" << ::std::endl;
		return false;
	}
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
        typedef double TFloat;
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
			::std::cerr << "\nNo session id given, using prefix 'userdef'"<<::std::endl;
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
			getAvgFromPrbDirectory(pm_options.fname[0],errorDistribution,pm_options);
		}

		::std::fstream file;
		//if(prefixCount)
		if(pm_options.doAllOneGapped) makeOneGappedStatsFile(errorDistribution,pm_options);
		if(pm_options.doSelectedGapped) makeSelectedStatsFile(errorDistribution,pm_options);
		if(pm_options.doUngapped) makeUngappedStatsFile(errorDistribution,pm_options);
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
        if (pm_options.verbose)
        {
               ::std::cerr << ::std::endl;
               ::std::cerr << "Read length      = " << pm_options.totalN << "bp\n";
               ::std::cerr << "Max num errors   = " << pm_options.totalK << "\n";
               ::std::cerr << "Recognition rate = " << 100.0*(1.0-pm_options.optionLossRate) << "%\n";
        }
	// parse loss rate file and find appropriate filter criterium
	if(pm_options.verbose) ::std::cerr << "\n--> Reading " <<  paramsfile.str()<<::std::endl;
	::std::fstream file;
	file.open(paramsfile.str().c_str(),::std::ios_base::in | ::std::ios_base::binary);
	if(!file.is_open())
	{
		if(pm_options.verbose)::std::cerr << "Couldn't open file "<<paramsfile.str()<<::std::endl;
		return false;
	}
	else
	{
		if(pm_options.doSelectedGapped || pm_options.doAllOneGapped) parseGappedParams(r_options,file,pm_options);
		else parseParams(r_options,file,pm_options);
		if(pm_options.verbose) ::std::cerr << "\n Choose \nshape: " << r_options.shape << "\n and \nthreshold: " << r_options.threshold<< "\n to achieve optimal performance for expected recognition rate >= " << (100.0-100.0*pm_options.optionLossRate) << "% (expected recognition = " << (100.0-pm_options.chosenLossRate*100.0) <<"%)\n\n";
		file.close();
	}

	return true;
}
}

#endif
