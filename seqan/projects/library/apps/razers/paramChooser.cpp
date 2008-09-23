#define USE_LOGVALUES		// this is recommended when using probability values
#define RUN_RAZERS
#define RUN_RAZERS_ONEGAPPED	
#define SEQAN_PROFILE
//#define LOSSRATE_VALIDATION	//generates output for loss rate validation (empirical vs computed), only for ungapped

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find_motif.h>

#include "recognitionRateDP.h"
#include "readSimulator.h"
#include "razers.h"


using namespace seqan;
using namespace std;

typedef long double TFloat;

static unsigned minThreshold = 1;					// minimum value for threshold parameter 
static unsigned maxWeight = 14;                                           // maximum value of q
static bool optionChooseOneGappedOnly = false;      // choose onegapped (or ungapped) shape (discard all other gapped shapes)


// global input parameters
static unsigned totalN = 32;				// sequence length
static unsigned totalK = 2;					// errors
static TFloat optionLossRate = 0.01;		// in
static TFloat chosenLossRate = 0.0;			// out
static TFloat optionErrorRate = 0.05;		// 
static bool optionHammingOnly = false;
#ifdef LOSSRATE_VALIDATION
static bool doUngapped = true;
static bool doAllOneGapped = false;
static bool doSelectedGapped = false;
#else
static bool doUngapped = false;
static bool doAllOneGapped = false;
static bool doSelectedGapped = true;
#endif

static TFloat optionProbINSERT = 0.0;
static TFloat optionProbDELETE = 0.0;

// output parameters
unsigned chosenQ = 4;			//
unsigned chosenThreshold = 2;		//
CharString chosenShape;			//

bool		fnameCount0 = 0;
bool		fnameCount1 = 0;
bool		prefixCount = 0;
const char	*fname[2] = { "","" };
const char	*fprefix[1] = { "" };
char		fparams[100];
char		fgparams[100];
bool		verbose = true;
bool		solexaQual = true;
char 		best_shape_folder[200];
bool 		best_shape_helpFolder = false;
String<bool> firstTimeK;

//////////////////////////////////////////////////////////////////////////////

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
inline TFloat
_parse_readEValue(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
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
inline TFloat
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

// ls in directory dir, store filenames in files
template<typename TPath, typename TFilenameString>
int 
getDir(TPath dir, TFilenameString &files)
{
    typedef typename Value<TFilenameString>::Type TFilename;
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(toCString(dir))) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
	TFilename name = (string(dirp->d_name)).c_str();
        appendValue(files,name);
//	cout <<  files[length(files)-1] << " ?\n";
    }
    closedir(dp);
    return 0;
}


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
template<typename TFile>
bool
parseGappedParams(TFile & file)
{
	String<CharString> shapes;
	resize(shapes,14); //best shape for each possible value of q
	String<unsigned> thresholds;
	resize(thresholds,14); //corresponding t
	String<unsigned> measure;
	resize(measure,14); //potential matches (or mincov if doAllOneGapped==true)
	String<TFloat> lossrates;
	resize(lossrates,14); //lossrates
	
	char c = _streamGet(file);
	if(_streamEOF(file)) cout << "Loss rate file is empty!\n";
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
                if(optionChooseOneGappedOnly && numGaps(currShape)>1)
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
		if(!doAllOneGapped) 
		{
			_parse_skipWhitespace(file,c);
			currMeasure = _parse_readNumber(file,c); //PM in the case of selectedGapped
		}
#endif
		if(currThreshold >= minThreshold && currLossrate <= optionLossRate /*&& val > bestSoFar*/)
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
				if((doAllOneGapped && currMeasure >= measure[weight-1]) || (!doAllOneGapped && currMeasure <= measure[weight-1]))
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
		if(verbose) cout << "\n!!! Something wrong with file? !!!\n";
		return false;
	}
	int i;
	for(i = maxWeight-1; i >= 0; --i )
		if(length(shapes[i]) > 0)  // if a shape of weight i+1 has been found
			break;
	chosenLossRate = lossrates[i];
	chosenShape = shapes[i];
	chosenThreshold = thresholds[i];
	// suggest a suitable combination of q and t

        return true;
//        return false;
}


//////////////////////////////////////////////////////////////////////////////
// Get parameters q and t optimal for given loss rate
template<typename TFile>
bool
parseParams(TFile & file)
{
        unsigned countQ = 0;
        unsigned countT = 0;
        TFloat bestSoFar = 0.0;
        unsigned bestQ = 1;
        unsigned bestT = 0;
	unsigned secondBestT = 0;
        //String<TFloat> loss;
        //reserve(loss, 20*readLen);

        char c = _streamGet(file);
	if(_streamEOF(file)) cout << "Loss rate file is empty!\n";

        while(!_streamEOF(file))
        {
                TFloat val = _parse_readEValue(file,c);
 //            cout << "\n"<<val;
//              appendValue(loss,val);  //t=0
                countT = 1;
                while(!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')))
                {
                        _parse_skipWhitespace(file,c);
                        val = _parse_readEValue(file,c);
  //                  cout << " " << val;
//                      appendValue(loss,val);
                        if( countT >= minThreshold && val <= optionLossRate /*&& val > bestSoFar*/)
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
		if(countQ>=maxWeight)
			break;
                _parse_skipWhitespace(file,c);
        }
	if(bestT<1) cout << "\n!!! Something wrong with file? !!!\n";
	chosenLossRate = bestSoFar;
        fill(chosenShape, bestQ, '1');
        chosenThreshold = bestT;

        return true;
}
//compute average position dependent error distribution (assumes solexa qualtiy values in prb.txt format)
template<typename TFile, typename TDistribution>
void
qualityDistributionFromPrbFile(TFile & file, TDistribution & avg)
{
	String<TFloat> qualitySum;
	String<int> count;
	fill(qualitySum,totalN,0);
	fill(count,totalN,0);

	if (_streamEOF(file)) return;

	char c = _streamGet(file);
	_parse_skipWhitespace(file, c);

	while (!_streamEOF(file))
	{
		for (unsigned pos = 0; (!_streamEOF(file)) && (pos < totalN); ++pos)
		{
			_parse_skipBlanks(file,c);
			int qualA = (int) _parse_readDouble(file,c);
			_parse_skipBlanks(file,c);
			int qualC = (int) _parse_readDouble(file,c);
			_parse_skipBlanks(file,c);
			int qualG = (int) _parse_readDouble(file,c);
			_parse_skipBlanks(file,c);
			int qualT = (int) _parse_readDouble(file,c);
			int qual = max(max(qualA, qualC), max(qualG, qualT));

//			cout << qual << " ";
//			f = _convertSolexaQual2ErrProb(f);

			qualitySum[pos] += _convertSolexaQual2ErrProb((TFloat)qual);
			++count[pos];
		}
//		cout << endl;
			
		_parse_skipLine2(file, c);
	}
	cout << " Readcount = " << count[0] << "\n";

	fill(avg,totalN,0.0);
	for(unsigned t = 0; t < totalN; ++t)
	{
		TFloat f = (TFloat) qualitySum[t] / (TFloat)count[t];
//		f = _convertSolexaQual2ErrProb(f);
		avg[t] = f;
	}
}


template<typename TFile, typename TDistribution>
void
qualityDistributionFromFastQFile(TFile & file, TDistribution & avg)
{
	String<int> qualitySum, count;
	fill(qualitySum,totalN,0);
	fill(count,totalN,0);

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
			if (++i == totalN) break;
		};
	}
	cout << " Readcount = " << count[0] << "\n";

	fill(avg,totalN,0.0);
	for(unsigned t = 0; t < totalN; ++t)
	{
		TFloat f = (TFloat) qualitySum[t] / (TFloat)count[t];
		if (solexaQual)	f = _convertSolexaQual2PhredQual(f);
		f = _convertPhredQual2ErrProb(f);
		avg[t] = f;
	}
}

template<typename TFile, typename TDistribution>
void
qualityDistributionFromFastQIntFile(TFile & file, TDistribution & avg)
{
	String<int> qualitySum, count;
	fill(qualitySum,totalN,0);
	fill(count,totalN,0);

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
			if (++i == totalN) break;
		};
	}
	cout << " Readcount = " << count[0] << "\n";

	fill(avg,totalN,0.0);
	for(unsigned t = 0; t < totalN; ++t)
	{
		TFloat f = (TFloat) qualitySum[t] / (TFloat)count[t];
		if (solexaQual)	f = _convertSolexaQual2PhredQual(f);
		f = _convertPhredQual2ErrProb(f);
		avg[t] = f;
	}
}


// find all *_prb.txt files in directory prbPath and compute average position dependent quality distribution
// compute average over all averages and store in errorDistribution
template<typename TPath, typename TError>
void
getAvgFromPrbDirectory(TPath prbPath, TError & errorDistribution)
{
	fill(errorDistribution,totalN,0.0);
	
	String<string> files;
	getDir(prbPath,files);
	unsigned countPrbs = 0;
	for (unsigned int i = 0;i < length(files);i++) 
	{
		if(suffix(files[i],length(files[i])-6) == ".fastq")
		{
			cout << "Processing "<< files[i] << "...\n";
			TError avg_act;
			resize(avg_act,totalN);
			fstream filestrm;
			stringstream sstrm;
			sstrm << prbPath << files[i];
			filestrm.open(sstrm.str().c_str(),ios_base::in);
			qualityDistributionFromFastQFile(filestrm,avg_act);
			filestrm.close();
			for(unsigned j=0; j < totalN; ++j)
			{
//				cout << " " << avg_act[j];
				errorDistribution[j] += avg_act[j];
			}
			++countPrbs;
			continue;
		}
		if(suffix(files[i],length(files[i])-9) == ".fastqint")
		{
			cout << "Processing "<< files[i] << "...\n";
			TError avg_act;
			resize(avg_act,totalN);
			fstream filestrm;
			stringstream sstrm;
			sstrm << prbPath << files[i];
			filestrm.open(sstrm.str().c_str(),ios_base::in);
			qualityDistributionFromFastQIntFile(filestrm,avg_act);
			filestrm.close();
			for(unsigned j=0; j < totalN; ++j)
			{
				cout << " " << avg_act[j];
				errorDistribution[j] += avg_act[j];
			}
			++countPrbs;
			continue;
		}
		if(suffix(files[i],length(files[i])-8) == "_prb.txt")
		{
			cout << "Processing "<< files[i] << "...\n";
			TError avg_act;
			resize(avg_act,totalN);
			fstream filestrm;
			stringstream sstrm;
			sstrm << prbPath << files[i];
			filestrm.open(sstrm.str().c_str(),ios_base::in);
			qualityDistributionFromPrbFile(filestrm,avg_act);
			filestrm.close();
			for(unsigned j=0; j < totalN; ++j)
			{
//				cout << " " << avg_act[j];
				errorDistribution[j] += avg_act[j];
			}
			++countPrbs;
			continue;
		}
	}
	for(unsigned j=0; j < totalN; ++j)
		errorDistribution[j] /= (TFloat)countPrbs;
	cout << "Writing average error probabilities to " << fprefix[0] << "_errorProb.dat\n";
	fstream out;
	stringstream avgOut;
	avgOut << fprefix[0] << "_errorProb.dat";
	out.open(avgOut.str().c_str(),ios_base::out);
	if(!out.is_open()) cout << "Couldn't write to file "<<avgOut.str()<<"\n";
	else
		for(unsigned j=0; j < totalN; ++j)
			out << errorDistribution[j] << "\n";
	out.close();

	
}


template<typename TError>
void
makeUngappedStatsFile(TError & errorDistr)
{
	unsigned maxErrors = (unsigned) totalN / 10;
	if(maxErrors<5)maxErrors=5;
	unsigned maxT = totalN;
	
	// prepare log error distribution 
	String<TFloat> logErrorDistribution;
	resize(logErrorDistribution, 4*totalN);

	// transformed probs for seeing 1s at positions 0...optionMaxN-1
	double remainingProb = 1.0 - optionProbINSERT - optionProbDELETE;
	for(unsigned j = 0; j < totalN; ++j) 
	{
		logErrorDistribution[SEQAN_MISMATCH*totalN+j] = _transform(errorDistr[j]);
		logErrorDistribution[SEQAN_INSERT*totalN+j]   = _transform(optionProbINSERT);
		logErrorDistribution[SEQAN_DELETE*totalN+j]   = _transform(optionProbDELETE);
		logErrorDistribution[SEQAN_MATCH*totalN+j]    = _transform(remainingProb - errorDistr[j]);
	}
	
	CharString shape;
	for(int qLen = 8; qLen < 15; ++qLen)
	{
		clear(shape);
		fill(shape, qLen, '1');
		
		String<TFloat> found;
		resize(found,maxT*maxErrors);
		
		String< State<TFloat> > states;
		initPatterns(states, shape, maxErrors-1, logErrorDistribution, optionHammingOnly);
		computeFilteringLoss(found, states, length(shape), maxT, maxErrors,  logErrorDistribution);
		
#ifdef LOSSRATE_VALIDATION	
	//for loss rate sanity check
		if(qLen > 7)
		{
			for(unsigned e=0; e < maxErrors; ++e)
			{
				stringstream filename;
				filename << fparams <<fprefix[0]<<"_QE0_N" << totalN<< "_E" << e<<".dat";
				ofstream outfile;
				if(qLen==8){
					outfile.open(filename.str().c_str(),ios_base::out);
					cout << "Creating file "<<filename.str() << "\n";
				}
				else outfile.open(filename.str().c_str(),ios_base::app);
				for(unsigned t = 20; t > 0; --t)
				{
					outfile.precision(10);
					//write(outfiles[e],1.0-exp(found[e*totalN+t]));
					//outfiles[e].write(1.0-exp(found[e*totalN+t]));
					outfile << (1.0 - _transformBack(found[e*maxT+t]));
					if(t>1)outfile << "\n";
				//      cout.precision(3);
				//      if (t > 0) cout << "\t";
				//      cout << (1.0-exp(found[e*totalN+t]));
				}
	//                     fprintf(outfiles[e],"%c",'\n');
	//      		               outfiles[e].write("\n");
				outfile << endl;
				outfile.close();
			}
		}
#endif
#ifndef LOSSRATE_VALIDATION	
		//regular output
		for(unsigned e=0; e < maxErrors; ++e)
		{
			stringstream filename;
			filename << fparams <<fprefix[0]<<"_QE0_N" << totalN<< "_E" << e<<".dat";
			ofstream outfile;
			if(qLen==1){
				outfile.open(filename.str().c_str(),ios_base::out);
				cout << "Creating file "<<filename.str() << "\n";
			}
			else outfile.open(filename.str().c_str(),ios_base::app);
			for(unsigned t = 0; t < maxT; ++t) 
			{
				if (t > 0) outfile << "\t";
				outfile.precision(8);
				//write(outfiles[e],1.0-_transformBack(found[e*maxT+t]));
				//outfiles[e].write(1.0-_transformBack(found[e*maxT+t]));
				outfile << (1.0-_transformBack(found[e*maxT+t]));
			//	cout.precision(3);
			//	if (t > 0) cout << "\t";
			//	cout << (1.0-_transformBack(found[e*maxT+t]));
			}
//			fprintf(outfiles[e],"%c",'\n');
//			outfiles[e].write("\n");
			outfile << endl;	
			outfile.close();
		}
//		cout << endl;	
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
makeSelectedStatsFile(TError & errorDistr)
{

	unsigned maxErrors = (unsigned) totalN / 10;
	if(maxErrors<5 && totalN > 30) maxErrors = 5;
	unsigned minQ = 7;
	unsigned maxT = totalN-minQ+1;
	unsigned minT = 0;//totalN-minQ+1;

	typedef typename Value<TError>::Type TErrorValue;
	String<TErrorValue> logErrorDistribution;
	
	String<CharString> shapeStrings;

	if(totalN < 32)
	{
	//q=6
	appendValue(shapeStrings,"111111");
	appendValue(shapeStrings,"1111100001");
	appendValue(shapeStrings,"11000000100100101");
	}
	
	if(totalN < 36)
	{
	//q=7
	appendValue(shapeStrings,"1111111");
	appendValue(shapeStrings,"1111100011");
	appendValue(shapeStrings,"10110000001100101");
	}

	if(totalN < 40)
	{
	//q=8
	appendValue(shapeStrings,"11111111");
	appendValue(shapeStrings,"11111100011");
	appendValue(shapeStrings,"101001111000101");  //median shape
	}
	
	if(totalN < 50)
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
	resize(logErrorDistribution, 4*totalN);
	// transformed probs for seeing 1s at positions 0...optionMaxN-1
	double remainingProb = 1.0 - optionProbINSERT - optionProbDELETE;
	for(unsigned j = 0; j < totalN; ++j) 
	{
		logErrorDistribution[SEQAN_MISMATCH*totalN+j] = _transform(errorDistr[j]);
		logErrorDistribution[SEQAN_INSERT*totalN+j]   = _transform(optionProbINSERT);
		logErrorDistribution[SEQAN_DELETE*totalN+j]   = _transform(optionProbDELETE);
		logErrorDistribution[SEQAN_MATCH*totalN+j]    = _transform(remainingProb - errorDistr[j]);
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
		50000, maxErrors-1, logErrorDistribution, 0.5);	// generate 50K reads
#endif


	
	for(int i = length(shapeStrings)-1; i >= 0; --i)
	{
		String<TFloat> found;
		resize(found,maxT*maxErrors);
		
		String< State<TFloat> > states;
		if(verbose)cout << "do DP\n";
		initPatterns(states, shapeStrings[i], maxErrors-1, logErrorDistribution, optionHammingOnly);
		computeFilteringLoss(found, states, length(shapeStrings[i]), maxT, maxErrors,  logErrorDistribution);
		
		for(unsigned e = 1; e < maxErrors; ++e) {
			bool highestOptimalFound = false;
			for(unsigned t = maxT-1; t > minT; --t) {
				TFloat lossrate = 1.0 - (TFloat) _transformBack(found[e*maxT+t]);
				typename map<TFloat, unsigned>::iterator it, itlow, itup, itend, itbegin;
				if(lossrate <= 0.0){
					if(highestOptimalFound) break;
					else highestOptimalFound = true;
				}
				if(lossrate > 0.2) continue;

				unsigned gminCov = getMinCov(weights[i], length(shapeStrings[i]), t);

				// create the whole file name
				stringstream datName;
				if(best_shape_helpFolder) datName << best_shape_folder;
				else datName << "gapped_params";
				datName << "/"<<fprefix[0]<<"_N" << totalN << "_E" << e << "_";
				if(!optionHammingOnly) datName << "L.dat";
				else datName <<"H.dat";
				
			
				// if datName-file doesnt exist, write the title on it
				if(firstTimeK[e]==true){
					firstTimeK[e] = false;
					ofstream fout(datName.str().c_str(), ios::out);
					fout << "shape\t\tt\t\tlossrate\t\tminCoverage";
					fout << "\tPM\t\truntime";
					fout << endl << endl;
					fout.close();
				}
				
#ifdef RUN_RAZERS
				// count verifications
				String<ReadMatch<unsigned> > matches;
				RazerSOptions<RazerSSpec<false, true> > razersOptions;
				razersOptions.errorRate = (double)e / (double)totalN;
				razersOptions.errorRate += 0.0000001;
				razersOptions.threshold = t;
				razersOptions._debugLevel = 2;

				assign(razersOptions.shape, shapeStrings[i]);
				mapReads(matches, testGenome, testReads, razersOptions);
#endif

				// write best shape with its properties on the file
				ofstream fout(datName.str().c_str(), ios::app | ios::out);
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
				fout << endl;
				fout.close();
				
			} // t-loop
		}

	
	}


}



template<typename TError>
void
makeOneGappedStatsFile(TError & errorDistr)
{

	unsigned maxE = (unsigned) totalN / 10;
	if(maxE<5) maxE = 5;
	unsigned minE = 0;				
	unsigned maxQ = 14;				// weights are considered from minQ..maxQ
	unsigned minQ = 10;
	unsigned minGap = 0; 
	unsigned maxGap = 6;				// spans are considered from minQ..maxS
	unsigned maxT = totalN-minQ+1;
	
	typedef typename Value<TError>::Type TErrorValue;
	String<TErrorValue> logErrorDistribution;


	resize(logErrorDistribution, 4*totalN);

	// transformed probs for seeing 1s at positions 0...optionMaxN-1
	double remainingProb = 1.0 - optionProbINSERT - optionProbDELETE;
	for(unsigned j = 0; j < totalN; ++j) 
	{
		logErrorDistribution[SEQAN_MISMATCH*totalN+j] = _transform(errorDistr[j]);
		logErrorDistribution[SEQAN_INSERT*totalN+j]   = _transform(optionProbINSERT);
		logErrorDistribution[SEQAN_DELETE*totalN+j]   = _transform(optionProbDELETE);
		logErrorDistribution[SEQAN_MATCH*totalN+j]    = _transform(remainingProb - errorDistr[j]);
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
		50000, maxE-1, logErrorDistribution, 0.5);	// generate 10M reads
#endif


	String<TErrorValue> found;
	resize(found,maxT*maxE);

	// for each weight q...
	for(unsigned q = minQ; q <= maxQ; ++q){

		// j = span of shape
		for(unsigned j = q+minGap; j <= q+maxGap /*|| j < q+q-1 */; ++j){

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
				
				if(verbose) cout << "doDP...\n";
				String< State<long double> > states;
				initPatterns(states, shapeString, maxE-1, logErrorDistribution, optionHammingOnly, true);
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
						//if(lossrate > 0.1)
							continue;
						
						unsigned gminCov = getMinCov(q, j, t);
	

						// create the whole file name
						stringstream datName;
						if(best_shape_helpFolder) datName << best_shape_folder;
						else datName << "gapped_params";
						datName << "/"<<fprefix[0]<<"_N" << totalN << "_E" << e << "_";
						if(!optionHammingOnly) datName << "L_";
						else datName <<"H_";
						datName << "onegapped.dat";
						//datName << q<<"_onegapped.dat";
					
						// if datName-file doesnt exist, write the title on it
						if(firstTimeK[e]==true){
							firstTimeK[e] = false;
							ofstream fout(datName.str().c_str(), ios::out);
							fout << "shape\t\tt\t\tlossrate\t\tminCoverage";
							fout << "\t\tPM\t\truntime";
							fout << endl << endl;
							fout.close();
						}
						
				
#ifdef RUN_RAZERS_ONEGAPPED
						// count verifications
						String<ReadMatch<unsigned> > matches;
						RazerSOptions<RazerSSpec<false, true> > razersOptions;
						razersOptions.errorRate = (double)e / (double)totalN;
						razersOptions.errorRate += 0.0000001;
						razersOptions.threshold = t;
						razersOptions._debugLevel = 2;
		
						assign(razersOptions.shape, shapeString);
						mapReads(matches, testGenome, testReads, razersOptions);
#endif
						// write shape with its properties into file
						ofstream fout(datName.str().c_str(), ios::app | ios::out);
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
						fout << endl;
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
getParamsFilename(TSStr & paramsfile)
{
	paramsfile.str("");
	if(doSelectedGapped || doAllOneGapped)
	{
		paramsfile << fgparams<< fprefix[0]<<"_N" << totalN << "_E" << totalK;
		//if(prefixCount) paramsfile << fgparams<< fprefix[0]<<"_N" << totalN << "_E" << totalK;
		//else paramsfile << fgparams<<"userdef_N" << totalN << "_E" << totalK;
		if(optionHammingOnly) paramsfile << "_H";
		else paramsfile << "_L";
		if(doAllOneGapped) paramsfile << "_onegapped.dat";
		else paramsfile << ".dat";
	}
	else
	{
		paramsfile << fparams<< fprefix[0]<<"_QE0_N" << totalN << "_E" << totalK << ".dat";
		//if(prefixCount) paramsfile << fparams<< fprefix[0]<<"_QE0_N" << totalN << "_E" << totalK << ".dat";
		//else paramsfile << fparams<<"userdef_QE0_N" << totalN << "_E" << totalK << ".dat";
	}
}


//////////////////////////////////////////////////////////////////////////////
// Print usage
void printHelp(int, const char *[], bool longHelp = false) 
{
	cerr << "*******************************************************************" << endl;
	cerr << "*** Calculate efficient filter parameters (shape and threshold) ***" << endl;
	cerr << "*******************************************************************" << endl << endl;
	cerr << "Usage: paramChooser [OPTIONS]... " << endl;
	if (longHelp) {
		cerr << endl << "Options:" << endl;
		cerr << "  -n,  --length NUM            \t" << "sequence length (32)" << endl;
		cerr << "  -i,  --percent-identity NUM  \t" << "set the percent identity threshold (95)" << endl;
		cerr << "  -r,  --recognition-rate NUM  \t" << "set the percent recognition rate (99.0)" << endl;
		cerr << "  -pf, --prb-folder STR        \t" << "directory of [_prb.txt|.fastq|fastqint] files containing qualitiy values (optional)" << endl;
		cerr << "  -pq, --phred-qualities       \t" << "fastq files contain Phred qualities (default: Solexa qualities)" << endl;
		cerr << "  -d,  --error-distribution    \t" << "file containing mismatch probabilities (must contain at least n values, one value per line)" << endl;
		cerr << "  -pi, --prob-insert           \t" << "probability of an insertion (" << optionProbINSERT << ")" << endl;
		cerr << "  -pd, --prob-delete           \t" << "probability of a deletion (" << optionProbDELETE << ")" << endl;
		cerr << "                               \t" << "(for hamming-only filters use -pi 0 -pd 0)" << endl;
		cerr << "  -p,  --prefix STR            \t" << "session identifier (prefix of computed files);\n\t\t\t\tif also option d or pd is specified the prefix will be used for file naming\n\t\t\t\tuserspecific settings can be accessed in later session without recomputing loss rates\n\t\t\t\tby specifying the session id, i.e. prefix" << endl;
		cerr << "  -ug, --ungapped              \t" << "only consider ungapped shapes (off)" << endl;
		cerr << "  -og, --one-gapped            \t" << "only consider shapes containing at most one gap (of arbitrary length)" << endl;
		cerr << "  -h,  --help                  \t" << "print this help" << endl;
	} else {
		cerr << "Try 'chooseFilteringParameters --help' for more information." << endl;
	}
}

int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line
#ifdef LOSSRATE_VALIDATION	
	strncpy ( fparams, argv[0], length(argv[0])-12 );
	strcat(fparams,"./");
#else
	strncpy ( fparams, argv[0], length(argv[0])-12 );
	strcat(fparams,"params/");
#endif
	strncpy ( fgparams, argv[0], length(argv[0])-12 );
	strcat(fgparams,"gapped_params/");
	static const TFloat epsilon = 0.0000001;	

	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option
			if (strcmp(argv[arg], "-n") == 0 || strcmp(argv[arg], "--length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> totalN;
					if (!istr.fail())
						if (totalN < 1 || totalN > 100)
							cerr << "sequence length must be a value between 1 and 100" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv);
				return 0;
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
			if (strcmp(argv[arg], "-r") == 0 || strcmp(argv[arg], "--recognition-rate") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionLossRate;
					optionLossRate = 100.0-optionLossRate;
					optionLossRate /= 100.0;
					if (!istr.fail())
						if (optionLossRate < 0.0 || optionLossRate > 1.0)
							cerr << "Loss rate must be a value between 0 and 100" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-pi") == 0 || strcmp(argv[arg], "--prob-insert") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionProbINSERT;
					if (!istr.fail())
						if (optionProbINSERT < 0 || optionProbINSERT > 1)
							cerr << "Insert probability must be a value between 0 and 1" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv);
				return 0;
			}

			if (strcmp(argv[arg], "-pd") == 0 || strcmp(argv[arg], "--prob-delete") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> optionProbDELETE;
					if (!istr.fail())
						if (optionProbDELETE < 0 || optionProbDELETE > 1)
							cerr << "Delete probability must be a value between 0 and 1" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv);
				return 0;
			}


			if (strcmp(argv[arg], "-pf") == 0 || strcmp(argv[arg], "--prb-folder") == 0) {
				if (arg + 1 < argc) {
					++arg;
					fnameCount0 = true;
					fstream file;
					fname[0] = argv[arg];
				}
				else 
				{
					printHelp(argc, argv);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-p") == 0 || strcmp(argv[arg], "--prefix") == 0) { //prefix for previously computed param files
				if (arg + 1 < argc) {
					++arg;
					prefixCount = true;
					fstream file;
					fprefix[0] = argv[arg];
//					cout << "Session id prefix specified\n";
				}
				else 
				{
					printHelp(argc, argv);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-d") == 0 || strcmp(argv[arg], "--distribution-file") == 0) { //should also support fastq files
				if (arg + 1 < argc) {
					++arg;
					fnameCount1 = true;
					fstream file;
					fname[1] = argv[arg];
				}
				else 
				{
					printHelp(argc, argv);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-ha") == 0 || strcmp(argv[arg], "--hamming") == 0) {
				optionHammingOnly = true;
				continue;
			}

			if (strcmp(argv[arg], "-pq") == 0 || strcmp(argv[arg], "--phred-qualities") == 0) {
				solexaQual = false;
				continue;
			}

			if (strcmp(argv[arg], "-og") == 0 || strcmp(argv[arg], "--one-gapped") == 0) {
				optionChooseOneGappedOnly = true;       //optionChooseOneGappedOnly chooses shape with at most one gap
                                doUngapped = true;
                                doAllOneGapped = true;
                                doSelectedGapped = false;
				continue;
			}

			if (strcmp(argv[arg], "-ug") == 0 || strcmp(argv[arg], "--ungapped") == 0) {
                                if(optionChooseOneGappedOnly) continue;     //if both ungapped and onegapped specified --> optionChooseOneGappedOnly 
                                doUngapped = true;
                                doAllOneGapped = false;
                                doSelectedGapped = false;
				continue;
			}
			if (strcmp(argv[arg], "-mt") == 0 || strcmp(argv[arg], "--min-threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> minThreshold;
					if (!istr.fail())
						if (minThreshold < 1 || minThreshold > 3)
							cerr << "minimum threshold should be a value between 1 and 3" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv);
				return 0;
			}
			if (strcmp(argv[arg], "-mq") == 0 || strcmp(argv[arg], "--max-weight") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> maxWeight;
					if (!istr.fail())
						if (maxWeight < 6 || maxWeight > 14)
							cerr << "maximum weight should be a value between 6 and 14" << endl << endl;
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
		}
	}

	optionErrorRate += epsilon;
	optionLossRate += epsilon;
	
	if(optionProbINSERT <= epsilon && optionProbDELETE <= epsilon)
		optionHammingOnly=true;

	fill(firstTimeK,10,true);//set maximal number of errors considered in parameter computation to <10

// compute data specific loss rates
	if (fnameCount0 || fnameCount1) 
	{
		if(!prefixCount)
		{
			fprefix[0] = "userdef";
			cout << "\nNo session id given, using prefix 'userdef'\n";
		}
		String<TFloat> errorDistribution;
		resize(errorDistribution,totalN);
		//error distribution given --> read file containing error distr and compute loss rates
		if(fnameCount1)
		{
			fstream file;
			file.open(fname[1],ios_base::in | ios_base::binary);
			if(!file.is_open())
			{
				cout << "Couldn't open file "<<fname[1]<<"\n";
				return 0;
			}
			unsigned count = 0;
			char c = _streamGet(file);
			while(!_streamEOF(file) && count < totalN)
			{
				_parse_skipWhitespace(file,c);
				errorDistribution[count] = _parse_readEValue(file,c);// + (TFloat) 1.0/maxN;
				++count;
			}
			file.close();
			if(count != totalN)
			{
				cerr << "Error distribution file must contain at least " << totalN << " probability values (one value per line).\n";
				return 0;
			}
		}
		else // read qualtiy files and compute position dependent avg error probabilites
		{
			getAvgFromPrbDirectory(fname[0],errorDistribution);
		}

		fstream file;
		//if(prefixCount)
		if(doAllOneGapped) makeOneGappedStatsFile(errorDistribution);
		if(doSelectedGapped) makeSelectedStatsFile(errorDistribution);
		if(doUngapped) makeUngappedStatsFile(errorDistribution);
	}
	else if(!prefixCount) fprefix[0] = "results";

	totalK = (int)(optionErrorRate * totalN);
	
	//prioritize
	if(doSelectedGapped)
	{
		doAllOneGapped = false;
		doUngapped = false;
	}
	if(doAllOneGapped) doUngapped = false;
	
	// decide on which loss rate file to parse
	stringstream paramsfile;
	getParamsFilename(paramsfile);

	cout << "\nRead length      = " << totalN << "bp\n";
	cout << "Max num errors   = " << totalK << "\n";
	cout << "Recognition rate = " <<  100.0*(1.0-optionLossRate) << "%\n";
			
	// parse loss rate file and find appropriate filter criterium
	if(verbose)cout << "\n--> Reading " <<  paramsfile.str()<<"\n";
	fstream file;
	file.open(paramsfile.str().c_str(),ios_base::in | ios_base::binary);
	if(!file.is_open())
	{
		cout << "Couldn't open file "<<paramsfile.str()<<"\n";
		return 0;
	}
	else
	{
		if(doSelectedGapped || doAllOneGapped) parseGappedParams(file);
		else parseParams(file);
		cout << "\n Choose \nshape: " << chosenShape << "\n and \nthreshold: " << chosenThreshold<< "\n to achieve optimal performance for expected recognition rate >= " << (100.0-100.0*optionLossRate) << "% (expected recognition = " << (100.0-chosenLossRate*100.0) <<"%)\n\n";
		file.close();
	}

	return 0;
}
