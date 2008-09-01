#define USE_LOGVALUES		// this is recommended when using probability values

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
#include "bestOneGapped.h"


using namespace seqan;
using namespace std;

typedef long double TFloat;

// global input parameters
static unsigned totalN = 32;				// sequence length
static unsigned totalK = 2;					// errors
static TFloat optionLossRate = 0.01;		// in
static TFloat chosenLossRate = 0.0;			// out
static TFloat optionErrorRate = 0.05;		// 
static bool optionHammingOnly = false;
static bool doUngapped = true;
static bool doAllOneGapped = false;
static bool doSelectedGapped = false;
static TFloat optionProbINSERT = 0.02;
static TFloat optionProbDELETE = 0.02;

// output parameters
unsigned qgramLen = 4;			//
unsigned threshold = 2;			//

bool		fnameCount0 = 0;
bool		fnameCount1 = 0;
bool		prefixCount = 0;
const char	*fname[2] = { "","" };
const char	*fprefix[1] = { "" };
char		fparams[100];

//////////////////////////////////////////////////////////////////////////////

void seqan::initCountBits();

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
	//return (((unsigned) c >=  48) && ((unsigned) c <=  57));
	return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
		    (c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
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
//	std::cout <<  files[length(files)-1] << " ?\n";
    }
    closedir(dp);
    return 0;
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
        unsigned bestT = 1;
	unsigned secondBestT = 0;
	unsigned qCutoff = 14; //depends on RAM
        //String<TFloat> loss;
        //reserve(loss, 20*readLen);

        char c = _streamGet(file);
	if(_streamEOF(file)) std::cout << "Loss rate file is empty!\n";

        while(!_streamEOF(file))
        {
                TFloat val = _parse_readEValue(file,c);
 //            std::cout << "\n"<<val;
//              appendValue(loss,val);  //t=0
                countT = 1;
                while(!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')))
                {
                        _parse_skipWhitespace(file,c);
                        val = _parse_readEValue(file,c);
  //                  std::cout << " " << val;
//                      appendValue(loss,val);
                        if(countT > 1 && val <= optionLossRate /*&& val > bestSoFar*/)
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
		if(countQ>=qCutoff)
			break;
                _parse_skipWhitespace(file,c);
        }
	if(bestT<2) std::cout << "\n!!! Something wrong with file? !!!\n";
	chosenLossRate = bestSoFar;
//      countT = length(loss)/countQ;
//      std::cout << "insge: Q="<<countQ << " T="<<countT<<"\n";
        qgramLen = bestQ;
        threshold = bestT;
	if(bestQ==14)
		std::cout << "\nIf RAM <= 1GB : Choose q = "<<bestQ-1 << " and t = "<<secondBestT<<"\nelse          : ";
	else std::cout <<"\n";
        return true;
}

//compute average position dependent error distribution (assumes solexa qualtiy values in prb.txt format)
template<typename TFile, typename TDistribution>
void
qualityDistributionFromPrbFile(TFile & file, TDistribution & avg)
{
	String<int> countMatrix;
	fill(countMatrix,totalN*46,0);

	char c;
	if(!_streamEOF(file))
	{
		c = _streamGet(file);
		_parse_skipWhitespace(file,c);
	}
	int read_count = 0;
	int qual,qual1;
	unsigned pos = 0;
	while(!_streamEOF(file))
	{
		qual = (int) _parse_readDouble(file,c);// + (TFloat) 1.0/maxN;
		_parse_skipWhitespace(file,c);
		qual1 = (int) _parse_readDouble(file,c);// + (TFloat) 1.0/maxN;
		int ind = (qual > qual1) ? 0 : 1;
		qual = (qual > qual1) ? qual : qual1;
		_parse_skipWhitespace(file,c);
		qual1 = (int) _parse_readDouble(file,c);// + (TFloat) 1.0/maxN;
		ind = (qual > qual1) ? ind : 2;
		qual = (qual > qual1) ? qual : qual1;
		_parse_skipWhitespace(file,c);
		qual1 = (int) _parse_readDouble(file,c);// + (TFloat) 1.0/maxN;
		ind = (qual > qual1) ? ind : 3;
		qual = (qual > qual1) ? qual : qual1;
		++countMatrix[totalN * (qual+5) + pos];
		++pos;
		if(pos==totalN) pos = 0;
		_parse_skipWhitespace(file,c);
		//std::cout << ind;
		if(pos==0) {
		//	std::cout << endl;
			++read_count;
		}
	}
//	std::cout <<"\nNumber of reads: "<< read_count << "\n";
	std::cout << " Readcount = "<< read_count << "\n";



	fill(avg,totalN,0.0);
//	TFloat sumavg = 0.0;
	for(unsigned t = 0; t < totalN; ++t)
	{
		long sum = 0;
		for(int q = 0; q<46;++q) 
			sum += countMatrix[q*totalN + t];
		for(int q = 0; q<46;++q) 
			avg[t] += (TFloat) (countMatrix[q*totalN + t] * (q-5)) / sum;
		avg[t] =  pow((TFloat)10.0,   (TFloat)(avg[t]/-10)) / (1+pow((TFloat)10.0,(TFloat)(avg[t]/-10) ));
	}
	TFloat avgsum = 0.0;

	for(unsigned t = 0; t < totalN; ++t)
		avgsum += avg[t];

/*	for(unsigned t = 0; t < totalN; ++t)
	{
		cout.precision(10);
		//std::cout << avg[t]/avgsum << "\t";
		std::cout << avg[t] << "\t";
	}
	std::cout << "\n";
*/
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
		if(suffix(files[i],length(files[i])-8) == "_prb.txt")
		{
			std::cout << "Processing "<< files[i] << "...\n";
			TError avg_act;
			resize(avg_act,totalN);
			std::fstream filestrm;
			stringstream sstrm;
			sstrm << prbPath << files[i];
			filestrm.open(sstrm.str().c_str(),ios_base::in);
			qualityDistributionFromPrbFile(filestrm,avg_act);
			filestrm.close();
			for(unsigned j=0; j < totalN; ++j)
			{
//				std::cout << " " << avg_act[j];
				errorDistribution[j] += avg_act[j];
			}
			++countPrbs;
		}
	}
	for(unsigned j=0; j < totalN; ++j)
		errorDistribution[j] /= (TFloat)countPrbs;
	std::cout << "Writing average error probabilities to " << fprefix[0] << "_errorProb.dat\n";
	fstream out;
	stringstream avgOut;
	avgOut << fprefix[0] << "_errorProb.dat";
	out.open(avgOut.str().c_str(),ios_base::out);
	if(!out.is_open()) std::cout << "Couldn't write to file "<<avgOut.str()<<"\n";
	else
		for(unsigned j=0; j < totalN; ++j)
			out << errorDistribution[j] << "\n";
	out.close();

	
}


template<typename TError>
void
makeUngappedStatsFile(TError & errorDistr)
{
	unsigned maxErrors = 5;//(unsigned) totalN / 5;
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
		
		//for loss rate sanity check
//                for(unsigned e=0; e < maxErrors; ++e)
//                {
//                        std::stringstream filename;
//                        filename << fparams <<fprefix[0]<<"_QE0_N" << totalN<< "_E" << e<<".dat";
//                        ofstream outfile;
//                        if(qLen==8){
//                                outfile.open(filename.str().c_str(),ios_base::out);
//                                std::cout << "Creating file "<<filename.str() << "\n";
//                        }
//                        else outfile.open(filename.str().c_str(),ios_base::app);
//                        for(unsigned t = 20; t > 0; --t)
//                        {
//                                outfile.precision(10);
//                                //write(outfiles[e],1.0-exp(found[e*totalN+t]));
//                                //outfiles[e].write(1.0-exp(found[e*totalN+t]));
//                                outfile << (1.0 - _transformBack(found[e*maxT+t]));
//                                if(t>1)outfile << "\n";
//                        //      cout.precision(3);
//                        //      if (t > 0) cout << "\t";
//                        //      cout << (1.0-exp(found[e*totalN+t]));
//                        }
// //                     fprintf(outfiles[e],"%c",'\n');
// //                     outfiles[e].write("\n");
//                        outfile << endl;
//                        outfile.close();
//                }

		//regular output
		for(unsigned e=0; e < maxErrors; ++e)
		{
			std::stringstream filename;
			filename << fparams <<fprefix[0]<<"_QE0_N" << totalN<< "_E" << e<<".dat";
			ofstream outfile;
			if(qLen==1){
				outfile.open(filename.str().c_str(),ios_base::out);
				std::cout << "Creating file "<<filename.str() << "\n";
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
	}
// 	for(unsigned e=0; e < maxErrors; ++e)optionHammingOnly
// 	{
// //		outfiles[e].close();
// 		fclose(outfiles[e]);
// 	}

}

template<typename TError>
void
makeSelectedStatsFile(TError & errorDistr)
{

//	unsigned maxE = (unsigned) totalN / 5;
	unsigned maxErrors = 5;
	unsigned minQ = 7;
	unsigned maxT = totalN-minQ+1;

	typedef typename Value<TError>::Type TErrorValue;
	String<TErrorValue> logErrorDistribution;
	
	String<CharString> shapeStrings;
	//q=7
	appendValue(shapeStrings,"1111111");
	appendValue(shapeStrings,"1111100011");
	appendValue(shapeStrings,"111110000011");
	//q=8
	appendValue(shapeStrings,"11111111");
	appendValue(shapeStrings,"11111000111");
	appendValue(shapeStrings,"1111110000011");
	appendValue(shapeStrings,"11111100000011");
	//q=9
	appendValue(shapeStrings,"111111111");
	appendValue(shapeStrings,"111111000111");
	appendValue(shapeStrings,"1111111000011");
	appendValue(shapeStrings,"111111000000111");
	//q=10
	appendValue(shapeStrings,"1111111111");
	appendValue(shapeStrings,"111111111001");
	appendValue(shapeStrings,"1111111100011");
	appendValue(shapeStrings,"11111110000111");
	appendValue(shapeStrings,"111111110000011");

	//q=11
	appendValue(shapeStrings,"11111111111");
	appendValue(shapeStrings,"11111111000111");
	appendValue(shapeStrings,"111111111000011");
	appendValue(shapeStrings,"1111111110000011");
	appendValue(shapeStrings,"11111111000000111");
	//q=12
	appendValue(shapeStrings,"111111111111");
	appendValue(shapeStrings,"1111111111000011");
	appendValue(shapeStrings,"11111111110000011");
	appendValue(shapeStrings,"11111111000001111");
	//q=13
	appendValue(shapeStrings,"1111111111111");
	appendValue(shapeStrings,"111111111100000111");
	appendValue(shapeStrings,"1111111111000000111");
	appendValue(shapeStrings,"111111111000001111");
	appendValue(shapeStrings,"1111111111100000011");
	//q=14
	appendValue(shapeStrings,"11111111111111");
	appendValue(shapeStrings,"1111111111100000111");
	appendValue(shapeStrings,"1111111111000001111");
	appendValue(shapeStrings,"11111111110000001111");
	appendValue(shapeStrings,"11111111111000000111");

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

	//loss rate buckets
	map<TErrorValue, unsigned> lossRateBuckets;
	String<TErrorValue> lossRatesProbe;
	resize(lossRatesProbe, 9);
	lossRatesProbe[0] = 0;
	lossRatesProbe[1] = 0.00001;
	lossRatesProbe[2] = 0.0001;
	lossRatesProbe[3] = 0.001;
	lossRatesProbe[4] = 0.01;
	lossRatesProbe[5] = 0.02;
	lossRatesProbe[6] = 0.05;
	lossRatesProbe[7] = 0.1;
	lossRatesProbe[8] = 0.2;

	
	for(unsigned j = 0; j <= length(lossRatesProbe); ++j){
		lossRateBuckets.insert ( pair<TErrorValue,unsigned>(lossRatesProbe[j],j) );
	}
	
	String<bool> firstTimeK;
	fill(firstTimeK,maxErrors*length(lossRatesProbe)+1,true);
	
	for(unsigned i = 0; i < length(shapeStrings); ++i)
	{
		
		String<TFloat> found;
		resize(found,maxT*maxErrors);
		
		String< State<TFloat> > states;
		std::cout << "do DP\n";
		initPatterns(states, shapeStrings[i], maxErrors-1, logErrorDistribution, optionHammingOnly);
		computeFilteringLoss(found, states, length(shapeStrings[i]), maxT, maxErrors,  logErrorDistribution);
		std::cout << "Printing\n";
		for(unsigned e = 0; e < maxErrors; ++e) {
			bool highestOptimalFound = false;
			for(unsigned t = maxT-1; t > 0; --t) {
				TFloat lossrate = 1.0 - (TFloat) _transformBack(found[e*maxT+t]);
				typename std::map<TFloat, unsigned>::iterator it, itlow, itup, itend, itbegin;
				if(lossrate <= 0.0){
					if(highestOptimalFound) continue;
					else highestOptimalFound = true;
				}
				// points to first item in loss rate file
				itbegin = lossRateBuckets.begin();
				if(lossrate < ((*itbegin).first)){ 
					continue;
				}
				itend = lossRateBuckets.end();
				// points to last item in loss rate file
				--itend;
				
				if(lossrate > ((*itend).first)){ 
					continue;
				}
				
				TFloat helpLoss = (*itend).first; 
				if(lossrate > helpLoss){
					continue;
				}
				
				
				itup = lossRateBuckets.upper_bound(lossrate);
				itlow = itup;
				--itlow;

				unsigned gminCov = getMinCov(weights[i], length(shapeStrings[i]), t);

				unsigned index = unsigned((*itlow).second);	
			
				// create the whole file name
				stringstream datName;
				if(best_shape_helpFolder) datName << best_shape_folder;
				else datName << "gapped_params";
				datName << "/"<<fprefix[0]<<"_" << totalN << "_" << e << "_";
				if(!optionHammingOnly) datName << "L_";
				else datName <<"H_";
				datName << (*itlow).first << "_" << (*itup).first << ".dat";
				
			
				// if datName-file doesnt exist, write the title on it
				if(firstTimeK[maxErrors*index + e]==true){
					firstTimeK[maxErrors*index + e] = false;
					ofstream fout(datName.str().c_str(), ios::out);
					fout << "shape\t\tt\t\tloss rate\t\tminCoverage\n\n";
					fout.close();
				}
				// write best shape with its properties on the file
				ofstream fout(datName.str().c_str(), ios::app | ios::out);
				fout << shapeStrings[i] << "\t\t";
				fout << t << "\t\t";
				fout << lossrate << "\t\t";
				fout << gminCov << endl; 
				fout.close();
				
			} // t-loop
		}

	
	}


}


template<typename TError>
void
makeOneGappedStatsFile(TError & errorDistr)
{

//	unsigned maxE = (unsigned) totalN / 5;
	unsigned maxE = 5;

	unsigned minE = 0;				
	unsigned maxQ = 14;				// weights are considered from minQ..maxQ
	unsigned minQ = 12;
	unsigned minGap = 4; 
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

	
		
	String<TErrorValue> found;
	resize(found,maxT*maxE);

	initCountBits();

	map<TErrorValue, unsigned> lossRateBuckets;
	String<TErrorValue> lossRatesProbe;
	resize(lossRatesProbe, 9);
	lossRatesProbe[0] = 0;
	lossRatesProbe[1] = 0.00001;
	lossRatesProbe[2] = 0.0001;
	lossRatesProbe[3] = 0.001;
	lossRatesProbe[4] = 0.01;
	lossRatesProbe[5] = 0.02;
	lossRatesProbe[6] = 0.05;
	lossRatesProbe[7] = 0.1;
	lossRatesProbe[8] = 0.2;

	
	for(unsigned j = 0; j <= length(lossRatesProbe); ++j){
		lossRateBuckets.insert ( pair<TErrorValue,unsigned>(lossRatesProbe[j],j) );
	}


//	best_shape(totalN, found, logErrorDistribution, lossRates, 12, 10, 3, 6, 4, 0, maxT);
	best_shape(totalN, found, logErrorDistribution, lossRateBuckets, maxQ, minQ, minGap, maxGap, maxE, minE, maxT,optionHammingOnly);

}


//////////////////////////////////////////////////////////////////////////////
// Print usage
void printHelp(int, const char *[], bool longHelp = false) 
{
	cerr << "*****************************************************" << endl;
	cerr << "*** Calculate efficient filter parameters q and t ***" << endl;
	cerr << "*****************************************************" << endl << endl;
	cerr << "Usage: paramChooser [OPTIONS]... " << endl;
	if (longHelp) {
		cerr << endl << "Options:" << endl;
		cerr << "  -n,  --length NUM            \t" << "sequence length (32)" << endl;
		cerr << "  -i,  --percent-identity NUM  \t" << "set the percent identity threshold (95)" << endl;
		cerr << "  -r,  --recognition-rate NUM  \t" << "set the percent recognition rate (99.0)" << endl;
		cerr << "  -pf, --prb-folder STR        \t" << "directory of _prb.txt files containing qualtiy values (optional)" << endl;
		cerr << "  -d,  --error-distribution    \t" << "file containing mismatch probabilities (must contain at least n values, one value per line)" << endl;
		cerr << "  -pi, --prob-insert           \t" << "probability of an insertion (" << optionProbINSERT << ")" << endl;
		cerr << "  -pd, --prob-delete           \t" << "probability of a deletion (" << optionProbDELETE << ")" << endl;
		cerr << "                               \t" << "(for hamming-only filters use -pi 0 -pd 0)" << endl;
		cerr << "  -p,  --prefix STR            \t" << "session identifier (prefix of computed files);\n\t\t\t\tif also option d or pd is specified the prefix will be used for file naming\n\t\t\t\tuserspecific settings can be accessed in later session without recomputing loss rates\n\t\t\t\tby specifying the session id, i.e. prefix" << endl;
		cerr << "  -h,  --help                  \t" << "print this help" << endl;
	} else {
		cerr << "Try 'chooseFilteringParameters --help' for more information." << endl;
	}
}

int main(int argc, const char *argv[]) 
{
	//////////////////////////////////////////////////////////////////////////////
	// Parse command line
	strncpy ( fparams, argv[0], length(argv[0])-12 );
	strcat(fparams,"params/");
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


			if (strcmp(argv[arg], "-pf") == 0 || strcmp(argv[arg], "--prb-folder") == 0) { //should also support fastq files
				if (arg + 1 < argc) {
					++arg;
					fnameCount0 = true;
					std::fstream file;
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
					std::fstream file;
					fprefix[0] = argv[arg];
//					std::cout << "Session id prefix specified\n";
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
					std::fstream file;
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
	{
		optionHammingOnly=true;
		cout << "hier";	
	}
	// compute data specific loss rates
	if (fnameCount0 || fnameCount1) 
	{
		if(!prefixCount)
		{
			fprefix[0] = "userdef";
			std::cout << "\nNo session id given, using prefix 'userdef'\n";
		}
		String<TFloat> errorDistribution;
		resize(errorDistribution,totalN);
		//error distribution given --> read file containing error distr and compute loss rates
		if(fnameCount1)
		{
			std::fstream file;
			file.open(fname[1],std::ios_base::in | std::ios_base::binary);
			if(!file.is_open())
			{
				std::cout << "Couldn't open file "<<fname[1]<<"\n";
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

		std::fstream file;
		//if(prefixCount)
		if(doAllOneGapped) makeOneGappedStatsFile(errorDistribution);
		if(doSelectedGapped) makeSelectedStatsFile(errorDistribution);
		if(doUngapped) makeUngappedStatsFile(errorDistribution);
	}
	
	totalK = (int)(optionErrorRate * totalN);
	
	// decide on which loss rate file to parse
	std::stringstream paramsfile;
	if (fnameCount0 > 0 || fnameCount1 > 0)
	{
		if(prefixCount) paramsfile << fparams<< fprefix[0]<<"_QE0_N" << totalN << "_E" << totalK << ".dat";
		else paramsfile << fparams<<"userdef_QE0_N" << totalN << "_E" << totalK << ".dat";
	}
	else{
		if(prefixCount) paramsfile << fparams<< fprefix[0]<<"_QE0_N" << totalN << "_E" << totalK << ".dat";
		else paramsfile << fparams<<"results_QE0_N" << totalN << "_E" << totalK << ".dat";
	}

	std::cout << "\nRead length      = " << totalN << "bp\n";
	std::cout << "Max num errors   = " << totalK << "\n";
	std::cout << "Recognition rate = " <<  100.0*(1.0-optionLossRate) << "%\n";
			
	// parse loss rate file and find appropriate filter criterium
	std::cout << "\n--> Reading " <<  paramsfile.str()<<"\n";
	fstream file;
	file.open(paramsfile.str().c_str(),ios_base::in | ios_base::binary);
	if(!file.is_open())
	{
		std::cout << "Couldn't open file "<<paramsfile.str()<<"\n";
		return 0;
	}
	else parseParams(file);
		
	// suggest a suitable combination of q and t
	std::cout << "Choose q = " << qgramLen << " and t = " << threshold << " to achieve optimal performance for expected recognition rate >= " << (100.0-100.0*optionLossRate) << "% (expected recognition = " << (100.0-chosenLossRate*100.0) <<"%)\n\n";

	return 0;
}
