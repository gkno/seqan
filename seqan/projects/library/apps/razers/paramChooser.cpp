#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;


// global input parameters
static unsigned totalN = 32;				// sequence length
static unsigned totalK = 2;				// errors
static double optionLossRate = 0.01;		// in
static double chosenLossRate = 0.0;		// out
static double optionErrorRate = 0.05;		// 

// output parameters
unsigned qgramLen = 4;			//
unsigned threshold = 2;			//

bool		fnameCount0 = 0;
bool		fnameCount1 = 0;
bool		prefixCount = 0;
const char	*fname[2] = { "","" };
const char	*fprefix[1] = { "" };
char	fparams[100];

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

//////////////////////////////////////////////////////////////////////////////
template<typename TFile, typename TChar>
inline double
_parse_readEValue(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
	// Read number
	String<char> str(c);
	bool e = false;
	double val1;
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
		return val1 * pow((double)10.0,(double)atof(toCString(str)));
	}	
 	else 
		return (double)atof(toCString(str));
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
        double bestSoFar = 0.0;
        unsigned bestQ = 1;
        unsigned bestT = 1;
	unsigned secondBestT = 0;
	unsigned qCutoff = 14; //depends on RAM
        //String<double> loss;
        //reserve(loss, 20*readLen);

        char c = _streamGet(file);
	if(_streamEOF(file)) std::cout << "Loss rate file is empty!\n";

        while(!_streamEOF(file))
        {
                double val = _parse_readEValue(file,c);
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
		qual = (int) _parse_readDouble(file,c);// + (double) 1.0/maxN;
		_parse_skipWhitespace(file,c);
		qual1 = (int) _parse_readDouble(file,c);// + (double) 1.0/maxN;
		int ind = (qual > qual1) ? 0 : 1;
		qual = (qual > qual1) ? qual : qual1;
		_parse_skipWhitespace(file,c);
		qual1 = (int) _parse_readDouble(file,c);// + (double) 1.0/maxN;
		ind = (qual > qual1) ? ind : 2;
		qual = (qual > qual1) ? qual : qual1;
		_parse_skipWhitespace(file,c);
		qual1 = (int) _parse_readDouble(file,c);// + (double) 1.0/maxN;
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
//	double sumavg = 0.0;
	for(unsigned t = 0; t < totalN; ++t)
	{
		long sum = 0;
		for(int q = 0; q<46;++q) 
			sum += countMatrix[q*totalN + t];
		for(int q = 0; q<46;++q) 
			avg[t] += (double) (countMatrix[q*totalN + t] * (q-5)) / sum;
		avg[t] =  pow((double)10.0,   (double)(avg[t]/-10)) / (1+pow((double)10.0,(double)(avg[t]/-10) ));
	}
	double avgsum = 0.0;

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
		errorDistribution[j] /= (double)countPrbs;
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

//////////////////////////////////////////////////////////////////////////////
// Returns the sum of two probability values in log space
template<typename TLog>
TLog
logAdd (TLog a, TLog b)
{
	if(isinf(a)) return b;
	if(isinf(b)) return a;
	if(isnan(a+log1p(exp(b-a)))) return a;
	return (a + log1p(exp(b-a)));
}



//////////////////////////////////////////////////////////////////////////////
// Returns log probability of q-gram-configuration q ending at position pos in sequence
template<typename TErrDistr>
typename Value<TErrDistr>::Type getProb(unsigned Q, unsigned q, unsigned pos, TErrDistr & logErrorDistr)
{
	typename Value<TErrDistr>::Type prob = 0.0;
	
	for(unsigned j = 0; j < Q; ++j)
	{
		if (q & 1) prob += (logErrorDistr[pos-Q+j+1]);
		else prob += (logErrorDistr[length(logErrorDistr)/2 + pos-Q+j+1]);
		q = q >> 1;
	}
	return prob;
	
}


//////////////////////////////////////////////////////////////////////////////
// DP recursion
template <typename TErrorDistr, typename TString>
void computeFilteringLoss(unsigned maxN, unsigned maxE, unsigned maxT, unsigned Q, TErrorDistr & errorDistr, TString & found)
{
////////////////////////////////////////////////
	//"global parameters"
	unsigned errorsPerQGram = 0;			// how many errors are allowed per qgram
	unsigned QPot = 1 << Q;
	unsigned char countBits[256];
	//taken from original main function

	// Count 1-bits in a byte
	for(unsigned i = 0; i < 256; ++i) {
		unsigned char bits = 0;
		for(unsigned char b = 0; b < 8; ++b)
			if ((i >> b) & 1) ++bits;
		countBits[i] = bits;
	}

	// Error probabilities
	// if we didn't get them from a file --> uniform distribution
	if(length(errorDistr) != maxN)
		fill(errorDistr,maxN,1.0/(long double) maxN);

	// prepare log error distribution 
	String<double> logErrorDistr;
	resize(logErrorDistr,2*maxN);

// 	for(unsigned j = 0; j < maxN; ++j)
// 		cout << errorDistr[j]<<"\t";
// 	cout << "\n";

	// log probs for seeing 1s at positions 0...maxN-1
	for(unsigned j = 0; j < maxN; ++j)
		logErrorDistr[j] = log(errorDistr[j]);

	// log probs for seeing 0s at positions 0...maxN-1
	for(unsigned j = 0; j < maxN; ++j)
		logErrorDistr[maxN+j] = log(1.0 - errorDistr[j]);
////////////////////////////////////////////////

	typedef String<long double> TMatrixCol;
	typedef typename Value<TErrorDistr>::Type TProbValue;

	// columns n-1 and n for recursion 
	TMatrixCol col0;
	TMatrixCol col1;
	fill(col0, maxE * QPot * maxT, log(0.0));
	resize(col1, maxE * QPot * maxT);

	// remembers sum of probabilities of all possible sequences having
	// exactly e errors (recursively updated during DP)
	TMatrixCol count_col;
	fill(count_col, maxE, log(0.0));

	//initialized with first position (can be either 0 or 1)
	count_col[0] = logErrorDistr[maxN];
	if(maxE > 1) count_col[1] = logErrorDistr[0];
	
	// same thing for last column i.e. sequence length maxN
	// (last column gets special treatment)
	TMatrixCol count_colFinal;
	fill(count_colFinal, maxE, log(0.0));

	// RECURSION BEGIN
	for(unsigned q = 0; q < QPot; ++q) 
	{
		col0[q*maxT] = log(1.0);
		// every bit in q set to 1 marks an error
		unsigned errors = 
			countBits[q & 255] + countBits[(q >> 8) & 255] + 
			countBits[(q >> 16) & 255] + countBits[(q >> 24) & 255];
		
		
		//for gapped q-grams only count those bits that are relevant positions
		
		// for n==0
		if (errors > errorsPerQGram) {
			// we miss 1 match for t>0 and q has more than <errorsPerQGram> errors
			// --> probability for finding the q-gram is 0
			for(unsigned t = 1; t < maxT; ++t) 
					col0[q*maxT+t] = log(0.0);
		} else {
			// we miss no match if q has no or one error
			// --> probability is 1.0
			col0[q*maxT+1] = log(1.0);
			for(unsigned t = 2; t < maxT; ++t)
				col0[q*maxT+t] = log(0.0);
			
		}
	}
//	dump(col0,0);

	// iterate over sequence length n
	TMatrixCol *col = &col1;
	TMatrixCol *colPrev = &col0;

	
/*	cout << "2:0";
	dump(col0, 0);
	cout << " :1";
	dump(col0, 1);
*/
	
	// RECURSION
	//
	// found(n,q,t,e) = (1-errorProb[n-Q]) * found(n-1,0|(q>>1),t-delta,e) delta=1/0 <-> q hat 0/>0 fehler
	//               + errorProb[n-Q] * found(n-1,1|(q>>1),t-delta,e-1)
	
	// rekursion (fuer q-gram matches <=1 fehler)
	// found(n,q,t,e) = (1-errorProb[n-Q]) * found(n-1,0|(q>>1),t-delta,e) delta=1/0 <-> q hat <=1/>1 fehler
	//               + errorProb[n-Q] * found(n-1,1|(q>>1),t-delta,e-1)
	
	for(unsigned n = Q; n < maxN; ++n)
	{
		if(n>Q) // if n-Q==0 count_col is already up to date (as this is how it was initialized)
		{	
			for(int e = maxE - 1; e > 0; --e)
			{
				// update with position n-Q: can be either 0 
				count_col[e] = ((logErrorDistr[maxN+n-Q]) + count_col[e]);
				// or 1
				count_col[e] = logAdd(count_col[e],logErrorDistr[n-Q] + count_col[e-1]);
				
			}
			// for 0 errors, updating with 0 is the only possibility
			count_col[0] += (logErrorDistr[maxN+n-Q]);
		}
		for(unsigned e = 0; e < maxE * QPot; e += QPot)
		{
			for(unsigned q = 0; q < QPot; ++q)
			{
				unsigned _q = (q << 1) & (QPot - 1);
				unsigned errors = 
					countBits[q & 255] + countBits[(q >> 8) & 255] + 
					countBits[(q >> 16) & 255] + countBits[(q >> 24) & 255];
				// again only count those bits that are relevant

				for(unsigned t = 0; t < maxT; ++t)
				{
					unsigned _t = t;
					if (_t > 0 && errors <= errorsPerQGram)
						--_t;
					
					long double recovered = ((logErrorDistr[maxN+n-Q]) + (*colPrev)[(e+(_q|0))*maxT+_t]);
					if (e > 0) recovered = logAdd(recovered,(logErrorDistr[n-Q] + (*colPrev)[((e-QPot)+(_q|1))*maxT+_t]));
					
					//if this is the last iteration (i.e. last columnof matrix), add log probability for the whole q gram
					if(n==maxN-1) recovered += getProb(Q,q,maxN-1,logErrorDistr);
					
					(*col)[(e+q)*maxT+t] = recovered;
				}
				
				// if this is the last iteration, add log probability of whole q-gram to the appropriate entry in count_col
				// könnte man auch erst später machen...
				if(e ==  ((maxE-1) * QPot) && n == maxN - 1)
				{
					long double prob = getProb(Q,q,maxN-1,logErrorDistr);
					
					//for gapped q-grams this needs to be changed! all positions are relevant, i.e. all errors must be counted!!!
					if(errors < maxE) 
						for(unsigned eSum = errors; eSum < maxE; ++eSum)
							count_colFinal[eSum] = logAdd(count_colFinal[eSum], (prob + count_col[eSum-errors]));
				}
			}

		}

		TMatrixCol *tmp = col;
		col = colPrev;
		colPrev = tmp;

/*		cout << n+1<<":0";
		dump(*colPrev, 0);
		cout << " :1";
		dump(*colPrev, 1);
		cout << " :2";http://news.google.de/
		dump(*colPrev, 2);
*/	}



	//cumulative sum for prob(found | num errors <=e)
//	for(unsigned eSum = 1; eSum < maxE; ++eSum)
//		count_colFinal[eSum] = logAdd(count_colFinal[eSum-1],count_colFinal[eSum]);

/*	for(unsigned i = 0; i < maxE; ++i)
		cout << exp(count_colFinal[i]) << "\t";
	cout << "\n\n";
*/
//0.0278743       0.138543

	// RECURSION END
	for(unsigned eSum = 0; eSum < maxE; ++eSum)
		for(unsigned t = 0; t < maxT; ++t) 
		{
			long double recovered = log(0.0);
			for(unsigned q = 0; q < QPot; ++q) 
			{
				//again: for gapped q-grams this needs to be changed! all positions are relevant, i.e. all errors must be counted!!!
				unsigned errors = 
					countBits[q & 255] + countBits[(q >> 8) & 255] + 
					countBits[(q >> 16) & 255] + countBits[(q >> 24) & 255];
				if (errors <= eSum) {
					unsigned e = eSum - errors;
					recovered = logAdd(recovered,(*colPrev)[(e*QPot+q)*maxT+t]);
				}
			}
			//divide by probabilitiy of seeing eSum errors
			found[eSum*maxT+t] = recovered - count_colFinal[eSum];
		}


}

/*
template<typename TError>
void
makeStatsFile(TError & errorDistr)
{
	int maxErrors = (int) totalN / 5;
	for(int e = 0; e < maxErrors; ++e)
	{
		std::fstream outfile;
		std::stringstream filename;
		filename << "./params/"<<fprefix[0]<<"_QE0_N" << totalN<< "_E" << e<<".dat";
		std::cout << "Creating file "<<filename.str() << "\n";
		outfile.open(filename.str().c_str(),ios_base::out);
		for(int qLen = 1; qLen < 16; ++qLen)
		{
			String<double> found;
			resize(found,totalN*(e+1));
			computeFilteringLoss(totalN,e+1,totalN,qLen,errorDistr,found);
			for(unsigned t = 0; t < totalN; ++t) 
			{
				if (t > 0) outfile << "\t";
				outfile.precision(8);
				outfile << (1.0-exp(found[e*totalN+t]));
			//	cout.precision(3);
			//	if (t > 0) cout << "\t";
			//	cout << (1.0-exp(found[e*totalN+t]));
			}
			outfile << endl;	
			//cout << endl;	
		}
	}

}*/


template<typename TError>
void
makeStatsFile(TError & errorDistr)
{
	unsigned maxErrors = (unsigned) totalN / 5;
//maxErrors = 1;
/*	vector<FILE*> outfiles;
	outfiles.resize(maxErrors);
	for(unsigned e=0; e < maxErrors; ++e)
	{
		std::stringstream filename;
		filename << "./params/"<<fprefix[0]<<"_QE0_N" << totalN<< "_E" << e<<".dat";
		std::cout << "Creating file "<<filename.str() << "\n";
		outfiles[e] = fopen (filename.str().c_str(),"r");
 		fprintf(outfiles[e],"%c",'t');
//		outfiles[e].open(filename.str().c_str(),ios_base::out);
	}*/
	for(int qLen = 1; qLen < 16; ++qLen)
	{
		String<double> found;
		resize(found,totalN*maxErrors);
		computeFilteringLoss(totalN,maxErrors,totalN,qLen,errorDistr,found);
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
			for(unsigned t = 0; t < totalN; ++t) 
			{
				if (t > 0) outfile << "\t";
				outfile.precision(8);
				//write(outfiles[e],1.0-exp(found[e*totalN+t]));
				//outfiles[e].write(1.0-exp(found[e*totalN+t]));
				outfile << (1.0-exp(found[e*totalN+t]));
			//	cout.precision(3);
			//	if (t > 0) cout << "\t";
			//	cout << (1.0-exp(found[e*totalN+t]));
			}
//			fprintf(outfiles[e],"%c",'\n');
//			outfiles[e].write("\n");
			outfile << endl;	
			outfile.close();
		}
//		cout << endl;	
	}
// 	for(unsigned e=0; e < maxErrors; ++e)
// 	{
// //		outfiles[e].close();
// 		fclose(outfiles[e]);
// 	}

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
		cerr << "  -pd, --prb-directory STR     \t" << "directory of _prb.txt files containing qualtiy values (optional)" << endl;
		cerr << "  -d,  --distribution-file STR \t" << "file containing user-defined (or precomputed) error probabilities" << endl;
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
	static const double epsilon = 0.0000001;	

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
			if (strcmp(argv[arg], "-pd") == 0 || strcmp(argv[arg], "--prb-directory") == 0) { //should also support fastq files
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
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, true);
				return 0;
			}
		}
	}

	optionErrorRate += epsilon;
	optionLossRate += epsilon;
	
	// compute data specific loss rates
	if (fnameCount0 || fnameCount1) 
	{
		if(!prefixCount)
		{
			fprefix[0] = "userdef";
			std::cout << "\nNo session id given, using prefix 'userdef'\n";
		}
		String<double> errorDistribution;
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
				errorDistribution[count] = _parse_readEValue(file,c);// + (double) 1.0/maxN;
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
		makeStatsFile(errorDistribution);
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
