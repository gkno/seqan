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
#include "paramChooser.h"



using namespace seqan;
using namespace std;

typedef double TFloat;



//////////////////////////////////////////////////////////////////////////////
// 
// template<typename TSStr>
// void
// getParamsFilename(TSStr & paramsfile)
// {
// 	paramsfile.str("");
// 	if(doSelectedGapped || doAllOneGapped)
// 	{
// 		paramsfile << pm_options.fgparams<< pm_options.fprefix[0]<<"_N" << pm_options.totalN << "_E" << totalK;
// 		//if(prefixCount) paramsfile << fgparams<< fprefix[0]<<"_N" << totalN << "_E" << totalK;
// 		//else paramsfile << fgparams<<"userdef_N" << totalN << "_E" << totalK;
// 		if(pm_options.optionHammingOnly) paramsfile << "_H";
// 		else paramsfile << "_L";
// 		if(pm_options.doAllOneGapped) paramsfile << "_onegapped.dat";
// 		else paramsfile << ".dat";
// 	}
// 	else
// 	{
// 		paramsfile << pm_options.fparams<< pm_options.fprefix[0]<<"_QE0_N" << pm_options.totalN << "_E" << totalK << ".dat";
// 		//if(prefixCount) paramsfile << fparams<< fprefix[0]<<"_QE0_N" << totalN << "_E" << totalK << ".dat";
// 		//else paramsfile << fparams<<"userdef_QE0_N" << totalN << "_E" << totalK << ".dat";
// 	}
// }


//////////////////////////////////////////////////////////////////////////////
// Print usage
void printHelp(int, const char *[], ParamChooserOptions & pm_options, bool longHelp = false) 
{
	cerr << "*******************************************************************" << endl;
	cerr << "*** Calculate efficient filter parameters (shape and threshold) ***" << endl;
	cerr << "*******************************************************************" << endl << endl;
	cerr << "Usage: paramChooser [OPTIONS]... " << endl;
	if (longHelp) {
		cerr << endl << "Options:" << endl;
		cerr << "  -n,  --length NUM            \t" << "sequence length ("<<pm_options.totalN<<")" << endl;
		cerr << "  -i,  --percent-identity NUM  \t" << "set the percent identity threshold (95)" << endl;
		cerr << "  -r,  --recognition-rate NUM  \t" << "set the percent recognition rate (99.0)" << endl;
		cerr << "  -pf, --prb-folder STR        \t" << "directory of [_prb.txt|.fastq|fastqint] files containing qualitiy values (optional)" << endl;
		cerr << "  -pq, --phred-qualities       \t" << "fastq files contain Phred qualities (default: Solexa qualities)" << endl;
		cerr << "  -d,  --error-distribution    \t" << "file containing mismatch probabilities (must contain at least n values, one value per line)" << endl;
		cerr << "  -pi, --prob-insert           \t" << "probability of an insertion (" << pm_options.optionProbINSERT << ")" << endl;
		cerr << "  -pd, --prob-delete           \t" << "probability of a deletion (" << pm_options.optionProbDELETE << ")" << endl;
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
	static const TFloat epsilon = 0.0000001;	

        RazerSOptions<> r_options;
        ParamChooserOptions pm_options;

	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option
			if (strcmp(argv[arg], "-n") == 0 || strcmp(argv[arg], "--length") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.totalN;
					if (!istr.fail())
						if (pm_options.totalN < 1 || pm_options.totalN > 100)
							cerr << "sequence length must be a value between 1 and 100" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-i") == 0 || strcmp(argv[arg], "--percent-identity") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionErrorRate;
					pm_options.optionErrorRate = (100.0 - pm_options.optionErrorRate) / 100.0;
					if (!istr.fail())
						if (pm_options.optionErrorRate < 0 || pm_options.optionErrorRate > 0.1)
							cerr << "Percent identity threshold must be a value between 90 and 100" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-r") == 0 || strcmp(argv[arg], "--recognition-rate") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionLossRate;
					pm_options.optionLossRate = 100.0 - pm_options.optionLossRate;
					pm_options.optionLossRate /= 100.0;
					if (!istr.fail())
						if (pm_options.optionLossRate < 0.0 || pm_options.optionLossRate > 1.0)
							cerr << "Loss rate must be a value between 0 and 100" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-pi") == 0 || strcmp(argv[arg], "--prob-insert") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionProbINSERT;
					if (!istr.fail())
						if (pm_options.optionProbINSERT < 0 || pm_options.optionProbINSERT > 1)
							cerr << "Insert probability must be a value between 0 and 1" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}

			if (strcmp(argv[arg], "-pd") == 0 || strcmp(argv[arg], "--prob-delete") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.optionProbDELETE;
					if (!istr.fail())
						if (pm_options.optionProbDELETE < 0 || pm_options.optionProbDELETE > 1)
							cerr << "Delete probability must be a value between 0 and 1" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}


			if (strcmp(argv[arg], "-pf") == 0 || strcmp(argv[arg], "--prb-folder") == 0) {
				if (arg + 1 < argc) {
					++arg;
					pm_options.fnameCount0 = true;
					fstream file;
					pm_options.fname[0] = argv[arg];
				}
				else 
				{
					printHelp(argc, argv, pm_options);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-p") == 0 || strcmp(argv[arg], "--prefix") == 0) { //prefix for previously computed param files
				if (arg + 1 < argc) {
					++arg;
					pm_options.prefixCount = true;
					fstream file;
					pm_options.fprefix[0] = argv[arg];
//					cout << "Session id prefix specified\n";
				}
				else 
				{
					printHelp(argc, argv, pm_options);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-d") == 0 || strcmp(argv[arg], "--distribution-file") == 0) { //should also support fastq files
				if (arg + 1 < argc) {
					++arg;
					pm_options.fnameCount1 = true;
					fstream file;
					pm_options.fname[1] = argv[arg];
				}
				else 
				{
					printHelp(argc, argv, pm_options);
					return 0;
				}
			}
			if (strcmp(argv[arg], "-ha") == 0 || strcmp(argv[arg], "--hamming") == 0) {
				pm_options.optionHammingOnly = true;
				continue;
			}

			if (strcmp(argv[arg], "-pq") == 0 || strcmp(argv[arg], "--phred-qualities") == 0) {
				pm_options.solexaQual = false;
				continue;
			}

			if (strcmp(argv[arg], "-og") == 0 || strcmp(argv[arg], "--one-gapped") == 0) {
				pm_options.optionChooseOneGappedOnly = true;       //optionChooseOneGappedOnly chooses shape with at most one gap
                                pm_options.doUngapped = true;
                                pm_options.doAllOneGapped = true;
                                pm_options.doSelectedGapped = false;
				continue;
			}

			if (strcmp(argv[arg], "-ug") == 0 || strcmp(argv[arg], "--ungapped") == 0) {
                                if(pm_options.optionChooseOneGappedOnly) continue;     //if both ungapped and onegapped specified --> optionChooseOneGappedOnly 
                                pm_options.doUngapped = true;
                                pm_options.doAllOneGapped = false;
                                pm_options.doSelectedGapped = false;
				continue;
			}
			if (strcmp(argv[arg], "-mt") == 0 || strcmp(argv[arg], "--min-threshold") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.minThreshold;
					if (!istr.fail())
						if (pm_options.minThreshold < 1 || pm_options.minThreshold > 3)
							cerr << "minimum threshold should be a value between 1 and 3" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}
			if (strcmp(argv[arg], "-mq") == 0 || strcmp(argv[arg], "--max-weight") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> pm_options.maxWeight;
					if (!istr.fail())
						if (pm_options.maxWeight < 6 || pm_options.maxWeight > 14)
							cerr << "maximum weight should be a value between 6 and 14" << endl << endl;
						else
							continue;
				}
				printHelp(argc, argv, pm_options);
				return 0;
			}

			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, pm_options, true);
				return 0;
			}
		}
	}
        pm_options.optionErrorRate += epsilon;
	pm_options.optionLossRate += epsilon;
	
	pm_options.verbose = true;
	r_options._debugLevel = 1;
	r_options.errorRate = pm_options.optionErrorRate;

	
	if(pm_options.optionProbINSERT <= epsilon && pm_options.optionProbDELETE <= epsilon)
		pm_options.optionHammingOnly=true;

	fill(pm_options.firstTimeK,20,true);//set maximal number of errors considered in parameter computation to <10

        chooseParams(r_options, pm_options);
/*
// compute data specific loss rates
	if (fnameCount0 || pm_options.fnameCount1) 
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
	file
	// decide on which loss rate file to parse
	stringstream paramsfile;
	getParamsFilename(paramsfile,options);

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
		if(doSelectedGapped || doAllOneGapped) parseGappedParams(file,pm_options);
		else parseParams(file);
		cout << "\n Choose \nshape: " << chosenShape << "\n and \nthreshold: " << chosenThreshold<< "\n to achieve optimal performance for expected recognition rate >= " << (100.0-100.0*optionLossRate) << "% (expected recognition = " << (100.0-chosenLossRate*100.0) <<"%)\n\n";
		file.close();
	}
*/
	return 0;
}
