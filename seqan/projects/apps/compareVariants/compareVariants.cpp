#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>


#include "seqan/platform.h"
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/misc/misc_parsing.h>
#include <seqan/refinement.h>

#include "compareVariants.h"

using namespace std;
using namespace seqan;




template<typename TChar>
inline bool
parse_isDigit(TChar const c)
{
	//return (((int ) c >=  48) && ((int) c <=  57));
	return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
		    (c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
}


template<typename TFile, typename TChar>
inline long double
parse_readDouble(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!parse_isDigit(c) && (c != '.')) break;
		append(str, c);
	}
 	return atof(toCString(str));
}

template<typename TFile, typename TChar>
inline void 
parse_skipWhitespace(TFile& file, TChar& c)
{
	if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) break;
	}
}

template<typename TFile, typename TChar, typename TString>
void
_parseReadWordUntilWhitespace(TFile& file, TString& str, TChar& c)
{
        append(str,c);
        if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
                c = _streamGet(file);
                return;
        }
        while (!_streamEOF(file)) {
                c = _streamGet(file);
                if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
                append(str, c);
        }
        return;
}   


template <typename TOptions>
bool loadGenomes(const char* fileName, 
		StringSet<Dna5String> &genomes,
		String<CharString> & genomeIDs, 
		::std::map<CharString,unsigned> &gIdStringToIdNumMap,
		TOptions &options)
{
	
	MultiFasta multiFasta;
	if (!open(multiFasta.concat,fileName,OPEN_RDONLY)) return false;
	split(multiFasta, Fasta());

	unsigned seqCount = length(multiFasta);
	resize(genomeIDs,seqCount);
	if(options.sequenceContext) resize(genomes,seqCount);
	
	for(unsigned i = 0; i < seqCount; ++i)
	{
		CharString temp;
		if(options.sequenceContext) assignSeq(genomes[i], multiFasta[i], Fasta());
		assignSeqId(temp, multiFasta[i], Fasta());
		for (unsigned pos = 0; pos < length(temp); ++pos)
		{
			if(temp[pos]=='\t' || temp[pos]=='\b' || temp[pos]==' ')
			{
				resize(temp,pos);
				break;
			}
		}
		genomeIDs[i] = temp;
		gIdStringToIdNumMap.insert(::std::make_pair<CharString,unsigned>(temp,i)); 
	}
	return (seqCount > 0);
}



/////////////////////////////////////////////////////////////
// read Gff input file containing indels
template <
	typename TIndelSet,
	typename TGenomeMap,
	typename TOptions
>
int readGFF(
	const char*				&filename,
	TIndelSet 				&indelSet,
	TGenomeMap				&gIdStringToIdNumMap,
	TOptions				&options)
{
	typedef typename Value<TIndelSet>::Type	TIndel;
	typedef int				TId;
	typedef int				TContigPos;
	
	
	::std::ifstream file;
	file.open(filename, ::std::ios_base::in | ::std::ios_base::binary);
	if (!file.is_open()) return 1;
		
	TIndel indel = {0,0,0,0,0,0};
	
	clear(indelSet);
	char c = _streamGet(file);
	while (!_streamEOF(file))
	{
		
		if(c == '#')
			_parse_skipLine(file,c);	
				
		// skip whitespaces just in case (actually there shouldnt be a whitespace at the beginning of a line)
		_parseSkipWhitespace(file, c);
	
		if(c == '#')
			_parse_skipLine(file,c);	
		// and read entry in column 1  --> genomeID
		CharString temp_str;
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		
		TId contigId;
		//check if the genomeID is in our map of relevant genomeIDs, otherwise skip match
		typename TGenomeMap::iterator it = gIdStringToIdNumMap.find(temp_str);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		if(it != gIdStringToIdNumMap.end()) contigId = it->second;
		else
		{
		//	std::cout  << "No1\n";
			_parse_skipLine(file,c);
			continue;
		}
		
		// skip whitespaces and read entry in column 2
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		
		// skip whitespaces and read entry in column 3
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		_parseReadWordUntilWhitespace(file,temp_str,c); 
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\t";
		
		// skip whitespaces and read entry in column 4  --> genomic begin position
		_parseSkipWhitespace(file, c);
		indel.originalPos = (TContigPos) _parseReadNumber(file,c) - 1;
		if(options._debugLevel > 1) 
			::std::cout << indel.originalPos << "\t";
		
		// skip whitespaces and read entry in column 5  --> genomic end position // not needed here
		_parseSkipWhitespace(file, c);
		_parseReadNumber(file,c);
		
		// skip whitespaces and read entry in column 6  --> score (percent identity or mapping quality) or a '.'
		int readSupport = 1000; //  --> no information about read support (reference indel)
		_parseSkipWhitespace(file, c);
		if(c=='.')
			c = _streamGet(file);               // 
		else 
			readSupport = (TContigPos) _parseReadDouble(file,c); // number of supporting reads
			
		if(options._debugLevel > 1) 
			::std::cout << readSupport << "\t";
		
		// skip whitespaces and read entry in column 7  --> strand information: '+' or '-' // not needed here
		_parseSkipWhitespace(file, c);
		c = _streamGet(file);
		
		// skip whitespaces and read entry in column 8  --> always '.' here
		_parseSkipWhitespace(file, c);
		c = _streamGet(file);
		
		// skip whitespaces and read entry in column 9  --> tags, extra information. first tag is always "ID"
		_parseSkipWhitespace(file, c);
		clear(temp_str);
		_parseReadIdentifier(file,temp_str,c);
		if(options._debugLevel > 1) 
			::std::cout << temp_str << "\n";
		if(temp_str!="ID") ::std::cout << "first feature field should be 'ID'"<<::std::endl;
		
		// skip the "="
		c = _streamGet(file);
		
		// read the ID
		clear(temp_str);
		clear(indel.idStr);
		CharString indelID;
		_parseReadIdentifier(file,indel.idStr,c);
		if(options._debugLevel > 1) 
			::std::cout << "myID = "<< indel.idStr << "\n";
		
		// process tags in a loop
		CharString current_tag;
		_parseSkipWhitespace(file,c); 
		while(!_streamEOF(file) && !(c == '\n' || (c == '\r' && _streamPeek(file) != '\n'))) // while in same line
		{
			// different tags are separated by ';'  
			while(c != ';')
			{
				if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
					break;
				c = _streamGet(file);
			}
			if(c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) // end of line
				break;
			
			// get the current tag
			clear(current_tag);
			c = _streamGet(file);
			_parseReadIdentifier(file,current_tag,c);
			if(options._debugLevel > 1) 
				::std::cout << current_tag << " in features\n";
			if(current_tag=="size")
			{
		//		if(c == '-') 
					c = _streamGet(file);
				indel.indelSize = _parseReadNumber(file,c); 
				if(options._debugLevel > 1) 
					::std::cout << indel.indelSize << " indel size\n";
			}
			else
			{
				if(current_tag=="duplication") indel.duplication = true;
		//		else _parse_skipLine(file,c);
			}
		}
		appendValue(indelSet,indel);
		_parseSkipWhitespace(file, c);
		
	}

	file.close();
	
	if(options._debugLevel > 0) ::std::cout << "Parsed "<<length(indelSet)<<" indels." << ::std::endl;
	
	return 0;
}



//////////////////////////////////////////////////////////////////////////////
// Print usage
template<typename TOptions>
void printHelp(int, const char *[],TOptions &, bool longHelp = false) 
{
	cerr << "***********************" << endl;
	cerr << "*** compareVariants ***" << endl;
	cerr << "***********************" << endl << endl;
	cerr << "Usage: compareVariants [OPTIONS]... <SOURCE SEQUENCE FILE>" << endl;
	cerr << "\n";
	if (longHelp) {
		cerr << "  -ip,  --input-predicted FILE     \t" << "input gff file containing predicted indels" << endl;
		cerr << "  -ir,  --input-reference FILE     \t" << "input gff file containing reference indels" << endl;
		cerr << "  -o,   --output FILE              \t" << "output filename" << endl;
		cerr << "  -pt,  --position-tolerance NUM   \t" << "position tolerance in bp" << endl;
		cerr << "  -st,  --size-tolerance NUM       \t" << "size tolerance in percent" << endl;
		cerr << "  -sc,  --sequence-context         \t" << "switch on sequence-context mode" << endl;
		cerr << "  -v,   --verbose                  \t" << "verbose mode" << endl;
//		cerr << "  -vv,  --very-verbose             \t" << "very verbose mode" << endl;
		cerr << "  -h,   --help                     \t" << "print this help" << endl;
	} else {
		cerr << "Try 'compareVariants --help' for more information." << endl;
	}
}


//////////////////////////////////////////////////////////////////////////////
// Main part
int main(int argc, const char *argv[])
{
	srand(time(NULL));

	unsigned fnameCount = 0;
	const char *fname[1] = {""};
	IndelCompareOptions<> options;

	// Command line parsing
	for(int arg = 1; arg < argc; ++arg) {
		if (argv[arg][0] == '-') {
			// parse option

			if (strcmp(argv[arg], "-ip") == 0 || strcmp(argv[arg], "--input-predicted") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.inputPredicted = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-ir") == 0 || strcmp(argv[arg], "--input-reference") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.inputReference = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-pt") == 0 || strcmp(argv[arg], "--position-tolerance") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.positionTolerance;
					if (!istr.fail())
					{
						if (options.positionTolerance < 0)
							cerr << "PositionTolerance must be a positive integer value" << endl << endl;
						else
							continue;
					}
				}
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-st") == 0 || strcmp(argv[arg], "--size-tolerance") == 0) {
				if (arg + 1 < argc) {
					++arg;
					istringstream istr(argv[arg]);
					istr >> options.sizeTolerance;
					if (!istr.fail())
					{
						if (options.sizeTolerance < 0 || options.sizeTolerance > 100)
							cerr << "SizeTolerance must be a value between 0 and 100" << endl << endl;
						else{
							options.sizeTolerance /= 100;
							continue;
						}
					}
				}
				printHelp(argc, argv, options);
				return 0;
			}
			if (strcmp(argv[arg], "-o") == 0 || strcmp(argv[arg], "--output") == 0) {
				if (arg + 1 == argc) {
					printHelp(argc, argv, options);
					return 0;
				}
				++arg;
				options.output = argv[arg];
				continue;
			}
			if (strcmp(argv[arg], "-sc") == 0 || strcmp(argv[arg], "--sequence-context") == 0) {
				options.sequenceContext = true;
				continue;
			}
			if (strcmp(argv[arg], "-v") == 0 || strcmp(argv[arg], "--verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 1);
				continue;
			}
			if (strcmp(argv[arg], "-vv") == 0 || strcmp(argv[arg], "--very-verbose") == 0) {
				options._debugLevel = max(options._debugLevel, 2);
				continue;
			}
			if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0) {
				// print help
				printHelp(argc, argv, options, true);
				return 0;
			}
		}
		else {
			// parse file name
			if (fnameCount == 1) {
				printHelp(argc, argv, options);
				return 0;
			}
			fname[fnameCount++] = argv[arg];
		}
	}
	if (fnameCount < 1) {
		printHelp(argc, argv, options);
		return 0;
	}
	
	StringSet<IndelInfo>		refIndels, predictedIndels;
	
	::std::map<CharString,unsigned> gIdStringToIdNumMap;
	String<CharString> genomeIDs;
	StringSet<Dna5String> genomes;
	
	loadGenomes(fname[0],genomes,genomeIDs,gIdStringToIdNumMap,options);
	
	if (readGFF(options.inputReference, refIndels, gIdStringToIdNumMap, options) > 0) 
	{
		cerr << "Reference indels " << options.inputReference << " can't be loaded." << endl;
		return 0;
	}
	if (readGFF(options.inputPredicted, predictedIndels, gIdStringToIdNumMap, options) > 0) 
	{
		cerr << "Predicted indels " << options.inputPredicted << " can't be loaded." << endl;
		return 0;
	}
	
	if(options._debugLevel > 0 )
	{
		::std::cout << "Number of reference indels: " << length(refIndels) << endl;
		::std::cout << "Number of predicted indels: " << length(predictedIndels) << endl;
	}
	
	int result = compareIndels(refIndels,predictedIndels,genomes,genomeIDs,options);
	if(result > 0)
	{
		cerr << "Something went wrong.. Exiting..\n";
		return 1;
	}

	return 0;
}
