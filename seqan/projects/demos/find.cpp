#include <iostream>
#include <fstream>
#include <seqan/sequence.h>
#include <seqan/score.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <time.h>

using namespace seqan;
using namespace std;

template <typename TAlgorithmSpec, typename TAlphabet>
void _findNeedle(String<TAlphabet>& haystack, String<TAlphabet>& needle, String<unsigned int>& pos)
{
	clear(pos);
	Finder<String<char> > finder(haystack);
	Pattern<String<TAlphabet>, TAlgorithmSpec> pattern(needle);
	while (find(finder, pattern)) {
		append(pos,position(finder));
	}
}

void _find(char const* path, String<char> ndl) {
	String<char> haystack;
	fstream strm_in;
	strm_in.open(path, ios_base::in | ios_base::binary);
	read(strm_in, haystack, Fasta());
	strm_in.close();
	
	String<unsigned int> posHorspool;
	String<unsigned int> posShiftAnd;
	String<unsigned int> posShiftOr;
	String<unsigned int> posBomAlgo;
	String<unsigned int> posBndmAlgo;

	_findNeedle<Horspool>(haystack, ndl, posHorspool);
	_findNeedle<ShiftAnd>(haystack, ndl, posShiftAnd);
	_findNeedle<ShiftOr>(haystack, ndl, posShiftOr);
	_findNeedle<BomAlgo>(haystack, ndl, posBomAlgo);
	_findNeedle<BndmAlgo>(haystack, ndl, posBndmAlgo);

	if ((posShiftAnd == posHorspool) &&
		(posShiftAnd == posShiftOr) &&
		(posShiftAnd == posBndmAlgo) &&
		(posShiftAnd == posBomAlgo))
	{
		std::cout << "All Finder Match!" << ::std::endl;
		std::cout << "Number of pattern occurences: " << length(posShiftAnd) << ::std::endl;
	} else {
		std::cout << "Mismatch!!!" << ::std::endl;
		exit(0);
	}
}

void testFinder(char const * path, char const * index_file_name, char const * pattern)
{
	char filepath[1024];
	strcpy(filepath, path);
	char * filename = filepath + strlen(path);

	strcpy(filename, index_file_name);

	fstream strm;
	strm.open(filepath, ios_base::in);
	while (!strm.eof())
	{
		strm.getline(filename, 512);
		cout << "______________________\n" << filename << "\n";
		_find(filepath, pattern);
	}
	strm.close();
}


template <typename TAlgorithmSpec, typename TAlphabet>
inline unsigned int
_measureTime(String<TAlphabet>& haystack, String<TAlphabet>& needle)
{
	unsigned int count = 0;
	Finder<String<TAlphabet> > finder(haystack);
	Pattern<String<TAlphabet>, TAlgorithmSpec>  pattern(needle);
	while (find(finder, pattern)) {
		++count;
	}
	return count;
}

void performanceTest(char const * path) {
	typedef Dna TAlphabet;

	
	String<TAlphabet> haystack;
	fstream strm_in;
	strm_in.open(path, ios_base::in | ios_base::binary);
	read(strm_in, haystack, Fasta());
	strm_in.close();

	String<TAlphabet> needle("AG");
	while(length(needle) < 128) {
		time_t startTime;
		time_t duration;
		unsigned int count;
		startTime = time(0);
		count = _measureTime<Horspool, TAlphabet>(haystack, needle); 
		duration = time(0) - startTime;
		std::cout << "Horspool;" << length(needle) << ";" << count << ";" << duration << ::std::endl;
		startTime = time(0);
		count = _measureTime<ShiftAnd, TAlphabet>(haystack, needle); 
		duration = time(0) - startTime;
		std::cout << "ShiftAnd;" << length(needle) << ";" << count << ";" << duration << ::std::endl;
		startTime = time(0);
		count = _measureTime<ShiftOr, TAlphabet>(haystack, needle); 
		duration = time(0) - startTime;
		std::cout << "ShiftOr;" << length(needle) << ";" << count << ";" << duration << ::std::endl;
		startTime = time(0);
		count = _measureTime<BndmAlgo, TAlphabet>(haystack, needle); 
		duration = time(0) - startTime;
		std::cout << "BndmAlgo;" << length(needle) << ";" << count << ";" << duration << ::std::endl;
		startTime = time(0);
		count = _measureTime<BomAlgo, TAlphabet>(haystack, needle); 
		duration = time(0) - startTime;
		std::cout << "BomAlgo;" << length(needle) << ";" << count << ";" << duration << ::std::endl;
		appendValue(needle, 'C');
		appendValue(needle, 'T');
	}
}


int main()
{
	testFinder("projects/apps/aav/sequences/", "dir.txt", "AA"); 
	//performanceTest("fungi.fasta2"); 
}

