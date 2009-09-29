#include <iostream>
#include <fstream>
#include <seqan/index.h>
#include "test_index_creation.h"

using namespace std;
using namespace seqan;

template < typename TSA, typename TText>
void print_sa(TSA &sa, TText &s) {
	typedef typename Iterator<TSA>::Type TIter;
	int i = 0;
	for (TIter it = begin(sa); it!=end(sa); ++it, ++i) {
// 		cout << i << " " << *it << " " << suffix(s,*it) << endl;
		cout << i << " " << *it << " " << infix(s,*it, min(length(s),*it+100)) << endl;
	}
	cout << endl;
}

template < typename TSA1, typename TSA2 >
bool compare_sa(TSA1 &sa1, TSA2 &sa2) {
	if (length(sa1) != length(sa2)) return false;
	typedef typename Iterator<TSA1>::Type TIter1;
	typedef typename Iterator<TSA2>::Type TIter2;
	TIter1 it1 = begin(sa1);
	TIter2 it2 = begin(sa2);
	int i = 0;
	for (; it1!=end(sa1); ++it1, ++it2, ++i) {
		if (*it1!=*it2) {
			cout << "Mismatch at position " << i << endl;
			return false;
		}
	}
	return true;
}

//template < typename TText, typename TTag, typename TValue, typename TAllowsFastRandomAccess >
//bool check_sa_algorithm(TText& text, TAllowsFastRandomAccess&);

// a suffix array of 'text' is computed using the algorithm indicated by 'tag'
// If 'saReference' is not empty, it is assumed to be the correct suffix array
// and compared against.
template < typename TText, typename TTag, typename TValue>
bool check_sa_algorithm(TText& text, const False&) {
	String<TValue> sa;
	resize(sa, length(text));
	createSuffixArray(sa, text, TTag());
	return isSuffixArray(sa, text);
}

// a suffix array of 'text' is computed using the algorithm indicated by 'tag'
// If 'saReference' is not empty, it is assumed to be the correct suffix array
// and compared against.
template < typename TText, typename TTag, typename TValue>
bool check_sa_algorithm(TText& text, const True&) {
	remove("sa.out");
	String<TValue, External<> > sa("sa.out");
	resize(sa, length(text));
	createSuffixArray(sa, text, TTag());
	bool ok = isSuffixArray(sa, text);
	remove("sa.out");
	return ok;
}

#define MYASSERT(tag, type, use64) if (!check_sa_algorithm<CharString, BwtWalk<tag>, type>(text, use64())) { cerr << "Assertion failed in line " << __LINE__ << ": " << #tag << " " << #type << " " << #use64 << endl; exit(1); }

void testBWTWalk(const char* inputFilename) 
{
#if defined(__GNUC__)
#  if defined(__OPTIMIZE__)
	std::cout << "Optimized build: Yes\n";
#  else
	std::cout << "Optimized build: No\n";
#  endif
#endif

	std::fstream inputFile(inputFilename);
	CharString text;
	read(inputFile, text, Raw());
	std::cout << "filesize: " << length(text) << std::endl;

	MYASSERT(BwtWalkFast, unsigned, False);
	MYASSERT(BwtWalkFast, unsigned, True);
	MYASSERT(BwtWalkFast, __uint64, False);
	MYASSERT(BwtWalkFast, __uint64, True);

	MYASSERT(BwtWalkInPlace, unsigned, False);
	MYASSERT(BwtWalkInPlace, unsigned, True);
	MYASSERT(BwtWalkInPlace, __uint64, False);
	MYASSERT(BwtWalkInPlace, __uint64, True);
}
