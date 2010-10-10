/*==========================================================================
 SeqAn - The Library for Sequence Analysis
 http://www.seqan.de 
 ============================================================================
 Copyright (C) 2007
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.
 
 ============================================================================
 Author: Marcel Martin <marcel.martin@tu-dortmund.de>
 Author: Tobias Marschall <tobias.marschall@tu-dortmund.de>
 Author: David Weese <david.weese@fu-berlin.de>
 ==========================================================================*/

#ifndef TESTS_INDEX_TEST_SA_BWTWALK_H
#define TESTS_INDEX_TEST_SA_BWTWALK_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

template < typename TSA, typename TText>
void print_sa(TSA &sa, TText &s) {
	typedef typename Iterator<TSA>::Type TIter;
	int i = 0;
	for (TIter it = begin(sa); it!=end(sa); ++it, ++i) {
// 		std::cout << i << " " << *it << " " << suffix(s,*it) << std::endl;
		std::cout << i << " " << *it << " " << infix(s,*it, min(length(s),*it+100)) << std::endl;
	}
	std::cout << std::endl;
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
			std::cout << "Mismatch at position " << i << std::endl;
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

#define MYASSERT(tag, type, use64) if (!check_sa_algorithm<CharString, BwtWalk<tag>, type>(text, use64())) { std::cerr << "Assertion failed in line " << __LINE__ << ": " << #tag << " " << #type << " " << #use64 << std::endl; exit(1); }

SEQAN_DEFINE_TEST(testBWTWalk)
{
#if defined(__GNUC__)
#  if defined(__OPTIMIZE__)
	std::cout << "Optimized build: Yes\n";
#  else
	std::cout << "Optimized build: No\n";
#  endif
#endif

    char buffer[1023];
    strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
    strcat(buffer, "/projects/tests/index/m_tuberculosis_h37rv.fa");
    
	std::fstream inputFile(buffer);
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

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
