#ifndef TESTS_ALIGN_TEST_ALIGN_GAPS_H_
#define TESTS_ALIGN_TEST_ALIGN_GAPS_H_

#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>


using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSpec>
void TestGapsBase()
{
	typedef Gaps<TSource, TSpec> TGaps;

	TGaps gaps1;					//default ctor
	TSource src1 = "hello";
	setSource(gaps1, src1);			//setSource
	SEQAN_TASSERT(source(gaps1) == src1)
	SEQAN_TASSERT(id(source(gaps1)) == id(src1))

	SEQAN_TASSERT(id(source(gaps1)) == id(gaps1))	//id

	assignSource(gaps1, "blabla");	//assignSource
	SEQAN_TASSERT(source(gaps1) == "blabla")
	SEQAN_TASSERT(src1 == "blabla")

	assignSource(gaps1, "abcdef", 1, 5);	//assignSource
	SEQAN_TASSERT(source(gaps1) == "abcdef")
	SEQAN_TASSERT(sourceSegment(gaps1) == "bcde") //sourceSegment

	moveSource(gaps1, "hullahulla");	//moveSource
	SEQAN_TASSERT(source(gaps1) == "hullahulla")

	moveSource(gaps1, "abcdef", 1, 5);	//moveSource
	SEQAN_TASSERT(source(gaps1) == "abcdef") //???Sollte das anders sein?
	SEQAN_TASSERT(sourceSegment(gaps1) == "bcde")

	detach(gaps1);			//detach, createSource
	SEQAN_TASSERT(source(gaps1) == src1)
	SEQAN_TASSERT(id(gaps1) != id(src1))

	TGaps gaps2(gaps1);				//copy ctor

//  it is a real copy
//	SEQAN_TASSERT(id(gaps1) == id(gaps2)) //(its not a real copy)

//____________________________________________________________________________

	setSource(gaps1, src1);
	src1 = "hello";
	SEQAN_TASSERT(id(gaps1) != id(gaps2))
	gaps2 = gaps1;					//operator =
	SEQAN_TASSERT(id(gaps1) == id(gaps2))
	SEQAN_TASSERT(id(gaps2) == id(src1))

	TGaps gaps3(src1);				//ctor with source
	TGaps const & c_gaps3 = gaps3;	//const version

	SEQAN_TASSERT(id(gaps3) == id(src1))
	SEQAN_TASSERT(id(c_gaps3) == id(src1))

	SEQAN_TASSERT(dependentSource(gaps3))	//dependentSource
	SEQAN_TASSERT(dependentSource(c_gaps3))

	SEQAN_TASSERT(length(gaps3) == length(src1))		//length
	SEQAN_TASSERT(sourceLength(gaps3) == length(src1))	//sourceLength

	SEQAN_TASSERT(sourceBeginPosition(gaps3) == 0)		//sourceBeginPosition
	SEQAN_TASSERT(sourceBegin(gaps3) == begin(source(gaps3)))		//sourceBegin
	SEQAN_TASSERT(sourceBegin(gaps3, Rooted()) == begin(source(gaps3), Rooted()))

	SEQAN_TASSERT(sourceEndPosition(gaps3) == length(src1))	//sourceEndPosition
	SEQAN_TASSERT(sourceEnd(gaps3) == end(source(gaps3)))		//sourceEnd
	SEQAN_TASSERT(sourceEnd(gaps3, Rooted()) == end(source(gaps3), Rooted()))

	SEQAN_TASSERT(*(--end(gaps3)) == 'o') //end
	SEQAN_TASSERT(*(--end(c_gaps3)) == 'o') //end

	setSourceBeginPosition(gaps3, 3); //"---lo"			//setSourceBeginPosition
	SEQAN_TASSERT(gaps3 == "lo");
	SEQAN_TASSERT(sourceBeginPosition(gaps3) == 3)		//sourceBeginPosition
	SEQAN_TASSERT(beginPosition(gaps3) == 3)			//beginPosition
	SEQAN_TASSERT(length(gaps3) == 2)					//length
	
	setSourceBeginPosition(gaps3, 1); //"-ello"			//setSourceBeginPosition
	SEQAN_TASSERT(gaps3 == "ello");
	SEQAN_TASSERT(sourceBeginPosition(gaps3) == 1)		//sourceBeginPosition
	SEQAN_TASSERT(beginPosition(gaps3) == 1)			//beginPosition
	SEQAN_TASSERT(length(gaps3) == 4)					//length
	SEQAN_TASSERT(isGap(gaps3, 0))						//isGap
	SEQAN_TASSERT(gaps3[1] == 'e')
	SEQAN_TASSERT(c_gaps3[1] == 'e')

	SEQAN_TASSERT(getValue(gaps3, 1) == 'e');			//getValue
	SEQAN_TASSERT(getValue(c_gaps3, 1) == 'e');

	setSourceEndPosition(gaps3, 3); //"-el"				//setSourceEndPosition
	SEQAN_TASSERT(gaps3 == "el")
	SEQAN_TASSERT(sourceEndPosition(gaps3) == 3)		//sourceEndPosition
	SEQAN_TASSERT(endPosition(gaps3) == 3)				//endPosition
	SEQAN_TASSERT(length(gaps3) == 2)

	setBeginPosition(gaps3, 0);	//"el"					//setBeginPosition
	SEQAN_TASSERT(gaps3[0] == 'e')
	SEQAN_TASSERT(beginPosition(gaps3) == 0)			
	SEQAN_TASSERT(sourceBeginPosition(gaps3) == 1)
	SEQAN_TASSERT(endPosition(gaps3) == 2)
	SEQAN_TASSERT(sourceEndPosition(gaps3) == 3)


//____________________________________________________________________________

	clear(gaps3);										//clear
	SEQAN_TASSERT(length(gaps3) == 0);	
	SEQAN_TASSERT(sourceBeginPosition(gaps3) == 0)
	SEQAN_TASSERT(sourceEndPosition(gaps3) == 0)

	setSourceEndPosition(gaps3, length(source(gaps3))); //reactivate after clear
	SEQAN_TASSERT(gaps3 == "hello")

	insertGaps(gaps3, 2, 3);
	setBeginPosition(gaps3, 2);
	SEQAN_TASSERT(gaps3 == "he---llo") //"--he---llo"
	SEQAN_TASSERT(beginPosition(gaps3) == 2)

	//toSourcePosition
	SEQAN_TASSERT(toSourcePosition(gaps3, 0) == 0)
	SEQAN_TASSERT(toSourcePosition(gaps3, 1) == 0)
	SEQAN_TASSERT(toSourcePosition(gaps3, 2) == 0)
	SEQAN_TASSERT(toSourcePosition(gaps3, 3) == 1)
	SEQAN_TASSERT(toSourcePosition(gaps3, 4) == 2)
	SEQAN_TASSERT(toSourcePosition(gaps3, 5) == 2)
	SEQAN_TASSERT(toSourcePosition(gaps3, 6) == 2)
	SEQAN_TASSERT(toSourcePosition(gaps3, 7) == 2)
	SEQAN_TASSERT(toSourcePosition(gaps3, 8) == 3)
	SEQAN_TASSERT(toSourcePosition(gaps3, 9) == 4)
	SEQAN_TASSERT(toSourcePosition(gaps3, 10) == 5)
	SEQAN_TASSERT(toSourcePosition(gaps3, 11) == 5)

	//toViewPosition
	SEQAN_TASSERT(toViewPosition(gaps3, 0) == 2)
	SEQAN_TASSERT(toViewPosition(gaps3, 1) == 3)
	SEQAN_TASSERT(toViewPosition(gaps3, 2) == 7)
	SEQAN_TASSERT(toViewPosition(gaps3, 3) == 8)
	SEQAN_TASSERT(toViewPosition(gaps3, 4) == 9)
	SEQAN_TASSERT(toViewPosition(gaps3, 5) == 10)

//____________________________________________________________________________

	SEQAN_TASSERT(gaps3 == "he---llo") //"--he---llo"
	SEQAN_TASSERT(beginPosition(gaps3) == 2)

	clearGaps(gaps3, 1, 5);
	SEQAN_TASSERT(gaps3 == "he--llo") //"-he--llo"
	SEQAN_TASSERT(beginPosition(gaps3) == 1)

	clearGaps(gaps3, 1, 5);
	SEQAN_TASSERT(gaps3 == "hello") //"-hello"
	SEQAN_TASSERT(beginPosition(gaps3) == 1)

	clearGaps(gaps3);
	SEQAN_TASSERT(gaps3 == "hello") //"hello"
	SEQAN_TASSERT(beginPosition(gaps3) == 0)

//____________________________________________________________________________

	SEQAN_TASSERT(gaps3 == "hello") //"hello"

	setSourceBeginPosition(gaps3, 1);
	setSourceEndPosition(gaps3, 3); //"el"
	SEQAN_TASSERT(gaps3 == "el")
	SEQAN_TASSERT(sourceSegment(gaps3) == "el") //sourceSegment
	SEQAN_TASSERT(sourceSegment(c_gaps3) == "el")


//____________________________________________________________________________
// Comparison Functions

	SEQAN_TASSERT(gaps3 == "el") //"hello"
	SEQAN_TASSERT(gaps3 != "ello")
	SEQAN_TASSERT(gaps3 <= "el")
	SEQAN_TASSERT(gaps3 < "ello")
	SEQAN_TASSERT(gaps3 > "a")
	SEQAN_TASSERT(gaps3 >= "el")


//____________________________________________________________________________
}

//////////////////////////////////////////////////////////////////////////////

// Iterator Functions
template <typename TSource, typename TSpec>
void TestGapsIterator()
{
//____________________________________________________________________________
	typedef Gaps<TSource, TSpec> TGaps;


	typedef typename Iterator<TGaps, Rooted>::Type TIterator;

	TSource src1 = "hello";
	TGaps gaps4(src1);

	TIterator it1 = begin(gaps4);						//begin
	TIterator const & c_it1 = it1; //const version
	SEQAN_TASSERT(*it1 == 'h');							//operator *
	SEQAN_TASSERT(source(it1) == begin(src1))			//source
	SEQAN_TASSERT(source(c_it1) == begin(src1))	
	SEQAN_TASSERT(atBegin(c_it1))						//atBegin

	++it1;												//operator ++
	SEQAN_TASSERT(*it1 == 'e');	

	--it1;												//operator --
	SEQAN_TASSERT(*it1 == 'h');	

	TIterator it2 = end(gaps4);							//end
	SEQAN_TASSERT(atEnd(it2))							//atEnd

	--it2;
	SEQAN_TASSERT(*it2 == 'o');
//____________________________________________________________________________


	TIterator it3;										//default ctor
	TIterator it4 = it1;								//copy ctor
	TIterator const & c_it4 = it4;	

	SEQAN_TASSERT(container(it4) == container(it1))
	SEQAN_TASSERT(*it4 == *it1)

	SEQAN_TASSERT(it4 == it1)							//operator ==
	SEQAN_TASSERT(it4 == c_it1)
	SEQAN_TASSERT(c_it4 == it1)
	SEQAN_TASSERT(c_it4 == c_it1)

	++it1;

	SEQAN_TASSERT(it4 != it1)							//operator !=
	SEQAN_TASSERT(it4 != c_it1)
	SEQAN_TASSERT(c_it4 != it1)
	SEQAN_TASSERT(c_it4 != c_it1)

	it4 = it2;											//operator =
	SEQAN_TASSERT(*it4 == *it2)

	TIterator it5(gaps4, 1, 1);							//special ctor
//____________________________________________________________________________



//____________________________________________________________________________
}

//////////////////////////////////////////////////////////////////////////////

// Manipulation of Gaps
template <typename TSource, typename TSpec>
void TestGapManipulation()
{
//____________________________________________________________________________
	typedef Gaps<TSource, TSpec> TGaps;

//inserting gaps
	TSource src1 = "hello";
	TGaps gaps5(src1);

	insertGaps(gaps5, 1, 2); //insert gap somewhere
	SEQAN_TASSERT(gaps5 != "hello");
	SEQAN_TASSERT(gaps5 == "h--ello");

	insertGaps(gaps5, 1, 1); //insert blank at the beginning of a gap
	SEQAN_TASSERT(gaps5 == "h---ello");

	insertGaps(gaps5, 4, 1); //insert blank at the end of a gap
	SEQAN_TASSERT(gaps5 == "h----ello");

	insertGap(gaps5, 8); //insert second gap
	SEQAN_TASSERT(gaps5 == "h----ell-o");

	insertGaps(gaps5, 0, 2); //insert at position 0
	SEQAN_TASSERT(gaps5 == "h----ell-o"); //note: leading and trailing gaps are not displayed
	SEQAN_TASSERT(beginPosition(gaps5) == 2);

	insertGaps(gaps5, 8, 2); //insert gap with beginPosition == 2
	SEQAN_TASSERT(gaps5 == "h----e--ll-o");
	SEQAN_TASSERT(length(gaps5) == 12);

	insertGaps(gaps5, 14, 1); //insert gap behind end. Nothing happens
	SEQAN_TASSERT(gaps5 == "h----e--ll-o");
	SEQAN_TASSERT(length(gaps5) == 12);

//counting gaps
	SEQAN_TASSERT(gaps5 == "h----e--ll-o"); // "--h----e--ll-o"
	SEQAN_TASSERT(beginPosition(gaps5) == 2);

	SEQAN_TASSERT(countGaps(gaps5, 0) == 2);
	SEQAN_TASSERT(countGaps(gaps5, 2) == 0);
	SEQAN_TASSERT(countGaps(gaps5, 3) == 4);
	SEQAN_TASSERT(countGaps(gaps5, 4) == 3);
	SEQAN_TASSERT(countGaps(gaps5, 6) == 1);
	SEQAN_TASSERT(countGaps(gaps5, 8) == 2);
	SEQAN_TASSERT(countGaps(gaps5, 20) == 0);

//removing gaps
	SEQAN_TASSERT(gaps5 == "h----e--ll-o"); // "--h----e--ll-o"
	SEQAN_TASSERT(beginPosition(gaps5) == 2);

	removeGap(gaps5, 4); //remove gap somewhere in gap area
	SEQAN_TASSERT(gaps5 == "h---e--ll-o");

	removeGaps(gaps5, 6, 1); //try to remove gap in non-gap area. Nothing happens
	SEQAN_TASSERT("h---e--ll-o" == gaps5);

	removeGaps(gaps5, 7, 2); //remove gap region completely
	SEQAN_TASSERT(gaps5 == "h---ell-o");

	removeGaps(gaps5, 4, 10); //remove rest of gap region
	SEQAN_TASSERT(gaps5 == "h-ell-o");
//____________________________________________________________________________

}
//////////////////////////////////////////////////////////////////////////////

template <typename TSource, typename TSpec>
void TestSequenceGapsBase()
{
	typedef Gaps<TSource, TSpec> TGaps;

	TGaps gaps3;
	assignSource(gaps3, "hello");

	insertGaps(gaps3, 2, 3);
	setBeginPosition(gaps3, 2);
	SEQAN_TASSERT(gaps3 == "he---llo") //"--he---llo"
	SEQAN_TASSERT(beginPosition(gaps3) == 2)

	//toSourcePosition
	SEQAN_TASSERT(toSourcePosition(gaps3, 0) == 0)
	SEQAN_TASSERT(toSourcePosition(gaps3, 1) == 0)
	SEQAN_TASSERT(toSourcePosition(gaps3, 2) == 0)
	SEQAN_TASSERT(toSourcePosition(gaps3, 3) == 1)
	SEQAN_TASSERT(toSourcePosition(gaps3, 4) == 2)
	SEQAN_TASSERT(toSourcePosition(gaps3, 5) == 2)
	SEQAN_TASSERT(toSourcePosition(gaps3, 6) == 2)
	SEQAN_TASSERT(toSourcePosition(gaps3, 7) == 2)
	SEQAN_TASSERT(toSourcePosition(gaps3, 8) == 3)
	SEQAN_TASSERT(toSourcePosition(gaps3, 9) == 4)
	SEQAN_TASSERT(toSourcePosition(gaps3, 10) == 5)
	SEQAN_TASSERT(toSourcePosition(gaps3, 11) == 5)

	//toViewPosition
	SEQAN_TASSERT(toViewPosition(gaps3, 0) == 2)
	SEQAN_TASSERT(toViewPosition(gaps3, 1) == 3)
	SEQAN_TASSERT(toViewPosition(gaps3, 2) == 7)
	SEQAN_TASSERT(toViewPosition(gaps3, 3) == 8)
	SEQAN_TASSERT(toViewPosition(gaps3, 4) == 9)
	SEQAN_TASSERT(toViewPosition(gaps3, 5) == 10)

//____________________________________________________________________________

	SEQAN_TASSERT(gaps3 == "he---llo") //"--he---llo"
	SEQAN_TASSERT(beginPosition(gaps3) == 2)

	clearGaps(gaps3, 1, 5);
	SEQAN_TASSERT(gaps3 == "he--llo") //"-he--llo"
	SEQAN_TASSERT(beginPosition(gaps3) == 1)

	clearGaps(gaps3, 1, 5);
	SEQAN_TASSERT(gaps3 == "hello") //"-hello"
	SEQAN_TASSERT(beginPosition(gaps3) == 1)

	clearGaps(gaps3);
	SEQAN_TASSERT(gaps3 == "hello") //"hello"
	SEQAN_TASSERT(beginPosition(gaps3) == 0)

	assign(gaps3, "el");
	SEQAN_TASSERT(gaps3 == "el") //"el"

//____________________________________________________________________________
// Comparison Functions

	SEQAN_TASSERT(gaps3 == "el") 
	SEQAN_TASSERT(gaps3 != "ello")
	SEQAN_TASSERT(gaps3 <= "el")
	SEQAN_TASSERT(gaps3 < "ello")
	SEQAN_TASSERT(gaps3 > "a")
	SEQAN_TASSERT(gaps3 >= "el")

}

//____________________________________________________________________________

template <typename TSource, typename TSpec>
void TestTrailingGaps()
{
	typedef Gaps<TSource, TSpec> TGaps;

	//inserting gaps
	TSource src1 = "hello";
	TGaps gaps6(src1);
	SEQAN_TASSERT(gaps6 == "hello");

	insertGaps(gaps6, 10, 2);//somewhere behind the last char
	SEQAN_TASSERT(gaps6 == "hello"); //nothing changes
	SEQAN_TASSERT(countGaps(gaps6, 5) == 0);
	SEQAN_TASSERT(countGaps(gaps6, 10) == 0);

	insertGaps(gaps6, 5, 3); //directly behind the last char: append intended trailing gap
	SEQAN_TASSERT(gaps6 == "hello"); //"hello---"
	SEQAN_TASSERT(countGaps(gaps6,5) == 3);
	SEQAN_TASSERT(countGaps(gaps6,6) == 2);

	insertGap(gaps6, 8);//expand trailing gaps
	SEQAN_TASSERT(gaps6 == "hello"); //"hello----"
	SEQAN_TASSERT(countGaps(gaps6,5) == 4);
	SEQAN_TASSERT(countGaps(gaps6,6) == 3);
	SEQAN_TASSERT(countGaps(gaps6,7) == 2);
	SEQAN_TASSERT(countGaps(gaps6,8) == 1);
	SEQAN_TASSERT(countGaps(gaps6,9) == 0);
	SEQAN_TASSERT(countGaps(gaps6,10) == 0);
	SEQAN_TASSERT(countGaps(gaps6,20) == 0);

	insertGaps(gaps6, 9, 2);//expand trailing gaps on last position
	SEQAN_TASSERT(gaps6 == "hello") //"hello------"
	SEQAN_TASSERT(countGaps(gaps6,5) == 6);

//removing gaps
	SEQAN_TASSERT(gaps6 == "hello"); // "hello------"
	SEQAN_TASSERT(beginPosition(gaps6) == 0);

	removeGap(gaps6, 7); //remove gap somewhere in gap area
	SEQAN_TASSERT(gaps6 == "hello");//"hello-----"
	SEQAN_TASSERT(countGaps(gaps6, 5) == 5);
	SEQAN_TASSERT(countGaps(gaps6, 7) == 3);
	SEQAN_TASSERT(countGaps(gaps6, 9) == 1);
	SEQAN_TASSERT(countGaps(gaps6, 10) == 0);

	removeGaps(gaps6, 8, 15); //remove rest of gap region
	SEQAN_TASSERT(gaps6 == "hello"); //"hello---"
	SEQAN_TASSERT(countGaps(gaps6, 5) == 3);
	SEQAN_TASSERT(countGaps(gaps6, 7) == 1);
	SEQAN_TASSERT(countGaps(gaps6, 8) == 0);

	removeGaps(gaps6, 5, 3); //remove gap region completely
	SEQAN_TASSERT(gaps6 == "hello"); //"hello"
	SEQAN_TASSERT(countGaps(gaps6, 5) == 0);

	insertGaps(gaps6, 5, 6);
	SEQAN_TASSERT(countGaps(gaps6,5) == 6);

	assignSource(gaps6, "new");		//clear trailing gaps when assign new source
	SEQAN_TASSERT(gaps6 == "new")   //"new"
	SEQAN_TASSERT(countGaps(gaps6, 3) == 0);
	insertGaps(gaps6, 3, 10);
	SEQAN_TASSERT(countGaps(gaps6, 3) == 10);

	TSource src2 = "hello";
	setSource(gaps6, src2);			//clear trailing gaps when set new source
	SEQAN_TASSERT(gaps6 == "hello") //"re"
	SEQAN_TASSERT(countGaps(gaps6, 5) == 0);

	insertGaps(gaps6, 5, 10);
	SEQAN_TASSERT(countGaps(gaps6, 5) == 10);

	TGaps gaps6a(gaps6); 		//copy all trailing gaps
	SEQAN_TASSERT(gaps6a == "hello");
	SEQAN_TASSERT(countGaps(gaps6a, 5) == 10);

	TGaps gaps6b = gaps6a;		//assign the trailing gaps
	insertGaps(gaps6b, 10, 5);
	SEQAN_TASSERT(gaps6b == "hello");
	SEQAN_TASSERT(countGaps(gaps6a, 5) == 10);
	SEQAN_TASSERT(countGaps(gaps6b, 5) == 15);

	clearGaps(gaps6a);				//remove all trailing gaps
	SEQAN_TASSERT(countGaps(gaps6a, 5) == 0);
	SEQAN_TASSERT(countGaps(gaps6b, 5) == 15);

	}

template <typename TSource, typename TGapSpec>
inline void
TestCountCharacters() {
	typedef Gaps<TSource, TGapSpec> TGap;

	TSource seq = "hello";
	TGap gap7(seq);

	SEQAN_TASSERT(gap7 == "hello");
	SEQAN_TASSERT(length(gap7) == 5);

	SEQAN_TASSERT(countCharacters(gap7, 0) == 5);
	SEQAN_TASSERT(countCharacters(gap7, 1) == 4);
	SEQAN_TASSERT(countCharacters(gap7, 0) == 5);
	SEQAN_TASSERT(countCharacters(gap7, 3) == 2);
	SEQAN_TASSERT(countCharacters(gap7, 5) == 0);

	insertGaps(gap7, 0, 2);
	SEQAN_TASSERT(gap7 == "hello"); //"--hello"
	SEQAN_TASSERT(length(gap7) == 5);

	SEQAN_TASSERT(countCharacters(gap7, 0) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 2) == 5);

	insertGap(gap7, 4);
	SEQAN_TASSERT(gap7 == "he-llo"); //"--he-llo"
	SEQAN_TASSERT(length(gap7) == 6);

	SEQAN_TASSERT(countCharacters(gap7, 0) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 2) == 2);
	SEQAN_TASSERT(countCharacters(gap7, 3) == 1);
	SEQAN_TASSERT(countCharacters(gap7, 4) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 5) == 3);
	SEQAN_TASSERT(countCharacters(gap7, 6) == 2);
	SEQAN_TASSERT(countCharacters(gap7, 7) == 1);
	SEQAN_TASSERT(countCharacters(gap7, 8) == 0);

	insertGaps(gap7, 6, 3);
	SEQAN_TASSERT(gap7 == "he-l---lo"); //"--he-l---lo"
	SEQAN_TASSERT(length(gap7) == 9);

	SEQAN_TASSERT(countCharacters(gap7, 0) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 2) == 2);
	SEQAN_TASSERT(countCharacters(gap7, 3) == 1);
	SEQAN_TASSERT(countCharacters(gap7, 4) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 5) == 1);
	SEQAN_TASSERT(countCharacters(gap7, 6) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 7) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 8) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 9) == 2);
	SEQAN_TASSERT(countCharacters(gap7, 10) == 1);

	insertGaps(gap7, 11, 5);
	SEQAN_TASSERT(gap7 == "he-l---lo"); //"--he-l---lo-----"
	SEQAN_TASSERT(length(gap7) == 9);
	SEQAN_TASSERT(countCharacters(gap7, 11) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 12) == 0);
	SEQAN_TASSERT(countCharacters(gap7, 30) == 0);
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_align_gaps_base_char_string_array_gaps) {
    TestGapsBase<String<char>, ArrayGaps>();
}


SEQAN_DEFINE_TEST(test_align_gaps_base_char_string_sumlist_gaps) {
    TestGapsBase<String<char>, SumlistGaps>();
}


SEQAN_DEFINE_TEST(test_align_gaps_test_gaps_iterator) {
    TestGapsIterator<String<char>, ArrayGaps>();
}


SEQAN_DEFINE_TEST(test_align_gaps_test_gap_manipulation_char_string_array_gaps) {
	TestGapManipulation<String<char>, ArrayGaps>();
}


SEQAN_DEFINE_TEST(test_align_gaps_test_gap_manipulation_char_string_sumlist_gaps) {
    TestGapManipulation<String<char>, SumlistGaps>(); 
}


SEQAN_DEFINE_TEST(test_align_gaps_test_sequence_gaps_base) {
    TestSequenceGapsBase<String<char>, SequenceGaps>();
}

SEQAN_DEFINE_TEST(test_align_gaps_test_trailing_gaps_char_string_array_gaps) {
	TestTrailingGaps<String<char>, ArrayGaps>();
}

SEQAN_DEFINE_TEST(test_align_gaps_test_count_characters_char_string_array_gaps) {
	TestCountCharacters<String<char>, ArrayGaps >();
}


#endif
