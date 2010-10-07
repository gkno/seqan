#ifndef TESTS_ALIGN_TEST_ALIGN_ALIGN_H_
#define TESTS_ALIGN_TEST_ALIGN_ALIGN_H_

#include <iostream>
#include <cstdio>
#include <vector>


#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>

using namespace std;
using namespace seqan;


template <typename TAlign>
void testAlignBasics()
{
	TAlign ali1;			//default C'tor
	TAlign const & c_ali1 = ali1;	//const version

	resize(rows(ali1), 2);	//rows

	//add some rows to ali1
	assignSource(row(ali1, 0), "abcdef");
	insertGaps(row(ali1, 0), 2, 3);
	assignSource(row(ali1, 1), "xyz");

	//test copy C'tor
	TAlign ali2(ali1);
	SEQAN_ASSERT_TRUE(row(ali1, 0) == row(ali2, 0));	//is the same as ali1
	SEQAN_ASSERT_TRUE(row(ali1, 1) == row(ali2, 1));

//  it is a real copy
//	SEQAN_ASSERT_TRUE(id(row(ali1, 0)); == id(row(ali2, 0));) //(its not a real copy)
//	SEQAN_ASSERT_TRUE(id(row(ali1, 1)) == id(row(ali2, 1));)

	TAlign ali3;
	ali3 = ali1;			//operator =
	SEQAN_ASSERT_TRUE(row(ali1, 0) == row(ali3, 0));	//is the same as ali1
	SEQAN_ASSERT_TRUE(row(ali1, 1) == row(ali3, 1));

//  it is a real copy
//	SEQAN_ASSERT_TRUE(id(row(ali1, 0)) == id(row(ali3, 0))) //(its not a real copy)
//	SEQAN_ASSERT_TRUE(id(row(ali1, 1)) == id(row(ali3, 1)))

	detach(ali3);			//detach
	SEQAN_ASSERT_TRUE(row(ali1, 0) == row(ali3, 0));	//is the same as ali1
	SEQAN_ASSERT_TRUE(row(ali1, 1) == row(ali3, 1));
//  it is a real copy
//	SEQAN_ASSERT_TRUE(id(row(ali1, 0)) != id(row(ali3, 0))) //ali3 is not dependent anymore
//	SEQAN_ASSERT_TRUE(id(row(ali1, 1)) != id(row(ali3, 1)))

	SEQAN_ASSERT_NOT(dependentSource(row(ali3, 0))); //dependent
	SEQAN_ASSERT_NOT(dependentSource(row(ali3, 1)));

//____________________________________________________________________________
// Test Column Access
	
	typedef typename Cols<TAlign>::Type TCols;
	typedef typename Cols<TAlign const>::Type TConstCols;
	TCols cols1 = cols(ali1);				//cols
	TConstCols cols2 = cols(c_ali1);		//const version

	typename Col<TAlign>::Type col1 = col(ali1, 0); //col
    SEQAN_ASSERT_TRUE(col1 == cols1[0]);

	typename Col<TAlign const>::Type c_col1 = col(c_ali1, 0);
    SEQAN_ASSERT_TRUE(c_col1 == col1);

//____________________________________________________________________________
// Test Rows Access

	typename Rows<TAlign>::Type rows1 = rows(ali1); //rows
	typename Rows<TAlign const>::Type rows2 = rows(c_ali1); //rows

	typename Row<TAlign>::Type row1 = row(ali1, 0); //row
	typename Row<TAlign>::Type row2 = rows1[0];		//rows::operator[]
	typename Row<TAlign const>::Type row3 = rows2[0];

//____________________________________________________________________________

}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign>
void testAlignColsBase()
{
//____________________________________________________________________________

	TAlign ali1;
	TAlign const & c_ali1 = ali1;
	resize(rows(ali1), 2);
	assignSource(row(ali1, 0), "abcdef");
	insertGaps(row(ali1, 0), 2, 3);
	assignSource(row(ali1, 1), "xyz");

	typedef typename Cols<TAlign>::Type TCols;
	TCols cols1(ali1);							//special c'tor
	TCols const & c_cols1 = cols1;				//(const version)
	TCols cols2;								//default c'tor
	TCols cols3 = cols1;						//copy c'tor
	cols2 = cols1;								//operator =

	typedef typename Cols<TAlign const>::Type TConstCols;
	TConstCols cols4(c_ali1);					//special c'tor

	typename Col<TAlign>::Type col1 = cols1[0];	//operator[]
	typename Col<TAlign>::Type const c_col1 = c_cols1[1];
	typename Col<TAlign const>::Type col2 = cols4[1];

	typedef typename Iterator<TCols, Rooted>::Type TColsIterator;
	typedef typename Iterator<TCols const, Rooted>::Type TConstColsIterator;

	TColsIterator col_it1 = begin(cols1);  //begin
	TConstColsIterator const_col_it1 = begin(c_cols1); 
    SEQAN_ASSERT_TRUE(value(col_it1, 0) == row(ali1, 0)[beginPosition(cols1)]);
    SEQAN_ASSERT_TRUE(value(col_it1, 1) == row(ali1, 1)[beginPosition(cols1)]);
    SEQAN_ASSERT_TRUE(value(const_col_it1, 0) == row(ali1, 0)[beginPosition(cols1)]);
    SEQAN_ASSERT_TRUE(value(const_col_it1, 1) == row(ali1, 1)[beginPosition(cols1)]);

	col_it1 = end(cols1);						//end
	const_col_it1 = end(c_cols1);
	--col_it1;									//operator --
	--const_col_it1;

    SEQAN_ASSERT_TRUE(value(col_it1, 0) == row(ali1, 0)[endPosition(cols1) - 1]);
    SEQAN_ASSERT_TRUE(value(col_it1, 1) == row(ali1, 1)[endPosition(cols1) - 1]);
    SEQAN_ASSERT_TRUE(value(const_col_it1, 0) == row(ali1, 0)[endPosition(c_cols1) - 1]);
    SEQAN_ASSERT_TRUE(value(const_col_it1, 1) == row(ali1, 1)[endPosition(c_cols1) - 1]);

	setBeginPosition(row(ali1, 0), 4);
	setBeginPosition(row(ali1, 1), 2);
	SEQAN_ASSERT_EQ(beginPosition(cols(ali1)), 2u); //beginPosition
//____________________________________________________________________________
// 	cols iterator

	TColsIterator col_it2;						//defaut c'tor
	TColsIterator col_it3(ali1);				//special c'tor

	TColsIterator const & c_col_it1 = col_it1;	// (const version, dont confuse with const_col_it1)
        // Test copying, but do not use further.
	TColsIterator const & c_col_it2 = col_it2;
        (void)c_col_it2;

	SEQAN_ASSERT_TRUE(&host(col_it1) == &ali1);		//host
	SEQAN_ASSERT_TRUE(&host(c_col_it1) == &ali1);

	setHost(col_it2, ali1);						//setHost
	SEQAN_ASSERT_TRUE(&host(col_it2) == &ali1);

	SEQAN_ASSERT_TRUE(&host(container(col_it1)) == &ali1); //container
	SEQAN_ASSERT_TRUE(&host(container(c_col_it1)) == &ali1); //container

	col_it1 = begin(cols1);						//begin
    SEQAN_ASSERT_TRUE(value(col_it1, 0) == row(ali1, 0)[beginPosition(cols1)]);
    SEQAN_ASSERT_TRUE(value(col_it1, 1) == row(ali1, 1)[beginPosition(cols1)]);

	col_it2 = ++col_it1;						//operator ++ (prefix)
    SEQAN_ASSERT_TRUE(value(col_it1, 0) == row(ali1, 0)[beginPosition(cols1) + 1]);
    SEQAN_ASSERT_TRUE(value(col_it1, 1) == row(ali1, 1)[beginPosition(cols1) + 1]);
    SEQAN_ASSERT_TRUE(value(col_it2, 0) == row(ali1, 0)[beginPosition(cols1) + 1]);
    SEQAN_ASSERT_TRUE(value(col_it2, 1) == row(ali1, 1)[beginPosition(cols1) + 1]);

	col_it2 = col_it1++;						//operator ++ (suffix)
    SEQAN_ASSERT_TRUE(value(col_it1, 0) == row(ali1, 0)[beginPosition(cols1) + 2]);
    SEQAN_ASSERT_TRUE(value(col_it1, 1) == row(ali1, 1)[beginPosition(cols1) + 2]);
    SEQAN_ASSERT_TRUE(value(col_it2, 0) == row(ali1, 0)[beginPosition(cols1) + 1]);
    SEQAN_ASSERT_TRUE(value(col_it2, 1) == row(ali1, 1)[beginPosition(cols1) + 1]);

	col_it2 = --col_it1;						//operator -- (prefix)
    SEQAN_ASSERT_TRUE(value(col_it1, 0) == row(ali1, 0)[beginPosition(cols1) + 1]);
    SEQAN_ASSERT_TRUE(value(col_it1, 1) == row(ali1, 1)[beginPosition(cols1) + 1]);
    SEQAN_ASSERT_TRUE(value(col_it2, 0) == row(ali1, 0)[beginPosition(cols1) + 1]);
    SEQAN_ASSERT_TRUE(value(col_it2, 1) == row(ali1, 1)[beginPosition(cols1) + 1]);

	col_it2 = col_it1--;						//operator -- (suffix)
    SEQAN_ASSERT_TRUE(value(col_it1, 0) == row(ali1, 0)[beginPosition(cols1)]);
    SEQAN_ASSERT_TRUE(value(col_it1, 1) == row(ali1, 1)[beginPosition(cols1)]);
    SEQAN_ASSERT_TRUE(value(col_it2, 0) == row(ali1, 0)[beginPosition(cols1) + 1]);
    SEQAN_ASSERT_TRUE(value(col_it2, 1) == row(ali1, 1)[beginPosition(cols1) + 1]);


	SEQAN_ASSERT_NOT(col_it1 == col_it2);			//operator !=
    SEQAN_ASSERT_NOT(c_col_it1 == col_it2) ;
    SEQAN_ASSERT_NOT(col_it2 == c_col_it1) ;
    SEQAN_ASSERT_NOT(c_col_it1 == c_col_it2)   ;
	
	--col_it2;

	SEQAN_ASSERT_TRUE(col_it1 == col_it2);			//operator ==
    SEQAN_ASSERT_TRUE(c_col_it1 == col_it2);
    SEQAN_ASSERT_TRUE(col_it2 == c_col_it1);
    SEQAN_ASSERT_TRUE(c_col_it1 == c_col_it2);


    SEQAN_ASSERT_TRUE(row(ali1, 0) == "ab---cdef");
    SEQAN_ASSERT_EQ(beginPosition(row(ali1, 0)), 4u);
	col_it1 = iter(cols1, 6);					//iter

	SEQAN_ASSERT_TRUE(value(col_it1, 0) == '-');		//value
    SEQAN_ASSERT_TRUE(value(c_col_it1, 0) == '-');

	--col_it1;
	SEQAN_ASSERT_TRUE(value(col_it1, 0) == 'b');		//value
    SEQAN_ASSERT_TRUE(value(c_col_it1, 0) == 'b');

	SEQAN_ASSERT_TRUE(getValue(col_it1, 0) == 'b');	//getValue
    SEQAN_ASSERT_TRUE(getValue(c_col_it1, 0) == 'b');

	char c = 'X';
	assignValue(col_it1, 0, 'x');				//assignValue
    SEQAN_ASSERT_TRUE(row(ali1, 0) == "ax---cdef");
	assignValue(col_it1, 0, c);				
    SEQAN_ASSERT_TRUE(row(ali1, 0) == "aX---cdef");

	c = 'Y';
	assignValue(c_col_it1, 0, 'y');
    SEQAN_ASSERT_TRUE(row(ali1, 0) == "ay---cdef");
	assignValue(c_col_it1, 0, c);
    SEQAN_ASSERT_TRUE(row(ali1, 0) == "aY---cdef");


	c = 'X';
	moveValue(col_it1, 0, 'x');					//moveValue
    SEQAN_ASSERT_TRUE(row(ali1, 0) == "ax---cdef");
	moveValue(col_it1, 0, c);				
    SEQAN_ASSERT_TRUE(row(ali1, 0) == "aX---cdef");

	c = 'Y';
	moveValue(c_col_it1, 0, 'y');
    SEQAN_ASSERT_TRUE(row(ali1, 0) == "ay---cdef");
	moveValue(c_col_it1, 0, c);
    SEQAN_ASSERT_TRUE(row(ali1, 0) == "aY---cdef");


	SEQAN_ASSERT_EQ(length(cols(ali1)), 11u);
}



template <typename TAlign>
void testGotohAlign()
{
	typedef typename Value<TAlign>::Type TValue;
	
	//align two sequences using Smith-Waterman-algorithm
	String<TValue> str0 = "atgt";
	String<TValue> str1 = "atagat";
	
	TAlign ali;
	resize(rows(ali), 2);
	assignSource(row(ali, 0), str0);
	assignSource(row(ali, 1), str1);

	Score<int> score_type = Score<int>(1,-1,-1,-1);
	int score = globalAlignment(ali,score_type);

	SEQAN_ASSERT_TRUE(score == 2);
	SEQAN_ASSERT_TRUE(row(ali,0) == "at-g-t" );
	SEQAN_ASSERT_TRUE(row(ali,1) == "atagat" );
	clearGaps(row(ali,0));
	clearGaps(row(ali,1));

	Score<int> score_type1 = Score<int>(1,-1,-1,-3) ;
	score = globalAlignment(ali,score_type1);

	SEQAN_ASSERT_TRUE(score == -2);
	SEQAN_ASSERT_TRUE((row(ali,0) == "at--gt") || (row(ali,0) == "atg--t" ));
	SEQAN_ASSERT_TRUE(row(ali,1) == "atagat" );
	
	

}


template <typename TAlign>
void testAlignBasics2()
{
	typedef typename Source<TAlign>::Type TSource;
	StringSet<TSource> ss;
	appendValue(ss, "accagtta");
	appendValue(ss, "ccactagggg");

	TAlign aa(ss);
	cout << ss[0]; 
    SEQAN_ASSERT_TRUE(row(aa, 0) == "accagtta");
    SEQAN_ASSERT_TRUE(id(row(aa, 0)) == id(value(ss, 0)));
    SEQAN_ASSERT_TRUE(row(aa, 1) == "ccactagggg");
    SEQAN_ASSERT_TRUE(id(row(aa, 1)) == id(value(ss, 1)));
}


template <typename TAlign>
void _TestAlign()
{
	testAlignBasics<TAlign>();
	testAlignBasics2<TAlign>();
	testAlignColsBase<TAlign>();
	testGotohAlign<TAlign>();
}


SEQAN_DEFINE_TEST(test_align_align_char_string_array_gaps) {
    _TestAlign<Align<String<char>, ArrayGaps> >();
}


SEQAN_DEFINE_TEST(test_align_align_dna_array_gaps) {
    _TestAlign<Align<String<Dna>, ArrayGaps> >();
}


SEQAN_DEFINE_TEST(test_align_align_char_string_sumlist_gaps) {
    _TestAlign<Align<String<char>, SumlistGaps> >();
}


SEQAN_DEFINE_TEST(test_align_align_dna_sumlist_gaps) {
	_TestAlign<Align<String<Dna>, SumlistGaps> >();
}


#endif
