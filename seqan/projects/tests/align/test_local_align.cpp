#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>


using namespace std;
using namespace seqan;


SEQAN_DEFINE_TEST(testLocalAlign) {
	//align two sequences using Smith-Waterman-algorithm
	String<char> str0 = "ataagcgtctcg";
	String<char> str1 = "tcatagagttgc";
	
	Align< String<char>, ArrayGaps> ali;
	resize(rows(ali), 2);
	setSource(row(ali, 0), str0);
	setSource(row(ali, 1), str1);

	Score<int> score_type = Score<int>(2,-1,-2,0) ;
	LocalAlignmentFinder<int> sw_finder = LocalAlignmentFinder<int>();
	
	int cutoff = 0;

	int score = _smith_waterman(ali,sw_finder,score_type,cutoff);
	
	SEQAN_TASSERT(score == 9);
	SEQAN_TASSERT(row(ali,0) == "ataagcgt");
	SEQAN_TASSERT(row(ali,1) == "ata-gagt");
	
	int i = 1;
	while (true){
		
		score = _smith_waterman_get_next(ali,sw_finder,score_type,cutoff);

		if(score==0){
	//		cout <<"No more alignments satisfying score > "<<cutoff<<"found.\n";		
			break;
		}
		if(i == 1){
			SEQAN_TASSERT(score == 5);
			SEQAN_TASSERT(row(ali,0) == "tc-tcg");
			SEQAN_TASSERT(row(ali,1) == "tcatag");
		}
		if(i == 2){
			SEQAN_TASSERT(score == 4);
			SEQAN_TASSERT(row(ali,0) == "tc");
			SEQAN_TASSERT(row(ali,1) == "tc");
		}
		if(i == 3){
			SEQAN_TASSERT(score == 4);
			SEQAN_TASSERT(row(ali,0) == "gc");
			SEQAN_TASSERT(row(ali,1) == "gc");
		}
		if(i == 4){
			SEQAN_TASSERT(score == 4);
			SEQAN_TASSERT(row(ali,0) == "ag");
			SEQAN_TASSERT(row(ali,1) == "ag");
		}
		if(i == 5){
			SEQAN_TASSERT(score == 4);
			SEQAN_TASSERT(row(ali,0) == "taagcgtctcg");
			SEQAN_TASSERT(row(ali,1) == "tcatagagttg");
		}
		if(i == 6){
			SEQAN_TASSERT(score == 3);
			SEQAN_TASSERT(row(ali,0) == "ata");
			SEQAN_TASSERT(row(ali,1) == "aga");
		}
		if(i == 7){
			SEQAN_TASSERT(score == 3);
			SEQAN_TASSERT(row(ali,0) == "cgt");
			SEQAN_TASSERT(row(ali,1) == "cat");
		}
		if(i == 8){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "g");
			SEQAN_TASSERT(row(ali,1) == "g");
		}
		if(i == 9){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "t");
			SEQAN_TASSERT(row(ali,1) == "t");
		}
		if(i == 10){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "a");
			SEQAN_TASSERT(row(ali,1) == "a");
		}
		if(i == 11){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "c");
			SEQAN_TASSERT(row(ali,1) == "c");
		}
		if(i == 12){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "a");
			SEQAN_TASSERT(row(ali,1) == "a");
		}
		if(i == 13){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "t");
			SEQAN_TASSERT(row(ali,1) == "t");
		}
		if(i == 14){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "c");
			SEQAN_TASSERT(row(ali,1) == "c");
		}
		if(i == 15){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "g");
			SEQAN_TASSERT(row(ali,1) == "g");
		}
		if(i == 16){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "t");
			SEQAN_TASSERT(row(ali,1) == "t");
		}
		if(i == 17){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "a");
			SEQAN_TASSERT(row(ali,1) == "a");
		}
		if(i == 18){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "t");
			SEQAN_TASSERT(row(ali,1) == "t");
		}
		++i;
	}
	
//test if every cell has been reduced to 0 
//only makes sense if cutoff=0
	if(cutoff==0){
		int str0len = length(str0) + 1;
		int str1len = length(str1) + 1;
		bool check = true;
		for(int i = 0; i <str1len; ++i){
			for(int j=0;j<str0len;++j){
				if(getValue(sw_finder.matrix,(i*str0len)+j)!=0){
					check = false;
				}
			}
		}

		SEQAN_TASSERT(check == true);

	}

//desweiteren nur so:
	push(sw_finder.pQ,LocalAlignmentFinder<int>::TPQEntry());
	SEQAN_TASSERT(empty(sw_finder.pQ) == false);
	clear(sw_finder.pQ);
	SEQAN_TASSERT(empty(sw_finder.pQ) == true);




}


SEQAN_DEFINE_TEST(testLocalAlign2) {
//new interface

	String<char> str0 = "ataagcgtctcg";
	String<char> str1 = "tcatagagttgc";
	
	Align< String<char>, ArrayGaps> ali;
	resize(rows(ali), 2);
	setSource(row(ali, 0), str0);
	setSource(row(ali, 1), str1);

	Score<int> score_type = Score<int>(2,-1,-2,0) ;
	LocalAlignmentFinder<int> sw_finder = LocalAlignmentFinder<int>();
	
	int score = localAlignment(ali, sw_finder, score_type, 5);
	SEQAN_TASSERT(score == 9);
	SEQAN_TASSERT(row(ali,0) == "ataagcgt");
	SEQAN_TASSERT(row(ali,1) == "ata-gagt");
	
	score = localAlignment(ali, sw_finder, score_type, 5);
	SEQAN_TASSERT(score == 5);
	SEQAN_TASSERT(row(ali,0) == "tc-tcg");
	SEQAN_TASSERT(row(ali,1) == "tcatag");

	score = localAlignment(ali, sw_finder, score_type, 5, WatermanEggert());
	SEQAN_ASSERT_EQ(score, 0);
}


SEQAN_DEFINE_TEST(testBandedLocalAlign) {
    typedef String<Dna> TString;
	TString str0("ggggcttaagcttgggg");
	TString str1("aaaacttagctctaaaa");

	Score<int> score_type = Score<int>(2,-1,-2,-2);
	
    Align<TString> align;
	resize(rows(align), 2);
	assignSource(row(align, 0), str0);
	assignSource(row(align, 1), str1);

	LocalAlignmentFinder<int> finder = LocalAlignmentFinder<int>();

	int score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 12);
    SEQAN_ASSERT_TRUE(row(align, 0) == "cttaagct");
    SEQAN_ASSERT_TRUE(row(align, 1) == "ctt-agct");

    score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 10);
    SEQAN_ASSERT_TRUE(row(align, 0) == "aagcttgg");
    SEQAN_ASSERT_TRUE(row(align, 1) == "aaacttag");

    score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 10);
    SEQAN_ASSERT_TRUE(row(align, 0) == "gct-taa");
    SEQAN_ASSERT_TRUE(row(align, 1) == "gctctaa");

    score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 5);
    SEQAN_ASSERT_TRUE(row(align, 0) == "aagcttg");
    SEQAN_ASSERT_TRUE(row(align, 1) == "aacttag");

    score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 0);
}
