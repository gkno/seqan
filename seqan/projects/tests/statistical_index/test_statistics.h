/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Jonathan Goeke <goeke@molgen.mpg.de>
  ===========================================================================
  Tests for the random number generation code in seqan/random.
  ===========================================================================
*/

#ifndef TEST_STATISTICAL_INDEX_TEST_STATISTICS_H_
#define TEST_STATISTICAL_INDEX_TEST_STATISTICS_H_

#include <seqan/statistical_index.h>

SEQAN_DEFINE_TEST(test_statistical_index_statistics)
{
	using namespace seqan;

	typedef Dna TAlphabet;
	typedef String<TAlphabet> TSequence;

	TSequence str2 = "CCCAAAGC";
	TSequence str3 = "CCCAAAGTAAATT";

    TSequence str1a =
            "AGAAGCCTCAGATGAGGAGGGTTTTGCTGTGTGCTGCAAGTATCAGGGAGAAAGCATTTCTGCCCTCTCT"
            "GGAACATGGTGTGAACTTCATCCCTGTAATGATATTGTTTGAATTTTCCATGAAAAATTGTCAGCATGAG"
            "AGTAAGAAAAGTGTACGATGGGAAAATATTGAACCAAACAGACAAAAATGGTAGAGTCACATGACCAGTT"
            "TACTCATTGGTAAAGTTAATGAGAGGGTGAGATTAAACAGAAATTGGTAAAGTTAATGAGAGGGTGAGAT"
            "TAAACAGAGGGTGAGATTAAACTTGGGAATGAGTTTGTCTGAGGAGTGAGGTGAAGCATCATTCCTCTGA"
            "TGCACAGGGTAAGGGTTTGTCTGTAAAGAGATAGCACAGGTGTCTGGAGAGCAGCGTGCATGGTAACCTG"
            "TCCTCCAGGCCAGTGGAGCTGTCTGTCTAACCTGGCCAAGGTACAGTCTTCATCAAAGGTCAGGATCCAG"
            "TCCATGCACAAGGGAGGAGCCATTTGCAGCAGAGCCCAGAAATGCCTCCTGCGACATCTTGTTTGTGTCA"
            "TTTACTAGAGTTGGCACTGTCTTAAGATGGGGGCATGGCTGACATTTTCAACTATCATCAGTGAGTCACT"
            "TGCCCAAATGAGGACCATGGTATTAATCTTGCATGTTTTTGGAACTGTTTAAAAAATGTCTGATTTTTGT"
            "TGTTTAGTGTCTGTTTTTGAATTTCCCCTTCTCTGCAGTTCTTGGTTTCTATCTCACTGAGTGCAGAGGA"
            "TTTTAATTGTTGCTGTCTATCTGTGCTTCGCAGCATGAGAGAGCAATGCCTACGGGCTCTTGTGGTGCTT"
            "TGGGGTTGACGGGTTTTATGTCTGAGCAAGCAGATGTCATAGTAGCCATGCTGGATTGCAGTAATAAATG"
            "TGTCCTTTTTTTCCTTCTGTAGCATTGAAAGCCGAAAAGAGAAGAAAGCTGACTCAGGGAAAGGTGTTGA"
            "CAGGGAGACTTGTCTATGACTCGATCTTCAATTTATTTTTTACATATATATGAGAAGAGTGTCACAATTA"
            "TTAATAAAACTGCTTTGATCATGTATTGTAAATTCTGTCCCTCAACCCAAATCCACCTTCATACTGTAAG"
            "TAGTGCAATACTTGTTTCATTTCTGTGTTTAAACTTCTGAGCAGTGAGACATCCCTGTGAGCAGATACAA"
            "TAGCCAATGCAAGAATCTGTGTGTTCCTTGCTGTACGTTAGACATTTGTAAACTGGATTCTGATTGTCAG"
            "TTTTATGAGAGCAATAGCTTCCTTAAAGAGATAAGTCATATTTACCTAGTTTGTATTTTCCTACTTTAGT"
            "GACCTGAAGATGCCTGATAATTTCATTCAGAAGAATTTTTGAAAGGTAGTCTTACTTCTTTTTAGTTTTT"
            "ATAGCTTAGCATTAGTGACTTATTTCAAAAGACCCAAATCAAAAAGTTAGTTTGAAAGCATTTTTTAATA"
            "ATTGTATTTATGCATTTCCTTGATTTAATATGATAAATTTAATACTTAACAATTTATATGTAACTAAAAC"
            "TTAAAGTCATTTGAAAAATATATAGAAACCTATTTACAACTTGTTAAGGACAATCAGACATAATGCAGAG"
            "TTAAGTAGTATTTGCTTAAAATTCAAGTTGTGACTAATGATCAAATACTAGGCTTGTACGAAATGCTTTA"
            "GAAAAACTTTGTAACAGTTTTGTGGGATTTTTCAATATAAACCTTTATCAGAAATATACTAAGTTTGTCT"
            "CCCACTGACAACAGATGTTTTCCAAATAAACATATTCTATACATACTTGTGGAATGCCACATGGTGAATC"
            "ATTGTATATGAAATTCCACTCCTGTACAGTTACTCTGCAGCTAATGGTCATGCACTGCTTAATGCTGGTC"
            "CTGAATCATGTTCTCATGTTAGACCAACAGCTCTCCAATTGTCATTTTTTTTCTGCAGAGTTTTTTTTTT"
            "CCACTTTTAAATTAAATGCATGTTGTGGAAAAACAGTCTTTTAAAATGAAATTTCAGATTCCATTTGAGA"
            "AGGTTCTGTAGATATTTCAGTCCATATAAAATAATACATCTTTACTAAACTTATATAAGGGGAGAGAAAG"
            "TTATGAAGTTTTGGACATTACTAAAAGTACAGTATTTGATTTCACTTTCAATGAATGTGAAGTTAATAAA"
            "ACTAAATCTCATAATGCTCTTGGTTCCTAAGAATGAGTAGTAATCATCAACTTTATAATACTCCAATATT"
            "CCGTTTTATAATAATTCAGAGCCCTGTGGCTTTTACACACCGTTAATTATGTACTCTGTTGGAAGTGCAC"
            "ATGAAAAGTGAAGAAAAGTTCCTCTTGTGATTAAACTAATGGGAGGAAATAAATCAACAAAGTCTCCATT"
            "AAGTTCTACATTTTGAGACCTTTTAAAAATTCCCCTCACAATTCTTTAAGGAGCCCCCCTTTTTATGGAA"
            "CATGAGCCTAAAAATTATAGAAAGAAGAATTTTAAGTTAATAAAGTTTGTATTTATAAATGCTGAAAAAA"
            "TACAGAAACTTTCTGTTCCAAATGTGTTGCTTTGTGTATTTTATAATACAGATACTACATTGTAAACATT"
            "TCCATTGTTTTATGATTTAGCCAGTGATTCCCCAAAGCAGCCTCTTAGTGTTTTAATATATTAATAACTG"
            "TTTTGTTAAAAATGATCATAGTGAATTTAAATCTTCACATGATCACCTATTTGAATAAGCAA";
	TSequence str1b =
            "TGAGGACCTCAGATGAGGAGGGTTTTGCTGTAAACAAGTATCAGGGAGAAAGCATTTCTGCCCTCTCTGC"
            "TGTGTGCTGCAAGTACTTCATCCCTGTAATGATATTGTTTGAATTTTCCATGAAAAATTGTCAGCATGAG"
            "AGTAAGAAAAGTGTACGATGGGAAAATATTGAACCAAACAGACAAAAATGGTAGAGTCACATGACCAGTT"
            "TACTCATTGGTAAAGTTAATGAGAGGGTGAGATTAAACAGAAATTGGTAAAGTTAATGAGAGGGTGAGAT"
            "TAAACAGAGGGTGAGATTAAACTTGGGAATGAGTTTGTCTGAGGAGTGAGGTGAAGCATCATTCCTCTGA"
            "TGCACAGGGTAAGGGTTTGTCTGTAAAGAGATAGCACAGGTGTCTGGAGAGCAGCGTGCATGGTAACCTG"
            "TCCTCCAGGCCAGTGGAGCTGTCTGTCTAACCTGGCCAAGGTACAGTCTTCATCAAAGGTCAGGATCCAG"
            "TCCATGCACAAGGGAGGAGCCATTTGCAGCAGAGCCCAGAAATGCCTCCTGCGACATCTTGTTTGTGTCA"
            "TTTACTAGAGTTGGCACTGTCTTAAGATGGGGGCATGGCTGACATTTTCAACTATCATCAGTGAGTCACT"
            "TGCCCAAATGAGGACCATGGTATTAATCTTGCATGTTTTTGGAACTGTTTAAAAAATGTCTGATTTTTGT"
            "TGTTTAGTGTCTGTTTTTGAATTTCCCCTTCTCTGCAGTTCTTGGTTTCTATCTCACTGAGTGCAGAGGA"
            "TTTTAATTGTTGCTGTCTATCTGTGCTTCGCAGCATGAGAGAGCAATGCCTACGGGCTCTTGTGGTGCTT"
            "TGGGGTTGACGGGTTTTATGTCTGAGCAAGCAGATGTCATAGTAGCCATGCTGGATTGCAGTAATAAATG"
            "TGTCCTTTTTTTCCTTCTGTAGCATTGAAAGCCGAAAAGAGAAGAAAGCTGACTCAGGGAAAGGTGTTGA"
            "CAGGGAGACTTGTCTATGACTCGATCTTCAATTTATTTTTTACATATATATGAGAAGAGTGTCACAATTA"
            "TTAATAAAACTGCTTTGATCATGTATTGTAAATTCTGTCCCTCAACCCAAATCCACCTTCATACTGTAAG"
            "TAGTGCAATACTTGTTTCATTTCTGTGTTTAAACTTCTGAGCAGTGAGACATCCCTGTGAGCATTTGCTG"
            "TAAAGCAAGAATCTGTGTGTTCCTTGCTGTACGTTAGACATTTGTAAACTGGATTCTGATTGTCAGTTTT"
            "ATGAGAGCAATAGCTTCCTTAAAGAGATAAGTCATATTTACCTAGTTTGTATTTTCCTACTTTAGTGACC"
            "TGAAGATGCCTGATAATTTCATTCAGAAGAATTTTTGAAAGGTAGTCTTACTTCTTTTTAGTTTTTATAG"
            "CTTAGCATTAGTGACTTATTTCAAAAGACCCAAATCAAAAAGTTAGTTTGAAAGCATTTTTTAATAATTG"
            "TATTTATGCATTTCCTTGATTTAATATGATAAATTTAATACTTAACAATTTATATGTAACTAAAACTTAA"
            "AGTCATTTGAAAAATATATAGAAACCTATTTACAACTTGTTAAGGACAATCAGACATAATGCAGAGTTAA"
            "GTAGTATTTGCTTAAAATTCAAGTTGTGACTAATGATCAAATACTAGGCTTGTACGAAATGCTTTAGAAA"
            "AACTTTGTAACAGTTTTGTGGGATTTTTCAATATAAACCTTTATCAGAAATATACTAAGTTTGTCTCCCA"
            "CTGACAACAGATGTTTTCCAAATAAACATATTCTATACATACTTGTGGAATGCCACATGGTGAATCATTG"
            "TATATGAAATTCCACTCCTGTACAGTTACTCTGCAGCTAATGGTCATGCACTGCTTAATGCTGGTCCTGA"
            "ATCATGTTCTCATGTTAGACCAACAGCTCTCCAATTGTCATTTTTTTTCTGCAGAGTTTTTTTTTTCCAC"
            "TTTTAAATTAAATGCATGTTGTGGAAAAACAGTCTTTTAAAATGAAATTTCAGATTCCATTTGAGAAGGT"
            "TCTGTAGATATTTCAGTCCATATAAAATAATACATCTTTACTAAACTTATATAAGGGGAGAGAAGTTTAC"
            "AAGGTAGTCTGGGATTACTAACAAAATAAACAAGAGCCTTTCTAGATAAATGTGTCCATATGCCAGTGCG"
            "GTTTAGGTCTTATTCAAGACACAAGTCATTACTT";
	TSequence str1c =
            "CATGGTGTGAACTTCATCCCTGTAATGATATTGTTTGAATTTTCCATGAAAAATTGTCAGCATGAGAGTA"
            "AGAAAAGTGTACGATGGGAAAATATTGAACCAAACAGACAAAAATGGTAGAGTCACATGACCAGTTTACT"
            "CATTGGTAAAGTTAATGAGAGGGTGAGATTAAACAGAAATTGGTAAAGTTAATGAGAGGGTGAGATTAAA"
            "CAGAGGGTGAGATTAAACTTGGGAATGAGTTTGTCTGAGGAGTGAGGTGAAGCATCATTCCTCTGATGCA"
            "CAGGGTAAGGGTTTGTCTGTAAAGAGATAGCACAGGTGTCTGGAGAGCAGCGTGCATGGTAACCTGTCCT"
            "CCAGGCCAGTGGAGCTGTCTGTCTAACCTGGCCAAGGTACAGTCTTCATCAAAGGTCAGGATCCAGTCCA"
            "TGCACAAGGGAGGAGCCATTTGCAGCAGAGCCCAGAAATGCCTCCTGCGAGTTTTGCTGTAAACAAGTAT"
            "CAGGGAGAAAGCATTTCTGCCCTCTCTGCTGTGTGCTGCAAGTACTTCATCCCTGTAATGATATTGTTTG"
            "AATTTTCCATGAAAAATTGTCAGCATGAGAGTAAGAAAAGTGTACGATGGGAAAATATTGAACCAAACAG"
            "ACAAAAATGGTAGAGTCACATGACCAGTTTACTCATTGGTAAAGTTAATGAGAGGGTGAGATTAAACAGA"
            "AATTGGTAAAGTTAATGAGAGGGTGAGATTAAACAGAGGGTGAGATTAAACTTGGGAATGAGTTTGTCTG"
            "AGGAGTGAGGTGAAGCATCATTCCTCTGATGCACAGGGTAAGGGTTTGTCTGTAAAGAGATAGCACAGGT"
            "GTCTGGAGAGCAGCGTGCATGGTAACCTGTCCTCCAGGCCAGTGGAGCTGTCTGTCTAACCTGGCCAAGG"
            "TACAGTCTTCATCAAAGGTCAGGATCCAGTCCATGCACAAGGGAGGAGCCATTTGCAGCAGAGCCCAGAA"
            "ATGCCTCCTGCGACATCTTGTTTGTGTCATTTACTAGAGTTGGCACTGTCTTAAGATGGGGGCATGGCTG"
            "ACATTTTCAACTATCATCAGTGAGTCACTTGCCCAAATGAGGACCATGGTATTAATCTTGCATGTTTTTG"
            "GAACTGTTTAAAAAATGTCTGATTTTTGTTGTTTAGTGTCTGTTTTTGAATTTCCCCTTCTCTGCAGTTC"
            "TTGGTTTCTATCTCACTGAGTGCAGAGGATTTTAATTGTTGCTGTCTATCTGTGCTTCGCAGCATGAGAG"
            "AGCAATGCCTACGGGCTCTTGTGGTGCTTTGGGGTTGACGGGTTTTATGTCTGAGCAAGCAGATGTCATA"
            "GTAGCCATGCTGGATTGCAGTAATAAATGTGTCCTTTTTTTCCTTCTGTAGCATTGAAAGCCGAAAAGAG"
            "AAGAAAGCTGACTCAGGGAAAGGTGTTGACAGGGAGACTTGTCTATGACTCGATCTTCAATTTATTTTTT"
            "ACATATATATGAGAAGAGTGTCACAATTATTAATAAAACTGCTTTGATCATGTATTGTAAATTCTGTCCC"
            "TCAACCCAAATCCACCTTCATACTGTAAGTAGTGCAATACTTGTTTCATTTCTGTGTTTAAACTTCTGAG"
            "CAGTGAGACATCCCTGTGAGCATTTGCTGTAAAGCAAGAATCTGTGTGTTCCTTGCTGTACGTTAGACAT"
            "TTGTAAACTGGATTCTGATTGTCAGTTTTATGAGAGCAATAGCTTCCTTAAAGAGATAAGTCATATTTAC"
            "CTAGTTTGTATTTTCCTACTTTAGTGACCTGAAGATGCCTGATAATTTCATTCAGAAGAATTTTTGAAAG"
            "GTAGTCTTACTTCTTTTTAGTTTTTATAGCTTAGCATTAGTGACTTATTTCAAAAGACCCAAATCAAAAA"
            "GTTAGTTTGAAAGCATTTTTTAATAATTGTATTTATGCATTTCCTTGATTTAATATGATAAATTTAATAC"
            "TTAACAATTTATATGTAACTAAAACTTAAAGTCATTTGAAAAATATATAGAAACCTATTTACAACTTGTT"
            "AAGGACAATCAGACATAATGCAGAGTTAAGTAGTATTTGCTTAAAATTCA"
            ;

	MarkovModel<TAlphabet> mm(3);
	MarkovModel<TAlphabet> mmNew(1);


	char buffer[1024];
	strcpy(buffer, SEQAN_PATH_TO_PROJECTS());
	strcat(buffer, "/projects/tests/statistical_index/zscore_human_mm.3");
	FILE *fd = fopen(buffer,"r");
	read(fd, mm);
	fclose(fd);
	StringSet<TSequence> X;
	appendValue(X, str1a);
    appendValue(X, str1b);
	appendValue(X, str1c);

	StringSet<TSequence> W;
	appendValue(W,str2);
	appendValue(W,str3);
	buildMarkovModel(mmNew,X);

    double x=zscore(W, X, mm, AhoCorasick());
    double x2=zscore(W, X, mmNew, AhoCorasick());

    SEQAN_ASSERT_LEQ(x, 5.9);
    SEQAN_ASSERT_GEQ(x, 5.8);
    SEQAN_ASSERT_LEQ(x2, 8.9);
    SEQAN_ASSERT_GEQ(x2, 8.8);
    double x3=zscore(W, X, mm, WuManber());
    SEQAN_ASSERT_EQ(x, x3);
    // TODO(holtgrew): Porting zscore program as a test with best effort. Nobody knowledgeable was around to verify this.
}

#endif  // TEST_STATISTICAL_INDEX_TEST_STATISTICS_H_
