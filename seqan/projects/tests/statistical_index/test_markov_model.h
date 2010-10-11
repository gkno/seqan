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
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  Tests for the probability distribution code in seqan/random.
  ===========================================================================
*/

#ifndef TEST_STATISTICAL_INDEX_TEST_MARKOV_MODEL_H_
#define TEST_STATISTICAL_INDEX_TEST_MARKOV_MODEL_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/statistical_index.h>

SEQAN_DEFINE_TEST(test_statistical_index_markov_model)
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
            "AAGGACAATCAGACATAATGCAGAGTTAAGTAGTATTTGCTTAAAATTC";

	StringSet<TSequence> X;
	appendValue(X, str1a);
    appendValue(X, str1b);
	appendValue(X, str1c);
	

	StringSet<TSequence> W;
	appendValue(W,str2);
	appendValue(W,str3);

 	MarkovModel<TAlphabet> mm(1);
	buildMarkovModel(mm,X);

    double x = emittedProbability(mm,W);
    SEQAN_ASSERT_IN_DELTA(x, 1.25542e-05, 1.25542e-07);
    double x2 = mm.emittedProbability(W);
    SEQAN_ASSERT_IN_DELTA(x2, 1.25542e-05, 1.25542e-07);

    ensureAuxMatrices(mm);
}

#endif  // TEST_STATISTICAL_INDEX_TEST_MARKOV_MODEL_H_
