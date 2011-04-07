// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jonathan Goeke <goeke@molgen.mpg.de>
// ==========================================================================
// Tests for the random number generation code in seqan/random.
// ==========================================================================

#ifndef TEST_STATISTICS_TEST_STATISTICS_H_
#define TEST_STATISTICS_TEST_STATISTICS_H_

#include <seqan/statistics.h>

SEQAN_DEFINE_TEST(test_statistics_statistics)
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
	strcat(buffer, "/projects/tests/statistics/zscore_human_mm.3");
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
    SEQAN_ASSERT_IN_DELTA(x, x3, 0.0001);
    // TODO(holtgrew): Porting zscore program as a test with best effort. Nobody knowledgeable was around to verify this.
}

#endif  // TEST_STATISTICS_TEST_STATISTICS_H_