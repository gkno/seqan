// ==========================================================================
//                               alignmentFree
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Jonathan Göke <goeke@molgen.mpg.de>
// Tests for all alignment free sequence comparison methods implemented in 
// Seqan
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/alignment_free.h>

template <typename TStringSet>
void alfTestHelperGetSequences(TStringSet & sequences)
{
    using namespace seqan;
    CharString seqIID1 =
        "TAGGTTTTCCGAAAAGGTAGCAACTTTACGTGATCAAACCTCTGACGGGGTTTTCCCCGTCGAAATTGGGTG"
        "TTTCTTGTCTTGTTCTCACTTGGGGCATCTCCGTCAAGCCAAGAAAGTGCTCCCTGGATTCTGTTGCTAACG"
        "AGTCTCCTCTGCATTCCTGCTTGACTGATTGGGCGGACGGGGTGTCCACCTGACGCTGAGTATCGCCGTCAC"
        "GGTGCCACATGTCTTATCTATTCAGGGATCAGAATTTATTCAGGAAATCAGGAGATGCTACACTTGGGTTAT"
        "CGAAGCTCCTTCCAAGGCGTAGCAAGGGCGACTGAGCGCGTAAGCTCTAGATCTCCTCGTGTTGCAACTACA"
        "CGCGCGGGTCACTCGAAACACATAGTATGAACTTAACGACTGCTCGTACTGAACAATGCTGAGGCAGAAGAT"
        "CGCAGACCAGGCATCCCACTGCTTGAAAAAACTATNNNNCTACCCGCCTTTTTATTATCTCATCAGATCAAG";
    CharString seqIID2 =
        "ACCGACGATTAGCTTTGTCCGAGTTACAACGGTTCAATAATACAAAGGATGGCATAAACCCATTTGTGTGAA"
        "AGTGCCCATCACATTATGATTCTGTCTACTATGGTTAATTCCCAATATACTCTCGAAAAGAGGGTATGCTCC"
        "CACGGCCATTTACGTCACTAAAAGATAAGATTGCTCAAANNNNNNNNNACTGCCAACTTGCTGGTAGCTTCA"
        "GGGGTTGTCCACAGCGGGGGGTCGTATGCCTTTGTGGTATACCTTACTAGCCGCGCCATGGTGCCTAAGAAT"
        "GAAGTAAAACAATTGATGTGAGACTCGACAGCCAGGCTTCGCGCTAAGGACGCAAAGAAATTCCCTACATCA"
        "GACGGCCGCGNNNAACGATGCTATCGGTTAGGACATTGTGCCCTAGTATGTACATGCCTAATACAATTGGAT"
        "CAAACGTTATTCCCACACACGGGTAGAAGAACNNNNATTACCCGTAGGCACTCCCCGATTCAAGTAGCCGCG";

    clear(sequences);
    appendValue(sequences, seqIID1);
    appendValue(sequences, seqIID2);
}

SEQAN_DEFINE_TEST(test_alignment_free_calculateProbability)
{
    using namespace seqan;
    double p = 0.0;
    DnaString word = "CCCAAGTTT";
    String<double> model;
    resize(model, 4);
    model[0] = 0.3;  // p(A)
    model[1] = 0.2;  // p(C)
    model[2] = 0.2;  // p(G)
    model[3] = 0.3;  // p(T)
    calculateProbability(p, word, model);
    SEQAN_ASSERT_IN_DELTA(p, 0.00000387, 0.000001);
}

SEQAN_DEFINE_TEST(test_alignment_free_calculateVariance)
{
    using namespace seqan;
    double var = 0.0;
    int n = 10000;
    DnaString word = "CAAGTC";
    String<double> model;
    resize(model, 4);
    model[0] = 0.3;  // p(A)
    model[1] = 0.2;  // p(C)
    model[2] = 0.2;  // p(G)
    model[3] = 0.3;  // p(T)
    calculateVariance(var, word, model, n);  // var = 2.16
    SEQAN_ASSERT_IN_DELTA(var, 2.15845, 0.001);
    StringSet<DnaString> sequences;
    appendValue(sequences, "CAGAAAAAAACACTGATTAACAGGAATAAGCAGTTTACTTATTTTGGGCCTGGGACCCGTGTCTCTAATTTAATTAGGTGATCCCTGCGAAGTTTCTCCA");
    MarkovModel<Dna, double> modelMM0(0);  // Bernoulli model
    modelMM0.build(sequences);
    calculateVariance(var, word, modelMM0, n);  // var = 2.16
    SEQAN_ASSERT_IN_DELTA(var, 2.15845, 0.001);
    MarkovModel<Dna, double> modelMM1(1);  // First order Markov model
    modelMM1.build(sequences);
    calculateVariance(var, word, modelMM1, n);  // var = 1.69716
    SEQAN_ASSERT_IN_DELTA(var, 1.69716, 0.001);
}

SEQAN_DEFINE_TEST(test_alignment_free_calculateCovariance)
{
    using namespace seqan;
    double covar = 0.0;
    int n = 10000;
    DnaString word1 = "ATATAT";
    DnaString word2 = "TATATA";
    String<double> model;
    resize(model, 4);
    model[0] = 0.3;  // p(A)
    model[1] = 0.2;  // p(C)
    model[2] = 0.2;  // p(G)
    model[3] = 0.3;  // p(T)
    calculateCovariance(covar, word1, word2, model, n);  // covar = 4.74
    SEQAN_ASSERT_IN_DELTA(covar, 4.741, 0.001);
    StringSet<DnaString> sequences;
    appendValue(sequences, "CAGCACTGATTAACAGGAATAAGCAGTTTACTTCTGTCAGAATATTGGGCATATATACTGGGACCCGTGTAATACTCTAATTTAATTAGGTGATCCCTGCGAAGTCTCCA");
    MarkovModel<Dna, double> modelMM0(0);  // Bernoulli model
    modelMM0.build(sequences);
    calculateCovariance(covar, word1, word2, modelMM0, n);  // covar = 4.74
    SEQAN_ASSERT_IN_DELTA(covar, 4.741, 0.001);
    MarkovModel<Dna, double> modelMM1(1);  // First order Markov model
    modelMM1.build(sequences);
    calculateCovariance(covar, word1, word2, modelMM1, n);  // covar = 4.74
    SEQAN_ASSERT_IN_DELTA(covar, 13.1541, 0.001);
}
SEQAN_DEFINE_TEST(test_alignment_free_d2_dna)
{
    using namespace seqan;
    StringSet<DnaString> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    // kmerSize = 3
    unsigned kmerSize = 3;
    AFScore<D2> myScoreD2(kmerSize);
    alignmentFreeComparison(myMatrix, sequences, myScoreD2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 4424.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 3965.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 3965.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 4814.0, 0.01);
}

SEQAN_DEFINE_TEST(test_alignment_free_d2_dna5)
{
    using namespace seqan;
    StringSet<Dna5String> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    // kmerSize = 3
    unsigned kmerSize = 3;
    AFScore<D2> myScoreD2(kmerSize);
    alignmentFreeComparison(myMatrix, sequences, myScoreD2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 4322.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 3652.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 3652.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 4030.0, 0.01);

    kmerSize=5;
    myScoreD2.kmerSize=kmerSize;
    alignmentFreeComparison(myMatrix, sequences, myScoreD2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 762.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 216.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 216.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 688.0, 0.01);
}

SEQAN_DEFINE_TEST(test_alignment_free_d2star_dna5)
{
    using namespace seqan;
    StringSet<Dna5String> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    // kmerSize = 3
    unsigned kmerSize = 3;
    unsigned bgModelOrder = 0;
    AFScore<D2Star> myScoreD2Star(kmerSize,bgModelOrder);
    alignmentFreeComparison(myMatrix, sequences, myScoreD2Star);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 58.2374, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -10.949, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -10.949, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 48.358, 0.01);

    bgModelOrder = 1;
    myScoreD2Star.bgModelOrder=bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences,myScoreD2Star);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 30.336, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -15.104, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -15.104, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 36.3913, 0.01);
}

SEQAN_DEFINE_TEST(test_alignment_free_d2z_dna5)
{
    using namespace seqan;
    StringSet<Dna5String> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    //kmerSize = 3
    unsigned kmerSize = 3;
    unsigned bgModelOrder = 0;
    AFScore<D2z> myScoreD2z(kmerSize,bgModelOrder);
    alignmentFreeComparison(myMatrix, sequences, myScoreD2z);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 5.13022, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.614828, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.614828, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 3.40064, 0.01);

    bgModelOrder=1;
    myScoreD2z.bgModelOrder=bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences,myScoreD2z);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.61939, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 0.218295, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 0.218295, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 2.47939, 0.01);
}
SEQAN_DEFINE_TEST(test_alignment_free_N2_dna5)
{
    using namespace seqan;
    StringSet<Dna5String> sequences;
    alfTestHelperGetSequences(sequences);

    typedef Matrix<double, 2> TMatrix;
    TMatrix myMatrix;

    unsigned kmerSize = 3;
    unsigned bgModelOrder = 0;
    String<char>  revCom = "";
    unsigned mismatches = 0;
    double mismatchWeight = 0.5;
    AFScore<N2> myScoreN2(kmerSize, bgModelOrder, revCom, mismatches, mismatchWeight);

    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.161242, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.161242, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    bgModelOrder = 1;
    myScoreN2.bgModelOrder = bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), 0.143021, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), 0.143021, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    myScoreN2.revCom = "bothStrands";
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.236594, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.236594, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    mismatches = 1;
    mismatchWeight = 0.5;
    myScoreN2.mismatches = mismatches;
    myScoreN2.mismatchWeight = mismatchWeight;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.382932, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.382932, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    bgModelOrder=0;
    myScoreN2.bgModelOrder = bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.675514, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.675514, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    kmerSize = 4;
    myScoreN2.kmerSize = kmerSize;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.479099, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.479099, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

    bgModelOrder = 2;
    myScoreN2.bgModelOrder = bgModelOrder;
    alignmentFreeComparison(myMatrix, sequences, myScoreN2);

    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 0), 1.0, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 0, 1), -0.0858019, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 0), -0.0858019, 0.01);
    SEQAN_ASSERT_IN_DELTA(value(myMatrix, 1, 1), 1.0, 0.01);

}

SEQAN_BEGIN_TESTSUITE(test_alignment_free)
{
    // Call tests.
    SEQAN_CALL_TEST(test_alignment_free_d2_dna);
    SEQAN_CALL_TEST(test_alignment_free_d2_dna5);
    SEQAN_CALL_TEST(test_alignment_free_d2star_dna5);
    SEQAN_CALL_TEST(test_alignment_free_d2z_dna5);
    SEQAN_CALL_TEST(test_alignment_free_N2_dna5);
    SEQAN_CALL_TEST(test_alignment_free_calculateProbability);
    SEQAN_CALL_TEST(test_alignment_free_calculateVariance);
    SEQAN_CALL_TEST(test_alignment_free_calculateCovariance);
}
SEQAN_END_TESTSUITE
