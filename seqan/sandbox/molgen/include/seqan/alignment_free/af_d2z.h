// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Jonathan Goeke <goeke@molgen.mpg.de>
// ==========================================================================
// This header contains the implementation of the D2z score for alignment free
// sequence comparison (see Kantorovitz et al. Bioinformatics 2007, Volume23, 
// Issue13, Pp. i249-i255).
// These functions can be called with alignmentFreeComparison().
// ==========================================================================
#ifndef SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2Z_H_
#define SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2Z_H_

namespace seqan {

/*
 * D2z
 */
template <typename TStringSet, typename TValue>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, AFScore<D2z> const & score)
{
    typedef typename Value<TStringSet>::Type TString;
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;
    typedef Matrix<TValue, 2> TMatrix;
    typedef typename Iterator<TStringSet const>::Type       TIteratorSet;
    typedef typename Iterator<StringSet<String<unsigned int> > >::Type  TIteratorSetUnsigned;
    typedef typename Iterator<StringSet<String<double> > >::Type        TIteratorSetDouble;

    typedef typename Iterator<String<MarkovModel<TUnmaskedAlphabet> > >::Type      TIteratorMarkovModel;

    unsigned int seqNumber = length(sequenceSet);

    // resize the scoreMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    // fill(scoreMatrix,(TValue) 0);
    resize(scoreMatrix, (TValue) 0);

    StringSet<String<unsigned int> > kmerCounts;
    resize(kmerCounts, seqNumber);
    if (score.bgModelOrder == 0)
    {
        StringSet<String<double> > backgroundFrequencies;
        resize(backgroundFrequencies, seqNumber);

        // Count all kmers and all background nucleotide frequencies and store them in StringSets
        TIteratorSetUnsigned itKmerCounts = begin(kmerCounts);
        TIteratorSetDouble itBackgroundFrequencies = begin(backgroundFrequencies);

        TIteratorSet itSeqSet = begin(sequenceSet);
        for (; itSeqSet < end(sequenceSet); ++itSeqSet)
        {
            countKmers(value(itKmerCounts), value(itBackgroundFrequencies), value(itSeqSet), score.kmerSize);
            ++itKmerCounts;
            ++itBackgroundFrequencies;
        }
        if(score.verbose)
	{
		std::cout << "\ncounted words";
	}
	// calculate all pairwise scores and store them in scoreMatrix
        for (unsigned int rowIndex = 0; rowIndex < (seqNumber); ++rowIndex) //(remove diagonal seqNumber-1)
        {
		if(score.verbose)
		{
			std::cout << "\nSequence number " << rowIndex;
		}
            for (unsigned int colIndex = rowIndex; colIndex < (seqNumber); ++colIndex) //(remove diagonal rowIndex+1)
            {
                alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), kmerCounts[rowIndex], backgroundFrequencies[rowIndex], kmerCounts[colIndex], backgroundFrequencies[colIndex], score);
                value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex); // Copy symmetric entries
            }
        }
    }
    else
    {
        String<MarkovModel<TUnmaskedAlphabet> > backgroundModels;
        resize(backgroundModels, seqNumber, MarkovModel<TUnmaskedAlphabet>(score.bgModelOrder));
        // Count all kmers and all background nucleotide frequencies and store them in StringSets

        TIteratorMarkovModel itMM = begin(backgroundModels);
        TIteratorSet itSeqSet = begin(sequenceSet);

        // Count all kmers and all background nucleotide frequencies and store them in StringSets
        TIteratorSetUnsigned itKmerCounts = begin(kmerCounts);

        for (; itSeqSet < end(sequenceSet); ++itSeqSet)
        {
            countKmers(value(itKmerCounts), value(itMM), value(itSeqSet), score.kmerSize);

            // countKmers(value(itSeqSet),value(itKmerCounts),value(itMM),k);
            ++itKmerCounts;
            if (itMM < end(backgroundModels))
            {
                value(itMM)._computeAuxiliaryMatrices();
                ++itMM;
            }

        }
	if(score.verbose)
	{
		std::cout << "\ncounted words";
	}
        // calculate all pairwise scores and store them in scoreMatrix
        for (unsigned int rowIndex = 0; rowIndex < (seqNumber); ++rowIndex) //(seqNumber-1)
        {
		if(score.verbose)
		{
			std::cout << "\nSequence number " << rowIndex;
		}
            for (unsigned int colIndex = rowIndex; colIndex < (seqNumber); ++colIndex)
            {
                alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), kmerCounts[rowIndex], backgroundModels[rowIndex], kmerCounts[colIndex], backgroundModels[colIndex], score);
                value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex); // Copy symmetric entries
            }
        }
    }
}

template <typename TValue, typename TStringBG>
void alignmentFreeCompareCounts(TValue & result, String<unsigned int> & kmerCounts1, TStringBG & backgroundFrequencies1, String<unsigned int> & kmerCounts2, TStringBG & backgroundFrequencies2, AFScore<D2z> const & score)
{
    typedef typename Value<TStringBG>::Type TValueBG;
    // sum : inner product /d2 distance
    TValueBG sum = 0;
    unsigned int len1 = score.kmerSize - 1;
    unsigned int len2 = score.kmerSize - 1;
    unsigned int nvals = length(kmerCounts1);

    for (unsigned int l = 0; l < nvals; l++)
    {
        len1 += kmerCounts1[l];
        len2 += kmerCounts2[l];

        sum += kmerCounts1[l] * kmerCounts2[l];
    }

    TValueBG q1[4];
    TValueBG q2[4];
    for (int l = 0; l < 4; l++)
    {
        q1[l] = backgroundFrequencies1[l];
        q2[l] = backgroundFrequencies2[l];
        // std::cout<<backgroundFrequencies1[l]<<"\n";
    }

    // Compute expected value and variance (IID)
    TValueBG E = computeExpectationD2(len1, len2, score.kmerSize, q1, q2);
    TValueBG var = computeVarianceD2(len1, len2, score.kmerSize, q1, q2);
    // std::cout<<"\nlength1:"<<len1<<"length2: "<<len2<<"Expected value: "<<E<<"Variance: "<<var<<"Sum: "<<sum<<"\n";

    if ((var <= 0))
    {
	if(score.verbose)
	{
		std::cout << "\nlength1:" << len1 << "length2: " << len2 << "Expected value: " << E << "Variance: " << var << "Sum: " << sum << "\n";
		std::cout << "Error: negative variance\n";
	}
        result = 0;
        return;
    }

    // Calculate z-score
    result = (TValue) (sum - E) / pow(var, 0.5);
}

template <typename TAlphabet, typename TValue, typename TSpec>
void alignmentFreeCompareCounts(TValue & result, String<unsigned int> & kmerCounts1, MarkovModel<TAlphabet, TValue, TSpec> & bgModel1, String<unsigned int> & kmerCounts2, MarkovModel<TAlphabet, TValue, TSpec> & bgModel2, AFScore<D2z> const & score)
{
    unsigned int nvals = length(kmerCounts1);

    // Compute D2
    int sum = 0;
    int sumCounts1 = score.kmerSize - 1;
    int sumCounts2 = score.kmerSize - 1;

    for (unsigned int l = 0; l < nvals; l++)
    {
        sumCounts1 += kmerCounts1[l];
        sumCounts2 += kmerCounts2[l];
        sum += value(kmerCounts1, l) * value(kmerCounts2, l);       // Inner product
    }
    // std::cout<<"\nSum: "<<sum<<"\n";

    TValue D2 = (TValue) sum;

    // compute mean and variance
    TValue indicatorexpecation = 0;

    TValue E = computeExpectationD2(sumCounts1, sumCounts2, score.kmerSize, bgModel1, bgModel2, indicatorexpecation);
    // std::cout<<"\nE: "<<E<<"\n"<<"indicatorexpectation: "<<indicatorexpecation<<"\n";
    TValue var = computeVarianceD2(sumCounts1, sumCounts2, score.kmerSize, bgModel1, bgModel2, indicatorexpecation);
    // std::cout<<"\nvar: "<<var<<"\n"<<indicatorexpecation<<"\n";

    if (var <= 0)
    {
        fprintf(stderr, "Error: negative variance\n");
        result = 0;
        return;
    }
    // return z-score of D2
    result = (D2 - E) / pow(var, 0.5);

}

template <typename TValue>
double computeExpectationD2(int len1, int len2, unsigned int k, TValue * q1, TValue * q2)
{
    TValue p2 = 0;
    for (int i = 0; i < 4; i++)
        p2 += q1[i] * q2[i];

    int nbar1 = len1 - k + 1;
    int nbar2 = len2 - k + 1;

    TValue retval = nbar1;
    retval *= nbar2;
    retval *= pow(p2, k);
    return retval;
}

template <typename TAlphabet, typename TValue, typename TSpec>
double computeExpectationD2(int slen1, int slen2, unsigned int k, MarkovModel<TAlphabet, TValue, TSpec> & bkg1, MarkovModel<TAlphabet, TValue, TSpec> & bkg2, TValue & indicatorexpectation)
{
    unsigned int mo = bkg1.order;

    if (mo >= k)
    {
        fprintf(stderr, "Error: Can't suppport markov order greater or equal to word length\n");
        exit(1);
    }

    long emo  = 1 << (2 * mo); // pow(4, mo);
    // std::cout<<"\nemo: "<<emo<<"\n";

    long ekmo = 1 << (2 * (k - mo)); // pow(4, k-mo);

    TValue mean = 0;
    for (long i = 0; i < emo; i++)
    {
        // i is w1
        // p^A(w1)*p^B(w1)
        // double term = bkg1->p[i]*bkg2->p[i];		// p is stationary distribution
        TValue term = value(bkg1.stationaryDistribution, i) * value(bkg2.stationaryDistribution, i);
        // std::cout<<"\nterm: "<<term<<"\n";
        TValue subterm = 0;

        for (long j = 0; j < ekmo; j++)
        {
            // p^A_*(w)p^B_*(w)
            TValue subprob1 = ComputeWordProbGivenPrefix(i, j, bkg1, k, mo);
            TValue subprob2 = ComputeWordProbGivenPrefix(i, j, bkg2, k, mo);
            subterm += subprob1 * subprob2;
        }
        mean += term * subterm;
    }

    indicatorexpectation = mean; // this is E[Y_ij]
    int nbar1 = slen1 - k + 1;      // Number of kmers in sequence1
    int nbar2 = slen2 - k + 1;      // Number of kmers in sequence2
    return mean * ((TValue) nbar1 * nbar2);
}

template <typename TValue>
double computeVarianceD2(int len1, int len2, unsigned int k, TValue * q1, TValue * q2)
{
    int nbar1 = len1 - k + 1; // Number of kmers in sequence 1
    int nbar2 = len2 - k + 1; // Number of kmers in sequence 2

    int qbar1 = len1 - 2 * k + 2; // Number of overlapping kmers in sequence 1
    int qbar2 = len2 - 2 * k + 2;

    TValue p2 = 0, p31 = 0, p32 = 0;
    for (int i = 0; i < 4; i++)
    {
        p2 += q1[i] * q2[i]; // match seq1 und seq2 g_1,1
        p31 += q1[i] * q2[i] * q1[i]; // match seq 1 und seq2 und overlap in seq 1 g_2,1
        p32 += q1[i] * q2[i] * q2[i]; // match seq 1 und seq2 und overlap in seq 2 g_1,2
    }

    TValue variance = 0;
    // crabgrass with l=0 (= complete overlap):
    TValue power1 = pow(p32, k) - pow(p2, 2 * k); // overlap for seq2, l=0
    power1 *= TValue(nbar1) * TValue(qbar2) * TValue(qbar2 - 1); // sum over all possible (i,j) and (s,t)
    TValue power2 = pow(p31, k) - pow(p2, 2 * k); // overlap for seq1, l=0
    power2 *= TValue(nbar2) * TValue(qbar1) * TValue(qbar1 - 1);
    variance += power1 + power2;

    // crabgrasses with l>0:
    for (unsigned int l = 1; l <= k - 1; l++)
        variance += 2 * TValue(nbar1 - l) * TValue(qbar2) * TValue(qbar2 - 1) * (pow(p2, 2 * l) * pow(p32, k - l) - pow(p2, 2 * k)) + 2 * TValue(nbar2 - l) * TValue(qbar1) * TValue(qbar1 - 1) * (pow(p2, 2 * l) * pow(p31, k - l) - pow(p2, 2 * k));

    // accordion main diagonal:
    variance += TValue(nbar1) * TValue(nbar2) * (pow(p2, k) - pow(p2, 2 * k)); // l=0 term
    for (unsigned int l = 1; l <= k - 1; l++)
        variance += 2 * TValue(nbar1 - l) * TValue(nbar2 - l) * (pow(p2, k + l) - pow(p2, 2 * k)); // l > 0 terms

    return variance;
}

template <typename TAlphabet, typename TValue, typename TSpec>
double computeVarianceD2(int slen1, int slen2, unsigned int k, MarkovModel<TAlphabet, TValue, TSpec> & bkg1, MarkovModel<TAlphabet, TValue, TSpec> & bkg2, TValue indicatorexpectation)
{
    unsigned int mo = bkg1.order;
    if (mo >= k)
    {
        fprintf(stderr, "Error: Can't support markov order greater or equal to word length\n");
        exit(1);
    }


    long emo  = 1 << (2 * mo); // pow(4, mo);
    long ekmo = 1 << (2 * (k - mo)); // pow(4, k-mo);

    int q1 = slen1 - 2 * k + 2;
    int q2 = slen2 - 2 * k + 2;
    int nbar1 = slen1 - k + 1;
    int nbar2 = slen2 - k + 1;

    // create space for Sbctilde and fill in the values.
    Matrix<TValue, 2> Sbctilde1;
    setLength(Sbctilde1, 0, emo);
    setLength(Sbctilde1, 1, emo);
    resize(Sbctilde1);

    Matrix<TValue, 2> Sbctilde2;
    setLength(Sbctilde2, 0, emo);
    setLength(Sbctilde2, 1, emo);
    resize(Sbctilde2);

    for (long i = 0; i < emo; i++)
    {
        for (long j = 0; j < emo; j++)
        {
            value(Sbctilde1, i, j) = (q1 * (q1 + 1) / 2) * value(bkg1.stationaryDistribution, j) - (q1 - 1) * value(bkg1._qppp, i, j) - value(bkg1._qppqpp, i, j);
            value(Sbctilde2, i, j) = (q2 * (q2 + 1) / 2) * value(bkg2.stationaryDistribution, j) - (q2 - 1) * value(bkg2._qppp, i, j) - value(bkg2._qppqpp, i, j);
        }
    }

    // compute sums of word probs and star probs of x-mers (x between mo+1 and k) conditional on last or first (resp) mo-word
    Matrix<TValue, 2> sump; // sump[x][wk] is total word prob of every x-word ending with wk (which is of length mo)
    setLength(sump, 0, k + 1);
    setLength(sump, 1, emo);
    // resize(sump);
    resize(sump, 0.0);

    Matrix<TValue, 2> sumpstar; // sump[x][wk] is total word prob of every x-word ending with wk (which is of length mo)
    setLength(sumpstar, 0, k + 1);
    setLength(sumpstar, 1, emo);
    resize(sumpstar, 0.0);


    for (unsigned int x = 0; x <= k; x++)
    {
        if (x <= mo)
        {

            continue;
        }
        // sump[x] = new double[emo]; // sump[x][wk] is total word prob of every x-word ending with wk (which is of length mo)
        // sumpstar[x] = new double[emo]; // sumpstar[u1] is star word prob of every x-word beginning with u1 (of length mo)
        for (long wk = 0; wk < emo; wk++)
        {
            TValue sum = 0;
            long exmo = 1 << (2 * (x - mo)); // pow(4,x-mo);
            for (long wpre = 0; wpre < exmo; wpre++)
            {
                // w = wpre wk
                sum += ComputeWordProb((wpre << (2 * mo)) + wk, bkg1, x, mo) * ComputeWordProb((wpre << (2 * mo)) + wk, bkg2, x, mo);
            }
            value(sump, x, wk) = sum;
        }

        for (long u1 = 0; u1 < emo; u1++)
        {
            TValue sum = 0;
            long exmo = 1 << (2 * (x - mo)); // pow(4,x-mo);
            for (long usuf = 0; usuf < exmo; usuf++)
            {
                // u = u1 usuf
                sum += ComputeWordProbGivenPrefix(u1, usuf, bkg1, x, mo) * ComputeWordProbGivenPrefix(u1, usuf, bkg2, x, mo);
            }
            value(sumpstar, x, u1) = sum;
        }
    }

    // handle the non overlap terms first
    TValue covnonoverlap = 0;
    for (long wk = 0; wk < emo; wk++)
    {
        for (long u1 = 0; u1 < emo; u1++)
        {
            TValue term = value(Sbctilde1, wk, u1) * value(Sbctilde2, wk, u1) * value(sump, k, wk) * value(sumpstar, k, u1);
            covnonoverlap += 4 * term;
        }
    }
    // std::cout<<"\nnonoverlap: "<<covnonoverlap;
    TValue subtractFromNonOverlap = (TValue)q1 * (q1 - 1) * q2 * (q2 - 1);
    subtractFromNonOverlap *= pow(indicatorexpectation, 2);
    covnonoverlap -= subtractFromNonOverlap;
    // std::cout<<" subtract:"<<subtractFromNonOverlap<<"\n";
    // compute crabgrass terms
    TValue covcrabgrass = 0;

    // case 1: overlap >= mo
    for (unsigned int m = 1; m <= k - mo; m++)
    {
        long ekm = 1 << (2 * (k - m)); // pow(4,k-m);
        long moonesflag = (1 << (2 * mo)) - 1;
        long kmmoonesflag = (1 << (2 * (k - m - mo))) - 1;
        for (long v = 0; v < ekm; v++) // all strings of length k-m
        {
            long vsuf = v & moonesflag; // last morder chars of v
            long vpre = v >> (2 * (k - m - mo)); // first morder chars of v, equivalent to v/pow(4,k-m-mo)
            long vsuf2 = v & kmmoonesflag; // remaining k-m-morder chars of v; equivalent to v%pow(4,k-m-mo)
            // note: in the following, if m=k-mo, ComputeWordProbGivenPrefix(.,.,.,k-m,mo) will return 1, as it should
            // overlap in A, separate in B
            TValue term1 = value(Sbctilde2, vsuf, vpre);
            term1 *= ComputeWordProbGivenPrefix(vpre, vsuf2, bkg1, k - m, mo);
            term1 *= pow(ComputeWordProbGivenPrefix(vpre, vsuf2, bkg2, k - m, mo), 2);
            term1 *= value(sump, m + mo, vpre);
            term1 *= value(sumpstar, m + mo, vsuf);
            // overlap in B, separate in A
            TValue term2 = value(Sbctilde1, vsuf, vpre);
            term2 *= ComputeWordProbGivenPrefix(vpre, vsuf2, bkg2, k - m, mo);
            term2 *= pow(ComputeWordProbGivenPrefix(vpre, vsuf2, bkg1, k - m, mo), 2);
            term2 *= value(sump, m + mo, vpre);
            term2 *= value(sumpstar, m + mo, vsuf);

            covcrabgrass += 4 * q1 * term1 + 4 * q2 * term2;
        }
    }

    // case 2: overlap < morder
    for (unsigned int m = k - mo + 1; m <= k - 1; m++)
    {
        int tlen = 2 * mo - (k - m);
        long et = 1 << (2 * tlen); // pow(4,tlen);
        long moonesflag = (1 << (2 * mo)) - 1;
        long tmoonesflag = (1 << (2 * (tlen - mo))) - 1;
        for (long t = 0; t < et; t++) // all strings of length tlen
        {
            long tsuf = t & moonesflag; // last morder chars of t; equivalent to t%emo;
            long tpre = t >> (2 * (tlen - mo)); // first morder chars of t, equivalent to t/pow(4,tlen-mo)
            long tsuf2 = t & tmoonesflag; // remaining tlen-morder chars of t, equivalent to t%etmo;
            // overlap in A, separate in B
            TValue term1 = value(Sbctilde2, tpre, tsuf);
            term1 *= ComputeWordProbGivenPrefix(tpre, tsuf2, bkg1, tlen, mo);
            term1 *= value(sump, k, tpre);
            term1 *= value(sumpstar, k, tsuf);
            // overlap in B, separate in A
            TValue term2 = value(Sbctilde1, tpre, tsuf);
            term2 *= ComputeWordProbGivenPrefix(tpre, tsuf2, bkg2, tlen, mo);
            term2 *= value(sump, k, tpre);
            term2 *= value(sumpstar, k, tsuf);

            covcrabgrass += 4 * q1 * term1 + 4 * q2 * term2;
        }
    }

    // case 3: m=0 (complete overlap)
    // long ek = 1<<(2*k); // pow(4,k); unused
    long moonesflag = (1 << (2 * mo)) - 1;
    for (long wpre = 0; wpre < emo; wpre++)
    {
        // first morder chars of w
        TValue term1 = 0;
        TValue term2 = 0;
        for (long wsuf2 = 0; wsuf2 < ekmo; wsuf2++)
        {
            // remaining k-morder chars of w
            long wsuf = ((wpre << (2 * (k - mo))) + wsuf2) & moonesflag; // last morder chars of w
            term1 += 2 * nbar1 * value(Sbctilde2, wsuf, wpre) * ComputeWordProbGivenPrefix(wpre, wsuf2, bkg1, k, mo) * pow(ComputeWordProbGivenPrefix(wpre, wsuf2, bkg2, k, mo), 2);
            term2 += 2 * nbar2 * value(Sbctilde1, wsuf, wpre) * ComputeWordProbGivenPrefix(wpre, wsuf2, bkg2, k, mo) * pow(ComputeWordProbGivenPrefix(wpre, wsuf2, bkg1, k, mo), 2);
        }
        covcrabgrass += (term1 + term2) * value(bkg1.stationaryDistribution, wpre) * value(bkg2.stationaryDistribution, wpre);
    }

    // IGNORED EDGE EFFECTS

    // finally, subtract expectations
    TValue subtractfromcrabgrass = (((2 * q1 * k - nbar1) * (double)q2 * double(q2 - 1)) + ((2 * q2 * k - nbar2) * (double)q1 * double(q1 - 1)));
    subtractfromcrabgrass *= pow(indicatorexpectation, 2);
    covcrabgrass -= subtractfromcrabgrass;

    // compute accordion main diagonal terms
    TValue covaccordiondiag = 0;

    for (unsigned int m = 1; m <= k - 1; m++)
    {
        long tlen;
        if (m <= mo - 1)
            tlen = m + mo;
        else
            tlen = 2 * mo - 1;
        long et = 1 << (2 * tlen); // pow(4,tlen);
        long moonesflag = (1 << (2 * mo)) - 1;
        long tmoonesflag = (1 << (2 * (tlen - mo))) - 1;

        TValue mterm = 0;
        for (long t = 0; t < et; t++) // all strings of length tlen
        {
            TValue term = 1;
            long tpre = t >> (2 * (tlen - mo)); // first morder chars of t, equivalent to t/pow(4,tlen-mo)
            long tsuf2 = t & tmoonesflag; // remaining tlen-morder chars of t, equivalent to t%etmo;
            term *= ComputeWordProbGivenPrefix(tpre, tsuf2, bkg1, tlen, mo) * ComputeWordProbGivenPrefix(tpre, tsuf2, bkg2, tlen, mo);
            term *= value(sump, k, tpre);
            if (m >= mo)
            {
                long tsuf = t & moonesflag; // last morder chars of t; equivalent to t%emo;
                term *= value(sumpstar, m + 1, tsuf);
            }
            mterm += term;
        }
        covaccordiondiag += 2 * (nbar1 - m) * (nbar2 - m) * mterm;
    }
    covaccordiondiag += nbar1 * nbar2 * indicatorexpectation; // complete overlap term
    covaccordiondiag -= (nbar1 * nbar2 * (2 * k - 1) - (nbar1 + nbar2) * k * (k - 1) + (k - 1) * k * (2 * k - 1) / 3) * pow(indicatorexpectation, 2);

    // IGNORED NON-DIAGONAL ACCORDION TERMS
/*
    // clean up
    for (int x=0; x<=k; x++)
    {
        if (x <= mo) continue;
        delete [] sump[x];
        delete [] sumpstar[x];
    }
    delete [] sump;
    delete [] sumpstar;
    for (long i=0; i<emo; i++) {
    delete [] Sbctilde1[i];
    delete [] Sbctilde2[i];
  }
  delete [] Sbctilde1;
  delete [] Sbctilde2;
*/
// return values
// fprintf(stderr,"%f\t%f\t%f\n",covnonoverlap,covcrabgrass,covaccordiondiag);
    TValue variance = covnonoverlap + covcrabgrass + covaccordiondiag;
    return variance;

    // return 0;
}

template <typename TAlphabet, typename TValue, typename TSpec>
double ComputeWordProb(long word, MarkovModel<TAlphabet, TValue, TSpec> & bkg, unsigned int k, int mo)
{
// need to know mo, since that decides the prefix
// need to know total length of word
    long prefix = word >> (2 * (k - mo)); // equivalent to word/pow(4,k-mo); gets first mo chars of word
    long suffix = word & ((1 << (2 * (k - mo))) - 1); // bit-and with a number that is all ones in the last 2*(k-mo) bits, gets last (k-mo) chars of word
    return value(bkg.stationaryDistribution, prefix) * ComputeWordProbGivenPrefix(prefix, suffix, bkg, k, mo);
    // return bkg->p[prefix]*ComputeWordProbGivenPrefix(prefix,suffix,bkg,k,mo);
}

template <typename TAlphabet, typename TValue, typename TSpec>
double ComputeWordProbGivenPrefix(long prefix, long suffix, MarkovModel<TAlphabet, TValue, TSpec> & bkg, unsigned int k, unsigned int mo)
{
// computes p_*(prefix suffix)
    // calculate only if k-mo>=1; otherwise return 1
    if (k - mo <= 0)
        return 1;

    TValue prob = 1;
    long pre = prefix; // assumed of length mo
    long jcopy = suffix; // assumed of length k-mo

    // need to iterate through successive positions of suffix
    for (unsigned int l = 0; l < k - mo; l++)
    {
        long jcopyfirst = jcopy >> (2 * (k - mo - 1)); // get first char of jcopy, which is of length k-mo; equivalent to jcopy/pow(4,k-mo-1)
        long suf = ((pre & ((1 << (2 * (mo - 1))) - 1)) << 2) + jcopyfirst; // last (mo-1) chars of pre, cat jcopyfirst; equivalent to (pre%pow(4,mo-1))*4+jcopyfirst;

        // prob *= bkg->P[pre][suf];	// P is transition matrix
        prob *= value(bkg.transition, pre, suf);
        // erase first char of jcopy, and preserve its length at k-mo by shifting one char to left; equivalent to (jcopy%pow(4,k-mo-1))*4;
        jcopy = (jcopy & ((1 << (2 * (k - mo - 1))) - 1)) << 2;
        pre = suf;
    }
    return prob;
}

}  // namespace seqan

#endif  // SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2Z_H_
