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
// This header contains the implementation of the D2star score for alignment 
// free sequence comparison.
// (see Reinert et al. J Comput Biol. 2009 Dec;16(12):1615-34.)
// These functions can be called with alignmentFreeComparison().
// ==========================================================================
#ifndef SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_ORIGINAL_H_
#define SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_ORIGINAL_H_

namespace seqan {

template <typename TValue, typename TStringSet>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, AFScore<D2Star> const & score)
{



    typedef typename Value<TStringSet>::Type TString;
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;
    typedef Matrix<TValue, 2> TMatrix;

    typedef typename Iterator<TStringSet const>::Type       TIteratorSet;
    typedef typename Iterator<StringSet<String<unsigned int> > >::Type      TIteratorSetInt;
    typedef typename Iterator<StringSet<String<double> > >::Type        TIteratorSetDouble;

    unsigned int seqNumber = length(sequenceSet);

    // resize the distMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    resize(scoreMatrix, (TValue) 0);

    StringSet<String<double> > standardisedKmerCounts;
    resize(standardisedKmerCounts, seqNumber);
    // Count all kmers and all background nucleotide frequencies and store them in StringSets
    TIteratorSetDouble itStandardisedKmerCounts = begin(standardisedKmerCounts);
    TIteratorSet itSeqSet = begin(sequenceSet);

    // calculate all pairwise scores and store them in scoreMatrix
    for (unsigned int rowIndex = 0; rowIndex < (seqNumber); ++rowIndex) //(remove diagonal values seqNumber-1)
    {
        if(score.verbose)
	{
	  std::cout << "\nSequence number " << rowIndex;
	}
	for (unsigned int colIndex = rowIndex; colIndex < (seqNumber); ++colIndex) // remove diagonal values rowIndex+1
        {
            d2star(value(scoreMatrix, rowIndex, colIndex), sequenceSet[rowIndex], sequenceSet[colIndex], score);
// std::cout<<"\noriginal: "<<value(scoreMatrix,rowIndex,colIndex)<<"\n";

            // alignmentFreeCompareCounts(value(scoreMatrix,rowIndex,colIndex), standardisedKmerCounts[rowIndex],standardisedKmerCounts[colIndex],score);
            value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  // Copy symmetric entries
        }
    }
}

/*
 * Calculate pairwise score given the counts of all kmers
 */
template <typename TValue, typename TString>
void
alignmentFreeCompareCounts(TValue & result, TString & kmerCounts1, TString & kmerCounts2, AFScore<D2Star> const & score)
{
    typedef typename Value<TString>::Type TStringValue;
    typedef typename Iterator<TString, Rooted>::Type        TIteratorTString;

    TIteratorTString it1 = begin(kmerCounts1);
    TIteratorTString it2 = begin(kmerCounts2);
    result = 0.0;
    // TValue resultRC=0.0;

    for (; it1 < end(kmerCounts1); ++it1)
    {
        // Multiply standardised counts
        result += (TValue)(value(it1) * value(it2));

        ++it2;
    }
}

template <typename TValue, typename TSequence>
void d2star(TValue & result, TSequence const & sequence1, TSequence const & sequence2, AFScore<D2Star> const & score)
{


    typedef typename Value<TSequence>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;
    // typedef typename Value<TString>::Type TValue;
    typedef typename Iterator<String<unsigned int>, Rooted>::Type       TIteratorUnsigned;
    // typedef typename Iterator<TString ,Rooted>::Type		TIteratorTString;

    TValue missing = -pow(10, 10);
    TSequence seq1seq2;
    append(seq1seq2, sequence1);
    append(seq1seq2, sequence2);
    result = 0.0;
    /*
    *--------ORDER 0 BACKGROUND MODEL----------
    */
    if (score.bgModelOrder == 0)
    {
        // String<unsigned int> kmerCounts;
        String<unsigned int> kmerCounts1;
        String<unsigned int> kmerCounts2;
        String<unsigned int> backgroundCounts;
        String<double> backgroundFrequencies;
        resize(backgroundFrequencies, 4, 0);
        countKmers(kmerCounts1, sequence1, score.kmerSize);
        countKmers(kmerCounts2, sequence2, score.kmerSize);
        countKmers(backgroundCounts, seq1seq2, 1);
        int sumBG = 0;
        for (unsigned i = 0; i < length(backgroundCounts); ++i)
        {
            sumBG += backgroundCounts[i];
        }
        for (unsigned i = 0; i < length(backgroundCounts); ++i)
        {
            backgroundFrequencies[i] = backgroundCounts[i] / ((double)sumBG);
        }
// int tmpSum=0;
// for(int tmp=0;tmp<length(kmerCounts);++tmp){std::cout<<kmerCounts[tmp]<<"\n";tmpSum+=kmerCounts[tmp];}
// std::cout<<"\ntmpSum: "<<tmpSum<<"\n";
        unsigned nvals = length(kmerCounts1); // Number of kmers
        int len1 = 0;
        int len2 = 0;

        for (unsigned l = 0; l < nvals; l++)
        {
            len1 += kmerCounts1[l];
            len2 += kmerCounts2[l];

        }
        // fill(standardisedCounts,nvals,(TValue) 0.0);

        // string of tvalue to store p_w
        String<TValue> probabilities;
        resize(probabilities, nvals, missing);

        // TIteratorUnsigned itCounts;
        // TIteratorTString itStandardisedCounts;

        // itCounts=begin(kmerCounts);
        // itStandardisedCounts=begin(standardisedCounts);
        //TValue normValue1 = 0.0;
        //TValue normValue2 = 0.0;

        for (unsigned i = 0; i < nvals; ++i)
        {
            TValue p_w = 1;   // Probability of kmer

            String<TUnmaskedAlphabet> w;
            unhash(w, i, score.kmerSize);
            calculateProbability(p_w, w, backgroundFrequencies);
            TValue variance1 = 0.0;
            TValue variance2 = 0.0;

            variance1 = pow(len1 * p_w, 0.5);
            variance2 = pow(len2 * p_w, 0.5);
            // TValue varTMP= pow(((TValue) len1)*((TValue) len1),0.5)*p_w;
// std::cout<<"\n1: "<<variance1<<"\n2: "<<variance2<<"\n1*2: "<<(variance1*variance2)<<"\nold: "<<varTMP;

            // test if variance not 0 or inf before dividing
            if ((variance1 > missing) && (variance1 < pow(10, 10)))
            {
                if (p_w > 0)
                {
                    TValue stCount1 = (kmerCounts1[i] - p_w * len1) / variance1;
                    TValue stCount2 = (kmerCounts2[i] - p_w * len2) / variance2;
                    //normValue1 += pow(stCount1, 2);
                    //normValue2 += pow(stCount2, 2);

                    // value(itStandardisedCounts)=(((TValue) ((TValue) kmerCounts1[i])-p_w*((TValue)len1))*((TValue) ((TValue) kmerCounts2[i])-p_w*((TValue)len2)))/((TValue) variance);
                    // result+=(((TValue) ((TValue) kmerCounts1[i])-p_w*((TValue)len1))*((TValue) ((TValue) kmerCounts2[i])-p_w*((TValue)len2)))/((TValue) variance);
                    result += stCount1 * stCount2;
                    // result+=((kmerCounts1[i]-p_w*len1)*(kmerCounts2[i]-p_w*len2))/varTMP;
// TValue normValue1=(((TValue) ((TValue) kmerCounts1[i])-p_w*((TValue)len1))*((TValue) ((TValue) kmerCounts2[i])-p_w*((TValue)len2)))/((TValue) variance);
// std::cout<<"\n1: "<<normValue1;
// TValue normValue2=((kmerCounts1[i]-p_w*len1)*(kmerCounts2[i]-p_w*len2))/variance;
// std::cout<<"\n2: "<<(normValue1)<<"="<<(normValue2);
                }
            }
// std::cout<<"\n"<<w<<"; counts:"<<value(itCounts)<<"; stCounts:"<<value(itStandardisedCounts)<<" p_w: "<<p_w<<" var:"<<variance;
        }

    }
    /*
    *--------HIGHER ORDER BACKGROUND MODEL----------
    */
    else
    {
        String<unsigned int> kmerCounts1;
        String<unsigned int> kmerCounts2;
        StringSet<String<TUnmaskedAlphabet> > bgSequences;
        stringToStringSet(bgSequences, seq1seq2); // create unmasked sequnces
        MarkovModel<TUnmaskedAlphabet, TValue> backgroundModel(score.bgModelOrder);
        buildMarkovModel(backgroundModel, bgSequences);
        countKmers(kmerCounts1, sequence1, score.kmerSize);
        countKmers(kmerCounts2, sequence2, score.kmerSize);


        // countKmers(sequence,kmerCounts,backgroundModel,k);
        unsigned nvals = length(kmerCounts1); // Number of kmers
        int len1 = 0;
        int len2 = 0;

        for (unsigned l = 0; l < nvals; l++)
        {
            len1 += kmerCounts1[l];
            len2 += kmerCounts2[l];

        }
        // fill(standardisedCounts,nvals,(TValue) 0.0);

        String<TValue> probabilities;

        resize(probabilities, nvals, missing);
        // TIteratorUnsigned itCounts;
        // TIteratorTString itStandardisedCounts;

        // itCounts=begin(kmerCounts);
        // itStandardisedCounts=begin(standardisedCounts);
        // double sumTMP=0;
//      for(;itCounts<end(kmerCounts);++itCounts)
//      {
        //TValue normValue1 = 0.0;
        //TValue normValue2 = 0.0;

        for (unsigned i = 0; i < nvals; ++i)
        {
            TValue p_w = 1.0; // Probability of kmer
            TValue variance = 0.0;
            String<TUnmaskedAlphabet> w;
            unhash(w, i, score.kmerSize);
            p_w = emittedProbability(backgroundModel, w);
            variance = ((TValue) pow(((TValue) len1 * len2), 0.5)) * p_w;
            TValue variance1 = 0.0;
            TValue variance2 = 0.0;

            variance1 = pow(len1 * p_w, 0.5);
            variance2 = pow(len2 * p_w, 0.5);

            // Calculate standardised kmer Count
            // std::cout<<"\nword:"<<w<<", p_w:"<<p_w<<", var: "<<((TValue) pow(((TValue) len1)*p_w,0.5));
            if ((variance > pow(10, -10)) && (variance < pow(10, 10)))
            {
                if (p_w > 0)
                {

                    TValue stCount1 = (kmerCounts1[i] - p_w * len1) / variance1;
                    TValue stCount2 = (kmerCounts2[i] - p_w * len2) / variance2;

                    //normValue1 += pow(stCount1, 2);
                   // normValue2 += pow(stCount2, 2);

                    result += stCount1 * stCount2;
                    // result+=(((TValue) ((TValue) kmerCounts1[i])-p_w*((TValue)len1))*((TValue) ((TValue) kmerCounts2[i])-p_w*((TValue)len2)))/((TValue) variance);
                    // value(itStandardisedCounts)=((TValue) ((TValue) value(itCounts)) -p_w*((TValue)len1))/((TValue) pow(((TValue) len1)*p_w,0.5));
                }
                //++itStandardisedCounts;
            }
        }

       // if (score.norm == true)
         //   result = result / pow(normValue1 * normValue2, 0.5);
    }
}

/*
 * count kmers and standardise count vectors for Dna5 and markov model background
 */
// template <typename TString, typename TSequence>
// void standardiseCounts(TString & standardisedCounts, TSequence const & sequence, AFScore<D2Star> const & score)
// {
// 
// 
//     typedef typename Value<TSequence>::Type TAlphabet;
//     typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;
//     typedef typename Value<TString>::Type TValue;
//     typedef typename Iterator<String<unsigned int>, Rooted>::Type       TIteratorUnsigned;
//     typedef typename Iterator<TString, Rooted>::Type        TIteratorTString;
// 
//     TValue missing = -pow(10, 10);
// 
//     /*
//     *--------ORDER 0 BACKGROUND MODEL----------
//     */
//     if (score.bgModelOrder == 0)
//     {
// 
//         String<unsigned int> kmerCounts;
//         String<double> backgroundFrequencies;
//         countKmers(kmerCounts, backgroundFrequencies, sequence, score.kmerSize);
// // int tmpSum=0;
// // for(int tmp=0;tmp<length(kmerCounts);++tmp){std::cout<<kmerCounts[tmp]<<"\n";tmpSum+=kmerCounts[tmp];}
// // std::cout<<"\ntmpSum: "<<tmpSum<<"\n";
//         int nvals = length(kmerCounts); // Number of kmers
//         int len1 = 0;
//         for (int l = 0; l < nvals; l++)
//         {
//             len1 += kmerCounts[l];
//         }
//         resize(standardisedCounts, nvals, (TValue) 0.0);
// 
// 
//         // string of tvalue to store p_w
//         String<TValue> probabilities;
//         resize(probabilities, nvals, missing);
// 
//         TIteratorUnsigned itCounts;
//         TIteratorTString itStandardisedCounts;
// 
//         itCounts = begin(kmerCounts);
//         itStandardisedCounts = begin(standardisedCounts);
// 
//         for (; itCounts < end(kmerCounts); ++itCounts)
//         {
//             TValue p_w = 1;   // Probability of kmer
// 
//             String<TUnmaskedAlphabet> w;
//             unhash(w, (unsigned)position(itCounts), score.kmerSize);
//             calculateProbability(p_w, w, backgroundFrequencies);
//             TValue variance = 0;
// 
//             variance = ((TValue) len1) * p_w;
// 
// 
//             // test if variance not 0 or inf before dividing
//             if ((variance > pow(10, -10)) && (variance < pow(10, 10)))
//             {
//                 if (p_w > 0)
//                     value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts)) - p_w * ((TValue)len1)) / ((TValue) pow(variance, 0.5));
//             }
// // std::cout<<"\n"<<w<<"; counts:"<<value(itCounts)<<"; stCounts:"<<value(itStandardisedCounts)<<" p_w: "<<p_w<<" var:"<<variance;
//             ++itStandardisedCounts;
//         }
//     }
//     /*
//     *--------HIGHER ORDER BACKGROUND MODEL----------
//     */
//     else
//     {
//         String<unsigned int> kmerCounts;
//         MarkovModel<TUnmaskedAlphabet, TValue> backgroundModel(score.bgModelOrder);
//         countKmers(kmerCounts, backgroundModel, sequence, score.kmerSize);
// 
//         // countKmers(sequence,kmerCounts,backgroundModel,k);
//         int nvals = length(kmerCounts); // Number of kmers
//         int len1 = 0;
//         for (int l = 0; l < nvals; l++)
//             len1 += kmerCounts[l];
//         resize(standardisedCounts, nvals, (TValue) 0.0);
// 
//         String<TValue> probabilities;
// 
//         resize(probabilities, nvals, missing);
// 
//         TIteratorUnsigned itCounts;
//         TIteratorTString itStandardisedCounts;
// 
//         itCounts = begin(kmerCounts);
//         itStandardisedCounts = begin(standardisedCounts);
//         // double sumTMP=0;
//         for (; itCounts < end(kmerCounts); ++itCounts)
//         {
//             TValue p_w = 1;   // Probability of kmer
//             TValue variance = 0;
//             String<TUnmaskedAlphabet> w;
//             unhash(w, (unsigned)position(itCounts), score.kmerSize);
//             p_w = emittedProbability(backgroundModel, w);
//             variance = ((TValue) pow(((TValue) len1) * p_w, 0.5));
//             // Calculate standardised kmer Count
//             // std::cout<<"\nword:"<<w<<", p_w:"<<p_w<<", var: "<<((TValue) pow(((TValue) len1)*p_w,0.5));
//             if ((variance > pow(10, -10)) && (variance < pow(10, 10)))
//             {
//                 if (p_w > 0)
//                 {
//                     value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts)) - p_w * ((TValue)len1)) / ((TValue) pow(((TValue) len1) * p_w, 0.5));
//                 }
//                 ++itStandardisedCounts;
//             }
//         }
//     }
// }

}  // namespace seqan

#endif  // SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_ORIGINAL_H_
