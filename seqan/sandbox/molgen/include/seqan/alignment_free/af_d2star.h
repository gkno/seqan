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

#ifndef SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_H_
#define SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_H_

namespace seqan {

String<unsigned> revComIndex_star;
void initialiseRevComIndex_star(unsigned k)
{
    unsigned myLength = (unsigned) pow(4, k);
    resize(revComIndex_star, myLength, 0);
    Shape<Dna, SimpleShape> myShape;
    //Declare variables
    //TShape myShape;		//Shape, length can be changed (kmer_length)
    resize(myShape, k);
    for (unsigned i = 0; i < myLength; ++i)
    {
        String<Dna> w;
        unhash(w, i, k);
        DnaStringReverseComplement wRC(w);
        unsigned hashValue = hash(myShape, begin(wRC));
        revComIndex_star[i] = hashValue;
    }

}

template <typename TValue, typename TStringSet>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, AF_Score<D2star> const & score)
{



    typedef typename Value<TStringSet>::Type TString;
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;
    typedef Matrix<TValue, 2> TMatrix;

    typedef typename Iterator<TStringSet const>::Type       TIteratorSet;
    typedef typename Iterator<StringSet<String<unsigned int> > >::Type      TIteratorSetInt;
    typedef typename Iterator<StringSet<String<double> > >::Type        TIteratorSetDouble;


    //Initialise reverse complement hash table
    initialiseRevComIndex_star(score.kmerSize);
    unsigned int seqNumber = length(sequenceSet);

    //resize the distMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    resize(scoreMatrix, (TValue) 0);


    StringSet<String<double> > standardisedKmerCounts;
    resize(standardisedKmerCounts, seqNumber);
    //Count all kmers and all background nucleotide frequencies and store them in StringSets
    TIteratorSetDouble itStandardisedKmerCounts = begin(standardisedKmerCounts);
    TIteratorSet itSeqSet = begin(sequenceSet);
    for (; itSeqSet < end(sequenceSet); ++itSeqSet)
    {
        standardisedCounts(value(itStandardisedKmerCounts), value(itSeqSet), score);
        std::cout << "\n" << position(itSeqSet);
        ++itStandardisedKmerCounts;
    }
    //output of pairwise kmer weight to file
    std::ofstream myfile;
    if (score.outputFile != "")
    {
        myfile.open(toCString(score.outputFile));
        for (unsigned i = 0; i < length(standardisedKmerCounts[0]); ++i)
        {
            String<TUnmaskedAlphabet> w;
            unhash(w, i, score.kmerSize);
            //std::ofstream myfile;

            myfile << "\t" << w;
        }
        myfile << "\n";
        for (unsigned int seqIndex = 0; seqIndex < seqNumber; ++seqIndex)
        {
            myfile << "Seq" << seqIndex;
            for (unsigned i = 0; i < length(standardisedKmerCounts[seqIndex]); ++i)
            {
                myfile << "\t" << standardisedKmerCounts[seqIndex][i];
            }
            myfile << "\n";
        }
        myfile.close();
    }

    std::cout << "\ncounted words";


    //calculate all pairwise scores and store them in scoreMatrix
    for (unsigned int rowIndex = 0; rowIndex < (seqNumber); ++rowIndex)
    {
        std::cout << "\nSequence number " << rowIndex;
        for (unsigned int colIndex = rowIndex; colIndex < (seqNumber); ++colIndex)
        {
            alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), standardisedKmerCounts[rowIndex], standardisedKmerCounts[colIndex], score);
            value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  //Copy symmetric entries
        }
    }
}

/*
 * Calculate pairwise score given the counts of all kmers
 */
template <typename TValue, typename TString>
void
alignmentFreeCompareCounts(TValue & result, TString & kmerCounts1, TString & kmerCounts2, AF_Score<D2star> const & score)
{
    typedef typename Value<TString>::Type TStringValue;
    typedef typename Iterator<TString, Rooted>::Type        TIteratorTString;

    TIteratorTString it1 = begin(kmerCounts1);
    TIteratorTString it2 = begin(kmerCounts2);
    result = 0.0;
    TValue resultRC = 0.0;

    for (; it1 < end(kmerCounts1); ++it1)
    {


        //Multiply standardised counts
        result += (TValue)(value(it1) * value(it2));


        //Computation of the reverse complement strand score
        if ((score.revCom != ""))  //"min" "max" "mean" ""
        {
            unsigned hashValue = revComIndex_star[position(it1)];
            resultRC += (TValue)(value(it1) * kmerCounts2[hashValue]);
        }
        ++it2;
    }

    //std::cout<<"\nresult"<<result<<"\n"<<"\nresultRC"<<resultRC<<"\n";

    //Compute reverse complements scores if revCom!=""
    if (score.revCom == "mean")
    {
        result = (TValue) (resultRC + result) / 2;
    }
    else if (score.revCom == "max")
    {
        result = max(resultRC, result);
    }
    else if (score.revCom == "min")
    {
        result = min(resultRC, result);
    }


}

/*
 * count kmers and standardise count vectors for Dna5 and markov model background
 */
template <typename TString, typename TSequence>
void standardisedCounts(TString & standardisedCounts, TSequence const & sequence, AF_Score<D2star> const & score)
{


    typedef typename Value<TSequence>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;
    typedef typename Value<TString>::Type TValue;
    typedef typename Iterator<String<unsigned int>, Rooted>::Type       TIteratorUnsigned;
    typedef typename Iterator<TString, Rooted>::Type        TIteratorTString;

    TValue missing = -pow(10, 10);

    /*
    *--------ORDER 0 BACKGROUND MODEL----------
    */
    if (score.bgModelOrder == 0)
    {

        String<unsigned int> kmerCounts;
        String<double> backgroundFrequencies;
        countKmers(kmerCounts, backgroundFrequencies, sequence, score.kmerSize);
//int tmpSum=0;
//for(int tmp=0;tmp<length(kmerCounts);++tmp){std::cout<<kmerCounts[tmp]<<"\n";tmpSum+=kmerCounts[tmp];}
//std::cout<<"\ntmpSum: "<<tmpSum<<"\n";
        int nvals = length(kmerCounts); //Number of kmers
        int len1 = 0;
        for (int l = 0; l < nvals; l++)
        {
            len1 += kmerCounts[l];
        }
        resize(standardisedCounts, nvals, (TValue) 0.0);


        //string of tvalue to store p_w
        String<TValue> probabilities;
        resize(probabilities, nvals, missing);

        TIteratorUnsigned itCounts;
        TIteratorTString itStandardisedCounts;

        itCounts = begin(kmerCounts);
        itStandardisedCounts = begin(standardisedCounts);

        for (; itCounts < end(kmerCounts); ++itCounts)
        {
            TValue p_w = 1;   //Probability of kmer

            String<TUnmaskedAlphabet> w;
            unhash(w, (unsigned)position(itCounts), score.kmerSize);
            calculateProbability(p_w, w, backgroundFrequencies);
            TValue variance = 0;

            variance = ((TValue) len1) * p_w;


            //test if variance not 0 or inf before dividing
            if ((variance > pow(10, -10)) && (variance < pow(10, 10)))
            {
                if (p_w > 0)
                {
                    value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts)) - p_w * ((TValue)len1)) / ((TValue) pow(variance, 0.5));

                }
            }
//std::cout<<"\n"<<w<<"; counts:"<<value(itCounts)<<"; stCounts:"<<value(itStandardisedCounts)<<" p_w: "<<p_w<<" var:"<<variance;
            ++itStandardisedCounts;
        }
    }
    /*
    *--------HIGHER ORDER BACKGROUND MODEL----------
    */
    else
    {
        String<unsigned int> kmerCounts;
        MarkovModel<TUnmaskedAlphabet, TValue> backgroundModel(score.bgModelOrder);
        countKmers(kmerCounts, backgroundModel, sequence, score.kmerSize);

        //countKmers(sequence,kmerCounts,backgroundModel,k);
        int nvals = length(kmerCounts); //Number of kmers
        int len1 = 0;
        for (int l = 0; l < nvals; l++)
        {
            len1 += kmerCounts[l];
        }
        resize(standardisedCounts, nvals, (TValue) 0.0);

        String<TValue> probabilities;

        resize(probabilities, nvals, missing);

        TIteratorUnsigned itCounts;
        TIteratorTString itStandardisedCounts;

        itCounts = begin(kmerCounts);
        itStandardisedCounts = begin(standardisedCounts);
        //double sumTMP=0;
        for (; itCounts < end(kmerCounts); ++itCounts)
        {
            TValue p_w = 1;   //Probability of kmer
            TValue variance = 0;
            String<TUnmaskedAlphabet> w;
            unhash(w, (unsigned)position(itCounts), score.kmerSize);
            p_w = emittedProbability(backgroundModel, w);
            variance = ((TValue) pow(((TValue) len1) * p_w, 0.5));
            //Calculate standardised kmer Count
            //std::cout<<"\nword:"<<w<<", p_w:"<<p_w<<", var: "<<((TValue) pow(((TValue) len1)*p_w,0.5));
            if ((variance > pow(10, -10)) && (variance < pow(10, 10)))
            {
                if (p_w > 0)
                {
                    value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts)) - p_w * ((TValue)len1)) / ((TValue) pow(((TValue) len1) * p_w, 0.5));
                }
                ++itStandardisedCounts;
            }
        }

    }
}

}  // namespace seqan

#endif  // SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_D2STAR_H_
