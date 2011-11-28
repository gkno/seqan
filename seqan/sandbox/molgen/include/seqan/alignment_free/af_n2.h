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
// This header contains the implementation of the N2 score for alignment free
// sequence comparison with word neighbourhood counts
// Goeke et al, to appear.
// These functions can be called with alignmentFreeComparison().
// ==========================================================================
#ifndef SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_N2_H_
#define SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_N2_H_

namespace seqan {

String<unsigned> revComIndex_N;

void initialiseRevComIndexN(unsigned k)
{
    unsigned myLength = (unsigned) pow(4, k);
    // resize(revComIndex_N,myLength);
    // fill(revComIndex_N,myLength, 0);
    resize(revComIndex_N, myLength, 0);
    Shape<Dna, SimpleShape> myShape;
    // Declare variables
    // TShape myShape;		// Shape, length can be changed (kmer_length)
    resize(myShape, k);
    for (unsigned i = 0; i < myLength; ++i)
    {
        String<Dna> w;
        unhash(w, i, k);
        DnaStringReverseComplement wRC(w);
        unsigned hashValue = hash(myShape, begin(wRC));
        revComIndex_N[i] = hashValue;
    }

}

// Currently only one mismatch maximum
StringSet<String<unsigned int> > kmerNeighbourhoodN;
// // template < typename TAlphabet, unsigned MISMATCHNUMBER>
// void initialiseKmerNeighbourhood_N2(unsigned const mismatchNumber, unsigned k)
// {
// typedef Dna TAlphabet;
// // std::cout<<"mismatches:"<<mismatchNumber;
// unsigned const tmp=value(mismatchNumber);
// typedef Enumerator< String<TAlphabet>, EditEnvironment<HammingDistance, (tmp)> > TEnumHamming;
//  // typedef Enumerator< DnaString, EditEnvironment<HammingDistance, MISMATCHNUMBER> > TEnumHamming;
// // typedef Enumerator< DnaString, TEditEnvronment> TEnumHamming;
//  unsigned myLength=(unsigned) pow(4,k);
//  Shape<Dna,SimpleShape >	myShape;
//  resize(myShape,k);
//  resize(kmerNeighbourhoodN,myLength);
//  for(unsigned i=0;i<1;++i)// myLength
//  {
//      String<TAlphabet> w;
//      unhash(w, i,k);
//      TEnumHamming enumHamming(w);
//          ::std::cout << w<<"\n";
//          // Iterator<TEnumHamming>::Type it = begin(enumHamming);
// // typename Iterator <TEnumHamming>::Type TIterator;
// // TIterator it;// = begin(enumHamming);
// // it= begin(enumHamming);
// //       // Iterator<TEnumHamming>::Type itEnd = end(enumHamming);
// // Iterator<TEnumHamming> itEnd = end(enumHamming);
// //       for (; it != itEnd; goNext(it))
// //       {
// //           ::std::cout << (*it) << " ";
// //           ::std::cout << ::std::endl;
// //       }
//  }
// }
void initialiseKmerNeighbourhoodN(unsigned k, bool revCom)
{
    unsigned myLength = (unsigned) pow(4, k);
    Shape<Dna, SimpleShape> myShape;
    resize(myShape, k);
    resize(kmerNeighbourhoodN, myLength);
    for (unsigned i = 0; i < myLength; ++i)
    {
        resize(kmerNeighbourhoodN[i], 1, i);

        String<Dna> w;
        unhash(w, i, k);
        // appendValue(kmerNeighbourhoodN[i],i);
        if ((revComIndex_N[i] != i) && (revCom == true))
        {
            appendValue(kmerNeighbourhoodN[i], revComIndex_N[i]);
        }
        for (unsigned j = 0; j < k; ++j)
        {

// std::cout<<"\n\n"<<w;
            for (unsigned l = 0; l < 4; ++l)
            {
                String<Dna> wTMP;
                wTMP = w;
                if (wTMP[j] != l)
                {
                    wTMP[j] = l;
                    // std::cout<<wTMP<<"\n";
                    unsigned hashValue = hash(myShape, begin(wTMP));
                    // check for double word occurences
                    bool duplicate = false;
                    if (revCom == true)
                    {

                        for (unsigned n = 0; n < length(kmerNeighbourhoodN[i]); ++n)
                        {
                            if ((hashValue) == kmerNeighbourhoodN[i][n])
                            {
                                duplicate = true;
                                // std::cout<<"\nbreak;"<<wTMP<<" w:"<<w;
                                break;
                            }
                        }

                    }

                    if (duplicate == false)
                    {
                        appendValue(kmerNeighbourhoodN[i], hashValue);
                        if (revCom == true)
                        {
                            if (revComIndex_N[hashValue] != hashValue)
                            {
                                appendValue(kmerNeighbourhoodN[i], revComIndex_N[hashValue]);
                            }
                        }
                    }
                }
            }
            // unsigned hashValue=hash(myShape,begin(w));
        }
    }
// // output of kmer neigbours
//  for(unsigned i = 0;i< length(kmerNeighbourhoodN[0]);++i)
//  {
//          String<Dna> w2;
//          unhash(w2, kmerNeighbourhoodN[0][i],k);
//          std::cout<<"\n"<<w2<<kmerNeighbourhoodN[0][i];
//  }
}

template <typename TValue, typename TStringSet>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, AFScore<N2> const & score)
{



    typedef typename Value<TStringSet>::Type TString;
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;
    typedef Matrix<TValue, 2> TMatrix;

    typedef typename Iterator<TStringSet const>::Type       TIteratorSet;
    typedef typename Iterator<StringSet<String<unsigned int> > >::Type      TIteratorSetInt;
    typedef typename Iterator<StringSet<String<double> > >::Type        TIteratorSetDouble;


    // Initialise reverse complement hash table
    initialiseRevComIndexN(score.kmerSize);
    if (score.revCom == "bothStrands")
    {
        initialiseKmerNeighbourhoodN(score.kmerSize, true);
    }
    else
    {
        initialiseKmerNeighbourhoodN(score.kmerSize, false);
    }
// unsigned myMismatches=score.mismatches;
// EditEnvironment<HammingDistance, 0> myEdit;
// initialiseKmerNeighbourhoodN2<Dna,1>(score.kmerSize);
    // typedef typename Iterator<String<MarkovModel<Dna> > > ::Type		TIteratorMarkovModel;

    unsigned int seqNumber = length(sequenceSet);

    // resize the distMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    // fill(scoreMatrix,(TValue) 0);
    resize(scoreMatrix, (TValue) 0);


    StringSet<String<double> > standardisedKmerCounts;
    resize(standardisedKmerCounts, seqNumber);
    // Count all kmers and all background nucleotide frequencies and store them in StringSets
    TIteratorSetDouble itStandardisedKmerCounts = begin(standardisedKmerCounts);
    TIteratorSet itSeqSet = begin(sequenceSet);
    for (; itSeqSet < end(sequenceSet); ++itSeqSet)
    {

        standardiseCounts(value(itStandardisedKmerCounts), value(itSeqSet), score);
	
        if(score.verbose)
	{
		std::cout << "\n" << position(itSeqSet);
	}
	++itStandardisedKmerCounts;
    }

    if (score.norm == true) // normalise score so that sequence-self-comparison is maximal (1)
    {
        itStandardisedKmerCounts = begin(standardisedKmerCounts);
        for (; itStandardisedKmerCounts < end(standardisedKmerCounts); ++itStandardisedKmerCounts)
        {
            TValue normValue = 0.0;
            for (unsigned i = 0; i < length(value(itStandardisedKmerCounts)); ++i)
            {
                normValue += value(itStandardisedKmerCounts)[i] * value(itStandardisedKmerCounts)[i];
            }
            for (unsigned i = 0; i < length(value(itStandardisedKmerCounts)); ++i)
            {
                value(itStandardisedKmerCounts)[i] /= sqrt(normValue);
            }

        }
    }
    // output of pairwise kmer weight to file, will be large very fast!
    std::ofstream myfile;
    if (score.outputFile != "")
    {
        myfile.open(toCString(score.outputFile));
        for (unsigned i = 0; i < length(standardisedKmerCounts[0]); ++i)
        {
            String<TUnmaskedAlphabet> w;
            unhash(w, i, score.kmerSize);
            // std::ofstream myfile;

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

        if(score.verbose)
	{
		std::cout << "\ncounted words";
	}

    // calculate all pairwise scores and store them in scoreMatrix
    for (unsigned int rowIndex = 0; rowIndex < (seqNumber); ++rowIndex) //(remove diagonal values seqNumber-1)
    {
	
        if(score.verbose)
	{
		std::cout << "\nSequence number " << rowIndex;
	}
        for (unsigned int colIndex = rowIndex; colIndex < (seqNumber); ++colIndex) // remove diagonal values rowIndex+1
        {

            alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), standardisedKmerCounts[rowIndex], standardisedKmerCounts[colIndex], score);
            value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  // Copy symmetric entries
        }
    }
}

/*
 * Calculate pairwise score given the counts of all kmers
 */
template <typename TValue, typename TString>
void
alignmentFreeCompareCounts(TValue & result, TString & kmerCounts1, TString & kmerCounts2, AFScore<N2> const & score)
{
    typedef typename Value<TString>::Type TStringValue;
    typedef typename Iterator<TString, Rooted>::Type        TIteratorTString;

    TIteratorTString it1 = begin(kmerCounts1);
    TIteratorTString it2 = begin(kmerCounts2);
    result = 0.0;
    TValue resultRC = 0.0;
    // output of pairwise comparison kmerWeights to file
    for (; it1 < end(kmerCounts1); ++it1)
    {

	result += (TValue)(value(it1) * value(it2));

// std::cout<<"\nval1:"<<value(it1)<<", val2:"<<value(it2)<<", prod:"<<(value(it1)*value(it2))<<", score: "<<result;

        // Computation of the reverse complement strand score
        if ((score.revCom != "") && (score.revCom != "bothStrands"))
        {
            unsigned hashValue = revComIndex_N[position(it1)];
            resultRC += (TValue)(value(it1) * kmerCounts2[hashValue]);
        }
        ++it2;
    }


    // std::cout<<"\nresult"<<result<<"\n"<<"\nresultRC"<<resultRC<<"\n";

    if (score.revCom == "mean")
    {
        result = (TValue) (resultRC + result) / 2;
    }
    else if (score.revCom == "max")
    {
        result = std::max(resultRC, result);
    }
    else if (score.revCom == "min")
    {
        result = std::min(resultRC, result);
    }


}

/*
 * count kmers and standardise count vectors for Dna5 and markov model background
 */
template <typename TString, typename TSequence>
void standardiseCounts(TString & standardisedCounts, TSequence const & sequence, AFScore<N2> const & score)
{


    typedef typename Value<TSequence>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;
    typedef typename Value<TString>::Type TValue;
    typedef typename Iterator<String<unsigned int>, Rooted>::Type       TIteratorUnsigned;
    typedef typename Iterator<TString, Rooted>::Type        TIteratorTString;
    unsigned int alphabetSize = ValueSize<TUnmaskedAlphabet>::VALUE;

    Matrix<TValue, 2> covarianceMatrix;
    TValue missing = -pow(10, 10);
    if (score.mismatches > 0)
    {
        setLength(covarianceMatrix, 0, pow(alphabetSize, score.kmerSize));
        setLength(covarianceMatrix, 1, pow(alphabetSize, score.kmerSize));
        // fill(covarianceMatrix, missing);
        resize(covarianceMatrix, missing);
    }

    /*
    *--------ORDER 0 BACKGROUND MODEL----------
    */
    if (score.bgModelOrder == 0)
    {

        String<unsigned int> kmerCounts;
        String<double> backgroundFrequencies;
        countKmers(kmerCounts, backgroundFrequencies, sequence, score.kmerSize);
// int tmpSum=0;for(int tmp=0;tmp<length(kmerCounts);++tmp){tmpSum+=kmerCounts[tmp];}std::cout<<"tmpSum\n"<<tmpSum;
        int nvals = length(kmerCounts); // Number of kmers
        int len1 = 0;
        for (int l = 0; l < nvals; l++)
        {
            len1 += kmerCounts[l];
        }
        // fill(standardisedCounts,nvals,(TValue) 0.0);
        resize(standardisedCounts, nvals, (TValue) 0.0);


        // string of tvalue to store p_w
        String<TValue> probabilities;
        // fill(probabilities,nvals,missing);
        resize(probabilities, nvals, missing);


        // temporary counter for mismatch kmer counting
        // unsigned counterTMP;

        TIteratorUnsigned itCounts;
        TIteratorTString itStandardisedCounts;

        itCounts = begin(kmerCounts);
        itStandardisedCounts = begin(standardisedCounts);

        for (; itCounts < end(kmerCounts); ++itCounts)
        {
// temporary counter for mismatch kmer counting
            TValue counterTMP = 0;
            TValue p_w = 1;   // Probability of kmer

            String<TUnmaskedAlphabet> w;
            unhash(w, (unsigned)position(itCounts), score.kmerSize);
            calculateProbability(p_w, w, backgroundFrequencies);
            TValue variance = 0;
/*
*	--------EXACT VARIANCE CALCULATION---------
*	used for mismatches, reversecomplement,
*/
            // if((score.exactVariance==true)||(score.revCom=="bothStrands")||(score.mismatches==1))
            //{
            if ((score.mismatches == 1)) // mismatch calculation
            {

                p_w = 0;
                // The first word in the kmerNeighbourhoodN is the kmer itself, it is weighted normally!
                // unsigned row=0;
                // counterTMP+=kmerCounts[kmerNeighbourhoodN[position(itCounts)][row]];
                // Sum of all entries in the covariance matrix. only once computed
                unsigned wordHash = position(itCounts);
                unsigned wordRCHash = revComIndex_N[wordHash];

                for (unsigned row = 0; row < length(kmerNeighbourhoodN[wordHash]); ++row)
                {
                    unsigned wordRowHash = kmerNeighbourhoodN[wordHash][row];
                    //// The kmer itself is weighted normally!
                    if (wordRowHash == wordHash) //(row==0)//(kmerNeighbourhoodN[position(itCounts)]==position(itCounts))//(row==0)// The first word in the kmerNeighbourhoodN is the kmer itself, it is weighted normally!
                    {
                        counterTMP += (TValue) kmerCounts[wordRowHash];
                    }
                    else if ((score.revCom == "bothStrands") && (wordRowHash == wordRCHash))
                    {
                        counterTMP += ((TValue) kmerCounts[wordRowHash]);
                    }
                    else
                    {
                        counterTMP += ((TValue) kmerCounts[wordRowHash]) * score.mismatchWeight;
                    }
                    String<Dna> wMM1;
                    unhash(wMM1, wordRowHash, score.kmerSize);
// std::cout<<"\n"<<wMM1;

                    for (unsigned col = row; col < length(kmerNeighbourhoodN[wordHash]); ++col)
                    {
                        unsigned wordColHash = kmerNeighbourhoodN[wordHash][col];
                        if (value(covarianceMatrix, wordColHash, wordRowHash) == missing)
                        {

                            String<Dna> wMM2;

                            unhash(wMM2, wordColHash, score.kmerSize);
                            calculateCovariance(value(covarianceMatrix, wordColHash, wordRowHash), wMM1, wMM2, backgroundFrequencies, (len1 + score.kmerSize - 1));
                            value(covarianceMatrix, wordRowHash, wordColHash) = value(covarianceMatrix, wordColHash, wordRowHash);
                        }

                        if (row == col) // Varaince of weighted variables
                        {
                            /*if((score.revCom=="bothStrands")&&(wordRowHash==wordRCHash)&&(wordHash==wordRCHash))
                            {
                            // std::cout<<"check";
                            variance+=(4.0)*value(covarianceMatrix,wordRowHash,wordColHash);
                            }// vierfache varianz für palindromische wörter da cov(x,x)=var(x)!
                            else
*/
                            if ((wordRowHash == wordHash) || (score.revCom == "bothStrands" && (wordRowHash == wordRCHash))) //(row==0)// variance of the kmer is counted full
                            {
                                variance += value(covarianceMatrix, wordRowHash, wordColHash);
                            }
                            else
                            {
                                variance += pow(score.mismatchWeight, 2) * value(covarianceMatrix, wordRowHash, wordColHash);
                            }
                        }
// the covariance of the kmer and the reverse complement is weighted full
                        else if ((score.revCom == "bothStrands") && (((wordRowHash == wordHash) && (wordColHash == wordRCHash)) || ((wordRowHash == wordRCHash) && (wordColHash == wordHash))))
                        {
                            variance += (2.0) * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                        else if ((wordRowHash == wordHash || wordColHash == wordHash) || (score.revCom == "bothStrands" && (wordRowHash == wordRCHash || wordColHash == wordRCHash))) // cov is weighted half
                        {
                            variance += (2.0) * score.mismatchWeight * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                        else    // covariance is weighted^2
                        {
                            variance += (2.0) * pow(score.mismatchWeight, 2) * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                    }

                    if (probabilities[wordRowHash] == missing)
                    {
                        calculateProbability(probabilities[wordRowHash], wMM1, backgroundFrequencies);
                    }

                    if (wordRowHash == wordHash) //(row==0)	// Weight the probabiliets /expected values, normal weight for the kmer itself!
                    {
                        p_w += probabilities[wordRowHash];
                    }
                    else if ((score.revCom == "bothStrands") && (wordRowHash == wordRCHash)) //(row==0)	// Weight the probabiliets /expected values, normal weight for the rev.com. kmer itself!
                    {
                        p_w += probabilities[wordRowHash];
                    }
                    else
                    {
                        p_w += score.mismatchWeight * probabilities[wordRowHash];
                    }
                }
                variance = pow(variance, 0.5);
// std::cout<<"\nw: "<<w<<", p: "<<p_w<<", var: "<<variance;
            }    // mismatch calculation end
            else if (score.revCom == "bothStrands")
            {

                TValue variance1;
                TValue variance2;
                TValue covariance;
                // DnaStringReverseComplement wRC(w);
                String<Dna> wRC;
// std::cout<<"\n";//<<revComIndex_N[(unsigned)position(itCounts)];
                unhash(wRC, (unsigned)  revComIndex_N[(unsigned)position(itCounts)], score.kmerSize);
                calculateVariance(variance1, w, backgroundFrequencies, (len1 + score.kmerSize - 1));
                calculateVariance(variance2, wRC, backgroundFrequencies, (len1 + score.kmerSize - 1));
                calculateCovariance(covariance, w, wRC, backgroundFrequencies, (len1 + score.kmerSize - 1));
                // variance=pow((variance1+pow(score.revComWeight,2)*variance2+(2.0)*score.revComWeight*covariance),0.5);
                variance = pow((variance1 + variance2 + (2.0) * covariance), 0.5);
                TValue p_wRC = 1;   // Probability of rev com kmer
                calculateProbability(p_wRC, wRC, backgroundFrequencies);
                p_w += p_wRC;
                // value(itCounts)+=kmerCounts[revComIndex_N[(unsigned)position(itCounts)]];
            }
            else
            {
                calculateVariance(variance, w, backgroundFrequencies, (len1 + score.kmerSize - 1));
                variance = pow(variance, 0.5);
            }
            //}
            // else
            //{
            //	variance=((TValue) pow(((TValue) len1)*p_w,0.5));
            //}

            if ((variance > pow(10, -10)) && (variance < pow(10, 10)))
            {
                if (p_w > 0)
                {

                    if (score.mismatches > 0)
                    {

                        value(itStandardisedCounts) = ((TValue) ((TValue) counterTMP) - p_w * ((TValue)len1)) / variance;
                    }
                    else if (score.revCom == "bothStrands")
                    {
                        // value(itStandardisedCounts)=((TValue) ((TValue) value(itCounts)+score.revComWeight*kmerCounts[revComIndex_N[(unsigned)position(itCounts)]])-p_w*((TValue)len1))/variance;
                        value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts) + kmerCounts[revComIndex_N[(unsigned)position(itCounts)]]) - p_w * ((TValue)len1)) / variance;

                        // value(itStandardisedCounts)=(((TValue) value(itCounts))-p_w*((TValue)len1))/((TValue) pow(((TValue) len1)*p_w,0.5));
                    }
                    else
                    {
                        value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts)) - p_w * ((TValue)len1)) / variance;
                    }
                }
            }
// std::cout<<"\nword:"<<w<<",count: "<<value(itCounts)<<",stcount: "<<value(itStandardisedCounts)<<", p_w:"<<p_w<<", var(approx): "<<((TValue) pow(((TValue) len1)*p_w,0.5))<<"var (ex)"<<variance;

// std::cout<<"\n"<<w<<"; counts:"<<value(itCounts)<<" counterTMP: "<<counterTMP<<"; stCounts: "<<value(itStandardisedCounts)<<" p_w: "<<p_w<<" var:"<<variance;
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

        // countKmers(sequence,kmerCounts,backgroundModel,k);
        int nvals = length(kmerCounts); // Number of kmers
        int len1 = 0;
        for (int l = 0; l < nvals; l++)
        {
            len1 += kmerCounts[l];
        }
        // fill(standardisedCounts,nvals,(TValue) 0.0);
        resize(standardisedCounts, nvals, (TValue) 0.0);

        String<TValue> probabilities;

        // fill(probabilities,nvals,missing);
        resize(probabilities, nvals, missing);

        TIteratorUnsigned itCounts;
        TIteratorTString itStandardisedCounts;

        itCounts = begin(kmerCounts);
        itStandardisedCounts = begin(standardisedCounts);
        // double sumTMP=0;
        for (; itCounts < end(kmerCounts); ++itCounts)
        {
            TValue p_w = 1; // Probability of kmer
            TValue variance = 0;
            String<TUnmaskedAlphabet> w;
            unhash(w, (unsigned)position(itCounts), score.kmerSize);
            p_w = emittedProbability(backgroundModel, w);

            TValue counterTMP = 0.0;
// if((score.exactVariance==true)||(score.revCom=="bothStrands")||(score.mismatches==1))
            //{

            if ((score.mismatches == 1)) // mismatch calculation
            {

                p_w = 0;
                // The first word in the kmerNeighbourhoodN is the kmer itself, it is weighted normally!
                // unsigned row=0;
                // counterTMP+=kmerCounts[kmerNeighbourhoodN[position(itCounts)][row]];
                // Sum of all entries in the covariance matrix. only once computed
                unsigned wordHash = position(itCounts);
                unsigned wordRCHash = revComIndex_N[wordHash];

                for (unsigned row = 0; row < length(kmerNeighbourhoodN[wordHash]); ++row)
                {
                    unsigned wordRowHash = kmerNeighbourhoodN[wordHash][row];
                    //// The kmer itself is weighted normally!
                    if (wordRowHash == wordHash) //(row==0)//(kmerNeighbourhoodN[position(itCounts)]==position(itCounts))//(row==0)// The first word in the kmerNeighbourhoodN is the kmer itself, it is weighted normally!
                    {
                        counterTMP += (TValue) kmerCounts[wordRowHash];
                    }
                    else if ((score.revCom == "bothStrands") && (wordRowHash == wordRCHash))
                    {
                        counterTMP += ((TValue) kmerCounts[wordRowHash]);
                    }
                    else
                    {
                        counterTMP += ((TValue) kmerCounts[wordRowHash]) * score.mismatchWeight;
                    }
                    String<Dna> wMM1;
                    unhash(wMM1, wordRowHash, score.kmerSize);
// std::cout<<"\n"<<wMM1;
                    for (unsigned col = row; col < length(kmerNeighbourhoodN[wordHash]); ++col)
                    {
                        unsigned wordColHash = kmerNeighbourhoodN[wordHash][col];
                        if (value(covarianceMatrix, wordColHash, wordRowHash) == missing)
                        {

                            String<Dna> wMM2;

                            unhash(wMM2, wordColHash, score.kmerSize);
                            calculateCovariance(value(covarianceMatrix, wordColHash, wordRowHash), wMM1, wMM2, backgroundModel, (len1 + score.kmerSize - 1));
                            value(covarianceMatrix, wordRowHash, wordColHash) = value(covarianceMatrix, wordColHash, wordRowHash);
                        }

                        if (row == col) // Varaince of weighted variables
                        {
                            if ((wordRowHash == wordHash) || (score.revCom == "bothStrands" && (wordRowHash == wordRCHash))) //(row==0)// variance of the kmer is counted full
                            {
                                variance += value(covarianceMatrix, wordRowHash, wordColHash);
                            }
                            else
                            {
                                variance += pow(score.mismatchWeight, 2) * value(covarianceMatrix, wordRowHash, wordColHash);
                            }
                        }
// the covariance of the kmer and the reverse complement is weighted full
                        else if ((score.revCom == "bothStrands") && (((wordRowHash == wordHash) && (wordColHash == wordRCHash)) || ((wordRowHash == wordRCHash) && (wordColHash == wordHash))))
                        {
                            variance += (2.0) * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                        else if ((wordRowHash == wordHash || wordColHash == wordHash) || (score.revCom == "bothStrands" && (wordRowHash == wordRCHash || wordColHash == wordRCHash))) // cov is weighted half
                        {
                            variance += (2.0) * score.mismatchWeight * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                        else    // covariance is weighted^2
                        {
                            variance += (2.0) * pow(score.mismatchWeight, 2) * value(covarianceMatrix, wordRowHash, wordColHash);
                        }
                    }
                    if (probabilities[wordRowHash] == missing)
                    {
                        probabilities[wordRowHash] = emittedProbability(backgroundModel, wMM1);
                    }

                    if (wordRowHash == wordHash) //(row==0)	// Weight the probabiliets /expected values, normal weight for the kmer itself!
                    {
                        p_w += probabilities[wordRowHash];
                    }
                    else if ((score.revCom == "bothStrands") && (wordRowHash == wordRCHash)) //(row==0)	// Weight the probabiliets /expected values, normal weight for the rev.com. kmer itself!
                    {
                        p_w += probabilities[wordRowHash];
                    }
                    else
                    {
                        p_w += score.mismatchWeight * probabilities[wordRowHash];
                    }
                }

                variance = pow(variance, 0.5);
// std::cout<<"\nw: "<<w<<", p: "<<p_w<<", var: "<<variance;
            }    // mismatch calculation end
            else if (score.revCom == "bothStrands")
            {

                TValue variance1;
                TValue variance2;
                TValue covariance;
                // DnaStringReverseComplement wRC(w);
                String<Dna> wRC;
                unhash(wRC, (unsigned)  revComIndex_N[(unsigned)position(itCounts)], score.kmerSize);
                calculateVariance(variance1, w, backgroundModel, (len1 + score.kmerSize - 1));
// std::cout<<"\nword:"<<w<<", wrc: "<<wRC;
                calculateVariance(variance2, wRC, backgroundModel, (len1 + score.kmerSize - 1));
                calculateCovariance(covariance, w, wRC, backgroundModel, (len1 + score.kmerSize - 1));
                // variance=pow((variance1+pow(score.revComWeight,2)*variance2+(2.0)*score.revComWeight*covariance),0.5);
                variance = pow((variance1 + variance2 + (2.0) * covariance), 0.5);
                TValue p_wRC = 1;   // Probability of rev com kmer
                p_wRC = emittedProbability(backgroundModel, wRC);
                // std::cout<<"\npw1:"<<p_w<<"\npwRC:"<<p_wRC;
                p_w += p_wRC;
                // std::cout<<"\npw2:"<<p_w;	// value(itCounts)+=kmerCounts[revComIndex_N[(unsigned)position(itCounts)]];
            }
            else
            {
// std::cout<<"test";
                calculateVariance(variance, w, backgroundModel, (len1 + score.kmerSize - 1));
                variance = pow(variance, 0.5);
            }
            //}
            // sumTMP+=((TValue) pow(((TValue) len1)*p_w,0.5));
            // Calculate standardised kmer Count
            // std::cout<<"\nword:"<<w<<", p_w:"<<p_w<<", var: "<<((TValue) pow(((TValue) len1)*p_w,0.5));
            if ((variance > pow(10, -10)) && (variance < pow(10, 10)))
            {
                if (p_w > 0)
                {

                    if (score.mismatches > 0)
                    {
                        value(itStandardisedCounts) = ((TValue) ((TValue) counterTMP) - p_w * ((TValue)len1)) / variance;
                    }
                    else if (score.revCom == "bothStrands")
                    {
                        // value(itStandardisedCounts)=((TValue) ((TValue) value(itCounts)+score.revComWeight*kmerCounts[revComIndex_N[(unsigned)position(itCounts)]])-p_w*((TValue)len1))/variance;
                        value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts) + kmerCounts[revComIndex_N[(unsigned)position(itCounts)]]) - p_w * ((TValue)len1)) / variance;
                        // value(itStandardisedCounts)=(((TValue) value(itCounts))-p_w*((TValue)len1))/((TValue) pow(((TValue) len1)*p_w,0.5));
                    }
                    else
                    {
                        value(itStandardisedCounts) = ((TValue) ((TValue) value(itCounts)) - p_w * ((TValue)len1)) / ((TValue) variance);
                    }

// std::cout<<"\nword:"<<w<<",count: "<<value(itCounts)<<",stcount: "<<value(itStandardisedCounts)<<", p_w:"<<p_w<<", var(approx): "<<((TValue) pow(((TValue) len1)*p_w,0.5))<<"var (ex)"<<variance;
// std::cout<<"\nstCounts:"<<value(itStandardisedCounts);
                }
                ++itStandardisedCounts;
            }
        }
        // std::cout<<"\nsum: "<<sumTMP;
    }
}

}  // namespace seqan

#endif  // SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_AF_N2_H_
