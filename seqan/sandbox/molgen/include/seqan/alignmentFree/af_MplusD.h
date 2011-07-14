/*
 * af_MPLUSD.h
 *
 *  Created on: 19 Nov 2009
 *      Author: goeke
 */

#ifndef AF_MPLUSD_H_
#define AF_MPLUSD_H_

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TValue, typename TStringSet>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, AF_Score<MplusD> const & score)
{
    typedef typename Value<TStringSet>::Type TString;
    typedef typename Value<TString>::Type TAlphabet;
    typedef Matrix<TValue, 2> TMatrix;
    typedef typename Size<TMatrix>::Type TSize;
    typedef typename Iterator<TStringSet const>::Type       TIteratorSet;
    typedef typename Iterator<StringSet<String<TValue> > >::Type        TIteratorSetTValue;

    TSize seqNumber = length(sequenceSet);

    //resize the distMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    //fill(scoreMatrix,(TValue) 0);
    resize(scoreMatrix, (TValue) 0);

    StringSet<String<TValue> > weightedCounts;
    resize(weightedCounts, seqNumber);

    TIteratorSetTValue itWeightedCounts = begin(weightedCounts);

    TIteratorSet itSeqSet = begin(sequenceSet);
    for (; itSeqSet < end(sequenceSet); ++itSeqSet)
    {
        weightedFrequencies(value(itWeightedCounts), value(itSeqSet), score);
        ++itWeightedCounts;
    }
    std::cout << "\ncounted words";
    //calculate all pairwise scores and store them in scoreMatrix
    for (unsigned int rowIndex = 0; rowIndex < (seqNumber - 1); ++rowIndex)
    {
        std::cout << "\nSequence number " << rowIndex;
        for (unsigned int colIndex = rowIndex + 1; colIndex < (seqNumber); ++colIndex)
        {
            alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), weightedCounts[rowIndex], weightedCounts[colIndex], score);

            //alignmentFreeCompareCounts(weightedCounts[rowIndex],weightedCounts[colIndex],value(scoreMatrix,rowIndex,colIndex), AF_kmerMM());
            value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  //Copy symmetric entries
        }
    }

}

/*
 * Calculate pairwise score given the counts of all kmers
 */
template <typename TValue, typename TString>
void
alignmentFreeCompareCounts(TValue & result, TString & kmerCounts1, TString & kmerCounts2, AF_Score<MplusD> const & score)
{
    typedef typename Value<TString>::Type TStringValue;
    typedef typename Iterator<TString>::Type       TStringIterator;
    TStringIterator it1 = begin(kmerCounts1);
    TStringIterator it2 = begin(kmerCounts2);
    int kmerNumber = length(kmerCounts1);
    result = 0;

    /*
    std::ofstream myfile;

    //Set this to 1 for output of kmer scores for every comparison
    bool singleOutput=1;
    if(singleOutput==1)
    {
        myfile.open ("model1Motif_afkmerMMKmerWeights.tab",std::ios::app);
    }
*/
    for (; it1 < end(kmerCounts1); ++it1)
    {
        TValue tmp = result;
        if (value(it1) != 0) //else definition as zero, nothing to add
        {
            result += (value(it1) * ((TValue) log((((TValue) 2.0) * value(it1)) / (value(it1) + value(it2)))));
        }
        if (value(it2) != 0) //else definition as zero, nothing to add
        {
            result += (value(it2) * ((TValue) log((((TValue) 2.0) * value(it2)) / (value(it2) + value(it1)))));
        }
        tmp = result - tmp;
        /*
        if(singleOutput==1)
        {
            myfile << tmp<<"\t";
        }
*/
        ++it2;
    }
    /*
    if(singleOutput==1)
    {
        myfile << "\n";
        myfile.close();
    }
    */
    result = result / kmerNumber + ((TValue) (2.0 * log(2)));
    //Correction term? calculated with R based on the mplusd program, error +-10^-14? Original publication is unclear here. Algorithm not reproducable
    //val=val*0.40959999987592+0.818468190977198;
}

/*
 * count kmers and standardise count vectors for Dna5 and markov model background
 */
template <typename TValue, typename TString>
void weightedFrequencies(String<TValue> & weightedFrequencies, TString const & sequence, AF_Score<MplusD> const & score)
{
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;

    typedef typename Iterator<String<unsigned int>, Rooted>::Type       TIteratorUnsigned;
    typedef typename Iterator<String<TValue>, Rooted>::Type     TIteratorTValue;

    String<unsigned int> kmerCounts;

    if (score.bgModelOrder == 0)
    {


        String<TValue> backgroundFrequencies;
        countKmers(kmerCounts, backgroundFrequencies, sequence, score.kmerSize);
        int nvals = length(kmerCounts); //Number of kmers
        int len1 = 0;
        for (int l = 0; l < nvals; l++)
        {
            len1 += kmerCounts[l];
        }
        //fill(weightedFrequencies,nvals,(TValue) 0.0);
        resize(weightedFrequencies, nvals, (TValue) 0.0);

        TIteratorUnsigned itCounts;
        TIteratorTValue itWeightedFrequencies;

        itCounts = begin(kmerCounts);
        itWeightedFrequencies = begin(weightedFrequencies);

        for (; itCounts < end(kmerCounts); ++itCounts)
        {
            TValue p_w = 1;   //Probability of kmer

            String<TUnmaskedAlphabet> w;
            unhash(w, (unsigned)position(itCounts), score.kmerSize);
            calculateProbability(p_w, w, backgroundFrequencies);
            if (p_w > 0)
            {
                value(itWeightedFrequencies) = ((TValue) value(itCounts)) / ((TValue) len1) * p_w;
            }
            ++itWeightedFrequencies;
        }
    }
    else
    {

        MarkovModel<TUnmaskedAlphabet, TValue> backgroundModel(score.bgModelOrder);
        countKmers(kmerCounts, backgroundModel, sequence, score.kmerSize);

        int nvals = length(kmerCounts); //Number of kmers
        int len1 = 0;
        for (int l = 0; l < nvals; l++)
        {
            len1 += kmerCounts[l];
        }
        //fill(weightedFrequencies,nvals,(TValue) 0.0);
        resize(weightedFrequencies, nvals, (TValue) 0.0);

        TIteratorUnsigned itCounts;
        TIteratorTValue itWeightedFrequencies;

        itCounts = begin(kmerCounts);
        itWeightedFrequencies = begin(weightedFrequencies);

        for (; itCounts < end(kmerCounts); ++itCounts)
        {
            TValue p_w = 1; //Probability of kmer

            String<TUnmaskedAlphabet> w;
            unhash(w, (unsigned)position(itCounts), score.kmerSize);
            p_w = emittedProbability(backgroundModel, w);
            if (p_w > 0)
            {
                value(itWeightedFrequencies) = ((TValue) value(itCounts)) / ((TValue) len1) * p_w;
            }
            ++itWeightedFrequencies;
        }
    }
}

}
#endif /* AF_MPLUSD_H_ */