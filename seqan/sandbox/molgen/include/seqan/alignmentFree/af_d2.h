/*
 * af_d2.h
 *
 *  Created on: 7 Dec 2009
 *      Author: goeke
 */

#ifndef AF_D2_H_
#define AF_D2_H_


namespace SEQAN_NAMESPACE_MAIN
{

template <typename TStringSet, typename TValue>
void _alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, AF_Score<D2> const & score) //AF_Score<D2> const & score
{

    typedef typename Value<TStringSet>::Type TString;
    typedef typename Value<TString>::Type TAlphabet;
    typedef Matrix<TValue, 2> TMatrix;
    //typedef typename Size<TMatrix>::Type TSize;
    typedef typename Iterator<TStringSet const>::Type       TIteratorSet;
    typedef typename Iterator<StringSet<String<unsigned int> > >::Type              TIteratorSetInt;

    unsigned seqNumber = length(sequenceSet);

    //resize the ScoreMatrix
    setLength(scoreMatrix, 0, seqNumber);
    setLength(scoreMatrix, 1, seqNumber);
    resize(scoreMatrix, (TValue) 0);

    StringSet<String<unsigned> > kmerCounts;
    //StringSet<String<int> > kmerCounts;
    resize(kmerCounts, seqNumber);

    //Count all kmers
    TIteratorSetInt itKmerCounts = begin(kmerCounts);
    TIteratorSet itSeqSet = begin(sequenceSet);

    for (; itSeqSet < end(sequenceSet); ++itSeqSet)
    {
        countKmers(value(itKmerCounts), value(itSeqSet), score.kmerSize);
        ++itKmerCounts;
    }
    std::cout << "\ncounted words";
    //calculate all pairwise scores and store them in scoreMatrix
    for (unsigned int rowIndex = 0; rowIndex < (seqNumber); ++rowIndex) //(remove diagonal: seqNumber-1)
    {
        std::cout << "\nSequence number " << rowIndex;
        for (unsigned int colIndex = rowIndex; colIndex < (seqNumber); ++colIndex) //(remove diagonal: rowIndex+1)
        {
            alignmentFreeCompareCounts(value(scoreMatrix, rowIndex, colIndex), kmerCounts[rowIndex], kmerCounts[colIndex], score);
            value(scoreMatrix, colIndex, rowIndex) = value(scoreMatrix, rowIndex, colIndex);  //Copy symmetric entries
        }
    }

}

/*
 * Calculate pairwise score given the counts of all kmers
 */
template <typename TValue>
void
alignmentFreeCompareCounts(TValue & result, String<unsigned int> & kmerCounts1, String<unsigned int> & kmerCounts2, AF_Score<D2> const & score)
{
    typedef typename Iterator<String<unsigned int> >::Type      TIteratorInt;

    TIteratorInt it1 = begin(kmerCounts1);
    TIteratorInt it2 = begin(kmerCounts2);

    result = 0;
    for (; it1 < end(kmerCounts1); ++it1)
    {
        result += value(it1) * value(it2);
        ++it2;
    }
}

}
#endif /* AF_D2_H_ */