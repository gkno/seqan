/*-------------------------------------------------------
 *-------------------------------------------------------
 *
 *This file contains functions for alignment free sequence comparisons
 *Can be used for pairwise sequence comparison
 *-------------------------------------------------------
 *-------------------------------------------------------
 */

#ifndef ALIGNMENTFREECOMPARISON_H_
#define ALIGNMENTFREECOMPARISON_H_

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Function.alignmentFreeComparison:
..summary:Computes the pairwise similarity scores or distances for a set of sequences
..cat:AlignmentFree
..signature:
alignmentFreeComparison(scoreMatrix, sequenceSet, score)
..param.scoreMatrix:A two-dimensional Matrix, used to store all pairwise scores
...type:Class.Matrix
..param.sequenceSet:StringSet containing all sequences for which pairwise scores will be computed
...type:Class.StringSet
..param.score:The score values to be used for computing the alignment.
...type:Class.AF_Score
*/
template <typename TStringSet, typename TValue, typename TComparisonMethod>
void alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, TComparisonMethod const & comparisonMethod)
{
    _alignmentFreeComparison(scoreMatrix, sequenceSet, comparisonMethod);

}

template <typename TStringSet, typename TValue, typename TComparisonMethod>
void alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, TComparisonMethod & comparisonMethod)
{
    _alignmentFreeComparison(scoreMatrix, sequenceSet, comparisonMethod);

}

}
#endif /* ALIGNMENTFREECOMPARISON_H_ */