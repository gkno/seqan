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
// This file contains the function that calls alignment free sequence 
// comparisons algorithms (see AFScore):
// alignmentFreeComparison()
// ==========================================================================

#ifndef SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_ALIGNMENT_FREE_COMPARISON_H_
#define SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_ALIGNMENT_FREE_COMPARISON_H_

namespace seqan {

/**
.Function.alignmentFreeComparison:
..summary:Computes the pairwise similarity scores for a set of sequences
..cat:Alignment Free
..signature:alignmentFreeComparison(scoreMatrix, sequenceSet, score)
..param.scoreMatrix:A two-dimensional Matrix, used to store all pairwise scores
...type:Class.Matrix
..param.sequenceSet:StringSet containing all sequences for which pairwise scores will be computed
...type:Class.StringSet
..param.score:The score values to be used for computing the alignment.
...type:Class.AFScore
*/

template <typename TStringSet, typename TValue, typename TComparisonMethod>
void alignmentFreeComparison(Matrix<TValue, 2> & scoreMatrix, TStringSet const & sequenceSet, TComparisonMethod const & comparisonMethod)
{
    _alignmentFreeComparison(scoreMatrix, sequenceSet, comparisonMethod);
}

}  // namespace seqan

#endif  // SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_COMPARISON_H_
