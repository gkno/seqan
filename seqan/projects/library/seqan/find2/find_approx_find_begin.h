/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007-2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
 ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ============================================================================
  A simple implementation of the algorithms and datastructures required for
  findBegin().  It is a "backwards" implementation of Sellers' algorithm for
  approxim prefix search.

  In order to support findBegin(), the approximate DP-based algorithms need
  to store more information than only for the forward search.  If findBegin()
  does not need to be supported, this information does not have to be stored.
  For this reason, the functionality is centralized in the struct template
  _ApproxFindBegin that is to be used as the base class for the approximate
  search patterns.

  If the parameter HasFindBeginSupport is set to True then the class gets the
  information and functionality for findBegin(), if it is set to False, it
  does not get this information.
 ============================================================================
*/
#ifndef SEQAN_FIND2_FIND_APPROX_DPSEARCH_
#define SEQAN_FIND2_FIND_APPROX_DPSEARCH_

namespace seqan {

// TODO(holtgrew): Add FindPrefix support.
// TODO(holtgrew): Add support for affine gap costs.
template <typename TNeedle, typename TScore, typename HasFindBeginSupport>
struct _ApproxFindBegin;


// Implementation for no findBegin() support.
template <typename TNeedle, typename TScore>
struct _ApproxFindBegin<TNeedle, TScore, False> {};


// Implementation for findBegin() support.
template <typename TNeedle, typename TScore>
struct _ApproxFindBegin<TNeedle, TScore, True> : public _FindState {
    // We use the position of the needle for the finder, too.  This is
    // the best we can do with the split-up interface.
    typedef typename Position<TNeedle>::Type TPosition;
    // We will use a matrix column of score values for the dynamic
    // programming and types for this.
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TScoreValue> TMatrixColumn;

    // The length of the match for findBegin().
    TPosition _findBeginMatchLength;

    // The score limit when searching backwards.
    int _findBeginScoreLimit;

    // The current score of the findBegin() call.
    int _findBeginCurrentScore;

    // The matrix column for the findBegin computation.
    TMatrixColumn _findBeginMatrixColumn;

    // The state of the begin search.
    TState _findBeginState;
};


// Called to initialize a _ApproxFindBegin datastructure.
template <typename TNeedle, typename TScore>
void _initFindBegin(_ApproxFindBegin<TNeedle, TScore, True> & findBeginStruct,
                    typename Position<TNeedle>::Type const & needleLength, TScore const & scoringScheme, int scoreLimit) {
    SEQAN_CHECKPOINT;
    typedef _ApproxFindBegin<TNeedle, TScore, True> TApproxFindBegin;
    typedef typename Position<TNeedle>::Type TPosition;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TScoreValue> TMatrixColumn;
    typedef typename Iterator<TMatrixColumn>::Type TIterator;

    SEQAN_ASSERT_EQ(scoreGapOpen(scoringScheme), scoreGapOpen(scoringScheme), "findBegin() only supports linear gap costs at the moment.");

    // Assumption: The end position must have been found already!  The
    // caller has to ensure this.
    findBeginStruct._findBeginState = TApproxFindBegin::STATE_FOUND;
    findBeginStruct._findBeginMatchLength = 0;
    findBeginStruct._findBeginScoreLimit = scoreLimit;
    // Re-initialize the matrix column.
    resize(findBeginStruct._findBeginMatrixColumn, needleLength);
    TScoreValue gapValue = scoreGap(scoringScheme);
    for (TIterator it = begin(findBeginStruct._findBeginMatrixColumn, Standard()); it != end(findBeginStruct._findBeginMatrixColumn, Standard()); ++it) {
        *it = gapValue;
        gapValue += scoreGap(scoringScheme);
    }
}


// Actual implementation of findBegin() through _ApproxFindBegin.
template <typename TNeedle, typename TScore, typename THaystack>
bool _findBeginImpl(_ApproxFindBegin<TNeedle, TScore, True> & findBeginStruct,
                    TScore const & scoringScheme,
                    typename Position<TNeedle>::Type const & endPosition,
                    THaystack const & haystack, TNeedle const & needle) {
    SEQAN_CHECKPOINT;
    typedef _ApproxFindBegin<TNeedle, TScore, True> TApproxFindBegin;
    typedef typename Position<TNeedle>::Type TPosition;
    typedef typename Value<TScore>::Type TScoreValue;
    SEQAN_ASSERT_EQ(scoreGapOpen(scoringScheme), scoreGapExtend(scoringScheme), "findBegin() only supports linear gap costs at the moment.");

    // Finding the start of the match follows Seller's algorithm but
    // works from the right to the left.

//     std::cout << "== start of findBegin() == " << std::endl;
//     for (unsigned i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i)
//         std::cout << findBeginStruct._findBeginMatrixColumn[i] << " ";
//     std::cout << std::endl;
//     std::cout << "=====" << std::endl;

    // Search for the next match.
    for (TPosition & j = findBeginStruct._findBeginMatchLength; j < endPosition; ++j) {
        // Compute best score if the begin of the match was j position
        // left of the end.
        findBeginStruct._findBeginCurrentScore = (j + 1) * scoreGap(scoringScheme);  // is v in SeqAn book
        TScoreValue d = (j) * scoreGap(scoringScheme);
        for (TPosition i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i) {
            TScoreValue h = findBeginStruct._findBeginMatrixColumn[i];
            findBeginStruct._findBeginCurrentScore = _max(d + score(scoringScheme, haystack[endPosition - j - 1], needle[length(needle) - 1 - i]),
                                                          _max(findBeginStruct._findBeginCurrentScore, h) + scoreGap(scoringScheme));
            findBeginStruct._findBeginMatrixColumn[i] = findBeginStruct._findBeginCurrentScore;
            d = h;
        }
//         std::cout << "== within of findBegin() == " << std::endl;
//         for (unsigned i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i)
//             std::cout << findBeginStruct._findBeginMatrixColumn[i] << " ";
//         std::cout << std::endl;
//         std::cout << "=====" << std::endl;
        if (findBeginStruct._findBeginCurrentScore >= findBeginStruct._findBeginScoreLimit) {
//             std::cout << "== end of findBegin() == " << std::endl;
//             for (unsigned i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i)
//                 std::cout << findBeginStruct._findBeginMatrixColumn[i] << " ";
//             std::cout << std::endl;
//             std::cout << "=====" << std::endl;
            findBeginStruct._findBeginMatchLength += 1;
            // Found a match, update the state and report the match.
            findBeginStruct._findBeginState = TApproxFindBegin::STATE_BEGIN_FOUND;
            return true;
        }
    }

//     std::cout << "== end of findBegin() == " << std::endl;
//     for (unsigned i = 0; i < length(findBeginStruct._findBeginMatrixColumn); ++i)
//         std::cout << findBeginStruct._findBeginMatrixColumn[i] << " ";
//     std::cout << std::endl;
//     std::cout << "=====" << std::endl;

    // No match found, update state accordingly.
    findBeginStruct._findBeginState = TApproxFindBegin::STATE_BEGIN_NOTFOUND;
	return false;
}

};

#endif  // SEQAN_FIND2_FIND_APPROX_DPSEARCH_
