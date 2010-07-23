/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2010

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

  Based on the code by Carsten Kemena <carsten.kemena@crg.es>, debugged
  by Birte Kehr <birte.kehr@fu-berlin.de>.
 ============================================================================
  Seed extension algorithms.

  The approach for gapped X-drop extension is based on the algorithm in
  Figure 2 from (Zhang et al., 2000).

    Zhang Z, Schwartz S, Wagner L, Miller W.  A greedy algorithm for aligning
    DNA sequences.  Journal of computational biologyA a journal of
    computational molecular cell biology.  2000;7(1-2):203-14.
    doi:10.1089/10665270050081478
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_EXTENSION_H_
#define SEQAN_SEEDS_SEEDS_EXTENSION_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Tag.Seed Extension
..cat:Seed Handling
..summary:The algorithms used to extend a seed.
..see:Function.extendSeed
..see:Function.extendSeeds
..see:Function.extendSeedScore
..see:Function.extendSeedsScore
..tag.MatchExtend:Extends a seed until a mismatch occurs.
..tag.UngappedXDrop:Ungapped extension of a seed until score drops below a Value.
..tag.GappedXDrop:Gapped extension of a seed until score drops below a Value. Only @Spec.SimpleSeed@s.
..include:seqan/seeds.h
*/
struct _MatchExtend;
typedef Tag<_MatchExtend> const MatchExtend;

struct _UnGappedXDrop;
typedef Tag<_UnGappedXDrop> const UnGappedXDrop;

struct _GappedXDrop;
typedef Tag<_GappedXDrop> const GappedXDrop;

enum ExtensionDirection
{
    EXTEND_LEFT,
    EXTEND_RIGHT,
    EXTEND_BOTH
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

/**
.Function.extendSeed
..summary:Extends a seed.
..cat:Seed Handling
..signature:extendSeed(seed, query, database, direction, MatchExtend)
..signature:extendSeed(seed, query, database, direction, scoreDropOff, scoreMatrix, {UngappedXDrop, GappedXDrop})
..param.seed: The seed to extend.
...type:Class.Seed
..param.query: The query sequence.
...type:Class.String
..param.query: The database sequence.
...type:Class.String
..param.direction: Defines the direction in which the seed should be extended. 0 = left, 1 = right, 2 = both
..param.scoreDropOff: The score drop after which the extension should stop. The extension stops if this value is exceeded.
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.scoreMatrix: The scoring scheme.
...type:Spec.Simple Score
...remarks:Only used for the algorithms @Tag.Seed Extension.UngappedXDrop@ and @Tag.Seed Extension.GappedXDrop@
..param.tag: The algorithm to use.
...type:Tag.Seed Extension.MatchExtend
...type:Tag.Seed Extension.UngappedXDrop
...type:Tag.Seed Extension.GappedXDrop
..include:seqan/seeds.h
*/

// We need one specialization for each combination of the extension
// variants and seeds.  It is not worth to extract the common parts
// for simple and chained seeds.

template <typename TConfig, typename TQuery, typename TDatabase>
inline void 
extendSeed(Seed<Simple, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
		   MatchExtend const &)
{
    // For match extension of Simple Seeds, we can simply update the
    // begin and end values in each dimension.
	SEQAN_CHECKPOINT;

    typedef Seed<Simple, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
		TPosition posDim0 = getBeginDim0(seed) ;
		TPosition posDim1 = getBeginDim1(seed);
		while (posDim0 >= 1 && posDim1 >= 1 && query[posDim0 - 1] == database[posDim1 - 1]) {
			--posDim0;
			--posDim1;
		}
		setBeginDim0(seed, posDim0);
		setBeginDim1(seed, posDim1);
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
		TPosition posDim0 = getEndDim0(seed) ;
		TPosition posDim1 = getEndDim1(seed);
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && query[posDim0] == database[posDim1]) {
			++posDim0;
			++posDim1;
		}
		setEndDim0(seed, posDim0);
		setEndDim1(seed, posDim1);
	}
}


template <typename TConfig, typename TQuery, typename TDatabase>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
		   MatchExtend const &)
{
    // For match extension of Chained Seeds, we extend the first and
    // the last Seed Diagonal.
	SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(seed), 0u);

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Value<TSeed>::Type TSeedDiagonal;
    typedef typename Position<TSeedDiagonal>::Type TPosition;
    typedef typename Size<TSeedDiagonal>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        TSeedDiagonal & diag = front(seed);
		TPosition posDim0 = diag.beginDim0;
		TPosition posDim1 = diag.beginDim1;
        TSize diagonalLength = diag.length;
		while (posDim0 >= 1 && posDim1 >= 1 && query[posDim0 - 1] == database[posDim1 - 1]) {
			--posDim0;
			--posDim1;
            ++diagonalLength;
		}
        diag.beginDim0 = posDim0;
        diag.beginDim1 = posDim1;
        diag.length = diagonalLength;
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
        TSeedDiagonal & diag = back(seed);
		TPosition posDim0 = diag.beginDim0 + diag.length;
		TPosition posDim1 = diag.beginDim1 + diag.length;
        TSize diagonalLength = diag.length;
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && query[posDim0] == database[posDim1]) {
			++posDim0;
			++posDim1;
            ++diagonalLength;
		}
        diag.length = diagonalLength;
	}
}


template <typename TConfig, typename TQuery, typename TDatabase, typename TScoreValue, typename TScoreSpec>
inline void 
extendSeed(Seed<Simple, TConfig> & seed,
           TQuery const & query,
           TDatabase const & database,
           ExtensionDirection direction,
           Score<TScoreValue, TScoreSpec> const & scoringScheme,
           TScoreValue scoreDropOff,
           UnGappedXDrop const &)
{
    // For ungapped X-drop extension of Simple Seeds, we can simply
    // update the begin and end values in each dimension.
	SEQAN_CHECKPOINT;

    scoreDropOff = -scoreDropOff;

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        TScoreValue tmpScore = 0;
		TPosition posDim0 = getBeginDim0(seed) ;
		TPosition posDim1 = getBeginDim1(seed);
        TPosition mismatchingSuffixLength = 0;
		while (posDim0 >= 1 && posDim1 >= 1 && tmpScore > scoreDropOff) {
            tmpScore += score(scoringScheme, posDim0, posDim1, query, database);
            if (query[posDim0 - 1] == database[posDim1 - 1]) {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            } else {
                mismatchingSuffixLength += 1;
            }
			--posDim0;
			--posDim1;
		}
		setBeginDim0(seed, posDim0 + mismatchingSuffixLength);
		setBeginDim1(seed, posDim1 + mismatchingSuffixLength);
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
        TScoreValue tmpScore = 0;
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
		TPosition posDim0 = getEndDim0(seed) ;
		TPosition posDim1 = getEndDim1(seed);
        TPosition mismatchingSuffixLength = 0;
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && tmpScore > scoreDropOff) {
            tmpScore += score(scoringScheme, posDim0, posDim1, query, database);
            if (query[posDim0] == database[posDim1]) {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            } else {
                mismatchingSuffixLength += 1;
            }
            ++posDim0;
            ++posDim1;
		}
		setEndDim0(seed, posDim0 - mismatchingSuffixLength);
		setEndDim1(seed, posDim1 - mismatchingSuffixLength);
    }

    // TODO(holtgrew): Update score?!
}


template <typename TConfig, typename TQuery, typename TDatabase, typename TScoreValue, typename TScoreSpec>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
           Score<TScoreValue, TScoreSpec> const & scoringScheme,
           TScoreValue scoreDropOff,
		   UnGappedXDrop const &)
{
    // For ungapped X-drop extension of Chained Seeds, we extend the
    // first and the last Seed Diagonal.
	SEQAN_CHECKPOINT;

    scoreDropOff = -scoreDropOff;

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Value<TSeed>::Type TSeedDiagonal;
    typedef typename Position<TSeedDiagonal>::Type TPosition;
    typedef typename Size<TSeedDiagonal>::Type TSize;
    
	// Extension to the left
	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        TScoreValue tmpScore = 0;
        TPosition mismatchingSuffixLength = 0;
        TSeedDiagonal & diag = front(seed);
		TPosition posDim0 = getBeginDim0(seed) ;
		TPosition posDim1 = getBeginDim1(seed);
        TSize diagonalLength = diag.length;
		while (posDim0 >= 1 && posDim1 >= 1 && tmpScore > scoreDropOff) {
            tmpScore += score(scoringScheme, posDim0, posDim1, query, database);
            if (query[posDim0 - 1] == database[posDim1 - 1]) {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            } else {
                mismatchingSuffixLength += 1;
            }
			--posDim0;
			--posDim1;
            ++diagonalLength;
		}
        diag.beginDim0 = posDim0 + mismatchingSuffixLength;
        diag.beginDim1 = posDim1 + mismatchingSuffixLength;
        diag.length = diagonalLength - mismatchingSuffixLength;
	}

	// Extension to the right
	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
        TScoreValue tmpScore = 0;
        TPosition mismatchingSuffixLength = 0;
		TSize lengthDim0 = length(query);
		TSize lengthDim1 = length(database);
        TSeedDiagonal & diag = back(seed);
		TPosition posDim0 = diag.beginDim0 + diag.length;
		TPosition posDim1 = diag.beginDim1 + diag.length;
        TSize diagonalLength = diag.length;
		while (posDim0 < lengthDim0 && posDim1 < lengthDim1 && tmpScore > scoreDropOff) {
            tmpScore += score(scoringScheme, posDim0, posDim1, query, database);
            if (query[posDim0] == database[posDim1]) {
                mismatchingSuffixLength = 0;
                if (tmpScore > static_cast<TScoreValue>(0))
                    tmpScore = 0;
            } else {
                mismatchingSuffixLength += 1;
            }
            ++posDim0;
            ++posDim1;
            ++diagonalLength;
		}
        diag.length = diagonalLength - mismatchingSuffixLength;
	}

    // TODO(holtgrew): Update score?!
}


// Helper function for gapped X-drop extension for Simple Seeds.
template <typename TConfig, typename TQuerySegment, typename TDatabaseSegment, typename TScoreValue>
inline void
_extendSeedGappedXDropOneDirection(
        Seed<Simple, TConfig> & seed,
        TQuerySegment const & querySeg,
        TDatabaseSegment const & databaseSeg,
        ExtensionDirection direction,
        Score<TScoreValue, Simple> const & scoringScheme,
        TScoreValue scoreDropOff)
{
    SEQAN_CHECKPOINT;
	
	typedef Seed<Simple, TConfig> TSeed;
	typedef __int64 TDiagonal;
	typedef __int64 TPosition;
	typedef __int64 TScore;
	
	
	// TODO(holtgrew): Adapt to new begin/end interface, requires a +1 on each right border.
    TScore gapCost = scoreGap(scoringScheme);
    // TODO(holtgrew): Why is this +1-gapCost here?
    TPosition infimum = infimumValue<TPosition>()+1-gapCost;
	
    // TODO(holtgrew): rename to upperDiagonal/lowerDiagonal? Initialize to seed's values?
    TPosition upperBound = 0;
    TPosition lowerBound = 0;
	
    // TODO(holtgrew): Rename to lengthDim{0,1}?
    TPosition xLength = length(querySeg);
    TPosition yLength = length(databaseSeg);
	
    // set variables for calculating the sequence positions
    int factor, xSummand, ySummand;
    if(direction == 0) {
        factor = -1;
        xSummand = xLength;
        ySummand = yLength;
    } else {
        factor = 1;
        xSummand = -1;
        ySummand = -1;
    }
	
    // antidiagonals of DP matrix
    String<TPosition> antiDiagonal1;
    String<TPosition> antiDiagonal2;
    String<TPosition> antiDiagonal3;
	
    fill(antiDiagonal1, 1, 0);
    fill(antiDiagonal2, 2, infimum);
    fill(antiDiagonal3, _min(3, xLength+1), infimum);
	
    // TODO(holtgrew): Is there a better way to do this? Probably not...
    String<TPosition> *antiDiag1 = &antiDiagonal1;	//smallest diagonal
	String<TPosition> *antiDiag2 = &antiDiagonal2;
	String<TPosition> *antiDiag3 = &antiDiagonal3;	//current diagonal
	String<TPosition> *tmpDiag;
	
	//Matrix initialization
	if (gapCost >= -scoreDropOff) {
		(*antiDiag2)[0] = gapCost;
		(*antiDiag2)[1] = gapCost;
	}
	if (2 * gapCost >= -scoreDropOff) {
		(*antiDiag3)[0] = 2*gapCost;
		(*antiDiag3)[2] = 2*gapCost;
	}
	
    // TODO(holtgrew): Probably rename to reflect the letters from (Zhang et al., 2000) and the SeqAn book?
	TPosition b = 1; // lower bound for i
	TPosition u = 0; // upper bound for i
	TPosition k = 1; // current antidiagonal
    // TODO(holtgrew): Maybe move declarations of the tmp variables to their usage scopes?
    // TODO(holtgrew): Probably rename these to best, M_i_j?
	TPosition tmp;   // for calculating the maximum for a single matrix entry
	TPosition tmpMax1 = 0; // maximum score without the current diagonal
	TPosition tmpMax2 = 0; // maximum score including the current diagonal
	
	//Extension as proposed by Zhang et al, following Figure 2.
	while(true) {
		++k;
		for (TPosition i = b; i <= u+1; ++i) {
            // calculate matrix entry
			tmp = infimum;  // TODO(holtgrew): Superflous? Move declaration of tmp in here.
			tmp = _max((*antiDiag2)[i-1], (*antiDiag2)[i]) + gapCost;
			tmp = _max(tmp, (*antiDiag1)[i-1] + score(scoringScheme, xSummand + factor*i, ySummand + factor*(k-i), querySeg, databaseSeg));
			
			tmpMax2 = _max(tmpMax2, tmp);
			if (tmp < tmpMax1-scoreDropOff)
				(*antiDiag3)[i] = infimum;
			else
				(*antiDiag3)[i] = tmp;
		}
        
        // narrow the relevant matrix region
		while ((b < (TPosition)length(*antiDiag3)-1) && ((*antiDiag3)[b] == infimum) && (*antiDiag2)[b-1] == infimum) {
			++b;
		}
		++u;
		while ((u >= 0) && ((*antiDiag3)[u] == infimum) && (*antiDiag2)[u] == infimum) {
			--u;
        }
		
		//borders for lower triangle of edit matrix
		b = _max(b, k-yLength);
		u = _min(u, xLength-1);
		
        if (b > u+1) break;
		
	    // Calculate upper/lower bound for diagonals.  The modulo operation is used here because of the "1/2 shift" between adjacent antidiagonals.  Zhang solves this by iterating in 1/2 steps.
        // TODO(holtgrew): Rename to lowerDiagonal/upperDiagonal?
        if (2*((k+1)/2 - b) - (k%2) > lowerBound) {
			lowerBound = 2*((k+1)/2 - b) - (k%2);
        }
        if (2*(u - k/2) - (k%2) > upperBound) {
            upperBound = 2*(u - k/2) - (k%2);
        }
		
        // swap diagonals
		tmpDiag = antiDiag1;
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag3 = tmpDiag;
		
        // extend last diagonal diagonal to be the new longest and initialize
		// TODO(holtgrew): Maybe limit diagonal size.  This requires an offset computation.
		int d = 0;
        // TODO(holtgrew): Replace the following while-loop with fill(*andiDiag3, _min(length(*antiDiag3) + 3, xLength), 0)?.
		while ((d < 3) && ((TPosition)length(*antiDiag3) <= xLength)) {
			appendValue(*antiDiag3, 0);
			++d;
		}
        // TODO(holtgrew): What about clear(), fill() or std::fill()?
		for (unsigned int eu = 0; eu < length(*antiDiag3); ++eu)
			(*antiDiag3)[eu] = infimum;
		
		if ((*antiDiag2)[0]+ gapCost >= tmpMax1-scoreDropOff)
			(*antiDiag3)[0] = (*antiDiag2)[0] + gapCost;
		if ((*antiDiag2)[length(*antiDiag2)-1] + gapCost >= tmpMax1-scoreDropOff)
			(*antiDiag3)[length(*antiDiag3)-1] = (*antiDiag2)[length(*antiDiag2)-1] + gapCost;
		
		tmpMax1 = tmpMax2;
	}
	
	//  // print anti diagonals
	//  for(int ii = length(*antiDiag1)-1; ii >= 0 ; --ii) {
	//      for(int jj = ii; jj > 0 ; --jj) std::cout << "  ";
	//      std::cout << " ";
	//      if ((*antiDiag1)[ii] == infimum ) std::cout << "i" << " ";
	//      else std::cout << (*antiDiag1)[ii] << " ";
	//      if (length(*antiDiag2) <= ii+1 ) std::cout << " ";
	//      else if ((*antiDiag2)[ii+1] == infimum ) std::cout << "i" << " ";
	//      else std::cout << (*antiDiag2)[ii+1] << " ";
	//if (length(*antiDiag3) <= ii+2 ) std::cout << std::endl;
	//      else if ((*antiDiag3)[ii+2] == infimum ) std::cout << "i" << std::endl;
	//      else std::cout << (*antiDiag3)[ii+2] << std::endl;
	//  }
	
	//Find seed start/end
	// TODO(holtgrew): With some more work at this point, this could maybe be simplified and made more robust.
	TPosition extLengthQuery = 0; // length of extension in query
    TPosition extLengthDatabase = 0; // length of extension in database
	TPosition tmpMax = infimum;
	if ((k >= xLength+yLength) && ((*antiDiag2)[u+1] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of both sequences
		extLengthQuery = xLength;
        extLengthDatabase = yLength;
		tmpMax = (*antiDiag2)[u+1];
    } else if ((b >= xLength) && b < (TPosition)length(*antiDiag2) && ((*antiDiag2)[b] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of query
        tmpMax = (*antiDiag2)[b];
        extLengthQuery = xLength;
        extLengthDatabase = k-(xLength+1);
    } else if ((k-u-1 >= yLength) && u >= 0 && ((*antiDiag2)[u] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of database
        tmpMax = (*antiDiag2)[u];
        extLengthQuery = u;
        extLengthDatabase = yLength;
    } else {
        for (unsigned int eu = 0; eu < length(*antiDiag1); ++eu) {
            if ((*antiDiag1)[eu] > tmpMax) {
                // extension ends with mismatch
		        tmpMax = (*antiDiag1)[eu];
		        extLengthQuery = _min(xLength, (TPosition)eu);
                extLengthDatabase = _min(yLength, (TPosition)(k - (eu+2)));
		    }
        }
	}	
	
	/*

    typedef Seed<Simple, TConfig> TSeed;
//    typedef typename Position<TSeed>::Type TPosition;
//    typedef typename Size<TSeed>::Type TSize;
//    typedef typename Diagonal<TSeed>::Type TDiagonal;
    typedef typename Diagonal<TSeed>::Type TPosition;
    typedef typename Diagonal<TSeed>::Type TSize;
    typedef typename Diagonal<TSeed>::Type TDiagonal;
	
    // Can either extend left or right.
    SEQAN_ASSERT_NEQ(direction, EXTEND_BOTH);

    // TODO(holtgrew): Adapt to new begin/end interface, requires a +1 on each right border.
    TScoreValue gapCost = scoreGap(scoringScheme);
    // TODO(holtgrew): Why is this +1-gapCost here?
    TScoreValue infimum = infimumValue<TScoreValue>() + 1 - gapCost;

    // TODO(holtgrew): rename to upperDiagonal/lowerDiagonal? Initialize to seed's values?
    TPosition upperBound = 0;
    TPosition lowerBound = 0;

    // TODO(holtgrew): Rename to lengthDim{0,1}?
    TPosition xLength = length(querySeg);
    TPosition yLength = length(databaseSeg);

    // Set variables for calculating the sequence positions depending
    // on the extension direction.  This way, we do not have to write
    // this function once for right and left.
    int factor, xSummand, ySummand;
    if (direction == EXTEND_LEFT) {
        factor = -1;
        xSummand = xLength;
        ySummand = yLength;
    } else {  // direction == EXTEND_RIGHT
        factor = 1;
        xSummand = -1;
        ySummand = -1;
    }

    // Declare and initialize the antidiagonals of DP matrix.
    String<TScoreValue> antiDiagonal1;
    String<TScoreValue> antiDiagonal2;
    String<TScoreValue> antiDiagonal3;
    fill(antiDiagonal1, 1, 0);
    fill(antiDiagonal2, 2, infimum);
    fill(antiDiagonal3, _min(3u, xLength + 1), infimum);

    // We keep pointers to the antidiagonals so we can rotate them
    // later on.  antiDiag1 is the smallest/leftmost one while
    // antiDiag3 is the current/rightmost one.
    String<TScoreValue> * antiDiag1 = &antiDiagonal1;
	String<TScoreValue> * antiDiag2 = &antiDiagonal2;
	String<TScoreValue> * antiDiag3 = &antiDiagonal3;
	String<TScoreValue> * tmpDiag;

	// Matrix initialization.
	if (gapCost >= -scoreDropOff) {
		(*antiDiag2)[0] = gapCost;
		(*antiDiag2)[1] = gapCost;
	}
	if (2 * gapCost >= -scoreDropOff) {
		(*antiDiag3)[0] = 2 * gapCost;
		(*antiDiag3)[2] = 2 * gapCost;
	}
		
    // TODO(holtgrew): Probably rename to reflect the letters from (Zhang et al., 2000) and the SeqAn book?
	TPosition b = 1; // lower bound for i
	TPosition u = 0; // upper bound for i
	TPosition k = 1; // current antidiagonal
    // TODO(holtgrew): Maybe move declarations of the tmp variables to their usage scopes?
    // TODO(holtgrew): Probably rename these to best, M_i_j?
	TScoreValue tmp;   // for calculating the maximum for a single matrix entry
	TScoreValue tmpMax1 = 0; // maximum score without the current diagonal
	TScoreValue tmpMax2 = 0; // maximum score including the current diagonal

	// Extension as proposed in Figure 2 of (Zhang et al., 2000).
	while(true) {
		++k;
		for (TPosition i = b; i <= u+1; ++i) {
            // calculate matrix entry
			tmp = infimum;  // TODO(holtgrew): Superflous? Move declaration of tmp in here.
			tmp = _max((*antiDiag2)[i-1], (*antiDiag2)[i]) + gapCost;
			tmp = _max(tmp, (*antiDiag1)[i-1] + score(scoringScheme, xSummand + factor*i, ySummand + factor*(k-i), querySeg, databaseSeg));

			tmpMax2 = _max(tmpMax2, tmp);
			if (tmp < tmpMax1-scoreDropOff)
				(*antiDiag3)[i] = infimum;
			else
				(*antiDiag3)[i] = tmp;
		}
        
        // Narrow the relevant matrix region.
		while ((b < (TPosition)length(*antiDiag3)-1) && ((*antiDiag3)[b] == infimum) && (*antiDiag2)[b-1] == infimum)
			++b;
		++u;
		while ((u >= 1) && ((*antiDiag3)[u] == infimum) && (*antiDiag2)[u] == infimum)
			--u;
			
		// Borders for lower triangle of edit matrix.
		b = _max(b, k > yLength ? k-yLength : 0);
		u = _min(u, xLength-1);

        if (b > u + 1)
            break;

	    // Calculate upper/lower bound for diagonals.  The modulo
	    // operation is used here because of the "1/2 shift" between
	    // adjacent antidiagonals.  Zhang et al. solve this by
	    // iterating in 1/2 steps.
        if (2 * ((k + 1) / 2 - b) - (k % 2) > lowerBound)
			lowerBound = 2 * ((k + 1) / 2 - b) - (k % 2);
        if (2 * (u - k / 2) - (k % 2) > upperBound)
            upperBound = 2 * (u - k / 2) - (k % 2);

        // Cycle-swap pointers to anti-diagonals.
		tmpDiag = antiDiag1;
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag3 = tmpDiag;

        // Extend last anti-diagonal to be the new longest and
        // initialize with infimum values.
        //
		// TODO(holtgrew): Maybe limit diagonal size.  This requires an offset computation.
		int d = 0;
        // TODO(holtgrew): Replace the following while-loop with fill(*andiDiag3, _min(length(*antiDiag3) + 3, xLength), 0)?.
		while ((d < 3) && ((TPosition)length(*antiDiag3) <= xLength)) {
			appendValue(*antiDiag3, 0);
			++d;
		}
        // TODO(holtgrew): What about clear(), fill() or std::fill()?
		for (unsigned int eu = 0; eu < length(*antiDiag3); ++eu)
			(*antiDiag3)[eu] = infimum;

		if ((*antiDiag2)[0]+ gapCost >= tmpMax1-scoreDropOff)
			(*antiDiag3)[0] = (*antiDiag2)[0] + gapCost;
		if ((*antiDiag2)[length(*antiDiag2)-1] + gapCost >= tmpMax1-scoreDropOff)
			(*antiDiag3)[length(*antiDiag3)-1] = (*antiDiag2)[length(*antiDiag2)-1] + gapCost;

		tmpMax1 = tmpMax2;
	}

	// Find start and end of the seed.
    //
	// TODO(holtgrew): With some more work at this point, this could maybe be simplified and made more robust.
	TPosition extLengthQuery = 0; // length of extension in query
    TPosition extLengthDatabase = 0; // length of extension in database
	TScoreValue tmpMax = infimum;
	if ((k >= xLength+yLength) && ((*antiDiag2)[u+1] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of both sequences
		extLengthQuery = xLength;
        extLengthDatabase = yLength;
		tmpMax = (*antiDiag2)[u+1];
    } else if ((b >= xLength) && b < (TPosition)length(*antiDiag2) && ((*antiDiag2)[b] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of query
        tmpMax = (*antiDiag2)[b];
        extLengthQuery = xLength;
        extLengthDatabase = k-(xLength+1);
    } else if ((k-u-1 >= yLength) && u >= 0 && ((*antiDiag2)[u] >= tmpMax1-scoreDropOff)) {
        // extension ends at end of database
        tmpMax = (*antiDiag2)[u];
        extLengthQuery = u;
        extLengthDatabase = yLength;
    } else {
        for (unsigned int eu = 0; eu < length(*antiDiag1); ++eu) {
            if ((*antiDiag1)[eu] > tmpMax) {
                // extension ends with mismatch
		        tmpMax = (*antiDiag1)[eu];
		        extLengthQuery = _min(xLength, (TPosition)eu);
                extLengthDatabase = _min(yLength, (TPosition)(k - (eu+2)));
		    }
        }
	}*/

    // If we reached a higher score then update the seed to reflect
    // the extension.
    if (tmpMax != infimum) {
        if (direction == EXTEND_LEFT) {
            // Set lower and upper diagonals.
            if (getLowerDiagonal(seed) < getStartDiagonal(seed) + static_cast<TDiagonal>(upperBound))
                setLowerDiagonal(seed, getStartDiagonal(seed) + upperBound);
            if (getUpperDiagonal(seed) > getStartDiagonal(seed) - static_cast<TDiagonal>(lowerBound))
                setUpperDiagonal(seed, getStartDiagonal(seed) - lowerBound);
            
            // Set new start position of seed.
            setBeginDim0(seed, getBeginDim0(seed) - extLengthQuery);
            setBeginDim1(seed, getBeginDim1(seed) - extLengthDatabase);
        } else {  // direction == EXTEND_RIGHT
            // Set new lower and upper diagonals.
            if (getUpperDiagonal(seed) > getEndDiagonal(seed) - static_cast<TDiagonal>(upperBound))
                setUpperDiagonal(seed, getEndDiagonal(seed) - upperBound);
            if (getLowerDiagonal(seed) < getEndDiagonal(seed) + static_cast<TDiagonal>(lowerBound))
                setLowerDiagonal(seed, getEndDiagonal(seed) + lowerBound);
            
            // Set new end position of seed.
            setEndDim0(seed, getEndDim0(seed) + extLengthQuery);
            setEndDim1(seed, getEndDim1(seed) + extLengthDatabase);
        }
    }
}


template <typename TConfig, typename TQuery, typename TDatabase, typename TScoreValue>
inline void 
extendSeed(Seed<Simple, TConfig> & seed, 
		   TQuery const & query,
		   TDatabase const & database,
		   ExtensionDirection direction,
           Score<TScoreValue, Simple> const & scoringScheme,
           TScoreValue scoreDropOff,
		   GappedXDrop const &)
{
    // For gapped X-drop extension of Simple Seeds, we can simply
    // update the begin and end values in each dimension.
	SEQAN_CHECKPOINT;

    typedef Seed<Simple, TConfig> TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;

    // The algorithm only works for linear gap scores < 0, mismatch scores < 0
    // and match scores > 0.
    SEQAN_ASSERT_GT(scoreMatch(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreMismatch(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreGapOpen(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreGapExtend(scoringScheme), 0);
    SEQAN_ASSERT_EQ(scoreGapExtend(scoringScheme), scoreGapOpen(scoringScheme));

	if (direction == EXTEND_LEFT || direction == EXTEND_BOTH) {
        // Do not extend to the left if we are already at the beginning of an
        // infix or the sequence itself.
        //
        // TODO(holtgrew): Can this be handled more elegantly, maybe in _extendSeedGappedXDropOneDirection?
        if (getBeginDim0(seed) != beginPosition(query) && getBeginDim1(seed) != beginPosition(database)) {
            typedef typename Prefix<TQuery const>::Type TQueryPrefix;
            typedef typename Prefix<TDatabase const>::Type TDatabasePrefix;

            TQueryPrefix queryPrefix = prefix(query, getBeginDim0(seed));
            TDatabasePrefix databasePrefix = prefix(database, getBeginDim1(seed));
            _extendSeedGappedXDropOneDirection(seed, queryPrefix, databasePrefix, EXTEND_LEFT, scoringScheme, scoreDropOff);
        }
    }

	if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH) {
        // Do not extend to the right if we are already at the beginning of an
        // infix or the sequence itself.
        //
        // TODO(holtgrew): Can this be handled more elegantly, maybe in _extendSeedGappedXDropOneDirection?
        if (getEndDim0(seed) < endPosition(query) && (getEndDim1(seed) < endPosition(database))) {
            typedef typename Suffix<TQuery const>::Type TQuerySuffix;
            typedef typename Suffix<TDatabase const>::Type TDatabaseSuffix;
            
            TQuerySuffix querySuffix = suffix(query, getEndDim0(seed));
            TDatabaseSuffix databaseSuffix = suffix(database, getEndDim1(seed));
			// std::cout << "database = " << database << std::endl;
			// std::cout << "database Suffix = " << databaseSuffix << std::endl;
			// std::cout << "query = " << query << std::endl;
			// std::cout << "query Suffix = " << querySuffix << std::endl;
            _extendSeedGappedXDropOneDirection(seed, querySuffix, databaseSuffix, EXTEND_RIGHT, scoringScheme, scoreDropOff);
        }
    }

    // TODO(holtgrew): Update seed's score?!
}


template <typename TConfig, typename TQuery, typename TDatabase, typename TScoreValue>
inline void 
extendSeed(Seed<ChainedSeed, TConfig> & seed,
		   TQuery const & /*query*/,
		   TDatabase const & /*database*/,
		   ExtensionDirection /*direction*/,
           Score<TScoreValue, Simple> const & scoringScheme,
           TScoreValue /*scoreDropOff*/,
		   GappedXDrop const &)
{
    // For ungapped X-drop extension of Chained Seeds, we have to append
    // diagonals to the front and end of the list of seed diagonals and modify
    // the first and last one of the current set of seed diagonals.
	SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(seed), 0u);

    typedef Seed<ChainedSeed, TConfig> TSeed;
    typedef typename Value<TSeed>::Type TSeedDiagonal;
    typedef typename Position<TSeedDiagonal>::Type TPosition;
    typedef typename Size<TSeedDiagonal>::Type TSize;
    
    // The algorithm only works for linear gap scores < 0, mismatch scores < 0
    // and match scores > 0.
    SEQAN_ASSERT_GT(scoreMatch(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreMismatch(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreGapOpen(scoringScheme), 0);
    SEQAN_ASSERT_LT(scoreGapExtend(scoringScheme), 0);
    SEQAN_ASSERT_EQ(scoreGapExtend(scoringScheme), scoreGapOpen(scoringScheme));

    SEQAN_ASSERT_FAIL("Write me!");

    // TODO(holtgrew): Update seed's score?!
    // TODO(holtgrew): For chained seeds, the code is similar to the code for simple seeds.  However, we need to store the whole matrix to compute the traceback, from this the edit script and from this compute how to add/update diagonals.
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_EXTENSION_H_
