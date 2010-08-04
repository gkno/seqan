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
 ============================================================================
  Global chaining algorithms.
 ==========================================================================*/

#ifndef SEQAN_SEEDS_SEEDS_GLOBAL_CHAINING_H_
#define SEQAN_SEEDS_SEEDS_GLOBAL_CHAINING_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Tag.Global Chaining
..cat:Seed Handling
..see:Function.chainSeeds
..tag:GusfieldChaining:
    Chaining as described in (Gusfield, 1997) section 13.3.
 */
struct _GusfieldChaining;
typedef Tag<_GusfieldChaining> GusfieldChaining;

// ===========================================================================
// Metafunctions
// ===========================================================================

// ===========================================================================
// Functions
// ===========================================================================

template <typename TTargetContainer, typename TSeedSpec, typename TSeedSetSpec, typename TSeedConfig>
void
chainSeedsGlobally(
        TTargetContainer & target,
        SeedSet<TSeedSpec, TSeedSetSpec, TSeedConfig> const & seedSet,
        GusfieldChaining const &)
{
    SEQAN_CHECKPOINT;

    typedef SeedSet<TSeedSpec, TSeedSetSpec, TSeedConfig> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;
    typedef typename Position<TSeed>::Type TPosition;
    typedef typename Size<TSeed>::Type TSize;
    typedef typename TSeedSet::THighQualitySeeds const THighQualitySeeds;
    typedef typename THighQualitySeeds::const_iterator THighQualitySeedsIterator;

    // -----------------------------------------------------------------------
    // Step 1: Generate the sorted list of interval points.
    // -----------------------------------------------------------------------
    // This list is I in Gusfield's description.  An interval point is
    // a triple of (dimension 0 border value, is start point, pointer
    // to seed it belongs to).
    typedef Triple<TPosition, bool, TSeed *> TIntervalPoint;
    typedef String<TIntervalPoint> TIntervalPoints;
    typedef typename Iterator<TIntervalPoints, Standard>::Type TIntervalPointsIterator;

    TIntervalPoints intervalPoints;
    for (THighQualitySeedsIterator it = seedSet._highQualitySeeds.begin(), itEnd = seedSet._highQualitySeeds.end(); it != itEnd; ++it) {
        appendValue(intervalPoints, TIntervalPoint(getBeginDim0(**it), true, *it));
        appendValue(intervalPoints, TIntervalPoint(getEndDim0(**it), false, *it));
    }
    std::sort(begin(intervalPoints, Standard()), end(intervalPoints, Standard()));
    for (unsigned i = 0; i < length(intervalPoints); ++i) {
        std::cout << "IP (" << intervalPoints[i].i1 << ", " << intervalPoints[i].i2 << ", " << intervalPoints[i].i3 << ")" << std::endl;
    }

    // -----------------------------------------------------------------------
    // Step 2: Build the chain.
    // -----------------------------------------------------------------------
    // We build a list of "intermediate solutions".  Each such
    // solution is represented by the triple (end position in dim1,
    // value of best chain so far, last seed of the chain).
    typedef Triple<TPosition, TSize, TSeed *> TIntermediateSolution;
    typedef std::multiset<TIntermediateSolution> TIntermediateSolutions;
    typedef typename TIntermediateSolutions::iterator TIntermediateSolutionsIterator;

    // For all interval points...
    TIntermediateSolutions intermediateSolutions;
    for (TIntervalPointsIterator it = begin(intervalPoints), itEnd = end(intervalPoints); it != itEnd; ++it) {
        // The seed belonging ot the interval point is seed k.
        TSeed const & seedK = value(value(it).i3);

        std::cout << "Hit (" << value(it).i1 << ", " << value(it).i2 << ", " << value(it).i3 << ")" << std::endl;
        if (value(it).i2) {  // Is is begin point.
            // Find the closest seed (in dimension 1) to seed k with an
            // entry in intermediateSolutions whose end coordinate in
            // dimension 1 is <= the begin coordinate in dimension 1
            // of seedK.
            //
            // STL gives us upper_bound which returns a pointer to the
            // *first* one that compares greater than the reference
            // one.  Searching for the next and decrementing the
            // result iterator gives the desired result.
            TIntermediateSolution referenceSolution(getBeginDim1(seedK), 0, 0);
            TIntermediateSolutionsIterator itJ = intermediateSolutions.upper_bound(referenceSolution);
            if (itJ == intermediateSolutions.begin())
                continue;
            SEQAN_ASSERT_GT(intermediateSolutions.size(), 0u);  // TODO(holtgrew): Remove this assertion?
            --itJ;
            // Now, we have found such a seed j.
            TSeed const & seedJ = value(itJ->i3);
            SEQAN_ASSERT_LEQ(getEndDim1(seedJ), getBeginDim1(seedK));
            // Update the intermediate solution value.
            TIntermediateSolution sol(*itJ);
            sol.i2 += getSeedSize(seedJ);
            std::cout << "V(k) = V(k) {=" << sol.i2 << "} + " << getSeedSize(seedJ) << std::endl;
            intermediateSolutions.erase(itJ);
            intermediateSolutions.insert(sol);
        } else {  // Is end point.
            // Search for the first triple in intermediateSolutions
            // where the end coordinate in dimension 1 is >= end
            // coordinate in dimension 1 for seed k.  The corresponding
            // seed is seed j.
            //
            // We work with upper_bound here which gives us the first
            // value that is > so we have to work around this to get
            // >= again...
            SEQAN_ASSERT_GT(getEndDim1(seedK), 0u);
            TIntermediateSolution referenceSolution(getEndDim1(seedK) - 1, 0, 0);
            TIntermediateSolutionsIterator itSol = intermediateSolutions.upper_bound(referenceSolution);
            if (itSol == intermediateSolutions.end()) {
                // None found.  Insert a new triple for seed k.
                TIntermediateSolution sol(getEndDim1(seedK) - 1, getSeedSize(seedK), value(it).i3);
                std::cout << "INSERT (" << sol.i1 << ", " << sol.i2 << ", " << sol.i3 << ")" << std::endl;
                intermediateSolutions.insert(sol);
            } else {
                // Found this intermediate solution.
                SEQAN_ASSERT_GEQ(itSol->i1, getEndDim1(seedK));
                TSeed const & seedJ = value(itSol->i3);
                // Possibly start a new chain at k if the end1 is
                // before the end1 of the chain ending in j or they
                // end at the same coordinate in dim1 but k already
                // has a higher quality than the whole chaing ending
                // at j.
                if (getEndDim1(seedJ) > getEndDim1(seedK) ||
                    (getEndDim1(seedJ) == getEndDim1(seedK) && getSeedSize(seedK) > itSol->i2)) {
                    TIntermediateSolution sol(getEndDim1(seedK), getSeedSize(seedK), value(it).i3);
                    std::cout << "INSERT (" << sol.i1 << ", " << sol.i2 << ", " << sol.i3 << ")" << std::endl;
                    intermediateSolutions.insert(sol);
                }
            }
            
            // Delete all intermediate solutions where end1 >= end1 of k and have a lower quality than k.
            TIntermediateSolutionsIterator itDel = intermediateSolutions.upper_bound(referenceSolution);
            TIntermediateSolutionsIterator itDelEnd = intermediateSolutions.end();
            while (itDel != itDelEnd) {
                TIntermediateSolutionsIterator ptr = itDel;
                ++itDel;
                if (ptr->i2 < V[k])
                    intermediateSolutions.erase(ptr);
            }
        }
    }

    std::cout << "Maximal quality: " << intermediateSolutions.rbegin()->i2 << std::endl;

    // -----------------------------------------------------------------------
    // Step 3: Write out the resulting chain.
    // -----------------------------------------------------------------------
    clear(target);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_SEEDS_GLOBAL_CHAINING_H_
