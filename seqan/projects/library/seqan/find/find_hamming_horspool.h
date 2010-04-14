/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  This header defines a Finder specialization that implements the Horspool
  algorithm for approximate string searching with the Hamming distance as
  described in:

    J. Tarhio and E. Ukkonen.  Boyer-moore approach to approximate
    string matching. SWAT '90, pages 348–359, 1990.
 ============================================================================
  Status: The code compiles but does NOT WORK CORRECTLY.
 ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ==========================================================================*/

#ifndef SEQAN_FIND_HAMMING_HORSPOOL_H_
#define SEQAN_FIND_HAMMING_HORSPOOL_H_

namespace SEQAN_NAMESPACE_MAIN {

/**
.Spec.HammingHorspool:
..summary:Hamming distance string matching using approximate Boyer-Moore-Horspool algorithm
..general:Class.Pattern
..cat:Pattern Matching
..signature:Pattern<TNeedle, HammingHorspool>
..param.TNeedle:The needle type.
...type:Class.String
..include:seqan/find.h
*/

///.Class.Pattern.param.TSpec.type:Spec.HammingHorspool

struct _HammingHorspool;
typedef Tag<_HammingHorspool> HammingHorspool;


template <typename TNeedle>
class Pattern<TNeedle, HammingHorspool> {
public:
    typedef typename Size<TNeedle>::Type TSize;

    // Holder of the pattern's needle.
    Holder<TNeedle> data_host;

    // Shift table d_k[i,a].
    String<TSize> shift_table;

    // The current end position
    TSize _currentEndPos;

    // The shift to apply in the next iteration;
    TSize _nextShift;

    // Maximal number of allowed mismatches.
    unsigned int k;

    // The distance at the current position, i.e. -score.
    int distance;

    Pattern() {
        SEQAN_CHECKPOINT;
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl, int k = -1) {
        SEQAN_CHECKPOINT;
        setHost(*this, ndl, k);
    }
};


template <typename TNeedle>
inline int score(const Pattern<TNeedle, HammingHorspool> &me) {
    SEQAN_CHECKPOINT;
    return -me.distance;
}


template <typename TNeedle>
inline int getScore(const Pattern<TNeedle, HammingHorspool> &me) {
    SEQAN_CHECKPOINT;
    return -me.distance;
}



template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, HammingHorspool> & me,
        TNeedle2 const & needle,
        int const k_) {
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_NOT(empty(needle));
    SEQAN_ASSERT_LEQ_MSG(k_, 0, "Are you confusing distances and scores?");

    typedef typename Position<TNeedle>::Type TPosition;
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TValue;
    typedef typename Iterator<TNeedle, Standard>::Type TNeedleIterator;

    // Initialize pattern.
    me.k = -k_;

    // Perform precomputation.  The notation roughly follows Algorithm
    // 3 from Tahio and Ukkonen, 1990.  However, me.k corresponds to k
    // and me.shift_table corresponds to d_k (we add a helper macro
    // SEQAN_d_k for this).
    const TSize &m = length(needle);  // Needle length.
    const TSize &c = ValueSize<TValue>::VALUE;  // Alphabet size.

    // Define and initialize the temporary "ready" table.  Note that
    // seqan::fill() also resizes containers.  This corresponds to lines
    // 2-4, but we fill me.shift_table a bit more than there.
    String<TSize> ready;
    fill(ready, c, m + 1);
    // Initialize the shift table (== d_k[i,a]).
    fill(me.shift_table, c * (me.k + 1), m);

    // Lines 5-8.  However, we use 0-based arrays.
    for (TPosition i = m - 1; i >= 1; --i) {
        const TPosition v = _max(i + 1, m - me.k);  // Shortcut.
        const TPosition p_i = ordValue(needle[i - 1]);
        for (TPosition j = ready[p_i] - 1; j >= v; --j) {
            int idx = p_i * (me.k + 1) + j - (m - me.k);
            me.shift_table[idx] = j - i;
        }
        ready[p_i] = v;
    }
    // Undefine helper macro.
    #undef SEQAN_d_k
    // Done with precomputation.

#if SEQAN_ENABLE_DEBUG
    // Validate precomputation.  All shifts should be > 0 and <= needle length.
    for (TPosition i = 0; i <= me.k; ++i) {
        std::cout << "i = " << i + (m - me.k) << " -- ";
        for (TPosition a = 0; a < c; ++a) {
            SEQAN_ASSERT_GT(me.shift_table[a * (me.k + 1) + i], 0u);
            SEQAN_ASSERT_LEQ(me.shift_table[a * (me.k + 1) + i], m);
            std::cout << me.shift_table[a * (me.k + 1) + i] << " ";
        }
        std::cout << std::endl;
    }
#endif  // SEQAN_DEBUG

    // In the beginning the first tentative end position is m.
    me._currentEndPos = m;
    // The shift is 0 then.
    me._nextShift = 0;

    // Actually assign the needle as the host.
    me.data_host = needle;
}


template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, HammingHorspool> & horsp, TNeedle2 & ndl, unsigned int k) {
    SEQAN_CHECKPOINT;
    setHost(horsp, reinterpret_cast<TNeedle2 const &>(ndl),k);
}


template <typename TNeedle>
inline void _finderInit(Pattern<TNeedle, HammingHorspool> & me) {
    SEQAN_CHECKPOINT;
    (void) me;  // Suppress unused variable warning.
}


template <typename TNeedle>
inline typename Host<Pattern<TNeedle, HammingHorspool> >::Type &
host(Pattern<TNeedle, HammingHorspool> & me)
{
    SEQAN_CHECKPOINT;
    return value(me.data_host);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, HammingHorspool> >::Type &
host(Pattern<TNeedle, HammingHorspool> const & me)
{
    SEQAN_CHECKPOINT;
    return value(me.data_host);
}


template <typename TFinder, typename TNeedle>
bool find(TFinder &finder, Pattern<TNeedle, HammingHorspool> &me) {
    SEQAN_CHECKPOINT;

    typedef typename Position<TFinder>::Type TPos;
    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Size<TNeedle>::Type TSize;

    THaystack &haystack = container(finder);
    TNeedle &needle = host(me);

    //std::cout << "haystack = " << haystack << std::endl;
    //std::cout << "needle = " << needle << std::endl;

    // Guard: Cannot find anything if needle longer than haystack.
    if (length(me) > length(haystack))
        return false;
    // Guard: If we are already at the end then return false.
    if (me._currentEndPos + me._nextShift > length(haystack))
        return false;

    // If the finder is empty then initialize it.
    if (empty(finder)) {
        _finderInit(me);
        _setFinderLength(finder, length(needle));
        _finderSetNonEmpty(finder);
    }

    // Actually run the finder.  The notation follows Algorithm 4 from
    // Tahio and Ukkonen, 1990.
    TPos n = length(haystack);
    TPos m = length(needle);
    // Helper macro d_k allows us to access me.shift_table the way it is
    // accessed in the paper.
    //const TSize offset = m - me.k;  // Offset for d_k[i,a]
    //std::cout << "m = " << m << ", me.k = " << me.k << std::endl;
    #define SEQAN_d_k(i, a) (me.shift_table[(a) * (me.k + 1) + (i) - m + me.k + 1])

    // Match with pattern ends at position j.
    TSize &j = me._currentEndPos;
    TPos &d = me._nextShift;
    while (j <= n) {
        // Algorithm 4, line 11 (need to emulate coroutines).
        j += d;  // Shift to the right.
        //std::cout << "beginning of loop" << std::endl;
        // Algorithm 4, lines 4-5.
        // h scans the text, i scans the text, both are indices.
        TPos h = j - 1;
        TPos i = m - 1;
        // neq is the number of mismatches.
        size_t neq = 0;
        // d is the current shift, set its initial value.
        d = m - me.k;
        // Algorithm 4, lines 6-9.
        // Compare needle and haystack from right to left, starting at h and i.
        // This could read "while (i >= 0 && neq <= me.k)" if i was signed.
        while (true) {  // We will break out of loop below.
            unsigned t_h = ordValue(haystack[h]);
            unsigned p_i = ordValue(needle[i]);
            if (i + 1 >= m - me.k) { // Minimize d over component shifts.
                d = _min(d, SEQAN_d_k(i, t_h));
            }
            //SEQAN_ASSERT_GT(d, 0);  // We have to make progress!
            //std::cout << "  h = " << h << ", i == " << i << std::endl;
            if (t_h != p_i) {  // Count mismatches.
                //std::cout << "  mismatch " << t_h << " != " << p_i << std::endl;
                neq += 1;
            } else {
                //std::cout << "  hit " << t_h << " == " << p_i << std::endl;
            }
            // We break out of this loop here because i is unsigned
            // and cannot become < 0.
            if (i == 0 || neq > me.k)
                break;
            // Proceed to the left.
            i -= 1;
            h -= 1;
        }
        // Algorithm 4, line 10.
        //std::cout << "  neq = " << neq << ", me.k = " << me.k << std::endl;
        if (neq <= me.k) {
            me.distance = neq;  // Set current score.
            // Report match at position j, begins at j - m.
            _setFinderEnd(finder, j);
            setPosition(finder, j - m);
            return true;
        }
    }
    // Undefine the helper macro again.
    #undef SEQAN_d_k

    return false;
}


template <typename TNeedle>
struct Host<Pattern<TNeedle, HammingHorspool> > {
    typedef TNeedle Type;
};


template <typename TNeedle>
struct Host<const Pattern<TNeedle, HammingHorspool> > {
    typedef TNeedle const Type;
};

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_FIND_HAMMING_HORSPOOL_H_
