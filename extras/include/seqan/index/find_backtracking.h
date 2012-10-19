// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Approximate string matching via backtracking on VSTrees
// ==========================================================================

#ifndef SANDBOX_ESIRAGUSA_INCLUDE_SEQAN_FIND_BACKTRACKING_H_
#define SANDBOX_ESIRAGUSA_INCLUDE_SEQAN_FIND_BACKTRACKING_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TDistance = HammingDistance, typename TSpec = void>
struct Backtracking {};

// ============================================================================

template <typename TPrefix, typename TDistance>
struct State {};

// ============================================================================

template <typename TSuffix, typename TDistance>
class SuffixAligner
{};

template <typename TSuffix>
class SuffixAligner<TSuffix, HammingDistance>
{
protected:
    typedef typename Size<TSuffix>::Type    TSize;

public:
//    Holder<TSuffix>     suffix;
    TSize               suffix_length;
    TSize               position;

    SuffixAligner() {}

    SuffixAligner(TSuffix & suffix) :
//        suffix(suffix),
        position(0)
    {}
};

// ============================================================================

template <typename TPrefix, typename TDistance>
class PrefixAligner
{};

template <typename TPrefix>
class PrefixAligner<TPrefix, HammingDistance>
{
protected:
    typedef typename Size<TPrefix>::Type    TSize;

public:
//    Holder<TPrefix>     prefix;
    TSize               prefix_length;
    TSize               position;
    unsigned            errors;

    PrefixAligner() {}

    PrefixAligner(TPrefix const & prefix) :
//        prefix(prefix),
        position(0)
    {}
};

// ============================================================================

/**
 .Spec.Backtracking:
 ..summary:Provides approximate string matching via backtracking on a substring index.
 ..general:Class.Pattern
 ..general:Class.Finder
 ..cat:Searching
 ..signature:Finder<TIndex, Backtracking<TDistance> >
 ..signature:Pattern<TNeedle, Backtracking<TDistance> >
 ..param.TIndex: An index of the sequence that should be searched.
 ...type:Spec.Index
 ..param.TNeedle: The type of the sequence(s) that should be searched for.
 ...type:Spec.String
 ...type:Spec.Index
 ..param.TDistance: Specifies the distance filter.
 ...type:Spec.HammingDistance
 ...type:Spec.EditDistance
 ..include:seqan/index.h
 ..remarks:
 The @Class.Pattern@ can be a sequence, or an index for a set of sequences. The tolerated backtracking distance must be given when @Function.find@ is called.
 */
///.Class.Pattern.param.TSpec.type:Spec.Backtracking
///.Class.Finder.param.TSpec.type:Spec.Backtracking

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec>
class Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> >
{
protected:
    typedef Index<TText, TSpec>                                         TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type                 TIndexIterator;
    typedef String<TIndexIterator, Block<> >                            TParentStack;

    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;
    typedef typename Size<TIndex>::Type                                 TSize;

    typedef typename Fibre<TIndex, EsaText>::Type                       TSAText;
    typedef typename Infix<TSAText const>::Type                         TSuffix;
    typedef SuffixAligner<TSuffix, TDistance>                           TSuffixAligner;

public:
    bool                index_iterator_at_root;

    Holder<TIndex>      index;
    TIndexIterator      index_iterator;
    TParentStack        index_parents;
    Pair<TSize>         index_range;

    TIterator           data_iterator;
    TSize               data_length;
    Pair<TIterator>     range;

    TSuffixAligner      suffix_aligner;

    Finder() {}

    Finder(TIndex & index) :
        index(index),
        index_iterator(index)
    {
        clear(*this);
    }

    Finder(TIndex const & index) :
        index(index),
        index_iterator(index)
    {
        clear(*this);
    }

    ~Finder()
    {
#ifdef SEQAN_DEBUG
        std::cout << "Finder Parents Height: " << length(this->index_parents) << std::endl;
#endif
    }

};

// ============================================================================

template <typename TNeedle, typename TDistance, typename TBacktrackingSpec>
class Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> >
{
protected:
    typedef TNeedle                                     TPrefix;
    typedef PrefixAligner<TPrefix, TDistance>           TPrefixAligner;

    typedef typename State<TPrefix, TDistance>::Type    TState;
    typedef String<TState, Block<> >                    TStateStack;

public:
    Holder<TNeedle>         data_host;
    TPrefixAligner          prefix_aligner;
    TStateStack             state;

    Pattern() {}

    Pattern(TNeedle const & needle)
    {
        setHost(*this, needle);
    }

    ~Pattern()
    {
        // Empty stack
        clear(this->state);
    }

};

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
class Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> >
{
protected:
    typedef Index<TNeedle, TSpec>                                       TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type                 TIndexIterator;
    typedef String<TIndexIterator, Block<> >                            TParentStack;

    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;
    typedef typename Size<TIndex>::Type                                 TSize;

    typedef typename Fibre<TIndex, EsaText>::Type                       TSAText;
    typedef typename Infix<TSAText const>::Type                         TPrefix;
    typedef PrefixAligner<TPrefix, TDistance>                           TPrefixAligner;

    typedef typename State<TPrefix, TDistance>::Type                    TState;
    typedef String<TState, Block<> >                                    TStateStack;

    typedef String<bool, Block<> >                                      TEndStack;

public:
    bool                index_iterator_at_root;

    Holder<TIndex>      data_host;
    TIndexIterator      index_iterator;
    TParentStack        index_parents;
    Pair<TSize>         index_range;

    TIterator           data_iterator;
    TSize               data_length;
    Pair<TIterator>     range;

    TPrefixAligner      prefix_aligner;
    TStateStack         state;
    TEndStack           atEnd;

    unsigned            exact;
    bool                search;

    // TODO(esiragusa): Remove depth, isLeaf(it) should return true when repLength(it) >= depth
    TSize               depth;

    Pattern() {}

    Pattern(TIndex & index, TSize depth) :
        data_host(index),
        index_iterator(index),
        depth(depth)
    {
        setHost(*this, index);
    }

    Pattern(TIndex const & index, TSize depth) :
        data_host(index),
        index_iterator(index),
        depth(depth)
    {
        setHost(*this, index);
    }

    ~Pattern()
    {
#ifdef SEQAN_DEBUG
        std::cout << "State Height: " << length(this->state) << std::endl;
        std::cout << "atEnd Height: " << length(this->atEnd) << std::endl;
        std::cout << "Pattern Parents Height: " << length(this->index_parents) << std::endl;
#endif

        // Empty stacks
        clear(this->index_parents);
        clear(this->state);
        clear(this->atEnd);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TPrefix>
struct State<TPrefix, HammingDistance>
{
    typedef typename Size<TPrefix>::Type            TSize;
    typedef Pair<TSize, unsigned>                   Type;
};

// ============================================================================

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
struct Position<Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > >:
    SAValue<Index<TNeedle, TSpec> >
{};

// ============================================================================
// Functions
// ============================================================================

template <typename TPrefix>
inline unsigned
getScore(PrefixAligner<TPrefix, HammingDistance> & me, bool isLeaf)
{
    return (atEnd(me) && isLeaf) ? me.errors : std::numeric_limits<unsigned>::max();
}

template <typename TPrefix>
inline unsigned
getMinScore(PrefixAligner<TPrefix, HammingDistance> & me)
{
    return me.errors;
}

template <typename TPrefix>
inline typename State<TPrefix, HammingDistance>::Type
getInitialState(PrefixAligner<TPrefix, HammingDistance> &)
{
    typedef typename State<TPrefix, HammingDistance>::Type  TState;
    return TState(0, 0);
}

template <typename TPrefix>
inline typename State<TPrefix, HammingDistance>::Type
getState(PrefixAligner<TPrefix, HammingDistance> & me)
{
    typedef typename State<TPrefix, HammingDistance>::Type  TState;
    return TState(me.position, me.errors);
}

template <typename TPrefix, typename TState>
inline void
setState(PrefixAligner<TPrefix, HammingDistance> & me, TState const & state)
{
    me.position = state.i1;
    me.errors = state.i2;
#ifdef SEQAN_DEBUG
    std::cout << "Old State:      " << "(" << state.i1 << ", " << state.i2 << ")" << std::endl;
#endif
}

template <typename TSuffix, typename TState>
inline void
setState(SuffixAligner<TSuffix, HammingDistance> & me, TState const & state)
{
    me.position = state.i1;
}

template <typename TSuffix, typename TSize>
inline void setSuffixLength(SuffixAligner<TSuffix, HammingDistance> & me, TSize suffix_length)
{
    me.suffix_length = suffix_length;
}

template <typename TPrefix, typename TSize>
inline void setPrefixLength(PrefixAligner<TPrefix, HammingDistance> & me, TSize prefix_length)
{
    me.prefix_length = prefix_length;
}

template <typename TSuffix, typename TPrefix, typename TErrors>
inline bool align(SuffixAligner<TSuffix, HammingDistance> & suffix_aligner,
                  PrefixAligner<TPrefix, HammingDistance> & prefix_aligner,
                  TSuffix & suffix,
                  TPrefix & prefix,
                  TErrors errors)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize min_length = std::min(suffix_aligner.suffix_length, prefix_aligner.prefix_length);

    // TODO(esiragusa): Use iterators instead of subscripts
    while ((prefix_aligner.position) < min_length)
    {
        if (suffix[suffix_aligner.position] != prefix[prefix_aligner.position])
            if (++prefix_aligner.errors > errors)
                return false;

        suffix_aligner.position++;
        prefix_aligner.position++;
    }

    return true;
}

template <typename TSuffix, typename TPrefix>
inline bool match(SuffixAligner<TSuffix, HammingDistance> & suffix_aligner,
                  PrefixAligner<TPrefix, HammingDistance> & prefix_aligner,
                  TSuffix & suffix,
                  TPrefix & prefix)
{
    typedef typename Size<TPrefix>::Type            TSize;

    TSize min_length = std::min(suffix_aligner.suffix_length, prefix_aligner.prefix_length);

    // TODO(esiragusa): Use iterators instead of subscripts
    while ((prefix_aligner.position) < min_length)
    {
        if (suffix[suffix_aligner.position] != prefix[prefix_aligner.position])
            return false;

        suffix_aligner.position++;
        prefix_aligner.position++;
    }

    return true;
}

template <typename TPrefix>
inline bool atEnd(PrefixAligner<TPrefix, HammingDistance> & me)
{
    return (me.position == me.prefix_length) ? true : false;
}

template <typename TSuffix>
inline bool atEnd(SuffixAligner<TSuffix, HammingDistance> & me)
{
    return (me.position == me.suffix_length) ? true : false;
}

// ============================================================================

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline void
clear(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    typedef Index<TText, TSpec>                                         TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type                 TIndexIterator;

    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;

    typedef typename Fibre<TIndex, EsaText>::Type                       TSAText;
    typedef typename Infix<TSAText const>::Type                         TSuffix;
    typedef SuffixAligner<TSuffix, TDistance>                           TSuffixAligner;

    // Init backtracking on root node
    me.index_iterator = TIndexIterator(host(me));
    me.index_iterator_at_root = true;

    // Empty parent stack
    clear(me.index_parents);

    // Init data iterator on empty range
    hostIterator(me) = begin(indexSA(host(me)), Standard());
    me.range.i1 = me.range.i2 = TIterator();
    me.data_length = 0;

    // Call SuffixAligner constructor
    me.suffix_aligner = TSuffixAligner();
}

// ============================================================================

template <typename TNeedle, typename TDistance, typename TBacktrackingSpec, typename TNewNeedle>
void setHost(Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & me, TNewNeedle const & needle)
{
    typedef TNeedle                                     TPrefix;
    typedef PrefixAligner<TPrefix, TDistance>           TPrefixAligner;
    typedef typename State<TPrefix, TDistance>::Type    TState;

    SEQAN_ASSERT_NOT(empty(needle));
    setValue(me.data_host, needle);

    // Call PrefixAligner constructor
    me.prefix_aligner = TPrefixAligner();

    // Init stack on empty word
    clear(me.state);
    appendValue(me.state, getInitialState(me.prefix_aligner));
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
void _moveIteratorAtRoot(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    typedef Index<TNeedle, TSpec>                                       TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type                 TIndexIterator;
    me.index_iterator = TIndexIterator(host(me));
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNewNeedle, typename TNewSpec>
void setHost(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me, Index<TNewNeedle, TNewSpec> const & index)
{
    typedef Index<TNeedle, TSpec>                                       TIndex;
//		typedef typename Iterator< TIndex, TopDown<> >::Type                TIndexIterator;
    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;

    typedef typename Fibre<TIndex, EsaText>::Type                       TSAText;
    typedef typename Infix<TSAText const>::Type                         TPrefix;
    typedef PrefixAligner<TPrefix, TDistance>                           TPrefixAligner;
    typedef typename State<TPrefix, TDistance>::Type                    TState;

    // TODO(esiragusa): Update index holder
    setValue(me.data_host, index);

    // Init backtracking on root node
    _moveIteratorAtRoot(me);
//		me.index_iterator = TIndexIterator(host(me));
    me.index_iterator_at_root = true;

    // Empty parent stack
    clear(me.index_parents);

    // Init data iterator on empty range
    hostIterator(me) = begin(indexSA(host(me)), Standard());
    me.range.i1 = me.range.i2 = TIterator();
    me.data_length = 0;

    // Init stack on empty word
    clear(me.state);
    appendValue(me.state, getInitialState(me.prefix_aligner));

    // Empty atEnd stack
    clear(me.atEnd);
//		appendValue(me.atEnd, true);

    // Empty exact search stack
    me.exact = 0;
    me.search = false;

    // Call PrefixAligner constructor
    me.prefix_aligner = TPrefixAligner();
}

template <typename TNeedle, typename TDistance, typename TBacktrackingSpec, typename TNewNeedle>
void setHost(Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & me, TNewNeedle & needle)
{
    setHost(me, reinterpret_cast<TNeedle const &>(needle));
}

// ============================================================================

template < typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec >
inline typename Iterator< typename Fibre<Index<TNeedle, TSpec>, FibreSA>::Type const, Standard>::Type const &
hostIterator(Pattern< Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > const & me)
{
    return me.data_iterator;
}

template < typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec >
inline typename Iterator< typename Fibre<Index<TNeedle, TSpec>, FibreSA>::Type const, Standard >::Type &
hostIterator(Pattern< Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return me.data_iterator;
}

// ============================================================================

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline bool
empty(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return me.range.i1 == me.range.i2;
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline bool
atBegin(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return (empty(me) || hostIterator(me) == me.range.i1);
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline bool
atEnd(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return (empty(me) || hostIterator(me) == me.range.i2);
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline void
goBegin(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    hostIterator(me) = me.range.i1;
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline void
goEnd(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    hostIterator(me) = me.range.i2;
}

// ============================================================================

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline typename Position<Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > >::Type
beginPosition(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    SEQAN_ASSERT_NOT(empty(me));
    return *me.data_iterator;
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline typename Position<Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > >::Type
endPosition(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return posAdd(beginPosition(me), me.data_length);
}

template <typename TNeedle, typename TSpec, typename TDistance, typename TBacktrackingSpec>
inline typename Position<Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > >::Type
position(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & me)
{
    return beginPosition(me);
}

// ============================================================================

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle>
inline bool
_resume(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
        Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    // Resume positively from root only the first time
    if (finder.index_iterator_at_root)
    {
        finder.index_iterator_at_root = false;
        return true;
    }

    return _cut(finder, pattern);
}

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TErrors>
inline bool
_backtrack(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
           Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
           TErrors errors)
{
    typedef Index<TText, TSpec>                             TIndex;
    typedef typename Iterator<TIndex, TopDown<> >::Type     TIndexIterator;
    typedef typename Fibre<TIndex, EsaText>::Type           TSA;
    typedef typename Infix<TSA const>::Type                 TSuffix;
    typedef typename Size<TIndex>::Type                     TSuffixSize;

    typedef typename State<TNeedle, TDistance>::Type        TState;
    
    setPrefixLength(pattern.prefix_aligner, length(needle(pattern)));

    do
    {
        SEQAN_ASSERT_NOT(empty(pattern.state));

#ifdef SEQAN_DEBUG
        std::cout << "Stack Height:   " << length(pattern.state) << std::endl;
        std::cout << "Suffix:         " <<
        prefix(representative(finder.index_iterator),
               std::min(repLength(finder.index_iterator), length(needle(pattern))))
        << std::endl;
#endif
        // Restore last state
        setState(finder.suffix_aligner, back(pattern.state));
        setState(pattern.prefix_aligner, back(pattern.state));

        // Update current suffix
        TSuffix suffixRepr = representative(finder.index_iterator);
        TSuffixSize suffixLength = length(suffixRepr);
        setSuffixLength(finder.suffix_aligner, suffixLength);

        // Align suffix with pattern
        align(finder.suffix_aligner, pattern.prefix_aligner, suffixRepr, needle(pattern), errors);

        // A complete match was found
        if (getScore(pattern.prefix_aligner, true) <= errors)
        {
            finder.index_range = range(finder.index_iterator);
            return true;
        }
        // Reduce to exact suffix search (speedup)
        else if (getMinScore(pattern.prefix_aligner) == errors)
        {
            TIndexIterator index_iterator(finder.index_iterator);

            // A complete match was found
            if (goDown(index_iterator, suffix(needle(pattern), pattern.prefix_aligner.position)))
            {
                // Move aligners to end of pattern
                finder.suffix_aligner.position = length(needle(pattern));
                pattern.prefix_aligner.position = length(needle(pattern));

                finder.index_range = range(index_iterator);
                return true;
            }
            else
                _cut(finder, pattern);
        }
        // Walk down text index only if an alignment is still possible
        else if (getMinScore(pattern.prefix_aligner) <= errors && !isLeaf(finder.index_iterator))
        {
            appendValue(finder.index_parents, finder.index_iterator);
            goDown(finder.index_iterator);
//            appendValue(pattern.state, getState(pattern.prefix_aligner));
            appendValue(pattern.state, TState(pattern.prefix_aligner.position, pattern.prefix_aligner.errors));
        }
        // Otherwise cut branch
        else
            _cut(finder, pattern);
    }
    while (!isRoot(finder.index_iterator));

    return false;
}

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle>
inline bool
_cut(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    // Walk to next node
    if (!goRight(finder.index_iterator))
        while (!isRoot(finder.index_iterator))
        {
            SEQAN_ASSERT_NOT(empty(finder.index_parents));
            finder.index_iterator = back(finder.index_parents);
            eraseBack(finder.index_parents);

            SEQAN_ASSERT_NOT(empty(pattern.state));
            eraseBack(pattern.state);
            if (goRight(finder.index_iterator))
                break;
        }

    return !isRoot(finder.index_iterator);
}

// ============================================================================

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec, typename TErrors>
inline bool
_resume(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
        Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
        TErrors errors)
{
    // Resume positively from root only the first time
    if (finder.index_iterator_at_root)
    {
        finder.index_iterator_at_root = false;
        return _backtrack(finder, pattern, errors);
    }

    if (pattern.search)
    {
        if (_cut_exact(finder, pattern) && _search(finder, pattern))
            return true;

        SEQAN_ASSERT_NOT(pattern.exact);

        SEQAN_ASSERT_NOT(empty(finder.index_parents));
        finder.index_iterator = back(finder.index_parents);
        eraseBack(finder.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.index_parents));
        pattern.index_iterator = back(pattern.index_parents);
        eraseBack(pattern.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.state));
        eraseBack(pattern.state);

        pattern.search = false;
    }

    return _cut(finder, pattern) && _backtrack(finder, pattern, errors);
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec, typename TErrors>
inline bool
_backtrack(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
           Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
           TErrors errors)
{
    typedef Index<TText, TTextSpec>                                     TTextIndex;
    typedef typename Iterator<TTextIndex, TopDown<> >::Type             TTextIndexIterator;
    typedef typename Fibre<TTextIndex, EsaText>::Type                   TSAText;
    typedef typename Infix<TSAText const>::Type                         TSuffix;
    typedef typename Size<TTextIndex>::Type                             TSuffixSize;

    typedef Index<TNeedle, TNeedleSpec>                                 TNeedleIndex;
    typedef typename Fibre<TNeedleIndex, EsaText>::Type                 TSANeedle;
    typedef typename Infix<TSANeedle const>::Type                       TPrefix;
    typedef typename Size<TNeedleIndex>::Type                           TPrefixSize;

    typedef typename State<TNeedle, TDistance>::Type                    TState;

    SEQAN_ASSERT_NOT(pattern.exact);

    do
    {
        SEQAN_ASSERT_NOT(empty(pattern.state));

#ifdef SEQAN_DEBUG
        std::cout << "Stack Height:   " << length(pattern.state) << std::endl;
//			std::cout << "Suffix:         "	<< representative(finder.index_iterator) << std::endl;
        std::cout << "Prefix:         " << representative(pattern.index_iterator) << std::endl;
#endif

        // Lookup last state
        setState(finder.suffix_aligner, back(pattern.state));
        setState(pattern.prefix_aligner, back(pattern.state));

        // Update current suffix and prefix
        TSuffix suffixRepr = representative(finder.index_iterator);
        TPrefix prefixRepr = representative(pattern.index_iterator);
        TSuffixSize suffixLength = length(suffixRepr);
        TPrefixSize prefixLength = length(prefixRepr);
        setSuffixLength(finder.suffix_aligner, suffixLength);
        setPrefixLength(pattern.prefix_aligner, std::min(prefixLength, pattern.depth));

        // Align suffix with prefix
        align(finder.suffix_aligner, pattern.prefix_aligner, suffixRepr, prefixRepr, errors);

        // A complete match was found
//			if (getScore(pattern.prefix_aligner, isLeaf(pattern.index_iterator)) <= errors)
        if (getScore(pattern.prefix_aligner, (pattern.prefix_aligner.position >= pattern.depth)) <= errors)
        {
            finder.index_range = range(finder.index_iterator);
            pattern.index_range = range(pattern.index_iterator);
            return true;
        }
        // Reduce to exact suffix search (speedup)
        else if (getMinScore(pattern.prefix_aligner) == errors)
        {
            pattern.search = true;

            appendValue(finder.index_parents, finder.index_iterator);
            appendValue(pattern.index_parents, pattern.index_iterator);
//            appendValue(pattern.state, getState(pattern.prefix_aligner));
            appendValue(pattern.state, TState(pattern.prefix_aligner.position, pattern.prefix_aligner.errors));

            if (_search(finder, pattern))
                return true;

            SEQAN_ASSERT_NOT(pattern.exact);

            SEQAN_ASSERT_NOT(empty(finder.index_parents));
            finder.index_iterator = back(finder.index_parents);
            eraseBack(finder.index_parents);
            SEQAN_ASSERT_NOT(empty(pattern.index_parents));
            pattern.index_iterator = back(pattern.index_parents);
            eraseBack(pattern.index_parents);
            SEQAN_ASSERT_NOT(empty(pattern.state));
            eraseBack(pattern.state);

            pattern.search = false;

            _cut(finder, pattern);
        }
        else if (getMinScore(pattern.prefix_aligner) <= errors)
        {
            if (atEnd(finder.suffix_aligner))
            {
                if (!isLeaf(finder.index_iterator))
                {
                    appendValue(finder.index_parents, finder.index_iterator);
                    appendValue(pattern.index_parents, pattern.index_iterator);
//                    appendValue(pattern.state, getState(pattern.prefix_aligner));
                    appendValue(pattern.state, TState(pattern.prefix_aligner.position, pattern.prefix_aligner.errors));
                    appendValue(pattern.atEnd, false);

                    goDown(finder.index_iterator);
                }
                else
                    _cut(finder, pattern);
            }
            else if (atEnd(pattern.prefix_aligner))
            {
                if (prefixLength < pattern.depth && !isLeaf(pattern.index_iterator))
                {
                    appendValue(finder.index_parents, finder.index_iterator);
                    appendValue(pattern.index_parents, pattern.index_iterator);
//                    appendValue(pattern.state, getState(pattern.prefix_aligner));
                    appendValue(pattern.state, TState(pattern.prefix_aligner.position, pattern.prefix_aligner.errors));
                    appendValue(pattern.atEnd, true);

                    goDown(pattern.index_iterator);
                }
                else
                    _cut(finder, pattern);
            }
        }
        else
            _cut(finder, pattern);
    }
    while (!(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator)));

    return false;
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec>
inline bool
_cut(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    SEQAN_ASSERT_NOT(pattern.exact);

    if (empty(pattern.atEnd))
        return false;

    do
    {
        SEQAN_ASSERT_NOT(empty(pattern.atEnd));

        if (!back(pattern.atEnd))
        {
            if (goRight(finder.index_iterator))
            {
                SEQAN_ASSERT_NOT(empty(pattern.index_parents));
                pattern.index_iterator = back(pattern.index_parents);
                break;
            }
        }
        else
        {
            if (goRight(pattern.index_iterator))
            {
                SEQAN_ASSERT_NOT(empty(finder.index_parents));
                finder.index_iterator = back(finder.index_parents);
                break;
            }
        }

        SEQAN_ASSERT_NOT(empty(finder.index_parents));
        finder.index_iterator = back(finder.index_parents);
        eraseBack(finder.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.index_parents));
        pattern.index_iterator = back(pattern.index_parents);
        eraseBack(pattern.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.state));
        eraseBack(pattern.state);
        SEQAN_ASSERT_NOT(empty(pattern.atEnd));
        eraseBack(pattern.atEnd);
    }
    while (!(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator)));

    return !(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator));
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec>
inline unsigned
_search(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
        Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    typedef Index<TText, TTextSpec>                                     TTextIndex;
    typedef typename Iterator<TTextIndex, TopDown<> >::Type             TTextIndexIterator;
    typedef typename Fibre<TTextIndex, EsaText>::Type                   TSAText;
    typedef typename Infix<TSAText const>::Type                         TSuffix;
    typedef typename Size<TTextIndex>::Type                             TSuffixSize;

    typedef Index<TNeedle, TNeedleSpec>                                 TNeedleIndex;
    typedef typename Fibre<TNeedleIndex, EsaText>::Type                 TSANeedle;
    typedef typename Infix<TSANeedle const>::Type                       TPrefix;
    typedef typename Size<TNeedleIndex>::Type                           TPrefixSize;

    typedef typename State<TNeedle, TDistance>::Type                    TState;
    
    do
    {
        SEQAN_ASSERT_NOT(empty(pattern.state));

#ifdef SEQAN_DEBUG
        std::cout << "Stack Height:   " << length(pattern.state) << std::endl;
        std::cout << "Exact Height:   " << pattern.exact << std::endl;
//			std::cout << "Suffix:         "	<< representative(finder.index_iterator) << std::endl;
        std::cout << "Prefix:         " << representative(pattern.index_iterator) << std::endl;
#endif

        // Lookup last state
        setState(finder.suffix_aligner, back(pattern.state));
        setState(pattern.prefix_aligner, back(pattern.state));

        // Update current suffix and prefix
        TSuffix suffixRepr = representative(finder.index_iterator);
        TPrefix prefixRepr = representative(pattern.index_iterator);
        TSuffixSize suffixLength = length(suffixRepr);
        TPrefixSize prefixLength = length(prefixRepr);
        setSuffixLength(finder.suffix_aligner, suffixLength);
        setPrefixLength(pattern.prefix_aligner, std::min(prefixLength, pattern.depth));

        // Align exactly suffix with prefix
        if (match(finder.suffix_aligner, pattern.prefix_aligner, suffixRepr, prefixRepr))
        {
            if (!atEnd(pattern.prefix_aligner) && atEnd(finder.suffix_aligner))
            {
                unsigned extension = pattern.prefix_aligner.prefix_length;

                TTextIndexIterator index_iterator(finder.index_iterator);
                if (!goDown(finder.index_iterator, infix(prefixRepr, pattern.prefix_aligner.position, extension)))
                {
                    finder.index_iterator = index_iterator;
                    if (!_cut_exact(finder, pattern))
                        break;

                    continue;
                }

                finder.suffix_aligner.position = extension;
                pattern.prefix_aligner.position = extension;

                TSuffix suffixRepr_ = representative(finder.index_iterator);
                TSuffixSize suffixLength_ = length(suffixRepr_);
                setSuffixLength(finder.suffix_aligner, suffixLength_);
            }

            // A complete match was found
            if (pattern.prefix_aligner.position == pattern.depth)
            {
                finder.index_range = range(finder.index_iterator);
                pattern.index_range = range(pattern.index_iterator);
                return true;
            }
            else if (prefixLength < pattern.depth && !isLeaf(pattern.index_iterator))
            {
                ++pattern.exact;

                appendValue(finder.index_parents, finder.index_iterator);
                appendValue(pattern.index_parents, pattern.index_iterator);
//                appendValue(pattern.state, getState(pattern.prefix_aligner));
                appendValue(pattern.state, TState(pattern.prefix_aligner.position, pattern.prefix_aligner.errors));

                goDown(pattern.index_iterator);
            }
        }
        else if (!_cut_exact(finder, pattern))
            break;

    }
    while (!(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator)));

    SEQAN_ASSERT_NOT(pattern.exact);

    return false;
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec>
inline bool
_cut_exact(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
           Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    do
    {
        if (!pattern.exact)
            return false;

        if (goRight(pattern.index_iterator))
        {
            SEQAN_ASSERT_NOT(empty(finder.index_parents));
            finder.index_iterator = back(finder.index_parents);
            break;
        }

        SEQAN_ASSERT_NOT(empty(finder.index_parents));
        finder.index_iterator = back(finder.index_parents);
        eraseBack(finder.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.index_parents));
        pattern.index_iterator = back(pattern.index_parents);
        eraseBack(pattern.index_parents);
        SEQAN_ASSERT_NOT(empty(pattern.state));
        eraseBack(pattern.state);

        --pattern.exact;
    }
    while (!(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator)));

    return !(isRoot(pattern.index_iterator) && isRoot(finder.index_iterator));
}

// ============================================================================

template <typename TText, typename TSpec, typename TBacktrackingSpec>
inline void
_setFinderLength(Finder<Index<TText, TSpec>, Backtracking<HammingDistance, TBacktrackingSpec> > & finder)
{
    _setFinderLength(finder, finder.suffix_aligner.position);
}

template <typename TNeedle, typename TBacktrackingSpec>
inline void
_setPatternLength(Pattern<TNeedle, Backtracking<HammingDistance, TBacktrackingSpec> > & pattern)
{
    pattern.data_length = pattern.prefix_aligner.position;
}

// ============================================================================

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TErrors>
inline bool
find(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
     TErrors errors)
{

    // Try to get another match from current node
    if (!atEnd(finder))
        goNext(hostIterator(finder));

    // Try to get more matches by backtracking some more
    if (atEnd(finder))
    {
        // Resume from last matching node and backtrack until next match
        if (_resume(finder, pattern) && _backtrack(finder, pattern, errors))
        {
            // Set data iterator range to the interval containing matches
//				Pair<TSize> rng = range(finder.index_iterator);
            hostIterator(finder) = begin(indexSA(host(finder)), Standard());
            finder.range.i1 = hostIterator(finder) + finder.index_range.i1;
            finder.range.i2 = hostIterator(finder) + finder.index_range.i2;
            hostIterator(finder) = finder.range.i1;

            // Set match length
            _setFinderLength(finder);
//				_setPatternLength(pattern);
        }
        // No more matches
        else
        {
            hostIterator(finder) = begin(indexSA(host(finder)), Standard());
            finder.range.i1 = hostIterator(finder);
            finder.range.i2 = hostIterator(finder);
        }
    }

    return !atEnd(finder);
}

template <typename TText, typename TTextSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle, typename TNeedleSpec, typename TErrors>
inline bool
find(Finder<Index<TText, TTextSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<Index<TNeedle, TNeedleSpec>, Backtracking<TDistance, TBacktrackingSpec> > & pattern,
     TErrors errors)
{

    // Try to get another match from current text node
    if (!atEnd(finder))
        goNext(hostIterator(finder));

    // Try to get another match from current pattern node
    if (atEnd(finder))
    {
        goNext(hostIterator(pattern));

        if (!atEnd(pattern))
        {
            // Set data iterator range to the interval containing matches
            hostIterator(finder) = begin(indexSA(host(finder)), Standard());
            finder.range.i1 = hostIterator(finder) + finder.index_range.i1;
            finder.range.i2 = hostIterator(finder) + finder.index_range.i2;
            hostIterator(finder) = finder.range.i1;
        }
        // Try to get more matches by backtracking some more
        else
        {
            // Resume from last matching node and backtrack until next match
//				if (_resume(finder, pattern) && _backtrack(finder, pattern, errors))
            if (_resume(finder, pattern, errors))
            {
                // Set data iterator range to the interval containing matches
                hostIterator(finder) = begin(indexSA(host(finder)), Standard());
                finder.range.i1 = hostIterator(finder) + finder.index_range.i1;
                finder.range.i2 = hostIterator(finder) + finder.index_range.i2;
                hostIterator(finder) = finder.range.i1;

				hostIterator(pattern) = begin(indexSA(host(pattern)), Standard());
				pattern.range.i1 = hostIterator(pattern) + pattern.index_range.i1;
				pattern.range.i2 = hostIterator(pattern) + pattern.index_range.i2;
				hostIterator(pattern) = pattern.range.i1;

                // Set match length
                _setFinderLength(finder);
                _setPatternLength(pattern);
            }
            // No more matches
            else
            {
                hostIterator(finder) = begin(indexSA(host(finder)), Standard());
                finder.range.i1 = hostIterator(finder);
                finder.range.i2 = hostIterator(finder);

				hostIterator(pattern) = begin(indexSA(host(pattern)), Standard());
				pattern.range.i1 = hostIterator(pattern);
				pattern.range.i2 = hostIterator(pattern);
            }
        }

    }

    return !(atEnd(finder) && atEnd(pattern));
}

template <typename TText, typename TSpec, typename TDistance, typename TBacktrackingSpec, typename TNeedle>
inline bool
find(Finder<Index<TText, TSpec>, Backtracking<TDistance, TBacktrackingSpec> > & finder,
     Pattern<TNeedle, Backtracking<TDistance, TBacktrackingSpec> > & pattern)
{
    return find(finder, pattern, 0u);
}

}

#endif  // #ifndef SANDBOX_ESIRAGUSA_INCLUDE_SEQAN_FIND_BACKTRACKING_H_
