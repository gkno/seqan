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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_SA_QSORT_H
#define SEQAN_HEADER_INDEX_SA_QSORT_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct SAQSort {};

	// compare two suffices of a given text
    template < typename TSAValue, typename TText >
	struct _SuffixLess : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TText, Standard>::Type TIter;
		TIter _begin, _end;

		_SuffixLess(TText &text): 
			_begin(begin(text, Standard())),
			_end(end(text, Standard())) {}

		// skip the first <offset> characters
		template <typename TSize>
		_SuffixLess(TText &text, TSize offset): 
			_begin(begin(text, Standard()) + offset),
			_end(end(text, Standard())) {}

		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
			if (a == b) return false;
			TIter itA = _begin + a;
			TIter itB = _begin + b;
			if (a <= b) {
				for(; itB != _end; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				for(; itA != _end; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};

	// compare two suffices of a given text
    template < typename TSAValue, typename TText, typename TSetSpec >
	struct _SuffixLess<TSAValue, StringSet<TText, TSetSpec> > : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typename Size<TText>::Type _offset;
		TText const &_text;

		_SuffixLess(TText &text): 
			_text(text) {}
			
		// skip the first <offset> characters
		template <typename TSize>
		_SuffixLess(TText &text, TSize offset):
			_text(text),
			_offset(offset) {}
		
		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
			typedef typename Iterator<TText>::Type TIter;
			if (a == b) return false;
			TIter itA = begin(getValue(_text, getSeqNo(a)), Standard()) + getSeqOffset(a);
			TIter itB = begin(getValue(_text, getSeqNo(b)), Standard()) + getSeqOffset(b);
			TIter itAEnd = end(getValue(_text, getSeqNo(a)), Standard());
			TIter itBEnd = end(getValue(_text, getSeqNo(b)), Standard());
			if (itAEnd - itA <= itBEnd - itB) {
				for(; itA != itAEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				for(; itB != itBEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}	
	};

	// compare two suffices of a given text and skip the first <lcp> characters
    template < typename TSAValue, typename TText >
	struct _SuffixLessOffset: _SuffixLess<TSAValue, TText> 
	{
		// skip the first <offset> characters
		template <typename TSize>
		_SuffixLessOffset(TText &text, TSize offset): 
			_SuffixLess<TSAValue, TText> (text, offset) {}
	};

		
		
	template < typename TSA, 
			typename TText,
			typename TSize>
	void _sortBucketQuickSort(
		TSA &sa,
		TText &text,
		TSize lcp)
	{
	SEQAN_CHECKPOINT
		// sort bucket with quicksort
		::std::sort(
			begin(sa, Standard()), 
			end(sa, Standard()), 
			_SuffixLessOffset<typename Value<TSA>::Type, TText>(text, lcp));
	}

	template < typename TSA,
               typename TText >
    inline void createSuffixArray(
		TSA &SA,
		TText &s,
		SAQSort const &)
	{
	SEQAN_CHECKPOINT
		typedef typename Size<TSA>::Type TSize;
		typedef typename Iterator<TSA, Standard>::Type TIter;

		// 1. Fill suffix array with a permutation (the identity)
		TIter it = begin(SA, Standard());
		TIter itEnd = end(SA, Standard());
		for(TSize i = 0; it != itEnd; ++it, ++i)
			*it = i;

		// 2. Sort suffix array with quicksort
		::std::sort(
			begin(SA, Standard()), 
			end(SA, Standard()), 
			_SuffixLess<typename Value<TSA>::Type, TText const>(s));
	}

	template < typename TSA,
               typename TText,
			   typename TSSetSpec >
    inline void createSuffixArray(
		TSA &SA,
		StringSet< TText, TSSetSpec > const &s,
		SAQSort const &)
	{
	SEQAN_CHECKPOINT
		typedef typename Size<TSA>::Type TSize;
		typedef typename Iterator<TSA, Standard>::Type TIter;

		// 1. Fill suffix array with a permutation (the identity)
		TIter it = begin(SA, Standard());
		TIter itEnd = end(SA, Standard());
		for(unsigned j = 0; j < length(s); ++j)
		{
			TSize len = length(s[j]);
			for(TSize i = 0; i < len; ++i)
				*it = Pair<unsigned, TSize>(j, i);
		}

		// 2. Sort suffix array with quicksort
		::std::sort(
			begin(SA, Standard()), 
			end(SA, Standard()), 
			_SuffixLess<typename Value<TSA>::Type, TText const>(s));
	}

    //////////////////////////////////////////////////////////////////////////////
    // suffix quicksort pipe
    template < typename TInput >
    struct Pipe< TInput, SAQSort >
    {
		typedef typename Value<TInput>::Type	TValue;
		typedef typename SAValue<TInput>::Type	TSAValue;

		typedef String<TValue, Alloc<> >		TText;
		typedef String<TSAValue, Alloc<> >		TSA;
		typedef Pipe<TSA, Source<> >			TSource;

		TSA		sa;
		TSource	in;

		Pipe(TInput &_textIn):
			in(sa)
		{
			TText text;
			text << _textIn;

			resize(sa, length(_textIn), Exact());
			createSuffixArray(sa, text, SAQSort());
		}

		inline typename Value<TSource>::Type const & operator*() {
            return *in;
        }
        
        inline Pipe& operator++() {
            ++in;
            return *this;
        }        
	};

}

#endif
