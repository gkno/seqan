/*
 *  index_sa_qsort.h
 *  SeqAn
 *
 *  Created by David Weese on 22.04.07.
 *
 */

#ifndef SEQAN_HEADER_INDEX_SA_QSORT_H
#define SEQAN_HEADER_INDEX_SA_QSORT_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct SAQSort {};

    template < typename TSAValue, typename TText >
	struct _SuffixLess : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TText, Standard>::Type TIter;
		TIter _begin, _end;

		_SuffixLess() {}
		_SuffixLess(TText &text): 
			_begin(begin(text, Standard())),
			_end(end(text, Standard())) {}

		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
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

    template < typename TSAValue, typename TText >
	struct _SuffixLessOffset : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TText, Standard>::Type	TIter;
		typedef typename Size<TText>::Type					TSize;

		TIter _begin, _end;

		_SuffixLessOffset() {}
		template <typename TSize>
		_SuffixLessOffset(TText &text, TSize lcp): 
			_begin(begin(text, Standard()) + lcp),
			_end(end(text, Standard())) {}

		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
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
		sort(
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
		sort(
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
