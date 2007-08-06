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

#ifndef SEQAN_HEADER_RADIX_H
#define SEQAN_HEADER_RADIX_H

namespace SEQAN_NAMESPACE_MAIN
{

    // map alphabet to unsigned ints
    // signed alphabet needs a function to map char to [0..n) not to [-128..128) like char does

	template <typename TValue>
	inline unsigned _ord(TValue const &c) {
		return (typename _MakeUnsigned<TValue>::Type const &)c;
	}

	template <typename TValue, typename TSpec>
	inline unsigned _ord(SimpleType<TValue,TSpec> const &c) {
		return c;
	}


	// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
	template <
		typename TSortedArray, 
		typename TUnsortedArray, 
		typename TCountArray, 
		typename TText >
	void radixPass(
		TSortedArray &b,			// result array (sorted by 1 character)
		TUnsortedArray const &a,	// source array (unsorted)
		TText const &r,				// text to compare with
		TCountArray &c,				// temp. counter array
		unsigned K)					// alphabet size
	{ // count occurrences
		typedef typename Value<TCountArray>::Type	TSize;
		typedef typename Value<TText>::Type			TValue;

		TSize i, sum = 0, n = length(a);								// for (i = 0;  i < K;  i++) c[i] = 0;
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		for (i = 0;  i < n;  i++)										// count occurences
			c[_ord(r[a[i]])]++;
		for (i = 0;  i < K;  i++) {										// exclusive prefix sums
			TSize t = c[i];  c[i] = sum;  sum += t;
		}
		for (i = 0;  i < n;  i++) {
			TSize j = a[i];												// sort
			b[c[_ord(r[j])]++] = j;
		}
	}

	template <
		typename TSortedArray, 
		typename TUnsortedArray, 
		typename TCountArray, 
		typename TText, 
		typename TShift >
	void radixPass(
		TSortedArray &b,			// result array (sorted by 1 character)
		TUnsortedArray const &a,	// source array (unsorted)
		TText const &r,				// text to compare with
		TCountArray &c,				// temp. counter array
		unsigned K,					// alphabet size
		TShift shift)				// shift value
	{ // count occurrences
		typedef typename Value<TCountArray>::Type	TSize;
		typedef typename Value<TText>::Type			TValue;

		TSize i, sum = 0, n = length(a), sn = length(r);
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		for (i = 0;  i < n;  i++) {
			TSize j = a[i] + shift;										// count occurences
			if (j < sn) c[_ord(r[j])]++;
			else        sum++;
		}
		for (i = 0;  i < K;  i++) {										// exclusive prefix sums
			TSize t = c[i];  c[i] = sum;  sum += t;
		}
		for (i = 0, sum = 0;  i < n;  i++) {							// sort
			TSize j = a[i] + shift;
			if (j < sn) b[c[_ord(r[j])]++] = a[i];	// On Exception: Make sure, you have resized your sufarray
			else        b[sum++    ] = a[i];		// before calling createSuffixArray(..) to length(text)?
		}
	}

    // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
	template <
		typename TSortedArray, 
		typename TUnsortedArray, 
		typename TCountArray, 
		typename TText >
	void radixExtend(
		TSortedArray &b,			// result array (sorted by 1 character)
		TUnsortedArray const &a,	// source array (unsorted)
		TText const &r,				// text to compare with
		TCountArray &c,				// temp. counter array
		unsigned K)					// alphabet size
	{ // count occurrences
		typedef typename Value<TCountArray>::Type	TSize;
		typedef typename Value<TText>::Type			TValue;

		const TValue *rp = begin(r, Standard()) - 1;
		TSize i, sum = 0, n = length(a);								// for (i = 0;  i < K;  i++) c[i] = 0;
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		for (i = 0;  i < n;  i++)										// count occurences
			c[_ord(rp[a[i]])]++;
		for (i = 0;  i < K;  i++) {										// exclusive prefix sums
			TSize t = c[i];  c[i] = sum;  sum += t;
		}
		for (i = 0;  i < n;  i++) {
			TSize j = a[i];												// sort
			b[c[_ord(rp[j])]++] = j - 1;
		}
	}

    // stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 from r
	template <
		typename TSortedArray, 
		typename TUnsortedArray, 
		typename TCountArray, 
		typename TText >
	void radixExtendClip(
		TSortedArray &b,			// result array (sorted by 1 character)
		TUnsortedArray const &a,	// source array (unsorted)
		TText const &r,				// text to compare with
		TCountArray &c,				// temp. counter array
		unsigned K)					// alphabet size
	{ // count occurrences
		typedef typename Value<TCountArray>::Type	TSize;
		typedef typename Value<TText>::Type			TValue;

		const TValue *rp = begin(r, Standard()) - 1;
        TSize i, sum = 0, n = length(a);								// for (i = 0;  i < K;  i++) c[i] = 0;
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		for (i = 0;  i < n;  i++) {										// count occurences
	        TSize j = a[i];
			if (j > 0) c[_ord(rp[j])]++;
		}
		for (i = 0;  i < K;  i++) {										// exclusive prefix sums
			TSize t = c[i];  c[i] = sum;  sum += t;
		}
		for (i = 0;  i < n;  i++) {
			TSize j = a[i];												// sort
			if (j > 0) b[c[_ord(rp[j])]++] = j - 1;
		}
	}

}

#endif
