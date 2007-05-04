/*
 *  radix.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_RADIX_H
#define SEQAN_HEADER_RADIX_H

namespace SEQAN_NAMESPACE_MAIN
{

    // map alphabets to unsigned ints
    // signed alphabet needs a function to map char to [0..n) not to [-128..128) like char does

    template < typename TValue >
    struct RadixMap: public ::std::unary_function<TValue, unsigned> {
        RadixMap(unsigned) {}
        inline TValue const & operator() (TValue const &x) const { return x; }
    };

    template <>
    struct RadixMap<char>: public ::std::unary_function<char, unsigned> {
        RadixMap(unsigned) {}
        inline unsigned char const & operator() (char const & x) const {
            return reinterpret_cast<unsigned char const &>(x); 
        }
    };

    template <typename TSpec>
    struct RadixMap<SimpleType<char,TSpec> >: public ::std::unary_function<SimpleType<char,TSpec>, unsigned> {
        RadixMap(unsigned) {}
        inline unsigned operator() (SimpleType<char,TSpec> const & x) const {
            return (unsigned)x;
        }
    };


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

        RadixMap<TValue> map(K);
		TSize i, sum = 0, n = length(a);						// for (i = 0;  i < K;  i++) c[i] = 0;
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		for (i = 0;  i < n;  i++)								// count occurences
			c[map(r[a[i]])]++;
		for (i = 0;  i < K;  i++) {								// exclusive prefix sums
			TSize t = c[i];  c[i] = sum;  sum += t;
		}
		for (i = 0;  i < n;  i++) {
			TSize j = a[i];										// sort
			b[c[map(r[j])]++] = j;
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

        RadixMap<TValue> map(K);
		TSize i, sum = 0, n = length(a), sn = length(r);
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		for (i = 0;  i < n;  i++) {
			TSize j = a[i] + shift;								// count occurences
			if (j < sn) c[map(r[j])]++;
			else        sum++;
		}
		for (i = 0;  i < K;  i++) {                             // exclusive prefix sums
			TSize t = c[i];  c[i] = sum;  sum += t;
		}
		for (i = 0, sum = 0;  i < n;  i++) {                    // sort
			TSize j = a[i] + shift;
			if (j < sn) b[c[map(r[j])]++] = a[i];	// On Exception: Make sure, you have resized your sufarray
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

        RadixMap<TValue> map(K);
		const TValue *rp = begin(r, Standard()) - 1;
		TSize i, sum = 0, n = length(a);						// for (i = 0;  i < K;  i++) c[i] = 0;
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		for (i = 0;  i < n;  i++)								// count occurences
			c[map(rp[a[i]])]++;
		for (i = 0;  i < K;  i++) {								// exclusive prefix sums
			TSize t = c[i];  c[i] = sum;  sum += t;
		}
		for (i = 0;  i < n;  i++) {
			TSize j = a[i];										// sort
			b[c[map(rp[j])]++] = j - 1;
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

        RadixMap<TValue> map(K);
		const TValue *rp = begin(r, Standard()) - 1;
        TSize i, sum = 0, n = length(a);						// for (i = 0;  i < K;  i++) c[i] = 0;
		arrayFill(begin(c, Standard()), begin(c, Standard()) + K, 0);	// reset counters
		for (i = 0;  i < n;  i++) {								// count occurences
	        TSize j = a[i];
			if (j > 0) c[map(rp[j])]++;
		}
		for (i = 0;  i < K;  i++) {								// exclusive prefix sums
			TSize t = c[i];  c[i] = sum;  sum += t;
		}
		for (i = 0;  i < n;  i++) {
			TSize j = a[i];										// sort
			if (j > 0) b[c[map(rp[j])]++] = j - 1;
		}
	}

}

#endif
