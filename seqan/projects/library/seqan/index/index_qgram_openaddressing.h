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

#ifndef SEQAN_HEADER_INDEX_QGRAM_OPENADRESSING_H
#define SEQAN_HEADER_INDEX_QGRAM_OPENADRESSING_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct OpenAddressing_;
	typedef Tag<OpenAddressing_> OpenAddressing;
	
	template <typename THashValue>
	struct BucketMap 
	{
		String<THashValue>	qgramHash;
		THashValue			prime;
	};

	// use the index value type as shape value type
	template < typename TObject, typename TShapeSpec >
	struct Fibre< Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >, FibreBucketMap>
	{
		typedef typename Fibre< Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >, FibreShape>::Type TShape;
		typedef typename Value<TShape>::Type	THashValue;
		typedef BucketMap<THashValue>			Type;
	};
	
/**
.Spec.OpenAddressing
..summary:An index based on an array of sorted q-grams.
..cat:Index
..general:Spec.IndexQGram
..signature:Index<TText, IndexQGram<TShapeSpec, OpenAddressing> >
..param.TText:The text type.
...type:Class.String
..param.TShapeSpec:The @Class.Shape@ specialization type.
...note:This can be either a $TSpec$ argument (e.g. $SimpleShape$) or a complete @Class.Shape@ class (e.g. Shape<Dna, SimpleShape>).
..remarks:This index uses a non-trivial hashing for mapping q-gram hash values to buckets.
This reduces the sizes of bucket directories (QGramDir, QGramCountsDir fibres) from |\Sigma|^q to min(\alpha*n,|\Sigma|^q), for a load factor \alpha>1.
..include:seqan/index.h
.Memvar.OpenAddressing#alpha
..summary:Load factor. Controls space/time-tradeoff and must be greater 1. Default value is 1.6.
..class:Spec.OpenAddressing
*/	

#ifdef PLATFORM_WINDOWS_VS
#pragma warning( push )
// Disable warning C4521 locally (multiple copy constructors).
#pragma warning( disable: 4521 )
#endif  // PLATFORM_WINDOWS_VS

	template < typename TObject, typename TShapeSpec >
	class Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >
	{
    private:
        static const double defaultAlpha;
	public:
		typedef typename Fibre<Index, QGramText>::Type			TText;
		typedef typename Fibre<Index, QGram_SA>::Type			TSA;
		typedef typename Fibre<Index, QGramDir>::Type			TDir;
		typedef typename Fibre<Index, QGramCounts>::Type		TCounts;
		typedef typename Fibre<Index, QGramCountsDir>::Type	TCountsDir;
		typedef typename Fibre<Index, QGramShape>::Type		TShape;
		typedef typename Fibre<Index, QGramBucketMap>::Type	TBucketMap;
		typedef typename Cargo<Index>::Type						TCargo;
		typedef typename Size<Index>::Type						TSize;

		Holder<TText>	text;		// underlying text
		TSA				sa;			// suffix array sorted by the first q chars
		TDir			dir;		// bucket directory
		TCounts			counts;		// counts each q-gram per sequence
		TCountsDir		countsDir;	// directory for count buckets
		TShape			shape;		// underlying shape
		TCargo			cargo;		// user-defined cargo
		TBucketMap		bucketMap;	// bucketMap table (used by open-addressing index)
		TSize			stepSize;	// store every <stepSize>'th q-gram in the index

		double			alpha;

		Index():
			stepSize(1),
			alpha(defaultAlpha) {}

		Index(Index &other):
			text(other.text),
			sa(other.sa),
			dir(other.dir),
			counts(other.counts),
			countsDir(other.countsDir),
			shape(other.shape),
			cargo(other.cargo),
			stepSize(1),
			alpha(defaultAlpha) {}

		Index(Index const &other):
			text(other.text),
			sa(other.sa),
			dir(other.dir),
			counts(other.counts),
			countsDir(other.countsDir),
			shape(other.shape),
			cargo(other.cargo),
			stepSize(1),
			alpha(defaultAlpha) {}

		template <typename _TText>
		Index(_TText &_text):
			text(_text),
			stepSize(1),
			alpha(defaultAlpha) {}

		template <typename _TText>
		Index(_TText const &_text):
			text(_text),
			stepSize(1),
			alpha(defaultAlpha) {}

		template <typename _TText, typename _TShape>
		Index(_TText &_text, _TShape const &_shape):
			text(_text),
			shape(_shape),
			stepSize(1),
			alpha(defaultAlpha) {}

		template <typename _TText, typename _TShape>
		Index(_TText const &_text, _TShape const &_shape):
			text(_text),
			shape(_shape),
			stepSize(1),
			alpha(defaultAlpha) {}
	};
#ifdef PLATFORM_WINDOWS_VS
// Enable warning C4521 again (multiple copy operators).
#pragma warning( pop )
#endif  // PLATFORM_WINDOWS_VS


    template < typename TObject, typename TShapeSpec >
    const double Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >::defaultAlpha = 1.6;

	//////////////////////////////////////////////////////////////////////////////
	// Counting sort - Step 1: Clear directory
	template < typename TDir, typename THashValue >
	inline void _qgramClearDir(TDir &dir, BucketMap<THashValue> &bucketMap)
	{
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);
		arrayFill(begin(bucketMap.qgramHash, Standard()), end(bucketMap.qgramHash, Standard()), (THashValue)-1);
	}

	// looks up the bucket for the hash
	// allocates a new bucket if non for this hash exsits yet
	// returns the position of the bucket
	template < typename THashValue, typename THashValue2 >
	inline THashValue
	requestBucket(BucketMap<THashValue> &bucketMap, THashValue2 hash)
	{
		typedef unsigned long TSize;
		// get size of the index
		TSize hlen = length(bucketMap.qgramHash) - 1;
		
		// check whether bucket map is disabled
		if (hlen == (TSize)-1)
			return hash;
		
		// where the hash should be found if no collision took place before
		TSize h1 = (TSize)(hash % hlen);
		
		// -1 is the undefiend value, hence the method works not for the largest word of length 32
		// if there is no collision with another hash value
		// the bucket is still empty
		if (bucketMap.qgramHash[h1] == (THashValue)-1)
		{
			bucketMap.qgramHash[h1] = hash;
			return h1;
		}
		// if there is a collision
		else
		{
			// the bucket for this hash was requestet before
			// return the same bucket
			if (bucketMap.qgramHash[h1] == hash)
				return h1;
			// another hash is occupying this bucket already
			else 
			{
				// TSize step = 1 + (hash % bucketMap.prime);
				TSize step = (TSize)bucketMap.prime;
				// look 'step' buckets further untill one is free or was requested by this hash earlier
				do {
					h1 = (h1 + step) % hlen;
				} while (bucketMap.qgramHash[h1] != (THashValue)-1 && bucketMap.qgramHash[h1] != hash);
				bucketMap.qgramHash[h1] = hash;
				return h1;
			}
		}
	}

	
	// looks up the bucket for the hash
	// returns the position of the bucket
	template < typename THashValue, typename THashValue2 >
	inline THashValue
	getBucket(BucketMap<THashValue> const &bucketMap, THashValue2 hash)
	{
		typedef unsigned long TSize;
		// get size of the index
		TSize hlen = length(bucketMap.qgramHash) - 1;

		// check whether bucket map is disabled
		if (hlen == (TSize)-1)
			return hash;
		
		// where the hash should be found if no collision took place before
		TSize h1 = hash % hlen;
		
		// -1 is the undefiend value, hence the method works not for the largest word of length 32
		// if there is no collision with another hash value
		// the bucket is still empty
		if (bucketMap.qgramHash[h1] == (THashValue)-1)
			return h1;
		// if there is a collision
		else
		{
			// the bucket for this hash was requestet before
			// return the same bucket
			if (bucketMap.qgramHash[h1] == hash)
				return h1;
			// another hash is occupying this bucket already
			else 
			{
				TSize step = bucketMap.prime;
				// look 'step' buckets further untill one is free or was requested by this hash earlier
				do {
					h1 = (h1 + step) % hlen;
				} while (bucketMap.qgramHash[h1] != (THashValue)-1 && bucketMap.qgramHash[h1] != hash);
				return h1;
			}
		}
	}

	template <typename TSize>
	inline unsigned coprimeTest(TSize hlen)
	{
		static unsigned _primes[42] = { 43, 47, 53, 59, 61, 67, 71, 
										73, 79, 83, 89, 97, 101, 103, 
										107, 109, 113, 127, 131, 137, 139, 
										149, 151, 157, 163, 167, 173, 179,
										181, 191, 193, 197, 199, 211, 223, 
										227, 229, 233, 239, 241, 251, 257 };

		for (int i = 0; i < 42; ++i)
			if ((hlen % _primes[i]) != 0)
				return _primes[i];
		
		return _primes[0];
	}

	template <typename TObject, typename TShapeSpec>
	inline __int64 _fullDirLength(Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> > const &index) 
	{
		typedef Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >	TIndex;
		typedef typename Fibre<TIndex, QGramDir>::Type						TDir;
		typedef typename Fibre<TIndex, FibreShape>::Type					TShape;
		typedef typename Host<TShape>::Type									TTextValue;
		typedef typename Value<TDir>::Type									TDirValue;
		typedef typename Value<TShape>::Type								THashValue;
		
		double num_qgrams = _qgramQGramCount(index) * index.alpha;
		double max_qgrams = pow((double)ValueSize<TTextValue>::VALUE, (double)weight(indexShape(index)));
		__int64 qgrams;
		
		// compare size of open adressing with 1-1 mapping and use the smaller one
		if (num_qgrams * (sizeof(TDirValue) + sizeof(THashValue)) < max_qgrams * sizeof(TDirValue))
		{
			qgrams = (__int64)ceil(num_qgrams);
			resize(const_cast<TIndex &>(index).bucketMap.qgramHash, qgrams + 1, Exact());
		} else
		{
			qgrams = (__int64)ceil(max_qgrams);
			clear(const_cast<TIndex &>(index).bucketMap.qgramHash);	// 1-1 mapping, no bucket map needed
		}
		
		const_cast<TIndex &>(index).bucketMap.prime = coprimeTest(qgrams);
		return qgrams + 1;
	}
	
}

#endif //#ifndef SEQAN_HEADER_...
