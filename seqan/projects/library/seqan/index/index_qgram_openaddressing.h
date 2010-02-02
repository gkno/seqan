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

	struct _OpenAddressing;
	typedef Tag<_OpenAddressing> OpenAddressing;
	
	template <typename THashValue>
	struct BucketMap 
	{
		String<THashValue>	qgramHash;
		THashValue			prime;
	};

	// use the index value type as shape value type
	template < typename TObject, typename TShapeSpec >
	struct Fibre< Index<TObject, Index_QGram<TShapeSpec, OpenAddressing> >, Fibre_BucketMap>
	{
		typedef typename Fibre< Index<TObject, Index_QGram<TShapeSpec, OpenAddressing> >, Fibre_Shape>::Type TShape;
		typedef typename Value<TShape>::Type	THashValue;
		typedef BucketMap<THashValue>			Type;
	};

	template < typename TObject, typename TShapeSpec >
	class Index<TObject, Index_QGram<TShapeSpec, OpenAddressing> >
	{
	public:
		typedef typename Fibre<Index, QGram_Text>::Type			TText;
		typedef typename Fibre<Index, QGram_SA>::Type			TSA;
		typedef typename Fibre<Index, QGram_Dir>::Type			TDir;
		typedef typename Fibre<Index, QGram_Counts>::Type		TCounts;
		typedef typename Fibre<Index, QGram_CountsDir>::Type	TCountsDir;
		typedef typename Fibre<Index, QGram_Shape>::Type		TShape;
		typedef typename Fibre<Index, QGram_BucketMap>::Type	TBucketMap;
		typedef typename Cargo<Index>::Type						TCargo;

		Holder<TText>	text;		// underlying text
		TSA				sa;			// suffix array sorted by the first q chars
		TDir			dir;		// bucket directory
		TCounts			counts;		// counts each q-gram per sequence
		TCountsDir		countsDir;	// directory for count buckets
		TShape			shape;		// underlying shape
		TCargo			cargo;		// user-defined cargo
		TBucketMap		bucketMap;	// bucketMap table (used by open-addressing index)

		double			alpha;

		Index():
			alpha(1.5) {}

		Index(Index &other):
			text(other.text),
			sa(other.sa),
			dir(other.dir),
			counts(other.counts),
			countsDir(other.countsDir),
			shape(other.shape),
			cargo(other.cargo),
			alpha(1.5) {}

		Index(Index const &other):
			text(other.text),
			sa(other.sa),
			dir(other.dir),
			counts(other.counts),
			countsDir(other.countsDir),
			shape(other.shape),
			cargo(other.cargo),
			alpha(1.5) {}

		template <typename _TText>
		Index(_TText &_text):
			text(_text),
			alpha(1.5) {}

		template <typename _TText>
		Index(_TText const &_text):
			text(_text),
			alpha(1.5) {}

		template <typename _TText, typename _TShape>
		Index(_TText &_text, _TShape const &_shape):
			text(_text),
			shape(_shape),
			alpha(1.5) {}

		template <typename _TText, typename _TShape>
		Index(_TText const &_text, _TShape const &_shape):
			text(_text),
			shape(_shape),
			alpha(1.5) {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// Counting sort - Step 1: Clear directory
	template < typename TDir, typename THashValue >
	inline void _qgramClearDir(TDir &dir, BucketMap<THashValue> &bucketMap)
	{
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);
		arrayFill(begin(bucketMap.qgramHash, Standard()), end(bucketMap.qgramHash, Standard()), (THashValue)-1);
	}

	template < typename THashValue, typename THashValue2 >
	inline THashValue
	requestBucket(BucketMap<THashValue> &bucketMap, THashValue2 hash)
	{
		long unsigned hlen = length(bucketMap.qgramHash) - 1;
		long unsigned h1 = hash % hlen;
		// -1 is the undefiend value, hence the method works not for the largest word of length 32
		if (bucketMap.qgramHash[h1] == (THashValue)-1)
		{
			bucketMap.qgramHash[h1] = hash;
			return h1;
		}
		else
		{
			if (bucketMap.qgramHash[h1] == hash)
				return h1;
			else 
			{
//				long unsigned step = 1 + (hash % bucketMap.prime);
				long unsigned step = bucketMap.prime;
				do {
					h1 = (h1 + step) % hlen;
				} while (bucketMap.qgramHash[h1] != (THashValue)-1 && bucketMap.qgramHash[h1] != hash);
				bucketMap.qgramHash[h1] = hash;
				return h1;
			}
		}
	}

	template < typename THashValue, typename THashValue2 >
	inline THashValue
	getBucket(BucketMap<THashValue> const &bucketMap, THashValue2 hash)
	{
		long unsigned hlen = length(bucketMap.qgramHash) - 1;
		long unsigned h1 = hash % hlen;
		// -1 is the undefiend value, hence the method works not for the largest word of length 32
		if (bucketMap.qgramHash[h1] == (THashValue)-1)
			return h1;
		else
		{
			if (bucketMap.qgramHash[h1] == hash)
				return h1;
			else 
			{
				long unsigned step = bucketMap.prime;
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
	inline int _fullDirLength(Index<TObject, Index_QGram<TShapeSpec, OpenAddressing> > const &index) 
	{
		typedef Index<TObject, Index_QGram<TShapeSpec, OpenAddressing> >	TIndex;
		typedef typename Fibre<TIndex, Fibre_Shape>::Type					TShape;
		typedef typename Host<TShape>::Type									TTextValue;
		
		unsigned long qgrams = _qgramQGramCount(index);
		double max_qgrams = pow((double)ValueSize<TTextValue>::VALUE, (double)weight(indexShape(index)));
		if (max_qgrams <= qgrams)
			qgrams = ceil(max_qgrams);
		
		qgrams = ceil(qgrams * index.alpha);
		
		const_cast<TIndex &>(index).bucketMap.prime = coprimeTest(qgrams);
		resize(const_cast<TIndex &>(index).bucketMap.qgramHash, qgrams + 1, Exact());
		return qgrams + 1;
	}
	
}

#endif //#ifndef SEQAN_HEADER_...
