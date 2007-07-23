/*
 *  index_qgram_nested.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_QGRAM_NESTED_H
#define SEQAN_HEADER_INDEX_QGRAM_NESTED_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// q-gram index fibres

	struct _Fibre_Dir_Nested;		// directory/hash table, contains start indices of buckets

	typedef Tag<_Fibre_Text>		QGram_Nested_Text;
	typedef Tag<_Fibre_RawText>		QGram_Nested_RawText;
	typedef Tag<_Fibre_SA>			QGram_Nested_SA;
	typedef Tag<_Fibre_RawSA>		QGram_Nested_RawSA;
	typedef Tag<_Fibre_Dir>			QGram_Nested_Dir;
	typedef Tag<_Fibre_SADir>		QGram_Nested_SADir;
	typedef Tag<_Fibre_Shape>		QGram_Nested_Shape;

//////////////////////////////////////////////////////////////////////////////
// nested q-gram table

	template <typename TPos>
	struct _NestedDirEntry {
		TPos	start;		// if the upper bit is set, position points to another directory
		TPos	end;
	};

	template <typename TSize>
	struct _QGramTreeNodeDescriptor {
		TSize	saBeginPos;
		TSize	saEndPos;
		TSize	dirBeginPos;
		TSize	dirEndPos;
		TSize	lcp;
		bool	rootNode;
	};


//////////////////////////////////////////////////////////////////////////////
// q-gram index

/**
.Spec.Index_QGram_Nested:
..summary:An index based on an array of sorted q-grams.
..cat:Index
..general:Class.Index
..signature:Index<TText, Index_QGram_Nested<> >
..param.TText:The text type.
...type:Class.String
..remarks:The fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a suffix array sorted by the first q characters (see @Tag.QGram_SA@) and a q-gram directory (see @Tag.QGram_Dir@).
*/

	template < typename TShapeSpec >
	struct Index_QGram_Nested;

/*
	template < typename TObject, typename TShapeSpec >
	struct Fibre< Index<TObject, Index_QGram_Nested<TShapeSpec> >, Tag<_Fibre_Dir> > 
	{
		typedef Index<TObject, Index_QGram_Nested<TShapeSpec> > TIndex;
		typedef String< 
			typename typename Size<TIndex>::Type,
			Alloc<>
		> Type;
	};
*/
	template < typename TObject, typename TShapeSpec >
	struct Fibre< Index<TObject, Index_QGram_Nested<TShapeSpec> >, Tag<_Fibre_Shape> > {
		typedef Index< TObject, Index_QGram_Nested<TShapeSpec> >	TIndex;
		typedef Shape< typename Value<TIndex>::Type, TShapeSpec >	Type;
	};

	
	template < typename TObject, typename TShapeSpec >
	class Index<TObject, Index_QGram_Nested<TShapeSpec> > {
	public:
		typedef typename Fibre<Index, QGram_Nested_Text>::Type		TText;
		typedef typename Fibre<Index, QGram_Nested_SA>::Type		TSA;
		typedef typename Fibre<Index, QGram_Nested_Dir>::Type		TDir;
		typedef typename Fibre<Index, QGram_Nested_Shape>::Type		TShape;

		typedef typename Size<TText>::Type							TSize;

		static TSize const SUBDIR        =  1 << (BitsPerValue<TSize>::VALUE - 1);
		static TSize const SORTED_SUBDIR =  1 << (BitsPerValue<TSize>::VALUE - 2);
		static TSize const SORTED_SA     =  1 << (BitsPerValue<TSize>::VALUE - 2);
		static TSize const BITMASK       = ~(SUBDIR | SORTED_SA);


		Holder<TText>	text;	// underlying text
		TSA				sa;		// suffix array sorted by the first q chars
		TDir			dir;	// bucket directory
		TShape			shape;	// underlying shape

		TSize			maxLcp;			// maximal qgram-tree depth
		TSize			minSufCount;	// minimal required sa entries per bucket
		TSize			minDirSize;		// if dir contains less entries it will be compressed

		String<typename Value<TSA>::Type, Alloc<> >		tempSA;
		String<typename Value<TShape>::Type, Alloc<> >	hashSet;
		
		Index():
			maxLcp(10),
			minSufCount(0*1024),
			minDirSize(0) {}

		template <typename _TText>
		Index(_TText &_text):
			text(_text),
			maxLcp(10),
			minSufCount(1024) {}

		template <typename _TText>
		Index(_TText const &_text):
			text(_text),
			maxLcp(10),
			minSufCount(1024) {}

		template <typename _TText, typename _TShape>
		Index(_TText &_text, _TShape const &_shape):
			text(_text),
			shape(_shape),
			maxLcp(10),
			minSufCount(1024) {}

		template <typename _TText, typename _TShape>
		Index(_TText const &_text, _TShape const &_shape):
			text(_text),
			shape(_shape),
			maxLcp(10),
			minSufCount(1024) {}
	};

    template < typename TText, typename TSpec >
    struct Value< Index<TText, Index_QGram_Nested<TSpec> > > {
		typedef typename Value< typename Fibre< Index<TText, Index_QGram_Nested<TSpec> >, QGram_RawText >::Type >::Type Type;
    };

	template < typename TText, typename TSpec >
    struct Size< Index<TText, Index_QGram_Nested<TSpec> > > {
		typedef typename Size< typename Fibre< Index<TText, Index_QGram_Nested<TSpec> >, QGram_RawText >::Type >::Type Type;
    };


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, Index_QGram_Nested<TSpec> >, Tag<_Fibre_SA> > {
        typedef Default Type;
    };

//////////////////////////////////////////////////////////////////////////////
/**
.Function.createQGramIndex:
..summary:Builds a q-gram index on a sequence. 
..cat:Index
..signature:createQGramIndex(sa, dir, text, shape)
..param.text:The sequence.
..param.shape:The shape to be used.
...type:Class.Shape
..param.sa:The resulting list in which all q-grams are sorted alphabetically.
..param.dir:The resulting array that indicates at which position in index the corresponding q-grams can be found.
..returns:Index contains the sorted list of qgrams. For each possible q-gram pos contains the first position in index that corresponds to this q-gram. 
*/

	template < typename TSA, 
			typename TDir,
			typename TText,
			typename TShape>
	void _sortFirstQGramBucket(
		TSA &sa,
		TDir &dir,
		TText const &text,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type			TTextIterator;
		typedef typename Iterator<TDir, Standard>::Type					TDirIterator;
		typedef typename Value<TShape>::Type							TValue;
		typedef typename Size<TText>::Type								TSize;
		
		// 1. clear counters and copy SA to temporary SA
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		// 2. count q-grams
		{
			TTextIterator itText = begin(text, Standard());
			TSize rest = length(text);
			++dir[hash2(shape, itText, rest)];
			if (rest > 1) {
				--rest;
				do {
					++itText;
					++dir[hash2Next(shape, itText, rest)];
				} while (rest-- != 0);
			}
		}

		// 3. cumulative sum
		{
			TDirIterator it = begin(dir, Standard());
			TDirIterator itEnd = end(dir, Standard());
			TSize diff = 0, diff_prev = 0, sum = 0;
			while (it != itEnd) {
				sum += diff_prev;
				diff_prev = diff;
				diff = *it;
				*it = sum;
				++it;
			}
			SEQAN_ASSERT(sum + diff_prev == (TSize)length(sa));
		}
		
		// 4. fill suffix array
		{
			TTextIterator itText = begin(text, Standard());
			TSize rest = length(text);
			sa[dir[hash2(shape, itText, rest) + 1]++] = 0;					// first hash
			if (rest > 1) {
				--rest;
				TSize i = 0;
				do {
					++itText; ++i;
					sa[dir[hash2Next(shape, itText, rest) + 1]++] = i;		// next hash
				} while (rest-- != 0);
			}
		}
	}

	template < typename TSA, 
			typename TDir,
			typename TText,
			typename TShape,
			typename TTempSA,
			typename TSize>
	TSize _sortQGramBucket(
		TSA &sa,
		TDir &dir,
		TText const &text,
		TShape &shape,
		TTempSA &temp_sa,
		TSize lcp,
		TSize saOffset)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type		TTextIterator;
		typedef typename Iterator<TTempSA const, Standard>::Type	TSAIterator;
		typedef typename Iterator<TDir, Standard>::Type				TDirIterator;
		typedef typename Value<TShape>::Type						TValue;
		typedef typename Size<TText>::Type							TTextSize;

		// 1. clear counters and copy SA to temporary SA
		TDirIterator dirBegin = begin(dir, Standard());

		arrayFill(dirBegin, end(dir, Standard()), 0);
		temp_sa = sa;

		// 2. count q-grams
		{
			TTextIterator itText = begin(text, Standard()) + lcp;
			TSAIterator itSA = begin((TTempSA const &)temp_sa, Standard());
			TSAIterator itSAEnd = end((TTempSA const &)temp_sa, Standard());
			TTextSize textRest = length(text) - lcp;
			for(; itSA != itSAEnd; ++itSA)
				if (textRest >= *itSA)
					++*(dirBegin + hash2(shape, itText + *itSA, textRest - *itSA));
		}

		TSize dirtyBuckets = 0;

		// 3. cumulative sum
		{
			TDirIterator it = begin(dir, Standard());
			TDirIterator itEnd = end(dir, Standard());
			TSize diff = 0, diff_prev = 0, sum = 0;
			while (it != itEnd) {
				sum += diff_prev;
				diff_prev = diff;
				if (diff = *it) ++dirtyBuckets;
				*it = sum;
				++it;
			}
			SEQAN_ASSERT(sum + diff_prev == (TSize)length(sa));
		}
		
		// 4. fill suffix array
		{
			TTextIterator itText = begin(text, Standard()) + lcp;
			TSAIterator itSA = begin((TTempSA const &)temp_sa, Standard());
			TSAIterator itSAEnd = end((TTempSA const &)temp_sa, Standard());
			TTextSize textRest = length(text) - lcp;
			for(; itSA != itSAEnd; ++itSA)
				if (textRest >= *itSA)
					sa[(*(dirBegin + hash2(shape, itText + *itSA, textRest - *itSA) + 1))++] = *itSA;
		}

		// 5. add saOffset to hash directory entries
		{
			TDirIterator it = begin(dir, Standard());
			TDirIterator itEnd = end(dir, Standard());
			for(; it != itEnd; ++it)
				*it += saOffset;
		}

		return dirtyBuckets;
	}


	template <typename TText, typename TShapeSpec, typename TPos>
	void _compressQGramBucket(
		Index<TText, Index_QGram_Nested<TShapeSpec> > &index, 
		TPos pos)
	{
	SEQAN_CHECKPOINT
		typedef Index<TText, Index_QGram_Nested<TShapeSpec> >			TIndex;
		typedef typename Fibre<TIndex, QGram_Nested_Dir>::Type			TDir;
		typedef typename Fibre<TIndex, QGram_Nested_Shape >::Type		TShape;
		typedef typename Iterator<TDir, Standard>::Type					TDirIterator;
		typedef typename Value<TDir>::Type								TValue;
		typedef typename Size<TIndex>::Type								TSize;

		typedef String<typename Value<TShape>::Type, Alloc<> >			THashSet;
		typedef typename Iterator<THashSet, Standard>::Type				THSIterator;
		
		// 1. clear set of hash values
		TSize dirLength = _fullDir2Length(index) + 1;
		resize(index.hashSet, index.minDirSize);
		arrayFill(begin(index.hashSet, Standard()), end(index.hashSet, Standard()), false);

		// 2. remove duplicate entries
		{
			TDirIterator itSrc = begin(indexDir(index), Standard()) + pos;	// original iterator
			TDirIterator itDst = itSrc + 1;						// compressed iterator
			TDirIterator itEnd = end(indexDir(index), Standard());
			THSIterator  itHS = begin(index.hashSet, Standard());

			// copy first src entry to second dst entry
			TValue old = *itSrc; ++itSrc;
			TSize size = 0;
			
			for (TSize i = 0; itSrc != itEnd; ++itSrc, ++i)
				if (old != *itSrc) {
					TValue tmp = *itSrc;
					*itDst = old; ++itDst;
					old = tmp;
					*itHS = i; ++itHS;
					++size;
				}

			indexDir(index)[pos] = size;

			THSIterator itHS2 = begin(index.hashSet, Standard());
			for (; itHS2 != itHS; ++itDst, ++itHS2)
				*itDst = *itHS2;
		}
	}


	template <typename TText, typename TShapeSpec, typename TPos>
	inline TPos
	_firstSuffixOfQGramBucket(
		Index<TText, Index_QGram_Nested<TShapeSpec> > &index, 
		TPos pos)
	{
		typedef Index<TText, Index_QGram_Nested<TShapeSpec> >		TIndex;
		typedef typename Fibre<TIndex, QGram_Nested_Dir>::Type		TDir;
		typedef typename Iterator<TDir, Standard>::Type				TDirIterator;

		TDirIterator dirBegin = begin(indexDir(index), Standard());
		pos = *(dirBegin + pos);
		while (pos & index.SUBDIR)
			pos = *(dirBegin + (pos & index.BITMASK));
		return pos & index.BITMASK;
	}

	template <typename TText, typename TShapeSpec, typename TPos, typename TSize>
	inline TPos 
	_expandQGramNode(
		Index<TText, Index_QGram_Nested<TShapeSpec> > &index, 
		TPos pos,
		TSize lcp)
	{
	SEQAN_CHECKPOINT
		typedef Index<TText, Index_QGram_Nested<TShapeSpec> >		TIndex;
		typedef typename Fibre<TIndex, QGram_Nested_Dir>::Type		TDir;
		typedef typename Fibre<TIndex, QGram_Nested_SA >::Type		TSA;
		typedef typename Fibre<TIndex, QGram_Nested_Shape >::Type	TShape;

		typedef typename Value<TDir>::Type							TDirValue;

		// expand directory table
		TSize sa1 = _firstSuffixOfQGramBucket(index, pos);
		TSize sa2 = _firstSuffixOfQGramBucket(index, pos + 1);

		// sub-directory contains too little suffices
		// maximal qgram-tree depth reached -> do a complete quicksort
		if (sa2 - sa1 <= index.minSufCount || lcp >= index.maxLcp) {
			::std::cout << "build bucket [" << sa1 << "," << sa2 << "] with qsort  lcp=" << lcp << ::std::flush;
			typename Infix<TSA>::Type	infSA = infix(indexSA(index), sa1, sa2);
			_sortBucketQuickSort(infSA, indexText(index), lcp);
			::std::cout << "." << ::std::endl;
			return indexDir(index)[pos] = (sa1 | index.SORTED_SA);
		}


		::std::cout << "build bucket " << pos << "  lcp=" << lcp << ::std::flush;

		TShape &shape = indexShape(index);

		// get space for new directory
		TDirValue dir1 = length(indexDir(index));
		TDirValue dir2 = dir1 + _fullDir2Length(index) + 1;
		resize(indexDir(index), dir2, Generous());

		// sort bucket
		typename Infix<TSA>::Type	infSA  = infix(indexSA(index), sa1, sa2);
		typename Infix<TDir>::Type	infDir = infix(indexDir(index), dir1, dir2);

		TDirValue size = _sortQGramBucket(
				infSA, infDir, indexText((TIndex const&)index),
				shape, index.tempSA, lcp, sa1);

		if (size < index.minDirSize) 
		{
			_compressQGramBucket(index, dir1);
			dir2 = dir1 + 2 * size + 1;
			resize(indexDir(index), dir2, Generous());
			dir1 = (dir1 + 1) | index.SORTED_SUBDIR;
		}


		::std::cout << "." << ::std::endl;
		// set link to sub-directory
		return indexDir(index)[pos] = (dir1 | index.SUBDIR);
	}

	template <typename TText, typename TShapeSpec>
	inline void
	_dump(Index<TText, Index_QGram_Nested<TShapeSpec> > &index)
	{
		typedef Index<TText, Index_QGram_Nested<TShapeSpec> >		TIndex;
		typedef typename Fibre<TIndex, QGram_Nested_Dir>::Type		TDir;
		typedef typename Fibre<TIndex, QGram_Nested_SA >::Type		TSA;
		typedef typename Fibre<TIndex, QGram_Nested_Shape >::Type	TShape;

		typedef typename Value<TDir>::Type							TDirValue;

		::std::cout << ::std::endl;
		for(int i=0; i < length(indexDir(index)); ++i) {
			TDirValue d = indexDir(index)[i];
			::std::cout << i << ":  " << (d & index.BITMASK);
			if (d & index.SUBDIR)    ::std::cout << "  (SubDir)";
			if (d & index.SORTED_SA) ::std::cout << "  (SubSA)";
			::std::cout << ::std::endl;
		}

		::std::cout << ::std::endl;
		for(int i=0; i < length(indexSA(index)); ++i)
			::std::cout << i << ":  " << indexSA(index)[i] << "  " << suffix(indexText(index), indexSA(index)[i]) << ::std::endl;

	}

	template <typename TText, typename TShapeSpec, typename TPos, typename THash>
	inline TPos
	_findHash(Index<TText, Index_QGram_Nested<TShapeSpec> > &index, TPos subDir, THash hash)
	{
		
	}


	template <typename TText, typename TShapeSpec, typename TPattern>
	inline Pair<unsigned>
	equalRange(Index<TText, Index_QGram_Nested<TShapeSpec> > &index, TPattern &p)
	{
		typedef Index<TText, Index_QGram_Nested<TShapeSpec> >		TIndex;
		typedef typename Fibre<TIndex, QGram_Nested_SA>::Type		TSA;
		typedef typename Fibre<TIndex, QGram_Nested_Dir>::Type		TDir;
		typedef typename Fibre<TIndex, QGram_Nested_Shape >::Type	TShape;

		typedef typename Iterator<TPattern, Standard>::Type			TPatternIterator;
		typedef typename Value<TDir>::Type							TDirValue;
		typedef typename Size<TIndex>::Type							TSize;
		typedef Pair<unsigned>										TPair;

		TShape &shape = indexShape(index);
		TDir &dir = indexDir(index);

		TPatternIterator it = begin(p, Standard());
		TSize rest = length(p);
		TSize lcp = 0;
		TSize dirOffset = 0;
		for(; rest > length(shape); rest -= length(shape), it += length(shape))
		{
			TDirValue pos = dirOffset + hash2(shape, it, rest);
			TDirValue subDir = dir[pos];
			if ((subDir & index.SUBDIR) == 0) {
				// no sub-directory -> create one
				lcp += length(shape);
				if ((subDir & index.SORTED_SA) == 0)
					subDir = _expandQGramNode(index, pos, lcp);

				_dump(index);
				if (subDir & index.SORTED_SA) {
					subDir &= index.BITMASK;
					TSize sa2 = _firstSuffixOfQGramBucket(index, pos + 1);
					typedef typename Infix<TSA>::Type TInfSA;
					typedef Pair<typename Iterator<TInfSA const>::Type> TItPair;
					TItPair itP = equalRangeSA(
						indexText(index), 
						infix(indexSA(index), subDir & index.BITMASK, sa2),
						p);
					TPair range;
					range.i1 = position(itP.i1) + subDir;
					range.i2 = position(itP.i2) + subDir;
					return range;
				}
			}
			dirOffset = subDir & index.BITMASK;
		}

		TPair range;
		TDirValue pos1 = dirOffset + hash2(shape, it, rest);
		TDirValue pos2 = dirOffset + hash2Upper(shape, it, rest);
		range.i1 = _firstSuffixOfQGramBucket(index, pos1);
		range.i2 = _firstSuffixOfQGramBucket(index, pos2);

		return range;
	}

//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation 

	template <typename TText, typename TShapeSpec>
	inline bool indexCreate(Index<TText, Index_QGram_Nested<TShapeSpec> > &index, QGram_Nested_SA const, Default const)
	{
		typedef Index<TText, Index_QGram_Nested<TShapeSpec> >	TIndex;
		resize(indexSA(index), length(indexRawText(index)) + 1);
		resize(indexDir(index), _fullDir2Length(index) + 1);

		index.minDirSize = length(indexDir(index)) / 4;

		_sortFirstQGramBucket(indexSA(index), indexDir(index), indexText(index), indexShape(index));
		return true;
	}

}


#endif //#ifndef SEQAN_HEADER_...
