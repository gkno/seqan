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
	struct NestedDirEntry {
		TPos	start;		// if the upper bit is set, position points to another directory
		TPos	end;
	};

/*	template < typename TObject, typename TSpec >
	struct Fibre< Index<TObject, TSpec>, Tag<_Fibre_Dir> > {
		typedef String< 
			typename Size< Index<TObject, TSpec> >::Type,
			Block<>
		> Type;
	};
*/
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

		Holder<TText>	text;	// underlying text
		TSA				sa;		// suffix array sorted by the first q chars
		TDir			dir;	// bucket directory
		TShape			shape;	// underlying shape

		typedef String<typename Value<TSA>::Type, Alloc<> >	tempDir;
		
		Index() {}

		template <typename _TText>
		Index(_TText &_text):
			text(_text) {}

		template <typename _TText>
		Index(_TText const &_text):
			text(_text) {}

		template <typename _TText, typename _TShape>
		Index(_TText &_text, _TShape const &_shape):
			text(_text),
			shape(_shape) {}

		template <typename _TText, typename _TShape>
		Index(_TText const &_text, _TShape const &_shape):
			text(_text),
			shape(_shape) {}
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
			typename TShape,
			typename TTempSA,
			typename TSize>
	void _sortQGramBucket(
		TSA &sa,
		TDir &dir,
		TText const &text,
		TShape &shape,
		TTempSA &temp_sa,
		TSize offset)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type			TTextIterator;
		typedef typename Iterator<TTempSA const, Standard>::Type		TSAIterator;
		typedef typename Iterator<TDir, Standard>::Type					TDirIterator;
		typedef typename Value<TShape>::Type							TValue;
		typedef typename _MakeSigned<typename Size<TText>::Type>::Type	TTextSize;
		
		// 1. clear counters and copy SA to temporary SA
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);
		temp_sa = sa;

		// 2. count q-grams
		{
			TTextIterator itText = begin(text, Standard());
			TSAIterator itSA = begin((TTempSA const &)temp_sa, Standard());
			TSAIterator itSAEnd = end((TTempSA const &)temp_sa, Standard());
			TSize textLength = length(text);
			for(; itSA != itSAEnd; ++itSA)
				++dir[hashCrossBorder(shape, itText + *itSA + offset, textLength - *itSA - offset)];
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
			TSAIterator itSA = begin((TTempSA const &)temp_sa, Standard());
			TSAIterator itSAEnd = end((TTempSA const &)temp_sa, Standard());
			TSize textLength = length(text);
			for(; itSA != itSAEnd; ++itSA)
				sa[dir[hashCrossBorder(shape, itText + *itSA + offset, textLength - *itSA - offset) + 1]++] = *itSA;
		}
	}

	template <typename TText, typename TShapeSpec, typename TPos>
	inline void _expandQGramNode(Index<TText, Index_QGram_Nested<TShapeSpec> > &index, TPos nodeID)
	{
	SEQAN_CHECKPOINT
		typedef Index<TText, Index_QGram_Nested<TShapeSpec> >		TIndex;
		typedef typename Fibre<TIndex, QGram_Nested_Dir>::Type		TDir;
		typedef typename Fibre<TIndex, QGram_Nested_SA >::Type		TSA;
		typedef typename Fibre<TIndex, QGram_Nested_Shape >::Type	TShape;

		String< typename Value<TSA>::Type, Alloc<> > tempSA;
		resize(tempSA, length(indexSA(index)));

		TShape &shape = indexShape(index);

		typename Size<TDir>::Type oldLength = length(indexDir(index));
		resize(indexDir(index), oldLength + _fullDirLength(index));
		_sortQGramBucket(indexSA(index), indexDir(index), indexText(index), shape, tempSA, 0);
	}

//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation 

	template <typename TText, typename TShapeSpec>
	inline bool indexCreate(Index<TText, Index_QGram_Nested<TShapeSpec> > &index, QGram_Nested_SA const, Default const)
	{
		typedef Index<TText, Index_QGram_Nested<TShapeSpec> >	TIndex;
		resize(indexSA(index), length(indexRawText(index)));
		clear(indexDir(index));

		_expandQGramNode(index, 0);

		return true;
	}

}


#endif //#ifndef SEQAN_HEADER_...
