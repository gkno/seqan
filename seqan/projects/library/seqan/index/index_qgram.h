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

#ifndef SEQAN_HEADER_INDEX_QGRAM_H
#define SEQAN_HEADER_INDEX_QGRAM_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// q-gram index fibres

	struct _Fibre_Dir;			// directory/hash table, contains start indices of buckets
	struct _Fibre_SADir;		// identifies algorithm to construct both SA and directory at once
	struct _Fibre_Shape;		// underlying shape
	struct _Fibre_Counts;		// counts each q-gram
	struct _Fibre_CountsDir;	// directory for counts buckets

	typedef Tag<_Fibre_Text> const		QGram_Text;
	typedef Tag<_Fibre_RawText> const	QGram_RawText;
	typedef Tag<_Fibre_SA> const		QGram_SA;
	typedef Tag<_Fibre_RawSA> const		QGram_RawSA;
	typedef Tag<_Fibre_Dir> const		QGram_Dir;
	typedef Tag<_Fibre_SADir> const		QGram_SADir;
	typedef Tag<_Fibre_Shape> const		QGram_Shape;
	typedef Tag<_Fibre_Counts> const	QGram_Counts;
	typedef Tag<_Fibre_CountsDir> const	QGram_CountsDir;

//////////////////////////////////////////////////////////////////////////////
// q-gram index

/**
.Spec.Index_QGram:
..summary:An index based on an array of sorted q-grams.
..cat:Index
..general:Class.Index
..signature:Index<TText, Index_QGram<> >
..param.TText:The text type.
...type:Class.String
..remarks:The fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a suffix array sorted by the first q characters (see @Tag.QGram_SA@) and a q-gram directory (see @Tag.QGram_Dir@).
*/

	template < typename TShapeSpec >
	struct Index_QGram;


	template < typename TObject, typename TShapeSpec >
	struct Fibre< Index<TObject, Index_QGram<TShapeSpec> >, Tag<_Fibre_Shape> const > {
		typedef Index< TObject, Index_QGram<TShapeSpec> >			TIndex;
		typedef Shape< typename Value<TIndex>::Type, TShapeSpec >	Type;
	};

	template < typename TObject, typename TShapeSpec >
	class Index<TObject, Index_QGram<TShapeSpec> > {
	public:
		typedef typename Fibre<Index, QGram_Text>::Type			TText;
		typedef typename Fibre<Index, QGram_SA>::Type			TSA;
		typedef typename Fibre<Index, QGram_Dir>::Type			TDir;
		typedef typename Fibre<Index, QGram_Counts>::Type		TCounts;
		typedef typename Fibre<Index, QGram_CountsDir>::Type	TCountsDir;
		typedef typename Fibre<Index, QGram_Shape>::Type		TShape;

		Holder<TText>	text;	// underlying text
		TSA				sa;		// suffix array sorted by the first q chars
		TDir			dir;	// bucket directory
		TCounts			counts;	// counts each q-gram per sequence
		TCountsDir		countsDir;	// directory for count buckets
		TShape			shape;	// underlying shape

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
    struct Value< Index<TText, Index_QGram<TSpec> > > {
		typedef typename Value< typename Fibre< Index<TText, Index_QGram<TSpec> >, QGram_RawText >::Type >::Type Type;
    };

	template < typename TText, typename TSpec >
    struct Size< Index<TText, Index_QGram<TSpec> > > {
		typedef typename Size< typename Fibre< Index<TText, Index_QGram<TSpec> >, QGram_RawText >::Type >::Type Type;
    };


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, Index_QGram<TSpec> >, Tag<_Fibre_SA> const > {
        typedef Default Type;
    };

//////////////////////////////////////////////////////////////////////////////
// counts array type

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, Tag<_Fibre_Counts> const > {
		typedef String<
				Pair<
					typename Size< TText >::Type,
					typename Size< Index<TText, TSpec> >::Type
				>,
				typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type 
		> Type;
	};


//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Tag<_Fibre_Dir> const >::Type & 
	getFibre(Index<TText, TSpec> &index, Tag<_Fibre_Dir> const) {
		return index.dir;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Tag<_Fibre_Dir> const >::Type & 
	getFibre(Index<TText, TSpec> const &index, Tag<_Fibre_Dir> const) {
		return index.dir;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Tag<_Fibre_Counts> const >::Type & 
	getFibre(Index<TText, TSpec> &index, Tag<_Fibre_Counts> const) {
		return index.counts;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Tag<_Fibre_Counts> const >::Type & 
	getFibre(Index<TText, TSpec> const &index, Tag<_Fibre_Counts> const) {
		return index.counts;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Tag<_Fibre_CountsDir> const >::Type & 
	getFibre(Index<TText, TSpec> &index, Tag<_Fibre_CountsDir> const) {
		return index.countsDir;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Tag<_Fibre_CountsDir> const >::Type & 
	getFibre(Index<TText, TSpec> const &index, Tag<_Fibre_CountsDir> const) {
		return index.countsDir;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Tag<_Fibre_Shape> const >::Type & 
	getFibre(Index<TText, TSpec> &index, Tag<_Fibre_Shape> const) {
		return index.shape;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Tag<_Fibre_Shape> const >::Type & 
	getFibre(Index<TText, TSpec> const &index, Tag<_Fibre_Shape> const) {
		return index.shape;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexDir:
..summary:Shortcut for $getFibre(.., QGram_Dir)$.
..cat:Index
..signature:indexDir(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_QGram
..returns:A reference to the @Tag.QGram_Dir@ fibre (q-gram directory).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Tag<_Fibre_Dir> const >::Type & 
	indexDir(Index<TText, TSpec> &index) { 
		return getFibre(index, Tag<_Fibre_Dir>()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Tag<_Fibre_Dir> const >::Type & 
	indexDir(Index<TText, TSpec> const &index) { 
		return getFibre(index, Tag<_Fibre_Dir>()); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.dirAt:
..summary:Shortcut for $value(indexDir(..), ..)$.
..cat:Index
..signature:dirAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_QGram
..returns:A reference or proxy to the value.
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, Tag<_Fibre_Dir> >::Type>::Type dirAt(TPos i, TIndex &index) {
		return value(getFibre(index, Tag<_Fibre_Dir> ()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, Tag<_Fibre_Dir> >::Type>::Type dirAt(TPos i, TIndex const &index) {
		return value(getFibre(index, Tag<_Fibre_Dir> ()), i);
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexShape:
..summary:Shortcut for $getFibre(.., QGram_Shape)$.
..cat:Index
..signature:indexShape(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_QGram
..returns:Returns a reference to the @Class.Shape@ object of a q-gram index.
Formally, this is a reference to the @Tag.QGram_Shape@ fibre.
...type:Class.Shape
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Tag<_Fibre_Counts> const >::Type & 
	indexCounts(Index<TText, TSpec> &index) {
		return getFibre(index, Tag<_Fibre_Counts>()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Tag<_Fibre_Counts> const >::Type & 
	indexCounts(Index<TText, TSpec> const &index) {
		return getFibre(index, Tag<_Fibre_Counts>()); 
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Tag<_Fibre_CountsDir> const >::Type & 
	indexCountsDir(Index<TText, TSpec> &index) {
		return getFibre(index, Tag<_Fibre_CountsDir>()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Tag<_Fibre_CountsDir> const >::Type & 
	indexCountsDir(Index<TText, TSpec> const &index) {
		return getFibre(index, Tag<_Fibre_CountsDir>()); 
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, Tag<_Fibre_Shape> const >::Type & 
	indexShape(Index<TText, TSpec> &index) { 
		return getFibre(index, Tag<_Fibre_Shape>()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, Tag<_Fibre_Shape> const >::Type & 
	indexShape(Index<TText, TSpec> const &index) { 
		return getFibre(index, Tag<_Fibre_Shape>()); 
	}

	template <typename TIndex>
	inline int _fullDirLength(TIndex const &index) 
	{
		typedef typename Fibre<TIndex, Tag<_Fibre_Shape> const >::Type	TShape;
		typedef typename Value<TIndex>::Type							TValue;
		return _intPow((unsigned)ValueSize<TValue>::VALUE, weight(indexShape(index))) + 1;
	}

	template <typename TIndex>
	inline int _fullDir2Length(TIndex const &index) 
	{
		typedef typename Fibre<TIndex, Tag<_Fibre_Shape> const >::Type	TShape;
		typedef typename Value<TIndex>::Type							TValue;
		return (_intPow(
					(unsigned)ValueSize<TValue>::VALUE,
					weight(indexShape(index)) + 1) - 1)
				/ ((unsigned)ValueSize<TValue>::VALUE - 1) + 1;
	}


    //////////////////////////////////////////////////////////////////////////////
    // QGramLess
	//
	// compare two q-grams of a given text (q-grams can be smaller than q)
    template < typename TSAValue, typename TText >
	struct _QGramLess : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TText, Standard>::Type TIter;
		TIter _begin, _end;
		typename Size<TText>::Type _q;

		template <typename TSize>
		_QGramLess(TText &text, TSize q): 
			_begin(begin(text, Standard())),
			_end(end(text, Standard())),
			_q(q) {}

		// skip the first <offset> characters
		template <typename TSize1, typename TSize2>
		_QGramLess(TText &text, TSize1 q, TSize2 offset): 
			_begin(begin(text, Standard()) + offset),
			_end(end(text, Standard())),
			_q(q) {}

		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
			if (a == b) return false;
			TIter itA = _begin + a;
			TIter itB = _begin + b;
			if (a <= b) {
				TIter itEnd = itB + _q;
				if (itEnd > _end)
					itEnd = _end;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd = itA + _q;
				if (itEnd > _end)
					itEnd = _end;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};

	template < typename TSAValue, typename TString, typename TSpec >
	struct _QGramLess<TSAValue, StringSet<TString, TSpec> const> : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TString, Standard>::Type TIter;
		StringSet<TString, TSpec> const &_stringSet;
		typename Size<TString>::Type _q;

		template <typename TSize>
		_QGramLess(StringSet<TString, TSpec> const &text, TSize q): 
			_stringSet(text),
			_q(q) {}

		inline bool operator() (TSAValue const &a, TSAValue const &b) const 
		{
			if (a == b) return false;
			TIter itA = begin(_stringSet[getValueI1(a)], Standard()) + getValueI2(a);
			TIter itB = begin(_stringSet[getValueI1(b)], Standard()) + getValueI2(b);
			if (suffixLength(a, _stringSet) > suffixLength(b, _stringSet)) {
				TIter _end = end(_stringSet[getValueI1(b)], Standard());
				TIter itEnd = itB + _q;
				if (itEnd > _end)
					itEnd = _end;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter _end = end(_stringSet[getValueI1(a)], Standard());
				TIter itEnd = itA + _q;
				if (itEnd > _end)
					itEnd = _end;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // QGramLessOffset
	//
	// compare two q-grams of a given text and skip the first <offset> characters
    template < typename TSAValue, typename TText >
	struct _QGramLessOffset :
		_QGramLess<TSAValue, TText>
	{
		template <typename TSize1, typename TSize2>
		_QGramLessOffset(TText &text, TSize1 q, TSize2 offset): 
			_QGramLess<TSAValue, TText> (text, q, offset) {}
	};

	template < typename TSAValue, typename TString, typename TSpec >
	struct _QGramLessOffset<TSAValue, StringSet<TString, TSpec> const> : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TString, Standard>::Type	TIter;
		typedef typename Size<TString>::Type				TSize;
		StringSet<TString, TSpec> const &_stringSet;
		typename Size<TString>::Type _q, _offset;

		template <typename TSize1, typename TSize2>
		_QGramLessOffset(StringSet<TString, TSpec> const &text, TSize1 q, TSize2 offset): 
			_stringSet(text),
			_q(q),
			_offset(offset) {}

		inline bool operator() (TSAValue const &a, TSAValue const &b) const 
		{
			if (a == b) return false;
			TString const &sA = _stringSet[getValueI1(a)];
			TString const &sB = _stringSet[getValueI1(b)];
			TIter itA = begin(sA, Standard()) + getValueI2(a) + _offset;
			TIter itB = begin(sB, Standard()) + getValueI2(b) + _offset;
			TSize restA = length(sA) - getValueI2(a);
			TSize restB = length(sB) - getValueI2(b);
			if (restA > restB) {
				TIter itEnd;
				if (restB >= _q)
					itEnd = itB + _q;
				else
					itEnd = end(sB, Standard());
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd;
				if (restA >= _q)
					itEnd = itA + _q;
				else
					itEnd = end(sA, Standard());
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // QGramLessNoCheck
	//
	// compare two q-grams of a given text (no check for q-grams smaller than q)
    template < typename TSAValue, typename TText >
	struct _QGramLessNoCheck : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TText, Standard>::Type TIter;
		TIter _begin;
		typename Size<TText>::Type _q;

		template <typename TSize>
		_QGramLessNoCheck(TText &text, TSize q): 
			_begin(begin(text, Standard())),
			_q(q) {}

		// skip the first <offset> characters
		template <typename TSize>
		_QGramLessNoCheck(TText &text, TSize q, TSize offset): 
			_begin(begin(text, Standard()) + offset),
			_q(q) {}

		inline bool operator() (TSAValue const a, TSAValue const b) const 
		{
			if (a == b) return false;
			TIter itA = _begin + a;
			TIter itB = _begin + b;
			if (a <= b) {
				TIter itEnd = itB + _q;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd = itA + _q;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};

    template < typename TSAValue, typename TString, typename TSpec >
	struct _QGramLessNoCheck<TSAValue, StringSet<TString, TSpec> const> : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TString, Standard>::Type TIter;
		StringSet<TString, TSpec> const &_stringSet;
		typename Size<TString>::Type _q;

		template <typename TSize>
		_QGramLessNoCheck(StringSet<TString, TSpec> const &text, TSize q): 
			_stringSet(text),
			_q(q) {}

		inline bool operator() (TSAValue const &a, TSAValue const &b) const 
		{
			if (a == b) return false;
			TIter itA = begin(_stringSet[getValueI1(a)], Standard()) + getValueI2(a);
			TIter itB = begin(_stringSet[getValueI1(b)], Standard()) + getValueI2(b);
			if (suffixLength(a, _stringSet) > suffixLength(b, _stringSet)) {
				TIter itEnd = itB + _q;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd = itA + _q;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // _QGramLessNoCheckOffset
	//
	// compare two q-grams of a given text and skip the first <offset> characters
    template < typename TSAValue, typename TText >
	struct _QGramLessNoCheckOffset: _QGramLessNoCheck<TSAValue, TText> 
	{
		template <typename TSize1, typename TSize2>
		_QGramLessNoCheckOffset(TText &text, TSize1 q, TSize2 offset): 
			_QGramLessNoCheck<TSAValue, TText> (text, q, offset) {}
	};

	template < typename TSAValue, typename TString, typename TSpec >
	struct _QGramLessNoCheckOffset<TSAValue, StringSet<TString, TSpec> const> : 
		public ::std::binary_function < TSAValue, TSAValue, bool >
    {
		typedef typename Iterator<TString, Standard>::Type TIter;
		StringSet<TString, TSpec> const &_stringSet;
		typename Size<TString>::Type _q, _offset;

		template <typename TSize1, typename TSize2>
		_QGramLessNoCheckOffset(StringSet<TString, TSpec> const &text, TSize1 q, TSize2 offset): 
			_stringSet(text),
			_q(q),
			_offset(offset) {}

		inline bool operator() (TSAValue const &a, TSAValue const &b) const 
		{
			if (a == b) return false;
			TIter itA = begin(_stringSet[getValueI1(a)], Standard()) + getValueI2(a) + _offset;
			TIter itB = begin(_stringSet[getValueI1(b)], Standard()) + getValueI2(b) + _offset;
			if (a <= b) {
				TIter itEnd = itB + _q;
				for(; itB != itEnd; ++itB, ++itA) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return false;
			} else {
				TIter itEnd = itA + _q;
				for(; itA != itEnd; ++itA, ++itB) {
					if (lexLess(*itA, *itB)) return true;
					if (lexLess(*itB, *itA)) return false;
				}
				return true;
			}
		}
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

	template <
		typename TSA, 
		typename TDir,
		typename TText,
		typename TShape >
	void createQGramIndex(
		TSA &sa,
		TDir &dir,
		TText const &text,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TIterator;
		typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
		typedef typename Value<TShape>::Type					TValue;
		typedef typename Size<TText>::Type						TSize;
		
		// 1. clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		// 2. count q-grams
		if (length(text) < length(shape)) return;
		TSize num_qgrams = length(text) - length(shape) + 1;

		TIterator itText = begin(text, Standard());
		++dir[hash(shape, itText)];
		for(TSize i = 1; i < num_qgrams; ++i)
		{
			++itText;
			++dir[hashNext(shape, itText)];
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
			SEQAN_ASSERT(sum + diff_prev == length(sa));
		}
		
		// 4. fill suffix array
		itText = begin(text, Standard());
		sa[dir[hash(shape, itText) + 1]++] = 0;				// first hash
		for(TSize i = 1; i < num_qgrams; ++i)
		{
			++itText;
			sa[dir[hashNext(shape, itText) + 1]++] = i;		// next hash
		}
	}

	template < 
		typename TSA, 
		typename TDir,
		typename TString,
		typename TSpec,
		typename TShape >
	void createQGramIndex(
		TSA &sa,
		TDir &dir,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Iterator<TDir, Standard>::Type				TDirIterator;
		typedef typename Value<TShape>::Type						TValue;
		typedef typename Size<TString>::Type						TSize;
		
		// 1. clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		// 2. count q-grams
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = length(sequence) - length(shape) + 1;

			TIterator itText = begin(sequence, Standard());
			++dir[hash(shape, itText)];
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				++dir[hashNext(shape, itText)];
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
			SEQAN_ASSERT(sum + diff_prev == length(sa));
		}
		
		// 4. fill suffix array
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = length(sequence) - length(shape) + 1;

			typename Value<TSA>::Type localPos;
			assignValueI1(localPos, seqNo);
			assignValueI2(localPos, 0);

			TIterator itText = begin(sequence, Standard());
			sa[dir[hash(shape, itText) + 1]++] = localPos;				// first hash
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				assignValueI2(localPos, i);
				sa[dir[hashNext(shape, itText) + 1]++] = localPos;		// next hash
			}
		}
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.createQGramIndexSAOnly:
..summary:Builds the suffix array of a q-gram index on a sequence. 
..cat:Index
..signature:createQGramIndexSAOnly(sa, text, shape)
..param.text:The sequence.
..param.shape:The shape to be used. q is the length of this shape
...type:Class.Shape
..param.sa:The resulting list in which all q-grams are sorted alphabetically.
*/

	template < 
		typename TSA, 
		typename TText,
		typename TShape >
	void createQGramIndexSAOnly(
		TSA &sa,
		TText const &text,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Size<TSA>::Type TSize;
		typedef typename Iterator<TSA, Standard>::Type TIter;

		// 1. Fill suffix array with a permutation (the identity)
		TIter it = begin(sa, Standard());
		TIter itEnd = end(sa, Standard());
		TSize i = 0;
		for(; it != itEnd; ++it, ++i)
			*it = i;

		// 2. Sort suffix array with quicksort
		TSize q = length(shape);
		if (i + q > length(text) + 1)
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				_QGramLess<typename Value<TSA>::Type, TText const>(text, q));
		else
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				_QGramLessNoCheck<typename Value<TSA>::Type, TText const>(text, q));
	}

	template < 
		typename TSA, 
		typename TString,
		typename TSpec,
		typename TShape >
	void createQGramIndexSAOnly(
		TSA &sa,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TSA, Standard>::Type	TIter;
		typedef typename Value<TSA>::Type				TValue;
		typedef typename Size<TString>::Type			TSize;
		typedef StringSet<TString, TSpec>				TStringSet;
		
		// 1. Fill suffix array with a permutation (the identity)
		TIter it = begin(sa, Standard());
		TIter itEnd = end(sa, Standard());
		TValue pair;
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			assignValueI1(pair, seqNo);
			TSize i = 0;
			for(; it != itEnd; ++it, ++i) {
				assignValueI2(pair, i);
				*it = pair;
			}
		}

		// 2. Sort suffix array with quicksort
		TSize q = length(shape);
		if (lengthSum(stringSet) == length(sa))
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				_QGramLess<typename Value<TSA>::Type, TStringSet const>(stringSet, q));
		else
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				_QGramLessNoCheck<typename Value<TSA>::Type, TStringSet const>(stringSet, q));
	}

	template < 
		typename TSA, 
		typename TDir,
		typename TText,
		typename TSize1,
		typename TSize2 >
	void _refineQGramIndex(
		TSA &sa,
		TDir &dir,
		TText const &text,
		TSize1 oldQ,
		TSize2 newQ)
	{
	SEQAN_CHECKPOINT
		typedef typename Size<TSA>::Type TSize;
		typedef typename Iterator<TSA, Standard>::Type		TIter;
		typedef typename Iterator<TDir, Standard>::Type		TDirIter;

		if (newQ <= (TSize2)oldQ) return;

		if (length(dir) < 2) {
			::std::sort(
				begin(sa, Standard()), 
				end(sa, Standard()), 
				_QGramLessOffset<typename Value<TSA>::Type, TText const>(text, newQ - oldQ, oldQ));
			return;
		}

		// 1. Sort each bucket with quicksort and compare substrings s[i+oldQ..i+newQ)
		TDirIter dirIt = begin(dir, Standard());
		TDirIter dirItEnd = end(dir, Standard());
		TIter itBegin = begin(sa, Standard());
		TIter itBktBegin = itBegin + *dirIt;
		TIter itBktEnd;
		++dirIt;
		for(; dirIt != dirItEnd; ++dirIt, itBktBegin = itBktEnd) {
			itBktEnd = itBegin + *dirIt;
			if (itBktEnd - itBktBegin < 2) continue;
			::std::sort(
				itBktBegin, 
				itBktEnd, 
				_QGramLessOffset<typename Value<TSA>::Type, TText const>(text, newQ - oldQ, oldQ));
		}
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.createQGramIndexDirOnly:
..summary:Builds the directory of a q-gram index on a sequence. 
..cat:Index
..signature:createQGramIndexDirOnly(sa, dir, text, shape)
..param.text:The sequence.
..param.shape:The shape to be used.
...type:Class.Shape
..param.sa:The resulting list in which all q-grams are sorted alphabetically.
..param.dir:The resulting array that indicates at which position in index the corresponding q-grams can be found.
..returns:Index contains the sorted list of qgrams. For each possible q-gram pos contains the first position in index that corresponds to this q-gram. 
*/

	template <
		typename TDir,
		typename TText,
		typename TShape >
	void createQGramIndexDirOnly(
		TDir &dir,
		TText const &text,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TIterator;
		typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
		typedef typename Value<TShape>::Type					TValue;
		typedef typename Size<TText>::Type						TSize;
		
		// 1. clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		// 2. count q-grams
		if (length(text) < length(shape)) return;
		TSize num_qgrams = length(text) - length(shape) + 1;

		TIterator itText = begin(text, Standard());
		++dir[hash(shape, itText)];
		for(TSize i = 1; i < num_qgrams; ++i)
		{
			++itText;
			++dir[hashNext(shape, itText)];
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
			SEQAN_ASSERT(sum + diff_prev == length(sa));
		}
	}

	template < 
		typename TDir,
		typename TString,
		typename TSpec,
		typename TShape >
	void createQGramIndexDirOnly(
		TDir &dir,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Iterator<TDir, Standard>::Type				TDirIterator;
		typedef typename Value<TShape>::Type						TValue;
		typedef typename Size<TString>::Type						TSize;
		
		// 1. clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		// 2. count q-grams
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = length(sequence) - length(shape) + 1;

			TIterator itText = begin(sequence, Standard());
			++dir[hash(shape, itText)];
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				++dir[hashNext(shape, itText)];
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
			SEQAN_ASSERT(sum + diff_prev == length(sa));
		}
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.createCountArray:
..summary:Builds an index on a StringSet storing how often a q-gram occurs in each sequence.
..cat:Index
..signature:createCountsArray(counts, dir, stringSet, shape)
..param.text:The StringSet.
...type:Class.StringSet
..param.shape:The shape to be used.
...type:Class.Shape
..param.count:The resulting list of pairs (seqNo,count).
..param.dir:The resulting array that indicates at which position in the count table the corresponding a certain q-gram can be found.
*/

	template < 
		typename TCounts, 
		typename TDir,
		typename TString,
		typename TSpec,
		typename TShape >
	void createCountsArray(
		TCounts &counts,
		TDir &dir,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Iterator<TDir, Standard>::Type				TDirIterator;
		typedef typename Value<TShape>::Type						TValue;
		typedef typename Size<TString>::Type						TSize;

		TDir lastSeqSeen;
		resize(lastSeqSeen, length(dir));
		
		// 1. clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);
		arrayFill(begin(lastSeqSeen, Standard()), end(lastSeqSeen, Standard()), -1);

		// 2. count distinct sequences for each q-gram
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = length(sequence) - length(shape) + 1;

			TIterator itText = begin(sequence, Standard());
			TValue hashValue = hash(shape, itText);
			lastSeqSeen[hashValue] = seqNo;
			++dir[hashValue];
			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				hashValue = hashNext(shape, itText);
				if (seqNo != lastSeqSeen[hashValue]) {
					lastSeqSeen[hashValue] = seqNo;
					++dir[hashValue];
				}
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
			resize(counts, sum + diff_prev);
		}

		// 4. fill count array
		arrayFill(begin(lastSeqSeen, Standard()), end(lastSeqSeen, Standard()), -1);
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);
			if (length(sequence) < length(shape)) continue;
			TSize num_qgrams = length(sequence) - length(shape) + 1;

			TIterator itText = begin(sequence, Standard());
			TValue hashValue = hash(shape, itText);
			lastSeqSeen[hashValue] = seqNo;
			TSize pos = dir[hashValue + 1]++;
			counts[pos].i1 = seqNo;				// first hash
			counts[pos].i2 = 1;

			for(TSize i = 1; i < num_qgrams; ++i)
			{
				++itText;
				hashValue = hashNext(shape, itText);
				if (seqNo == lastSeqSeen[hashValue])
					++(counts[dir[hashValue + 1] - 1].i2);
				else {
					lastSeqSeen[hashValue] = seqNo;
					pos = dir[hashValue + 1]++;
					counts[pos].i1 = seqNo;				// next hash
					counts[pos].i2 = 1;
				}
			}
		}
	}


//////////////////////////////////////////////////////////////////////////////
// create q-gram index of *one* sequence in external memory 

	// *** COMPARATORS & MAPS ***
        
    template <typename InType, typename Result = int>
    struct _qgram_comp : public ::std::binary_function<InType,InType,Result> {
        inline Result operator()(InType const &a, InType const &b) const
        {
			typedef typename Value<InType, 2>::Type TQGram;
			typename Value<TQGram>::Type const *sa = a.i2.i;
            typename Value<TQGram>::Type const *sb = b.i2.i;
			typename Value<TQGram>::Type const *saEnd = sa + length(a.i2);

            for(; sa != saEnd; ++sa, ++sb) {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }
			return posCompare(a.i1, b.i1);
        }
    };

    // optimized for bitvectors
    template <typename T1, typename TValue, unsigned _size, typename Result>
    struct _qgram_comp< Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >, Result > :
        public ::std::binary_function<
            Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >,
            Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >,
            Result> {       
        inline Result operator()(
            const Pair<T1, Tuple<TValue, _size, Compressed>, Compressed > &a,
            const Pair<T1, Tuple<TValue, _size, Compressed>, Compressed > &b) const
        {
            if (a.i2 < b.i2) return -1;
            if (a.i2 > b.i2) return 1;
			return posCompare(a.i1, b.i1);
        }
    };


    template <typename TValue, typename TResult = unsigned>
    struct _qgram_hash : public ::std::unary_function<TValue, TResult> {
        inline TResult operator()(TValue const &a) const
        {
			typedef typename Value<TValue, 2>::Type	TQGram;
			TResult hash = 0;
			unsigned len = length(a.i2);
            for(unsigned i = 0; i < len; ++i) {
				hash *= ValueSize< typename Value<TQGram>::Type >::VALUE;
				hash += (TResult)a.i2[i];
            }
            return hash;
        }
    };

	// TODO: replace fixed tuple size of 6 with q and add q to Shape template arguments
	template < 
		typename TSA, 
		typename TDir,
		typename TText,
		typename TShape >
	void createQGramIndexExt(
		TSA &suffixArray,
		TDir &dir,
		TText &text,
		TShape &shape)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename _MakeUnsigned< typename Value<TText>::Type >::Type TUValue;

        // *** SPECIALIZATION ***

		typedef Pipe< TText, Source<> >				TSource;
        typedef Pipe< TSource, Caster<TUValue> >    TUnsigner;
		typedef Pipe< TUnsigner, Tupler<7> >	    TTupler;
						                typedef _qgram_comp<_TypeOf(TTupler)> qcomp_t;
        typedef Pool< 
					_TypeOf(TTupler), 
					SorterSpec< SorterConfigSize<qcomp_t, _TSizeOf(TTupler) > > 
				> TSortTuples;
										typedef _qgram_hash<_TypeOf(TTupler), typename Size<TDir>::Type> qhash_t;

        // *** INSTANTIATION ***

		TSource			src(text);
        TUnsigner		unsigner(src);
		TTupler			tupler(unsigner);
		TSortTuples		sorter;

		// sort q-grams
		sorter << tupler;

		// fill sa and dir
		if (!beginRead(sorter)) return;

		typename Iterator<TSA>::Type itSA = begin(suffixArray);
		typename Iterator<TDir>::Type itDir = begin(dir);

		qcomp_t	qcomp;
		qhash_t qhash;

		typename Value<TSortTuples>::Type	old_qgram;
		typename Size<TDir>::Type			hash, old_hash = 0;
        typename Size<TSortTuples>::Type	leftToRead = length(sorter);
		bool first = true;

		while (leftToRead) {
			// copy occurence position
			*itSA = (*sorter).i1;
			if (first || qcomp(old_qgram, *sorter) != 0) {
				old_qgram = *sorter;
				hash = qhash(old_qgram);

				SEQAN_ASSERT(old_hash < hash);

				// copy bucket begin
				typename Size<TSortTuples>::Type i = length(sorter) - leftToRead;
				for(; old_hash < hash; ++old_hash, ++itDir)
					*itDir = i;
				first = false;
			}
			++src;
			--leftToRead;
		}

		// fill bucket table
		typename Size<TSortTuples>::Type i = length(sorter);
		hash = length(dir);
		for(; old_hash < hash; ++old_hash, ++itDir)
			*itDir = i;

		endRead(sorter);
	}


//////////////////////////////////////////////////////////////////////////////
// create q-gram index of *multiple* sequences in external memory 

	template < 
		typename TSA, 
		typename TDir,
		typename TString,
		typename TSpec,
		typename TShape,
		typename TLimitsString >
	void createQGramIndexExt(
		TSA &suffixArray,
		TDir &dir,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape,
		TLimitsString &limits)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type			TConcat;
        typedef typename _MakeUnsigned< typename Value<TConcat>::Type >::Type	TUValue;
		typedef Multi<
			Tupler<7, true, Compressed>, 
			typename Value<TSA>::Type,
			typename StringSetLimits< StringSet<TString, TSpec> >::Type >		TTuplerSpec;

        // *** SPECIALIZATION ***

		typedef Pipe< TConcat, Source<> >			TSource;
        typedef Pipe< TSource, Caster<TUValue> >    TUnsigner;
		typedef Pipe< TUnsigner, TTuplerSpec >	    TTupler;
						                typedef _qgram_comp<_TypeOf(TTupler)> qcomp_t;
        typedef Pool< 
					_TypeOf(TTupler), 
					SorterSpec< SorterConfigSize<qcomp_t, _TSizeOf(TTupler) > > 
				> TSortTuples;
										typedef _qgram_hash<_TypeOf(TTupler), typename Size<TDir>::Type> qhash_t;

        // *** INSTANTIATION ***

		TSource			src(concat(stringSet));
        TUnsigner		unsigner(src);
		TTupler			tupler(unsigner, limits);
		TSortTuples		sorter;

		// sort q-grams
		sorter << tupler;

		// fill sa and dir
		if (!beginRead(sorter)) return;

		typename Iterator<TSA>::Type itSA = begin(suffixArray);
		typename Iterator<TDir>::Type itDir = begin(dir);

		qcomp_t	qcomp;
		qhash_t qhash;

		typename Value<TSortTuples>::Type	old_qgram;
		typename Size<TDir>::Type			hash, old_hash = 0;
        typename Size<TSortTuples>::Type	leftToRead = length(sorter);
		bool first = true;

		while (leftToRead) {
			// copy occurence position
			*itSA = (*sorter).i1;
			if (first || qcomp(old_qgram, *sorter) != 0) {
				old_qgram = *sorter;
				hash = qhash(old_qgram);

				SEQAN_ASSERT(old_hash < hash);

				// copy bucket begin
				typename Size<TSortTuples>::Type i = length(sorter) - leftToRead;
				for(; old_hash < hash; ++old_hash, ++itDir)
					*itDir = i;
				first = false;
			}
			++src;
			--leftToRead;
		}

		// fill bucket table
		typename Size<TSortTuples>::Type i = length(sorter);
		hash = length(dir);
		for(; old_hash < hash; ++old_hash, ++itDir)
			*itDir = i;

		endRead(sorter);
	}



//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation 

	template <typename TText, typename TShapeSpec>
	inline bool indexCreate(Index<TText, Index_QGram<TShapeSpec> > &index, Tag<_Fibre_SADir> const, Default const) 
	{
		typedef Index<TText, Index_QGram<TShapeSpec> >			TIndex;
		typedef Shape<typename Value<TIndex>::Type, TShapeSpec>	TShape;

		TShape &shape = indexShape(index);

		// count all overlapping q-grams
		typename Size<TIndex>::Type qgram_count = 0;
		for(unsigned i = 0; i < countSequences(index); ++i)
			if (sequenceLength(i, index) >= length(shape))
				qgram_count += sequenceLength(i, index) - (length(shape) - 1);

		resize(indexSA(index), qgram_count, Exact());
		resize(indexDir(index), _fullDirLength(index), Exact());
		createQGramIndex(indexSA(index), indexDir(index), indexText(index), indexShape(index));
		return true;
	}

	template <typename TText, typename TSpec>
	inline bool indexSupplied(Index<TText, TSpec> &index, Tag<_Fibre_SADir> const) {
		return !(empty(getFibre(index, Tag<_Fibre_SA>())) || empty(getFibre(index, Tag<_Fibre_Dir>())));
	}

	template <typename TText, typename TShapeSpec>
	inline bool indexCreate(Index<TText, Index_QGram<TShapeSpec> > &index, Tag<_Fibre_SA> const, Default const alg)
	{
		typedef Index<TText, Index_QGram<TShapeSpec> >			TIndex;
		typedef Shape<typename Value<TIndex>::Type, TShapeSpec>	TShape;

		TShape &shape = indexShape(index);

		// count all overlapping q-grams
		typename Size<TIndex>::Type qgram_count = 0;
		for(unsigned i = 0; i < countSequences(index); ++i)
			if (sequenceLength(i, index) >= length(shape))
				qgram_count += sequenceLength(i, index) - (length(shape) - 1);

		resize(indexSA(index), qgram_count, Exact());
		createQGramIndexSAOnly(indexSA(index), indexText(index), indexShape(index));
		return true;
	}

	template <typename TText, typename TShapeSpec>
	inline bool indexCreate(Index<TText, Index_QGram<TShapeSpec> > &index, Tag<_Fibre_Counts> const, Default const) 
	{
		resize(indexCountsDir(index), _fullDirLength(index), Exact());
		createCountsArray(indexCounts(index), indexCountsDir(index), indexText(index), indexShape(index));
		return true;
	}

/*
	template <typename TText, typename TShapeSpec>
	inline bool indexCreate(Index<TText, Index_QGram<TShapeSpec> > &index, Tag<_Fibre_Dir> const, Default const alg)
	{
		typedef Index<TText, Index_QGram<TShapeSpec> >			TIndex;

		resize(indexDir(index), _fullDirLength(index), Exact());
		createQGramIndexDirOnly(indexDir(index), indexText(index), indexShape(index));
		return true;
	}
*/

//////////////////////////////////////////////////////////////////////////////
// getKmerSimilarityMatrix
//     for the whole StringSet

	template < typename TObject, typename TSpec, typename TDistMatrix >
	inline void getKmerSimilarityMatrix(
		Index< TObject, Index_QGram<TSpec> > &index, 
		TDistMatrix &distMat)
	{
		typedef Index< TObject, Index_QGram<TSpec> >			TIndex;
		typedef typename Size<TIndex>::Type						TSize;
		typedef typename Size<TDistMatrix>::Type				TSizeMat;
		typedef typename Value<TDistMatrix>::Type				TValueMat;

		typedef typename Fibre<TIndex, QGram_CountsDir>::Type	TCountsDir;
		typedef typename Iterator<TCountsDir, Standard>::Type	TIterCountsDir;
		typedef typename Fibre<TIndex, QGram_Counts>::Type		TCounts;
		typedef typename Iterator<TCounts, Standard>::Type		TIterCounts;

		// declare requirements
		indexRequire(index, QGram_Counts());

		// initialize distance matrix
		TSizeMat seqNoLength = countSequences(index);
		clear(distMat);
		resize(distMat, seqNoLength * seqNoLength);
		arrayFill(begin(distMat, Standard()), end(distMat, Standard()), 0);

		TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
		TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
		TIterCounts itCountsBegin = begin(indexCounts(index), Standard());

		// for each bucket count common q-grams for each sequence pair
		TSize bucketBegin = *itCountsDir;
		for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir) 
		{
			TSize bucketEnd = *itCountsDir;

			// q-gram must occur in at least 2 different sequences
			if (bucketBegin != bucketEnd) 
			{
				TIterCounts itA = itCountsBegin + bucketBegin;
				TIterCounts itEnd = itCountsBegin + bucketEnd;
				for(; itA != itEnd; ++itA) 
				{
					TSizeMat ofs = (*itA).i1 * seqNoLength;
					TSize countA = (*itA).i2;
					TIterCounts itB = itA;

					for(; itB != itEnd; ++itB) 
					{
						TSize countB = (*itB).i2;
						if (countA < countB)
							distMat[ofs + (*itB).i1] += countA;
						else
							distMat[ofs + (*itB).i1] += countB;
					}
				}
			}
			bucketBegin = bucketEnd;
		}

		// copy upper triangle to lower triangle and scale
		for(TSizeMat row = 0; row < seqNoLength; ++row) 
		{
			TValueMat maxValRow = distMat[row * (seqNoLength + 1)];
			for(TSizeMat col = row + 1; col < seqNoLength; ++col)
			{
				// fractional common kmer count
				TValueMat maxValCol = distMat[col * (seqNoLength + 1)];
				TValueMat val = distMat[row * seqNoLength + col];

				// number of common q-grams / Number of possible common q-grams
				if (maxValRow < maxValCol) {
					if (maxValRow != 0)
						val /= maxValRow;
					distMat[col * seqNoLength + row] = val;
					distMat[row * seqNoLength + col] = val;
				} else {
					if (maxValCol != 0)
						val /= maxValCol;
					distMat[col * seqNoLength + row] = val;
					distMat[row * seqNoLength + col] = val;
				}
			}
		}

		// set diagonal to 1
		for(TSizeMat i = 0; i < seqNoLength; ++i)
			distMat[i * (seqNoLength + 1)] = 1;
	}


//////////////////////////////////////////////////////////////////////////////
// getKmerSimilarityMatrix 
//     for a subset of the StringSet given as a sorted string of sequence numbers

	template < typename TObject, typename TSpec, typename TDistMatrix, typename TSeqNoString >
	inline void getKmerSimilarityMatrix(
		Index< TObject, Index_QGram<TSpec> > &index, 
		TDistMatrix &distMat,
		TSeqNoString const &seqNo)
	{
		typedef Index< TObject, Index_QGram<TSpec> >			TIndex;
		typedef typename Size<TIndex>::Type						TSize;
		typedef typename Size<TDistMatrix>::Type				TSizeMat;
		typedef typename Value<TDistMatrix>::Type				TValueMat;

		typedef typename Fibre<TIndex, QGram_CountsDir>::Type	TCountsDir;
		typedef typename Iterator<TCountsDir, Standard>::Type	TIterCountsDir;
		typedef typename Fibre<TIndex, QGram_Counts>::Type		TCounts;
		typedef typename Iterator<TCounts, Standard>::Type		TIterCounts;
		typedef typename Iterator<TSeqNoString, Standard>::Type	TIterSeqNo;

		// declare requirements
		indexRequire(index, QGram_Counts());

		// initialize distance matrix
		TSizeMat seqNoLength = length(seqNo);
		clear(distMat);
		resize(distMat, seqNoLength * seqNoLength);
		arrayFill(begin(distMat, Standard()), end(distMat, Standard()), 0);

		TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
		TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
		TIterCounts itCountsBegin = begin(indexCounts(index), Standard());
		TIterSeqNo itSetBegin = begin(seqNo, Standard());
		TIterSeqNo itSetEnd = end(seqNo, Standard());

		// for each bucket count common q-grams for each sequence pair
		TSize bucketBegin = *itCountsDir;
		for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir) 
		{
			TSize bucketEnd = *itCountsDir;

			// q-gram must occur in at least 2 different sequences
			if (bucketBegin != bucketEnd) 
			{
				TIterCounts itA = itCountsBegin + bucketBegin;
				TIterCounts itEnd = itCountsBegin + bucketEnd;
				TIterSeqNo itSetA = itSetBegin;

				while (itA != itEnd && itSetA != itSetEnd)
				{
					if ((*itA).i1 < *itSetA)
						++itA;
					else if ((*itA).i1 > *itSetA)
						++itSetA;
					else 
					{
						TSizeMat ofs = (*itA).i1 * seqNoLength;
						TSize countA = (*itA).i2;
						TIterCounts itB = itA;
						TIterSeqNo itSetB = itSetA;

						while (itB != itEnd && itSetB != itSetEnd)
						{
							if ((*itB).i1 < *itSetB)
								++itB;
							else if ((*itB).i1 > *itSetB)
								++itSetB;
							else 
							{
								TSize countB = (*itB).i2;
								if (countA < countB)
									distMat[ofs + (*itB).i1] += countA;
								else
									distMat[ofs + (*itB).i1] += countB;
								++itB;
								++itSetB;
							}
						}
						++itA;
						++itSetA;
					}
				}
			}
			bucketBegin = bucketEnd;
		}

		// copy upper triangle to lower triangle and scale
		for(TSizeMat row = 0; row < seqNoLength; ++row) 
		{
			TValueMat maxValRow = distMat[row * (seqNoLength + 1)];
			for(TSizeMat col = row + 1; col < seqNoLength; ++col)
			{
				// fractional common kmer count
				TValueMat maxValCol = distMat[col * (seqNoLength + 1)];
				TValueMat val = distMat[row * seqNoLength + col];

				// number of common q-grams / Number of possible common q-grams
				if (maxValRow < maxValCol) {
					if (maxValRow != 0)
						val /= maxValRow;
					distMat[col * seqNoLength + row] = val;
					distMat[row * seqNoLength + col] = val;
				} else {
					if (maxValCol != 0)
						val /= maxValCol;
					distMat[col * seqNoLength + row] = val;
					distMat[row * seqNoLength + col] = val;
				}
			}
		}

		// set diagonal to 1
		for(TSizeMat i = 0; i < seqNoLength; ++i)
			distMat[i * (seqNoLength + 1)] = 1;
	}


//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_QGram<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!open(getFibre(index, QGram_Text()), toCString(name), openMode)) && 
			(!open(getFibre(index, QGram_Text()), fileName), openMode)) return false;

		name = fileName;	append(name, ".sa");	open(getFibre(index, QGram_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".dir");	open(getFibre(index, QGram_Dir()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_QGram<TSpec> > &index, 
		const char *fileName) 
	{
		return open(index, fileName, OPEN_RDONLY);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_QGram<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!save(getFibre(index, QGram_Text()), toCString(name), openMode)) && 
			(!save(getFibre(index, QGram_Text()), fileName), openMode)) return false;

		name = fileName;	append(name, ".sa");	save(getFibre(index, QGram_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".dir");	save(getFibre(index, QGram_Dir()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_QGram<TSpec> > &index, 
		const char *fileName)
	{
		return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

}

#endif //#ifndef SEQAN_HEADER_...
