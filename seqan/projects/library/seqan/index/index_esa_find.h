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

#ifndef SEQAN_HEADER_INDEX_FIND_H
#define SEQAN_HEADER_INDEX_FIND_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct SortedList {};
	struct LeftCompleteTree {};

	template < unsigned BlockSize = 4096 >
	struct BTree {};

	template < typename TString, typename TSpec >
	class SearchTreeIterator {};

	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search tree like a real tree
	//

	template <typename TString>
	class SearchTreeIterator< TString, SortedList >
	{
	public:

		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		inline SearchTreeIterator(TString &string):
			_First(begin(string, Standard())),
			_Count(length(string))
		{
			_Count2 = _Count / 2;
			_Mid = _First;
			goFurther(_Mid, _Count2);
		}
			
        inline const TValue& operator*() const {
			return *_Mid;
		}

        inline const TValue* operator->() const {
			return &*_Mid;
		}

		inline TSize mid() {
			return _Count2;
		}

		// descend left
		inline SearchTreeIterator & left()
		{
			_Count = _Count2;
			_Count2 /= 2;
			_Mid = _First;
			goFurther(_Mid, _Count2);
			return *this;
		}

		// descend right
        inline SearchTreeIterator & right()
		{
			_First = ++_Mid, _Count -= _Count2 + 1;
			_Count2 = _Count / 2;
			goFurther(_Mid, _Count2);
			return *this;
		}

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline bool eof() {
            return !_Count;
        }

		inline operator TIterator & () {
			return _Mid;
		}

	private:
		TIterator	_First, _Mid;
		TSize		_Count, _Count2;
	};


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search tree like a real tree
	//

	template <typename TString>
	class SearchTreeIterator< TString, LeftCompleteTree > 
	{
	public:
		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		TIterator	it;
		TSize		size;

		inline SearchTreeIterator(TString &string):
			it(begin(string, Standard())),
			size(length(string))
		{
			_left = 0;
			_lSize = 1;
            _iSize = size;
			for(_xSize = 1; _xSize < size; _xSize <<= 1);
            if (!size) _xSize = 0;
		}

		inline SearchTreeIterator():
			it(),
            size(0),
            _xSize(0) {}

        inline const TValue& operator*() const {
			return *it;
		}

        inline const TValue* operator->() const {
			return &*it;
		}

		inline TSize mid() {
			return _xSize >> 1;
		}

		inline TSize leftSize() {
			return mid();
		}

		inline TSize rightSize() {
			return _iSize - mid();
		}

		// descend left
        inline SearchTreeIterator & left()
		{
            if (_xSize <= 1) {
                _xSize = 0;
                return *this;
            }
            _descendLeft();
			_iSize = _xSize;	    // = mid();
            return *this;
        }

		// descend right
        inline SearchTreeIterator & right()
		{
            if (_xSize <= 1) {
                _xSize = 0;
                return *this;
            }
			_iSize -= mid();
            SEQAN_ASSERT(_iSize != 0);    // _xSize/2 is less than _iSize by invariant

            // step down right
            _descendRight();

            // step down left until we have two sons or are a leaf
            while (_iSize <= (_xSize >> 1))
                _descendLeft();

            return *this;
        }

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline SearchTreeIterator& operator--() {
            --it;
            return *this;
        }

        inline SearchTreeIterator operator--(int) {
            SearchTreeIterator before = *this;
            --it;
            return before;
        }

        inline SearchTreeIterator & operator++() {
            ++it;
            return *this;
        }

        inline SearchTreeIterator operator++(int) {
            SearchTreeIterator before = *this;
            ++it;
            return before;
        }

        inline bool operator==(SearchTreeIterator const &I) {
            return (_xSize == I._xSize) && (_xSize == 0 || it == I.it);
        }

        //operator FlatFwdIt() {
        //    return it;
        //}

        inline bool eof() {
            return !_xSize;
        }

		inline operator TIterator & () {
			return it;
		}

	private:
		TSize _left;		// left iterator offset of current interval
		TSize _lSize;		// iterator elements of current level
		TSize _xSize;		// max interval size of current level
		TSize _iSize;		// current interval size

        inline void _descendLeft() 
		{
			goFurther(it, _left + _lSize);
			_left <<= 1;
			_xSize >>= 1;
			_lSize = (size + _xSize - 1) / _xSize;
        }

        inline void _descendRight() 
		{
			goFurther(it, _left + 1 + _lSize);
			_left = (_left << 1) + 1;
			_xSize >>= 1;
			_lSize = (size + _xSize - 1) / _xSize;
        }
    };


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search b-tree like a real tree
	//

	template <typename TString, unsigned BlockSize>
	class SearchTreeIterator< TString, BTree<BlockSize> > 
	{
	public:
		typedef typename Value<TString>::Type				TValue;
		typedef typename Size<TString>::Type				TSize;
		typedef typename Iterator<TString, Standard>::Type	TIterator;

		enum { BlockHeight = Log2Floor<BlockSize>::VALUE };
		enum { BlockElements = (1 << BlockHeight) - 1 };
		enum { BlockInnerElements = (1 << (BlockHeight - 1)) - 1 };
		enum { BlockLeafs = 1 << (BlockHeight - 1) };

		TIterator	it;
		TSize		size;

		inline SearchTreeIterator(TString &string):
			it(begin(string, Standard())),
			size(length(string))
		{
			//_left = 0;
			//_lSize = 0;
   //         _iSize = size;

			_heightHigh = 1;
			_stepSizeLow = 1;
			for(TSize _xSizeLow = 2; _xSizeLow <= size; _xSizeLow <<= 1) {
				if (_stepSizeLow == BlockLeafs) {
					_stepSizeLow = 1;
					++_heightHigh;
				} else
					_stepSizeLow <<= 1;
			}
			
			_stepSizeLow >>= 1;
			for(_xSizeHigh = 1; _xSizeHigh * BlockSize <= size; _xSizeHigh *= BlockSize);

			_leftLow = (_stepSizeLow << 1) - 1;		// point to the middle
			_leftHigh = 0;
			_lSizeHigh = 1;

			it += _leftLow;

			_updateElements();

			//if (!size) _xSizeLow = 0;
   //         if (!size) _xSizeHigh = 0;
		}

		inline SearchTreeIterator():
			it() {}

        inline const TValue& operator*() const {
			return *it;
		}

		inline const TValue* operator->() const {
			return &*it;
		}

		// descend left
        inline SearchTreeIterator & left() 
		{
            if (_heightHigh == 1 && !_stepSizeLow) {
                _heightHigh = 0;
                return *this;
            }
            _descendLeft();
            return *this;
        }

		// descend right
        inline SearchTreeIterator & right() 
		{
            if (_heightHigh == 1 && !_stepSizeLow) {
				++it;
                _heightHigh = 0;
                return *this;
            }

            // step down right
            _descendRight();

            // step down left until we have two sons or are a leaf
            while (_elements <= _leftLow)
                _descendLeft();

            return *this;
        }

        inline SearchTreeIterator leftChild() const {
            return SearchTreeIterator(*this).left();
        }

        inline SearchTreeIterator rightChild() const {
            return SearchTreeIterator(*this).right();
        }

        inline SearchTreeIterator& operator--() {
            --it;
            return *this;
        }

        inline SearchTreeIterator operator--(int) {
            SearchTreeIterator before = *this;
            --it;
            return before;
        }

        inline SearchTreeIterator& operator++() {
            ++it;
            return *this;
        }

        inline SearchTreeIterator operator++(int) {
            SearchTreeIterator before = *this;
            ++it;
            return before;
        }

        inline bool operator==(SearchTreeIterator const &I) {
			return it == I.it || (eof() && I.eof());
        }

        //operator FlatFwdIt() {
        //    return it;
        //}

        inline bool eof() {
            return !_heightHigh;
        }

		inline operator TIterator & () {
			return it;
		}

	private:
		TSize _heightHigh;	// height measured in BBlocks
		unsigned _elements;		// elements in current BBlock

		unsigned _leftLow;		// left iterator offset of current interval
		unsigned _stepSizeLow;	// left/right step size of current level

		TSize _leftHigh;		// left BBlock offset of current interval
		TSize _lSizeHigh;	// BBlocks of current level
		TSize _xSizeHigh;	// max BBlocks of current level

		inline void _descendLeft() 
		{
			if (_stepSizeLow) {
				it -= _stepSizeLow;
				_leftLow -= _stepSizeLow;
				_stepSizeLow >>= 1;
			} else
				if (--_heightHigh) {
					_descendBBlock(_leftLow);

					_leftLow = BlockInnerElements;		// point to the middle
					_stepSizeLow = BlockLeafs / 2;		// with step size of half an interval length
				}
        }

		inline void _descendRight() 
		{
			if (_stepSizeLow) {
				it += _stepSizeLow;
				_leftLow += _stepSizeLow;
				_stepSizeLow >>= 1;
			} else
				if (--_heightHigh) {
					_descendBBlock(_leftLow + 1);

					_leftLow = BlockInnerElements;		// point to the middle
					_stepSizeLow = BlockLeafs / 2;		// with step size of half an interval length
				}
        }

		inline void _updateElements() 
		{
			TSize firstElement = (1 + _leftHigh * BlockSize) * _xSizeHigh - 1;
			TSize lastElement = (1 + (_leftHigh + 1) * BlockSize) * _xSizeHigh - 2;

			if (lastElement >= size)
				_elements = (size - firstElement) / _xSizeHigh;
			else
				_elements = BlockElements;
		}

        inline void _descendBBlock(TSize _childIndex) 
		{
			// goFurther to the begin of the current BBlock and further to the destination BBlock
			goFurther(it, BlockSize * (_leftHigh * (BlockSize - 1) + _childIndex + _lSizeHigh) + BlockInnerElements - _leftLow);

			_leftHigh = _leftHigh * BlockSize + _childIndex;
			_xSizeHigh /= BlockSize;
			_lSizeHigh = (size / _xSizeHigh + BlockSize - 1) / BlockSize;

			_updateElements();
        }
    };


	//////////////////////////////////////////////////////////////////////////////
	// substring search with SA table and w/o LCP-table or enhancement
	//

	template <
		typename TText,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	_lowerBoundSA(
		TText &text,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query)
	{	// find first element not before query, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;

		TDiff lcpLower = 0;
		TDiff lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.eof(); )
		{	// divide and conquer, find half that contains answer

			TSuffix		suf = suffix(text, *treeIter);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TDiff		lcp = Min(lcpLower, lcpUpper);

			goFurther(t, lcp);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && *t == *q) {
				++t;
				++q;
				++lcp;
			}
			
            // is text < query ?
			if (q != qEnd && (t == tEnd || *t < *q)) {
				treeIter.right();
				lcpLower = lcp;
			} else {
				treeIter.left();
				lcpUpper = lcp;
			}
		}
		return treeIter;
	}

	//////////////////////////////////////////////////////////////////////////////

	template <
		typename TText,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline typename Iterator<TSA, Standard>::Type
	_upperBoundSA(
		TText &text,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query)
	{	// find first element that query is before, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;

		TDiff lcpLower = 0;
		TDiff lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.eof(); )
		{	// divide and conquer, find half that contains answer

			TSuffix		suf = suffix(text, *treeIter);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TDiff		lcp = Min(lcpLower, lcpUpper);

			goFurther(t, lcp);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && *t == *q) {
				++t;
				++q;
				++lcp;
			}
			
            // is text <= query ?
			if (q == qEnd || t == tEnd || !(*q < *t)) {
				treeIter.right();
				lcpLower = lcp;
			} else {
				treeIter.left();
				lcpUpper = lcp;
			}
		}
		return treeIter;
	}

	//////////////////////////////////////////////////////////////////////////////

	template <
		typename TText,
		typename TSA,
		typename TSpec,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA, Standard>::Type >
	_equalRangeSA(
		TText &text,
		SearchTreeIterator< TSA, TSpec > treeIter,
		TQuery &query)
	{	// find range equivalent to query, using operator<
		typedef typename Difference<TText>::Type			TDiff;
		typedef typename Suffix<TText>::Type				TSuffix;
		typedef typename Iterator<TSuffix, Standard>::Type	TTextIter;
		typedef typename Iterator<TQuery, Standard>::Type	TQueryIter;
		typedef typename Iterator<TSA, Standard>::Type		TSAIter;

		TDiff lcpLower = 0;
		TDiff lcpUpper = 0;

		TQueryIter qBegin = begin(query, Standard());
		TQueryIter qEnd = end(query, Standard());

		for (; !treeIter.eof(); )
		{	// divide and conquer, check midpoint

			TSuffix		suf = suffix(text, *treeIter);
			TTextIter	t = begin(suf, Standard());
			TTextIter	tEnd = end(suf, Standard());
			TQueryIter	q = qBegin;
            TDiff		lcp = Min(lcpLower, lcpUpper);

			goFurther(t, lcp);
			goFurther(q, lcp);
			while (t != tEnd && q != qEnd && *t == *q) {
				++t;
				++q;
				++lcp;
			}
			
            // is text < query ?
			if (q != qEnd && (t == tEnd || *t < *q))
			{	// range begins above _Mid, loop
				treeIter.right();
				lcpLower = lcp;
			}
            // is text > query ?
			else if (q != qEnd && (t != tEnd && *q < *t))
			{	// range in first half, loop
				treeIter.left();
				lcpUpper = lcp;
			} else
            // is text == query ?
			{	// range straddles mid, find each end and return
				return Pair<TSAIter> (
					_lowerBoundSA(text, treeIter.leftChild(), query),
					_upperBoundSA(text, treeIter.rightChild(), query)
				);
			}
		}
		return Pair<TSAIter> (treeIter, treeIter);	// empty range
	}


	//////////////////////////////////////////////////////////////////////////////
	// Finder wrappers (return iterators instead of positions)

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	lowerBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	upperBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA const, Standard>::Type >
	equalRangeSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find range equivalent to query, using operator<
		return _equalRangeSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}


	//////////////////////////////////////////////////////////////////////////////
	// workarounds for the Visual Studio array problem

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	lowerBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Iterator<TSA const, Standard>::Type
	upperBoundSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Iterator<TSA const, Standard>::Type >
	equalRangeSAIterator(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find range equivalent to query, using operator<
		return _equalRangeSA(text, SearchTreeIterator<TSA const, SortedList>(sa), query);
	}


	//////////////////////////////////////////////////////////////////////////////
	// wrappers with nice interfaces

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find first element not before query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element not before query, using operator<
		return lowerBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find first element that query is before, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element that query is before, using operator<
		return upperBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query,
		TFlatTreeSpec const)
	{	// find range equivalent to query, using operator<
		Pair< typename Iterator<TSA, Standard>::Type > itPair = 
			_equalRangeSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query);
		return Pair< typename Position<TSA>::Type >
			(itPair.i1 - begin(sa, Standard()), itPair.i2 - begin(sa, Standard()));
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery const &query)
	{	// find first element that query is before, using operator<
		return equalRangeSA(text, sa, query, SortedList());
	}


	//////////////////////////////////////////////////////////////////////////////
	// workarounds for the Visual Studio array problem

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find first element not before query, using operator<
		return _lowerBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element not before query, using operator<
		return lowerBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find first element that query is before, using operator<
		return _upperBoundSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query) - begin(sa, Standard());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline typename Position<TSA>::Type 
	upperBoundSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element that query is before, using operator<
		return upperBoundSA(text, sa, query, SortedList());
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery,
		typename TFlatTreeSpec
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery *query,
		TFlatTreeSpec const)
	{	// find range equivalent to query, using operator<
		Pair< typename Iterator<TSA, Standard>::Type > itPair = 
			_equalRangeSA(text, SearchTreeIterator<TSA const, TFlatTreeSpec>(sa), query);
		return Pair< typename Position<TSA>::Type >
			(itPair.i1 - begin(sa, Standard()), itPair.i2 - begin(sa, Standard()));
	}

	template <
		typename TText,
		typename TSA,
		typename TQuery
	>
	inline Pair< typename Position<TSA>::Type >
	equalRangeSA(
		TText const &text,
		TSA const &sa,
		TQuery *query)
	{	// find first element that query is before, using operator<
		return equalRangeSA(text, sa, query, SortedList());
	}

/*
	//////////////////////////////////////////////////////////////////////////////
	// substring search with enhanced LCP-table
	//

	template <
		typename TTextIter,
		typename TSAIter,
		typename LCPTreeIt,
		typename TQueryIter
	>
	inline TSAIter _Lower_bound_lcp_enhanced(
		TTextIter tBegin,
		TTextIter tEnd,
		TSAIter _First,
		TSAIter _Last,
		LCPTreeIt _LCPTop,
		TQueryIter qBegin,
		TQueryIter qEnd,
		typename Difference<TTextIter>::Type lcpLower,
		typename Difference<TTextIter>::Type lcpUpper)
	{	// find first element not before query, using operator<
		typedef typename Difference<TTextIter>::Type TDiff;
        TDiff delta = difference(_First, _Last) - 1;
		TDiff lcp;
		#ifdef SEQAN_PROFILE_LCPEFIND
			TDiff skippedCompares = 0;	// difference of char compares related to xxx_bound_sa
		#endif

        // binary search with intervals >= 3 elements
		for (; 1 < delta; )
		{	// divide and conquer, find half that contains answer
			TDiff _Delta2 = _LCPTop.leftSize();
			TSAIter _Mid = _First;
			goFurther(_Mid, _Delta2);

			if (lcpLower > lcpUpper) {
                LCPTreeIt leftChild = _LCPTop;
                leftChild.left();
				TDiff _lcpMidLower = *leftChild;

				if (_lcpMidLower > lcpLower) {
					// second half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += lcpLower - lcpUpper;
					#endif
					_First = _Mid;
					_LCPTop.right();
					if ((delta -= _Delta2) == 1) {
						++_First;
						delta = 0;
					}
					continue;
				} else	if (_lcpMidLower < lcpLower) {
					// first half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += _lcpMidLower - lcpUpper;
					#endif
					lcpUpper = _lcpMidLower;
					_LCPTop = leftChild;
					if ((delta = _Delta2) == 1)
						delta = 0;
					continue;
				}
				lcp = lcpLower;
			} else if (lcpLower < lcpUpper) {
                LCPTreeIt rightChild = _LCPTop;
                rightChild.right();
				TDiff _lcpMidUpper = *rightChild;

				if (_lcpMidUpper > lcpUpper) {
					// first half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += lcpUpper - lcpLower;
					#endif
					_LCPTop.left();
					if ((delta = _Delta2) == 1)
						delta = 0;
					continue;
				} else if (_lcpMidUpper < lcpUpper) {
					// second half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += _lcpMidUpper - lcpLower;
					#endif
					lcpLower = _lcpMidUpper;
					_First = _Mid;
                    _LCPTop = rightChild;
					if ((delta -= _Delta2) == 1) {
						++_First;
						delta = 0;
					}
					continue;
				}
				lcp = lcpUpper;
			} else
				lcp = lcpUpper;

			TTextIter t = tBegin;
			TQueryIter q = qBegin;
            goFurther(t, *_Mid);

			// lcp search changes MIN to MAX here
//			TDiff lcp = Max(lcpLower, lcpUpper);
			#ifdef SEQAN_PROFILE_LCPEFIND
				skippedCompares += lcp - Min(lcpLower, lcpUpper);
			#endif
			goFurther(t, lcp);
			goFurther(q, lcp);
			for(TDiff i = Min(difference(t, tEnd), difference(q, qEnd)); 
				i && *t == *q;
				--i, ++t, ++q, ++lcp);

            // is text < query ?
			if (q != qEnd && (t == tEnd || *t < *q)) {
				// second half
				lcpLower = lcp;
				_First = _Mid;
				_LCPTop.right();
				if ((delta -= _Delta2) == 1) {
					++_First;
					delta = 0;
				}
			} else {
				// first half
				lcpUpper = lcp;
				_LCPTop.left();
				if ((delta = _Delta2) == 1)
					delta = 0;
			}
		}

		#ifdef SEQAN_PROFILE_LCPEFIND
			SEQAN_PROADD(SEQAN_PROEXTRA3, skippedCompares);
		#endif

        // binary search for intervals of 2 or less elements
        lcp = Min(lcpLower, lcpUpper);

        TQueryIter q = qBegin;
		goFurther(q, lcp);

        while (true) {
			TTextIter t = tBegin;
			goFurther(t, *_First + lcp);

			for(TDiff i = Min(difference(t, tEnd), difference(q, qEnd)); 
            	i && *t == *q;
            	 --i, ++t, ++q);

            // is text < query ?
			if (q != qEnd && (t == tEnd || *t < *q)) {
				// second half
				++_First;
				if (!delta) return _First;
                --delta;
			} else {
				// first half -> end
				return _First;
			}
        }
	}

	template <
		typename TTextIter,
		typename TSAIter,
		typename LCPFwdIt,
		typename TQueryIter
	>
	inline TSAIter lower_bound(
		TTextIter tBegin,
		TTextIter tEnd,
		TSAIter _First,
		TSAIter _Last,
		SearchTreeIterator<LCPFwdIt, LeftCompleteTree> _LCPTop,
		TQueryIter qBegin,
		TQueryIter qEnd)
	{	// find first element not before query, using operator<
		return _Lower_bound_lcp_enhanced(
			tBegin, tEnd,
			_First, _Last,
			_LCPTop,
			qBegin, qEnd,
			0, 0);
	}


	//////////////////////////////////////////////////////////////////////////////

	template <
		typename TTextIter,
		typename TSAIter,
		typename LCPTreeIt,
		typename TQueryIter
	>
	inline TSAIter _Upper_bound_lcp_enhanced(
		TTextIter tBegin,
		TTextIter tEnd,
		TSAIter _First,
		TSAIter _Last,
		LCPTreeIt _LCPTop,
		TQueryIter qBegin,
		TQueryIter qEnd,
		typename Difference<TTextIter>::Type lcpLower,
		typename Difference<TTextIter>::Type lcpUpper)
	{	// find first element not before query, using operator<
		typedef typename Difference<TTextIter>::Type TDiff;
        TDiff delta = difference(_First, _Last) - 1;

        // binaray search with intervals >= 3 elements
		for (; 1 < delta; )
		{	// divide and conquer, find half that contains answer
			TDiff _Delta2 = _LCPTop.leftSize();
			TSAIter _Mid = _First;
			goFurther(_Mid, _Delta2);

			if (lcpLower > lcpUpper) {
                LCPTreeIt leftChild = _LCPTop;
                leftChild.left();
				TDiff _lcpMidLower = *leftChild;

				if (_lcpMidLower > lcpLower) {
					// second half
					_First = _Mid;
					_LCPTop.right();
					if ((delta -= _Delta2) == 1) {
						++_First;
						delta = 0;
					}
					continue;
				} else	if (_lcpMidLower < lcpLower) {
					// first half
					lcpUpper = _lcpMidLower;
					_LCPTop = leftChild;
					if ((delta = _Delta2) == 1)
						delta = 0;
					continue;
				}
			} else if (lcpLower < lcpUpper) {
                LCPTreeIt rightChild = _LCPTop;
                rightChild.right();
				TDiff _lcpMidUpper = *rightChild;

				if (_lcpMidUpper > lcpUpper) {
					// first half
					_LCPTop.left();
					if ((delta = _Delta2) == 1)
						delta = 0;
					continue;
				} else if (_lcpMidUpper < lcpUpper) {
					// second half
					lcpLower = _lcpMidUpper;
					_First = _Mid;
                    _LCPTop = rightChild;
					if ((delta -= _Delta2) == 1) {
						++_First;
						delta = 0;
					}
					continue;
				}
			}

			TTextIter t = tBegin;
			TQueryIter q = qBegin;
            goFurther(t, *_Mid);

			// lcp search changes MIN to MAX here
			TDiff lcp = Max(lcpLower, lcpUpper);
			goFurther(t, lcp);
			goFurther(q, lcp);
			TDiff max = Min(difference(t, tEnd), difference(q, qEnd));
            TDiff i = max;
			for(; i && *t == *q; --i, ++t, ++q);
			lcp += max - i;

            // is text <= query ?
			if (q == qEnd || t == tEnd || !(*q < *t)) {
				// second half
				lcpLower = lcp;
				_First = _Mid;
				_LCPTop.right();
				if ((delta -= _Delta2) == 1) {
					++_First;
					delta = 0;
				}
			} else {
				// first half
				lcpUpper = lcp;
				_LCPTop.left();
				if ((delta = _Delta2) == 1)
					delta = 0;
			}
		}

        // binary search for intervals of 2 or less elements
        TDiff lcp = Min(lcpLower, lcpUpper);

        TQueryIter q = qBegin;
		goFurther(q, lcp);

        while (true) {
			TTextIter t = tBegin;
			goFurther(t, *_First + lcp);

			TDiff i = Min(difference(t, tEnd), difference(q, qEnd));
            for(; i && *t == *q; --i, ++t, ++q);

            // is text <= query ?
			if (q == qEnd || t == tEnd || !(*q < *t)) {
				// second half
				++_First;
				if (!delta) return _First;
                --delta;
			} else {
				// first half -> end
				return _First;
			}
        }
	}

	template <
		typename TTextIter,
		typename TSAIter,
		typename LCPFwdIt,
		typename TQueryIter
	>
	inline TSAIter upper_bound(
		TTextIter tBegin,
		TTextIter tEnd,
		TSAIter _First,
		TSAIter _Last,
		SearchTreeIterator<LCPFwdIt, LeftCompleteTree> _LCPTop,
		TQueryIter qBegin,
		TQueryIter qEnd)
	{	// find first element not before query, using operator<
		return _Upper_bound_lcp_enhanced(
			tBegin, tEnd,
			_First, _Last,
			_LCPTop,
			qBegin, qEnd,
			0, 0);
	}


	//////////////////////////////////////////////////////////////////////////////

	template <
		typename TTextIter,
		typename TSAIter,
		typename LCPTreeIt,
		typename TQueryIter
	>
	inline Pair<TSAIter> _Equal_range_lcp_enhanced(
		TTextIter tBegin,
		TTextIter tEnd,
		TSAIter _First,
		TSAIter _Last,
		LCPTreeIt _LCPTop,
		TQueryIter qBegin,
		TQueryIter qEnd)
	{	// find first element not before query, using operator<
		typedef typename Difference<TTextIter>::Type TDiff;
		TDiff lcpLower = 0;
		TDiff lcpUpper = 0;
        TDiff delta = difference(_First, _Last) - 1;

        // binaray search with intervals >= 3 elements
		for (; 1 < delta; )
		{	// divide and conquer, find half that contains answer
			TDiff _Delta2 = _LCPTop.leftSize();
			TSAIter _Mid = _First;
			goFurther(_Mid, _Delta2);

			if (lcpLower > lcpUpper) {
                LCPTreeIt leftChild = _LCPTop;
                leftChild.left();
				TDiff _lcpMidLower = *leftChild;

				if (_lcpMidLower > lcpLower) {
					// second half
					_First = _Mid;
					_LCPTop.right();
					if ((delta -= _Delta2) == 1) {
						++_First;
						delta = 0;
					}
					continue;
				} else	if (_lcpMidLower < lcpLower) {
					// first half
					lcpUpper = _lcpMidLower;
					_LCPTop = leftChild;
					if ((delta = _Delta2) == 1)
						delta = 0;
					continue;
				}
			} else if (lcpLower < lcpUpper) {
                LCPTreeIt rightChild = _LCPTop;
                rightChild.right();
				TDiff _lcpMidUpper = *rightChild;

				if (_lcpMidUpper > lcpUpper) {
					// first half
					_LCPTop.left();
					if ((delta = _Delta2) == 1)
						delta = 0;
					continue;
				} else if (_lcpMidUpper < lcpUpper) {
					// second half
					lcpLower = _lcpMidUpper;
					_First = _Mid;
                    _LCPTop = rightChild;
					if ((delta -= _Delta2) == 1) {
						++_First;
						delta = 0;
					}
					continue;
				}
			}

			TTextIter t = tBegin;
			TQueryIter q = qBegin;
            goFurther(t, *_Mid);

			// lcp search changes MIN to MAX here
			TDiff lcp = Max(lcpLower, lcpUpper);
			goFurther(t, lcp);
			goFurther(q, lcp);
			TDiff max = Min(difference(t, tEnd), difference(q, qEnd));
            TDiff i = max;
			for(; i && *t == *q; --i, ++t, ++q);
			lcp += max - i;

            // is text < query ?
			if (q != qEnd && (t == tEnd || !(*q < *t))) {
				// second half
				lcpLower = lcp;
				_First = _Mid;
				_LCPTop.right();
				if ((delta -= _Delta2) == 1) {
					++_First;
					delta = 0;
				}
			} else
            // is text > query ?
			if (q != qEnd && t != tEnd && (*q < *t)) {
				// first half
				lcpUpper = lcp;
				_LCPTop.left();
				if ((delta = _Delta2) == 1)
					delta = 0;
			} else
			{	// range straddles mid, find each end and return
				TSAIter _First2 = _Lower_bound_lcp_enhanced(
					tBegin,	tEnd,
        			_First,	_Mid,
					_LCPTop.leftChild(),
					qBegin, qEnd,
					lcpLower, lcp);
				_LCPTop.right();
				if ((delta -= _Delta2) == 1) {
					delta = 0;
					++_First;
				}
				TSAIter _Last2 = _Upper_bound_lcp_enhanced(
					tBegin,	tEnd,
        			_Mid, _Last,
					_LCPTop,
					qBegin, qEnd,
					lcp, lcpUpper);
				return Pair<TSAIter> (_First2, _Last2);
			}
		}

		// range straddles mid, find each end and return
		return Pair<TSAIter> (
			_Lower_bound_lcp_enhanced(
				tBegin,	tEnd,
    			_First,	_Last,
				_LCPTop,
				qBegin, qEnd,
				lcpLower, lcpUpper),
			_Upper_bound_lcp_enhanced(
				tBegin,	tEnd,
    			_First,	_Last,
				_LCPTop,
				qBegin, qEnd,
				lcpLower, lcpUpper));
	}

	template <
		typename TTextIter,
		typename TSAIter,
		typename LCPFwdIt,
		typename TQueryIter
	>
	inline Pair<TSAIter> equal_range(
		TTextIter tBegin,
		TTextIter tEnd,
		TSAIter _SAFirst,
		TSAIter _SALast,
		LCPFwdIt _LCPTop,
		TQueryIter qBegin,
		TQueryIter qEnd)
	{	// find first element not before query, using operator<
		return _Equal_range_lcp_enhanced(
			tBegin,	tEnd,
			_SAFirst,	_SALast,
			SearchTreeIterator<LCPFwdIt, LeftCompleteTree>(_LCPTop),
			qBegin, qEnd);
	}
*/

    //////////////////////////////////////////////////////////////////////////////
	// little helpers (not used)
/*
	template < typename LCPFwdIt, typename TSize >
	TSize lcp(LCPFwdIt _First, TSize _Count) {
		if (_Count > 1) {
			TSize lcp = *_First;
			++_First;
			_Count-=2;
			while (_Count) {
				if (lcp < *_First) lcp = *_First;
				++_First;
				--_Count;
			}
			return lcp;
		} else
			return 0;
	}
*/

    //////////////////////////////////////////////////////////////////////////////
	// wrappers with nice interfaces
/*
	template <
		typename TText,
		typename TSA,
		typename TLCPE,
		typename TSubText >
	inline Pair<typename Iterator<TSA const>::Type> equalRangeLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TSubText const &subtext)
	{
		return _Equal_range_lcp_enhanced(
			begin(text), end(text),
			begin(sa), end(sa),
            SearchTreeIterator<typename Iterator<TLCPE const>::Type, LeftCompleteTree>(begin(lcpe), (length(text)>1)?length(text)-1:0),
			begin(subtext), end(subtext));
	}
*/
}
#endif
