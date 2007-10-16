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

	template < class FlatFwdIt, typename TSpec >
	struct TreeIterator {};

	template < class FlatFwdIt, typename TSize, typename TSpec >
    TreeIterator<FlatFwdIt, TSpec> toTreeIterator(FlatFwdIt const &_begin, TSize _size, const TSpec)
    {
        return TreeIterator<FlatFwdIt, TSpec> (_begin, _size);
	}


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search tree like a real tree
	//

	template <class FlatFwdIt>
	struct TreeIterator< FlatFwdIt, SortedList > {

		typedef typename Value<FlatFwdIt>::Type	Type;
		typedef typename Size<FlatFwdIt>::Type	SizeType;

        typedef TreeIterator iterator;

		inline TreeIterator(FlatFwdIt const &_begin, SizeType _size):
			_First(_begin),
			_Count(_size)
		{
			_Count2 = _Count / 2;
			_Mid = _First;
			goFurther(_Mid, _Count2);
		}
			
        inline const Type& operator*() const {
			return *_Mid;
		}

        inline const Type* operator->() const {
			return &*_Mid;
		}

		inline SizeType mid() {
			return _Count2;
		}

		// descend left
		inline iterator& left() {
			_Count = _Count2;
			_Count2 /= 2;
			_Mid = _First;
			goFurther(_Mid, _Count2);
			return *this;
		}

		// descend right
        inline iterator& right() {
			_First = ++_Mid, _Count -= _Count2 + 1;
			_Count2 = _Count / 2;
			goFurther(_Mid, _Count2);
			return *this;
		}

        inline iterator leftChild() const {
            return iterator(*this).left();
        }

        inline iterator rightChild() const {
            return iterator(*this).right();
        }

        inline bool eof() {
            return !_Count;
        }

		inline operator FlatFwdIt() {
			return _Mid;
		}

	private:
		FlatFwdIt _First, _Mid;
		SizeType _Count, _Count2;
	};


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search tree like a real tree
	//

	template <class FlatFwdIt>
	struct TreeIterator< FlatFwdIt, LeftCompleteTree > {

		typedef typename Value<FlatFwdIt>::Type	Type;
		typedef typename Size<FlatFwdIt>::Type	SizeType;

        typedef TreeIterator iterator;

		FlatFwdIt   it;
		SizeType    size;

		inline TreeIterator():
			it(),
            size(0),
            _xSize(0) {}

		inline TreeIterator(FlatFwdIt const &_begin, SizeType _size):
			it(_begin),
			size(_size)
		{
			_left = 0;
			_lSize = 1;
            _iSize = _size;
			for(_xSize = 1; _xSize < _size; _xSize <<= 1);
            if (!_size) _xSize = 0;
		}

        inline const Type& operator*() const {
			return *it;
		}

        inline const Type* operator->() const {
			return &*it;
		}

		inline SizeType mid() {
			return _xSize >> 1;
		}

		inline SizeType leftSize() {
			return mid();
		}

		inline SizeType rightSize() {
			return _iSize - mid();
		}

		// descend left
        inline iterator& left() {
            if (_xSize <= 1) {
                _xSize = 0;
                return *this;
            }
            _descendLeft();
			_iSize = _xSize;	    // = mid();
            return *this;
        }

		// descend right
        inline iterator& right() {
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

        inline iterator leftChild() const {
            return iterator(*this).left();
        }

        inline iterator rightChild() const {
            return iterator(*this).right();
        }

        inline iterator& operator--() {
            --it;
            return *this;
        }

        inline iterator operator--(int) {
            iterator before = *this;
            --it;
            return before;
        }

        inline iterator& operator++() {
            ++it;
            return *this;
        }

        inline iterator operator++(int) {
            iterator before = *this;
            ++it;
            return before;
        }

        inline bool operator==(iterator const &I) {
            return (_xSize == I._xSize) && (_xSize == 0 || it == I.it);
        }

        //operator FlatFwdIt() {
        //    return it;
        //}

        inline bool eof() {
            return !_xSize;
        }

		inline operator FlatFwdIt() {
			return it;
		}

	private:
		SizeType _left;			// left iterator offset of current interval
		SizeType _lSize;		// iterator elements of current level
		SizeType _xSize;		// max interval size of current level
		SizeType _iSize;		// current interval size

        inline void _descendLeft() {
			goFurther(it, _left + _lSize);
			_left <<= 1;
			_xSize >>= 1;
			_lSize = (size + _xSize - 1) / _xSize;
        }

        inline void _descendRight() {
			goFurther(it, _left + 1 + _lSize);
			_left = (_left << 1) + 1;
			_xSize >>= 1;
			_lSize = (size + _xSize - 1) / _xSize;
        }
    };


	//////////////////////////////////////////////////////////////////////////////
	// class to access a flat search b-tree like a real tree
	//

	template < class FlatFwdIt,
			   unsigned BlockSize >
	struct TreeIterator<FlatFwdIt, BTree<BlockSize> > {

		typedef typename Value<FlatFwdIt>::Type	Type;
		typedef typename Size<FlatFwdIt>::Type	SizeType;

		enum { BlockHeight = Log2Floor<BlockSize>::VALUE };
		enum { BlockElements = (1 << BlockHeight) - 1 };
		enum { BlockInnerElements = (1 << (BlockHeight - 1)) - 1 };
		enum { BlockLeafs = 1 << (BlockHeight - 1) };

        typedef TreeIterator iterator;

		FlatFwdIt   it;
		SizeType    size;

		inline TreeIterator():
			it() {}

		inline TreeIterator(FlatFwdIt const &_begin, SizeType _size):
			it(_begin),
			size(_size)
		{
			//_left = 0;
			//_lSize = 0;
   //         _iSize = _size;

			_heightHigh = 1;
			_stepSizeLow = 1;
			for(SizeType _xSizeLow = 2; _xSizeLow <= _size; _xSizeLow <<= 1) {
				if (_stepSizeLow == BlockLeafs) {
					_stepSizeLow = 1;
					++_heightHigh;
				} else
					_stepSizeLow <<= 1;
			}
			
			_stepSizeLow >>= 1;
			for(_xSizeHigh = 1; _xSizeHigh * BlockSize <= _size; _xSizeHigh *= BlockSize);

			_leftLow = (_stepSizeLow << 1) - 1;		// point to the middle
			_leftHigh = 0;
			_lSizeHigh = 1;

			it += _leftLow;

			_updateElements();

			//if (!_size) _xSizeLow = 0;
   //         if (!_size) _xSizeHigh = 0;
		}

        inline const Type& operator*() const {
			return *it;
		}

        inline const Type* operator->() const {
			return &*it;
		}

		// descend left
        inline iterator& left() {
            if (_heightHigh == 1 && !_stepSizeLow) {
                _heightHigh = 0;
                return *this;
            }
            _descendLeft();
            return *this;
        }

		// descend right
        inline iterator& right() {
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

        inline iterator leftChild() const {
            return iterator(*this).left();
        }

        inline iterator rightChild() const {
            return iterator(*this).right();
        }

        inline iterator& operator--() {
            --it;
            return *this;
        }

        inline iterator operator--(int) {
            iterator before = *this;
            --it;
            return before;
        }

        inline iterator& operator++() {
            ++it;
            return *this;
        }

        inline iterator operator++(int) {
            iterator before = *this;
            ++it;
            return before;
        }

        inline bool operator==(iterator const &I) {
			return it == I.it || (eof() && I.eof());
        }

        //operator FlatFwdIt() {
        //    return it;
        //}

        inline bool eof() {
            return !_heightHigh;
        }

		inline operator FlatFwdIt&() {
			return it;
		}

	private:
		SizeType _heightHigh;	// height measured in BBlocks
		unsigned _elements;		// elements in current BBlock

		unsigned _leftLow;		// left iterator offset of current interval
		unsigned _stepSizeLow;	// left/right step size of current level

		SizeType _leftHigh;		// left BBlock offset of current interval
		SizeType _lSizeHigh;	// BBlocks of current level
		SizeType _xSizeHigh;	// max BBlocks of current level

		inline void _descendLeft() {
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

		inline void _descendRight() {
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

		inline void _updateElements() {
			SizeType firstElement = (1 + _leftHigh * BlockSize) * _xSizeHigh - 1;
			SizeType lastElement = (1 + (_leftHigh + 1) * BlockSize) * _xSizeHigh - 2;

			if (lastElement >= size)
				_elements = (size - firstElement) / _xSizeHigh;
			else
				_elements = BlockElements;
		}

        inline void _descendBBlock(SizeType _childIndex) {
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
		class TextRndIt,
		class SAFwdIt,
		class SubTextFwdIt,
		class TSpec
	>
	inline SAFwdIt _Lower_bound_sa(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		TreeIterator< SAFwdIt, TSpec > _Tree,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast)
	{	// find first element not before _Val, using operator<
		typedef typename Difference<TextRndIt>::Type _Diff;
		_Diff _lcpLower = 0;
		_Diff _lcpUpper = 0;

		for (; !_Tree.eof(); )
		{	// divide and conquer, find half that contains answer

			TextRndIt _T = _TextFirst;
			SubTextFwdIt _S = _SubTextFirst;
            _Diff _lcp = Min(_lcpLower, _lcpUpper);

			goFurther(_T, (*_Tree) + _lcp);
			goFurther(_S, _lcp);
			while (_T != _TextLast && _S != _SubTextLast && *_T == *_S) {
				++_T;
				++_S;
				++_lcp;
			}
			
            // is T < S ?
			if (_S != _SubTextLast && (_T == _TextLast || *_T < *_S)) {
				_Tree.right();
				_lcpLower = _lcp;
			} else {
				_Tree.left();
				_lcpUpper = _lcp;
			}
		}
		return _Tree;
	}

	template <
		class TextRndIt,
		class SAFwdIt,
		class SubTextFwdIt,
		class TSpec
	>
	inline SAFwdIt lower_bound(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _First,
		SAFwdIt _Last,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast,
		const TSpec)
	{	// find first element not before _Val, using operator<
		return _Lower_bound_sa(
			_TextFirst, _TextLast,
			toTreeIterator(_First, difference(_First, _Last), TSpec()),
			_SubTextFirst, _SubTextLast);
	}


	//////////////////////////////////////////////////////////////////////////////

	template <
		class TextRndIt,
		class SAFwdIt,
		class SubTextFwdIt,
		class TSpec
	>
	inline SAFwdIt _Upper_bound_sa(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		TreeIterator< SAFwdIt, TSpec > _Tree,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast)
	{	// find first element not before _Val, using operator<
		typedef typename Difference<TextRndIt>::Type _Diff;
		_Diff _lcpLower = 0;
		_Diff _lcpUpper = 0;

		for (; !_Tree.eof(); )
		{	// divide and conquer, find half that contains answer

			TextRndIt _T = _TextFirst;
			SubTextFwdIt _S = _SubTextFirst;
            _Diff _lcp = Min(_lcpLower, _lcpUpper);

			goFurther(_T, (*_Tree) + _lcp);
            goFurther(_S, _lcp);
			while (_T != _TextLast && _S != _SubTextLast && *_T == *_S) {
				++_T;
				++_S;
				++_lcp;
			}
			
            // is T <= S ?
			if (_S == _SubTextLast || _T == _TextLast || !(*_S < *_T)) {
				_Tree.right();
				_lcpLower = _lcp;
			} else {
				_Tree.left();
				_lcpUpper = _lcp;
			}
		}
		return _Tree;
	}

	template <
		class TextRndIt,
		class SAFwdIt,
		class SubTextFwdIt,
		class TSpec
	>
	inline SAFwdIt upper_bound(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _First,
		SAFwdIt _Last,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast,
		const TSpec)
	{	// find first element that _Val is before, using operator<
		return _Upper_bound_sa(
			_TextFirst, _TextLast,
			toTreeIterator(_First, difference(_First, _Last), TSpec()),
			_SubTextFirst, _SubTextLast);
	}


	//////////////////////////////////////////////////////////////////////////////

	template <
		class TextRndIt,
		class SAFwdIt,
		class SubTextFwdIt,
		class TSpec
	>
	inline Pair<SAFwdIt> _Equal_range_sa(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		TreeIterator< SAFwdIt, TSpec > _Tree,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast)
	{	// find range equivalent to _Val, using operator<
		typedef typename Difference<TextRndIt>::Type _Diff;
		_Diff _lcpLower = 0;
		_Diff _lcpUpper = 0;

		for (; !_Tree.eof(); )
		{	// divide and conquer, check midpoint

			TextRndIt _T = _TextFirst;
			SubTextFwdIt _S = _SubTextFirst;
			_Diff _lcp = Min(_lcpLower, _lcpUpper);
			goFurther(_T, (*_Tree) + _lcp);
			goFurther(_S, _lcp);
			while (_T != _TextLast && _S != _SubTextLast && *_T == *_S) {
				++_T;
				++_S;
				++_lcp;
			}
			
            // is T < S ?
			if (_S != _SubTextLast && (_T == _TextLast || *_T < *_S))
			{	// range begins above _Mid, loop
				_Tree.right();
				_lcpLower = _lcp;
			}
            // is T > S ?
			else if (_S != _SubTextLast && (_T != _TextLast && *_S < *_T))
			{	// range in first half, loop
				_Tree.left();
				_lcpUpper = _lcp;
			} else
            // is T == S ?
			{	// range straddles mid, find each end and return
				SAFwdIt _First2 = _Lower_bound_sa(
					_TextFirst,	_TextLast,
					_Tree.leftChild(),
					_SubTextFirst, _SubTextLast);
				SAFwdIt _Last2 = _Upper_bound_sa(
					_TextFirst,	_TextLast,
					_Tree.rightChild(),
					_SubTextFirst, _SubTextLast);
				return Pair<SAFwdIt> (_First2, _Last2);
			}
		}
		return Pair<SAFwdIt> (_Tree, _Tree);	// empty range
	}

	template <
		class TextRndIt,
		class SAFwdIt,
		class SubTextFwdIt,
		class TSpec
	>
	inline Pair<SAFwdIt> equal_range(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _First,
		SAFwdIt _Last,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast)
	{	// find range equivalent to _Val, using operator<
		return _Equal_range_sa(
			_TextFirst, _TextLast,
			toTreeIterator(_First, difference(_First, _Last), TSpec()),
			_SubTextFirst, _SubTextLast);
	}


	//////////////////////////////////////////////////////////////////////////////
	// substring search with enhanced LCP-table
	//

	template <
		class TextRndIt,
		class SAFwdIt,
		class LCPTreeIt,
		class SubTextFwdIt
	>
	inline SAFwdIt _Lower_bound_lcp_enhanced(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _First,
		SAFwdIt _Last,
		LCPTreeIt _LCPTop,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast,
		typename Difference<TextRndIt>::Type _lcpLower,
		typename Difference<TextRndIt>::Type _lcpUpper)
	{	// find first element not before _Val, using operator<
		typedef typename Difference<TextRndIt>::Type _Diff;
        _Diff _Delta = difference(_First, _Last) - 1;
		_Diff _lcp;
		#ifdef SEQAN_PROFILE_LCPEFIND
			_Diff skippedCompares = 0;	// difference of char compares related to xxx_bound_sa
		#endif

        // binary search with intervals >= 3 elements
		for (; 1 < _Delta; )
		{	// divide and conquer, find half that contains answer
			_Diff _Delta2 = _LCPTop.leftSize();
			SAFwdIt _Mid = _First;
			goFurther(_Mid, _Delta2);

			if (_lcpLower > _lcpUpper) {
                LCPTreeIt leftChild = _LCPTop;
                leftChild.left();
				_Diff _lcpMidLower = *leftChild;

				if (_lcpMidLower > _lcpLower) {
					// second half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += _lcpLower - _lcpUpper;
					#endif
					_First = _Mid;
					_LCPTop.right();
					if ((_Delta -= _Delta2) == 1) {
						++_First;
						_Delta = 0;
					}
					continue;
				} else	if (_lcpMidLower < _lcpLower) {
					// first half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += _lcpMidLower - _lcpUpper;
					#endif
					_lcpUpper = _lcpMidLower;
					_LCPTop = leftChild;
					if ((_Delta = _Delta2) == 1)
						_Delta = 0;
					continue;
				}
				_lcp = _lcpLower;
			} else if (_lcpLower < _lcpUpper) {
                LCPTreeIt rightChild = _LCPTop;
                rightChild.right();
				_Diff _lcpMidUpper = *rightChild;

				if (_lcpMidUpper > _lcpUpper) {
					// first half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += _lcpUpper - _lcpLower;
					#endif
					_LCPTop.left();
					if ((_Delta = _Delta2) == 1)
						_Delta = 0;
					continue;
				} else if (_lcpMidUpper < _lcpUpper) {
					// second half
					#ifdef SEQAN_PROFILE_LCPEFIND
						skippedCompares += _lcpMidUpper - _lcpLower;
					#endif
					_lcpLower = _lcpMidUpper;
					_First = _Mid;
                    _LCPTop = rightChild;
					if ((_Delta -= _Delta2) == 1) {
						++_First;
						_Delta = 0;
					}
					continue;
				}
				_lcp = _lcpUpper;
			} else
				_lcp = _lcpUpper;

			TextRndIt _T = _TextFirst;
			SubTextFwdIt _S = _SubTextFirst;
            goFurther(_T, *_Mid);

			// lcp search changes MIN to MAX here
//			_Diff _lcp = Max(_lcpLower, _lcpUpper);
			#ifdef SEQAN_PROFILE_LCPEFIND
				skippedCompares += _lcp - Min(_lcpLower, _lcpUpper);
			#endif
			goFurther(_T, _lcp);
			goFurther(_S, _lcp);
			for(_Diff i = Min(difference(_T, _TextLast), difference(_S, _SubTextLast)); 
				i && *_T == *_S;
				--i, ++_T, ++_S, ++_lcp);

            // is T < S ?
			if (_S != _SubTextLast && (_T == _TextLast || *_T < *_S)) {
				// second half
				_lcpLower = _lcp;
				_First = _Mid;
				_LCPTop.right();
				if ((_Delta -= _Delta2) == 1) {
					++_First;
					_Delta = 0;
				}
			} else {
				// first half
				_lcpUpper = _lcp;
				_LCPTop.left();
				if ((_Delta = _Delta2) == 1)
					_Delta = 0;
			}
		}

		#ifdef SEQAN_PROFILE_LCPEFIND
			SEQAN_PROADD(SEQAN_PROEXTRA3, skippedCompares);
		#endif

        // binary search for intervals of 2 or less elements
        /*_Diff*/ _lcp = Min(_lcpLower, _lcpUpper);

        SubTextFwdIt _S = _SubTextFirst;
		goFurther(_S, _lcp);

        while (true) {
			TextRndIt _T = _TextFirst;
			goFurther(_T, *_First + _lcp);

			for(_Diff i = Min(difference(_T, _TextLast), difference(_S, _SubTextLast)); 
            	i && *_T == *_S;
            	 --i, ++_T, ++_S);

            // is T < S ?
			if (_S != _SubTextLast && (_T == _TextLast || *_T < *_S)) {
				// second half
				++_First;
				if (!_Delta) return _First;
                --_Delta;
			} else {
				// first half -> end
				return _First;
			}
        }
	}

	template <
		class TextRndIt,
		class SAFwdIt,
		class LCPFwdIt,
		class SubTextFwdIt
	>
	inline SAFwdIt lower_bound(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _First,
		SAFwdIt _Last,
		TreeIterator<LCPFwdIt, LeftCompleteTree> _LCPTop,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast)
	{	// find first element not before _Val, using operator<
		return _Lower_bound_lcp_enhanced(
			_TextFirst, _TextLast,
			_First, _Last,
			_LCPTop,
			_SubTextFirst, _SubTextLast,
			0, 0);
	}


	//////////////////////////////////////////////////////////////////////////////

	template <
		class TextRndIt,
		class SAFwdIt,
		class LCPTreeIt,
		class SubTextFwdIt
	>
	inline SAFwdIt _Upper_bound_lcp_enhanced(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _First,
		SAFwdIt _Last,
		LCPTreeIt _LCPTop,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast,
		typename Difference<TextRndIt>::Type _lcpLower,
		typename Difference<TextRndIt>::Type _lcpUpper)
	{	// find first element not before _Val, using operator<
		typedef typename Difference<TextRndIt>::Type _Diff;
        _Diff _Delta = difference(_First, _Last) - 1;

        // binaray search with intervals >= 3 elements
		for (; 1 < _Delta; )
		{	// divide and conquer, find half that contains answer
			_Diff _Delta2 = _LCPTop.leftSize();
			SAFwdIt _Mid = _First;
			goFurther(_Mid, _Delta2);

			if (_lcpLower > _lcpUpper) {
                LCPTreeIt leftChild = _LCPTop;
                leftChild.left();
				_Diff _lcpMidLower = *leftChild;

				if (_lcpMidLower > _lcpLower) {
					// second half
					_First = _Mid;
					_LCPTop.right();
					if ((_Delta -= _Delta2) == 1) {
						++_First;
						_Delta = 0;
					}
					continue;
				} else	if (_lcpMidLower < _lcpLower) {
					// first half
					_lcpUpper = _lcpMidLower;
					_LCPTop = leftChild;
					if ((_Delta = _Delta2) == 1)
						_Delta = 0;
					continue;
				}
			} else if (_lcpLower < _lcpUpper) {
                LCPTreeIt rightChild = _LCPTop;
                rightChild.right();
				_Diff _lcpMidUpper = *rightChild;

				if (_lcpMidUpper > _lcpUpper) {
					// first half
					_LCPTop.left();
					if ((_Delta = _Delta2) == 1)
						_Delta = 0;
					continue;
				} else if (_lcpMidUpper < _lcpUpper) {
					// second half
					_lcpLower = _lcpMidUpper;
					_First = _Mid;
                    _LCPTop = rightChild;
					if ((_Delta -= _Delta2) == 1) {
						++_First;
						_Delta = 0;
					}
					continue;
				}
			}

			TextRndIt _T = _TextFirst;
			SubTextFwdIt _S = _SubTextFirst;
            goFurther(_T, *_Mid);

			// lcp search changes MIN to MAX here
			_Diff _lcp = Max(_lcpLower, _lcpUpper);
			goFurther(_T, _lcp);
			goFurther(_S, _lcp);
			_Diff max = Min(difference(_T, _TextLast), difference(_S, _SubTextLast));
            _Diff i = max;
			for(; i && *_T == *_S; --i, ++_T, ++_S);
			_lcp += max - i;

            // is T <= S ?
			if (_S == _SubTextLast || _T == _TextLast || !(*_S < *_T)) {
				// second half
				_lcpLower = _lcp;
				_First = _Mid;
				_LCPTop.right();
				if ((_Delta -= _Delta2) == 1) {
					++_First;
					_Delta = 0;
				}
			} else {
				// first half
				_lcpUpper = _lcp;
				_LCPTop.left();
				if ((_Delta = _Delta2) == 1)
					_Delta = 0;
			}
		}

        // binary search for intervals of 2 or less elements
        _Diff _lcp = Min(_lcpLower, _lcpUpper);

        SubTextFwdIt _S = _SubTextFirst;
		goFurther(_S, _lcp);

        while (true) {
			TextRndIt _T = _TextFirst;
			goFurther(_T, *_First + _lcp);

			_Diff i = Min(difference(_T, _TextLast), difference(_S, _SubTextLast));
            for(; i && *_T == *_S; --i, ++_T, ++_S);

            // is T <= S ?
			if (_S == _SubTextLast || _T == _TextLast || !(*_S < *_T)) {
				// second half
				++_First;
				if (!_Delta) return _First;
                --_Delta;
			} else {
				// first half -> end
				return _First;
			}
        }
	}

	template <
		class TextRndIt,
		class SAFwdIt,
		class LCPFwdIt,
		class SubTextFwdIt
	>
	inline SAFwdIt upper_bound(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _First,
		SAFwdIt _Last,
		TreeIterator<LCPFwdIt, LeftCompleteTree> _LCPTop,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast)
	{	// find first element not before _Val, using operator<
		return _Upper_bound_lcp_enhanced(
			_TextFirst, _TextLast,
			_First, _Last,
			_LCPTop,
			_SubTextFirst, _SubTextLast,
			0, 0);
	}


	//////////////////////////////////////////////////////////////////////////////

	template <
		class TextRndIt,
		class SAFwdIt,
		class LCPTreeIt,
		class SubTextFwdIt
	>
	inline Pair<SAFwdIt> _Equal_range_lcp_enhanced(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _First,
		SAFwdIt _Last,
		LCPTreeIt _LCPTop,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast)
	{	// find first element not before _Val, using operator<
		typedef typename Difference<TextRndIt>::Type _Diff;
		_Diff _lcpLower = 0;
		_Diff _lcpUpper = 0;
        _Diff _Delta = difference(_First, _Last) - 1;

        // binaray search with intervals >= 3 elements
		for (; 1 < _Delta; )
		{	// divide and conquer, find half that contains answer
			_Diff _Delta2 = _LCPTop.leftSize();
			SAFwdIt _Mid = _First;
			goFurther(_Mid, _Delta2);

			if (_lcpLower > _lcpUpper) {
                LCPTreeIt leftChild = _LCPTop;
                leftChild.left();
				_Diff _lcpMidLower = *leftChild;

				if (_lcpMidLower > _lcpLower) {
					// second half
					_First = _Mid;
					_LCPTop.right();
					if ((_Delta -= _Delta2) == 1) {
						++_First;
						_Delta = 0;
					}
					continue;
				} else	if (_lcpMidLower < _lcpLower) {
					// first half
					_lcpUpper = _lcpMidLower;
					_LCPTop = leftChild;
					if ((_Delta = _Delta2) == 1)
						_Delta = 0;
					continue;
				}
			} else if (_lcpLower < _lcpUpper) {
                LCPTreeIt rightChild = _LCPTop;
                rightChild.right();
				_Diff _lcpMidUpper = *rightChild;

				if (_lcpMidUpper > _lcpUpper) {
					// first half
					_LCPTop.left();
					if ((_Delta = _Delta2) == 1)
						_Delta = 0;
					continue;
				} else if (_lcpMidUpper < _lcpUpper) {
					// second half
					_lcpLower = _lcpMidUpper;
					_First = _Mid;
                    _LCPTop = rightChild;
					if ((_Delta -= _Delta2) == 1) {
						++_First;
						_Delta = 0;
					}
					continue;
				}
			}

			TextRndIt _T = _TextFirst;
			SubTextFwdIt _S = _SubTextFirst;
            goFurther(_T, *_Mid);

			// lcp search changes MIN to MAX here
			_Diff _lcp = Max(_lcpLower, _lcpUpper);
			goFurther(_T, _lcp);
			goFurther(_S, _lcp);
			_Diff max = Min(difference(_T, _TextLast), difference(_S, _SubTextLast));
            _Diff i = max;
			for(; i && *_T == *_S; --i, ++_T, ++_S);
			_lcp += max - i;

            // is T < S ?
			if (_S != _SubTextLast && (_T == _TextLast || !(*_S < *_T))) {
				// second half
				_lcpLower = _lcp;
				_First = _Mid;
				_LCPTop.right();
				if ((_Delta -= _Delta2) == 1) {
					++_First;
					_Delta = 0;
				}
			} else
            // is T > S ?
			if (_S != _SubTextLast && _T != _TextLast && (*_S < *_T)) {
				// first half
				_lcpUpper = _lcp;
				_LCPTop.left();
				if ((_Delta = _Delta2) == 1)
					_Delta = 0;
			} else
			{	// range straddles mid, find each end and return
				SAFwdIt _First2 = _Lower_bound_lcp_enhanced(
					_TextFirst,	_TextLast,
        			_First,	_Mid,
					_LCPTop.leftChild(),
					_SubTextFirst, _SubTextLast,
					_lcpLower, _lcp);
				_LCPTop.right();
				if ((_Delta -= _Delta2) == 1) {
					_Delta = 0;
					++_First;
				}
				SAFwdIt _Last2 = _Upper_bound_lcp_enhanced(
					_TextFirst,	_TextLast,
        			_Mid, _Last,
					_LCPTop,
					_SubTextFirst, _SubTextLast,
					_lcp, _lcpUpper);
				return Pair<SAFwdIt> (_First2, _Last2);
			}
		}

		// range straddles mid, find each end and return
		return Pair<SAFwdIt> (
			_Lower_bound_lcp_enhanced(
				_TextFirst,	_TextLast,
    			_First,	_Last,
				_LCPTop,
				_SubTextFirst, _SubTextLast,
				_lcpLower, _lcpUpper),
			_Upper_bound_lcp_enhanced(
				_TextFirst,	_TextLast,
    			_First,	_Last,
				_LCPTop,
				_SubTextFirst, _SubTextLast,
				_lcpLower, _lcpUpper));
	}

	template <
		class TextRndIt,
		class SAFwdIt,
		class LCPFwdIt,
		class SubTextFwdIt
	>
	inline Pair<SAFwdIt> equal_range(
		TextRndIt _TextFirst,
		TextRndIt _TextLast,
		SAFwdIt _SAFirst,
		SAFwdIt _SALast,
		LCPFwdIt _LCPTop,
		SubTextFwdIt _SubTextFirst,
		SubTextFwdIt _SubTextLast)
	{	// find first element not before _Val, using operator<
		return _Equal_range_lcp_enhanced(
			_TextFirst,	_TextLast,
			_SAFirst,	_SALast,
			TreeIterator<LCPFwdIt, LeftCompleteTree>(_LCPTop),
			_SubTextFirst, _SubTextLast);
	}


    //////////////////////////////////////////////////////////////////////////////
	// little helpers (not used)
/*
	template < class LCPFwdIt, typename SizeType >
	SizeType lcp(LCPFwdIt _First, SizeType _Count) {
		if (_Count > 1) {
			SizeType lcp = *_First;
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
	// find first element not before _Val, using operator<

	template <
		class TText,
		class TSA,
		class TSubText >
	inline typename Iterator<TSA const>::Type lowerBoundSA(
		TText const &text,
		TSA const &sa,
		TSubText const &subtext)
	{
		return _Lower_bound_sa(
			begin(text), end(text),
			toTreeIterator(begin(sa), length(sa), SortedList()),
			begin(subtext), end(subtext));
	}

    //////////////////////////////////////////////////////////////////////////////
	// wrappers with nice interfaces

	template <
		class TText,
		class TSA,
		class TSubText >
	inline Pair<typename Iterator<TSA const>::Type> equalRangeSA(
		TText const &text,
		TSA const &sa,
		TSubText const &subtext)
	{
		return _Equal_range_sa(
			begin(text), end(text),
			toTreeIterator(begin(sa), length(sa), SortedList()),
			begin(subtext), end(subtext));
	}

	template <
		class TText,
		class TSA,
		class TLCPE,
		class TSubText >
	inline Pair<typename Iterator<TSA const>::Type> equalRangeLCPE(
		TText const &text,
		TSA const &sa,
		TLCPE const &lcpe,
		TSubText const &subtext)
	{
		return _Equal_range_lcp_enhanced(
			begin(text), end(text),
			begin(sa), end(sa),
            TreeIterator<typename Iterator<TLCPE const>::Type, LeftCompleteTree>(begin(lcpe), (length(text)>1)?length(text)-1:0),
			begin(subtext), end(subtext));
	}

}
#endif
