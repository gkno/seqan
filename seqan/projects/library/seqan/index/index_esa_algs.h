/*
 *  index_esa_algs.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_ESA_ALGS_H
#define SEQAN_HEADER_INDEX_ESA_ALGS_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	//////////////////////////////////////////////////////////////////////////////
	// more sophisticated algorithms on suffix trees / enhanced suffix arrays
	//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// MUMs - generalized suffix tree version
	//////////////////////////////////////////////////////////////////////////////

	template < typename TSTree >
	struct GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<MUMs> > > > {
		typedef PostorderEmptyEdges	Type;
	};

	template < typename TSTree >
	class Iter< TSTree, VSTree< BottomUp<MUMs> > >:
		public Iter< TSTree, VSTree< BottomUp<> > > 
	{
	public:
		typedef Iter< TSTree, VSTree< BottomUp<> > >	TBase;
		typedef typename Size<TSTree>::Type				TSize;
//____________________________________________________________________________

		TSize						minLength;
		TSize						seqCount;
		VectorSet<TSize, Alloc<> >	seqSet;
//____________________________________________________________________________

		Iter(TSTree &_tree):
			TBase(_tree),
			minLength(1),
			seqCount(countSequences(_tree)),
			seqSet(countSequences(_tree))
		{
			indexRequire(_tree, ESA_BWT());
			goNext(*this);	// the iterator starts in a suffix, i.e. not a MUM node (length(occ)<2<=seqCount)
		}

		Iter(TSTree &_tree, TSize _minLength):
			TBase(_tree),
			minLength(_minLength),
			seqCount(countSequences(_tree)),
			seqSet(countSequences(_tree))
		{
			indexRequire(_tree, ESA_BWT());
			goNext(*this);	// the iterator starts in a suffix, i.e. not a MUM node (length(occ)<2<=seqCount)
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			seqCount(countSequences(container(_origin))),
			seqSet(countSequences(container(_origin))) {}
	};

	template < typename TSTree >
	inline void goNext(Iter< TSTree, VSTree< BottomUp<MUMs> > > &it) {
		do {
			goNext(it, PostorderEmptyEdges());
		} while (!atEnd(it) && 
			     !(	(countOccurences(it) == it.seqCount) && 
					(repLength(it) >= it.minLength) &&
					isUnique(it, it.seqSet) && 
					isLeftMaximal(it)) );
	}


	//////////////////////////////////////////////////////////////////////////////
	// super-maximal repeats - suffix tree version
	//////////////////////////////////////////////////////////////////////////////

	template < typename TSTree >
	struct _GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > > > > {
		typedef PostorderEmptyEdges	Type;
	};

	template < typename TSTree >
	class Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > >:
		public Iter< TSTree, VSTree< BottomUp<> > > 
	{
	public:
		typedef Iter< TSTree, VSTree< BottomUp<> > >	TBase;
		typedef typename Size<TSTree>::Type				TSize;
		typedef typename Value<TSTree>::Type			TValue;
//____________________________________________________________________________

		TSize						minLength;
		typename Set<TValue>::Type	charSet;
//____________________________________________________________________________

		Iter(TSTree &_tree):
			TBase(_tree),
			minLength(1)
		{
			indexRequire(_tree, ESA_ChildTab());
			indexRequire(_tree, ESA_BWT());
			goNext(*this);	// the iterator starts in a suffix, i.e. not a repeat
		}

		Iter(TSTree &_tree, TSize _minLength):
			TBase(_tree),
			minLength(_minLength)
		{
			indexRequire(_tree, ESA_BWT());
			goNext(*this);	// the iterator starts in a suffix, i.e. not a repeat
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			charSet(_origin.charSet) {}
	};

	template < typename TSTree >
	inline void goNext(Iter< TSTree, VSTree< BottomUp<SuperMaxRepeats> > > &it) {
		do {
			goNext(it, PostorderEmptyEdges());
		} while (!atEnd(it) && 
			     !(	childrenAreLeaves(it) && 
					(repLength(it) >= it.minLength) &&
					!isPartiallyLeftExtensible(it, it.charSet)) );
	}
		

	//////////////////////////////////////////////////////////////////////////////
	// supermaximal repeats - specialized for Enhanced Suffix Arrays
	//////////////////////////////////////////////////////////////////////////////

	template < typename TSTree >
	struct _GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<SuperMaxRepeatsFast> > > > {
		typedef PostorderEmptyEdges	Type;
	};

	template < typename TText, typename TSpec >
	class Iter< Index<TText, Index_ESA<TSpec> >, VSTree< BottomUp<SuperMaxRepeatsFast> > >:
		public Iter< Index<TText, Index_ESA<TSpec> >, VSTree< BottomUp<> > > 
	{
	public:
		typedef Index< TText, Index_ESA<TSpec> >		TIndex;
		typedef Iter< TIndex, VSTree< BottomUp<> > >	TBase;
		typedef typename Size<TIndex>::Type				TSize;
		typedef typename Value<TIndex>::Type			TValue;

		typedef typename Iterator<typename Fibre<TIndex, ESA_LCP>::Type const>::Type	TLCPIter;
//____________________________________________________________________________

		TSize		minLength;
		TLCPIter	lIter, lEnd;	// lcp table iterators (optimization)
		TSize		lValueLast;		// current l-value of interval
		bool 		rising;			// is the left interval border valid
		typename Set<TValue>::Type	charSet;
//____________________________________________________________________________

		Iter(TIndex &_index):
			TBase(_index),
			minLength(1),
			lValueLast(0),
			rising(true)
		{
			this->vDesc.i1 = Pair<TSize>(0,0);
			indexRequire(_index, ESA_BWT());
			lIter = begin(indexLCP(this->index));
			lEnd  = end(indexLCP(this->index));
			goNext(*this);
		}

		Iter(TIndex &_index, TSize _minLength):
			TBase(_index),
			minLength(_minLength),
			lValueLast(0),
			rising(true)
		{
			this->vDesc.i1 = Pair<TSize>(0,0);
			indexRequire(_index, ESA_BWT());
			lIter = begin(indexLCP(this->index));
			lEnd  = end(indexLCP(this->index));
			goNext(*this);
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			lIter(_origin.lIter),
			lEnd(_origin.lEnd),
			lValueLast(_origin.lValueLast),
			rising(_origin.rising) {}
	};

	template < typename TText, typename TSpec >
	inline void goNext(Iter< Index<TText, Index_ESA<TSpec> >, VSTree< BottomUp<SuperMaxRepeatsFast> > > &it) 
	{
		typedef Index<TText, Index_ESA<TSpec> >		TIndex;
		typename Size<TIndex>::Type					lcp;

		while (it.lIter != it.lEnd) {
			lcp = *it.lIter;

			if (lcp < it.lValueLast) {
				if (it.rising) {
					if (it.lValueLast > it.minLength) {
						_dfsLCP(it) = it.lValueLast;
						++_dfsRange(it).i2;
						++it.lIter;
						it.lValueLast = lcp;
						if (!isPartiallyLeftExtensible(it, it.charSet)) return;
						continue;
					}
					it.rising = false;
				}
			} else
			if (lcp > it.lValueLast) {
				_dfsRange(it).i1 = _dfsRange(it).i2;
				it.rising = true;
			}

			++_dfsRange(it).i2;
			++it.lIter;
			it.lValueLast = lcp;
		}
		_dfsRange(it).i2 = 0;
		return;
	}

	
	//////////////////////////////////////////////////////////////////////////////
	// maximal repeats - suffix tree version
	//////////////////////////////////////////////////////////////////////////////

	// contains a list of indices of the same bwt value (fraction)
	template <typename TSize>
	struct _FractionHeader {
		TSize	begin, end;
		TSize	size;
		_FractionHeader() {}
		_FractionHeader(TSize _begin, TSize _end, TSize _size):
			begin(_begin), end(_end), size(_size) {}
	};

	// contains a set of fractions (one for each bwt value) 
	// and a fraction for the undefined bwt value (for the virtual character at position -1)
	template <typename TValue, typename TSize>
	struct _FractionCompound {
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;	// TFraction = (c,(begin,end))	c..char, begin/end indices in posList
		typedef typename Set<TFraction>::Type	TSet;

		TSet			set;
		TFractionHeader	leftmost;

		_FractionCompound():
			leftmost(0,0,0) {}
	};

	template < typename TSTree >
	struct _GetVSTreeIteratorTraits< Iter< TSTree, VSTree< BottomUp<MaxRepeats> > > > {
		typedef PostorderEmptyEdges	Type;
	};

	template < typename TSTree >
	class Iter< TSTree, VSTree< BottomUp<MaxRepeats> > >:
		public Iter< TSTree, VSTree< BottomUp<> > > 
	{
	public:
		typedef Iter< TSTree, VSTree< BottomUp<> > >	TBase;
		typedef typename Size<TSTree>::Type				TSize;
		typedef typename Value<TSTree>::Type			TValue;

		typedef _FractionCompound<TValue, TSize>	TFractionCompound;
		typedef String<TFractionCompound, Block<> >	TSetStack;
		typedef String<TSize>						TPositionList;
		
		typedef typename TFractionCompound::TSet	TSet;
		typedef typename Iterator<TSet>::Type		TSetIterator;

		typedef typename TBase::TStackEntry			TStackEntry;

//____________________________________________________________________________

		TSize			minLength;
		TSetStack		setStack;
		TPositionList	posList;	// this list is indexed just as SA is and contains the next entry's index
		bool			canMerge;	// is false, if parent node appears after its first child on stack
//____________________________________________________________________________

		Iter(TSTree &_index):
			TBase(_index, MinimalCtor()),
			minLength(1),
			canMerge(true)
		{
			indexRequire(_index, ESA_SA());
			indexRequire(_index, ESA_LCP());
			indexRequire(_index, ESA_BWT());
			resize(posList, length(_index));

			if (!empty(indexSA(_index))) 
			{
				_dfsOnPush(*this, TStackEntry(0,0));
				goNext(*this);
			}
		}

		Iter(TSTree &_index, TSize _minLength):
			TBase(_index, MinimalCtor()),
			minLength(_minLength),
			canMerge(true)
		{
			indexRequire(_index, ESA_SA());
			indexRequire(_index, ESA_LCP());
			indexRequire(_index, ESA_BWT());
			resize(posList, length(_index));

			if (!empty(indexSA(_index))) 
			{
				_dfsOnPush(*this, TStackEntry(0,0));
				goNext(*this);
			}
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			setStack(_origin.setStack),
			posList(_origin.posList),
			canMerge(_origin.canMerge) {}

//____________________________________________________________________________

		inline bool hasRepeats() 
		{
			if (length(setStack) < 2) return false;

			TFractionCompound &child  = top(setStack);
			TFractionCompound &parent = topPrev(setStack);

			TSize cs = length(child.set), ps = length(parent.set);

			if (child.leftmost.size  > 0) ++cs;
			if (parent.leftmost.size > 0) ++ps;

			if (cs == 0 || ps == 0) return false;
			if (cs  > 1 || ps  > 1) return true;

			if (child.leftmost.size > 0 || parent.leftmost.size > 0)
				return true;

			return keyOf(begin(child.set)) != keyOf(begin(parent.set));
		}

		inline TSize countRepeats() 
		{
			if (length(setStack) < 2) return 0;

			TFractionCompound &child  = top(setStack);
			TFractionCompound &parent = topPrev(setStack);

			TSetIterator child_fraction		= begin(child.set);
			TSetIterator child_end			= end(child.set);
			TSetIterator parent_fraction	= begin(parent.set);
			TSetIterator parent_end			= end(parent.set);

			TSize sum = 0;
			for(; child_fraction != child_end; ++child_fraction) {
				for(; parent_fraction != parent_end; ++parent_fraction) {
					if (keyOf(child_fraction) != keyOf(parent_fraction))
						sum += (*child_fraction).size * (*parent_fraction).size;

					sum += child.leftmost.size * (*parent_fraction).size;
				}
				sum += (*child_fraction).size * parent.leftmost.size;
			}
			sum += child.leftmost.size * parent.leftmost.size;
			return sum;
		}
//____________________________________________________________________________

		inline void _dump() const {
			::std::cerr << "SETSTACK of " << representative(*this) << ":" << ::std::endl;
			typename Iterator<TSetStack const>::Type it = begin(setStack), itEnd = end(setStack);
			while (it != itEnd) {
				TSet const &set = (*it).set;
				typename Iterator<TSet const>::Type sit = begin(set), sitEnd = end(set);

				while (sit != sitEnd) {
					::std::cerr << keyOf(sit) << "::";
					typename TFractionCompound::TFractionHeader head = objectOf(sit);
					TSize i = head.begin;
					while (!_isSizeInval(i)) {
						::std::cerr << i << "  ";
						i = posList[i];
					}
					::std::cerr << ::std::endl;
					++sit;
				}

				::std::cerr << "_________________________" << ::std::endl;
				++it;
			}
		}
	};

	template < typename TSTree >
	inline typename VertexDescriptor<TSTree>::Type 
	value(Iter< TSTree, VSTree< BottomUp<MaxRepeats> > > const &it) 
	{
		if (empty(it.history))
			return it.vDesc;
		typedef typename VertexDescriptor<TSTree>::Type TDesc;
		typedef typename Value<TDesc, 1>::Type			TRange;
		return TDesc(TRange(top(it.history).i1, it.vDesc.i1.i2), 0);
	}

	template < typename TSTree >
	inline typename Size<TSTree>::Type 
	repLength(Iter< TSTree, VSTree< BottomUp<MaxRepeats> > > const &it) 
	{
		return top(it.history).i2;
	}

	// add bwt partitions of child to parent node
	template < typename TSTree, typename TFractionCompound >
	inline void maxRepeatMerge(
		Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > &it, 
		TFractionCompound &parent,
		TFractionCompound &child)
	{
		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;
		typedef typename Set<TFraction>::Type	TSet;
		typedef typename Iterator<TSet>::Type	TSetIterator;

		TSetIterator _end = end(child.set);
		for(TSetIterator i = begin(child.set); i != _end; ++i) {
			if (in(keyOf(i), parent.set)) {	// append child fraction to parent's fraction
				TFractionHeader &parent_header = objectOf(find(keyOf(i), parent.set));
				TFractionHeader &child_header = objectOf(i);
				it.posList[parent_header.end] = child_header.begin;
				parent_header.end = child_header.end;
				parent_header.size += child_header.size;
			} else
				insert(TFraction(keyOf(i), objectOf(i)), parent.set);	// insert child fraction in parent's set
		}
		if (parent.leftmost.size > 0) {
			if (child.leftmost.size > 0) {
				it.posList[parent.leftmost.end] = child.leftmost.begin;
				parent.leftmost.end = child.leftmost.end;
				parent.leftmost.size += child.leftmost.size;
			}
		} else
			parent.leftmost = child.leftmost;
	}

	// maximal repeat push/leaf handlers of lcp-dfs-traversal
	template < typename TSTree, typename TElement >
	inline void _dfsOnPush(Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > &it, TElement const &e) 
	{
		typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
		_dfsOnPush((TBase&)it, e);

		if (it.canMerge)
			push(it.setStack);
/*
		::std::cerr << "PUSH ";
		_dumpHistoryStack(it);
		it._dump();
*/	}

	template < typename TSTree >
	inline void _dfsOnLeaf(Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > &it) 
	{
		typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
		_dfsOnLeaf((TBase&)it);

		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;

		push(it.setStack);

		TSize index = _dfsRange(it).i1;
		if (!posAtFirstLocal(saAt(index, container(it)), stringSetLimits(it))) 
			insert(
				TFraction(
					bwtAt(index, container(it)),
					TFractionHeader(index, index, 1)), 
				top(it.setStack).set);
		else
			top(it.setStack).leftmost = TFractionHeader(index, index, 1);

		_setSizeInval(it.posList[index]);
/*
		::std::cerr << "LEAF ";
		_dumpHistoryStack(it);
		it._dump();
*/	}

	template < typename TSTree >
	inline void goNext(Iter< TSTree, VSTree< BottomUp<MaxRepeats> > > &it) {
		do {
			if (it.canMerge && length(it.setStack) >= 2) {
				maxRepeatMerge(it, topPrev(it.setStack), top(it.setStack));
				pop(it.setStack);
			}
			goNext(it, PostorderEmptyEdges());
			if (empty(it.history))
				it.canMerge = false;
			else
				it.canMerge = !_dfsReversedOrder(it);
		} while (!eof(it) && !(it.canMerge && (repLength(it) >= it.minLength) && it.hasRepeats()));
	}





	//////////////////////////////////////////////////////////////////////////////
	// maximal repeat representation
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSTree>
	struct MaxRepeat {
		Iter< TSTree, VSTree<BottomUp<MaxRepeats> > > &it;
	};

	template <typename TSTree>
	struct Value< MaxRepeat<TSTree> > {
		typedef Pair< typename SAValue<TSTree>::Type > Type;
	};

	template <typename TSTree>
	struct Size< MaxRepeat<TSTree> > {
		typedef typename Size<TSTree>::Type Type;
	};


	template <typename TSTree>
	inline typename Size< MaxRepeat<TSTree> >::Type 
	length(MaxRepeat<TSTree> const &repeat) {
		return repeat.it.countRepeats();
	}
/*
	template <typename TSTree>
	inline typename Iterator< MaxRepeat<TSTree> >::Type 
	begin(MaxRepeat<TSTree> &repeat) {
		return Iterator< MaxRepeat<TSTree> >::Type(repeat.it);
	}

	template <typename TSTree>
	inline typename Iterator< MaxRepeat<TSTree> const >::Type 
	begin(MaxRepeat<TSTree> const &repeat) {
		return Iterator< MaxRepeat<TSTree> >::Type(repeat.it);
	}
*/


	template <typename TSTree>
	class Iter< MaxRepeat<TSTree>, MaxRepeatOccurences > {
	public:

		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef	Pair<TSize>						TPair;

		typedef _FractionCompound<TValue, TSize>	TFractionCompound;
		typedef typename TFractionCompound::TSet	TSet;
		typedef typename Iterator<TSet const>::Type	TSetIterator;

		TSize			child_ptr, parent_ptr;
		TSetIterator	child_fraction, child_end;
		TSetIterator	parent_fraction, parent_end;
		bool			_atEnd;
		TPair			tmp;
		bool			leftmost_child, leftmost_parent;
		
		Iter<TSTree, VSTree<BottomUp<MaxRepeats> > > const *maxIt;

		inline Iter(Iter<TSTree, VSTree<BottomUp<MaxRepeats> > > const &_maxIt):
			maxIt(&_maxIt)
		{
			_init();
		}
		
		inline bool _innerStep() {
			if (_isSizeInval(child_ptr = maxIt->posList[child_ptr])) {
				if (_isSizeInval(parent_ptr = maxIt->posList[parent_ptr])) return false;
				child_ptr = objectOf(child_fraction).begin;
			}
			return true;
		}

		inline void _firstParentFraction() {
			TFractionCompound const &parent = topPrev(maxIt->setStack);

			parent_fraction	= begin(parent.set);
			parent_end		= end(parent.set);

			if (parent_fraction != parent_end) {
				leftmost_parent = false;
				parent_ptr = objectOf(parent_fraction).begin;
			} else {
				leftmost_parent = true;
				parent_ptr = parent.leftmost.begin;
			}
		}

		inline void _firstChildFraction() {
			TFractionCompound const &child = top(maxIt->setStack);

			child_fraction	= begin(child.set);
			child_end		= end(child.set);

			if (child_fraction != child_end) {
				leftmost_child = false;
				child_ptr = objectOf(child_fraction).begin;
			} else {
				leftmost_child = true;
				child_ptr = child.leftmost.begin;
			}
		}

		inline bool _nextParentFraction() {
			if (leftmost_parent)
				return false;

			if (++parent_fraction == parent_end) {
				if (topPrev(maxIt->setStack).leftmost.size > 0) {
					leftmost_parent = true;
					parent_ptr = topPrev(maxIt->setStack).leftmost.begin;
				} else
					return false;
			} else
				parent_ptr = objectOf(parent_fraction).begin;

			return true;
		}

		inline bool _nextChildFraction() {
			if (leftmost_child)
				return false;

			if (++child_fraction == child_end) {
				if (top(maxIt->setStack).leftmost.size > 0) {
					leftmost_child = true;
					child_ptr = top(maxIt->setStack).leftmost.begin;
				} else
					return false;
			} else
				child_ptr = objectOf(child_fraction).begin;

			return true;
		}

		inline bool _outerStep() {
			do {
				if (!_nextChildFraction()) {
					_firstChildFraction();
					if (!_nextParentFraction()) {
						_atEnd = true;
						return false;
					}
				}
				if (leftmost_child || leftmost_parent) break;
			} while (keyOf(child_fraction) == keyOf(parent_fraction));		// ignore occurences with equal bwt entries
			return true;
		}

		inline void _init() 
		{
			if (length(maxIt->setStack) < 2) {
				_atEnd = true;
				return;
			}

			_firstChildFraction();
			_firstParentFraction();

			if (!leftmost_child && !leftmost_parent &&
				(keyOf(child_fraction) == keyOf(parent_fraction)))
				_atEnd = !_outerStep();
			else
				_atEnd = false;

			if (!_atEnd) {
				tmp.i1 = saAt(parent_ptr, container(*maxIt));
				tmp.i2 = saAt(child_ptr, container(*maxIt));
			}
		}
	};


	template < typename TRepeat >
	inline typename Value< Iter<TRepeat, MaxRepeatOccurences> >::Type &
	value(Iter<TRepeat, MaxRepeatOccurences> const &it)  {
		return it.tmp;
	}

	template < typename TRepeat >
	inline typename Value< Iter<TRepeat, MaxRepeatOccurences> >::Type &
	value(Iter<TRepeat, MaxRepeatOccurences> &it)  {
		return it.tmp;
	}

	template < typename TRepeat >
	inline Iter<TRepeat, MaxRepeatOccurences> &
	goNext(Iter<TRepeat, MaxRepeatOccurences> &it)  {
		if (it._innerStep()) {
			it.tmp.i1 = saAt(it.parent_ptr, container(*it.maxIt));
			it.tmp.i2 = saAt(it.child_ptr, container(*it.maxIt));
			return it;
		}
		if (it._outerStep()) {
			it.tmp.i1 = saAt(it.parent_ptr, container(*it.maxIt));
			it.tmp.i2 = saAt(it.child_ptr, container(*it.maxIt));
		}
		return it;
	}

	template < typename TRepeat >
	inline bool atEnd(Iter<TRepeat, MaxRepeatOccurences> const &it) {
		return it._atEnd;
	}

	template < typename TRepeat >
	inline bool atEnd(Iter<TRepeat, MaxRepeatOccurences> &it) {
		return it._atEnd;
	}


/*
	template <typename TSTree>
	struct ItValue< Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > > {
		typedef MaxRepeat<TSTree> Type;
	};
*/
	template <typename TSTree>
	struct Size< Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > > {
		typedef typename Size<TSTree>::Type Type;
	};



	template <typename TSTree>
	struct Iterator< MaxRepeat<TSTree> > {
		typedef Iter<MaxRepeat<TSTree>, MaxRepeatOccurences> Type;
	};

	template <typename TSTree>
	struct Size< Iter<MaxRepeat<TSTree>, MaxRepeatOccurences> > {
		typedef typename Size<TSTree>::Type Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// Iterator wrappers
	//////////////////////////////////////////////////////////////////////////////

	template <typename TObject, typename TSpec>
	struct Iterator< TObject, BottomUp<TSpec> > {
		typedef Iter< TObject, VSTree< BottomUp<TSpec> > > Type;
	};

	template <typename TObject, typename TSpec>
	struct Iterator< TObject, TopDown<TSpec> > {
		typedef Iter< TObject, VSTree< TopDown<TSpec> > > Type;
	};

	template <typename TObject>
	struct Iterator< TObject, MaxRepeats > {
		typedef Iter< TObject, VSTree< BottomUp<MaxRepeats> > > Type;
	};

	template <typename TObject>
	struct Iterator< TObject, SuperMaxRepeats > {
		typedef Iter< TObject, VSTree< BottomUp<SuperMaxRepeats> > > Type;
	};

	template <typename TObject>
	struct Iterator< TObject, SuperMaxRepeatsFast > {
		typedef Iter< TObject, VSTree< BottomUp<SuperMaxRepeatsFast> > > Type;
	};

	template <typename TObject>
	struct Iterator< TObject, MUMs > {
		typedef Iter< TObject, VSTree< BottomUp<MUMs> > > Type;
	};



//}

}

#endif
