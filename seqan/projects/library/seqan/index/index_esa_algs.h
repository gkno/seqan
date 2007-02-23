/*
 *  index_esa_algs.h
 *  genindex
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
			minLength(0),
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
			goNext(it, Postorder());
		} while (!atEnd(it) && 
			     !(	(countOccurences(it) == it.seqCount) && 
					(repLength(it) >= it.minLength) &&
					isUnique(it) && 
					isLeftMaximal(it)) );
	}


	//////////////////////////////////////////////////////////////////////////////
	// super-maximal repeats - suffix tree version
	//////////////////////////////////////////////////////////////////////////////

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
			minLength(0)
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
			goNext(it, Postorder());
		} while (!atEnd(it) && 
			     !(	childrenAreLeaves(it) && 
					(repLength(it) >= it.minLength) &&
					!isPartiallyLeftExtensible(it, it.charSet)) );
	}
		

	//////////////////////////////////////////////////////////////////////////////
	// supermaximal repeats - specialized for Enhanced Suffix Arrays
	//////////////////////////////////////////////////////////////////////////////

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
		TSize		lValue;			// current l-value of interval
		bool 		rising;			// is the left interval border valid
		typename Set<TValue>::Type	charSet;
//____________________________________________________________________________

		Iter(TIndex &__index):
			TBase(__index),
			minLength(0),
			lValue(0),
			rising(true)
		{
			this->_range = Pair<TSize>(0,0);
			indexRequire(__index, ESA_BWT());
			lIter = begin(indexLCP(this->_index));
			lEnd  = end(indexLCP(this->_index));
			goNext(*this);
		}

		Iter(TIndex &__index, TSize _minLength):
			TBase(__index),
			minLength(_minLength),
			lValue(0),
			rising(true)
		{
			this->_range = Pair<TSize>(0,0);
			indexRequire(__index, ESA_BWT());
			lIter = begin(indexLCP(this->_index));
			lEnd  = end(indexLCP(this->_index));
			goNext(*this);
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			lIter(_origin.lIter),
			lEnd(_origin.lEnd),
			lValue(_origin.lValue),
			rising(_origin.rising) {}
	};

	template < typename TText, typename TSpec >
	inline void goNext(Iter< Index<TText, Index_ESA<TSpec> >, VSTree< BottomUp<SuperMaxRepeatsFast> > > &it) 
	{
		typedef Index<TText, Index_ESA<TSpec> >		TIndex;
		typename Size<TIndex>::Type					lcp;

		while (it.lIter != it.lEnd) {
			lcp = *it.lIter;

			if (lcp < it.lValue) {
				if (it.rising) {
					if (it.lValue > it.minLength) {
						it._lcp = it.lValue;
						++it._range.i2;
						++it.lIter;
						it.lValue = lcp;
						if (!isPartiallyLeftExtensible(it, it.charSet)) return;
						continue;
					}
					it.rising = false;
				}
			} else
			if (lcp > it.lValue) {
				it._range.i1 = it._range.i2;
				it.rising = true;
			}

			++it._range.i2;
			++it.lIter;
			it.lValue = lcp;
		}
		it._range.i2 = 0;
		return;
	}

	
	//////////////////////////////////////////////////////////////////////////////
	// maximal repeats - suffix tree version
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSize>
	struct _FractionHeader {
		TSize	begin, end;
		TSize	size;
		_FractionHeader() {}
		_FractionHeader(TSize _begin, TSize _end, TSize _size):
			begin(_begin), end(_end), size(_size) {}
	};

	template < typename TSTree >
	class Iter< TSTree, VSTree< BottomUp<MaxRepeats> > >:
		public Iter< TSTree, VSTree< BottomUp<> > > 
	{
	public:
		typedef Iter< TSTree, VSTree< BottomUp<> > >	TBase;
		typedef typename Size<TSTree>::Type				TSize;
		typedef typename Value<TSTree>::Type			TValue;

		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;	// TFraction = (c,(begin,end))	c..char, begin/end indices in posList
		typedef typename Set<TFraction>::Type	TSet;
		typedef typename Iterator<TSet>::Type	TSetIterator;
		typedef String<TSet, Block<> >			TSetStack;
		typedef String<TSize>					TPositionList;
//____________________________________________________________________________

		TSize			minLength;
		TSetStack		setStack;
		TPositionList	posList;	// this list is indexed just as SA is and contains the next entry's index
//____________________________________________________________________________

		Iter(TSTree &_tree):
			TBase(_tree),
			minLength(0)
		{
			indexRequire(_tree, ESA_BWT());
			push(setStack, TSet());
			resize(posList, length(_tree));
			goNext(*this);
		}

		Iter(TSTree &_tree, TSize _minLength):
			TBase(_tree),
			minLength(_minLength)
		{
			indexRequire(_tree, ESA_BWT());
			push(setStack, TSet());
			resize(posList, length(_tree));
			goNext(*this);
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			setStack(_origin.setStack),
			posList(_origin.posList) {}

		inline bool hasRepeats() 
		{
			if (length(setStack) < 2) return false;

			TSet &child  = top(setStack);
			TSet &parent = *(end(setStack) - 2);

			TSize cs = length(child), ps = length(parent);

			if (cs == 0 || ps == 0) return false;
			if (cs  > 1 || ps  > 1) return true;

			return keyOf(begin(child)) == keyOf(begin(parent));
		}

		inline TSize countRepeats() 
		{
			if (length(setStack) < 2) return 0;

			TSet &child  = top(setStack);
			TSet &parent = *(end(setStack) - 2);

			TSetIterator child_fraction		= begin(child);
			TSetIterator child_end			= end(child);
			TSetIterator parent_fraction	= begin(parent);
			TSetIterator parent_end			= end(parent);

			TSize sum = 0;
			for(; child_fraction != child_end; ++child_fraction)
				for(; parent_fraction != parent_end; ++parent_fraction) {
					if (keyOf(child_fraction) != keyOf(parent_fraction))
						sum += (*child_fraction).size * (*parent_fraction).size;
				}
			return sum;
		}
	};

	template < typename TSTree >
	inline typename VertexDescriptor<TSTree>::Type 
	value(Iter< TSTree, VSTree< BottomUp<MaxRepeats> > > const &it) {
		return typename VertexDescriptor<TSTree>::Type ((*(end(it.history) - 2)).i1, it._range.i2);
	}

	template < typename TSTree >
	inline typename Size<TSTree>::Type repLength(Iter< TSTree, VSTree< BottomUp<MaxRepeats> > > const &it) {
		return top(it.history).i2;
	}

	// merge bwt partitions of current and parent node
	template < typename TSTree >
	inline void maxRepeatMerge(Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > &it) 
	{
		if (length(it.setStack) < 2) return;

		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;
		typedef typename Set<TFraction>::Type	TSet;
		typedef typename Iterator<TSet>::Type	TSetIterator;

		TSet &child  = top(it.setStack);
		TSet &parent = *(end(it.setStack) - 2);

		TSetIterator _end = end(child);
		for(TSetIterator i = begin(child); i != _end; ++i) {
			if (in(keyOf(i), parent)) {	// append child fraction to parent's fraction
				TFractionHeader &parent_header = objectOf(find(keyOf(i), parent));
				TFractionHeader &child_header = objectOf(i);
				it.posList[parent_header.end] = child_header.begin;
				parent_header.end = child_header.end;
				parent_header.size += child_header.size;
			} else
				insert(TFraction(keyOf(i), objectOf(i)), parent);	// insert child fraction in parent's set
		}
	}

	// maximal repeat push/pop handlers of lcp-dfs-traversal
	template < typename TSTree >
	inline void _postorderPop(Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > &it) 
	{
		typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
		_postorderPop((TBase&)it);
	
		maxRepeatMerge(it);
		pop(it.setStack);
	}

	template < typename TSTree, typename TElement >
	inline void _postorderPush(Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > &it, TElement const &e) 
	{
		typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
		_postorderPush((TBase&)it, e);
	
		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;
		typedef typename Set<TFraction>::Type	TSet;

		push(it.setStack, TSet());
	}

	template < typename TSTree >
	inline void _postorderLeaf(Iter<TSTree, VSTree< BottomUp<MaxRepeats> > > &it) 
	{
		typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
		_postorderLeaf((TBase&)it);

		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;

		TSize index = _dfsRange(it).i1;
		if (!posAtFirstLocal(saAt(index, container(it)), stringSetLimits(it)))
			insert(TFraction(bwtAt(index, container(it)), TFractionHeader(index, index, 1)), top(it.setStack));
		_setSizeInval(it.posList[index]);
	}

	template < typename TSTree >
	inline void goNext(Iter< TSTree, VSTree< BottomUp<MaxRepeats> > > &it) {
		do {
			goNext(it, Postorder());
		} while (!eof(it) && !it.hasRepeats());
	}





	//////////////////////////////////////////////////////////////////////////////
	// maximal repeat representation
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSTree>
	struct MaxRepeat {
		Iter< TSTree, VSTree<MaxRepeats> > &it;
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
	inline typename Size< MaxRepeat<TSTree> >::Type length(MaxRepeat<TSTree> &repeat) {
		return repeat.it.countRepeats();
	}



	template <typename TSTree>
	class Iter< MaxRepeat<TSTree>, MaxRepeatOccurences> {
	public:

		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef	Pair<TSize>						TPair;
		typedef _FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;	// TFraction = (c,(begin,end))	c..char, begin/end indices in posList

		typedef typename Set<TFraction>::Type	TSet;
		typedef typename Iterator<TSet>::Type	TSetIterator;

		TSize			child_ptr, parent_ptr;
		TSetIterator	child_fraction, child_end;
		TSetIterator	parent_fraction, parent_end;
		
		Iter<TSTree, VSTree<MaxRepeats> >	*repeat;

		Iter(Iter<TSTree, VSTree<MaxRepeats> > &_repeat):
			repeat(_repeat)
		{
			_update();
		}
		
		inline TPair operator*() const {
			return TPair(saAt(child_ptr, container(*repeat)), saAt(parent_ptr, container(*repeat)));
		}

		inline Iter operator++() {
			TSize next = repeat->posList[child_ptr];

			do {
				if (!_innerStep()) {
					_outerStep();
				}
			} while (parent_fraction != parent_end);
			return true;

			return *this;
		}

	private:

		inline bool _innerStep() {
			if (_isSizeInval(child_ptr = repeat->posList[child_ptr])) {
				if (_isSizeInval(parent_ptr = repeat->posList[parent_ptr])) return false;
				child_ptr = begin(repeat->_set);
			}
			return true;
		}

		inline bool _outerStep() {
			do {
				if (++child_fraction == child_end) {
					if (++parent_fraction == parent_end) return false;
					child_fraction = begin(top(repeat->setStack));
				}
			} while (keyOf(child_fraction) == keyOf(parent_fraction));		// ignore occurences with equal bwt entries
			child_ptr = child_fraction.begin;
			parent_ptr = parent_fraction.begin;
			return true;
		}

		inline void _update() 
		{
			if (length(repeat->setStack) < 2) return false;

			TSet &child  = top(repeat->setStack);
			TSet &parent = *(end(repeat->setStack) - 2);

			TSetIterator child_fraction		= begin(child);
			TSetIterator child_end			= end(child);
			TSetIterator parent_fraction	= begin(parent);
			TSetIterator parent_end			= end(parent);

			if (keyOf(child_fraction) == keyOf(parent_fraction))
				_outerStep();
			else {
				child_ptr = child_fraction.begin;
				parent_ptr = parent_fraction.begin;
			}
		}
	};
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
