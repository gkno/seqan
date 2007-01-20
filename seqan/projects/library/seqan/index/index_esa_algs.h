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
	// supermaximal repeats - enhanced suffix array version (needs bwt)

	template < typename TObject, typename TSpec >
	class Iter< Index<TObject, TSpec>, VSTree< BottomUp<SuperMaxRepeats> > > {
	public:

		typedef typename Value< Index<TObject, TSpec> >::Type	TValue;
		typedef typename Size< Index<TObject, TSpec> >::Type	TSize;
		typedef	Pair<TSize>								TPair;
		typedef ::std::stack<TPair>						TStack;
		typedef Iter									iterator;

		typedef typename Iterator<typename Fibre<Index<TObject, TSpec>, ESA_LCP>::Type const>::Type			LCPIter;
		typedef typename Iterator<typename Fibre<Index<TObject, TSpec>, ESA_BWT>::Type::TText const>::Type	BWTIter;

		Index<TObject, TSpec> const	&_index;	// container of all necessary tables
		TPair						_range;		// current interval in suffix array

		Iter(Index<TObject, TSpec> const &__index):
			_index(__index),
			_range(TPair(0,0)),
			lastLcp(0),
			leftMax(true),
			unique(true),
			minLength(0)
		{
			lIter = begin(indexLCP(_index));
			bIter = begin(indexBWT(_index).tab);
			undefined = indexBWT(_index).undefined;
			++*this;
		}

		Iter(Index<TObject, TSpec> const &__index, TSize _minLength):
			_index(__index),
			_range(TPair(0,0)),
			lastLcp(0),
			leftMax(true),
			unique(true),
			minLength(_minLength)
		{
			lIter = begin(getFibre(_index, ESA_LCP()));
			lEnd = end(getFibre(_index, ESA_LCP()));
			bIter = begin(getFibre(_index, ESA_BWT()).tab);
			undefined = getFibre(_index, ESA_BWT()).undefined;
			++*this;
		}

		Iter(Iter const &_origin):
			_index(container(_origin)),
			_range(_origin._range),
			lastLcp(_origin.lastLcp),
			leftMax(_origin.leftMax) {}

		template < typename _TText, typename _TSpec >
		friend Iter< Index<_TText, _TSpec>, VSTree< BottomUp<SuperMaxRepeats> > >& operator++(Iter< Index<_TText, _TSpec>, VSTree< BottomUp<SuperMaxRepeats> > > &it);

	private:
		TSize		minLength;
		LCPIter		lIter, lEnd;
		BWTIter		bIter;
		TSize		lastLcp;
		TSize		undefined;
		bool 		leftMax;	// is the current l-interval left-maximal?
		bool		unique;		// is the maximal repeat left-maximal and not part of another (== supermaximal)
		typename Set<TValue>::Type charset;
	};

	template < typename TIndex, class TSpec >
	inline Pair< typename Size<TIndex>::Type > operator*(Iter< TIndex, VSTree< BottomUp<SuperMaxRepeats> > > const &it) {
		return it.value;
	}

	template < typename TText, typename TSpec >
	finline Iter< Index<TText, TSpec>, VSTree< BottomUp<SuperMaxRepeats> > >& operator++(Iter< Index<TText, TSpec>, VSTree< BottomUp<SuperMaxRepeats> > > &it) {
		typename Size< Index<TText, TSpec> >::Type lcp;

		while (it.lIter != it.lEnd) {
			lcp = *it.lIter;

			if (lcp < it.lastLcp || lcp < it.minLength) {
				if (it.leftMax) {
					if (it._range.i2 != it.undefined && it.unique)
						if (in(*it.bIter, it.charset))
							it.unique = false;
						else
							insert(*it.bIter, it.charset);

					it.leftMax = false;
					++it._range.i2;
					++it.lIter;
					++it.bIter;
					it.lastLcp = lcp;
					if (it.unique) return it;
					continue;
				}
			} else
			if (lcp > it.lastLcp) {
				it._range.i1 = it._range.i2;
				it.leftMax = true;
				clear(it.charset);
				if (it._range.i2 != it.undefined)
					insert(*it.bIter, it.charset);
				it.unique = true;
			} else
				if (it._range.i2 != it.undefined && it.unique)
					if (in(*it.bIter, it.charset))
						it.unique = false;
					else
						insert(*it.bIter, it.charset);

			++it._range.i2;
			++it.lIter;
			++it.bIter;
			it.lastLcp = lcp;
		}
		it._range.i2 = 0;
		return it;
	}

/*

	//////////////////////////////////////////////////////////////////////////////
	// supermaximal repeats - suffix tree version

	template < typename TTree >
	class Iter< TTree, VSTree< BottomUp<SuperMaxRepeats> > > {
	public:

		typedef typename Value<TTree>::Type	TValue;
		typedef typename Size<TTree>::Type	TSize;
		typedef	Pair<TSize>					TPair;
		typedef ::std::stack<TPair>			TStack;
		typedef Iter						iterator;

		TTree		&_index;	// container of all necessary tables
		TStack		history;	// contains all left borders of current l-intervals (== left borders of history intervals)
		TPair		_range;		// current interval in suffix array
		TSize		_lcp;		// current l-value

		Iter(TTree const &__index):
			_index(__index),
			_range(TPair(0,1)),
			minLength(0)
		{
			_lcp = length(_index) - saAt(0, _index);
			history.push(TPair(0,0));
			++*this;
		}

		Iter(TTree const &__index, TSize _minLength):
			_index(__index),
			_range(TPair(0,1)),
			minLength(_minLength)
		{
			_lcp = length(_index) - saAt(0, _index);
			history.push(TPair(0,0));
			++*this;
		}

		Iter(Iter const &_origin):
			_index(container(_origin)),
			_range(_dfsRange(_origin)),
			_lcp(_dfsLCP(_origin)),
			history(_origin.history) {}

		template < typename TSpec >
		friend Iter<TTree, VSTree<SuperMaxRepeats> >& operator++(Iter<TTree, VSTree<SuperMaxRepeats> > &it);

	public:
		TSize						minLength;
		typename Set<TValue>::Type	charset;
	};

	template < typename TTree >
	inline Iter<TTree, VSTree< BottomUp<SuperMaxRepeats> > >& operator++(Iter<TTree, VSTree< BottomUp<SuperMaxRepeats> > > &it) {
		typedef typename Fibre<TTree, ESA_RawText>::Type	TText;
		typedef typename Value<TTree>::Type					TValue;
		typedef typename Size<TTree>::Type					TSize;

		typedef typename Infix<typename Fibre<TTree, ESA_SA>::Type>::Type	TOccs;
		typedef typename Iterator<TOccs>::Type								TOccIter;

		TOccIter i, e;
		TValue c;
		goNextPostorder(it, it.history);
		while (!eof(it)) {
			if (childrenAreLeaves(it)) {
				TOccs occs = getOccurences(it);
				clear(it.charset);
				e = end(occs);
				for(i = begin(occs); i != e; ++i) {
					if (bwtValidAt(*i, container(it))) {
						c = bwtAt(*i, container(it));
						if (in(c, it.charset))
							break;
					}
					insert(c, it.charset);
				}
				if (i == e && lcpAt(occs[0], index) >= it.minLength) return it;
			}
			goNextPostorder(it, it.history);
		}
		return it;
	}

*/

	//////////////////////////////////////////////////////////////////////////////
	// maximal repeats - suffix tree version

	template <typename TSize>
	struct FractionHeader {
		TSize	begin, end;
		TSize	size;
		FractionHeader() {}
		FractionHeader(TSize _begin, TSize _end, TSize _size):
			begin(_begin), end(_end), size(_size) {}
	};

	template < typename TTree >
	class Iter< TTree, VSTree< BottomUp<MaxRepeats> > > {
	public:

		typedef typename Value<TTree>::Type		TValue;
		typedef typename Size<TTree>::Type		TSize;
		typedef	Pair<TSize>						TPair;
		typedef FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;	// TFraction = (c,(begin,end))	c..char, begin/end indices in nextPos

		typedef typename Set<TFraction>::Type	TSet;
		typedef typename Iterator<TSet>::Type	TSetIterator;
		typedef Triple<TSize, TSize, TSet*>		TStackEntry;
		typedef ::std::vector<TStackEntry>		TStack;
		typedef Iter							iterator;

		typedef String<TSize>					TPositionList;

		typedef Allocator< SinglePool<sizeof(TSet)> >						TAllocator;
		typedef typename Infix<typename Fibre<TTree, ESA_SA>::Type>::Type	TOccs;
		typedef typename Iterator<TOccs>::Type								TOccIter;

		TTree const	&_index;		// container of all necessary tables
		TStack		history;		// contains all left borders of current l-intervals (== left borders of history intervals)
		TPair		_range;			// current interval (bottom-up traversal)
		TSize		_lcp;			// current l-value
		TSet		*_set;			// current set of repeat occurences

		TAllocator	alloc;

		Iter(TTree const &__index):
			_index(__index),
			_range(TPair(0,1)),
			minLength(0)
		{
			_lcp = length(_index) - saAt(0, _index);
			history.push_back(TStackEntry(0,0, NULL));
			allocate(alloc, _atomicTop(history).i3, 1, *this);
			++*this;
		}

		Iter(TTree const &__index, TSize _minLength):
			_index(__index),
			_range(TPair(0,1)),
			minLength(_minLength)
		{
			_lcp = length(_index) - saAt(0, _index);
			_atomicPush(history, TStackEntry(0,0, NULL));
			allocate(alloc, _atomicTop(history).i3, 1, *this);
			++*this;
		}

		Iter(Iter const &_origin):
			_index(container(_origin)),
			_range(_dfsRange(_origin)),
			_lcp(_dfsLCP(_origin)),
			history(_origin.history) {}

		template < typename TSpec >
		friend Iter<TTree, VSTree<MaxRepeats> >& operator++(Iter<TTree, VSTree<MaxRepeats> > &it);

		inline bool noRepeats() {
			TStackEntry &child = history.back();
			TStackEntry &parent = *(history.end() - 2);

			TSize cs = length(child.i3), ps = length(parent.i3);

			if (cs == 0 || ps == 0) return true;
			if (cs > 1 || ps > 1) return false;

			return keyOf(begin(child.i3)) == keyOf(begin(parent.i3));
		}

		inline TSize countRepeats() {
			TStackEntry &child = history.back();
			TStackEntry &parent = *(history.end() - 2);

			TSetIterator child_fraction		= begin(child.i3);
			TSetIterator child_end			= end(child.i3);
			TSetIterator parent_fraction	= begin(parent.i3);
			TSetIterator parent_end			= end(parent.i3);

			TSize sum = 0;
			for(; child_fraction != child_end; ++child_fraction)
				for(; parent_fraction != parent_end; ++parent_fraction) {
					if (keyOf(child_fraction) != keyOf(parent_fraction))
						sum += (*child_fraction).size * (*parent_fraction).size;
				}
			return sum;
		}

	public:
		TSize			minLength;
		TPositionList	nextPos;


		static inline void beforeUp(iterator &it) {
			TSet &parentSet = it.history.back().i3;

			TSetIterator _end = end(it._set);
			for(TSetIterator i = begin(it._set); i != _end; ++i) {
				if (in((*i).i1, parentSet)) {
					TFractionHeader &parent_header = objectOf(find(keyOf(i), parentSet));
					TFractionHeader &child_header = objectOf(i);
					it.nextPos[parent_header.end] = child_header.begin;
					parent_header.end = child_header.end;
					parent_header.size += child_header.size;
				} else {
					TSize index = getOccurence(it);
					insert(TFraction(keyOf(i), TFractionHeader(index, index, 1)), parentSet);	// create new fraction in parent node
					_setSizeInval(it.nextPos[index]);											// begin und end pointers are the same (== lastIndex)
				}
			}
		}
		static inline void afterUp(iterator const &) {}
		static inline void beforeDown(iterator const &) {}

		static inline void afterDown(iterator &it) {
			if (isLeaf(_dfsRange(it)) && bwtValidAt(_dfsRange(it).i1, container(it))) 
				insert(TFraction(bwtAt(_dfsRange(it).i1, container(it)), TFractionHeader(_dfsRange(it).i1, _dfsRange(it).i1, 1)), it._set); 
		}
	};

	template < typename TIndex >
	inline Pair<typename Size<TIndex>::Type> value(Iter< TIndex, VSTree< BottomUp<MaxRepeats> > > const &it) {
		return Pair<typename Size<TIndex>::Type> ((*(it.history.end() - 2)).i1, it._range.i2);
	}

	template < typename TIndex >
	inline typename Size<TIndex>::Type length(Iter< TIndex, VSTree< BottomUp<MaxRepeats> > > const &it) {
		return it.history.back().i2;
	}



	// merge bwt partitions of current and father node
	template < typename TTree >
	inline void maxRepeatMerge(Iter<TTree, VSTree< BottomUp<MaxRepeats> > > &it) {
		typedef typename Value<TTree>::Type		TValue;
		typedef typename Size<TTree>::Type		TSize;
		typedef FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;	// TFraction = (c,(begin,end))	c..char, begin/end indices in nextPos
		typedef typename Set<TFraction>::Type	TSet;
		typedef typename Iterator<TSet>::Type	TSetIterator;

		TSet &parentSet = it.history.back().i3;

		TSetIterator _end = end(it._set);
		for(TSetIterator i = begin(it._set); i != _end; ++i) {
			if (in((*i).i1, parentSet)) {
				TFractionHeader &parent_header = objectOf(find(keyOf(i), parentSet));
				TFractionHeader &child_header = objectOf(i);
				it.nextPos[parent_header.end] = child_header.begin;
				parent_header.end = child_header.end;
				parent_header.size += child_header.size;
			} else {
				TSize index = getOccurence(it);
				insert(TFraction(keyOf(i), TFractionHeader(index, index, 1)), parentSet);	// create new fraction in parent node
				_setSizeInval(it.nextPos[index]);											// begin und end pointers are the same (== lastIndex)
			}
		}
	}


	// maximal repeat push/pop handlers of lcp-dfs-traversal
	template < typename TIndex >
	inline void _postorderPop(Iter<TIndex, VSTree< BottomUp<MaxRepeats> > > &it) {
        _dfsRange(it).i1 = _atomicTop(it.history).i1;
		_dfsLCP(it) = _atomicTop(it.history).i2;
		
		if (it._set) {
			maxRepeatMerge(it);
			deallocate(it.alloc, it._set, 1, it);
		}
		it._set = _atomicTop(it.history)._set;
		_atomicPop(it.history);
	}

	template < typename TIndex, typename TElement >
	inline void _postorderPush(Iter<TIndex, VSTree< BottomUp<MaxRepeats> > > &it, TElement &e) {
		allocate(it.alloc, e._set, 1, it);
		_atomicPush(it.history, e);
	}

	template < typename TIndex, typename TElement >
	inline void _postorderLeaf(Iter<TIndex, VSTree< BottomUp<MaxRepeats> > > &it) {
		typedef typename Value<TIndex>::Type	TValue;
		typedef typename Size<TIndex>::Type	TSize;
		typedef FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;	// TFraction = (c,(begin,end))	c..char, begin/end indices in nextPos

		_setSizeInval(_dfsLCP(it));
		//length(container(it)) - saAt(_dfsRange(it).i1, container(it));
		if (it._set) {
			deallocate(it.alloc, it._set, 1, it);
			it._set = NULL;
		}

		if (bwtValidAt(_dfsRange(it).i1, container(it)))
			insert(TFraction(bwtAt(_dfsRange(it).i1, container(it)), TFractionHeader(_dfsRange(it).i1, _dfsRange(it).i1, 1)), it._set); 
	}



	template < typename TTree >
	inline Iter<TTree, VSTree< BottomUp<MaxRepeats> > >& operator++(Iter<TTree, VSTree< BottomUp<MaxRepeats> > > &it) {
		do {
			goNextPostorder(it, it.history, it);
		} while (!eof(it) && it.noRepeats());
		return it;
	}



	template <typename TTree>
	struct MaxRepeat {
		Iter< TTree, VSTree<MaxRepeats> > &it;
	};

	template <typename TTree>
	struct Value< MaxRepeat<TTree> > {
		typedef Pair< typename Size<TTree>::Type > Type;
	};

	template <typename TTree>
	struct Size< MaxRepeat<TTree> > {
		typedef typename Size<TTree>::Type Type;
	};


	template <typename TTree>
	inline typename Size< MaxRepeat<TTree> >::Type length(MaxRepeat<TTree> &repeat) {
		return repeat.it.countRepeats();
	}



	template <typename TTree>
	class Iter< MaxRepeat<TTree>, MaxRepeatOccurences> {
	public:

		typedef typename Value<TTree>::Type		TValue;
		typedef typename Size<TTree>::Type		TSize;
		typedef	Pair<TSize>						TPair;
		typedef FractionHeader<TSize>			TFractionHeader;
		typedef Pair<TValue, TFractionHeader>	TFraction;	// TFraction = (c,(begin,end))	c..char, begin/end indices in nextPos

		typedef typename Set<TFraction>::Type	TSet;
		typedef typename Iterator<TSet>::Type	TSetIterator;
		typedef Triple<TSize, TSize, TSet>		TStackEntry;

		TSize			child_ptr, parent_ptr;
		TSetIterator	child_fraction, child_end;
		TSetIterator	parent_fraction, parent_end;
		
		Iter<TTree, VSTree<MaxRepeats> >	*repeat;

		Iter(Iter<TTree, VSTree<MaxRepeats> > &_repeat):
			repeat(_repeat)
		{
			_update();
		}
		
		inline TPair operator*() const {
			return TPair(child_ptr, parent_ptr);
		}

		inline Iter operator++() {
			TSize next = repeat->nextPos[child_ptr];
			//if (next != 

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
			if (_isSizeInval(child_ptr = repeat->nextPos[child_ptr])) {
				if (_isSizeInval(parent_ptr = repeat->nextPos[parent_ptr])) return false;
				child_ptr = begin(repeat->_set);
			}
			return true;
		}

		inline bool _outerStep() {
			do {
				if (++child_fraction == child_end) {
					if (++parent_fraction == parent_end) return false;
					child_fraction = begin(repeat->history.back().i3);
				}
			} while (keyOf(child_fraction) == keyOf(parent_fraction));		// ignore occurences with equal bwt entries
			child_ptr = child_fraction.begin;
			parent_ptr = parent_fraction.begin;
			return true;
		}

		inline void _update() {
			TStackEntry &child = repeat->history.back();
			TStackEntry &parent = *(repeat->history.end() - 2);

			TSetIterator child_fraction		= begin(child.i3);
			TSetIterator child_end			= end(child.i3);
			TSetIterator parent_fraction	= begin(parent.i3);
			TSetIterator parent_end			= end(parent.i3);

			if (keyOf(child_fraction) == keyOf(parent_fraction))
				_outerStep();
			else {
				child_ptr = child_fraction.begin;
				parent_ptr = parent_fraction.begin;
			}
		}
	};
/*
	template <typename TTree>
	struct ItValue< Iter<TTree, VSTree< BottomUp<MaxRepeats> > > > {
		typedef MaxRepeat<TTree> Type;
	};
*/
	template <typename TTree>
	struct Size< Iter<TTree, VSTree< BottomUp<MaxRepeats> > > > {
		typedef typename Size<TTree>::Type Type;
	};



	template <typename TTree>
	struct Iterator< MaxRepeat<TTree> > {
		typedef Iter<MaxRepeat<TTree>, MaxRepeatOccurences> Type;
	};

	template <typename TTree>
	struct Size< Iter<MaxRepeat<TTree>, MaxRepeatOccurences> > {
		typedef typename Size<TTree>::Type Type;
	};


//////////////////////////////////////////////////////////////////////////////
// Iterator wrappers

	template <typename TObject, typename TSpec>
	struct Iterator< Index<TObject, Index_ESA<TSpec> >, MaxRepeats > {
		typedef Iter< Index<TObject, Index_ESA<TSpec> >, VSTree< BottomUp<SuperMaxRepeats> > > Type;
	};

	template <typename TObject, typename TSpec>
	struct Iterator< Index<TObject, Index_ESA<TSpec> >, SuperMaxRepeats > {
		typedef Iter< Index<TObject, Index_ESA<TSpec> >, VSTree< BottomUp<SuperMaxRepeats> > > Type;
	};

	template <typename TObject, typename TSpec>
	struct Iterator< Index<TObject, Index_ESA<TSpec> >, MUMs > {
		typedef Iter< Index<TObject, Index_ESA<TSpec> >, VSTree< BottomUp<MUMs> > > Type;
	};



//}

}

#endif
