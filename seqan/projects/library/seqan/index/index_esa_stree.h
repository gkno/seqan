/*
 *  index_esa_stree.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_ESA_STREE_H
#define SEQAN_HEADER_INDEX_ESA_STREE_H

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Spec.VSTree Iterator:
..cat:Iterators
..summary:Abstract iterator for Suffix Trees.
..signature:Iter<TContainer, VSTree<TSpec> >
..param.TContainer:Type of the container that can be iterated.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
..param.TSpec:The specialization type.
..remarks:This iterator is a pointer to a node in the Suffix Tree (given by the Enhanced Suffix Array @Spec.Index_ESA@).
Every node can uniquely be mapped to an interval of the Suffix Array containing all suffixes of the node's subtree.
This interval is the @Function.value@ of the iterator.
*/

	template < typename TIndex, typename TSpec >
    struct Value< Iter< TIndex, VSTree<TSpec> > > {
		typedef Pair< typename Size<TIndex>::Type > Type;
	};
 
	template < typename TIndex, typename TSpec >
	struct Size< Iter< TIndex, VSTree<TSpec> > > {
		typedef typename Size<TIndex>::Type Type;
	};
 
	template < typename TIndex, typename TSpec >
	struct Position< Iter< TIndex, VSTree<TSpec> > > {
		typedef typename Position<TIndex>::Type Type;
	};
 

/**
.Spec.TopDown Iterator:
..cat:Iterators
..general:Spec.VSTree Iterator
..summary:Iterator for Suffix Trees that can go down and right beginning from the root.
..signature:Iter<TContainer, VSTree< TopDown<TSpec> > >
..param.TContainer:Type of the container that can be iterated.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
..param.TSpec:The specialization type.
*/

	template < typename TIndex, class TSpec >
	class Iter< TIndex, VSTree< TopDown<TSpec> > > {
	public:

		typedef typename Size<TIndex>::Type	TSize;
		typedef	Pair<TSize>					TPair;
		typedef Iter						iterator;

		TIndex const	&_index;	// container of all necessary tables
		TPair			_range;		// current interval in suffix array
		TSize			topRight;	// right border of parent interval (needed in goRight)

		Iter(TIndex &__index):
			_index(__index)
		{
			indexRequire(__index, ESA_SA());
			indexRequire(__index, ESA_LCP());
			indexRequire(__index, ESA_ChildTab());

			_range.i1 = 0;
			_setSizeInval(_range.i2);
		}

		Iter(Iter const &_origin):
			_index(container(_origin)),
			_range(value(_origin)),
			topRight(_origin.topRight) {}

		Iter(Iter const &_origin, TPair const &_childRange):
			_index(container(_origin)),
			_range(_childRange),
			topRight(value(_origin).i2) {}

		template <typename _TSpec>
		Iter(Iter<TIndex, VSTree< TopDown< ParentLinks<_TSpec> > > > const &_origin):
			_index(container(_origin)),
			_range(value(_origin)) 
		{
			if (!empty(_origin.history))
				topRight = top(_origin.history).i2;
		}

	};


/**
.Spec.TopDownHistory Iterator:
..cat:Iterators
..general:Spec.TopDown Iterator
..summary:Iterator for Suffix Trees that can go down, right, and up.
..signature:Iter<TContainer, VSTree< TopDown< ParentLinks<TSpec> > > >
..param.TContainer:Type of the container that can be iterated.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
..implements:Concept.Iterator
..param.TSpec:The specialization type.
*/

	template < typename TIndex, class TOrder >
	class Iter< TIndex, VSTree< TopDown< ParentLinks<TOrder> > > > {
	public:

		typedef	Pair<typename Size<TIndex>::Type>	TPair, TStackEntry;
		typedef String<TStackEntry, Block<> >		TStack;
		typedef Iter								iterator;

		TIndex	const	&_index;	// container of all necessary tables
		TPair			_range;		// current interval in suffix array
		TStack			history;	// contains all previously visited intervals (allows to go up)

		Iter(TIndex &__index):
			_index(__index)
		{
			indexRequire(__index, ESA_SA());
			indexRequire(__index, ESA_LCP());
			indexRequire(__index, ESA_ChildTab());

			_range.i1 = 0;
			_setSizeInval(_range.i2);
		}

		Iter(Iter const &_origin):
			_index(container(_origin)),
			_range(value(_origin)),
			history(_origin.history) {}

		Iter(Iter const &_origin, TPair const &_childRange):
			_index(container(_origin)),
			_range(_childRange),
			history(_origin.history)
		{
			push(history, value(_origin));
		}
	};


/**
.Spec.BottomUp Iterator:
..cat:Iterators
..general:Spec.VSTree Iterator
..summary:Iterator for an efficient depth-first search in a Suffix Tree.
..signature:Iter<TContainer, VSTree< BottomUp<TSpec> > >
..param.TContainer:Type of the container that can be iterated.
...type:Spec.Index_ESA
...metafunction:Metafunction.Container
..implements:Concept.Iterator
..param.TSpec:The specialization type.
*/

	template < typename TIndex, typename TSpec >
	class Iter< TIndex, VSTree< BottomUp<TSpec> > > {
	public:

		typedef typename Size<TIndex>::Type	TSize;
		typedef	Pair<TSize>						TPair, TStackEntry;
		typedef String<TStackEntry, Block<> >	TStack;
		typedef Iter							iterator;

		TIndex	const	&_index;	// container of all necessary tables
		TPair			_range;		// current interval in suffix array
		TSize			_lcp;		// current l-value
		TStack			history;	// contains all left borders of current l-intervals (== left borders of history intervals)

		Iter(TIndex &__index):
			_index(__index),
			_range(TPair(0,1))
		{
			indexRequire(__index, ESA_SA());
			indexRequire(__index, ESA_LCP());

			if (empty(indexSA(__index))) {
				_lcp = 0;
				return;
			}
			_lcp = suffixLength(getOccurence(*this), container(*this));
			push(history, TStackEntry(0,0));
		}

		Iter(Iter const &_origin):
			_index(container(_origin)),
			_range(_dfsRange(_origin)),
			_lcp(_dfsLCP(_origin)),
			history(_origin.history) {}
	};


	template < typename TIndex, typename TSpec >
	inline void _dumpHistoryStack(Iter<TIndex, VSTree<TSpec> > &it) {
		for(typename Size<TIndex>::Type i = 0; i < length(it.history); ++i)
			::std::cout << it.history[i] << '\t';
		::std::cout << value(it) << ::std::endl;		
	}

	// standard push/pop handlers of lcp-dfs-traversal
	template < typename TIndex, typename TSpec >
	inline void _postorderPop(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it) {
        _dfsRange(it).i1 = top(it.history).i1;
		_dfsLCP(it) = top(it.history).i2;
		pop(it.history);
//		_dumpHistoryStack(it);
	}

	template < typename TIndex, typename TSpec, typename TElement >
	inline void _postorderPush(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, TElement const &e) {
		push(it.history, e);
//		_dumpHistoryStack(it);
	}

	template < typename TIndex, typename TSpec >
	inline void _postorderLeaf(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it) {
		_setSizeInval(_dfsLCP(it));
//		_dumpHistoryStack(it);
		//length(container(it)) - saAt(_dfsRange(it).i1, container(it));
	}


	template < typename TIndex, typename TSpec >
	inline void goNextImpl(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, Postorder const) {
		TIndex const &index = container(it);
		do {
			// postorder dfs via lcp-table
			if (isRoot(it)) {
				_dfsClear(it);
				return;
			}

			typedef typename Size<TIndex>::Type TSize;
			typedef typename Iter<TIndex, VSTree< BottomUp<TSpec> > >::TStackEntry TStackEntry;
			TStackEntry	_top_ = top(it.history);
			TSize		lcp_i = lcpAt(_dfsRange(it).i2 - 1, index);

			if (lcp_i < _top_.i2) {
				_postorderPop(it);			// go up
				return;
			}

			if (lcp_i > _top_.i2) {
				_top_.i1 = _dfsRange(it).i1;
				_top_.i2 = lcp_i;
				_postorderPush(it, _top_);		// go down
			}

			// last lcp entry (== 0) causes removal of toplevel interval
			if ((_dfsRange(it).i1 = _dfsRange(it).i2++) == length(index)) {
				_postorderPop(it);
				_dfsRange(it).i2 = _dfsRange(it).i1;
			} else {
				// skip $ leafs
				if (suffixLength(saAt(_dfsRange(it).i1, index), index) == lcpAt(_dfsRange(it).i1, index))
					continue;
				_postorderLeaf(it);
			}
			return;
		} while (true);
	}

	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type repLength(Iter< TIndex, VSTree<BottomUp<TSpec> > > const &it) {
		typename Size<TIndex>::Type lcp;
		if (!_isSizeInval(lcp = it._lcp))
			return lcp;
		else
			return suffixLength(getOccurence(it), container(it));
	}


	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type repLength(Iter< TIndex, VSTree<TopDown<TSpec> > > const &it) {
		return repLength(container(it), value(it));
	}

	template < typename TIndex, class TSize >
	inline TSize repLength(TIndex const &index, Pair<TSize> const &value) {
		if (_isLeaf(value)) return suffixLength(saAt(value.i1, index), index);
		if (_isRoot(value)) return 0;

		TSize lval = _getUp(value.i2, index);
		if (!(value.i1 < lval && lval < value.i2))
			lval = _getDown(value.i1, index);
		
		return lcpAt(lval - 1, index);
	}

/**
.Function.lca:
..summary:Returns the last common ancestor of two tree nodes.
..cat:Index
..signature:bool lca(a, b, result)
..param.a:The first node.
...type:Spec.TopDownHistory Iterator
..param.b:The second node.
...type:Spec.TopDownHistory Iterator
..param.result:A reference to the resulting lca node.
...type:Spec.TopDownHistory Iterator
..returns:$false$ if the lca of $a$ and $b$ is the root node, otherwise $true$.
*/

	template < typename TIndex, class TSpec1, class TSpec2 >
	inline bool lca(
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > > &a, 
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec2> > > > &b, 
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > > &_lca)
	{
		typename Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > >::TStack::const_iterator iA;
		typename Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec2> > > >::TStack::const_iterator iB;

		typedef typename Size<TIndex>::Type TSize;

		// push current intervals
		push(a.history, value(a));
		push(b.history, value(b));

		TSize s = min(a.history.size(), b.history.size()), i0 = 0;
		
		while (s) {
			TSize m = s / 2;
			iA = a.history.begin() + i0 + m;
			iB = b.history.begin() + i0 + m;
			if ((*iA).i1 == (*iB).i1 && (*iA).i2 == (*iB).i2) {
				i0 += m + 1;
				s -= m + 1;
			} else
				s = m;
		}

		_lca.history.resize(i0);
		copy(a.history.begin(), a.history.begin() + i0, _lca.history.begin());

		// pop current intervals
		pop(a.history);
		pop(b.history);
		if (i0) {
			value(_lca) = top(_lca.history);
			pop(_lca.history);
		}

		return i0;
	}

/**
.Function.lcp:
..summary:Returns the length of the longest-common-prefix of two Suffix Tree nodes.
..cat:Index
..signature:lcp(a, b)
..param.a:The first node.
...type:Spec.TopDownHistory Iterator
..param.b:The second node.
...type:Spec.TopDownHistory Iterator
..returns:The lcp-length of $a$ and $b$.
*/

	// return the lcp of a and b by seeking the lca of them
	template < typename TIndex, class TSpec1, class TSpec2 >
	inline typename Size<TIndex>::Type lcp(
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > > &a, 
		Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec2> > > > &b) 
	{
		typename Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec1> > > >::TStack::const_iterator iA;
		typename Iter<TIndex, VSTree< TopDown< ParentLinks<TSpec2> > > >::TStack::const_iterator iB;

		typedef typename Size<TIndex>::Type TSize;

		// push current intervals
		push(a.history, value(a));
		push(b.history, value(b));

		TSize s = min(a.history.size(), b.history.size()), i0 = 0;
		
		while (s) {
			TSize m = s / 2;
			iA = a.history.begin() + i0 + m;
			iB = b.history.begin() + i0 + m;
			if ((*iA).i1 == (*iB).i1 && (*iA).i2 == (*iB).i2) {
				i0 += m + 1;
				s -= m + 1;
			} else
				s = m;
		}

		TSize _lcp = (i0 > 0)? repLength(container(a), a.history[i0 - 1]): 0;

		// pop current intervals
		pop(a.history);
		pop(b.history);

		return _lcp;
	}

///.Function.container.param.iterator.type:Spec.VSTree Iterator

	template < typename TIndex, class TSpec >
	inline TIndex const & container(Iter< TIndex, VSTree<TSpec> > const &it) { return it._index; }

	template < typename TIndex, class TSpec >
	inline TIndex & container(Iter< TIndex, VSTree<TSpec> > &it) { return const_cast<TIndex&>(it._index); }


///.Function.value.param.iterator.type:Spec.VSTree Iterator

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> const & value(Iter< TIndex, VSTree<TSpec> > &it) {
		return it._range; 
	}

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> const & value(Iter< TIndex, VSTree<TSpec> > const &it) { 
		return it._range; 
	}

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> & value(Iter< TIndex, VSTree< TopDown<TSpec> > > &it) { 
		return it._range; 
	}

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> const & value(Iter< TIndex, VSTree< TopDown<TSpec> > > const &it) { 
		return it._range; 
	}


/**
.Function.getOccurence:
..summary:Returns an occurence of the @Function.representative@ substring in the index text.
..cat:Index
..signature:getOccurence(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:A position where the @Function.representative@ of $iterator$ occurs in the text (see @Tag.ESA_Text@).
...metafunction:Metafunction.SAValue.<TIndex>
*/

	template < typename TIndex, class TSpec >
	inline typename SAValue<TIndex>::Type getOccurence(Iter< TIndex, VSTree<TSpec> > const &it) {
		return saAt(value(it).i1, container(it));
	}


/**
.Function.countOccurences:
..summary:Returns the number of occurences of @Function.representative@ in the index text.
..cat:Index
..signature:countOccurences(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:The number of positions where the @Function.representative@ of $iterator$ occurs in the text (see @Tag.ESA_Text@).
If $iterator$'s container type is $TIndex$ the return type is $Size<TIndex>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type countOccurences(Iter< TIndex, VSTree<TSpec> > const &it) {
		if (_isSizeInval(value(it).i2))
			return length(indexSA(container(it))) - value(it).i1;
		else
			return value(it).i2 - value(it).i1;
	}

/**
.Function.getOccurences:
..summary:Returns all occurences of the @Function.representative@ substring in the index text.
..cat:Index
..signature:getOccurences(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:All positions where the @Function.representative@ of $iterator$ occurs in the text (see @Tag.ESA_Text@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_SA>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type getOccurences(Iter< TIndex, VSTree<TSpec> > const &it) {
		if (_isSizeInval(value(it).i2))
			return infix(indexSA(container(it)), value(it).i1, length(indexSA(container(it))));
		else
			return infix(indexSA(container(it)), value(it).i1, value(it).i2);
	}

/**
.Function.alignment:
..summary:Returns an alignment of the occurences of the @Function.representative@ substring in the index text.
..cat:Index
..signature:alignment(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:A local alignment corresponding to the seed of the $iterator$.
..remarks:The @Function.representative@ must uniquely occur in every sequence (e.g. in MUMs), 
otherwise the seed returned is one many.
*/
/*
	template < typename TString, typename TSSetSpec, typename TIndexSpec, class TSpec >
	inline typename Align<TString const, ArrayGaps>
	alignment(Iter< Index< StringSet<TString, TSSetSpec>, TIndexSpec >, VSTree<TSpec> > const &it) 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TIndexSpec > TIndex;
		typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typename Iterator<TOccs>::Type TIter;

		Align<TString const, ArrayGaps> align;
		TIndex const &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs), occEnd = end(occs);
		while (occ != occEnd) {
			typename Size<TIndex>::Type seqNo = getSeqNo(*occ, stringSetLimits(index));
			typename Size<TIndex>::Type seqOfs = getSeqOffset(*occ, stringSetLimits(index));
			setSource(row(align, seqNo), value(indexText(index), seqNo), seqOfs, seqOfs + repLen);
			++occ;
		}
		return align;
	}

	template < typename TString, typename TConcSpec, typename TIndexSpec, class TSpec >
	inline typename Align<TString const, ArrayGaps>
	alignment(Iter< Index< StringSet<TString, ConcatDirect<TConcSpec> >, TIndexSpec >, VSTree<TSpec> > const &it) 
	{
		typedef Index< StringSet<TString, ConcatDirect<TConcSpec> >, TIndexSpec > TIndex;
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs>::Type TIter;

		Align<TString const, ArrayGaps> align;
		TIndex const &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs), occEnd = end(occs);
		while (occ != occEnd) {
			typename Size<TIndex>::Type seqNo = getSeqNo(*occ, stringSetLimits(index));
			typename Size<TIndex>::Type globOfs = posGlobalize(*occ, stringSetLimits(index));
			setSource(row(align, seqNo), concat(indexText(index)), globOfs, globOfs + repLen);
			++occ;
		}
		return align;
	}
*/
	template < typename TString, typename TConcSpec, typename TIndexSpec, class TSpec >
	inline Align<TString, ArrayGaps>
	alignment(Iter< Index< StringSet<TString, ConcatDirect<TConcSpec> >, TIndexSpec >, VSTree<TSpec> > &it) 
	{
		typedef Index< StringSet<TString, ConcatDirect<TConcSpec> >, TIndexSpec > TIndex;
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs>::Type TIter;

		Align<TString, ArrayGaps> align;
		TIndex &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs), occEnd = end(occs);
		while (occ != occEnd) {
			typename Size<TIndex>::Type seqNo = getSeqNo(*occ, stringSetLimits((TIndex const&)index));
			typename Size<TIndex>::Type globOfs = posGlobalize(*occ, stringSetLimits((TIndex const&)index));
			setSource(row(align, seqNo), concat(indexText(index)), globOfs, globOfs + repLen);
			++occ;
		}
		return align;
	}

/**
.Function.getOccurencesBWT:
..summary:Returns the characters left beside all occurence of the @Function.representative@ substring in the index text.
..cat:Index
..signature:getOccurencesBWT(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:All positions where the @Function.representative@ of $iterator$ occurs in the text (see @Tag.ESA_Text@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_BWT>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, ESA_BWT>::Type const >::Type getOccurencesBWT(Iter< TIndex, VSTree<TSpec> > const &it) {
		if (_isSizeInval(value(it).i2))
			return infix(indexBWT(container(it)), value(it).i1, length(indexSA(container(it))));
		else
			return infix(indexBWT(container(it)), value(it).i1, value(it).i2);
	}

/**
.Function.representative:
..summary:Returns a substring representing the path from root to $iterator$ node.
..cat:Index
..signature:representative(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:An @Spec.InfixSegment@ of the raw text of an index (see @Tag.ESA_RawText@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_RawText>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, ESA_RawText>::Type const >::Type 
	representative(Iter< TIndex, VSTree<TSpec> > const &it) {
		typedef typename Size<TIndex>::Type TSize;
		TSize occ = posGlobalize(getOccurence(it), stringSetLimits(container(it)));
		TSize len = repLength(it);
		return infix(indexRawText(container(it)), occ, occ + len);
	}


	// deprecated (use goDown and goRight instead)
	// generic find method
	// tests functor F on every child interval and aborts if one call returns true (= found)
	template < typename TIndex, class TSpec, class TFunctor >
	inline bool _processChildren(Iter< TIndex, VSTree< TopDown<TSpec> > > const &it, TFunctor &F) {
		if (isLeaf(it)) return false;

		typedef	Pair<typename Size<TIndex>::Type>			TPair;
		typedef Iter< TIndex, VSTree< TopDown<TSpec> > >	iterator;

		TPair child(value(it).i1, _getUp(value(it).i2, container(it)));
		if (!(value(it).i1 < child.i2 && child.i2 < value(it).i2))
			child.i2 = _getDown(value(it).i1, container(it));

		if (F(iterator(it, child))) return true;
		child.i1 = child.i2;
		while (_isNextl(child.i2, container(it))) {
			child.i2 = _getNextl(child.i2, container(it));
			if (F(iterator(it, child))) return true;
			child.i1 = child.i2;
		}
		if (!isRoot(it)) {
			child.i2 = value(it).i2;
			return F(iterator(it, child));
		} else
			return false;
	}

/**
.Function.countChildren:
..summary:Count the number of children of a tree node.
..cat:Index
..signature:countChildren(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:The number of children of a tree node.
If $iterator$'s container type is $TIndex$, the return type is $Size<TIndex>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type countChildren(Iter< TIndex, VSTree<TSpec> > const &it) {
		if (isLeaf(it)) return 0;

		typedef typename Size<TIndex>::Type TSize;

		TSize i = _getUp(value(it).i2, container(it));
		if (!(value(it).i1 < i && i < value(it).i2))
			i = _getDown(value(it).i1, container(it));

		TSize result = (isRoot(it))? 1: 2;
		while (_isNextl(i, container(it))) {
			i = _getNextl(i, container(it));
			++result;
		}
		return result;
	}

	// get the interval of SA of the subtree under the edge beginning with character c
	template < typename TIndex, class TSpec, typename TValue, typename TPair >
	inline bool _getInterval(Iter< TIndex, VSTree<TSpec> > const &it, TValue c, TPair &found) {
		if (isLeaf(it)) return false;

		typedef typename Size<TIndex>::Type	TSize;
		typedef Iter< TIndex, VSTree<TSpec> >	iterator;

		TPair child(value(it).i1, _getUp(value(it).i2, container(it)));
		if (!(value(it).i1 < child.i2 && child.i2 < value(it).i2))
			child.i2 = _getDown(value(it).i1, container(it));

		TSize _lcp = lcpAt(child.i2 - 1, container(it));
		if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
			found = child;
			return true;
		}
		child.i1 = child.i2;
		while (_isNextl(child.i2, container(it))) {
			child.i2 = _getNextl(child.i2, container(it));
		if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
				found = child;
				return true;
			}
			child.i1 = child.i2;
		}
		if (!isRoot(it)) {
			if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
				found.i1 = child.i1;
				found.i2 = value(it).i2;
				return true;
			}
		}
		return false;
	}


///.Function.goNext.param.iterator.type:Spec.BottomUp Iterator
///.Function.goNext.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TIndex, typename TSpec >
	inline void goNext(Iter<TIndex, VSTree<TSpec> > &it) {
		goNext(it, typename DefaultDFSOrder< Iter<TIndex, VSTree<TSpec> > >::Type());
	}

	template < typename TIndex, typename TSpec, typename TSpecTag >
	inline void goNext(Iter<TIndex, VSTree<TSpec> > &it, Tag<TSpecTag> const tagOrder) {
		goNextImpl(it, tagOrder);
	}


/**
.Function.goDown:
..summary:Iterates down one edge or a path in a tree.
..cat:Index
..signature:bool goDown(iterator)
..signature:bool goDown(iterator, char)
..signature:bool goDown(iterator, text[, lcp])
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDown Iterator
..param.char:$iterator$ goes down the edge beginning with $char$.
..param.text:$iterator$ goes down the path representing $text$. If $text$ ends within an edge, $iterator$ will point to the child-end of this edge.
..param.lcp:A reference of a size type. When $goDown$ returns, $lcp$ contains the length of the longest-common-prefix of $text$ and a path beginning at the $iterator$ node.
...type:Class.String
...type:Class.Segment
..remarks:$goDown(iterator)$ goes down the leftmost edge in the Suffix Tree, i.e. the edge beginning with the lexicographically smallest character.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
*/

    //////////////////////////////////////////////////////////////////////////////
	// unified history stack access for goDown(..)
       
	template < typename TIndex, class TSpec, typename TStackEntry >
	inline void _historyPush(Iter< TIndex, VSTree<TSpec> > &it, TStackEntry _value) {
		it.topRight = value(it).i2;
	}

	template < typename TIndex, class TSpec, typename TStackEntry >
	inline void _historyPush(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, TStackEntry _value) {
		push(it.history, _value);
	}


	// go down the leftmost edge
	template < typename TIndex, class TSpec >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<TSpec> > > &it) {
		if (isLeaf(it)) return false;
		_historyPush(it, value(it));

		typename Size<TIndex>::Type i = _getUp(value(it).i2, container(it));
		if (value(it).i1 < i && i < value(it).i2)
			value(it).i2 = i;
		else
			value(it).i2 = _getDown(value(it).i1, container(it));
		return true;
	}

	// go down the edge beginning with c (returns false iff this edge doesn't exists)
	template < typename TIndex, class TSpec, typename TValue >
	inline bool _goDownChar(Iter< TIndex, VSTree< TopDown<TSpec> > > &it, TValue c) {
		Pair<typename Size<TIndex>::Type> prev = value(it);
		if (_getInterval(it, c, value(it))) {
			_historyPush(it, prev);
			return true;
		}
		return false;
	}

	// go down the path corresponding to pattern
	// lcp is the longest prefix of pattern and path
	template < typename TText, class TSpec, typename TString, typename TSize >
	inline bool _goDownString(Iter< Index<TText, Index_ESA<void> >, VSTree< TopDown<TSpec> > > &node, TString const &pattern, TSize &lcp) {
		typedef typename Infix< typename Fibre<Index<TText, Index_ESA<void> >, ESA_RawText>::Type const >::Type	TInfix;
		typedef typename Iterator<TInfix>::Type											IText;
		typedef typename Iterator<TString const>::Type									IPattern;
		
		IPattern p_begin = begin(pattern), p_end = end(pattern);
		IText t_begin, t_end;

		if (p_begin == p_end) {
			lcp = 0;
			return true;
		}

		TSize c = 0;
		while (goDown(node, *p_begin)) {
			TInfix t = representative(node);
			t_begin = begin(t) + c;
			t_end = end(t);

			while (t_begin != t_end && p_begin != p_end && *p_begin == *t_begin) {
				++t_begin;
				++p_begin;
			}
			if (p_begin == p_end) {
				lcp = length(pattern);
				return true;
			}
			c = length(t);
		}
		lcp = p_begin - begin(pattern);
		return false;
	}

	template < typename TIndex, typename TSpecIter, typename TValue >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, TValue const &c) {
		return _goDownChar(it, c);
	}

	template < typename TIndex, typename TSpecIter, typename TValue, typename TSpec >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, String<TValue, TSpec> const &pattern) {
		typename Size<TIndex>::Type dummy;
		return _goDownString(it, pattern, dummy);
	}

	template < typename TIndex, typename TSpecIter, typename THost, typename TSpec >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, Segment<THost, TSpec> const &pattern) {
		typename Size<TIndex>::Type dummy;
		return _goDownString(it, pattern, dummy);
	}


	template < typename TIndex, typename TSpecIter, typename TValue, typename TSpec, typename TSize >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, String<TValue, TSpec> const &pattern, TSize lcp) {
		return _goDownString(it, pattern, lcp);
	}

	template < typename TIndex, typename TSpecIter, typename THost, typename TSpec, typename TSize  >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, Segment<THost, TSpec> const &pattern, TSize lcp) {
		return _goDownString(it, pattern, lcp);
	}

		
/**
.Function.goUp:
..summary:Iterates up one edge to the parent in a tree.
..cat:Index
..signature:goUp(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDownHistory Iterator
*/

	// go up one edge (returns false if in root node)
	template < typename TIndex, class TSpec >
	inline bool goUp(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		if (!isRoot(it)) {
			value(it) = top(it.history);
			pop(it.history);
			return true;
		}
		return false;
	}

/**
.Function.goRight:
..summary:Iterates to the next sibling in a tree.
..cat:Index
..signature:goRight(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDown Iterator
*/

	// go right to the lexic. next sibling
	template < typename TIndex, class TSpec >
	inline bool goRight(Iter< TIndex, VSTree< TopDown<TSpec> > > &it) {
		if (value(it).i2 == length(container(it))) return false;		// right-most intervals have no right (would cause trouble with i2=\infty

		if (!isRoot(it)) {
			if (value(it).i2 < it.topRight) {
				if (_isNextl(value(it).i2, container(it))) {
					value(it).i1 = value(it).i2;
					value(it).i2 = _getNextl(value(it).i2, container(it));
					return true;
				}
				value(it).i1 = value(it).i2;
				value(it).i2 = it.topRight;
				return true;
			}
		}
		return false;
	}

	template < typename TIndex, class TSpec >
	inline bool goRight(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		if (value(it).i2 == length(container(it))) return false;		// right-most intervals have no right (would cause trouble with i2=\infty

		if (!isRoot(it)) {
			Pair<typename Size<TIndex>::Type> _top_ = top(it.history);
			if (value(it).i2 < _top_.i2) {
				if (_isNextl(value(it).i2, container(it))) {
					value(it).i1 = value(it).i2;
					value(it).i2 = _getNextl(value(it).i2, container(it));
					return true;
				}
				value(it).i1 = value(it).i2;
				value(it).i2 = _top_.i2;
				return true;
			}
		}
		return false;
	}

/**
.Function.parentEdgeLabel:
..summary:Returns a substring representing the edge from an $iterator$ node to its parent.
..cat:Index
..signature:parentEdgeLabel(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.TopDownHistory Iterator
..returns:An @Spec.InfixSegment@ of the raw text of an index (see @Tag.ESA_RawText@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_RawText>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, ESA_RawText>::Type const >::Type 
	parentEdgeLabel(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > const &it) {
		typedef typename Size<TIndex>::Type TSize;
		TSize occ = posGlobalize(getOccurence(it), stringSetLimits(container(it)));
		TSize last_len, len;

		if (isRoot(it)) {
			last_len = 0;
			len = 0;
		} else {
			last_len = repLength(container(it), top(it.history));
			len = repLength(it);
		}
		return infix(indexRawText(container(it)), occ + last_len, occ + len);
	}


///.Function.clear.param.iterator.type:Spec.BottomUp Iterator
///.Function.clear.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TIndex, class TSpec >
	inline void clear(Iter<TIndex, VSTree<TSpec> > &it) {
		value(it) = Pair<typename Size<TIndex>::Type>(0, 0);
    }

	template < typename TIndex, class TSpec >
	inline void _dfsClear(Iter<TIndex, VSTree<TSpec> > &it) {
		_dfsRange(it) = Pair<typename Size<TIndex>::Type>(0, 0);
    }


    //////////////////////////////////////////////////////////////////////////////
	// dfs traversal for ParentLink iterators

	template < typename TIndex, typename TSpec >
	inline void goNextImpl(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, Preorder const) {
		// preorder dfs
		if (!goDown(it) && !goRight(it))
			while (goUp(it) && !goRight(it));
		if (isRoot(it)) clear(it);
	}

	template < typename TIndex, typename TSpec >
	inline void goNextImpl(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, Postorder const) {
		// postorder dfs
		if (goRight(it))
			while (goDown(it));
		else
			if (!goUp(it)) clear(it);
	}


    //////////////////////////////////////////////////////////////////////////////
	// boolean functions

	template < typename TIndex, class TSpec >
	inline bool eof(Iter<TIndex, VSTree<TSpec> > &it) {
		return !value(it).i2;
	}

	template < typename TIndex, class TSpec >
	inline bool eof(Iter<TIndex, VSTree<TSpec> > const &it) {
		return !value(it).i2;
	}

	template < typename TIndex, class TSpec >
	inline bool empty(Iter<TIndex, VSTree<TSpec> > &it) {
		return !value(it).i2;
	}

	template < typename TIndex, class TSpec >
	inline bool empty(Iter<TIndex, VSTree<TSpec> > const &it) {
		return !value(it).i2;
	}

///.Function.atEnd.param.iterator.type:Spec.BottomUp Iterator
///.Function.atEnd.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TIndex, class TSpec >
	inline bool atEnd(Iter<TIndex, VSTree<TSpec> > &it) {
		return !value(it).i2;
	}

	template < typename TIndex, class TSpec >
	inline bool atEnd(Iter<TIndex, VSTree<TSpec> > const &it) {
		return !value(it).i2;
	}

/**
.Function.isRoot:
..summary:Test whether iterator points to the root node.
..cat:Index
..signature:bool isRoot(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to the root of the tree, otherwise $false$.
*/

	template < typename TIndex, class TSpec >
	inline bool isRoot(Iter<TIndex, VSTree< BottomUp<TSpec> > > const &it) {
		return empty(it.history);
	}

	template < typename TIndex, class TSpec >
	inline bool isRoot(Iter<TIndex, VSTree<TSpec> > const &it) {
		return _isRoot(value(it));
	}

	template < typename TSize >
	inline bool _isRoot(Pair<TSize> const &value) {
		return _isSizeInval(value.i2);
	}

/**
.Function.isRightTerminal:
..summary:Test whether iterator points to a suffix.
..cat:Index
..signature:bool isRightTerminal(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to the node representing a suffix, otherwise $false$.
..remarks:Every leaf is also a right terminal (see @Function.isLeaf@), but not vice versa.
*/

	template < typename TIndex, class TSpec >
	inline bool isRightTerminal(Iter<TIndex, VSTree<TSpec> > const &it) {
		// do we reach a leaf in a suffix tree with trailing '$'
		typename SAValue<TIndex>::Type pos = getOccurence(it);
		TIndex const &index = container(it);
		
		return (getSeqOffset(pos, stringSetLimits(index)) + repLength(it) 
			== sequenceLength(getSeqNo(pos, stringSetLimits(index)), container(it)));
	}

/**
.Function.isLeftMaximal:
..summary:Test whether the occurences of an iterator's @Function.representative@ mutually differ in the character left of the hits.
..cat:Index
..signature:bool isLeftMaximal(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if there are at least two different characters left of the occurences, otherwise $false$.
..see:@Function.getOccurences@
*/

	template < typename TIndex, class TSpec >
	inline bool isLeftMaximal(Iter<TIndex, VSTree<TSpec> > const &it)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type	TOccs;
		typedef typename Infix< typename Fibre<TIndex, ESA_BWT>::Type const >::Type	TOccsBWT;
		typedef typename Value< typename Fibre<TIndex, ESA_BWT>::Type const >::Type	TValue;

		typedef typename Iterator<TOccs>::Type TIter;
		typedef typename Iterator<TOccsBWT>::Type TIterBWT;
		
		TIndex const &index = container(it);
		typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

		TOccs occs = getOccurences(it);
		TOccsBWT bwts = getOccurencesBWT(it);

		TIter oc = begin(occs), ocEnd = end(occs);
		TIterBWT bw = begin(bwts);

		if (oc == ocEnd) return true;
		if (posAtFirstLocal(*oc, limits)) return true;

		TValue seen = *bw;
		++oc; 
		++bw;
		if (oc == ocEnd) return true;

		do {
			if (posAtFirstLocal(*oc, limits)) return true;
			if (seen != *bw) return true;
			++oc;
			++bw;
		} while (oc != ocEnd);

		return false;
	}

/**
.Function.isPartiallyLeftExtensible:
..summary:Test whether the characters left of the two occurences of @Function.representative@ are equal.
..cat:Index
..signature:bool isPartiallyLeftExtensible(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if there are at least two different characters left of the occurences, otherwise $false$.
..see:@Function.getOccurences@
*/

	template < typename TIndex, class TSpec, typename TSet >
	inline bool isPartiallyLeftExtensible(Iter<TIndex, VSTree<TSpec> > const &it, TSet &charSet)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type	TOccs;
		typedef typename Infix< typename Fibre<TIndex, ESA_BWT>::Type const >::Type	TOccsBWT;
		typedef typename Value< typename Fibre<TIndex, ESA_BWT>::Type const >::Type	TValue;

		typedef typename Iterator<TOccs>::Type TIter;
		typedef typename Iterator<TOccsBWT>::Type TIterBWT;
		
		TIndex const &index = container(it);
		typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

		clear(charSet);

		TOccs occs = getOccurences(it);
		TOccsBWT bwts = getOccurencesBWT(it);

		TIter oc = begin(occs), ocEnd = end(occs);
		TIterBWT bw = begin(bwts);

		while (oc != ocEnd) {
			if (!posAtFirstLocal(*oc, limits)) {
				TValue c = *bw;
				if (in(c, charSet)) return true;
				insert(c, charSet);
			}
			++oc;
			++bw;
		}

		return false;
	}

	template < typename TIndex, class TSpec >
	inline bool isPartiallyLeftExtensible(Iter<TIndex, VSTree<TSpec> > const &it)
	{
		typename Set<typename Value<TIndex>::Type>::Type set;
		return isPartiallyLeftExtensible(it, set);
	}

/**
.Function.isUnique:
..summary:Test whether the @Function.representative@ occurs only once in every sequence.
..cat:Index
..signature:bool isUnique(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if there are at least two different characters left of the occurences, otherwise $false$.
..see:@Function.getOccurences@
*/

	template < typename TIndex, class TSpec, typename TSet >
	inline bool isUnique(Iter<TIndex, VSTree<TSpec> > const &it, TSet &set)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs>::Type TIter;
		typedef typename Size<TIndex>::Type TSize;

		TIndex const &index = container(it);

		clear(set);

		TOccs occs = getOccurences(it);
		TIter oc = begin(occs), ocEnd = end(occs);

		while (oc != ocEnd) {
			TSize seqNo = getSeqNo(*oc, stringSetLimits(index));
			if (in(seqNo, set)) return false;
			insert(seqNo, set);
			++oc;
		}

		return true;
	}

	template < typename TIndex, class TSpec >
	inline bool isUnique(Iter<TIndex, VSTree<TSpec> > const &it) {
		VectorSet<typename Size<TIndex>::Type, Alloc<> > set(countSequences(container(it)));
		return isUnique(it, set);
	}

/**
.Function.childrenAreLeaves:
..summary:Test whether iterator points to a node with only leaf-children.
..cat:Index
..signature:bool childrenAreLeaves(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to an inner node of the tree, whose children are leaves. Otherwise it is $false$.
*/

	template < typename TIndex, class TSpec >
	inline bool childrenAreLeaves(Iter<TIndex, VSTree<TSpec> > const &it) {
		return countChildren(it) == countOccurences(it);
	}

/**
.Function.isLeaf:
..summary:Test whether iterator points to a leaf.
..cat:Index
..signature:bool isLeaf(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to a leaf of the tree, otherwise $false$.
*/

	template < typename TIndex, class TSpec >
	inline bool isLeaf(Iter<TIndex, VSTree<TSpec> > const &it) {
		return _isLeaf(value(it));
	}

	template < typename TSize >
	inline bool _isLeaf(Pair<TSize> const &value) {
		// do we reach a leaf?
		return value.i1 + 1 >= value.i2;
	}


	//////////////////////////////////////////////////////////////////////////////
	// (more or less) internal functions for accessing the childtab

	template < typename TSize, typename TIndex >
	inline bool _isNextl(TSize i, TIndex const &index) {
		if (i >= length(index)) return false;
		TSize j = childAt(i, index);
		return (j > i) && lcpAt(j - 1, index) == lcpAt(i - 1, index);
	}

	template < typename TSize, typename TIndex >
	inline bool _isUp(TSize i, TIndex const &index) {
		if (i >= length(index)) return false;
		TSize j = childAt(i, index);
		return (j <= i) && lcpAt(j - 1, index) > lcpAt(i - 1, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getNextl(TSize i, TIndex const &index) {
		return childAt(i, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getUp(TSize i, TIndex const &index) {
		if (!_isSizeInval(i))
			return childAt(i - 1, index);
		else
			return childAt(0, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getDown(TSize i, TIndex const &index) {
		return childAt(i, index);
	}


	//////////////////////////////////////////////////////////////////////////////
	// depth-first search 

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> & _dfsRange(Iter< TIndex, VSTree< BottomUp<TSpec> > > &it) {
		return it._range; 
	}

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> const & _dfsRange(Iter< TIndex, VSTree< BottomUp<TSpec> > > const &it) {
		return it._range; 
	}

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type & _dfsLCP(Iter< TIndex, VSTree< BottomUp<TSpec> > > &it) { 
		return it._lcp; 
	}

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type _dfsLCP(Iter< TIndex, VSTree< BottomUp<TSpec> > > const &it) { 
		return it._lcp; 
	}


}

#endif
