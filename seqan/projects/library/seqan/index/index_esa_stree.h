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
		typedef typename VertexDescriptor<TIndex>::Type Type;
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
	class Iter< TIndex, VSTree< TopDown<TSpec> > > 
	{
	public:

		typedef typename VertexDescriptor<TIndex>::Type	TVertexDesc;
		typedef Iter									iterator;

		TIndex const	&index;		// container of all necessary tables
		TVertexDesc		vDesc;		// current interval in suffix array and
									// right border of parent interval (needed in goRight)

		Iter(TIndex &_index):
			index(_index)
		{
			indexRequire(_index, ESA_SA());
			indexRequire(_index, ESA_LCP());
			indexRequire(_index, ESA_ChildTab());
			
			vDesc.i1.i1 = 0;				// start in root node with range (0,infty)
			_setSizeInval(vDesc.i1.i2);		// infty is equivalent to length(index) and faster to compare
		}

		Iter(TIndex &_index, TVertexDesc const &_vDesc):
			index(_index),
			vDesc(_vDesc)
		{
			indexRequire(_index, ESA_SA());
			indexRequire(_index, ESA_LCP());
			indexRequire(_index, ESA_ChildTab());
		}

		Iter(Iter const &_origin):
			index(container(_origin)),
			vDesc(value(_origin)) {}

		Iter(Iter const &_origin, typename Value<TVertexDesc, 1>::Type const &_childRange):
			index(container(_origin)),
			vDesc(_childRange, value(_origin).i1.i2) {}	// _origin will be parent of this node
/*
		template <typename _TSpec>
		Iter(Iter<TIndex, VSTree< TopDown< ParentLinks<_TSpec> > > > const &_origin):
			_index(container(_origin)),
			_range(value(_origin)) 
		{
			if (!empty(_origin.history))
				topRight = top(_origin.history).i2;
		}
*/
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

	template < typename TIndex, class TSpec >
	class Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > >:
		public Iter< TIndex, VSTree< TopDown<> > >
	{
	public:

		typedef Iter< TIndex, VSTree< TopDown<> > >	TBase;
		typedef	Pair<typename Size<TIndex>::Type>	TStackEntry, TRange;
		typedef String<TStackEntry, Block<> >		TStack;
		typedef Iter								iterator;

		TStack			history;	// contains all previously visited intervals (allows to go up)

		Iter(TIndex &_index):
			TBase(_index) {}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			history(_origin.history) {}

 		Iter(Iter const &_origin, TRange const &_childRange):
			TBase((TBase const &)_origin, _childRange),
			history(_origin.history)
		{
			push(history, value(_origin.i1));
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
	class Iter< TIndex, VSTree< BottomUp<TSpec> > > 
	{
	public:

		typedef typename VertexDescriptor<TIndex>::Type	TVertexDesc;
		typedef typename Size<TIndex>::Type				TSize;
		typedef	Pair<TSize>								TStackEntry;
		typedef String<TStackEntry, Block<> >			TStack;
		typedef Iter									iterator;

		TIndex	const	&index;			// container of all necessary tables
		TVertexDesc		vDesc;			// current interval in suffix array and
										// right border of parent interval (unused here)
		TSize			lValue;			// current l-value
		TStack			history;		// contains all left borders of current l-intervals (== left borders of history intervals)

		Iter(TIndex &_index):
			index(_index),
			vDesc(TStackEntry(0,0), 0),
			lValue(0)
		{
			indexRequire(_index, ESA_SA());
			indexRequire(_index, ESA_LCP());

			if (!empty(indexSA(_index))) 
			{
				_dfsOnPush(*this, TStackEntry(0,0));
				goNextImpl(*this, typename GetVSTreeIteratorTraits< Iter >::Type());
			}
		}

		Iter(TIndex &_index, MinimalCtor):
			index(_index),
			vDesc(TStackEntry(0,0), 0),
			lValue(0) {}

		Iter(Iter const &_origin):
			index(container(_origin)),
			vDesc(value(_origin)),
			lValue(_dfsLCP(_origin)),
			history(_origin.history) {}
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




	template < typename TIndex, typename TSpec >
	inline void _dumpHistoryStack(Iter<TIndex, VSTree<TSpec> > &it) {
		for(typename Size<TIndex>::Type i = 0; i < length(it.history); ++i)
			::std::cerr << it.history[i] << '\t';
		::std::cerr << value(it) << ::std::endl;
	}

	template < typename TIndex, typename TSpec >
	inline bool _dfsReversedOrder(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it) {
        return lcpAt(_dfsRange(it).i2 - 1, container(it)) > top(it.history).i2;
	}

	// standard push/pop handlers of lcp-dfs-traversal
	template < typename TIndex, typename TSpec, typename TSize >
	inline void _dfsOnPop(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, TSize const) {
        _dfsRange(it).i1 = top(it.history).i1;
		_dfsLCP(it) = top(it.history).i2;
		pop(it.history);
	}

	template < typename TIndex, typename TSpec, typename TElement >
	inline void _dfsOnPush(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, TElement const &e) {
		push(it.history, e);
	}

	template < typename TIndex, typename TSpec >
	inline void _dfsOnLeaf(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it) {
		_setSizeInval(_dfsLCP(it));
	}


	// postorder bottom up iterator (dfs)
	template < typename TIndex, typename TSpec, typename THideEmptyEdges >
	inline void goNextImpl(
		Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, 
		VSTreeIteratorTraits<_Postorder, THideEmptyEdges> const) 
	{
		TIndex const &index = container(it);
		do {
			// postorder dfs via lcp-table
			if (isRoot(it)) {
				_dfsClear(it);
				return;
			}

			if (_dfsRange(it).i2)
			{
				typedef typename Size<TIndex>::Type TSize;
				typedef typename Iter<TIndex, VSTree< BottomUp<TSpec> > >::TStackEntry TStackEntry;
				TStackEntry	_top_ = top(it.history);
				TSize		lcp_i = lcpAt(_dfsRange(it).i2 - 1, index);

				if (lcp_i < _top_.i2) {
					_dfsOnPop(it, lcp_i);
					return;
				}

				if (lcp_i > _top_.i2) {
					_top_.i1 = _dfsRange(it).i1;
					_top_.i2 = lcp_i;
					_dfsOnPush(it, _top_);
				}

	// innerer Knoten:
	// wenn kein Pop, aber Push -> begehe mind. 2. Teilbaum irgendeines Vorfahrs
	// wenn kein Pop, kein Push -> verlasse mind. zweites Kindblatt
	// wenn Pop, aber kein Push -> verlasse Wurzel des mind. 2.Teilbaums
	// wenn Pop und Push        -> verlasse ersten Teilbaum (sieht Vater zum ersten Mal und pusht jenen)

	// wenn nach Pop ein Pop folgen würde	-> Vater ist Top of Stack
	// wenn nach Pop ein Push folgen würde	-> Vater erst beim Push auf Stack (-> zwischenspeichern)
			}

			// last lcp entry (== 0) causes removal of toplevel interval
			if ((_dfsRange(it).i1 = _dfsRange(it).i2++) == length(index)) {
				_dfsOnPop(it, 0);
				_dfsRange(it).i2 = _dfsRange(it).i1;
			} else {
				// skip $ leafs (empty edges)
				if (THideEmptyEdges::VALUE &&
					suffixLength(saAt(_dfsRange(it).i1, index), index) == lcpAt(_dfsRange(it).i1, index))
					continue;

				_dfsOnLeaf(it);
	// Blatt:
	// wenn danach kein Pop, aber Push -> Vater wird erst noch gepusht
	// wenn danach Pop				   -> Vater ist Top des Stack
			}
			return;
		} while (true);
	}

	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type 
	repLength(Iter< TIndex, VSTree<BottomUp<TSpec> > > const &it) 
	{
		typename Size<TIndex>::Type lcp;
		if (!_isSizeInval(lcp = it.lValue))
			return lcp;
		else
			return suffixLength(getOccurence(it), container(it));
	}


	template < typename TIndex, typename TRange, typename TPos >
	inline typename Size<TIndex>::Type
	repLength(TIndex const &index, Pair<TRange, TPos> const &vDesc) 
	{
		if (_isLeaf(vDesc)) return suffixLength(saAt(vDesc.i1.i1, index), index);
		if (_isRoot(vDesc)) return 0;

		typename Size<TIndex>::Type lval = _getUp(vDesc.i1.i2, index);
		if (!(vDesc.i1.i1 < lval && lval < vDesc.i1.i2))
			lval = _getDown(vDesc.i1.i1, index);
		
		return lcpAt(lval - 1, index);
	}

	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type
	repLength(Iter< TIndex, VSTree<TopDown<TSpec> > > const &it) 
	{
		return repLength(container(it), value(it));
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
		push(a.history, value(a).i1);
		push(b.history, value(b).i1);

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
		goUp(_lca);

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

		typedef typename Size<TIndex>::Type				TSize;
		typedef typename VertexDescriptor<TIndex>::Type	TDesc;

		// push current intervals
		push(a.history, value(a).i1);
		push(b.history, value(b).i1);

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

		TSize _lcp = (i0 > 0)? repLength(container(a), TDesc(a.history[i0 - 1], 0)): 0;

		// pop current intervals
		pop(a.history);
		pop(b.history);

		return _lcp;
	}

///.Function.container.param.iterator.type:Spec.VSTree Iterator

	template < typename TIndex, class TSpec >
	inline TIndex const & container(Iter< TIndex, VSTree<TSpec> > const &it) { return it.index; }

	template < typename TIndex, class TSpec >
	inline TIndex & container(Iter< TIndex, VSTree<TSpec> > &it) { return const_cast<TIndex&>(it.index); }


///.Function.value.param.iterator.type:Spec.VSTree Iterator

	template < typename TIndex, class TSpec >
	inline typename VertexDescriptor<TIndex>::Type & 
	value(Iter< TIndex, VSTree<TSpec> > &it) { 
		return it.vDesc;
	}

	template < typename TIndex, class TSpec >
	inline typename VertexDescriptor<TIndex>::Type const & 
	value(Iter< TIndex, VSTree<TSpec> > const &it) { 
		return it.vDesc;
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
	inline typename SAValue<TIndex>::Type 
	getOccurence(Iter< TIndex, VSTree<TSpec> > const &it) {
		return saAt(value(it).i1.i1, container(it));
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
	inline typename Size<TIndex>::Type 
	countOccurences(Iter< TIndex, VSTree<TSpec> > const &it) {
		if (_isSizeInval(value(it).i1.i2))
			return length(indexSA(container(it))) - value(it).i1.i1;
		else
			return value(it).i1.i2 - value(it).i1.i1;
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
	inline typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type 
	getOccurences(Iter< TIndex, VSTree<TSpec> > const &it) {
		if (_isSizeInval(value(it).i1.i2))
			return infix(indexSA(container(it)), value(it).i1.i1, length(indexSA(container(it))));
		else
			return infix(indexSA(container(it)), value(it).i1.i1, value(it).i1.i2);
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

	template < typename TString, typename TSSetSpec, typename TIndexSpec, class TSpec >
	inline Align<TString, ArrayGaps>
	alignment(Iter< Index< StringSet<TString, TSSetSpec>, TIndexSpec >, VSTree<TSpec> > &it) 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TIndexSpec > TIndex;
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;

		Align<TString, ArrayGaps> align;
		TIndex &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs, Standard()), occEnd = end(occs, Standard());
		while (occ != occEnd) {
			typename Size<TIndex>::Type seqNo = getSeqNo(*occ, stringSetLimits((TIndex const&)index));
			typename Size<TIndex>::Type seqOfs = getSeqOffset(*occ, stringSetLimits((TIndex const&)index));
			setSource(row(align, seqNo), value(indexText(index), seqNo), seqOfs, seqOfs + repLen);
			++occ;
		}
		return align;
	}
/*
	template < typename TString, typename TConcSpec, typename TIndexSpec, class TSpec >
	inline typename Align<TString const, ArrayGaps>
	alignment(Iter< Index< StringSet<TString, Owner<ConcatDirect<TConcSpec> > >, TIndexSpec >, VSTree<TSpec> > const &it) 
	{
		typedef Index< StringSet<TString, Owner<ConcatDirect<TConcSpec> > >, TIndexSpec > TIndex;
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;

		Align<TString const, ArrayGaps> align;
		TIndex const &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs, Standard()), occEnd = end(occs, Standard());
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
	alignment(Iter< Index< StringSet<TString, Owner<ConcatDirect<TConcSpec> > >, TIndexSpec >, VSTree<TSpec> > &it) 
	{
		typedef Index< StringSet<TString, Owner<ConcatDirect<TConcSpec> > >, TIndexSpec > TIndex;
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;

		Align<TString, ArrayGaps> align;
		TIndex &index = container(it);
		resize(rows(align), length(indexText(index)));	// resize alignment to number of sequences
		TOccs occs = getOccurences(it);
		typename Size<TIndex>::Type repLen = repLength(it);
		TIter occ = begin(occs, Standard()), occEnd = end(occs, Standard());
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
	inline typename Infix< typename Fibre<TIndex, ESA_BWT>::Type const >::Type 
	getOccurencesBWT(Iter< TIndex, VSTree<TSpec> > const &it) {
		if (_isSizeInval(value(it).i1.i2))
			return infix(indexBWT(container(it)), value(it).i1.i1, length(indexSA(container(it))));
		else
			return infix(indexBWT(container(it)), value(it).i1.i1, value(it).i1.i2);
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
	inline bool 
	_processChildren(Iter< TIndex, VSTree< TopDown<TSpec> > > const &it, TFunctor &F)
	{
		if (isLeaf(it)) return false;

		typedef	Pair<typename Size<TIndex>::Type>			TPair;
		typedef Iter< TIndex, VSTree< TopDown<TSpec> > >	iterator;

		TPair child(value(it).i1.i1, _getUp(value(it).i1.i2, container(it)));
		if (!(value(it).i1.i1 < child.i2 && child.i2 < value(it).i1.i2))
			child.i2 = _getDown(value(it).i1.i1, container(it));

		if (F(iterator(it, child))) return true;
		child.i1 = child.i2;
		while (_isNextl(child.i2, container(it))) {
			child.i2 = _getNextl(child.i2, container(it));
			if (F(iterator(it, child))) return true;
			child.i1 = child.i2;
		}
		if (!isRoot(it)) {
			child.i2 = value(it).i1.i2;
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
	inline typename Size<TIndex>::Type 
	countChildren(Iter< TIndex, VSTree<TSpec> > const &it) {
		if (isLeaf(it)) return 0;

		typedef typename Size<TIndex>::Type TSize;

		TSize i = _getUp(value(it).i1.i2, container(it));
		if (!(value(it).i1.i1 < i && i < value(it).i1.i2))
			i = _getDown(value(it).i1.i1, container(it));

		TSize result = (isRoot(it))? 1: 2;
		while (_isNextl(i, container(it))) {
			i = _getNextl(i, container(it));
			++result;
		}
		return result;
	}

	// get the interval of SA of the subtree under the edge beginning with character c
	template < typename TIndex, class TSpec, typename TValue >
	inline bool 
	_getNodeByChar(
		Iter< TIndex, VSTree<TSpec> > const &it, 
		TValue c, 
		typename VertexDescriptor<TIndex>::Type &childDesc)
	{
		if (isLeaf(it)) return false;

		typedef typename Size<TIndex>::Type		TSize;
		typedef Iter< TIndex, VSTree<TSpec> >	iterator;

		Pair<TSize> child(value(it).i1.i1, _getUp(value(it).i1.i2, container(it)));
		if (!(value(it).i1.i1 < child.i2 && child.i2 < value(it).i1.i2))
			child.i2 = _getDown(value(it).i1.i1, container(it));

		TSize _lcp = lcpAt(child.i2 - 1, container(it));
		if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
			childDesc.i1 = child;
			childDesc.i2 = value(it).i1.i2;
			return true;
		}
		child.i1 = child.i2;
		while (_isNextl(child.i2, container(it))) {
			child.i2 = _getNextl(child.i2, container(it));
		if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
				childDesc.i1 = child;
				childDesc.i2 = value(it).i1.i2;
				return true;
			}
			child.i1 = child.i2;
		}
		if (!isRoot(it)) {
			if (textAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
				childDesc.i1.i1 = child.i1;
				childDesc.i1.i2 = childDesc.i2 = value(it).i1.i2;
				return true;
			}
		}
		return false;
	}


///.Function.goNext.param.iterator.type:Spec.BottomUp Iterator
///.Function.goNext.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TIndex, typename TSpec >
	inline void goNext(Iter<TIndex, VSTree<TSpec> > &it) {
		goNext(it, typename GetVSTreeIteratorTraits< Iter<TIndex, VSTree<TSpec> > >::Type());
	}

	template < typename TIndex, typename TSpec, typename TTraits >
	inline void goNext(Iter<TIndex, VSTree<TSpec> > &it, TTraits const traits) {
		goNextImpl(it, traits);
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
	inline void _historyPush(Iter< TIndex, VSTree<TSpec> > &it, TStackEntry) {
		value(it).i2 = value(it).i1.i2;
	}

	template < typename TIndex, class TSpec, typename TStackEntry >
	inline void _historyPush(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, TStackEntry range) {
		value(it).i2 = value(it).i1.i2;
		push(it.history, range);
	}


    //////////////////////////////////////////////////////////////////////////////
	// goDown

	// go down the leftmost edge (including empty $-edges)
	template < typename TIndex, class TSpec, typename TDFSOrder >
	inline bool _goDown(
		Iter< TIndex, VSTree< TopDown<TSpec> > > &it,
		VSTreeIteratorTraits<TDFSOrder, False> const)
	{
		if (isLeaf(it)) return false;
		_historyPush(it, value(it).i1);

		TIndex const &index = container(it);

		typename Size<TIndex>::Type lval = _getUp(value(it).i1.i2, index);
		if (!(value(it).i1.i1 < lval && lval < value(it).i1.i2))
			lval = _getDown(value(it).i1.i1, index);
		value(it).i1.i2 = lval;
		return true;
	}

	// go down the leftmost edge (skip empty $-edges) 
	// without calling goUp(..)
	template < typename TIndex, class TSpec, typename TDFSOrder >
	inline bool _goDown(
		Iter< TIndex, VSTree< TopDown<TSpec> > > &it,
		VSTreeIteratorTraits<TDFSOrder, True> const)
	{
		typedef Iter<TIndex, VSTree< TopDown<TSpec> > >	TIter;
		typedef typename Value<TIter>::Type				TDesc;

		if (isLeaf(it)) return false;

		TDesc desc = value(it);			// save descriptor of the current node
		_historyPush(it, desc.i1);

		TIndex const &index = container(it);

		typename Size<TIndex>::Type lval = _getUp(value(it).i1.i2, index);
		if (!(value(it).i1.i1 < lval && lval < value(it).i1.i2))
			lval = _getDown(value(it).i1.i1, index);
		value(it).i1.i2 = lval;

		typename Size<TIndex>::Type lcp = lcpAt(lval - 1, index);
		//typename typename StringSetLimits<TIndex const>::Type &limits = stringSetLimits(index);
		while (isLeaf(it)) {
			typename SAValue<TIndex>::Type pos = getOccurence(it);
			if (getSeqOffset(pos, stringSetLimits(index)) + lcp == sequenceLength(getSeqNo(pos, stringSetLimits(index)), index)) {
				if (!goRight(it)) {
					value(it) = desc;	// restore descriptor
					return false;
				}
			} else
				break;
		}
		return true;
	}

	// go down the leftmost edge (skip empty $-edges)
	template < typename TIndex, class TSpec, typename TDFSOrder >
	inline bool _goDown(
		Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it,
		VSTreeIteratorTraits<TDFSOrder, True> const)
	{
		if (isLeaf(it)) return false;
		_historyPush(it, value(it).i1);

		TIndex const &index = container(it);

		typename Size<TIndex>::Type lval = _getUp(value(it).i1.i2, index);
		if (!(value(it).i1.i1 < lval && lval < value(it).i1.i2))
			lval = _getDown(value(it).i1.i1, index);
		value(it).i1.i2 = lval;

		typename Size<TIndex>::Type lcp = lcpAt(lval - 1, index);
		//typename typename StringSetLimits<TIndex const>::Type &limits = stringSetLimits(index);
		while (isLeaf(it)) {
			typename SAValue<TIndex>::Type pos = getOccurence(it);
			if (getSeqOffset(pos, stringSetLimits(index)) + lcp == sequenceLength(getSeqNo(pos, stringSetLimits(index)), index)) {
				if (!goRight(it)) {
					goUp(it);
					return false;
				}
			} else
				break;
		}
		return true;
	}

	// go down the leftmost edge
	template < typename TIndex, class TSpec >
	inline bool goDown(Iter< TIndex, VSTree< TopDown<TSpec> > > &it) {
		return _goDown(it, typename GetVSTreeIteratorTraits< Iter<TIndex, VSTree< TopDown<TSpec> > > >::Type());
	}


    //////////////////////////////////////////////////////////////////////////////
	// goDown a specific edge (chosen by the first character)

	// go down the edge beginning with c (returns false iff this edge doesn't exists)
	template < typename TIndex, class TSpec, typename TValue >
	inline bool _goDownChar(Iter< TIndex, VSTree< TopDown<TSpec> > > &it, TValue c) 
	{
		Pair<typename Size<TIndex>::Type> oldRange = value(it).i1;
		if (_getNodeByChar(it, c, value(it))) {
			_historyPush(it, oldRange);
			return true;
		}
		return false;
	}

	// go down the path corresponding to pattern
	// lcp is the longest prefix of pattern and path
	template < typename TText, class TSpec, typename TString, typename TSize >
	inline bool
	_goDownString(
		Iter< Index<TText, Index_ESA<void> >, 
		VSTree< TopDown<TSpec> > > &node, 
		TString const &pattern, 
		TSize &lcp) 
	{
		typedef typename Infix< 
					typename Fibre<Index<TText, Index_ESA<void> >, ESA_RawText>::Type const 
				>::Type	TInfix;
		typedef typename Iterator<TInfix, Standard>::Type			IText;
		typedef typename Iterator<TString const, Standard>::Type	IPattern;
		
		IPattern p_begin = begin(pattern, Standard()), p_end = end(pattern, Standard());
		IText t_begin, t_end;

		if (p_begin == p_end) {
			lcp = 0;
			return true;
		}

		TSize c = 0;
		while (goDown(node, *p_begin)) {
			TInfix t = representative(node);
			t_begin = begin(t, Standard()) + c;
			t_end = end(t, Standard());

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
		lcp = p_begin - begin(pattern, Standard());
		return false;
	}

	template < typename TIndex, typename TSpecIter, typename TValue >
	inline bool 
	goDown(
		Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, 
		TValue const &c) 
	{
		return _goDownChar(it, c);
	}

	template < typename TIndex, typename TSpecIter, typename TValue, typename TSpec >
	inline bool 
	goDown(
		Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, 
		String<TValue, TSpec> const &pattern) 
	{
		typename Size<TIndex>::Type dummy;
		return _goDownString(it, pattern, dummy);
	}

	template < typename TIndex, typename TSpecIter, typename THost, typename TSpec >
	inline bool 
	goDown(
		Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, 
		Segment<THost, TSpec> const &pattern) 
	{
		typename Size<TIndex>::Type dummy;
		return _goDownString(it, pattern, dummy);
	}


	template < typename TIndex, typename TSpecIter, typename TValue, typename TSpec, typename TSize >
	inline bool 
	goDown(
		Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, 
		String<TValue, TSpec> const &pattern, 
		TSize lcp) 
	{
		return _goDownString(it, pattern, lcp);
	}

	template < typename TIndex, typename TSpecIter, typename THost, typename TSpec, typename TSize  >
	inline bool 
	goDown(
		Iter< TIndex, VSTree< TopDown<TSpecIter> > > &it, 
		Segment<THost, TSpec> const &pattern, 
		TSize lcp) 
	{
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
	inline bool 
	goUp(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) 
	{
		if (!empty(it.history)) {
			value(it).i1 = top(it.history);
			pop(it.history);
			if (!empty(it.history))
				value(it).i2 = top(it.history).i2;	// copy right boundary of parent's range
			return true;
		}
		return false;
	}

	// return vertex descriptor of parent's node
	template < typename TIndex, class TSpec >
	inline typename VertexDescriptor<TIndex>::Type
	nodeUp(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > const &it) 
	{
		if (!empty(it.history)) {
			typename Size<TIndex>::Type parentRight = 0;
			if (length(it.history) > 2)
				parentRight = topPrev(it.history).i2;
			return typename VertexDescriptor<TIndex>::Type(top(it.history), parentRight);
		} else
			return value(it);
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
	inline bool goRight(Iter< TIndex, VSTree< TopDown<TSpec> > > &it) 
	{
		if (value(it).i1.i2 == length(container(it)))
			return false;		// right-most intervals have no right (would cause trouble with i2=\infty)

		if (!isRoot(it)) {
			if (value(it).i1.i2 < value(it).i2)			// not the right-most child?
			{
				if (_isNextl(value(it).i1.i2, container(it))) 
				{
					value(it).i1.i1 = value(it).i1.i2;	// go right
					value(it).i1.i2 = _getNextl(value(it).i1.i2, container(it));
					return true;
				}
				value(it).i1.i1 = value(it).i1.i2;		// now it is the right-most child
				value(it).i1.i2 = value(it).i2;
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
	parentEdgeLabel(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > const &it) 
	{
		typedef typename Size<TIndex>::Type TSize;

		if (isRoot(it))
			return infix(indexRawText(container(it)), 0, 0);
		else {
			TSize occ = posGlobalize(getOccurence(it), stringSetLimits(container(it)));
			return infix(
				indexRawText(container(it)), 
				occ + repLength(container(it), nodeUp(it)), 
				occ + repLength(it));
		}
	}


///.Function.clear.param.iterator.type:Spec.BottomUp Iterator
///.Function.clear.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TIndex, class TSpec >
	inline void clear(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		value(it).i1 = Pair<typename Size<TIndex>::Type>(0, 0);
		value(it).i2 = 0;
    }

	template < typename TIndex, class TSpec >
	inline void _dfsClear(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		_dfsRange(it) = Pair<typename Size<TIndex>::Type>(0, 0);
    }


    //////////////////////////////////////////////////////////////////////////////
	// dfs traversal for ParentLink iterators

	template < typename TIndex, typename TSpec, typename THideEmptyEdges >
	inline void goNextImpl(
		Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, 
		VSTreeIteratorTraits<_Preorder, THideEmptyEdges> const)
	{
		// preorder dfs
		if (!goDown(it) && !goRight(it))
			while (goUp(it) && !goRight(it));
		if (isRoot(it)) clear(it);
	}

	template < typename TIndex, typename TSpec, typename THideEmptyEdges >
	inline void goNextImpl(
		Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, 
		VSTreeIteratorTraits<_Postorder, THideEmptyEdges> const)
	{
		// postorder dfs
		if (goRight(it))
			while (goDown(it));
		else
			if (!goUp(it)) clear(it);
	}


    //////////////////////////////////////////////////////////////////////////////
	// boolean functions

	template < typename TIndex, class TSpec >
	inline bool eof(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		return !value(it).i1.i2;
	}

	template < typename TIndex, class TSpec >
	inline bool eof(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return !value(it).i1.i2;
	}

	template < typename TIndex, class TSpec >
	inline bool empty(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		return !value(it).i1.i2;
	}

	template < typename TIndex, class TSpec >
	inline bool empty(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return !value(it).i1.i2;
	}

///.Function.atEnd.param.iterator.type:Spec.BottomUp Iterator
///.Function.atEnd.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TIndex, class TSpec >
	inline bool atEnd(Iter<TIndex, VSTree<TSpec> > &it) 
	{
		return !value(it).i1.i2;
	}

	template < typename TIndex, class TSpec >
	inline bool atEnd(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return !value(it).i1.i2;
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
	inline bool isRoot(Iter<TIndex, VSTree< BottomUp<TSpec> > > const &it) 
	{
		return empty(it.history);
	}

	template < typename TIndex, class TSpec >
	inline bool isRoot(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return _isRoot(value(it));
	}

	template < typename TRange, typename TPos >
	inline bool _isRoot(Pair<TRange, TPos> const &value) 
	{
		return _isSizeInval(value.i1.i2);
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
	inline bool isRightTerminal(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		// do we reach a leaf in a suffix tree with trailing '$'
		typename SAValue<TIndex>::Type pos = getOccurence(it);
		TIndex const &index = container(it);
		typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

		return (getSeqOffset(pos, limits) + repLength(it) 
			== sequenceLength(getSeqNo(pos, limits), index));
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

		typedef typename Iterator<TOccs, Standard>::Type	TIter;
		typedef typename Iterator<TOccsBWT, Standard>::Type TIterBWT;
		
		TIndex const &index = container(it);
		typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

		TOccs occs = getOccurences(it);
		TOccsBWT bwts = getOccurencesBWT(it);

		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());
		TIterBWT bw = begin(bwts, Standard());

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

		typedef typename Iterator<TOccs, Standard>::Type	TIter;
		typedef typename Iterator<TOccsBWT, Standard>::Type TIterBWT;
		
		TIndex const &index = container(it);
		typename StringSetLimits<typename Host<TIndex>::Type const>::Type &limits = stringSetLimits(index);

		clear(charSet);

		TOccs occs = getOccurences(it);
		TOccsBWT bwts = getOccurencesBWT(it);

		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());
		TIterBWT bw = begin(bwts, Standard());

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
		typedef typename Iterator<TOccs, Standard>::Type TIter;
		typedef typename Size<TIndex>::Type TSize;

		TIndex const &index = container(it);

		clear(set);

		TOccs occs = getOccurences(it);
		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());

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
		VectorSet<
			typename Size<TIndex>::Type,
			Alloc<> 
		> set(countSequences(container(it)));
		return isUnique(it, set);
	}

/**
.Function.getFrequency:
..summary:Returns the number of sequences, which contain the @Function.representative@ as a substring.
..cat:Index
..signature:int getFrequency(iterator)
..param.iterator:An iterator of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:The number of different sequences containing the @Function.representative@.
..see:@Function.getOccurences@
*/

	template < typename TIndex, class TSpec, typename TSet >
	inline int getFrequency(Iter<TIndex, VSTree<TSpec> > const &it, TSet &set)
	{
		typedef typename Infix< typename Fibre<TIndex, ESA_SA>::Type const >::Type TOccs;
		typedef typename Iterator<TOccs, Standard>::Type TIter;
		typedef typename Size<TIndex>::Type TSize;

		TIndex const &index = container(it);

		clear(set);

		TOccs occs = getOccurences(it);
		TIter oc = begin(occs, Standard()), ocEnd = end(occs, Standard());

		int counter = 0;
		while (oc != ocEnd) {
			TSize seqNo = getSeqNo(*oc, stringSetLimits(index));
			if (!in(seqNo, set)) {
				++counter;
				insert(seqNo, set);
			}
			++oc;
		}

		return counter;
	}

	template < typename TIndex, class TSpec >
	inline int getFrequency(Iter<TIndex, VSTree<TSpec> > const &it) {
		VectorSet<
			typename Size<TIndex>::Type,
			Alloc<> 
		> set(countSequences(container(it)));
		return getFrequency(it, set);
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
	inline bool childrenAreLeaves(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
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
	inline bool isLeaf(Iter<TIndex, VSTree<TSpec> > const &it) 
	{
		return _isLeaf(value(it));
	}

	template < typename TRange, typename TPos >
	inline bool _isLeaf(Pair<TRange, TPos> const &vDesc) 
	{
		// do we reach a leaf?
		return vDesc.i1.i1 + 1 >= vDesc.i1.i2;
	}


	//////////////////////////////////////////////////////////////////////////////
	// (more or less) internal functions for accessing the childtab

	template < typename TSize, typename TIndex >
	inline bool _isNextl(TSize i, TIndex const &index) 
	{
		if (i >= length(index)) return false;
		TSize j = childAt(i, index);
		return (j > i) && lcpAt(j - 1, index) == lcpAt(i - 1, index);
	}

	template < typename TSize, typename TIndex >
	inline bool _isUp(TSize i, TIndex const &index) 
	{
		if (i >= length(index)) return false;
		TSize j = childAt(i, index);
		return (j <= i) && lcpAt(j - 1, index) > lcpAt(i - 1, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getNextl(TSize i, TIndex const &index) 
	{
		return childAt(i, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getUp(TSize i, TIndex const &index) 
	{
		if (!_isSizeInval(i))
			return childAt(i - 1, index);
		else
			return childAt(0, index);
	}

	template < typename TSize, typename TIndex >
	inline TSize _getDown(TSize i, TIndex const &index) 
	{
		return childAt(i, index);
	}


	//////////////////////////////////////////////////////////////////////////////
	// depth-first search 

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> &
	_dfsRange(Iter< TIndex, VSTree< BottomUp<TSpec> > > &it)
	{
		return value(it).i1;
	}

	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> const & 
	_dfsRange(Iter< TIndex, VSTree< BottomUp<TSpec> > > const &it) 
	{
		return value(it).i1;
	}

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type & _dfsLCP(Iter< TIndex, VSTree< BottomUp<TSpec> > > &it)
	{
		return it.lValue;
	}

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type _dfsLCP(Iter< TIndex, VSTree< BottomUp<TSpec> > > const &it)
	{
		return it.lValue;
	}


}

#endif
