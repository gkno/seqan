/*
 *  index_esa_stree.h
 *  genindex
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
		typedef ::std::vector<TStackEntry>			TStack;
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
			_historyPush(*this, value(_origin));
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
		typedef ::std::stack<TStackEntry>		TStack;
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

			_lcp = length(_index) - saAt(0, _index);
			history.push(TStackEntry(0,0));
		}

		Iter(Iter const &_origin):
			_index(container(_origin)),
			_range(_dfsRange(_origin)),
			_lcp(_dfsLCP(_origin)),
			history(_origin.history) {}
	};


	// standard push/pop handlers of lcp-dfs-traversal
	template < typename TIndex, typename TSpec >
	inline void _postorderPop(Iter<TIndex, TSpec> &it) {
        _dfsRange(it).i1 = _atomicTop(it.history).i1;
		_dfsLCP(it) = _atomicTop(it.history).i2;
		_atomicPop(it.history);
	}

	template < typename TIndex, typename TSpec, typename TElement >
	inline void _postorderPush(Iter<TIndex, TSpec> &it, TElement const &e) {
		_atomicPush(it.history, e);
	}

	template < typename TIndex, typename TSpec >
	inline void _postorderLeaf(Iter<TIndex, TSpec> &it) {
		_setSizeInval(_dfsLCP(it));
		//length(container(it)) - saAt(_dfsRange(it).i1, container(it));
	}


	template < typename TIndex, typename TSpec >
	inline void goNextImpl(Iter<TIndex, VSTree< BottomUp<TSpec> > > &it, Postorder const) {
		// postorder dfs via lcp-table
		if (isRoot(it)) {
			_dfsClear(it);
			return;
		}

		typedef typename Size<TIndex>::Type TSize;
		typedef typename Iter<TIndex, VSTree< BottomUp<TSpec> > >::TStackEntry TStackEntry;
		TStackEntry	top = _atomicTop(it.history);
		TSize		lcp_i = lcpAt(_dfsRange(it).i2 - 1, container(it));

        if (lcp_i < top.i2) {
			_postorderPop(it);			// go up
			return;
		}

		if (lcp_i > top.i2) {
			top.i1 = _dfsRange(it).i1;
            top.i2 = lcp_i;
			_postorderPush(it, top);		// go down
        }

		// last lcp entry (== 0) causes removal of toplevel interval
		if ((_dfsRange(it).i1 = _dfsRange(it).i2++) == length(container(it))) {
			_postorderPop(it);
			_dfsRange(it).i2 = _dfsRange(it).i1;
		} else
			_postorderLeaf(it);
	}
/*
	// asumes the history container to be a vector
	template < typename TIndex, typename TSpec, typename TEntry, typename THandler >
	inline void goNextImpl(Iter<TIndex, TSpec> &it, ::std::vector<TEntry> &history, THandler &handler) {
		// postorder dfs via lcp-table
		if (isRoot(it)) {
			_dfsClear(it);
			return;
		}

		typedef typename Size<TIndex>::Type TSize;
		typedef typename Iter<TIndex, TSpec>::TStackEntry TStackEntry;
		TStackEntry	top = history.back();
		TSize		lcp_i = lcpAt(_dfsRange(it).i2 - 1, container(it));

        if (lcp_i < top.i2) {
			handler.beforeUp(it);
            _dfsRange(it).i1 = history.back().i1;
			_dfsLCP(it) = history.back().i2;
			history.pop_back();		// go up
			handler.afterUp(it);
			return;
		}

		if (lcp_i > top.i2) {
			handler.beforeDown(it);
			top.i1 = _dfsRange(it).i1;
            top.i2 = lcp_i;
			history.push_back(top);	// go down
			handler.afterDown(it);
        }

		// last lcp entry (== 0) causes removal of toplevel interval
		if ((_dfsRange(it).i1 = _dfsRange(it).i2++) == length(container(it))) {
			handler.beforeUp(it);
			_dfsRange(it).i2 = _dfsRange(it).i1;
            _dfsRange(it).i1 = 0;
            history.pop_back();
			_dfsLCP(it) = history.back().i2;
			handler.afterUp(it);
		} else
			_setSizeInval(_dfsLCP(it))
			//length(container(it)) - saAt(_dfsRange(it).i1, container(it));
	}
*/

	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type _repLength(Iter< TIndex, VSTree<BottomUp<TSpec> > > const &it) {
		typename Size<TIndex>::Type lcp;
		if (!_isSizeInval(lcp = _dfsLCP(it)))
			return lcp;
		else
			return length(container(it)) - saAt(_dfsRange(it).i1, container(it));
	}


	template < typename TIndex, typename TSpec >
	inline typename Size<TIndex>::Type _repLength(Iter< TIndex, VSTree<TopDown<TSpec> > > const &it) {
		return _repLength(container(it), value(it));
	}

	template < typename TIndex, class TSize >
	inline TSize _repLength(TIndex const &index, Pair<TSize> const &value) {
		if (_isLeaf(value)) return length(index) - saAt(value.i1, index);
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
		_historyPush(a, value(a));
		_historyPush(b, value(b));

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
		_historyPop(a);
		_historyPop(b);
		if (i0) _historyPop(_lca, value(_lca));

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
		_historyPush(a, value(a));
		_historyPush(b, value(b));

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

		TSize _lcp = (i0 > 0)? _repLength(container(a), a.history[i0 - 1]): 0;

		// pop current intervals
		_historyPop(a);
		_historyPop(b);

		return _lcp;
	}

///.Function.container.param.iterator.type:Spec.VSTree Iterator

	template < typename TIndex, class TSpec >
	inline TIndex const & container(Iter< TIndex, VSTree<TSpec> > const &it) { return it._index; }


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
..param.iterator:An iterater of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:A position where the @Function.representative@ of $iterator$ occurs in the raw text (see @Tag.ESA_RawText@).
*/

	template < typename TIndex, class TSpec >
	inline typename Size<TIndex>::Type getOccurence(Iter< TIndex, VSTree<TSpec> > const &it) {
		return saAt(value(it).i1, container(it));
	}

/**
.Function.getOccurences:
..summary:Returns all occurence of the @Function.representative@ substring in the index text.
..cat:Index
..signature:getOccurences(iterator)
..param.iterator:An iterater of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:All positions where the @Function.representative@ of $iterator$ occurs in the raw text (see @Tag.ESA_RawText@).
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
.Function.representative:
..summary:Returns a substring representing the path from root to $iterator$ node.
..cat:Index
..signature:representative(iterator)
..param.iterator:An iterater of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:An @Spec.InfixSegment@ of the raw text of an index (see @Tag.ESA_RawText@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, ESA_RawText>::Type>::Type$.
*/

	template < typename TIndex, class TSpec >
	inline typename Infix< typename Fibre<TIndex, ESA_RawText>::Type const >::Type 
	representative(Iter< TIndex, VSTree<TSpec> >  &it) {
		typename Size<TIndex>::Type occ = getOccurence(it), len = _repLength(it);
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
..param.iterator:An iterater of a Suffix Tree.
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
		if (rawtextAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
			found = child;
			return true;
		}
		child.i1 = child.i2;
		while (_isNextl(child.i2, container(it))) {
			child.i2 = _getNextl(child.i2, container(it));
		if (rawtextAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
				found = child;
				return true;
			}
			child.i1 = child.i2;
		}
		if (!isRoot(it)) {
			if (rawtextAt(saAt(child.i1, container(it)) + _lcp, container(it)) == c) {
				found.i1 = child.i1;
				found.i2 = value(it).i2;
				return true;
			}
		}
		return false;
	}


    //////////////////////////////////////////////////////////////////////////////
	// unified access to stacks and vectors

	template < typename TEntry, typename TStackEntry >
	inline void _atomicPush(::std::stack<TEntry> &stack, TStackEntry const &entry)	{	stack.push(entry);		}

	template < typename TEntry >
	inline void _atomicPop(::std::stack<TEntry> &stack)								{	stack.pop();			}

	template < typename TEntry, typename TStackEntry >
	inline void _atomicTop(::std::stack<TEntry> const &stack, TStackEntry &entry)	{	entry = stack.top();	}

	template < typename TEntry >
	inline TEntry const & _atomicTop(::std::stack<TEntry> const &stack)				{	return stack.top();		}


	template < typename TEntry, typename TStackEntry >
	inline void _atomicPush(::std::vector<TEntry> &stack, TStackEntry const &entry)	{	stack.push_back(entry);	}

	template < typename TEntry >
	inline void _atomicPop(::std::vector<TEntry> &stack)							{	stack.pop_back();		}

	template < typename TEntry, typename TStackEntry >
	inline void _atomicTop(::std::vector<TEntry> const &stack, TStackEntry &entry)	{	entry = stack.back();	}

	template < typename TEntry >
	inline TEntry const & _atomicTop(::std::vector<TEntry> const &stack)			{	return stack.back();	}


        
    //////////////////////////////////////////////////////////////////////////////
	// navigate through virtual suffix tree
       
	template < typename TIndex, class TSpec, typename TStackEntry >
	inline void _historyPush(Iter< TIndex, VSTree<TSpec> > &it, TStackEntry value) {}

	template < typename TIndex, class TSpec, typename TStackEntry >
	inline void _historyPush(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, TStackEntry value) {
		_atomicPush(it.history, value);
	}

	template < typename TIndex, class TSpec, typename TStackEntry >
	inline void _historyTop(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > const &it, TStackEntry &value) {
		value = _atomicTop(it.history);
	}

	template < typename TIndex, class TSpec, typename TStackEntry >
	inline void _historyPop(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it, TStackEntry &value) {
		value = _atomicTop(it.history);
		_atomicPop(it.history);
	}

	template < typename TIndex, class TSpec >
	inline void _historyPop(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		_atomicPop(it.history);
	}

/*
	template < typename TIndex, class TSpec >
	inline Pair<typename Size<TIndex>::Type> _historyPop(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		Pair<typename Size<TIndex>::Type> p = _atomicTop(it);
		_atomicPop(it.history);
		return p;
	}
*/


///.Function.goNext.param.iterator.type:Spec.BottomUp Iterator
///.Function.goNext.param.iterator.type:Spec.TopDownHistory Iterator

	template < typename TText, typename TSpecIndex, typename TSpecIter >
	inline void goNext(Iter<Index<TText, TSpecIndex>, TSpecIter> &it) {
		goNext(it, typename DefaultDFSOrder< Iter<Index<TText, TSpecIndex>, TSpecIter> >::Type());
	}

	template < typename TText, typename TSpecIndex, typename TSpecIter, typename TSpecTag >
	inline void goNext(Iter<Index<TText, TSpecIndex>, TSpecIter> &it, Tag<TSpecTag> const tagOrder) {
		goNextImpl(it, tagOrder);
	}


/**
.Function.goDown:
..summary:Iterates down one edge or a path in a tree.
..cat:Index
..signature:bool goDown(iterator)
..signature:bool goDown(iterator, char)
..signature:bool goDown(iterator, text[, lcp])
..param.iterator:An iterater of a Suffix Tree.
...type:Spec.TopDown Iterator
..param.char:$iterator$ goes down the edge beginning with $char$.
..param.text:$iterator$ goes down the path representing $text$. If $text$ ends within an edge, $iterator$ will point to the child-end of this edge.
..param.lcp:A reference of a size type. When $goDown$ returns, $lcp$ contains the length of the longest-common-prefix of $text$ and a path beginning at the $iterator$ node.
...type:Class.String
...type:Class.Segment
..remarks:$goDown(iterator)$ goes down the leftmost edge in the Suffix Tree, i.e. the edge beginning with the lexicographically smallest character.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
*/

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
	inline bool _goDownString(Iter< Index<TText, Index_ESA<> >, VSTree< TopDown<TSpec> > > &node, TString const &pattern, TSize &lcp) {
		typedef typename Infix< typename Fibre<TIndex, ESA_RawText>::Type const >::Type	TInfix;
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
..param.iterator:An iterater of a Suffix Tree.
...type:Spec.TopDownHistory Iterator
*/

	// go up one edge (returns false if in root node)
	template < typename TIndex, class TSpec >
	inline bool goUp(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		if (!isRoot(it)) {
			_historyPop(it, value(it));
			return true;
		}
		return false;
	}

/**
.Function.goRight:
..summary:Iterates to the next sibling in a tree.
..cat:Index
..signature:goRight(iterator)
..param.iterator:An iterater of a Suffix Tree.
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
				value(it).i2 = it.topRight.i2;
				return true;
			}
		}
		return false;
	}

	template < typename TIndex, class TSpec >
	inline bool goRight(Iter< TIndex, VSTree< TopDown< ParentLinks<TSpec> > > > &it) {
		if (value(it).i2 == length(container(it))) return false;		// right-most intervals have no right (would cause trouble with i2=\infty

		if (!isRoot(it)) {
			Pair<typename Size<TIndex>::Type> top;
			_historyTop(it, top);
			if (value(it).i2 < top.i2) {
				if (_isNextl(value(it).i2, container(it))) {
					value(it).i1 = value(it).i2;
					value(it).i2 = _getNextl(value(it).i2, container(it));
					return true;
				}
				value(it).i1 = value(it).i2;
				value(it).i2 = top.i2;
				return true;
			}
		}
		return false;
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
..param.iterator:An iterater of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to the root of the tree, otherwise $false$.
*/

	template < typename TIndex, class TSpec >
	inline bool isRoot(Iter<TIndex, VSTree< BottomUp<TSpec> > > const &it) {
		return it.history.empty();
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
..param.iterator:An iterater of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to the node representing a suffix, otherwise $false$.
..remarks:Every leaf is also a right terminal (see @Function.isLeaf@), but not vice versa.
*/

	template < typename TIndex, class TSpec >
	inline bool isRightTerminal(Iter<TIndex, VSTree<TSpec> > const &it) {
		// do we reach a leaf in a suffix tree with trailing '$'
		return (*it).i2 == length(container(it));
	}

/**
.Function.childrenAreLeaves:
..summary:Test whether iterator points to a node with only leaf-children.
..cat:Index
..signature:bool childrenAreLeaves(iterator)
..param.iterator:An iterater of a Suffix Tree.
...type:Spec.VSTree Iterator
..returns:$true$ if $iterator$ points to an inner node of the tree, whose children are leaves. Otherwise it is $false$.
*/

	template < typename TIndex, class TSpec >
	inline bool childrenAreLeaves(Iter<TIndex, VSTree<TSpec> > const &it) {
		return countChildren(it) == (value(it).i2 - value(it).i1);
	}

/**
.Function.isLeaf:
..summary:Test whether iterator points to a leaf.
..cat:Index
..signature:bool isLeaf(iterator)
..param.iterator:An iterater of a Suffix Tree.
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
