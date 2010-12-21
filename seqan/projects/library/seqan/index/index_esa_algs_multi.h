// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_ESA_ALGS_MULTI_H
#define SEQAN_HEADER_INDEX_ESA_ALGS_MULTI_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	//////////////////////////////////////////////////////////////////////////////
	// more sophisticated algorithms on suffix trees of
	// multiple sequences (generalized suffix tree)
	//////////////////////////////////////////////////////////////////////////////


/**
.Spec.MUMs Iterator:
..cat:Index
..general:Spec.BottomUp Iterator
..summary:Iterator to search for all maximum unique matches.
..signature:Iterator<TContainer, MUMs>::Type
..signature:Iter<TContainer, VSTree< BottomUp<MUMs> > >
..param.TContainer:Type of an index that can be iterated with a bottom-up iterator.
...type:Spec.IndexEsa
...metafunction:Metafunction.Container
..include:seqan/index.h

.Memfunc.MUMs Iterator#Iterator
..class:Spec.MUMs Iterator
..summary:Constructor
..signature:Iterator(index[, minLength])
..signature:Iterator(iterator)
..param.index:An index object.
..param.minLength:Minimum length of the maximal unique matches, default value is 1.
..param.iterator:Another MUMs iterator. (copy constructor)
...type:Spec.MUMs Iterator
*/

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
		typedef VectorSet_<TSize, Alloc<> >				TSeqSet;
//____________________________________________________________________________

		TSize		minLength;
		TSize		seqCount;
		TSeqSet		seqSet;
//____________________________________________________________________________

		Iter(TSTree &_tree):
			TBase(_tree),
			minLength(1),
			seqCount(countSequences(_tree)),
			seqSet(countSequences(_tree))
		{
			indexRequire(_tree, EsaBwt());
			goNext(*this);	// the iterator starts in a suffix, i.e. not a MUM node (length(occ)<2<=seqCount)
		}

		Iter(TSTree &_tree, MinimalCtor):
			TBase(_tree, MinimalCtor()) {}

		Iter(TSTree &_tree, TSize _minLength):
			TBase(_tree),
			minLength(_minLength),
			seqCount(countSequences(_tree)),
			seqSet(countSequences(_tree))
		{
			indexRequire(_tree, EsaBwt());
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
			     !(	(countOccurrences(it) == it.seqCount) && 
					(repLength(it) >= it.minLength) &&
					isUnique(it, it.seqSet) && 
					isLeftMaximal(it)) );
	}



	
/**
.Spec.MultiMEMs Iterator:
..cat:Index
..general:Spec.BottomUp Iterator
..summary:Iterator to search for MultiMEMs.
..signature:Iterator<TContainer, MultiMEMs>::Type
..signature:Iter<TContainer, VSTree< BottomUp<MultiMEMs> > >
..param.TContainer:Type of an index that can be iterated with a bottom-up iterator.
...type:Spec.IndexEsa
...metafunction:Metafunction.Container
..include:seqan/index.h

.Memfunc.MultiMEMs Iterator#Iterator
..class:Spec.MultiMEMs Iterator
..summary:Constructor
..signature:Iterator(index[, minLength])
..signature:Iterator(iterator)
..param.index:An index object.
..param.minLength:Minimum length of the multiple maximal exact matches, default value is 1.
..param.iterator:Another MultiMEMs iterator. (copy constructor)
...type:Spec.MultiMEMs Iterator
*/

	//////////////////////////////////////////////////////////////////////////////
	// MultiMEMs
	//////////////////////////////////////////////////////////////////////////////

	// contains a set of fraction compounds
	// one compound for each sequence
	template <typename TValue, typename TSize>
	struct FractionMultiCompound_ {
		typedef FractionCompound_<TValue, TSize>	TCompound;
		typedef String<TCompound>					TSet;	// seqNo..unsigned, compound: bwt character->suffixes

		TSet	set;
	};


	template < typename TSTree >
	class Iter< TSTree, VSTree< BottomUp<MultiMEMs> > >:
		public Iter< TSTree, VSTree< BottomUp<> > > 
	{
	public:
		typedef Iter< TSTree, VSTree< BottomUp<> > >	TBase;
		typedef typename Value<TSTree>::Type			TValue;
		typedef typename Size<TSTree>::Type				TSize;

		typedef FractionMultiCompound_<TValue, TSize>	TMultiCompound;
		typedef String<TMultiCompound, Block<> >		TSetStack;
		typedef String<TSize>							TPositionList;
		
		typedef typename TMultiCompound::TSet			TSet;
		typedef typename Iterator<TSet>::Type			TSetIterator;

		typedef typename TBase::TStackEntry				TStackEntry;

//____________________________________________________________________________

		TSize			minLength;
		unsigned		minSupport;	// the support is the number of distinct sequences
		unsigned		maxSupport;	// a repeat/match occurs in
		TSetStack		setStack;
		TPositionList	posList;	// this list is indexed just as SA is and contains the next entry's index
		bool			canMerge;	// is false, if parent node appears after its first child on stack
//____________________________________________________________________________

		Iter(TSTree &_index):
			TBase(_index, MinimalCtor()),
			minSupport(countSequences(_index)),
			maxSupport(countSequences(_index)),
			canMerge(true)
		{
			indexRequire(_index, EsaSA());
			indexRequire(_index, EsaLcp());
			indexRequire(_index, EsaBwt());
			resize(posList, length(_index));

			if (!empty(indexSA(_index))) 
			{
				TStackEntry e;
				e.range.i1 = 0;
				e.range.i2 = 0;
				_dfsOnPush(*this, e);
				goNext(*this);
			}
		}

		Iter(TSTree &_tree, MinimalCtor):
			TBase(_tree, MinimalCtor()) {}

		Iter(TSTree &_index, TSize _minLength):
			TBase(_index, MinimalCtor()),
			minLength(_minLength),
			minSupport(countSequences(_index)),
			maxSupport(countSequences(_index)),
			canMerge(true)
		{
			indexRequire(_index, EsaSA());
			indexRequire(_index, EsaLcp());
			indexRequire(_index, EsaBwt());
			resize(posList, length(_index));

			if (!empty(indexSA(_index))) 
			{
				TStackEntry e;
				e.range.i1 = 0;
				e.range.i2 = 0;
				_dfsOnPush(*this, e);
				goNext(*this);
			}
		}

		Iter(Iter const &_origin):
			TBase((TBase const &)_origin),
			minLength(_origin.minLength),
			minSupport(_origin.minSupport),
			maxSupport(_origin.maxSupport),
			setStack(_origin.setStack),
			posList(_origin.posList),
			canMerge(_origin.canMerge) {}

//____________________________________________________________________________

		inline bool hasRepeats() 
		{
			if (length(setStack) < 2) return false;

			TMultiCompound &child  = top(setStack);
			TMultiCompound &parent = topPrev(setStack);

			TValue prevKey = TValue();
			TValue equalKey = TValue();

			unsigned distinctSeqs = 0;
			bool distinctKeys = false;

			TSetIterator parentCompound	= begin(parent.set);
			TSetIterator parentEnd		= end(parent.set);
			TSetIterator childCompound	= begin(child.set);
			TSetIterator childEnd		= end(child.set);

			while (childCompound != childEnd && parentCompound != parentEnd)
			{
				int result = _haveMaximalRepeats(*childCompound, *parentCompound, equalKey);
				if (result > 0) {
					if (!distinctKeys && result == 1) {
						if (distinctSeqs > 0 && prevKey != equalKey)
							distinctKeys = true;						// there is a left maximal repeat
						prevKey = equalKey;
					} else
						distinctKeys = true;

					if (++distinctSeqs > minSupport && distinctKeys)	// if it is also a  repeat in at least
						return true;									// minSequences distinct sequences then
				}														// we have at least one repeat
				++childCompound;
				++parentCompound;
			}
			return false;
		}
/*
		inline TSize countRepeats() 
		{
			if (length(setStack) < 2) return 0;

			TFractionCompound &child  = top(setStack);
			TFractionCompound &parent = topPrev(setStack);

			TSetIterator childFraction	= begin(child.set);
			TSetIterator childEnd		= end(child.set);
			TSetIterator parentFraction	= begin(parent.set);
			TSetIterator parentEnd		= end(parent.set);

			TSize sum = 0;
			for(; childFraction != childEnd; ++childFraction) {
				for(; parentFraction != parentEnd; ++parentFraction) {
					if (keyOf(childFraction) != keyOf(parentFraction))
						sum += (*childFraction).size * (*parentFraction).size;

					sum += child.leftmost.size * (*parentFraction).size;
				}
				sum += (*childFraction).size * parent.leftmost.size;
			}
			sum += child.leftmost.size * parent.leftmost.size;
			return sum;
		}
*/
//____________________________________________________________________________
/*
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
*/
	};

	// add bwt partitions of child to parent node
	template < typename TSTree, typename TSpec, typename TValue, typename TSize >
	inline void _fractionMerge(
		Iter<TSTree, VSTree< BottomUp<TSpec> > > &it, 
		FractionMultiCompound_<TValue, TSize> &parent,
		FractionMultiCompound_<TValue, TSize> &child)
	{
		typedef FractionMultiCompound_<TValue, TSize>	TCompound;
		typedef typename TCompound::TSet				TSet;
		typedef typename Iterator<TSet, Standard>::Type	TSetIterator;

		TSetIterator parentCompound	= begin(parent.set, Standard());
		TSetIterator parentEnd		= end(parent.set, Standard());
		TSetIterator childCompound	= begin(child.set, Standard());
		TSetIterator childEnd		= end(child.set, Standard());

		while (childCompound != childEnd && parentCompound != parentEnd)
		{
			// append child compound to parent compound
			_fractionMerge(it, *parentCompound, *childCompound);
			++childCompound;
			++parentCompound;
		}
	}

	template < typename TSTree >
	inline void _dfsOnLeaf(Iter<TSTree, VSTree< BottomUp<MultiMEMs> > > &it) 
	{
		typedef Iter<TSTree, VSTree< BottomUp<> > > TBase;
		_dfsOnLeaf((TBase&)it);

		typedef typename Value<TSTree>::Type		TValue;
		typedef typename Size<TSTree>::Type			TSize;
		typedef typename SAValue<TSTree>::Type		TSAValue;
		typedef FractionHeader_<TSize>				TFractionHeader;
		typedef Pair<TValue, TFractionHeader>		TFraction;
		typedef typename Set<TFraction>::Type		TFractionSet;

		typedef FractionCompound_<TValue, TSize>	TCompound;
		typedef Pair<unsigned, TCompound>			TCompoundPair;
		typedef typename Set<TCompoundPair>::Type	TSet;

		push(it.setStack);

		TSTree &index = container(it);

		TSize		gPos = posGlobalize(_dfsRange(it).i1, stringSetLimits(index));
		TSAValue	lPos;
		posLocalize(lPos, _dfsRange(it).i1, stringSetLimits(index));

		TCompound &compound = top(it.setStack).set[getValueI1(lPos)];
		if (!posAtFirstLocal(lPos))
			insert(
				TFraction(
					bwtAt(gPos, container(it)),
					TFractionHeader(gPos, gPos, 1)), 
				compound.set);
		else
			compound.leftmost = TFractionHeader(gPos, gPos, 1);

		_setSizeInval(it.posList[gPos]);
/*
		::std::cerr << "LEAF ";
		_dumpHistoryStack(it);
		it._dump();
*/	}



	//////////////////////////////////////////////////////////////////////////////
	// maximal repeat representation
	//////////////////////////////////////////////////////////////////////////////

	template <typename TSTree>
	struct MultiMEM {
//		Iter< TSTree, VSTree<BottomUp<MultiMEMs> > > &it;
	};

	template <typename TSTree>
	struct Value< MultiMEM<TSTree> > {
		typedef Pair< typename SAValue<TSTree>::Type > Type;
	};

	template <typename TSTree>
	struct Size< MultiMEM<TSTree> > {
		typedef typename Size<TSTree>::Type Type;
	};


	template <typename TSTree>
	inline typename Size< MultiMEM<TSTree> >::Type 
	length(MultiMEM<TSTree> const &repeat) {
		return repeat.it.countRepeats();
	}

/*
	template <typename TSTree>
	inline typename Iterator< MultiMEM<TSTree> >::Type 
	begin(MultiMEM<TSTree> &repeat) {
		return Iterator< MultiMEM<TSTree> >::Type(repeat.it);
	}

	template <typename TSTree>
	inline typename Iterator< MultiMEM<TSTree> const >::Type 
	begin(MultiMEM<TSTree> const &repeat) {
		return Iterator< MultiMEM<TSTree> >::Type(repeat.it);
	}
*/


	template <typename TSTree>
	class Iter< MultiMEM<TSTree>, MultiMEMOccurences > {
	public:

		typedef typename Value<TSTree>::Type	TValue;
		typedef typename Size<TSTree>::Type		TSize;
		typedef typename SAValue<TSTree>::Type	TSAValue;
		typedef	Pair<TSAValue>					TPair;

		typedef FractionCompound_<TValue, TSize> const	TFractionCompound;
		typedef typename TFractionCompound::TSet const	TSet;
		typedef typename Iterator<TSet>::Type			TSetIterator;

		typedef Iter<TSTree, VSTree<BottomUp<MultiMEMs> > >	const	TIterator;
		typedef typename TIterator::TPositionList const				TPositionList;

		TIterator	*mmemIt;
		bool		_atEnd;
		unsigned	seqCount;
		TPair		tmp;
	
		
		// for every sequence there is a SubState stucture
		// storing the current enumeration state
		// that is necessary to enumerate every combination
		// of left maximal multi match

		struct SubState 
		{
			SubState			*prevState;
			TPositionList		*posList;
			TFractionCompound	*child, *parent;
			TSize				childPtr, parentPtr;			// per seq. suffix iterators
			TSetIterator		childFraction,  childBegin,  childEnd;
			TSetIterator		parentFraction, parentBegin, parentEnd;	
			bool				leftmostChild, leftmostParent;	// use undef. bwt (leftmost) set
			TValue				leftChar;						// are the keys of seq[1..i] equal?
			bool				charsEqual;

			inline void _updateLeftMaximality() {
				if (prevState)
					charsEqual = prevState->charsEqual && (leftChar == keyOf(childFraction));
			}

			inline bool _innerStep()
			{
				if (_isSizeInval(childPtr = ((*posList)[childPtr]))) {
					if (_isSizeInval(parentPtr = (*posList)[parentPtr])) {
						parentPtr = objectOf(parentFraction).begin;
						return false;
					}
					childPtr = objectOf(childFraction).begin;
				}
				return true;
			}

			inline void _firstParentFraction()
			{
				parentBegin		= parentFraction	= begin(parent->set);
				parentEnd		= end(parent->set);

				if (parentFraction != parentEnd) {
					leftmostParent = false;
					parentPtr = objectOf(parentFraction).begin;
				} else {
					leftmostParent = true;
					parentPtr = parent->leftmost.begin;
				}

				leftChar = keyOf(parentFraction);
			}

			inline void _firstChildFraction() 
			{
				childBegin		= childFraction		= begin(child->set);
				childEnd		= end(child->set);

				if (childFraction != childEnd) {
					leftmostChild = false;
					childPtr = objectOf(childFraction).begin;
				} else {
					leftmostChild = true;
					childPtr = child->leftmost.begin;
				}
			}

			inline bool _nextParentFraction() 
			{
				if (leftmostParent) return false;

				if (++parentFraction == parentEnd) {
					if (parent->leftmost.size > 0) {
						leftmostParent = true;
						parentPtr = parent->leftmost.begin;
					} else
						return false;
				} else
					parentPtr = objectOf(parentFraction).begin;

				return true;
			}

			inline bool _nextChildFraction() 
			{
				if (leftmostChild) return false;

				if (++childFraction == childEnd) {
					if (child->leftmost.size > 0) {
						leftmostChild = true;
						childPtr = child->leftmost.begin;
					} else
						return false;
				} else
					childPtr = objectOf(childFraction).begin;

				return true;
			}

			// single per-sequence enumeration step
			inline bool _outerStep() 
			{
				if (!_nextChildFraction()) {
					_firstChildFraction();
					if (!_nextParentFraction()) {
						_firstParentFraction();
						return false;
					}
				}
				return true;
			}

		};

		String<SubState>	subState;


		inline Iter(Iter<TSTree, VSTree<BottomUp<MultiMEMs> > > const &_maxIt):
			mmemIt(&_maxIt),
			seqCount(countSequences(container(_maxIt)))
		{
			_init();
		}

		inline bool _isLeftMaximal() {
			return !subState[seqCount - 1].charsEqual;
		}
		
		// inner enumeration
		inline bool _innerStep() {
			for(unsigned seq = 0; seq < seqCount; ++seq)
				if (subState[seq]._innerStep()) return true;
			return false;
		}

		inline bool _outerStepNoCheck() {
			for(unsigned seq = 0; seq < seqCount; ++seq)
				if (subState[seq]._outerStep()) return true;
			return false;
		}

		// outer enumeration
		inline bool _outerStep() {
			while (!((_atEnd = !_outerStepNoCheck()) || _isLeftMaximal())) ;
			return !_atEnd;
		}

		inline void _init() 
		{
			if (length(mmemIt->setStack) < 2) {
				_atEnd = true;
				return;
			}

			resize(subState, seqCount);
			SubState *prev = NULL;
			for(unsigned seq = 0; seq < seqCount; ++seq)
			{
				SubState &state = subState[seq];

				state.prevState = prev;
				state.posList = &(mmemIt->posList);
				state.parent = &(topPrev(mmemIt->setStack).set[seq]);
				state.child = &(top(mmemIt->setStack).set[seq]);

				state._firstChildFraction();
				state._firstParentFraction();

				prev = &state;
			}

			_atEnd = false;
			while (!(_isLeftMaximal() || (_atEnd = !_outerStep()))) ;
		}
	};


	template < typename TRepeat >
	inline typename Value< Iter<TRepeat, MultiMEMOccurences> >::Type &
	value(Iter<TRepeat, MultiMEMOccurences> const &it)  {
		return it.tmp;
	}

	template < typename TRepeat >
	inline typename Value< Iter<TRepeat, MultiMEMOccurences> >::Type &
	value(Iter<TRepeat, MultiMEMOccurences> &it)  {
		return it.tmp;
	}

//TODO:fix me
	template < typename TRepeat >
	inline Iter<TRepeat, MultiMEMOccurences> &
	goNext(Iter<TRepeat, MultiMEMOccurences> &it)  {
		if (it._innerStep()) {
//			it.tmp.i1 = saAt(it.subState.parentPtr, container(*it.mmemIt));
//			it.tmp.i2 = saAt(it.subState.childPtr, container(*it.mmemIt));
			return it;
		}
		if (it._outerStep()) {
//			it.tmp.i1 = saAt(it.subState.parentPtr, container(*it.mmemIt));
//			it.tmp.i2 = saAt(it.subState.childPtr, container(*it.mmemIt));
		}
		return it;
	}

	template < typename TRepeat >
	inline bool atEnd(Iter<TRepeat, MultiMEMOccurences> const &it) {
		return it._atEnd;
	}

	template < typename TRepeat >
	inline bool atEnd(Iter<TRepeat, MultiMEMOccurences> &it) {
		return it._atEnd;
	}


	template <typename TSTree>
	struct Iterator< MultiMEM<TSTree> > {
		typedef Iter<MultiMEM<TSTree>, MultiMEMOccurences> Type;
	};

	template <typename TSTree>
	struct Size< Iter<MultiMEM<TSTree>, MultiMEMOccurences> > {
		typedef typename Size<TSTree>::Type Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// Iterator wrappers
	//////////////////////////////////////////////////////////////////////////////

	template <typename TObject>
	struct Iterator< TObject, MUMs > {
		typedef Iter< TObject, VSTree< BottomUp<MUMs> > > Type;
	};

//}

}

#endif
