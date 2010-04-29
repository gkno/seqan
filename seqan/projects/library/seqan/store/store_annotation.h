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

#ifndef SEQAN_HEADER_STORE_ANNOTATION_H
#define SEQAN_HEADER_STORE_ANNOTATION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Annotation Store
//////////////////////////////////////////////////////////////////////////////

template <typename _TPos, typename TSpec = void>
struct AnnotationStoreElement
{
	typedef typename Id<AnnotationStoreElement>::Type		TId;
	typedef _TPos											TPos;
	typedef StringSet<CharString, Owner< ConcatDirect<> > >	TValues;

	static const TId  INVALID_ID;
	static const TPos INVALID_POS;
	
	TId					parentId;
	TId					contigId;
	TId					countId;
	TId					typeId;			// gene, intron, ...

	TPos				beginPos;		// begin position of the gapped sequence in gapped contig sequence
	TPos				endPos;			// end position of ..., for reverse aligned reads holds end < begin
	
	TId					lastChildId;	// generated back links to child
	TId					nextSiblingId;	// and sibling
	
	TValues				values;			// stores values for each keyId of (key,value) pairs

	AnnotationStoreElement() : 
		parentId(INVALID_ID), contigId(INVALID_ID), countId(INVALID_ID), typeId(INVALID_ID), 
		beginPos(INVALID_POS), endPos(INVALID_POS),
		lastChildId(INVALID_ID), nextSiblingId(INVALID_ID) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TSpec>
const typename Id<AnnotationStoreElement<TPos, TSpec> >::Type
AnnotationStoreElement<TPos, TSpec>::INVALID_ID = SupremumValue<typename Id<AnnotationStoreElement<TPos, TSpec> >::Type>::VALUE;

template <typename TPos, typename TSpec>
const TPos
AnnotationStoreElement<TPos, TSpec>::INVALID_POS = SupremumValue<TPos>::VALUE;

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct AnnotationTree {};

template <typename TFragmentStore, typename TSpec>
class Iter<TFragmentStore, AnnotationTree<TSpec> >
{
public:
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TFragmentStore *store;
	TId id;
	
	Iter():
		store(NULL),
		id(TAnnotation::INVALID_ID) {}

	Iter(TFragmentStore &_store):
		store(&_store),
		id(0) {}

	Iter(TFragmentStore &_store, MinimalCtor):
		store(&_store),
		id(TAnnotation::INVALID_ID) {}

	inline Iter const &
	operator = (Iter const &_origin)
	{
		store = &container(_origin);
		id = _origin.id;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
struct Iterator< TFragmentStore, AnnotationTree<TSpec> > {
	typedef Iter< TFragmentStore, AnnotationTree<TSpec> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
struct Value< Iter< TFragmentStore, AnnotationTree<TSpec> > >:
	VertexDescriptor<TFragmentStore> {};

template <typename TFragmentStore, typename TSpec>
struct Size< Iter< TFragmentStore, AnnotationTree<TSpec> > > :
	Size<TFragmentStore> {};

template <typename TFragmentStore, typename TSpec>
struct Position< Iter< TFragmentStore, AnnotationTree<TSpec> > > :
	Position<TFragmentStore> {};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline typename VertexDescriptor<TFragmentStore>::Type & 
value(Iter< TFragmentStore, AnnotationTree<TSpec> > &it) { 
	return it.id;
}

template <typename TFragmentStore, typename TSpec>
inline typename VertexDescriptor<TFragmentStore>::Type const & 
value(Iter< TFragmentStore, AnnotationTree<TSpec> > const &it) { 
	return it.id;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline TFragmentStore & 
container(Iter< TFragmentStore, AnnotationTree<TSpec> > &it) { 
	return *it.store;
}

template <typename TFragmentStore, typename TSpec>
inline TFragmentStore & 
container(Iter< TFragmentStore, AnnotationTree<TSpec> > const &it) { 
	return *it.store;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline typename Reference<typename TFragmentStore::TAnnotationStore>::Type
getAnnotation(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return it.store->annotationStore[it.id];
}

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
getName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoName(*it.store, it.id);
}

template <typename TFragmentStore, typename TSpec, typename TName>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
setName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it, TName & name)
{
	if (length(it.store->annotationNameStore) <= it.id)
		resize(it.store->annotationNameStore, it.id + 1);
	it.store->annotationNameStore[it.id] = name;
}

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
getParentName(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TId parentId = it.store->annotationStore[it.id].parentId;
	if (parentId == TAnnotation::INVALID_ID) parentId = it.id;
	return getAnnoName(*it.store, parentId);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationTypeStore>::Type
getType(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoType(*it.store, it.id);
}

template <typename TFragmentStore, typename TSpec, typename TTypeName>
inline void
setType(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it, TTypeName & typeName)
{
	_storeAppendType(*it.store, getAnnotation(it).typeId, typeName);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline CharString
getUniqueName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoUniqueName(*it.store, it.id);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void 
clearValues(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	clear(getAnnotation(it).values);
}

template <typename TFragmentStore, typename TSpec, typename TKey, typename TValue>
inline void 
assignValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > & it,
	TKey const & key,
	TValue const & value)
{
	annotationAssignValueByKey(*it.store, getAnnotation(it), key, value);
}

template <typename TFragmentStore, typename TSpec, typename TKey, typename TValue>
inline bool 
getValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > const & it,
	TKey const & key,
	TValue & value)
{
	return annotationGetValueByKey(*it.store, getAnnotation(it), key, value);
}

template <typename TFragmentStore, typename TSpec, typename TKey>
inline CharString
getValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > const & it,
	TKey const & key)
{
	return annotationGetValueByKey(*it.store, getAnnotation(it), key);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void
goBegin(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	it.id = 0;
}

template <typename TFragmentStore, typename TSpec>
inline void
goEnd(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	it.id = TAnnotation::INVALID_ID;
}

template <typename TFragmentStore, typename TSpec>
inline void
clear(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	it.id = TAnnotation::INVALID_ID;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline bool
atBegin(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	return it.id == 0;
}

template <typename TFragmentStore, typename TSpec>
inline bool
atEnd(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	return it.id == TAnnotation::INVALID_ID;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void
goNext(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	if (!goDown(it) && !goRight(it))
		while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

template <typename TFragmentStore, typename TSpec>
inline void
goNextRight(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	if (!goRight(it))
		while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

template <typename TFragmentStore, typename TSpec>
inline void
goNextUp(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline void
goRoot(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	it.id = 0;
}

template <typename TFragmentStore, typename TSpec>
inline bool
goUp(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TId parentId = getAnnotation(it).parentId;
	if (parentId != TAnnotation::INVALID_ID)
	{
		it.id = parentId;
		return true;
	}
	return false;
}

template <typename TFragmentStore, typename TSpec>
inline bool
goDown(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TId lastChildId = getAnnotation(it).lastChildId;
	if (lastChildId != TAnnotation::INVALID_ID)
	{
		it.id = it.store->annotationStore[lastChildId].nextSiblingId;
		return true;
	}
	return false;
}

template <typename TFragmentStore, typename TSpec>
inline bool
goRight(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TAnnotation &anno = getAnnotation(it);
	TId nextSiblingId = anno.nextSiblingId;
	if (nextSiblingId != TAnnotation::INVALID_ID)
	{
		TId lastChildId = it.store->annotationStore[anno.parentId].lastChildId;
		if (it.id != lastChildId)
		{
			it.id = nextSiblingId;
			return true;
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeUp(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goUp(tmp);
	return tmp;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeDown(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goDown(tmp);
	return tmp;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeRight(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goRight(tmp);
	return tmp;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createLeftChild(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &anno = getAnnotation(it);

	TId childId = length(it.store->annotationStore) - 1;
	TAnnotation &childAnno = it.store->annotationStore[childId];
	
	TId firstChildId;
	TId lastChildId = anno.lastChildId;
	
	if (lastChildId != TAnnotation::INVALID_ID)
	{
		TAnnotation &lastChild = it.store->annotationStore[lastChildId];
		firstChildId = lastChild.nextSiblingId;
		lastChild.nextSiblingId = childId;
	} else
		firstChildId = childId;

	childAnno.nextSiblingId = firstChildId;
	childAnno.parentId = it.id;
	childAnno.lastChildId = TAnnotation::INVALID_ID;
	it.store->annotationStore[lastChildId].nextSiblingId = childId;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > childIter(it);
	childIter.id = childId;
	return childIter;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createRightChild(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &parentAnno = getAnnotation(it);

	TId childId = length(it.store->annotationStore) - 1;
	TAnnotation &childAnno = it.store->annotationStore[childId];
	
	TId firstChildId;
	TId lastChildId = parentAnno.lastChildId;
	
	if (lastChildId != TAnnotation::INVALID_ID)
	{
		TAnnotation &lastChild = it.store->annotationStore[lastChildId];
		firstChildId = lastChild.nextSiblingId;
		lastChild.nextSiblingId = childId;
	} else
		firstChildId = childId;

	childAnno.nextSiblingId = firstChildId;
	childAnno.parentId = it.id;
	childAnno.lastChildId = TAnnotation::INVALID_ID;
	parentAnno.lastChildId = childId;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > childIter(it);
	childIter.id = childId;
	return childIter;
}

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createSibling(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &anno = getAnnotation(it);

	TId siblingId = length(it.store->annotationStore) - 1;

	TAnnotation &parentAnno = it.store->annotationStore[anno.parentId];
	if (parentAnno.lastChildId == it.id)
		parentAnno.lastChildId = siblingId;

	TAnnotation &siblingAnno = it.store->annotationStore[siblingId];
	siblingAnno.nextSiblingId = anno.nextSiblingId;
	siblingAnno.parentId = anno.parentId;
	siblingAnno.lastChildId = TAnnotation::INVALID_ID;
	anno.nextSiblingId = siblingId;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > siblingIter(it);
	siblingIter.id = siblingId;
	return siblingIter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline bool
isRoot(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

//	if (it.id >= length(it.store->annotationStore)) return false;
	return it.store->annotationStore[it.id].parentId == TAnnotation::INVALID_ID;
}

template <typename TFragmentStore, typename TSpec>
inline bool
isLeaf(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

//	if (it.id >= length(it.store->annotationStore)) return false;
	return it.store->annotationStore[it.id].lastChildId == TAnnotation::INVALID_ID;
}

template <typename TFragmentStore, typename TSpec>
inline bool
isLastChild(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TAnnotation &anno = getAnnotation(it);
	TId nextSiblingId = anno.nextSiblingId;
	if (nextSiblingId != TAnnotation::INVALID_ID)
	{
		TId lastChildId = it.store->annotationStore[anno.parentId].lastChildId;
		return it.id == lastChildId;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAnnotationStore>
inline void
_storeClearAnnoBackLinks(TAnnotationStore & me)
{
	typedef typename Value<TAnnotationStore>::Type				TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type TAnnoIter;

	TAnnoIter it = begin(me, Standard());
	TAnnoIter itEnd = end(me, Standard());
	
	for (; it != itEnd; ++it)
	{
		(*it).lastChildId = TAnnotation::INVALID_ID;
		(*it).nextSiblingId = TAnnotation::INVALID_ID;
	}
}

template <typename TAnnotationStore>
inline void
_storeCreateAnnoBackLinks(TAnnotationStore & me)
{
	typedef typename Value<TAnnotationStore>::Type				TAnnotation;
	typedef typename TAnnotation::TId							TId;
	typedef typename Iterator<TAnnotationStore, Standard>::Type TAnnoIter;
	
	TAnnoIter itBegin = begin(me, Standard());
	TAnnoIter itEnd = end(me, Standard());
	TId id = (itEnd - itBegin) - 1;
	TAnnoIter it = itBegin + id;
	
	for (; itBegin <= it; --it, --id)
	{
		if ((*it).parentId != TAnnotation::INVALID_ID)
		{
			TAnnoIter parent = itBegin + (*it).parentId;
			if ((*parent).lastChildId == TAnnotation::INVALID_ID)
			{
				(*parent).lastChildId = id;
				(*it).nextSiblingId = id;
			}

			if ((*it).nextSiblingId == TAnnotation::INVALID_ID)
			{
				TAnnoIter lastChild = itBegin + (*parent).lastChildId;
				(*it).nextSiblingId = (*lastChild).nextSiblingId;
				(*lastChild).nextSiblingId = id;
			}
		}
		else 
			(*it).nextSiblingId = TAnnotation::INVALID_ID;
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
