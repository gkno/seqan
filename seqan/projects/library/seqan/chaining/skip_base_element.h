/*	Copyright (c) 2006 Hendrik Woehrle
*	All rights reserved.
*
*	Deferred Skip List Datastructure
*
*	Elements in the base layer of the Skip List
*
*	Specializations:
*
*	TModus:
*		* Dynamic: contains predecessor/successor pointers to have the properties of an double linked list in the base layer
*		* Static: comes without these pointers to reduce size
*
*	TSpec:
*		* Complete: Not using Deferred Data Structuring
*		* Deferred: For use with Deferred Data Structuring
*
*/

#ifndef SEQAN_HEADER_SKIP_BASE_ELEMENT_H
#define SEQAN_HEADER_SKIP_BASE_ELEMENT_H

namespace seqan
{


//////////////////////////////////////////////////////////////////////////////
//
//		SkipBaseElement
//
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//
//	Dynamic <-> Static and Complete <-> Deferred member wrapper structs
//
//	depending on this specializations, the elements contain different members
//
//////////////////////////////////////////////////////////////////////////////

// special members of static/dynamic elements

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct
_DynamicStruct
{};

template< typename TObject, typename TSpec, typename TStructuring >
struct
_DynamicStruct< TObject, Dynamic, TSpec, TStructuring >
{
		// predecessing element in the skip list
	SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > * _pred;
		// succeeding element
	SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > * _succ;

	_DynamicStruct()
		: _pred( NULL )
		, _succ( NULL )
	{}

	_DynamicStruct( SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > * goPrevious, 
					SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > * succ )
		: _pred( goPrevious )
		, _succ( succ )
	{}
};

//////////////////////////////////////////////////////////////////////////////

// special members of complete/deferred elements

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct
_DeferredStruct
{
	union{
			// number of unsorted elements to the right of this
		typename Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type _count;
			// height of related tower
		typename Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type _height;
	};

	_DeferredStruct()
		: _height( 0 )
	{}

	_DeferredStruct( typename Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type h )
		: _height( h )
	{}
};

template< typename TObject, typename TSpec >
struct
_DeferredStruct< TObject, Dynamic, TSpec, Deferred >
{
		// next sorted element on the left side
	SkipBaseElement< TObject, Dynamic, TSpec, Deferred > * _left;
		// next sorted element on the right side
	SkipBaseElement< TObject, Dynamic, TSpec, Deferred > * _right;
		// number of unsorted elements to the right of this
	typename Size< SkipBaseElement< TObject, Dynamic, TSpec, Deferred > >::Type _count;
		// height of related tower
	typename Size< SkipBaseElement< TObject, Dynamic, TSpec, Deferred > >::Type _height;

	_DeferredStruct()
		: _left( NULL )
		, _right( NULL )
		, _count( 0 )
		, _height( 0 )
	{}

	_DeferredStruct( SkipBaseElement< TObject, Dynamic, TSpec, Deferred > * left, 
						SkipBaseElement< TObject, Dynamic, TSpec, Deferred > * right )
		: _left( left )
		, _right( right )
		, _count( 0 )
		, _height( 0 )
	{}

	template< typename TSize >
	_DeferredStruct( SkipBaseElement< TObject, Dynamic, TSpec, Deferred > * left, 
						SkipBaseElement< TObject, Dynamic, TSpec, Deferred > * right,
						TSize _count,
						TSize _height )
		: _left( left )
		, _right( right )
		, _count( _count )
		, _height( _height )
	{}
};


template< typename TObject, typename TSpec >
struct
_DeferredStruct< TObject, Static, TSpec, Deferred >
{
	//	// next sorted element on the left side
	//SkipBaseElement< TObject, Static, TSpec, Deferred > * _left;
		// number of unsorted elements to the right of this
	typename Size< SkipBaseElement< TObject, Static, TSpec, Deferred > >::Type _count;
		// height of related tower
	typename Size< SkipBaseElement< TObject, Static, TSpec, Deferred > >::Type _height;

	_DeferredStruct()
		: _count( 0 )
		, _height( 0 )
	{}

};

//////////////////////////////////////////////////////////////////////////////
//
// struct SkipBaseElement
//
//////////////////////////////////////////////////////////////////////////////

/*
.Internal.SkipBaseElement:
..summary:Type of elements in the base layer of a SkipList
..cat:SkipList
..signature:SkipBaseElement< TObject, TModus, TSpec, TStructuring >
.signature:SkipElement< TObject, TModus, TSpec, TStructuring >
..param.TObject:Type of the stored object.
..param.TModus:The TModus parameter of the SkipList.
..param.TSpec:The TSpec parameter of the SkipList.
..param.TStructuring:The TStructuring parameter of the SkipList.
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct SkipBaseElement
{
	
	_DynamicStruct< TObject, TModus, TSpec, TStructuring > _dynStruct;

	_DeferredStruct< TObject, TModus, TSpec, TStructuring > _defStruct;
	
		// pointer to tower in the upper layers
	SkipElement< TObject, TModus, TSpec, TStructuring > * _up;

		// the actual key
	typename Key< TObject >::Type _key;

		// saved key/value pair
	TObject * _obj;
		

	typename Cargo< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type _cargo;

	
	friend inline
	typename Cargo< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type * 
	cargo( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &me._cargo;
	}

	template< typename TCargo > friend inline
	void
	_setCargo( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
					TCargo & cargo )
	{
		SEQAN_CHECKPOINT
		me._cargo = cargo;
	}
	
public:

			// standard constructor
	SkipBaseElement(void)
		: _up(NULL)
		, _obj(NULL)
	{
		_initCargo( this, _cargo );
	}

	template< typename TKey >
	SkipBaseElement( TObject * obj, TKey key )
		: _up(NULL)
		, _key( key )
		, _obj(obj)
	{
		_initCargo( this, _cargo );
	}
	
	~SkipBaseElement(void)
	{
	}
  
};	// struct SkipBaseElement

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
inline typename Reference< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type
value( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me)
{
SEQAN_CHECKPOINT
	return *me._obj;
}

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
inline typename Reference< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type
value( SkipBaseElement< TObject, TModus, TSpec, TStructuring > const & me)
{
SEQAN_CHECKPOINT
	return *me._obj;
}

//____________________________________________________________________________

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
inline typename GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >::Type
getValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
{
SEQAN_CHECKPOINT
	return *me._obj;
}

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
inline typename GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >::Type
getValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > const & me )
{
SEQAN_CHECKPOINT
	return *me._obj;
}

//____________________________________________________________________________

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue2 >
inline void
assignValue(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
			TValue2 & val)
{
SEQAN_CHECKPOINT
	*me._obj = val;
}

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue2 >
inline void
assignValue(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
			TValue2 const & val)
{
SEQAN_CHECKPOINT
	*me._obj = val;
}

//____________________________________________________________________________

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue2 >
inline void
moveValue(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
		  TValue2 & val)
{
SEQAN_CHECKPOINT
	move(*me._obj, val);
}

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue2 >
inline void
moveValue(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
		  TValue2 const & val)
{
SEQAN_CHECKPOINT
	move(*me._obj, val);
}

//____________________________________________________________________________

/*
.Internal.getObject:
..summary:Get the saved object which stores key-value information.
..cat:SkipList
..signature:getObject(element)
..param.element:The desired element.
...type:Class.SkipBaseElement
*/


template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
TObject * 
getObject( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me )
{
	SEQAN_CHECKPOINT
	return me->_obj;
}

//____________________________________________________________________________

/*
.Internal._setObject:
..summary:Set the saved object of a SkipBaseElement which stores key-value information.
..remarks: For Internal use only. To insert an object into a dynamic skip list, use insert.
..cat:SkipList
..signature:getObject(element, object)
..param.element:The desired element.
...type:Class.SkipBaseElement
..param.element:The desired element.
...type:Class.SkipBaseElement
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
void 
_setObject( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
			TObject * obj )
{
	SEQAN_CHECKPOINT
	me._obj = obj;
}


template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
void 
_setObject( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * data )
{
	SEQAN_CHECKPOINT
	me._obj = getObject( data );
}

//____________________________________________________________________________

/*
.Internal._setUp:
..summary:Set the pointer to the lowest SkipElement of the related tower.
..cat:SkipList
..signature:_setUp(element, up)
..param.element:The desired object.
...type:Class.SkipBaseElement
..param.up:The SkipElement of the tower.
...type:Class.SkipElement
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
void 
_setUp( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me, 
		SkipElement< TObject, TModus, TSpec, TStructuring > & up )
{
	SEQAN_CHECKPOINT
	me._up = &up;
}

//____________________________________________________________________________

/*
.Internal._getUp:
..summary:_get a pointer to the lowest SkipElement of the related tower.
..cat:SkipList
..signature:_getUp(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipElement.$SkipElement*$@ pointing to the lowest SkipElement in the related tower. NULL if no tower is related to $element$.
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
SkipElement< TObject, TModus, TSpec, TStructuring > & 
_getUp( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
{
	SEQAN_CHECKPOINT
	return *me._up;
}

//____________________________________________________________________________

/*
.Internal._getSucc:
..summary:Get a pointer to the succeeding SkipBaseElement.
..cat:SkipList
..signature:_getSucc(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to $element1$ succeeding $element$. NULL if no such $element1$ exists.
..remarks:If the containing SkipList is Deferred and not completely sorted, it is unlikely that the ordering of the keys is already established. 
..So $key( _getSucc( element ) ) >= key( element )$ does not hold for sure.
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
SkipBaseElement< TObject, TModus, TSpec, TStructuring > *
_getSucc( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
{
	SEQAN_CHECKPOINT
	return ( &me + 1 );
}

template< typename TObject, typename TSpec, typename TStructuring > inline
SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > *
_getSucc( SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > & me )
{
	SEQAN_CHECKPOINT
	return me._dynStruct._succ;
}
//____________________________________________________________________________

/*
.Internal._setSucc:
..summary:Set the pointer to the succeeding SkipBaseElement.
..cat:SkipList
..signature:_setSucc(element, succ)
..param.element:The element.
...type:Class.SkipBaseElement
..param.succ:The succeeding SkipBaseElement.
...type:Class.SkipBaseElement
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
void 
_setSucc( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me, 
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * succ )
{
	SEQAN_ASSERT2( false, "No dynamic related members in default mode" )
}

template< typename TObject, typename TSpec, typename TStructuring > inline
void 
_setSucc( SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > & me, 
			SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > * succ )
{
	SEQAN_CHECKPOINT
	me._dynStruct._succ = succ;
}

//____________________________________________________________________________

/*
.Function._getPred:
..summary:Get a pointer to the preceding SkipBaseElement.
..cat:SkipList
..signature:_getPred(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to $element1$ preceding $element$. NULL if no such $element1$ exists.
..remarks:If the containing SkipList is Deferred and not completely sorted, it is unlikely that the ordering of the keys is already established. 
..So $key( _getPred( element ) ) >= key( element )$ does not hold for sure.
*/
template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
_getPred( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me );

template< typename TObject, typename TSpec, typename TStructuring > inline
SkipBaseElement< TObject, Static, TSpec, TStructuring > * 
_getPred( SkipBaseElement< TObject, Static, TSpec, TStructuring > & me )
{
	SEQAN_CHECKPOINT
	return ( &me - 1 );
}

template< typename TObject, typename TSpec, typename TStructuring > inline
SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > * 
_getPred( SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > & me )
{
	SEQAN_CHECKPOINT
	return me._dynStruct._pred;
}

//____________________________________________________________________________

/*
.Internal._setPred:
..summary:Set the pointer to the preceding SkipBaseElement.
..cat:SkipList
..signature:_setPred(element, goPrevious)
..param.element:The element.
...type:Class.SkipBaseElement
..param.goPrevious:The succeeding SkipBaseElement.
...type:Class.SkipBaseElement
*/
template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
void 
_setPred(	SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me, 
			SkipBaseElement< TObject, TModus, TSpec, TStructuring > * goPrevious )
{
	SEQAN_ASSERT2( false, "No dynamic related members in default mode" )
}

template< typename TObject, typename TSpec, typename TStructuring > inline
void 
_setPred(	SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > & me, 
			SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > * goPrevious )
{
	SEQAN_CHECKPOINT
	me._dynStruct._pred = goPrevious;
}

//____________________________________________________________________________

/*
.Internal._getLeft:
..summary:Get a pointer to the next sorted SkipBaseElement on the left side.
..cat:SkipList
..signature:_getLeft(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to the next sorted element on the left side.
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
_getLeft( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
{
	SEQAN_CHECKPOINT
	return _getPred( me );
}

template< typename TObject, typename TModus, typename TSpec > inline 
SkipBaseElement< TObject, TModus, TSpec, Deferred > * 
_getLeft( SkipBaseElement< TObject, TModus, TSpec, Deferred > & me )
{
	SEQAN_CHECKPOINT
	return me._defStruct._left;
}

//____________________________________________________________________________

/*
.Internal._setLeft:
..summary:Get a pointer to the next sorted SkipBaseElement on the left side.
..cat:SkipList
..signature:_getLeft(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to the next sorted element on the left side.
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
void 
_setLeft( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
						SkipBaseElement< TObject, TModus, TSpec, TStructuring > * left )
{
	SEQAN_ASSERT2( false, "No deferred related members in default mode" )
}

template< typename TObject, typename TModus, typename TSpec > inline 
void 
_setLeft( SkipBaseElement< TObject, TModus, TSpec, Deferred > & me,
						SkipBaseElement< TObject, TModus, TSpec, Deferred > * left )
{
	SEQAN_CHECKPOINT
	me._defStruct._left = left;
}

//____________________________________________________________________________

/*
.Internal._getRight:
..summary:Get a pointer to the next sorted SkipBaseElement on the right side.
..cat:SkipList
..signature:_getRight(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to the next sorted element on the right side.
*/


template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
_getRight( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
{
	SEQAN_CHECKPOINT
	return _getSucc( me );
}


template< typename TObject, typename TSpec > inline 
SkipBaseElement< TObject, Dynamic, TSpec, Deferred > * 
_getRight( SkipBaseElement< TObject, Dynamic, TSpec, Deferred > & me )
{
	SEQAN_CHECKPOINT
	return me._defStruct._right;
}

template< typename TObject, typename TSpec > inline 
SkipBaseElement< TObject, Static, TSpec, Deferred > * 
_getRight( SkipBaseElement< TObject, Static, TSpec, Deferred > & me )
{
	SEQAN_CHECKPOINT
	return &me +( _getCount( me ) + 1 );
}

//____________________________________________________________________________

//right: durch _getRight ersetzen???

template< typename TObject, typename TSpec, typename TStructuring > inline
void
right( SkipBaseElement< TObject, Static, TSpec, TStructuring > *& me )
{
	SEQAN_CHECKPOINT
	me += ( _getCount( *me ) + 1 );
}

template< typename TObject, typename TSpec, typename TStructuring > inline
void
right( SkipBaseElement< TObject, Dynamic, TSpec, TStructuring > *& me )
{
	SEQAN_CHECKPOINT
	me = _getRight( *me );
}


//____________________________________________________________________________

/*
.Internal._setRight:
..summary:Get a pointer to the next sorted SkipBaseElement on the right side.
..cat:SkipList
..signature:_getRight(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:A @Class.SkipBaseElement.$SkipBaseElement*$@, pointing to the next sorted element on the right side.
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
void 
_setRight(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
		  SkipBaseElement< TObject, TModus, TSpec, TStructuring > * right )
{
}

template< typename TObject, typename TSpec > inline 
void 
_setRight(SkipBaseElement< TObject, Static, TSpec, Deferred > & me,
		  SkipBaseElement< TObject, Static, TSpec, Deferred > * right )
{
	// do nothing
}

template< typename TObject, typename TSpec > inline 
void 
_setRight(SkipBaseElement< TObject, Dynamic, TSpec, Deferred > & me,
		  SkipBaseElement< TObject, Dynamic, TSpec, Deferred > * right )
{
	SEQAN_CHECKPOINT
	me._defStruct._right = right;
}

//____________________________________________________________________________

/*
.Internal._getHeight:
..summary:Get the height of the related tower.
..cat:SkipList
..signature:_getHeight(element)
..param.element:The element.
...type:Class.SkipBaseElement
..returns:Height of the associated tower.
...type:
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
_getHeight(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
{
	SEQAN_CHECKPOINT
	return me._defStruct._height;	
}

//____________________________________________________________________________

/*
.Internal._setHeight:
..summary:Set the height of the related tower.
..cat:SkipList
..signature:_setHeight(element, height)
..param.element:The desired object.
...type:Class.SkipBaseElement
..param.height:The height.
...type:@Metafunction.Height.$Size< element-type >::Type$@
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
void 
_setHeight(	SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
			TSize height )
{
	SEQAN_CHECKPOINT
	me._defStruct._height = height;
}

//____________________________________________________________________________

// key

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > inline 
void 
assignKey( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
		TKey theKey)
{
SEQAN_CHECKPOINT
	me._key= theKey;
}
/* falsch einsortiert
template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > inline 
void 
assignKey( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me,
		TKey theKey)
{
SEQAN_CHECKPOINT
	me->_key= theKey;
}
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
typename Key< TObject >::Type 
key( SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
{
SEQAN_CHECKPOINT
	return me._key;
}

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
typename Key< TObject >::Type 
key(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
	TParam & p )
{
SEQAN_CHECKPOINT
	return key(me);
}


//____________________________________________________________________________

/*
.Internal._getCount:
..summary:_get the count-value of the SkipBaseElement.
..cat:SkipList
..signature:_getCount(element)
..param.element:The desired object.
...type:Class.SkipBaseElement
..returns:The count-value of element. Type is $TKey$.
..remarks:The count value is defined as the number of unsorted elements between $element$ and $_getRight(element)$
*/
template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
_getCount(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me )
{
SEQAN_CHECKPOINT
	return me._defStruct._count;
}

//____________________________________________________________________________

/*
.Internal._setCount:
..summary:Set the count value.
..cat:SkipList
..signature:_setCount(element, count)
..param.element:The desired object.
...type:Class.SkipBaseElement
..param.count:The new count-value of element.
*/

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
void 
_setCount(SkipBaseElement< TObject, TModus, TSpec, TStructuring > & me,
		  TSize count )
{
SEQAN_CHECKPOINT
	me._defStruct._count = count;
}

//////////////////////////////////////////////////////////////////////////////

	
} // namespace seqan
#endif //_SKIP_BASE_ELEMENT_H


