#ifndef SEQAN_HEADER_SKIP_LIST_TYPE_H
#define SEQAN_HEADER_SKIP_LIST_TYPE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Keys
//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Key:
..summary:Type of the key attribute of an object. 
..signature:Key<T>::Type
..param.T:Type for which the key type is determined.
..returns.param.Type:Key type of $T$.
..remarks.text:The the key type of an object is used for sorting and searching.
*/

template <typename T>
struct Key
{
	typedef T Type;
};

//////////////////////////////////////////////////////////////////////////////
// Pair: first argument is key

template < typename _T1, typename _T2, typename _TCompression>
struct Pair;

//____________________________________________________________________________

template <typename TKey, typename TVal, typename _TCompression>
struct Key<Pair<TKey, TVal, _TCompression> >
{
	typedef TKey Type;
};
template <typename TKey, typename TVal, typename _TCompression>
struct Key<Pair<TKey, TVal, _TCompression> const>
{
	typedef TKey Type;
};

//____________________________________________________________________________

/**
.Function.key:
..summary:Get the the key of the element.
..cat:SkipList
..signature:key(element)
..param.element:The desired object.
..returns:The the key of the element. Type is @Metafunction.Key.$Key< "element-type" >::Type$@.
*/

template <typename TKey, typename TVal, typename _TCompression> 
inline TKey &
key(Pair<TKey, TVal, _TCompression> & me)
{
	return me.i1;
}
template <typename TKey, typename TVal, typename _TCompression> 
inline TKey &
key(Pair<TKey, TVal, _TCompression> const & me)
{
	return me.i1;
}

//____________________________________________________________________________

template <typename TKey, typename TVal, typename _TCompression> 
inline TKey
getKey(Pair<TKey, TVal, _TCompression> & me)
{
	return getValueI1(me);
}
template <typename TKey, typename TVal, typename _TCompression> 
inline TKey
getKey(Pair<TKey, TVal, _TCompression> const & me)
{
	return getValueI1(me);
}

//____________________________________________________________________________

/**
.Function.assignKey:
..summary:Set the theKeyattribute of the element.
..cat:SkipList
..signature:assignKey(element, theKey)
..param.element:The desired object.
...type:Class.SkipElement
..param.key:The desired object.
...type:Metafunction.Key.$Key< "element-type" >::Type
*/

template <typename TKey, typename TVal, typename _TCompression, typename TKey2> 
inline void
assignKey(Pair<TKey, TVal, _TCompression> & me,
		  TKey2 & _key)
{
	assignValueI1(me, _key);
}
template <typename TKey, typename TVal, typename _TCompression, typename TKey2> 
inline void
assignKey(Pair<TKey, TVal, _TCompression> & me,
		  TKey2 const & _key)
{
	assignValueI1(me, _key);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

struct Dynamic;

struct Static;

struct Complete;

struct Deferred;

//////////////////////////////////////////////////////////////////////////////
// Forward declarations
//////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus = Dynamic, typename TSpec = Default, typename TStructuring = Complete >
struct SkipElement;

template< typename TObject, typename TModus = Dynamic, typename TSpec = Default, typename TStructuring = Complete >
struct SkipBaseElement;

template< typename TObject, typename TModus = Dynamic, typename TSpec = Default, typename TStructuring = Complete >
struct SkipList;


//////////////////////////////////////////////////////////////////////////////
//
// METAFUNCTIONS
//
//////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef size_t Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Size< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

///////////////////////////////////////////////////////////////////////////////////////////////
//		Position Type


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Position< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef SkipBaseElement< TObject, TModus, TSpec, TStructuring > * Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//		GetValue Type


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct GetValue< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename GetValue< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

///////////////////////////////////////////////////////////////////////////////////////////////
//		Value Type


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef TObject Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Value< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Value< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//		Key Type


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< TObject >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipList< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Key< SkipBaseElement< TObject, TModus, TSpec, TStructuring > const >
{
	typedef typename Key< SkipList< TObject, TModus, TSpec, TStructuring > >::Type const Type;
};


///////////////////////////////////////////////////////////////////////////////////////////////
//		Cargo Type


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >
{
	typedef Nothing Type;	// default: no cargo
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};

template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
struct Cargo< SkipBaseElement< TObject, TModus, TSpec, TStructuring > >
{
	typedef typename Cargo< SkipList< TObject, TModus, TSpec, TStructuring > >::Type Type;
};


template< typename TTag, typename TCargo > inline
void
_initCargo( TTag * tag, TCargo & _cargo )
{}


}

#endif // SEQAN_HEADER_SKIPLIST_TYPE_H
