#ifndef SEQAN_HEADER_SKIP_LIST_ITERATOR_H
#define SEQAN_HEADER_SKIP_LIST_ITERATOR_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Iterator:
// 
//	bidirectional iterator for dynamic skip list (base layer is a list)
//	random access iterator for static skip list (base layer is an array)
//
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Metafunctions

// the standard iterator is a pointer to the entrys in the base layer
template <typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TIteratorSpec>
struct Iterator<SkipList<TObject, TModus, TSpec, TStructuring>, TIteratorSpec>
{
	typedef SkipBaseElement<TObject, TModus, TSpec, TStructuring> * Type;
};

//____________________________________________________________________________
// Value of iterator is value of SkipBaseElement

template <typename TObject, typename TModus, typename TSpec, typename TStructuring>
struct Value<SkipBaseElement<TObject, TModus, TSpec, TStructuring> *>:
	Value<TObject>
{
};

//____________________________________________________________________________

template <typename TObject, typename TModus, typename TSpec, typename TStructuring>
struct GetValue<SkipBaseElement<TObject, TModus, TSpec, TStructuring> *>:
	GetValue<TObject>
{
};

//____________________________________________________________________________

template <typename TObject, typename TModus, typename TSpec, typename TStructuring>
struct Reference<SkipBaseElement< TObject, TModus, TSpec, TStructuring> *>:
	Reference<TObject>
{
};

//////////////////////////////////////////////////////////////////////////////
// Functions

//____________________________________________________________________________

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
void
goNext(SkipBaseElement<TObject, TModus, TSpec, TStructuring> * & me )
{
SEQAN_CHECKPOINT
	me = _getSucc( *me );
}

//____________________________________________________________________________

template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
void
goPrevious(SkipBaseElement< TObject, TModus, TSpec, TStructuring > * & me )
{
SEQAN_CHECKPOINT
	me = _getPred( *me );
}

//____________________________________________________________________________

template < typename TObject, typename TModus, typename TSpec, typename TStructuring>
inline typename Value<SkipBaseElement<TObject, TModus, TSpec, TStructuring> *>::Type
value(SkipBaseElement<TObject, TModus, TSpec, TStructuring> * me)
{
SEQAN_CHECKPOINT
	return value(*me);
} 
template < typename TObject, typename TModus, typename TSpec, typename TStructuring>
inline typename Value<SkipBaseElement<TObject, TModus, TSpec, TStructuring> const *>::Type
value(SkipBaseElement<TObject, TModus, TSpec, TStructuring> const * me)
{
SEQAN_CHECKPOINT
	return value(*me);
} 

//____________________________________________________________________________

template < typename TObject, typename TModus, typename TSpec, typename TStructuring>
inline typename GetValue<SkipBaseElement<TObject, TModus, TSpec, TStructuring> *>::Type
getValue(SkipBaseElement<TObject, TModus, TSpec, TStructuring> * me)
{
SEQAN_CHECKPOINT
	return getValue(*me);
} 
template < typename TObject, typename TModus, typename TSpec, typename TStructuring>
inline typename GetValue<SkipBaseElement<TObject, TModus, TSpec, TStructuring> const *>::Type
getValue(SkipBaseElement<TObject, TModus, TSpec, TStructuring> const * me)
{
SEQAN_CHECKPOINT
	return getValue(*me);
} 

//____________________________________________________________________________

template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue2>
inline void
assignValue(SkipBaseElement<TObject, TModus, TSpec, TStructuring> * me,
			TValue2 & val)
{
SEQAN_CHECKPOINT
	return assignValue(*me, val);
} 
template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue2>
inline void
assignValue(SkipBaseElement<TObject, TModus, TSpec, TStructuring> * me,
			TValue2 const & val)
{
SEQAN_CHECKPOINT
	return assignValue(*me, val);
} 

//____________________________________________________________________________

template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue2>
inline void
moveValue(SkipBaseElement<TObject, TModus, TSpec, TStructuring> * me,
		  TValue2 & val)
{
SEQAN_CHECKPOINT
	return moveValue(*me, val);
} 
template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue2>
inline void
moveValue(SkipBaseElement<TObject, TModus, TSpec, TStructuring> * me,
		  TValue2 const & val)
{
SEQAN_CHECKPOINT
	return moveValue(*me, val);
} 

//____________________________________________________________________________
/*
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
inline void
valueConstruct( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * it )
{
SEQAN_CHECKPOINT
	new( it ) SkipBaseElement< TObject, TModus, TSpec, TStructuring >;
}

template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
inline void
valueConstruct( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * it,
				TObject * obj )
{
SEQAN_CHECKPOINT
	new( it ) SkipBaseElement< TObject, TModus, TSpec, TStructuring >( obj );
}
*/
template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey >
inline void
valueConstruct( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * it,
				TObject * obj,
				TKey key )
{
SEQAN_CHECKPOINT
	new( it ) SkipBaseElement< TObject, TModus, TSpec, TStructuring >( obj, key );
}

//____________________________________________________________________________
/*
template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
inline void
valueDestruct( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * it )
{
SEQAN_CHECKPOINT
	it->~SkipBaseElement< TObject, TModus, TSpec, TStructuring >();
}
*/

//____________________________________________________________________________
//key

template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TKey > inline 
void 
assignKey(SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me,
		  TKey theKey)
{
SEQAN_CHECKPOINT
	me->_key= theKey;
}


template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
typename Key< TObject >::Type 
key(SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me )
{
SEQAN_CHECKPOINT
	return me->_key;
}

//____________________________________________________________________________

//////////////////////////////////////////////////////////////////////////////
/*
template< typename TType, typename TIterType > inline
typename Key< typename Value< typename Iterator< TType >::Type >::Type >::Type
key( Iter< TType, TIterType > & it )
{
	return key( value( it ) );
}

template< typename TType, typename TIterType, typename TParam > inline
typename Key< typename Value< typename Iterator< TType >::Type >::Type >::Type
key( Iter< TType, TIterType > & it,
		TParam & param )
{
	return key( value( it ), param );
}
*/


//////////////////////////////////////////////////////////////////////////////
//Krimskrams fuer Full Iterator - braucht man vermutlich gar nicht
/* 
template <typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec>
inline typename Position<Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const>::Type 
position(Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & me)
{
SEQAN_CHECKPOINT
	typename Position<Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const>::Type pos = 0;
	Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > it = begin( container(me) );
	while( it != me )
	{
		goNext( it );
		++pos;
	}
	return pos;
}

template <typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TContainer2>
inline typename Position<Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const>::Type 
position(Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & me,
			TContainer2 const &)
{
SEQAN_CHECKPOINT
	return position( me );
}

template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
inline Iter< SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> >  
operator + (Iter< SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	Iter< SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > buffer(container(left), hostIterator(left) );
	while( right > 0 )
	{
		goNext( buffer );
		--right;
	}
	return buffer;
}
template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> >  
operator + (TIntegral left,
			Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & right)
{
SEQAN_CHECKPOINT
	Iter< SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > buffer(container(left), hostIterator(left) );
	while( right > 0 )
	{
		goNext( buffer );
		--right;
	}
	return buffer;
}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > &
operator += (Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > & left,
				TIntegral right)
{
SEQAN_CHECKPOINT
	while( right > 0 )
	{
		goNext( left );
		--right;
	}
	return left;
}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> >  
operator - (Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & left,
			TIntegral right)
{
SEQAN_CHECKPOINT
	Iter< SkipList< TObject, Dynamic, TSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > buffer(container(left), hostIterator(left) );
	while( right > 0 )
	{
		goPrevious( buffer );
		--right;
	}
	return buffer;
}


//____________________________________________________________________________

// ???
//template < typename TObject, typename Dynamic, typename TSpec, typename TStructuring, typename TIterator, typename TSpec>
//inline typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type  
//operator - (Iter<SkipList< TObject, Dynamic, TSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & left,
//			Iter<SkipList< TObject, Dynamic, TSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & right)
//{
//SEQAN_CHECKPOINT
//	return hostIterator(left) - hostIterator(right);
//}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > &
operator -= (Iter<SkipList< TObject, Dynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > & left,
				TIntegral right)
{
SEQAN_CHECKPOINT
	while( right > 0 )
	{
		goPrevious( left );
		--right;
	}
	return left;;
}
*/


} // namespace ...

//////////////////////////////////////////////////////////////////////////////


#endif //SEQAN_HEADER_SKIP_LIST_ITER_H
