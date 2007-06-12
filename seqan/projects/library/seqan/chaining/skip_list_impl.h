/*	Copyright (c) 2006 Hendrik Woehrle
*	All rights reserved.
*
*	Contains class SkipList< TObject, TModus, TSpec, TStructuring > implementation
*
*/

#ifndef SEQAN_HEADER_SKIP_LIST_IMPL_H
#define SEQAN_HEADER_SKIP_LIST_IMPL_H

namespace seqan
{

/**
.Class.SkipList:
..cat:SkipList
..summary:A SkipList is a randomized alternative for rooted trees. Offers fast searching, insertion and deletion operations
objects. Saved objects are sorted with respect to their "key"-attribute.
..signature:SkipList< TObject, [ TModus, TSpec, TStructuring] >
..param.TObject:Type of stored objects.
..param.TModus:Modus of operation of a SkipList. A SkipList can either be dynamic or static. 
Dynamic Skip Lists admit insertion and deletion operations of elements, but the construction time is higher 
compared to a static SkipList. Default is Dynamic.
..param.TSpec:Specialization of the SkipList.
..param.TStructuring:Parameter to specify whether the SkipList uses Deferred Data Structuring or not.
..remarks:The object given to the SkipList should offer the following functions:
..remarks:$key( obj )$: returns the key of the object.
..remarks:$assignKey( obj, k )$: set the key of the object to k.
..remarks:In contrast to STL-like containers, the objects are not cloned by the SkipList. It only supports searching operations on a set of objects. This set must be handled by the user.
*/


///////////////////////////////////////////////////////////////////////////////
//
//								 class SkipList
//
///////////////////////////////////////////////////////////////////////////////


template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
class SkipList
{	
public:
	typedef SkipElement< TObject, TModus, TSpec, TStructuring > TSkipElement;
	typedef SkipBaseElement< TObject, TModus, TSpec, TStructuring > TSkipBaseElement;


		// container for elements in the lowest layer
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * _baseStore;

		// container for elements on the left side
		// search operations start from this bording elements
	SkipElement< TObject, TModus, TSpec, TStructuring > * _leftSideStore;

		// number of elements in the list
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type _numOfElements;

		// pool allocators
	Allocator< SinglePool2< TSkipBaseElement > > _baseAlloc;
	Allocator< ChunkPool2< TSkipElement > > _elementAlloc;
		
		// search path
	_SearchPath< TObject, TModus, TSpec, TStructuring > _sp;
		
		// border element for the right side
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * _rightBorder;

		// border objects
	TObject l_border_obj;
	TObject r_border_obj;

		// 
	bool _initialState;
	
//*************************************** internal help functions: ************************************************


private:

	SkipList & operator=( const SkipList & old )
	{}

	SkipList( const SkipList & old )
		: _numOfElements( length( old ) ) // bording elements are not included in Entries array => numOfElements + 2 elements in base layer	
		, _baseAlloc( _numOfElements + 2 )	// guessing needed space: a basic layer must be possible
		, _elementAlloc( _numOfElements )
		, _sp( _getMaximalSLTowerHeight( _numOfElements ) )
		, _initialState( true )
	{
			// construct bording elements
		assignKey( l_border_obj, infimumValue< typename Key< TObject >::Type >() );
		assignKey( r_border_obj, supremumValue< typename Key< TObject >::Type >() );

		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base_right;
		
		_initBases( *this, _baseStore, base_right, begin( old ), end( old ), _numOfElements ); 

		_initSL( *this, _baseStore, base_right, _numOfElements );
		_rightBorder = base_right;
	}

public:
				
	SkipList(void)
		: _numOfElements( 0 ) // bording elements are not included in Entries array => numOfElements + 2 elements in base layer	
		, _baseAlloc( 100  )	// guessing needed space: a basic layer must be possible
		, _elementAlloc( 100 )
		, _sp( 10 )
		, _initialState( true )
	{
			// construct bording elements
		assignKey( l_border_obj, infimumValue< typename Key< TObject >::Type >() );
		assignKey( r_border_obj, supremumValue< typename Key< TObject >::Type >() );
			
	}

	template< typename TContainer >
	SkipList( TContainer & cont )
		: _numOfElements( length( cont ) ) // bording elements are not included in Entries array => numOfElements + 2 elements in base layer	
		, _baseAlloc( _numOfElements + 2 )	// guessing needed space: a basic layer must be possible
		, _elementAlloc( _numOfElements )
		, _sp( _getMaximalSLTowerHeight( _numOfElements ) )
		, _initialState( true )
	{
			// construct bording elements
		assignKey( l_border_obj, infimumValue< typename Key< TObject >::Type >() );
		assignKey( r_border_obj, supremumValue< typename Key< TObject >::Type >() );

		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base_right;
		
		_initBases( *this, _baseStore, base_right, begin( cont ), end( cont ), _numOfElements ); 

		_initSL( *this, _baseStore, base_right, _numOfElements );
		_rightBorder = base_right;
	}

		// Constructor
		// Needs a range, defined by to iterators

		// TODO: für allgmeinen iterator verfügbar machen
	template< typename TIterator >
	SkipList( TIterator beg, TIterator end )
		: _numOfElements( end - beg ) // bording elements are not included in Entries array => numOfElements + 2 elements in base layer	
		, _baseAlloc( _numOfElements + 2 )	// guessing needed space: a basic layer must be possible
		, _elementAlloc( _numOfElements )
		, _sp( _getMaximalSLTowerHeight( _numOfElements ) )
		, _initialState( true )
	{			
		assignKey( l_border_obj, infimumValue< typename Key< TObject >::Type >() );
		assignKey( r_border_obj, supremumValue< typename Key< TObject >::Type >() );

			// construct bording elements
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base_right;

		_initBases( *this, _baseStore, base_right, beg, end, _numOfElements ); 

		_initSL( *this, _baseStore, base_right, _numOfElements );
		_rightBorder = base_right;
	}
	

	~SkipList(void)
	{
		_clearSearchPath( this->_sp, _getMaximalSLTowerHeight( _numOfElements ) );
	}

	// Debug print methods
private:
/*
	template< typename TSize1, typename TSize2 > friend
	void 
	printLayer(	SkipList< TObject, TModus, TSpec, TStructuring > & me,
				TSize1 layer,
				TSize2 column )
	{
		if( layer == 0 )
		{
			for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type j = 0; j < 11; ++j )
			{
				std::cout<< "______";
			}
			std::cout<<std::endl;
			typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type temp = begin( me );
			goPrevious( temp );
			goFurther( temp, column );
			typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type tempEnd = temp;
			goFurther( tempEnd, 11 );
			while( temp != end( me ) && temp != tempEnd )
			{
				std::cout.width(7);
				if( key( temp ) == infimumValue< typename Key< TObject >::Type >( ) )
					std::cout << std::left << "L";
				else
					std::cout << std::left << key( temp );
				goNext( temp );
			}
			//std::cout<<std::endl;
			//printCounts( me );
			//std::cout<<std::endl;
			//printSorting( me );
			//std::cout<<std::endl;
			//printHeights( me );
			//std::cout<<std::endl;
			//printValues( me );
		}
		else if( layer <= _getCurrentLayer( me ) && !_getInitialState( me ) )
		{
			typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type temp = begin( me );
			goPrevious( temp );
			goFurther( temp, column );
			typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > >::Type tempEnd = temp;
			goFurther( tempEnd, 11 );
			while( temp != end( me ) && temp != tempEnd )
			{
				if( _getHeight( *hostIterator( temp ) ) >= layer )
				{
					if( _getRight( *( &_getUp( *hostIterator( temp ) ) + layer - 1) ) )
					{
						std::stringstream s;
						if( key( temp ) == infimumValue< typename Key< TObject >::Type >( ) )
							s << std::left << "L";
						else
							s << key( temp );
						s << ">";
						std::cout.width(7);
						std::cout << std::left << s.str();
					}
					else
						dump( *( &_getUp( *hostIterator( temp ) ) + layer - 1 ) );
				}
				else 
				{
					std::cout.width(7);
					std::cout << " ";
				}
				goNext( temp );
			} 
		}
		std::cout<<std::endl;
	}

	friend
	void 
	dump( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		int column = 0;

		while( column < length( me ) )
		{
			typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type j = _getMaximalSLTowerHeight( me._numOfElements );
			while( j > 0 ){
				printLayer( me, --j, column );
			}
			std::cout<< std::endl;
			column += 14;
		}
	}

	friend
	void 
	printCounts( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * temp = &me._baseStore[0];
		for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < me._numOfElements + 1; ++i ){
			std::cout.width(5);
			std::cout << std::left << _getCount( *temp );
			goNext( temp );
		}
	}

	friend
	void 
	printSorting( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SkipBaseElement< TObject, TModus, TSpec, TStructuring >* temp = &me._baseStore[0];
		for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < me._numOfElements + 1; ++i ){
			std::cout.width(5);
			if( _getRight( *temp ) != NULL )
				std::cout << std::left << "x";
			else
				std::cout << std::left << "o";
			goNext( temp );
		}
	}

	friend
	void
	printHeights( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SkipBaseElement< TObject, TModus, TSpec, TStructuring >* temp = &me._baseStore[0];
		for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < me._numOfElements + 1; ++i ){
			std::cout.width(5);
			std::cout << std::left << _getHeight( *temp );
			goNext( temp );
		}
	}

	friend
	void
	printValues( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		//SkipBaseElement< TObject, TModus, TSpec, TStructuring >* temp = &me._baseStore[0];
		//for( typename seqan::Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type i = 0; i < me._numOfElements + 1; ++i ){
		//	std::cout.width(5);
		//	std::cout << std::left << getValue( getObject( *temp ) );
		//	temp = _getSucc( *temp );
		//}
	}
*/	

}; // class SkipList



/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////



	
		// no append operations on skip lists
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	void
	append(SkipList< TObject, TModus, TSpec, TStructuring > & me,
		TParam & param)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam1, typename TParam2 > inline 
	void
	append( SkipList< TObject, TModus, TSpec, TStructuring > & me ,
			TParam1 & param1,
			TParam2 & param2 )
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	void
	appendValue( SkipList< TObject, TModus, TSpec, TStructuring > & me ,
					TParam & param)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam1, typename TParam2 > inline 
	void
	appendValue( SkipList< TObject, TModus, TSpec, TStructuring > & me ,
				TParam1 & param1,
				TParam2 & param2 )
	{
		// do nothing
	}

		// no assign operations on skip lists
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	void
	assign( SkipList< TObject, TModus, TSpec, TStructuring > & me,
			TParam & param)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam1, typename TParam2 > inline 
	void
	assign( SkipList< TObject, TModus, TSpec, TStructuring > & me,
			TParam1 & param1,
			TParam2 & param2 )
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam > inline 
	void
	assignValue( SkipList< TObject, TModus, TSpec, TStructuring > & me ,
					TParam & param)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TParam1, typename TParam2 > inline 
	void
	assignValue( SkipList< TObject, TModus, TSpec, TStructuring > & me ,
				TParam1 & param1,
				TParam2 & param2 )
	{
		// do nothing
	}


	// get the beginning of the skiplist
		// i.e. the left bording element with score - infinity
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > , Standard >::Type 
	_begin_default( SkipList< TObject, TModus, TSpec, TStructuring > & me,
				   Standard)
	{
	SEQAN_CHECKPOINT
		return _getSucc( *me._baseStore );
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > , Standard >::Type 
	_begin_default( SkipList< TObject, TModus, TSpec, TStructuring > const & me,
				   Standard)
	{
	SEQAN_CHECKPOINT
		return _getSucc( *me._baseStore );
	}


		// capacity
		// static case: the size of the list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	capacity( SkipList< TObject, TModus, TSpec, TStructuring > & me)
	{
	SEQAN_CHECKPOINT
		return length(me);
	}

		// dynamic case: the list can hold an infinite number of objevts
	template< typename TObject, typename TSpec, typename TStructuring > inline 
	typename Size< SkipList< TObject, Dynamic, TSpec, TStructuring > >::Type
	capacity( SkipList< TObject, Dynamic, TSpec, TStructuring > & me)
	{
	SEQAN_CHECKPOINT
		return supremumValue< typename Size< SkipList< TObject, Dynamic, TSpec, TStructuring > >::Type >();
	}


	template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring >, Standard>::Type 
	_end_default( SkipList< TObject, TModus, TSpec, TStructuring > & me,
					Standard)
	{
	SEQAN_CHECKPOINT
		return _getRightBorder( me );
	}
	template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
	inline typename Iterator< SkipList< TObject, TModus, TSpec, TStructuring > const, Standard>::Type 
	_end_default( SkipList< TObject, TModus, TSpec, TStructuring > const & me,
					Standard)
	{
	SEQAN_CHECKPOINT
		return _getRightBorder( me );
	}

		// the length == size
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	length( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
	SEQAN_CHECKPOINT
		return me._numOfElements;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline
	void 
	_setLength( SkipList< TObject, TModus, TSpec, TStructuring > & me,
				TSize num_elems )
	{
		SEQAN_CHECKPOINT
		me._numOfElements = num_elems;
	}

		// no moving objects in skip lists
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	moveValue( SkipList< TObject, TModus, TSpec, TStructuring > & me, 
				TPos & pos)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	moveValue( SkipList< TObject, TModus, TSpec, TStructuring > const & me, 
				TPos & pos)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	moveValue( SkipList< TObject, TModus, TSpec, TStructuring > & me, 
				TPos const & pos)
	{
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	moveValue( SkipList< TObject, TModus, TSpec, TStructuring > const & me, 
				TPos const & pos)
	{
		// do nothing
	}


			// no replacements in skip lists
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	replace( SkipList< TObject, TModus, TSpec, TStructuring > & me, 
				TPos & pos)
	{
	SEQAN_CHECKPOINT
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	replace( SkipList< TObject, TModus, TSpec, TStructuring >  const & me, 
				TPos & pos)
	{
	SEQAN_CHECKPOINT
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos >
	void
	replace( SkipList< TObject, TModus, TSpec, TStructuring > & me, 
				TPos const & pos)
	{
	SEQAN_CHECKPOINT
		// do nothing
	}


		// the value at a given position
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TPos>
	inline typename Reference< SkipList< TObject, TModus, TSpec, TStructuring > >::Type
	value( SkipList< TObject, TModus, TSpec, TStructuring > & me, 
			TPos pos)
	{
	SEQAN_CHECKPOINT
		SkipBaseElement< TObject, TModus, TSpec, TStructuring > * base = me._baseStore;
		goNext( base );
		while( pos > 0 )
		{
			goNext( base );
			--pos;
		}
		return value( *base );
	}


	template< typename TObject, typename TSpec, typename TStructuring, typename TPos>
	inline typename Reference< SkipList< TObject, Static, TSpec, TStructuring > >::Type
	value( SkipList< TObject, Static, TSpec, TStructuring > & me, 
			TPos pos)
	{
	SEQAN_CHECKPOINT
		return value( me._baseStore + pos + 1 );
	} 

///////////////////////////////////////////////////////
//
//	member accessors (internal)
//
///////////////////////////////////////////////////////
	

		// get root node
		// i.e. the skip element in left bording tower in the current layer
		// current layer == 0 at the beginning
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	SkipElement< TObject, TModus, TSpec, TStructuring > * 
	_getRoot( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getLeftSideStore( me )+ _getCurrentLayer( me ) - 1;
	}

//____________________________________________________________________________

		// get the element allocator of a skiplist
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ChunkPool2< SkipElement<TObject, TModus, TSpec, TStructuring> > > & 
	_getElementAlloc( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT2( &me._elementAlloc != NULL, "Allocator for layer elements not initialized" )
		return me._elementAlloc;
	}

		// get the base element allocator of skiplist
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< SinglePool2< SkipBaseElement<TObject, TModus, TSpec, TStructuring> > > & 
	_getBaseAlloc( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT2( &me._baseAlloc != NULL, "Allocator for base elements not initialized" )
		return me._baseAlloc;
	}

//____________________________________________________________________________

		//get the right bording element
		//i.e. the lement with theKey== + infinity
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > * 
	_getRightBorder( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT2( me._rightBorder != NULL, "Right border not created" )
		return me._rightBorder;
	}

		// get the left side store of the skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipElement< TObject, TModus, TSpec, TStructuring > * 
	_getLeftSideStore( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT2( me._leftSideStore != NULL, "_leftSideStore not initialized" );
		return me._leftSideStore;
	}

		// get the base layer of the skip list
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipBaseElement< TObject, TModus, TSpec, TStructuring > *
	_getBaseStore( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		SEQAN_ASSERT2( me._baseStore != NULL, "_baseStore not initialized" );
		return me._baseStore;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipElement< TObject, TModus, TSpec, TStructuring > ** 
	_getSearchPath( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getSearchPath( me._sp );
	}

//____________________________________________________________________________

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, TSpec, TStructuring > >::Type 
	_getCurrentLayer( SkipList< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return _getHeight( *me._baseStore );
	}

		// renew current layer
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TSize > inline 
	void 
	_setCurrentLayer(	SkipList< TObject, TModus, TSpec, TStructuring > & me,
						TSize layer )
	{
		SEQAN_CHECKPOINT
		SEQAN_CHECK2( layer > _getCurrentLayer( me ), "Setting smaller value for _currentLayer" )
		_setHeight( *me._baseStore, layer );
	}
/////////////////////////////////////////////////////////////////////////////////////////

} // namespace seqan
#endif //_SKIP_LIST_STATIC_H

