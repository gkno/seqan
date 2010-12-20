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

#ifndef SEQAN_H_CHAIN_WRAPPER_POINT
#define SEQAN_H_CHAIN_WRAPPER_POINT

namespace seqan{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class WrapperPoint_
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
.Class.WrapperPoint_:
..summary:Data structure to represent a begin or end point of a fragment
..cat:Chaining
..signature:WrapperPoint_<TBorder>
..param.TBorder:Type of the class that represents the limits (multidimensional point) of an fragment.
..include:seqan/chaining.h
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			general functions
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

			// get the Key of the bording points of the related fragment depends on _end 
	template< typename TFragType > inline 
	typename Key< TFragType >::Type
	key( WrapperPoint_< TFragType > & me )
	{
		return me._key;
	}

	template< typename TFragType > inline 
	typename Key< TFragType >::Type
	key( const WrapperPoint_< TFragType > & me )
	{
		return me._key;
	}


		// get the related fragment
	
	template< typename TFragType > inline
	TFragType & 
	_getFrag( WrapperPoint_< TFragType > & me )
	{
		SEQAN_CHECK( me._frag != NULL )
		SEQAN_CHECK( me._meta != NULL )
		SEQAN_CHECK( &_getFrag( *me._meta ) == me._frag )
		return *me._frag;
	}

	template< typename TFragType > inline
	TFragType & 
	_getFrag( const WrapperPoint_< TFragType > & me )
	{
		SEQAN_CHECK( me._frag != NULL )
		SEQAN_CHECK( me._meta != NULL )
		SEQAN_CHECK( _getFrag( *me._meta ) == me._frag )
		return *me._frag;
	}

	template< typename TFragType > inline
	MetaFragment_< TFragType > & 
	_meta( WrapperPoint_< TFragType > & me )
	{
		SEQAN_CHECK( me._frag != NULL )
		SEQAN_CHECK( me._meta != NULL )
		return *me._meta;
	}

	template< typename TFragType > inline
	void
	_setMeta( WrapperPoint_< TFragType > & me,
				MetaFragment_< TFragType > & meta )
	{
		me._meta = &meta;
		me._frag = _getFrag( meta );
	}

	
	template< typename TFragType > inline
	bool 
	_isEnd( WrapperPoint_< TFragType > & me )
	{
		return me._end;
	}

	template< typename TFragType > inline 
	bool 
	_isEnd( const WrapperPoint_< TFragType > & me )
	{
		return me._end;
	}
	
	template< typename TFragType > inline 
	bool
	_isBegin( WrapperPoint_< TFragType > & me )
	{
		return !me._end;
	}

	template< typename TFragType > inline 
	bool
	_isBegin( const WrapperPoint_< TFragType > & me )
	{
		return !me._end;
	}

	template< typename TFragType > inline 
	typename Size< TFragType >::Type
	dimension( const WrapperPoint_< TFragType > & me )
	{
		return dimension( *me._frag );
	}

	
	template< typename TFragType >
	struct WrapperPoint_
	{
			// related triple, which stores the preceding frgament, chain value...
		TFragType * _frag;
			// the related key
		typename Key< TFragType >::Type _key;
			// meta information struct for that point
		MetaFragment_< TFragType > * _meta;
			// the point is either the end ( _end == true ) or the beginning of a fragment
		bool _end;
			
		
	#ifdef _SEQAN_CHAIN_DEBUG
		friend inline
		void 
		dump( WrapperPoint_< TFragType > & me )
		{
			std::cout << "[ ";
			typename Size< TFragType >::Type dim = 0;
			if( me._end )
				std::cout << rightPosition( *me._frag, dim );
			else
				std::cout << leftPosition( *me._frag, dim );
			++dim;
			while( dim != dimension( *me._frag ) )
			{
				if( me._end )
					std::cout << " , " << rightPosition( *me._frag, dim );
				else
					std::cout << " , " << leftPosition( *me._frag, dim );
				++dim;
			}
			std::cout << " ]" << std::endl;
		}

	#endif // _SEQAN_CHAIN_DEBUG

		WrapperPoint_( )
			: _frag( NULL )
			, _key( 0 )
			, _meta( NULL )
			, _end( false )
		{
		}


		WrapperPoint_( const WrapperPoint_ & old )
			: _frag( old._frag )
			, _key( old._key )
			, _meta( old._meta )
			, _end( old._end )
		{
		}


		WrapperPoint_( MetaFragment_< TFragType > & meta, 
						bool end )
			: _frag( &_getFrag( meta ) )
			, _key( end ? rightPosition( _getFrag( meta ), dimension( _getFrag( meta ) ) - 1 ) : leftPosition( _getFrag( meta ), dimension( _getFrag( meta ) ) - 1 ) )
			, _meta( & meta )
			, _end( end )
		{}

		template< typename TKey >
		WrapperPoint_( MetaFragment_< TFragType > & meta,
						TKey key,
						bool end )
			: _frag( &_getFrag( meta ) )
			, _key( key )
			, _meta( & meta )
			, _end( end )
		{}


		WrapperPoint_ & operator=( const WrapperPoint_ & old )
		{
			if ( this == &old ) 
				return *this;
			_key = old._key;
			_end = old._end;
			_frag = old._frag;
			_meta = old._meta;
			return *this;
		}

	};

}
#endif // SEQAN_H_...
