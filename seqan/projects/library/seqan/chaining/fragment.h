#ifndef SEQAN_HEADER_FRAGMENT_H
#define SEQAN_HEADER_FRAGMENT_H


namespace seqan
{



/**
.Class.Fragment:
..summary:Data structure which represents an fragment
..cat:Chaining
..signature:Fragment<TBorder, TSpec>
..param.TBorder:Type of the class which represents the limits (multidimensional point) of an fragment.
..param.TSpec:Specialization tag
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class Fragment implementation
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	template< typename TBorder, typename TSpec > inline
	typename Size< Fragment< TBorder, TSpec > >::Type
	dimension( Fragment< TBorder, TSpec > & me )
	{
		return me._dim;
	}
	
	template< typename TBorder, typename TSpec, typename TSize > inline 
	TBorder 
	leftPosition( Fragment< TBorder, TSpec >  & me, 
					TSize dim )
	{
		return me._left[ dim ];
	}

	template< typename TBorder, typename TSpec, typename TSize > inline 
	TBorder 
	leftPosition( const Fragment< TBorder, TSpec >  & me, 
					TSize dim )
	{
		return me._left[ dim ];
	}

	template< typename TBorder, typename TSpec, typename TSize > inline 
	TBorder 
	rightPosition( Fragment< TBorder, TSpec >  & me, 
					TSize dim )
	{
		return me._right[dim];
	}

	template< typename TBorder, typename TSpec, typename TSize > inline 
	TBorder 
	rightPosition( const Fragment< TBorder, TSpec >  & me,
					TSize dim )
	{
SEQAN_ASSERT( me._left->size() == me._right->size() )
		return me._right[dim];
	}

	template< typename TBorder, typename TSpec > inline 
	typename Weight< Fragment< TBorder, TSpec > >::Type 
	weight( Fragment< TBorder, TSpec > & me )
	{
		return me._weight;
	}

	template< typename TBorder, typename TSpec > inline 
	typename Weight< Fragment< TBorder, TSpec > >::Type 
	weight( const Fragment< TBorder, TSpec > & me )
	{
		return me._weight;
	}


	template< typename TBorder, typename TSpec, typename TWeight > inline 
	typename Weight< Fragment< TBorder, TSpec > >::Type 
	setWeight( Fragment< TBorder, TSpec > & me,
				TWeight weight )
	{
		return me._weight = weight;;
	}

	template< typename TBorder, typename TSpec, typename TSize, typename TPosition > inline 
	void 
	_setLeftPosition( Fragment< TBorder, TSpec > & me,
						TSize dim,
						TPosition value )
	{
		SEQAN_CHECK2( dim >= 0 && dim < me._dim, "Dimension index out of bounds" );
		me._left[ dim ] = value;
	}

	template< typename TBorder, typename TSpec, typename TSize, typename TPosition > inline 
	void 
	_setRightPosition( Fragment< TBorder, TSpec > & me,
						TSize dim, 
						TPosition value )
	{
		SEQAN_CHECK2( dim >= 0 && dim < me._dim, "Dimension index out of bounds" );
		me._right[ dim ] = value;
	}


	template< typename TBorder, typename TSpec >
	struct Fragment
	{
		typename Key< Fragment< TBorder, TSpec > >::Type * _left;
		typename Key< Fragment< TBorder, TSpec > >::Type * _right;
		typename Size< Fragment< TBorder, TSpec > >::Type _dim;
		typename Weight< Fragment< TBorder, TSpec > >::Type _weight;

	#ifdef _SEQAN_CHAIN_DEBUG // some debugging variables to identify fragments while debugging
		char _id;

		friend inline 
		char 
		_getID( Fragment & me )
		{
			return me._id;
		}

		friend inline 
		char 
		_getID( const Fragment & me )
		{
			return me._id;
		}

	#endif // _SEQAN_CHAIN_DEBUG

	#ifdef _SEQAN_CHAIN_DEBUG
		static int _frag_counter;
	#endif // _SEQAN_CHAIN_DEBUG

		Fragment( )
			: _left( NULL )
			, _right( NULL )
			, _dim( 0 )
			, _weight( 0 )
		{}

		Fragment( typename Size< Fragment< TBorder, TSpec > >::Type dim )
			: _left( NULL )
			, _right( NULL )
			, _dim( dim )
			, _weight( 0 )
		{
			allocate( *this, _left, _dim );
			allocate( *this, _right, _dim );
		}

		Fragment( TBorder * left,
					TBorder * right,
					typename Size< Fragment< TBorder, TSpec > >::Type dim )
			: _left( NULL )
			, _right( NULL )
			, _dim( dim )
			, _weight( 0 )
		{
			allocate( *this, _left, _dim );
			allocate( *this, _right, _dim );
			_weight = 0;
			for( typename Size< Fragment< TBorder, TSpec > >::Type i = 0; i < dim; ++i )
			{
				_left[ i ] = left[ i ];
				_right[ i ] = right[ i ];
			}
		}

		Fragment( TBorder * left,
					TBorder * right,
					typename Size< Fragment< TBorder, TSpec > >::Type dim,
					typename Weight< Fragment< TBorder, TSpec > >::Type weight )
			: _left( NULL )
			, _right( NULL )
			, _dim( dim )
			, _weight( weight )
		{
			allocate( *this, _left, _dim );
			allocate( *this, _right, _dim );
			for( typename Size< Fragment< TBorder, TSpec > >::Type i = 0; i < dim; ++i )
			{
				_left[ i ] = left[ i ];
				_right[ i ] = right[ i ];
			}
		}

		
		Fragment & operator=( const Fragment & old )
		{
			if( this == &old) 
				return *this;
			if( _left )
				deallocate( *this, _left, _dim );
			if( _right )
				deallocate( *this, _right, _dim );
			_dim = old._dim;
			_weight =  old._weight;
			allocate( *this, _left, _dim );
			allocate( *this, _right, _dim );
			for( typename Size< Fragment< TBorder, TSpec > >::Type i = 0; i < _dim; ++i )
			{
				_left[ i ] = old._left[ i ];
				_right[ i ] = old._right[ i ];
			}			
		#ifdef _SEQAN_CHAIN_DEBUG
			_id = old._id;
		#endif
			return *this;
		}

		Fragment( const Fragment & old )
		{
			_dim = old._dim;
			_weight =  old._weight;
			allocate( *this, _left, _dim );
			allocate( *this, _right, _dim );
			for( typename Size< Fragment< TBorder, TSpec > >::Type i = 0; i < _dim; ++i )
			{
				_left[ i ] = old._left[ i ];
				_right[ i ] = old._right[ i ];
			}			
		#ifdef _SEQAN_CHAIN_DEBUG
			_id = old._id;
		#endif
		}

		~Fragment()
		{
			deallocate( *this, _left, _dim );
			deallocate( *this, _right, _dim );
			_dim = 0;
			_weight = 0;
		}

		friend inline
		void 
		dump( Fragment & me )
		{
			if( me._left != NULL )
			{
				std::cout << "[ ";
				typename Size< Fragment< TBorder, TSpec > >::Type dim = 0;
				std::cout << leftPosition( me, dim );
				++dim;
				while( dim != me._dim )
				{
					std::cout << " , " << leftPosition( me, dim );
					++dim;
				}
				std::cout << " ] * [ ";
				dim = 0;
				std::cout << rightPosition( me, dim );
				++dim;
				while( dim != me._dim )
				{
					std::cout << " , "<< rightPosition( me, dim );
					++dim;
				}
				std::cout << " ] "  << weight( me ) << std::endl;
			}
		}

		#ifdef _SEQAN_CHAIN_DEBUG
		friend inline
		void 
		dump( Fragment & me, std::ostream & os )
		{
			if( me._left != NULL )
			{
				typename Size< Fragment< TBorder, TSpec > >::Type dim = 0;
				os << "# " << weight( me )/10 << std::endl;
				while( dim != me._dim )
				{
					os << "[" << leftPosition( me, dim ) << "," << rightPosition( me, dim ) << "]";
					++dim;
				}
				os << std::endl;
			}
		}
		#endif // _SEQAN_CHAIN_DEBUG

	};


}; // namespace seqan

#endif // SEQAN_HEADER_FRAGMENT_H

