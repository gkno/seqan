#ifndef SEQAN_RT_BASE_H
#define SEQAN_RT_BASE_H



namespace seqan{

/////////////////////////////////////////////////////////////////////////////////////////
//
//	declarations
//
/////////////////////////////////////////////////////////////////////////////////////////
		
		// standard tag struct for range tree
	template< typename TSpec = Default >
	struct RT
	{};

		// tag structs for grade of deferredness of the range tree
	struct SemiDeferred
	{};

	struct Complete
	{};

	struct Deferred
	{};

		// main classes
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct SkipList;
	
	template< typename TObject, typename TModus, typename TSpec = Default, typename TStructuring = Complete >
	class RangeTree;


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >
	{
		typedef _Empty Type;
	};

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Cargo< SkipBaseElement< TObject, TModus, RT< TSpec >, TStructuring > >
	{
		typedef typename Cargo< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type Type;
	};

/////////////////////////////////////////////////////////////////////////////////////////
//
//	dependent members
//
/////////////////////////////////////////////////////////////////////////////////////////

		// the memory allocators of the range tree
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct
	_RTreeAllocators;

	template< typename TObject, typename TSpec, typename TStructuring >
	struct
	_RTreeAllocators< TObject, Static, TSpec, TStructuring >
	{
		
		Allocator< ClassPool< SkipElement< TObject, Static, TSpec, TStructuring >, Limited > > _elementAlloc;
		Allocator< ClassPool< SkipList< TObject, Static, TSpec, TStructuring >, Unlimited > > _listAlloc;
		Allocator< SimpleAlloc<> > _baseAlloc;

		_RTreeAllocators()
			: _elementAlloc( NULL )
		{}

		template< typename TSize >
		_RTreeAllocators( TSize size1, TSize size2 )
			:_elementAlloc( size1 )
			,_listAlloc( size2 )
		{
			SEQAN_CHECKPOINT
		}

		~_RTreeAllocators( )
		{
			SEQAN_CHECKPOINT
		}
	};



/////////////////////////////////////////////////////////////////////////////////////////
//
//	utilities
//
/////////////////////////////////////////////////////////////////////////////////////////
		

	const size_t _rt_thresh = 16;

		
	template< typename TObject, typename TSpec, typename TStructuring > inline
	bool
	_checkAssocThresh( SkipBaseElement< TObject, Static, TSpec, TStructuring > * first,
						SkipBaseElement< TObject, Static, TSpec, TStructuring > * second )
	{
		SEQAN_CHECKPOINT
		return ( ( second - first ) > _rt_thresh );
	}

	template< typename TTarget, typename TSource > inline
	void
	_pushBack( TTarget & target, TSource const & source )
	{
		SEQAN_CHECKPOINT
#ifdef RTTIMETEST
		volatile TSource temp = source;
#else
		appendValue( target, source );
#endif
	}


			// setting the element to - infinity
	template< typename TObject, typename TSize > inline 
	void 
	_setMinInfty(	TObject & me,  
					TSize dim )
	{
		SEQAN_CHECKPOINT
		me = TObject( dim );
		typename Key< TObject>::Type infValue = infimumValue< typename Key< TObject>::Type >();
		for( typename Size< TObject >::Type i = 0; i < dimension( me ); ++i )
		{
			setKey( me, i, infValue );
		}
	}

		// setting the element to + infinity
	template< typename TObject, typename TSize > inline 
	void 
	_setMaxInfty(	TObject & me,  
					TSize dim )
	{
		
		SEQAN_CHECKPOINT
		me = TObject( dim );
		typename Key< TObject>::Type supValue = supremumValue< typename Key< TObject>::Type >();
		for( typename Size< TObject >::Type i = 0; i < dimension( me ); ++i )
		{
			setKey( me, i, supValue );
		}
	}

/////////////////////////////////////////////////////////////////////////////////////////
//
//	basic accessor functions
//
/////////////////////////////////////////////////////////////////////////////////////////
	

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ClassPool< SkipElement< TObject, TModus, RT< TSpec >, TStructuring >, Limited > > & 
	_getElementAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._elementAlloc;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< ClassPool< SkipList< TObject, TModus, RT< TSpec >, TStructuring >, Unlimited > > & 
	_getListAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._listAlloc;
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	Allocator< SimpleAlloc<> > &
	_getBaseAlloc( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._allocs._baseAlloc;
	}

	
		// accessor f�r grenzobjekte
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline 
	TObject * 
	_getLBorderObj( RangeTree< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &me._LBorderObj;
	}


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	TObject * 
	_getRBorderObj( RangeTree< TObject, TModus, TSpec, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return &me._RBorderObj;
	}

		// skip list of the main tree
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	SkipList< TObject, TModus, RT< TSpec >, TStructuring > *
	_getList( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & me )
	{
		SEQAN_CHECKPOINT
		return me._list;
	}



	
/////////////////////////////////////////////////////////////////////////////////////////
//
//	algorithms
//
/////////////////////////////////////////////////////////////////////////////////////////
	
	
/**
.Function.rangeQuery:
..summary:Get the object with maximum priority in the RMT in a given intervall
..cat:Range Tree
..signature:rangeMaxQuery(tree, lower_border, upper_border, dest)
..param.tree:A Range Tree.
...type:RangeTree
..param.lower_border:The object that stores the lower borders for all dimensions, i.e. $key( lower_border ) <= key( point in range )$
...type:TObject.
..param.lower_border:The object that stores the upper borders for all dimensions, i.e. $key( point in range ) <= key( upper_border )$
...type:TObject.
..param.dest:A container to save the objects.
..remarks:The size of $dest$ should be sufficient.
*/

		// perform a range query
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TResultSet >
	void
	rangeQuery( RangeTree< TObject, TModus, TSpec, TStructuring > & me, 
				TObject & first, 
				TObject & second,
				TResultSet & results )
	{
		SEQAN_CHECKPOINT
		if( dimension( me ) > 1 )
			_fingerSearch( _getList( me ), &first, &second, dimension( me ) - 1, results );
		else
			_bottomSearch( _getList( me ), &first, &second, results );
	}

	template< typename TObject, typename TSize > inline
	bool 
	_testBruteForce( TObject & elem,
					 TObject & first,
					 TObject & second,
					 TSize dim )
	{
		SEQAN_CHECKPOINT
		bool in_range = true;
		typename Size< TObject >::Type _dim = 0;
		while( in_range && _dim <= dim )
		{
			in_range = ( ( key( first, _dim ) <= key( elem, _dim ) ) && ( key( elem, _dim ) <= key( second, _dim ) ) );
			++_dim;
		}
		return in_range;
	}

		// test if an element is in range
		// from dim to dim - 1 to 0
	template< typename TObject, typename TSize > inline
	bool 
	_testRange(	TObject & elem,
				TObject & first,
				TObject & second,
				TSize dim )
	{
		SEQAN_CHECKPOINT
		bool in_range = true;
		typename Key< TObject >::Type theKey;
		while ( in_range && dim > 0 )
		{
			theKey = key( elem, dim );
			in_range = ( ( key( first, dim ) <= theKey ) && ( theKey <= key( second, dim ) ) );
			--dim;
		}
		theKey = key( elem, dim );
		in_range &= ( ( key( first, dim ) <= theKey ) && ( theKey <= key( second, dim ) ) );
		return in_range;
	}

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	typename Size< SkipList< TObject, TModus, RT< TSpec >, TStructuring > >::Type
	_getMaximalSLTowerHeight( RangeTree< TObject, TModus, RT< TSpec >, TStructuring > & rt )
	{
		SEQAN_CHECKPOINT
		return log2( length( rt ) ) + 2;
	}


	template< typename TObject, typename TSpec, typename TStructuring, typename TSize >
	void 
	printLayerScores(	SkipList< TObject, Static, RT< TSpec >, TStructuring > * list,
						TSize layer,
						TSize _dim )
	{}

}

#endif // SEQAN_RT_BASE
