#ifndef SEQAN_HEADER_TREECHAIN_H
#define SEQAN_HEADER_TREECHAIN_H

#include <algorithm>
#include <vector>

namespace seqan{

		// compute chain spec for G0 and G1 cost metric
	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring, typename TCostModell, typename TSpec >
	TScoreValue
	_compute_chain( TSource & source, 
					TDest & dest, 
					TCostModell cost, 
					Score< TScoreValue, TScoreType > const & score_, 
					Chainer tag,
					TStructuring structuring,
					TSpec spec )
	{
		SEQAN_CHECK( !empty( source ) )
			// define some basic types
		typedef typename Value< TSource >::Type FragType;
		typedef typename Weight< FragType >::Type WeightType;
		typedef typename Key< FragType >::Type PositionType;
		typedef typename Size< FragType >::Type SizeType;
		typedef typename TSpec::Type SpecType;

		SizeType dim = dimension( value( begin( source ) ) );
		
			// construct containers for classes
		String< _MetaFragment< FragType > > metas;
		reserve( metas, length( source ) + 2 );

		std::vector< _WrapperPoint< FragType > > points;
		points.reserve( 2 * ( length( source ) + 2 ) );
		
		//String< _WrapperPoint< FragType > > points;
		//reserve( points, 2 * ( length( source ) + 2 ) );

		String< _ChainPoint< FragType, SpecType > > end_points;
		reserve( end_points, length( source ) + 2 );
		
			// define origin and terminus fragment
		FragType startingFrag( dim );		
		FragType endFrag( dim );

			// build the environment (construct wrapper points, chain point, get coordinates of the terminus)
		_build_chain_environment( source, metas, points, end_points, startingFrag, endFrag, spec );

		typename Iterator< String< _MetaFragment< FragType > > >::Type lastMeta = end( metas );
		goPrevious( lastMeta );

			// set the score of the origin from -infinity to 0
		setScore( value( begin( metas ) ), 0 );

			// sort the wrapper points to apply the line sweep paradigma
		std::sort( points.begin(), points.end(), _ChainSorter< _WrapperPoint< FragType > >( ) );
	
			// build the RMT
		RangeTree< _ChainPoint< FragType, SpecType >, Static, RT< MaxTree< > >, TStructuring > tree( end_points, dim-1 );

			// algorithm main loop
			// traverse wrapper points
		typename std::vector< _WrapperPoint< FragType > >::iterator pointIt = points.begin();
		while( pointIt != points.end() )
		{
				// actual point is the beginning of a frag
				// => search for preceding fragment
			_MetaFragment< FragType > & meta = _meta(*pointIt );
			if( !_isEnd( *pointIt ) )
			{
				_ChainPoint< FragType, SpecType > buffer( meta, dim - 1, true );
				_ChainPoint< FragType, SpecType > * result = rangeMaxQuery( tree, buffer );

				SEQAN_CHECK( result != NULL )

				_setPred( meta, _meta( *result ) );
				setScore( meta, _maxPriority( value( lastMeta ), meta, *result, cost, score_, dim ) );
			}
			else{
					// point is the end of a frag
					// => activate it
				size_t offset = &meta - &_meta( *points.begin() );
				_ChainPoint< FragType, SpecType > & point = end_points[ offset ];

				setPriority( point, score( meta ) - _activatePriority( value( lastMeta ), point, cost, score_, dim ) );
				activate( tree, point );
			}
			goNext( pointIt );
		}
			// perform backtracking
		return _chain_trace( dest, metas );
	};
	
}

#endif

