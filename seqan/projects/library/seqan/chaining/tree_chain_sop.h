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

#ifndef SEQAN_HEADER_TREECHAIN_SOP_H
#define SEQAN_HEADER_TREECHAIN_SOP_H

#include <algorithm>
#include <vector>

namespace seqan
{
		// compute chain spec for G_SoP cost metric
	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring, typename TSpec >
	TScoreValue
	_computeChain( TSource & source, 
					TDest & dest, 
					GSumOfPairCost, 
					Score< TScoreValue, TScoreType > const & score_,
					TStructuring,
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
		SizeType facValue = fac( dim );
		
			// construct containers for classes
			// meta information
		String< MetaFragment_< FragType > > meta_inf;
		reserve( meta_inf, length( source ) + 2 );
//		typename Iterator< String< MetaFragment_< FragType > > >::Type metaIt = begin( meta_inf );
		
			// wrapper points for the line sweep paradigma
		std::vector< WrapperPoint_< FragType > > points;
		points.reserve( 2 * ( length( source ) + 2 ) );
		
		//String< WrapperPoint_< FragType > > points;
		//reserve( points, 2 * ( length( source ) + 2 ) );
		//typename Iterator< String< WrapperPoint_< FragType > > >::Type pointIt = begin( points );

			// points for the rmt with transformated coodinates
		String< String< ChainPoint_< FragType, SpecType > > > trans_points;
		resize( trans_points, fac( dim ) );
		typename Iterator< String< String< ChainPoint_< FragType, SpecType > > > >::Type transIt = begin( trans_points );
	
			// define origin and terminus fragment
		FragType startingFrag( dim );		
		FragType endFrag( dim );

			// get permutations
		String< SizeType > permutation;
		reserve( permutation, dim );
		_initPerm( permutation, dim );
		
		String< SizeType > bufferPermutation;
		reserve( bufferPermutation, dim );
		_initPerm( bufferPermutation, dim );

			// build the environment (construct wrapper points, chain point, get coordinates of the terminus)
		_buildChainEnvironment( source, meta_inf, points, trans_points, startingFrag, endFrag, permutation, facValue, spec );
		typename Iterator< String< MetaFragment_< FragType > > >::Type lastMeta = end( meta_inf );
		goPrevious( lastMeta );
	
			// sort the wrapper points to apply the line sweep paradigma
		std::sort( points.begin(), points.end(), ChainSorter_< WrapperPoint_< FragType > >( ) );

			// build the RMT's for all points
		String< RangeTree< ChainPoint_< FragType, SpecType >, SkipListStatic, RT< MaxTree< > >, TStructuring > * > trees;
		resize( trees, facValue );
		_resetPerm( permutation );
		_buildChainTrees( trees, trans_points, dim, facValue );
			
			// some buffer values
		WeightType * weights;
		allocate( weights, weights, dim );
		MetaFragment_< FragType > & firstMeta = value( begin( meta_inf ) );
		setScore( firstMeta, 0 );

			// algorithm main loop
			// traverse wrapper points
		typename std::vector< WrapperPoint_< FragType > >::iterator pointIt = points.begin();
		while( pointIt != points.end() )
		{
			MetaFragment_< FragType > & meta = _meta( *pointIt );
			typename Iterator< String< RangeTree< ChainPoint_< FragType, SpecType >, SkipListStatic, RT< MaxTree< > >, TStructuring > * >, Rooted >::Type treeIt = begin( trees );
	
				// actual point is the beginning of a frag
				// => search for preceding fragment
			if( !_isEnd( *pointIt ) )
			{
				
				typename Weight< FragType >::Type maxWeight = minValue< typename Weight< FragType >::Type >();
				typename Weight< FragType >::Type tempWeight = 0;
				MetaFragment_< FragType > * bestPoint = NULL; //&_meta( *points.begin() );
				ChainPoint_< FragType, SpecType > origin( meta, dim, true );
				ChainPoint_< FragType, SpecType > buffer = origin;

					// choose the best preceding point
				while( treeIt!= end( trees ) )
				{
						// transform coordinates for current permutation
					_chainTransformCoordsSearch( origin, buffer, permutation );
						
						// perform RMQ
					ChainPoint_< FragType, SpecType > * result = rangeMaxQuery( *value( treeIt ), buffer );
						// choose the best point of current hyper corner,
						// if one exists for the current permutation
					if( key( *result, 0 ) != minValue< typename Key< ChainPoint_< FragType, SpecType > >::Type >() )
					{
						tempWeight = score( _meta( *result ) );
						typename Weight< FragType >::Type cost_for_rest = _costGSoP( meta, _meta( *result ), score_, begin( permutation ), end( permutation ), dim );
						tempWeight -= cost_for_rest;
					}
					else
					{		// no active point in hyper corner
						std::next_permutation( begin( permutation ), end( permutation ) );
						goNext( treeIt );
						continue;
					}		
						// get point with maximal priority
					if( tempWeight > maxWeight )
					{
						maxWeight = tempWeight;
						bestPoint = &_meta( *result );
					}
					std::next_permutation( begin( permutation ), end( permutation ) );
					goNext( treeIt );
				}
				_setPred( meta, *bestPoint );
				setScore( meta, maxWeight + weight( meta ) );
			}
			else{
					// activate point with geometric costs
				size_t offset = &meta - &firstMeta;
				_getPermDifference( weights, dim, _getFrag( value( lastMeta ) ), _getFrag( meta ) );
                _getPerm( weights, begin( bufferPermutation ), dim );
				setPriority( meta, score( meta ) - _costGSoP( value( lastMeta ), meta, score_, begin( bufferPermutation ), end( bufferPermutation ), dim ) );
				transIt = begin( trans_points );

				_resetPerm( bufferPermutation );

				while( treeIt!= end( trees ) )
				{
						// activate in all trees
					ChainPoint_< FragType, SpecType > & point = value( transIt )[ offset ];

					typename Weight< FragType >::Type new_prio = score( meta ) - _costGSoP( value( lastMeta ), meta, score_, begin( bufferPermutation ), end( bufferPermutation ), dim );
					setPriority( point, new_prio );
					activate( *value( treeIt ), point );
					goNext( transIt );
					goNext( treeIt );
					std::next_permutation( begin( bufferPermutation ), end( bufferPermutation ) );
				}
			}
			goNext( pointIt );
		}
			// delete the rmt's
		_deleteChainTrees( trees, facValue );
		deallocate( weights, weights, dim );
			// backtracking
		return _chainTrace( dest, meta_inf );
	}

}

#endif
