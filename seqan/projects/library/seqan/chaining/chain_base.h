#ifndef SEQAN_HEADER_CHAIN_BASE_H
#define SEQAN_HEADER_CHAIN_BASE_H

/*
 *  chain_base.h
 *  chaining
 *
 *  Created by Hendrik Woehrle
 *
 *	Basis declarations & definitions for chaining algorithms
 *
 */

namespace seqan
{
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			class forward declarations
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		// basic data structure for fragments
	template< typename TBorder, typename TSpec = Default >
	struct Fragment;

		// wrapper data structure for point in the RMT
	template< typename T, typename TSpec >
	struct _ChainPoint;

		// wrapper for points for the line sweep paradigma
	template< typename T >
	struct _WrapperPoint;

		// structure for metainformation about fragments in chains
	template< typename T >
	struct _MetaFragment;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			tag classes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
		
	struct Chainer
	{};
	//typedef Tag<Chainer_> const Chainer;

	struct G_0_Cost
	{};
	//typedef Tag<ZeroCost_> const G_0_Cost;

	struct G_1_Cost
	{};
	//typedef Tag<G_1_Cost_> const G_1_Cost;

	struct G_Inf_Cost
	{};
	//typedef Tag<G_Inf_Cost_> const G_Inf_Cost;

	struct G_SoP_Cost
	{};
	//typedef Tag<F_SoP_Cost_> const G_Inf_Cost;

	template< typename Spec >
	struct _ChainSpecType
	{
		typedef Spec Type;
	};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			special metatype declarations
//
///////////////////////////////////////////////////////////////////////////////////////////////////////
// sizes

	template< typename TBorder, typename TSpec >
	struct Size< Fragment< TBorder, TSpec > >
	{
		typedef size_t Type;
	};

	template< typename T, typename TSpec >
	struct Size< _ChainPoint< T, TSpec > >
	{
		typedef size_t Type;
	};

	template< typename T >
	struct Size< _WrapperPoint< T > >
	{
		typedef size_t Type;
	};

	template< typename T >
	struct Size< _MetaFragment< T > >
	{
		typedef size_t Type;
	};


//////////////////////////////////////////////////////////////////////////
// key types

	template< typename TBorder, typename TSpec >
	struct Key< Fragment< TBorder, TSpec > >
	{
		typedef TBorder Type;
	};

	template< typename T, typename TSpec >
	struct Key< _ChainPoint< T, TSpec > >
	{
		typedef typename Key< T >::Type Type;
	};

	template< typename T >
	struct Key< _WrapperPoint< T > >
	{
		typedef typename Key< T >::Type Type;
	};

	template< typename T >
	struct Key< _MetaFragment< T > >
	{
		typedef typename Key< T >::Type Type;
	};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			functions
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Function.compute_chain:
..summary:Computes the chain on a set of fragments.
..cat:Chaining
..signature:compute_chain(source, dest, score, algorithm, structuring )
..param.source:The set of fragments
..param.dest:A destination conainer.
..param.score:The penalties for gaps between fragments.
...remarks:Values should be positive integers.
..param.algorithm:A tag that identifies the algorithm which is used for chaining.
..param.structuring:Which type of RMT should be used.
*/

	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring, typename TAlgorithm >
	TScoreValue
	compute_chain( TSource & source, TDest & dest, Score< TScoreValue, TScoreType > const & score_, TAlgorithm tag, TStructuring structuring );

	template< typename TSource, typename TDest, typename TCostModel, typename TScoreValue, typename TScoreType, typename TStructuring, typename TAlgorithm, typename TSpec >
	TScoreValue
	_compute_chain( TSource & source, TDest & dest, TCostModel type, Score< TScoreValue, TScoreType > const & score_, TAlgorithm tag, TStructuring structuring, TSpec spec );


	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring >
	TScoreValue
	compute_chain( TSource & source, TDest & dest, Score< TScoreValue, TScoreType > const & score_, Chainer tag, TStructuring structuring )
	{

		SEQAN_CHECK2( scoreGapOpen( score_ ) == scoreGapExtend( score_ ), "Chaining only supports linear gap costs" )
		SEQAN_CHECK2( scoreGapOpen( score_ ) >= 0 && scoreGapExtend( score_ ) >= 0 && scoreMismatch( score_ ) >= 0, "Scores should be positive" )
		switch( dimension( value( begin( source ) ) ) )
		{
			case 1: SEQAN_REPORT("One dimensional chaining not supported")
					return 0;
			case 2: if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, tag, structuring, _ChainSpecType< Array< 1 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, tag, structuring, _ChainSpecType< Array< 1 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, tag, structuring,_ChainSpecType< Array< 2 > >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );
					break;
			case 3: if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, tag, structuring, _ChainSpecType< Array< 2 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, tag, structuring, _ChainSpecType< Array< 2 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, tag, structuring,_ChainSpecType< Array< 3 > >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );
					break;
			case 4:	if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, tag, structuring, _ChainSpecType< Array< 3 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, tag, structuring, _ChainSpecType< Array< 3 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, tag, structuring,_ChainSpecType< Array< 4 > >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );
					break;
			case 5:	if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, tag, structuring, _ChainSpecType< Array< 4 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, tag, structuring, _ChainSpecType< Array< 4 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, tag, structuring,_ChainSpecType< Array< 5 > >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );
					break;
			default:if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, tag, structuring, _ChainSpecType< Default >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, tag, structuring, _ChainSpecType< Default >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, tag, structuring, _ChainSpecType< Default >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );

		}
		return infimumValue< TScoreValue >();
	}
	
}

#endif	// SEQAN_HEADER_CHAIN_BASE_H
