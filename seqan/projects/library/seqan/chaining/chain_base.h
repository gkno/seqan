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
		
	struct _GenericChaining;
	typedef Tag<_GenericChaining> const GenericChaining;

	struct _RangetreeChaining;
	typedef Tag<_RangetreeChaining> const RangetreeChaining;


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


/*
	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring, typename TAlgorithm >
	TScoreValue
	compute_chain( TSource & source, TDest & dest, Score< TScoreValue, TScoreType > const & score_, TAlgorithm tag, TStructuring structuring );

	template< typename TSource, typename TDest, typename TCostModel, typename TScoreValue, typename TScoreType, typename TStructuring, typename TAlgorithm, typename TSpec >
	TScoreValue
	_compute_chain( TSource & source, TDest & dest, TCostModel type, Score< TScoreValue, TScoreType > const & score_, TAlgorithm tag, TStructuring structuring, TSpec spec );

*/
	template< typename TSource, typename TDest, typename TScoreValue, typename TScoreType, typename TStructuring >
	TScoreValue
	compute_chain( TSource & source, TDest & dest, Score< TScoreValue, TScoreType > const & score_, TStructuring structuring )
	{

//		SEQAN_CHECK2( scoreGapOpen( score_ ) == scoreGapExtend( score_ ), "Chaining only supports linear gap costs" )
//		SEQAN_CHECK2( scoreGapOpen( score_ ) >= 0 && scoreGapExtend( score_ ) >= 0 && scoreMismatch( score_ ) >= 0, "Scores should be positive" )
		switch( dimension( value( begin( source ) ) ) )
		{
			case 1: SEQAN_REPORT("One dimensional chaining not supported")
					return 0;
			case 2: if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					{
						return _compute_chain( source, dest, G_0_Cost(), score_, structuring, _ChainSpecType< Array< 1 > >() );
					}
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, structuring, _ChainSpecType< Array< 1 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, structuring,_ChainSpecType< Array< 2 > >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );
					break;
			case 3: if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, structuring, _ChainSpecType< Array< 2 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, structuring, _ChainSpecType< Array< 2 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, structuring,_ChainSpecType< Array< 3 > >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );
					break;
			case 4:	if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, structuring, _ChainSpecType< Array< 3 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, structuring, _ChainSpecType< Array< 3 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, structuring,_ChainSpecType< Array< 4 > >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );
					break;
			case 5:	if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, structuring, _ChainSpecType< Array< 4 > >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, structuring, _ChainSpecType< Array< 4 > >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, structuring,_ChainSpecType< Array< 5 > >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );
					break;
			default:if( scoreMismatch( score_ ) == 0 && scoreGapExtend( score_ ) == 0 )
					return _compute_chain( source, dest, G_0_Cost(), score_, structuring, _ChainSpecType< Default >() );
					else if( scoreMismatch( score_ ) > 0 && scoreGapExtend( score_ ) > 0 )
					{
						if( scoreMismatch( score_ ) == scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_1_Cost(), score_, structuring, _ChainSpecType< Default >() );
						else //if( scoreMismatch( score_ ) > 2 * scoreGapExtend( score_ ) )
							return _compute_chain( source, dest, G_SoP_Cost(), score_, structuring, _ChainSpecType< Default >() );
					}
					SEQAN_ASSERT2( false, "Gap/mismatch model not supported" );

		}
		return infimumValue< TScoreValue >();
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename TSource, typename TDest, typename TScoring>
inline typename Value<TScoring>::Type
chain(TSource & source,
	  TDest & dest,
	  TScoring const & scoring,
	  RangetreeChaining)
{
	typedef typename Iterator<TSource, Standard>::Type TSourceIterator;
	typedef typename Iterator<TDest, Standard>::Type TDestIterator;

	SEQAN_ASSERT(length(source))

	//transform coordinates to old style ("end is last item")
	unsigned int dim = dimension(source[0]);
	for (TSourceIterator it = begin(source, Standard()); it < end(source, Standard()); ++it)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			_setLeftPosition(*it, i, leftPosition(*it, i) + 1);
		}
	}

	//compute chain
	typename Value<TScoring>::Type ret_value = compute_chain(source, dest, scoring, SemiDeferred()); //Deferred(), Complete()

	//retransform coordinates to new style ("end is behind last item")
	for (TSourceIterator it = begin(source, Standard()); it < end(source, Standard()); ++it)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			_setLeftPosition(*it, i, leftPosition(*it, i) - 1);
		}
	}
	for (TDestIterator it = begin(dest, Standard()) + 1; it < end(dest, Standard()); ++it)
	{
		for (unsigned int i = 0; i < dim; ++i)
		{
			_setLeftPosition(*it, i, leftPosition(*it, i) - 1);
		}
	}

	//return score
	return ret_value;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Function.chain:
..summary:Computes the chain on a set of fragments.
..cat:Chaining
..signature:chain(source, dest, score [, algorithm])
..param.source:The set of fragments
..param.dest:A destination container.
..param.score:The penalties for gaps between fragments.
...remarks:Values should be positive integers.
..param.algorithm:A tag that identifies the algorithm which is used for chaining.
...default:$Default$
...value:$Default$: Compiler selects best algorithm.
...value:$GenericChaining$: A simple generic chaining algorithm.
...value:$RangetreeChaining$: An elaborated chaining algorithm for @Spec.Score Zero@, @Spec.Score Manhattan@, and @Score ChainSoP@ scoring schemes.
*/



//implementation for Default
template< typename TSource, typename TDest, typename TScoring>
inline typename Value<TScoring>::Type
chain(TSource & source, 
	  TDest & dest, 
	  TScoring const & scoring,
	  Default)
{
	return chain(source, dest, scoring, GenericChaining()); //default is GenericChaining
}
template< typename TSource, typename TDest, typename TValue>
inline TValue
chain(TSource & source, 
	  TDest & dest, 
	  Score<TValue, Zero> const & scoring,
	  Default)
{
	return chain(source, dest, scoring, RangetreeChaining());
}
template< typename TSource, typename TDest, typename TValue>
inline TValue
chain(TSource & source, 
	  TDest & dest, 
	  Score<TValue, Manhattan> const & scoring,
	  Default)
{
	return chain(source, dest, scoring, RangetreeChaining());
}
template< typename TSource, typename TDest, typename TValue>
inline TValue
chain(TSource & source, 
	  TDest & dest, 
	  Score<TValue, ChainSoP> const & scoring,
	  Default)
{
	return chain(source, dest, scoring, RangetreeChaining());
}


//chain(3) => chain(4)
template< typename TSource, typename TDest, typename TScoring>
inline typename Value<TScoring>::Type
chain(TSource & source, 
	  TDest & dest, 
	  TScoring const & scoring)
{
	return compute_chain(source, dest, scoring, Default());
}


///////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif	// SEQAN_HEADER_CHAIN_BASE_H
