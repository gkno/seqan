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

#ifndef SEQAN_HEADER_META_FRAGMENT
#define SEQAN_HEADER_META_FRAGMENT

namespace seqan{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Class MetaFragment_
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*DISABLED
.Class.MetaFragment_:
..summary:Basic data which associates fragments with prededing fragments and stores chain score informations
..cat:Chaining
..signature:MetaFragment_< TFragType >
..param.TFragType:Type of the fragment
..include:seqan/chaining.h
*/

		// get/set the weight of the related fragment
	template< typename TFragType > inline
	typename Weight< TFragType >::Type
	weight( MetaFragment_< TFragType > & me )
	{
		return weight( *me._frag );
	}

	template< typename TFragType, typename TWeight> inline
	void
	setWeight( MetaFragment_< TFragType > & me,
				TWeight weight )
	{
		setWeight( *me._frag, weight );
	}

		// get/set the score of the chain
	template< typename TFragType > inline
	typename Weight< TFragType >::Type 
	score( MetaFragment_< TFragType > & me )
	{
		return me._score;
	}

	template< typename TFragType, typename TWeight > inline
	void
	setScore( MetaFragment_< TFragType > & me,
						TWeight score )
	{
		me._score = score;
	}

		// get/set the priority
	template< typename TFragType > inline
	typename Weight< TFragType >::Type 
	priority( MetaFragment_< TFragType > & me )
	{
		return me._priority;
	}

	template< typename TFragType, typename TWeight > inline
	void
	setPriority( MetaFragment_< TFragType > & me,
					TWeight prio )
	{
		me._priority = prio;
	}

		// get the associated fragment
	template< typename TFragType > inline
	TFragType & 
	_getFrag( MetaFragment_< TFragType > & me )
	{
		return *me._frag;
	}

	template< typename TFragType > inline
	TFragType & 
	_getFrag( const MetaFragment_< TFragType > & me )
	{
		return *me._frag;
	}

		// get preceding fragment
	template< typename TFragType > inline
	MetaFragment_< TFragType > & 
	_getPred( MetaFragment_< TFragType > & me )
	{
		return *me._pred;
	}

	template< typename TFragType > inline 
	MetaFragment_< TFragType > & 
	_getPred( const MetaFragment_< TFragType > & me )
	{
		return *me._pred;
	}

		// set preceding fragment
	template< typename TFragType > inline
	void
	_setPred( MetaFragment_< TFragType > & me, 
				MetaFragment_< TFragType > & pred )
	{
		me._pred = &pred;
	}

	template< typename TFragType > inline
	void
	_setPred( const MetaFragment_< TFragType > & me, 
				MetaFragment_< TFragType > & pred )
	{
		me._pred = &pred;
	}

	template< typename TFragType > inline
	void 
	dump( MetaFragment_< TFragType > & me )
	{
		if( me._frag )
			dump( *me._frag );
		std::cout << me._priority << " " << me._score << std::endl;
	}


	template< typename TFragType >
	struct MetaFragment_
	{
		TFragType * _frag;
			// preceding element in a chain
		typename Weight< TFragType >::Type _priority;
		typename Weight< TFragType >::Type _score;
		MetaFragment_< TFragType > * _pred;

		MetaFragment_()		
			: _frag( NULL )
			, _priority( minValue< typename Weight< TFragType >::Type >() )
			, _score( minValue< typename Weight< TFragType >::Type >() )
			, _pred( NULL )
		{}

		MetaFragment_( TFragType & frag )
			: _frag( &frag )
			, _priority( minValue< typename Weight< TFragType >::Type >() )
			, _score( minValue< typename Weight< TFragType >::Type >() )
			, _pred( NULL )
		{}

		MetaFragment_( const MetaFragment_ & old )
			: _frag( old._frag)
			, _priority( old._priority )
			, _score( old._score )
			, _pred( old._pred )
		{}

		MetaFragment_ & operator=( const MetaFragment_ & old )
		{
			if ( this == &old ) 
				return *this;
			_frag = old._frag;
			_priority = old._priority;
			_score = old._score;
			_pred = old._pred;
			return *this;
		}

	};


}

#endif // SEQAN_HEADER_META_FRAGMENT
