/*
 *  modifier_base.h
 *  SeqAn
 *
 *  Created by David Weese on 24.01.07.
 *
 */

#ifndef SEQAN_HEADER_MODIFIER_BASE_H
#define SEQAN_HEADER_MODIFIER_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////


	template < typename THost, typename TSpec = void >
	struct Modifier;


//////////////////////////////////////////////////////////////////////////////


	template < typename THost, typename TSpec >
	struct Value< Modifier<THost, TSpec> > {
		typedef typename Value<THost>::Type Type;
	}

	template < typename THost, typename TSpec >
	struct Size< Modifier<THost, TSpec> > {
		typedef typename Size<THost>::Type Type;
	}

	template < typename THost, typename TSpec >
	struct Position< Modifier<THost, TSpec> > {
		typedef typename Position<THost>::Type Type;
	}

	template < typename THost, typename TSpec >
	struct Host< Modifier<THost, TSpec> > {
		typedef THost Type;
	}

	template < typename THost, typename TSpec >
	struct Iterator< Modifier<THost, TSpec> > {
		typedef typename Iterator<THost>::Type Type;
	}


//////////////////////////////////////////////////////////////////////////////


	template < typename THost, typename TSpec >
	inline typename Size< Modifier<THost, TSpec> >::Type 
	length(Modifier<THost, TSpec> const & me) {
		return length(host(me));
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< Modifier<THost, TSpec> const >::Type 
	begin(Modifier<THost, TSpec> const & me) {
		return begin(host(me));
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< Modifier<THost, TSpec> >::Type 
	begin(Modifier<THost, TSpec> & me) {
		return begin(host(me));
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< Modifier<THost, TSpec> const >::Type 
	end(Modifier<THost, TSpec> const & me) {
		return end(host(me));
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< Modifier<THost, TSpec> >::Type 
	end(Modifier<THost, TSpec> & me) {
		return end(host(me));
	}


}

#endif
