/*
 *  modifier_string.h
 *  genindex
 *
 *  Created by David Weese on 26.01.06.
 *
 */

#ifndef SEQAN_HEADER_MODIFIER_STRING_H
#define SEQAN_HEADER_MODIFIER_STRING_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////


	template < typename THost, typename TSpec = void >
	struct ModifiedString;

	template < typename THost, typename TSpec >
	struct Spec< ModifiedString<THost, TSpec> > {
		typedef TSpec Type;
	}


	template < typename THost, typename TSpec >
	struct Value< ModifiedString<THost, TSpec> >:
		Value<THost> {};

	template < typename THost, typename TSpec >
	struct Value< ModifiedString<THost, TSpec> const >:
		Value<THost const> {};


	template < typename THost, typename TSpec >
	struct GetValue< ModifiedString<THost, TSpec> >:
		GetValue<THost> {};

	template < typename THost, typename TSpec >
	struct GetValue< ModifiedString<THost, TSpec> const >:
		GetValue<THost const> {};


	template < typename THost, typename TSpec >
	struct Reference< ModifiedString<THost, TSpec> >:
		Reference<THost> {};

	template < typename THost, typename TSpec >
	struct Reference< ModifiedString<THost, TSpec> const >:
		Reference<THost const> {};


	template < typename THost, typename TSpec >
	struct Size< ModifiedString<THost, TSpec> >:
		Size<THost> {};


	template < typename THost, typename TSpec >
	struct Iterator< ModifiedString<THost, TSpec> >:
		Iterator<THost> {};

	template < typename THost, typename TSpec >
	struct Iterator< ModifiedString<THost, TSpec> const >:
		Iterator<THost const> {};


	template < typename THost, typename TSpec >
	struct Host< ModifiedString<THost, TSpec> > {
		typedef THost Type;
	}

	template < typename THost, typename TSpec >
	struct Host< ModifiedString<THost, TSpec> const > {
		typedef THost const Type;
	}


	//////////////////////////////////////////////////////////////////////////////
	// value
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TPos>
	inline typename Reference<ModifiedString<THost, TSpec> >::Type 
	value(ModifiedString<THost, TSpec> & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		return *(begin(me) + pos);
	}

	template <typename THost, typename TSpec, typename TPos>
	inline typename Reference<ModifiedString<THost, TSpec> const >::Type 
	value(ModifiedString<THost, TSpec> const & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		return *(begin(me) + pos);
	}

	//////////////////////////////////////////////////////////////////////////////
	// length
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TSpec >
	inline typename Size< ModifiedString<THost, TSpec> >::Type 
	length(ModifiedString<THost, TSpec> const & me) {
		return length(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// begin
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const >::Type 
	begin(ModifiedString<THost, TSpec> const & me) {
		return begin(host(me));
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
	begin(ModifiedString<THost, TSpec> & me) {
		return begin(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// end
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const >::Type 
	end(ModifiedString<THost, TSpec> const & me) {
		return end(host(me));
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
	end(ModifiedString<THost, TSpec> & me) {
		return end(host(me));
	}

}

#endif
