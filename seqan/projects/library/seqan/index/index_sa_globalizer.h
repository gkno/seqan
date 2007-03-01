/*
 *  index_sa_globalizer.h
 *  SeqAn
 *
 *  Created by David Weese on 26.01.06.
 *
 */

#ifndef SEQAN_HEADER_INDEX_SA_GLOBALIZER_H
#define SEQAN_HEADER_INDEX_SA_GLOBALIZER_H

namespace SEQAN_NAMESPACE_MAIN
{

	template <typename TLimitsString>
	struct ModGlobalizer {};


	template <typename THost, typename TLimitsString>
	class ModifiedIterator<THost, ModGlobalizer<TLimitsString> > {
	public:
		Holder<THost, Simple>	data_host;
		Holder<TLimitsString>	limits;

		mutable typename Value<TLimitsString>::Type globalPos;

		ModifiedIterator(THost &_host, TLimitsString &_limits):
			data_host(_host),
			limits(_limits) {}
	};

	template <typename THost, typename TLimitsString>
	struct Value< ModifiedIterator<THost, ModGlobalizer<TLimitsString> > >:
		Value<TLimitsString> {};

	template <typename THost, typename TLimitsString>
	struct GetValue< ModifiedIterator<THost, ModGlobalizer<TLimitsString> > >:
		Value<TLimitsString> {};

	template <typename THost, typename TLimitsString>
	struct Reference< ModifiedIterator<THost, ModGlobalizer<TLimitsString> > > {
		typedef typename Value<TLimitsString>::Type & Type;
	};



	template <typename THost, typename TLimitsString>
	class ModifiedString<THost, ModGlobalizer<TLimitsString> > {
	public:
		Holder<THost>			data_host;
		Holder<TLimitsString>	limits;

		ModifiedString(THost &_host, TLimitsString &_limits):
			data_host(_host),
			limits(_limits) {}

		template <typename TPos>
		inline typename Reference<ModifiedString>::Type 
		operator [] (TPos pos)
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<ModifiedString const>::Type 
		operator [] (TPos pos) const
		{
		SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};

	template <typename THost, typename TLimitsString>
	struct Value< ModifiedString<THost, ModGlobalizer<TLimitsString> > >:
		Value<TLimitsString> {};

	template <typename THost, typename TLimitsString>
	struct Value< ModifiedString<THost, ModGlobalizer<TLimitsString> > const> {
		typedef typename Value<TLimitsString>::Type const Type;
	};
/*
	template < typename THost, typename TSpec >
	struct GetValue< ModifiedString<THost, TSpec> >:
		GetValue<THost> {};

	template < typename THost, typename TSpec >
	struct GetValue< ModifiedString<THost, TSpec> const >:
		GetValue<THost const> {};
*/

	template <typename THost, typename TLimitsString>
	struct Reference< ModifiedString<THost, ModGlobalizer<TLimitsString> > > {
		typedef typename Value<TLimitsString>::Type & Type;
	};

	template <typename THost, typename TLimitsString>
	struct Reference< ModifiedString<THost, ModGlobalizer<TLimitsString> > const> {
		typedef typename Value<TLimitsString>::Type const & Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// begin
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TLimitsString >
	inline typename Iterator< ModifiedString<THost, ModGlobalizer<TLimitsString> > const >::Type 
	begin(ModifiedString<THost, ModGlobalizer<TLimitsString> > const & me) {
		return typename Iterator< ModifiedString<THost, ModGlobalizer<TLimitsString> > const >::Type
			(begin(host(me)), value(me.limits));
	}

	template < typename THost, typename TLimitsString >
	inline typename Iterator< ModifiedString<THost, ModGlobalizer<TLimitsString> > >::Type 
	begin(ModifiedString<THost, ModGlobalizer<TLimitsString> > & me) {
		return typename Iterator< ModifiedString<THost, ModGlobalizer<TLimitsString> > >::Type
			(begin(host(me)), value(me.limits));
	}

	//////////////////////////////////////////////////////////////////////////////
	// end
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TLimitsString >
	inline typename Iterator< ModifiedString<THost, ModGlobalizer<TLimitsString> > const >::Type 
	end(ModifiedString<THost, ModGlobalizer<TLimitsString> > const & me) {
		return typename Iterator< ModifiedString<THost, ModGlobalizer<TLimitsString> > const >::Type
			(end(host(me)), me.limits);
	}

	template < typename THost, typename TLimitsString >
	inline typename Iterator< ModifiedString<THost, ModGlobalizer<TLimitsString> > >::Type 
	end(ModifiedString<THost, ModGlobalizer<TLimitsString> > & me) {
		return typename Iterator< ModifiedString<THost, ModGlobalizer<TLimitsString> > >::Type
			(end(host(me)), me.limits);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator *
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TLimitsString>
	inline typename Reference<ModifiedIterator<THost, ModGlobalizer<TLimitsString> > >::Type 
	value(ModifiedIterator<THost, ModGlobalizer<TLimitsString> > & me)
	{
	SEQAN_CHECKPOINT
		return posGlobalize(value(host(me)), value(me.limits));
	}

	template <typename THost, typename TLimitsString>
	inline typename Reference<ModifiedIterator<THost, ModGlobalizer<TLimitsString> > const>::Type 
	value(ModifiedIterator<THost, ModGlobalizer<TLimitsString> > const & me)
	{
	SEQAN_CHECKPOINT
		return posGlobalize(value(host(me)), value(me.limits));
	}

}

#endif
