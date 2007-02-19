/*
 *  modifier_view.h
 *  SeqAn
 *
 *  Created by David Weese on 26.01.06.
 *
 */

#ifndef SEQAN_HEADER_MODIFIER_VIEW_H
#define SEQAN_HEADER_MODIFIER_VIEW_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Spec.ModView:
..summary:Transforms the characters of the $THost$ string/iterator using a custom function.
..cat:Modifiers
..general:Class.ModifiedIterator
..general:Class.ModifiedString
..signature:ModifiedIterator<THost, ModView<TFunctor> >
..signature:ModifiedString<THost, ModView<TFunctor> >
..param.THost:Original string/iterator.
...type:Concept.Iterator
..param.TFunctor:A unary function (see STL's $unary_function$).
...remarks:The argument type of $TFunctor$ must be $VALUE<THost>::Type$.
..remarks:The @Metafunction.Value@ type of this modifier is the result type of $TFunctor$.
*/

	template <typename TFunctor>
	struct ModView {};

	template <typename TFunctor>
	struct ModViewCargo {
		TFunctor	func;
	};


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// view iterator
	//////////////////////////////////////////////////////////////////////////////


	template <typename THost, typename TFunctor>
	struct Cargo< ModifiedIterator<THost, ModView<TFunctor> > > {
		typedef ModViewCargo<TFunctor>	Type;
	};

	template <typename THost, typename TFunctor>
	class ModifiedIterator<THost, ModView<TFunctor> > {
	public:
		Holder<THost, Simple>					data_host;
		typename Cargo<ModifiedIterator>::Type	data_cargo;

		mutable typename Value<ModifiedIterator>::Type	tmp_value;

		ModifiedIterator() {}

		explicit ModifiedIterator(TFunctor &_func) {
			assignModViewFunctor(*this, _func);
		}

		explicit ModifiedIterator(TFunctor const &_func) {
			assignModViewFunctor(*this, _func);
		}

		ModifiedIterator(ModifiedIterator &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		ModifiedIterator(ModifiedIterator const &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		template <typename T>
		ModifiedIterator(T & _origin) {
			assign(*this, _origin);
		}

		template <typename T>
		ModifiedIterator(T const & _origin) {
			assign(*this, _origin);
		}
//____________________________________________________________________________

		template <typename T>
		inline ModifiedIterator const &
		operator = (T & _origin) {
			assign(*this, _origin);
			return *this;
		}

		template <typename T>
		inline ModifiedIterator const &
		operator = (T const & _origin) {
			assign(*this, _origin);
			return *this;
		}
	};

	template <typename THost, typename TFunctor>
	struct Value< ModifiedIterator<THost, ModView<TFunctor> > > {
		typedef typename TFunctor::result_type	Type;
	};

	template <typename THost, typename TFunctor>
	struct GetValue< ModifiedIterator<THost, ModView<TFunctor> > >:
		Value< ModifiedIterator<THost, ModView<TFunctor> > > {};

	template <typename THost, typename TFunctor>
	struct Reference< ModifiedIterator<THost, ModView<TFunctor> > > {
		typedef typename Value< ModifiedIterator<THost, ModView<TFunctor> > >::Type & Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// operator *
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TFunctor>
	inline typename Reference<ModifiedIterator<THost, ModView<TFunctor> > >::Type 
	value(ModifiedIterator<THost, ModView<TFunctor> > & me)
	{
	SEQAN_CHECKPOINT
		me.tmp_value = cargo(me).func(value(host(me)));
		return me.tmp_value;
	}

	template <typename THost, typename TFunctor>
	inline typename Reference<ModifiedIterator<THost, ModView<TFunctor> > const>::Type 
	value(ModifiedIterator<THost, ModView<TFunctor> > const & me)
	{
	SEQAN_CHECKPOINT
		me.tmp_value = cargo(me).func(value(host(me)));
		return me.tmp_value;
	}

	//////////////////////////////////////////////////////////////////////////////
	// assignModViewFunctor
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TFunctor>
	inline void
	assignModViewFunctor(ModifiedIterator<THost, ModView<TFunctor> > & me, TFunctor const & _func) 
	{
	SEQAN_CHECKPOINT
		cargo(me).func = _func;
	}


//////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	// view string
	//////////////////////////////////////////////////////////////////////////////


	template <typename THost, typename TFunctor>
	struct Cargo< ModifiedString<THost, ModView<TFunctor> > > {
		typedef ModViewCargo<TFunctor>	Type;
	};

	template <typename THost, typename TFunctor>
	class ModifiedString<THost, ModView<TFunctor> > {
	public:
		Holder<THost>							data_host;
		typename Cargo<ModifiedString>::Type	data_cargo;

		ModifiedString() {}

		explicit ModifiedString(TFunctor &_func) {
			cargo(*this).func = _func;
		}

		explicit ModifiedString(TFunctor const &_func) {
			cargo(*this).func = _func;
		}

		explicit ModifiedString(ModifiedString const &_origin, TFunctor const &_func):
			data_host(_origin.data_host)
		{
			cargo(*this).func = _func;
		}

		ModifiedString(ModifiedString &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		ModifiedString(ModifiedString const &_origin):
			data_host(_origin.data_host),
			data_cargo(_origin.data_cargo) {}

		template <typename T>
		ModifiedString(T & _origin) {
			assign(*this, _origin);
		}

		template <typename T>
		inline ModifiedString const &
		operator = (T & _origin) {
			assign(*this, _origin);
			return *this;
		}

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


	//////////////////////////////////////////////////////////////////////////////
	// assignModViewFunctor
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TFunctor>
	inline void
	assignModViewFunctor(ModifiedString<THost, ModView<TFunctor> > & me, TFunctor const & _func)
	{
	SEQAN_CHECKPOINT
		cargo(me).func = _func;
	}


}

#endif
