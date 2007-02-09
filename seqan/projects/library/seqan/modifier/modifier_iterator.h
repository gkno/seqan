/*
 *  modifier_iterator.h
 *  SeqAn
 *
 *  Created by David Weese on 26.01.06.
 *
 */

#ifndef SEQAN_HEADER_MODIFIER_ITERATOR_H
#define SEQAN_HEADER_MODIFIER_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TSpec = void >
	class ModifiedIterator {
	public:
		Holder<THost, Simple>					data_host;
		typename Cargo<ModifiedIterator>::Type	data_cargo;

		ModifiedIterator() {}

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

	template < typename THost, typename TSpec >
	struct Spec< ModifiedIterator<THost, TSpec> > {
		typedef TSpec Type;
	};

	template < typename THost, typename TSpec >
	struct Spec< ModifiedIterator<THost, TSpec> const > {
		typedef TSpec Type;
	};


	// an iterator is not the owner of the values pointing at
	// it can be constant while
	// - pointing to an alterable object
	// - returning an non-constant value
	// - being an iterator of an alterable container

	template < typename THost, typename TSpec >
	struct Value< ModifiedIterator<THost, TSpec> >:
		Value<THost> {};

	template < typename THost, typename TSpec >
	struct Value< ModifiedIterator<THost, TSpec> const >:
		Value< ModifiedIterator<THost, TSpec> > {};


	template < typename THost, typename TSpec >
	struct GetValue< ModifiedIterator<THost, TSpec> >:
		GetValue<THost> {};

	template < typename THost, typename TSpec >
	struct GetValue< ModifiedIterator<THost, TSpec> const >:
		GetValue< ModifiedIterator<THost, TSpec> > {};


	template < typename THost, typename TSpec >
	struct Reference< ModifiedIterator<THost, TSpec> >:
		Reference<THost> {};

	template < typename THost, typename TSpec >
	struct Reference< ModifiedIterator<THost, TSpec> const >:
		Reference< ModifiedIterator<THost, TSpec> > {};

	template < typename THost, typename TSpec >
	struct Size< ModifiedIterator<THost, TSpec> >:
		Size<THost> {};

	template < typename THost, typename TSpec >
	struct Position< ModifiedIterator<THost, TSpec> >:
		Position<THost> {};

	template < typename THost, typename TSpec >
	struct Difference< ModifiedIterator<THost, TSpec> >:
		Difference<THost> {};


	template < typename THost, typename TSpec >
	struct Host< ModifiedIterator<THost, TSpec> > {
		typedef THost Type;
	};

	template < typename THost, typename TSpec >
	struct Host< ModifiedIterator<THost, TSpec> const > {
		typedef THost const Type;
	};


	template < typename THost, typename TSpec >
	struct Container< ModifiedIterator<THost, TSpec> >:
		Container<THost> {};

	template < typename THost, typename TSpec >
	struct Container< ModifiedIterator<THost, TSpec> const >:
		Container< ModifiedIterator<THost, TSpec> > {};


	//////////////////////////////////////////////////////////////////////////////
	// host interface
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline Holder<THost, Simple> &
	_dataHost(ModifiedIterator<THost, TSpec> & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_host;
	}
	
	template <typename THost, typename TSpec>
	inline Holder<THost, Simple> const &
	_dataHost(ModifiedIterator<THost, TSpec> const & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_host;
	}

	template <typename THost, typename TSpec>
	inline typename Reference< typename Cargo<ModifiedIterator<THost, TSpec> >::Type >::Type
	cargo(ModifiedIterator<THost, TSpec> & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_cargo;
	}

	template <typename THost, typename TSpec>
	inline typename Reference< typename Cargo<ModifiedIterator<THost, TSpec> const>::Type >::Type
	cargo(ModifiedIterator<THost, TSpec> const & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_cargo;
	}

	//////////////////////////////////////////////////////////////////////////////
	// assign
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedIterator<THost, TSpec> const &
	assign(ModifiedIterator<THost, TSpec> & me, ModifiedIterator<THost2, TSpec> & _origin) {
		host(me) = host(_origin);
		cargo(me) = cargo(_origin);
		return me;
	}

	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedIterator<THost, TSpec> const &
	assign(ModifiedIterator<THost, TSpec> & me, ModifiedIterator<THost2, TSpec> const & _origin) {
		host(me) = host(_origin);
		cargo(me) = cargo(_origin);
		return me;
	}

	template <typename THost, typename TSpec, typename T>
	inline ModifiedIterator<THost, TSpec> const &
	assign(ModifiedIterator<THost, TSpec> & me, T & _origin) {
		host(me) = _origin;
		return me;
	}

	template <typename THost, typename TSpec, typename T>
	inline ModifiedIterator<THost, TSpec> const &
	assign(ModifiedIterator<THost, TSpec> & me, T const & _origin) {
		host(me) = _origin;
		return me;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator *
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline typename Reference<ModifiedIterator<THost, TSpec> >::Type 
	value(ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		return value(host(me));
	}

	template <typename THost, typename TSpec>
	inline typename Reference<ModifiedIterator<THost, TSpec> const>::Type 
	operator * (ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		return value(me);
	}


	//////////////////////////////////////////////////////////////////////////////
	// operator ++
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline void
	goNext(ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		goNext(host(me));
	}

	template <typename THost, typename TSpec>
	inline ModifiedIterator<THost, TSpec> const &
	operator ++ (ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		goNext(me);
		return me;
	}

	template <typename THost, typename TSpec>
	inline ModifiedIterator<THost, TSpec> const &
	operator ++ (ModifiedIterator<THost, TSpec> & me, int)
	{
	SEQAN_CHECKPOINT
		ModifiedIterator<THost, TSpec> temp_(me);
		goNext(me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator --
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline void
	goPrevious(ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		goPrevious(host(me));
	}

	template <typename THost, typename TSpec>
	inline ModifiedIterator<THost, TSpec> const &
	operator -- (ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		goPrevious(me);
		return me;
	}

	template <typename THost, typename TSpec>
	inline ModifiedIterator<THost, TSpec> const &
	operator -- (ModifiedIterator<THost, TSpec> & me, int)
	{
	SEQAN_CHECKPOINT
		ModifiedIterator<THost, TSpec> temp_(me);
		goPrevious(me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator +
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec>
	operator + (ModifiedIterator<THost, TSpec> const & me, TDelta delta) {
		ModifiedIterator<THost, TSpec> temp_(me);
		host(temp_) = host(me) + delta;
		return temp_;
	}

	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec>
	operator += (ModifiedIterator<THost, TSpec> & me, TDelta delta) {
		host(me) += delta;
		return me;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator -
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec>
	operator - (ModifiedIterator<THost, TSpec> const & me, TDelta delta) {
		ModifiedIterator<THost, TSpec> temp_(me);
		host(temp_) = host(me) - delta;
		return temp_;
	}

	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec>
	operator -= (ModifiedIterator<THost, TSpec> & me, TDelta delta) {
		host(me) -= delta;
		return me;
	}

	template <typename THost, typename TSpec>
	inline typename Difference< ModifiedIterator<THost, TSpec> >::Type
	operator - (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
		return host(a) - host(b);
	}

	//////////////////////////////////////////////////////////////////////////////
	// goBegin
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TContainer>
	inline void
	goBegin(ModifiedIterator<THost, TSpec> & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		host(me) = begin(container);
	}

	template <typename THost, typename TSpec>
	inline void
	goBegin(ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		goBegin(me, container(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// goEnd
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TContainer>
	inline void
	goEnd(ModifiedIterator<THost, TSpec> & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		host(me) = end(container);
	}

	template <typename THost, typename TSpec>
	inline void
	goEnd(ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		goEnd(me, container(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// position
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline typename Position<ModifiedIterator<THost, TSpec> const>::Type 
	position(ModifiedIterator<THost, TSpec> const & me)
	{
	SEQAN_CHECKPOINT
		return position(host(me));
	}

	template <typename THost, typename TSpec, typename TContainer>
	inline typename Position<ModifiedIterator<THost, TSpec> const>::Type 
	position(ModifiedIterator<THost, TSpec> const & me, TContainer const &cont)
	{
	SEQAN_CHECKPOINT
		return position(host(me), cont);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator ==
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline bool
	operator == (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
		return host(a) == host(b);
	}

	template <typename THost, typename TSpec>
	inline bool
	operator != (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
		return !(a == b);
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator <
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline bool
	operator < (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
		return host(a) < host(b);
	}

	template <typename THost, typename TSpec>
	inline bool
	operator > (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b) {
		return host(a) > host(b);
	}

	//////////////////////////////////////////////////////////////////////////////
	// atBegin
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TContainer>
	inline bool
	atBegin(ModifiedIterator<THost, TSpec> & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		return atBegin(host(me), container);
	}

	template <typename THost, typename TSpec, typename TContainer>
	inline bool
	atBegin(ModifiedIterator<THost, TSpec> const & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		return atBegin(host(me), container);
	}

	template <typename THost, typename TSpec>
	inline bool
	atBegin(ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		return atBegin(host(me));
	}

	template <typename THost, typename TSpec>
	inline bool
	atBegin(ModifiedIterator<THost, TSpec> const & me)
	{
	SEQAN_CHECKPOINT
		return atBegin(host(me));
	}

	//////////////////////////////////////////////////////////////////////////////
	// atEnd
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TContainer>
	inline bool
	atEnd(ModifiedIterator<THost, TSpec> & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		return atEnd(host(me), container);
	}

	template <typename THost, typename TSpec, typename TContainer>
	inline bool
	atEnd(ModifiedIterator<THost, TSpec> const & me,
			TContainer const & container)
	{
	SEQAN_CHECKPOINT
		return atEnd(host(me), container);
	}

	template <typename THost, typename TSpec>
	inline bool
	atEnd(ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		return atEnd(host(me));
	}

	template <typename THost, typename TSpec>
	inline bool
	atEnd(ModifiedIterator<THost, TSpec> const & me)
	{
	SEQAN_CHECKPOINT
		return atEnd(host(me));
	}

}

#endif
