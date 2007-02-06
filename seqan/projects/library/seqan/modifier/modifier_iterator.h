/*
 *  modifier_iterator.h
 *  genindex
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
	struct ModifiedIterator;

	template < typename THost, typename TSpec >
	struct Spec< ModifiedIterator<THost, TSpec> > {
		typedef TSpec Type;
	}


	template < typename THost, typename TSpec >
	struct Value< ModifiedIterator<THost, TSpec> >:
		Value<THost> {};

	template < typename THost, typename TSpec >
	struct Value< ModifiedIterator<THost, TSpec> const >:
		Value<THost const> {};


	template < typename THost, typename TSpec >
	struct GetValue< ModifiedIterator<THost, TSpec> >:
		GetValue<THost> {};

	template < typename THost, typename TSpec >
	struct GetValue< ModifiedIterator<THost, TSpec> const >:
		GetValue<THost const> {};


	template < typename THost, typename TSpec >
	struct Reference< ModifiedIterator<THost, TSpec> >:
		Reference<THost> {};

	template < typename THost, typename TSpec >
	struct Reference< ModifiedIterator<THost, TSpec> const >:
		Reference<THost const> {};


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
	}

	template < typename THost, typename TSpec >
	struct Host< ModifiedIterator<THost, TSpec> const > {
		typedef THost const Type;
	}


	template < typename THost, typename TSpec >
	struct Container< ModifiedIterator<THost, TSpec> >:
		Container<THost> {};

	template < typename THost, typename TSpec >
	struct Container< ModifiedIterator<THost, TSpec> const >:
		Container<THost const> {};


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
	inline typename Reference<ModifiedIterator<THost, TSpec> >::Type 
	operator * (ModifiedIterator<THost, TSpec> & me)
	{
	SEQAN_CHECKPOINT
		return value(me);
	}

	template <typename THost, typename TSpec>
	inline typename Reference<ModifiedIterator<THost, TSpec> const>::Type 
	operator * (ModifiedIterator<THost, TSpec> const & me)
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
		return host(me) + delta;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator -
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename TDelta>
	inline ModifiedIterator<THost, TSpec>
	operator - (ModifiedIterator<THost, TSpec> const & me, TDelta delta) {
		return host(me) - delta;
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
	atBegin(ModifiedIterator<THost, TSpec> & const me,
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
	atBegin(ModifiedIterator<THost, TSpec> & const me)
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
	atEnd(ModifiedIterator<THost, TSpec> & const me,
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
	atEnd(ModifiedIterator<THost, TSpec> & const me)
	{
	SEQAN_CHECKPOINT
		return atEnd(host(me));
	}

}

#endif
