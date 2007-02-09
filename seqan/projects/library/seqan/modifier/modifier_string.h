/*
 *  modifier_string.h
 *  SeqAn
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
	class ModifiedString {
	public:
		Holder<THost>							data_host;
		typename Cargo<ModifiedString>::Type	data_cargo;

		ModifiedString() {}

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
	};

	template < typename THost, typename TSpec >
	struct Spec< ModifiedString<THost, TSpec> > {
		typedef TSpec Type;
	};

	template < typename THost, typename TSpec >
	struct Spec< ModifiedString<THost, TSpec> const > {
		typedef TSpec Type;
	};

/*
	template < typename THost, typename TSpec >
	struct Value< ModifiedString<THost, TSpec> >:
		Value<THost> {};

	template < typename THost, typename TSpec >
	struct Value< ModifiedString<THost, TSpec> const>:
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
*/


	// use Value, GetValue, Reference, Size, ... from corresponding iterator
	template < typename THost, typename TSpec >
	struct Value< ModifiedString<THost, TSpec> >:
		Value< typename Iterator< ModifiedString<THost, TSpec> >::Type > {};

	template < typename THost, typename TSpec >
	struct Value< ModifiedString<THost, TSpec> const >:
		Value< typename Iterator< ModifiedString<THost, TSpec> const >::Type > {};


	template < typename THost, typename TSpec >
	struct GetValue< ModifiedString<THost, TSpec> >:
		GetValue< typename Iterator< ModifiedString<THost, TSpec> >::Type > {};

	template < typename THost, typename TSpec >
	struct GetValue< ModifiedString<THost, TSpec> const >:
		GetValue< typename Iterator< ModifiedString<THost, TSpec> const >::Type > {};


	template < typename THost, typename TSpec >
	struct Reference< ModifiedString<THost, TSpec> >:
		Reference< typename Iterator< ModifiedString<THost, TSpec> >::Type > {};

	template < typename THost, typename TSpec >
	struct Reference< ModifiedString<THost, TSpec> const >:
		Reference< typename Iterator< ModifiedString<THost, TSpec> const >::Type > {};

	template < typename THost, typename TSpec >
	struct Size< ModifiedString<THost, TSpec> >:
		Size< typename Iterator< ModifiedString<THost, TSpec> >::Type > {};

	template < typename THost, typename TSpec >
	struct Position< ModifiedString<THost, TSpec> >:
		Position< typename Iterator< ModifiedString<THost, TSpec> >::Type > {};

	template < typename THost, typename TSpec >
	struct Difference< ModifiedString<THost, TSpec> >:
		Difference< typename Iterator< ModifiedString<THost, TSpec> >::Type > {};




	template <typename THost, typename TSpec>
	struct Iterator< ModifiedString<THost, TSpec>, Standard > {
		typedef ModifiedIterator<typename Iterator<THost, Standard>::Type, TSpec> Type;
	};

	template <typename THost, typename TSpec >
	struct Iterator< ModifiedString<THost, TSpec> const, Standard > {
		typedef ModifiedIterator<typename Iterator<THost const, Standard>::Type, TSpec> Type;
	};

	template <typename THost, typename TSpec>
	struct Iterator< ModifiedString<THost, TSpec>, Rooted > {
		typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, TSpec> Type;
	};

	template <typename THost, typename TSpec >
	struct Iterator< ModifiedString<THost, TSpec> const, Rooted > {
		typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, TSpec> Type;
	};


	template < typename THost, typename TSpec >
	struct Host< ModifiedString<THost, TSpec> > {
		typedef THost Type;
	};

	template < typename THost, typename TSpec >
	struct Host< ModifiedString<THost, TSpec> const > {
		typedef THost const Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// host interface
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec>
	inline Holder<THost> &
	_dataHost(ModifiedString<THost, TSpec> & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_host;
	}
	
	template <typename THost, typename TSpec>
	inline Holder<THost> const &
	_dataHost(ModifiedString<THost, TSpec> const & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_host;
	}

	template <typename THost, typename TSpec>
	inline typename Reference< typename Cargo<ModifiedString<THost, TSpec> >::Type >::Type
	cargo(ModifiedString<THost, TSpec> & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_cargo;
	}

	template <typename THost, typename TSpec>
	inline typename Reference< typename Cargo<ModifiedString<THost, TSpec> const>::Type >::Type
	cargo(ModifiedString<THost, TSpec> const & me) 
	{
	SEQAN_CHECKPOINT
		return me.data_cargo;
	}

	//////////////////////////////////////////////////////////////////////////////
	// assign
	//////////////////////////////////////////////////////////////////////////////

	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedString<THost, TSpec> const &
	assign(ModifiedString<THost, TSpec> & me, ModifiedString<THost2, TSpec> & _origin) {
		host(me) = host(_origin);
		cargo(me) = cargo(_origin);
		return me;
	}

	template <typename THost, typename TSpec, typename THost2>
	inline ModifiedString<THost, TSpec> const &
	assign(ModifiedString<THost, TSpec> & me, ModifiedString<THost2, TSpec> const & _origin) {
		host(me) = host(_origin);
		cargo(me) = cargo(_origin);
		return me;
	}

	template <typename THost, typename TSpec, typename T>
	inline ModifiedString<THost, TSpec> const &
	assign(ModifiedString<THost, TSpec> & me, T & _origin) {
		host(me) = _origin;
		return me;
	}

	template <typename THost, typename TSpec, typename T>
	inline ModifiedString<THost, TSpec> const &
	assign(ModifiedString<THost, TSpec> & me, T const & _origin) {
		host(me) = _origin;
		return me;
	}

	//////////////////////////////////////////////////////////////////////////////
	// value
	//////////////////////////////////////////////////////////////////////////////
/*
class type_info {
public:
    _CRTIMP virtual ~type_info();
    _CRTIMP bool operator==(const type_info& rhs) const;
    _CRTIMP bool operator!=(const type_info& rhs) const;
    _CRTIMP int before(const type_info& rhs) const;
    _CRTIMP const char* name() const;
    _CRTIMP const char* raw_name() const;
};
	template <typename THost, typename TSpec, typename TPos>
	inline typename Reference<ModifiedString<THost, TSpec> >::Type 
	value(ModifiedString<THost, TSpec> & me, TPos pos)
	{
	SEQAN_CHECKPOINT
//		return typename Reference<ModifiedString<THost, TSpec> >::Type = *(begin(me) + pos);
	typedef typename Iterator<ModifiedString<THost, TSpec> >::Type TIter;
	typedef ModifiedIterator<typename Iterator<THost>::Type, TSpec> TIterB;
	typedef typename Reference<TIter>::Type TIter1;
	typedef typename Reference<TIterB>::Type TIter2;
	typedef typename Reference<ModifiedString<THost, TSpec> >::Type TIter3;
	typedef unsigned & TRef1;
	typedef Reference<unsigned>::Type TRef2;
		TIter it = begin(me);
		it = it + pos;
//		::std::cout << value(it);
		::std::cout << typeid(TIter).name() << ::std::endl;
		::std::cout << typeid(TIterB).name() << ::std::endl;
		::std::cout << typeid(TIter1).name() << ::std::endl;
		::std::cout << typeid(TIter2).name() << ::std::endl;
		::std::cout << typeid(TIter3).name() << ::std::endl;
		::std::cout << typeid(TRef1).name() << ::std::endl;
		::std::cout << typeid(TRef2).name() << ::std::endl;
//		typename Reference<TIter>::Type ref=value(it);
//	static typename Value<TIter>::Type ref = 0;

//		return *it ;
		return *(begin(me) + pos);
	}
*/
	template <typename THost, typename TSpec, typename TPos>
	inline typename Reference<ModifiedString<THost, TSpec> >::Type 
	value(ModifiedString<THost, TSpec> & me, TPos pos)
	{
		return value(begin(me) + pos);
	}
/*
	template <typename THost, typename TSpec, typename TPos>
	inline typename Reference<ModifiedString<THost, TSpec> const >::Type 
	value(ModifiedString<THost, TSpec> const & me, TPos pos)
	{
	SEQAN_CHECKPOINT
		return value(begin(me) + pos);
	}
*/
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
		typename Iterator< ModifiedString<THost, TSpec> const >::Type temp_(begin(host(me)));
		cargo(temp_) = cargo(me);
		return temp_;
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
	begin(ModifiedString<THost, TSpec> & me) {
		typename Iterator< ModifiedString<THost, TSpec> >::Type temp_(begin(host(me)));
		cargo(temp_) = cargo(me);
		return temp_;
	}

	template < typename THost, typename TSpec, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type 
	begin(ModifiedString<THost, TSpec> const & me, Tag<TTagSpec> const tag_) {
		typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type temp_(begin(host(me), tag_));
		cargo(temp_) = cargo(me);
		return temp_;
	}

	template < typename THost, typename TSpec, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type 
	begin(ModifiedString<THost, TSpec> & me, Tag<TTagSpec> const tag_) {
		typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type temp_(begin(host(me), tag_));
		cargo(temp_) = cargo(me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// end
	//////////////////////////////////////////////////////////////////////////////

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const >::Type 
	end(ModifiedString<THost, TSpec> const & me) {
		typename Iterator< ModifiedString<THost, TSpec> const >::Type temp_(end(host(me)));
		cargo(temp_) = cargo(me);
		return temp_;
	}

	template < typename THost, typename TSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
	end(ModifiedString<THost, TSpec> & me) {
		typename Iterator< ModifiedString<THost, TSpec> >::Type temp_(end(host(me)));
		cargo(temp_) = cargo(me);
		return temp_;
	}

	template < typename THost, typename TSpec, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type 
	end(ModifiedString<THost, TSpec> const & me, Tag<TTagSpec> const tag_) {
		typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type temp_(end(host(me), tag_));
		cargo(temp_) = cargo(me);
		return temp_;
	}

	template < typename THost, typename TSpec, typename TTagSpec >
	inline typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type 
	end(ModifiedString<THost, TSpec> & me, Tag<TTagSpec> const tag_) {
		typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type temp_(end(host(me), tag_));
		cargo(temp_) = cargo(me);
		return temp_;
	}

	//////////////////////////////////////////////////////////////////////////////
	// stream operators
	//////////////////////////////////////////////////////////////////////////////

	template < typename TStream, typename THost, typename TSpec >
	inline TStream &
	operator << (TStream & target, ModifiedString<THost, TSpec> const & source)
	{
	SEQAN_CHECKPOINT
		write(target, source);
		return target;
	}

	//////////////////////////////////////////////////////////////////////////////

	template < typename TStream, typename THost, typename TSpec >
	inline TStream &
	operator >> (TStream & source, ModifiedString<THost, TSpec> & target)
	{
	SEQAN_CHECKPOINT
		read(source, target);
		return source;
	}


}

#endif
