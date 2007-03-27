/*
 *  string_set.h
 *  SeqAn
 *
 *  Created by David Weese on 26.01.06.
 *
 */

#ifndef SEQAN_HEADER_STRING_SET_H
#define SEQAN_HEADER_STRING_SET_H

namespace SEQAN_NAMESPACE_MAIN
{

	template < typename TSpec = void >
	struct ConcatVirtual;

	template < typename TLimiter = void >
	struct ConcatDirect;

    template < typename TStringSet >
    class StringSetIterator;


//	template < typename TString, typename TSpec = ConcatVirtual<> >
	template < typename TString, typename TSpec = ConcatDirect<> >
	class StringSet;

    template <typename TObject>
	struct Concatenator {
		typedef TObject Type;
	};

    template <typename TObject>
	struct Concatenator<TObject const> {
		typedef typename Concatenator<TObject>::Type const Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// StringSet limits
	//////////////////////////////////////////////////////////////////////////////

	template <typename TString>
	struct StringSetLimits {
		typedef Nothing Type;
	};

	template <typename TString>
	struct StringSetLimits<TString const> {
		typedef typename StringSetLimits<TString>::Type const Type;
	};

	template <typename TString>
	struct StringSetPosition {
		typedef typename Size<TString>::Type	Type;
	};

	template <typename TString, typename TSpec>
	struct StringSetLimits< StringSet<TString, TSpec> > {
		typedef typename Size<TString>::Type	TSize;
		typedef String<TSize>					Type;
	};

	template <typename TString, typename TSpec>
	struct StringSetPosition< StringSet<TString, TSpec> > {
		typedef typename Size<TString>::Type	TSize;
		typedef Pair<TSize>						Type;
	};



	//////////////////////////////////////////////////////////////////////////////
	// get StringSet limits
	//////////////////////////////////////////////////////////////////////////////

    template <typename TStringSet>
	inline typename StringSetLimits<TStringSet>::Type
	stringSetLimits(TStringSet &stringSet) {
		return Nothing();
	}

    template <typename TStringSet>
	inline typename StringSetLimits<TStringSet const>::Type
	stringSetLimits(TStringSet const &stringSet) {
		return Nothing();
	}

    template <typename TString, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSpec> >::Type & 
	stringSetLimits(StringSet<TString, TSpec> &stringSet) {
		if (!_validStringSetLimits(stringSet))
			_refreshStringSetLimits(stringSet);
		return stringSet.limits;
	}

    template <typename TString, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSpec> const>::Type & 
	stringSetLimits(StringSet<TString, TSpec> const &stringSet) {
		if (!_validStringSetLimits(stringSet))
			_refreshStringSetLimits(const_cast< StringSet<TString, TSpec>& >(stringSet));
		return stringSet.limits;
	}

	//////////////////////////////////////////////////////////////////////////////
	// 1 sequence
	//////////////////////////////////////////////////////////////////////////////

	template <typename TPosition>
	inline TPosition getSeqNo(TPosition const &pos, Nothing const &) {
		return 0;
	}

	template <typename TPosition>
	inline TPosition getSeqOffset(TPosition const &pos, Nothing const &) {
		return pos;
	}

	//////////////////////////////////////////////////////////////////////////////
	// n sequences (position type is Pair)
	//////////////////////////////////////////////////////////////////////////////

	template <typename T1, typename T2, typename TCompression, typename TLimitsString>
	inline T1 getSeqNo(Pair<T1, T2, TCompression> const &pos, TLimitsString const &) {
		return getValueI1(pos);
	}

	template <typename T1, typename T2, typename TCompression, typename TLimitsString>
	inline T2 getSeqOffset(Pair<T1, T2, TCompression> const &pos, TLimitsString const &) {
		return getValueI2(pos);
	}

	//////////////////////////////////////////////////////////////////////////////
	// n sequences (position type is an integral type)
	//////////////////////////////////////////////////////////////////////////////

	template <typename TPos, typename TLimitsString>
	inline TPos getSeqNo(TPos const &pos, TLimitsString const &limits) {
		typedef typename Iterator<TLimitsString const, Standard>::Type TIter;
        TIter _begin = begin(limits, Standard());
        TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), pos) - 1;
        return distance(_begin, _upper);
	}

	template <typename TPos, typename TLimitsString>
	inline TPos getSeqOffset(TPos const &pos, TLimitsString const &limits) {
		typedef typename Iterator<TLimitsString const, Standard>::Type TIter;
        TIter _begin = begin(limits, Standard());
        TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), pos) - 1;
        return pos - *_upper;
	}

	//////////////////////////////////////////////////////////////////////////////
	// local -> global position
	//////////////////////////////////////////////////////////////////////////////

	template <typename TPosition>
	inline TPosition posGlobalize(TPosition const &pos, Nothing const &) {
		return pos;
	}

	template <typename TLimitsString, typename T1, typename T2, typename TCompression>
	inline typename Value<TLimitsString>::Type 
	posGlobalize(Pair<T1, T2, TCompression> const &pos, TLimitsString const &limits) {
		return limits[getSeqNo(pos, limits)] + getSeqOffset(pos, limits);
	}

	//////////////////////////////////////////////////////////////////////////////
	// global -> local position
	//////////////////////////////////////////////////////////////////////////////

	template <typename TResult, typename TPosition>
	inline void posLocalize(TResult &result, TPosition const &pos, Nothing const &) {
		result = pos;
	}

	template <typename TResult, typename TSize, typename TSpec, typename TPosition>
	inline void posLocalize(TResult &result, TPosition const &pos, String<TSize, TSpec> const &limits) {
		typedef typename Iterator<String<TSize> const, Standard>::Type TIter;
        TIter _begin = begin(limits, Standard());
		TIter _upper = ::std::upper_bound(_begin, end(limits, Standard()), pos) - 1;
        result.i1 = distance(_begin, _upper);
        result.i2 = pos - *_upper;
	}

	//////////////////////////////////////////////////////////////////////////////
	// position arithmetics
	//////////////////////////////////////////////////////////////////////////////

	// posAtFirstLocal
	template <typename TPos, typename TLimitsString>
	inline bool posAtFirstLocal(TPos pos, TLimitsString const &limits) {
		return getSeqOffset(pos, limits) == 0;
	}

	// posPrev
	template <typename TPos>
	inline TPos posPrev(TPos pos) {
		return pos - 1;
	}

	template <typename T1, typename T2, typename TCompression>
	inline Pair<T1, T2, TCompression> posPrev(Pair<T1, T2, TCompression> const &pos) {
		return Pair<T1, T2, TCompression>(getValueI1(pos), getValueI2(pos) - 1);
	}

	// posNext
	template <typename TPos>
	inline TPos posNext(TPos pos) {
		return pos + 1;
	}

	template <typename T1, typename T2, typename TCompression>
	inline Pair<T1, T2, TCompression> 
	posNext(Pair<T1, T2, TCompression> const &pos) {
		return Pair<T1, T2, TCompression>(getValueI1(pos), getValueI2(pos) + 1);
	}

	// posAdd
	template <typename TPos, typename TDelta>
	inline TPos posAdd(TPos pos, TDelta delta) {
		return pos + delta;
	}

	template <typename T1, typename T2, typename TCompression, typename TDelta>
	inline Pair<T1, T2, TCompression> 
	posAdd(Pair<T1, T2, TCompression> const &pos, TDelta delta) {
		return Pair<T1, T2, TCompression>(getValueI1(pos), getValueI2(pos) + delta);
	}


	//////////////////////////////////////////////////////////////////////////////
    // StringSet Container
    //////////////////////////////////////////////////////////////////////////////

    template < typename TString, typename TSpec >
    class StringSet< TString, ConcatVirtual<TSpec> >
    {
	public:

        typedef String<TString>										TStrings;
		typedef typename StringSetLimits<StringSet>::Type			TLimits;
		typedef typename Iterator<TLimits const, Standard>::Type	TLimitIterator;

		typedef typename Value<TString>::Type			Type;
        typedef typename Size<TString>::Type			TSize;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator types

        friend class StringSetIterator<StringSet>;

		typedef StringSetIterator<StringSet>	StringSetConstIterator;
		typedef StringSetIterator<StringSet>	StringSetIterator;

		// used by Iter<.., ConcatVirtual<..> >
		typedef StringSetIterator		iterator;
		typedef StringSetConstIterator	const_iterator;
        typedef Type					value_type;
		typedef Type&					reference;
		typedef Type const &			const_reference;
		typedef Type*					pointer;
		typedef Type const *			const_pointer;
		typedef TSize				size_type;
		typedef TSize				difference_type;

        typedef typename Iterator<TString, Standard>::Type			obj_iterator;
        typedef typename Iterator<TString const, Standard>::Type	const_obj_iterator;

		TStrings	strings;
        TLimits		limits;
		bool		limitsValid;		// is true if limits contains the cumulative sum of the sequence lengths

		typename Concatenator<StringSet>::Type	concat;

		StringSet():
			limitsValid(true)
		{
			appendValue(limits, 0);
			concat.set = this;
		};
//____________________________________________________________________________

		template <typename TPos>
		inline typename Reference<StringSet>::Type
		operator [] (TPos pos)
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<StringSet const>::Type 
		operator [] (TPos pos) const
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};

    template < typename TString, typename TLimiter >
    class StringSet< TString, ConcatDirect<TLimiter> >
    {
	public:
		typedef typename StringSetLimits<StringSet>::Type			TLimits;
        typedef typename Iterator<TLimits const, Standard>::Type	TLimitIterator;

		typedef typename Value<TString>::Type			TValue;
        typedef typename Size<TString>::Type			TSize;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator types

        friend class StringSetIterator<StringSet>;

		typedef StringSetIterator<StringSet>	StringSetConstIterator;
		typedef StringSetIterator<StringSet>	StringSetIterator;

		// used by Iter<.., ConcatVirtual<..> >
        typedef TValue					value_type;
		typedef TValue&					reference;
		typedef TValue*					pointer;
		typedef TSize					size_type;

		// used by Iter<.., ConcatVirtual<..> >
        typedef typename Iterator<TString, Standard>::Type			obj_iterator;
        typedef typename Iterator<TString const, Standard>::Type	const_obj_iterator;

        TLimits		limits;

		typename Concatenator<StringSet>::Type	concat;

		StringSet()	{
			appendValue(limits, 0);
		}

//____________________________________________________________________________

		template <typename TPos>
		inline typename Reference<StringSet>::Type
		operator [] (TPos pos)
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<StringSet const>::Type 
		operator [] (TPos pos) const
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}
	};


    template < typename TString, typename TSpec >
    struct Value< StringSet< TString, TSpec > > {
        typedef TString Type;
    };

    template < typename TString, typename TSpec >
    struct Size< StringSet< TString, TSpec > > {
		typedef typename Size< typename StringSetLimits< StringSet<TString, TSpec> >::Type >::Type Type;
    };

	// direct concatenation
	template < typename TString, typename TSpec >
    struct Value< StringSet< TString, ConcatDirect<TSpec> > > {
		typedef typename Infix<TString>::Type Type;
    };

    template < typename TString, typename TSpec >
    struct GetValue< StringSet< TString, ConcatDirect<TSpec> > > {
		typedef typename Infix<TString>::Type Type;
    };

    template < typename TString, typename TSpec >
    struct GetValue< StringSet< TString, ConcatDirect<TSpec> > const > {
		typedef typename Infix<TString const>::Type Type;
    };

    template < typename TString, typename TSpec >
    struct Reference< StringSet< TString, ConcatDirect<TSpec> > > {
		typedef typename Infix<TString>::Type Type;
    };

    template < typename TString, typename TSpec >
    struct Reference< StringSet< TString, ConcatDirect<TSpec> > const > {
		typedef typename Infix<TString const>::Type Type;
    };
//____________________________________________________________________________

	template < typename T >
    inline bool _validStringSetLimits(T const &me) {
        return true;
    }

	template < typename TString, typename TSpec >
    inline bool _validStringSetLimits(StringSet< TString, ConcatVirtual<TSpec> > const &me) {
        return me.limitsValid;
    }
//____________________________________________________________________________

	template < typename T >
    inline void _refreshStringSetLimits(T &me) {
	}

	template < typename TString, typename TSpec >
    inline void _refreshStringSetLimits(StringSet< TString, ConcatVirtual<TSpec> > &me) 
	{
		typedef StringSet< TString, ConcatVirtual<TSpec> >	TStringSet;
		typedef typename StringSetLimits<TStringSet>::Type	TLimits;

		typename Value<TLimits>::Type	sum = 0;
		typename Size<TStringSet>::Type	i = 0;

		for(; i < length(me); ++i) {
			me.limits[i] = sum;
			sum += length(me[i]);
		}
		me.limits[i] = sum;
		me.limitsValid = true;
    }
//____________________________________________________________________________

	template < typename TString, typename TSpec >
    inline typename Size<TString>::Type lengthSum(StringSet< TString, TSpec > const &me) {
		if (!_validStringSetLimits(me))
			_refreshStringSetLimits(const_cast< StringSet<TString, TSpec>& >(me));
        return back(me.limits);
    }
//____________________________________________________________________________

	template < typename TString, typename TSpec, typename TString2, typename TExpand >
    inline void appendValue(
		StringSet< TString, ConcatVirtual<TSpec> > &me, 
		TString2 const &obj,
		Tag<TExpand> const tag) 
	{
        appendValue(me.strings, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }

	template < typename TString, typename TSpec, typename TString2, typename TExpand >
    inline void appendValue(
		StringSet< TString, ConcatVirtual<TSpec> > &me, 
		TString2 &obj,
		Tag<TExpand> const tag) 
	{
        appendValue(me.strings, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }
/*
	template < typename TString, typename TSpec, typename TString2, typename TExpand >
    inline void appendValue(
		StringSet< TString, ConcatVirtual<TSpec> > &me, 
		TString2 &obj,
		Tag<TExpand> const tag) 
	{
        appendValue(me.strings, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }

	template < typename TString, typename TSpec >
    inline void appendValue(StringSet< TString, ConcatVirtual<TSpec> > &me, TString &obj) {
        appendValue(me.strings, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }
*/

	template < typename TString, typename TString2, typename TExpand >
    inline void appendValue(
		StringSet< TString, ConcatDirect<void> > &me, 
		TString2 &obj,
		Tag<TExpand> const tag) 
	{
        append(me.concat, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }

    template < typename TString, typename TLimiter, typename TString2, typename TExpand >
    inline void appendValue(
		StringSet< TString, ConcatDirect<TLimiter> > &me, 
		TString2 &obj,
		Tag<TExpand> const tag) 
	{
        append(me.concat_string, obj);
        appendValue(me.concat_string, TLimiter());
        appendValue(me.limits, lengthSum(me) + length(obj) + 1);
    }

  
/*
    inline void append(TString *_objs[], unsigned count) {
        for(unsigned i = 0; i < count; ++i)
            add(_objs[i]);
    }
*/

//____________________________________________________________________________

    template < typename TString, typename TSpec >
    inline typename Size<TString>::Type 
	clear(StringSet< TString, ConcatVirtual<TSpec> > &me) {
		clear(me.strings);
		resize(me.limits, 1);
		me.limits = 0;
		me.valid = true;
    }

    template < typename TString, typename TLimiter >
    inline typename Size<TString>::Type 
	clear(StringSet< TString, ConcatDirect<TLimiter> > &me) {
		clear(me.concat_string);
		resize(me.limits, 1);
		me.limits = 0;
    }
//____________________________________________________________________________

    template < typename TString, typename TSpec >
    inline typename Size< StringSet< TString, TSpec > >::Type 
	length(StringSet< TString, TSpec > const &me) {
        return length(me.limits) - 1;
    }
//____________________________________________________________________________

	template < typename TString, typename TSpec, typename TSize >
    inline typename Size< StringSet< TString, ConcatVirtual<TSpec> > >::Type 
	resize(StringSet< TString, ConcatVirtual<TSpec> > &me, TSize new_size) {
		resize(me.limits, new_size + 1);
		me.limitsValid = (new_size == 0);
		return resize(me.strings, new_size);
    }

	template < typename TString, typename TSpec, typename TSize >
    inline typename Size< StringSet< TString, ConcatDirect<TSpec> > >::Type 
	resize(StringSet< TString, ConcatDirect<TSpec> > &me, TSize new_size) {
		return resize(me.limits, new_size + 1) - 1;
    }
//____________________________________________________________________________
	
	template < typename TString, typename TSpec, typename TPos >
	inline typename Reference< StringSet< TString, ConcatVirtual<TSpec> > >::Type
	value(StringSet< TString, ConcatVirtual<TSpec> > & me, TPos pos)
	{
		return me.strings[pos];
	}

	template < typename TString, typename TSpec, typename TPos >
	inline typename Reference< StringSet< TString, ConcatVirtual<TSpec> > const >::Type
	value(StringSet< TString, ConcatVirtual<TSpec> > const & me, TPos pos)
	{
		return me.strings[pos];
	}

	template < typename TString, typename TSpec, typename TPos >
	inline typename Infix<TString>::Type
	value(StringSet< TString, ConcatDirect<TSpec> > & me, TPos pos)
	{
		return infix(me.concat, me.limits[pos], me.limits[pos + 1]);
	} 

	template < typename TString, typename TSpec, typename TPos >
	inline typename Infix<TString const>::Type
	value(StringSet< TString, ConcatDirect<TSpec> > const & me, TPos pos)
	{
		return infix(me.concat, me.limits[pos], me.limits[pos + 1]);
	} 

//____________________________________________________________________________



	//////////////////////////////////////////////////////////////////////////////
	// ConcatenatorNto1 - a StringSet to String converter
	//////////////////////////////////////////////////////////////////////////////

    template <typename TStringSet>
	struct ConcatenatorNto1 {
		TStringSet *set;
		ConcatenatorNto1 () {}
		ConcatenatorNto1 (TStringSet &_set): set(&_set) {}
	};

    template <typename TStringSet>
	struct Value< ConcatenatorNto1<TStringSet> > {
		typedef typename Value< typename Value<TStringSet>::Type >::Type Type;
	};

    template <typename TStringSet>
	struct Value< ConcatenatorNto1<TStringSet> const > {
		typedef typename Value< typename Value<TStringSet>::Type >::Type Type;
	};

    template <typename TStringSet>
    struct Size< ConcatenatorNto1<TStringSet> > {
		typedef typename Size< typename Value<TStringSet>::Type >::Type Type;
    };

    template <typename TString, typename TSpec >
    struct Iterator< ConcatenatorNto1< StringSet<TString, TSpec> > > {
        typedef Iter<StringSet<TString, TSpec>, TSpec> Type;
    };

    template <typename TString, typename TSpec >
    struct Iterator< ConcatenatorNto1< StringSet<TString, TSpec> const > > {
        typedef Iter<StringSet<TString, TSpec> const, TSpec> Type;
    };

    template <typename TStringSet >
    struct Iterator< ConcatenatorNto1<TStringSet> const > {
		typedef typename Iterator< ConcatenatorNto1<TStringSet> >::Type Type;
    };

//////////////////////////////////////////////////////////////////////////////
// length

	template < typename TString, typename TSpec >
    inline typename Size< ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > > >::Type 
	length(ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > > const &concat) {
        return lengthSum(*concat.set);
    }

//////////////////////////////////////////////////////////////////////////////
// begin

	template < typename TString, typename TSpec >
	inline typename Iterator< ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > const> >::Type
	begin(ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > const> concat)
	{
		return Iter< StringSet<TString, ConcatVirtual<TSpec> > const, ConcatVirtual<TSpec> > (*concat.set);
	}

	template < typename TString, typename TSpec >
	inline typename Iterator< ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > > >::Type
	begin(ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > > concat)
	{
		return Iter< StringSet<TString, ConcatVirtual<TSpec> >, ConcatVirtual<TSpec> > (*concat.set);
	}

//////////////////////////////////////////////////////////////////////////////
// end

	template < typename TString, typename TSpec >
	inline typename Iterator< ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > const> >::Type
	end(ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > const> concat) {
		return Iter< StringSet<TString, ConcatVirtual<TSpec> > const, ConcatVirtual<TSpec> > 
			(*concat.set, length(*concat.set), 0);
	}

	template < typename TString, typename TSpec >
	inline typename Iterator< ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > > >::Type
	end(ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > > concat) {
		return Iter< StringSet<TString, ConcatVirtual<TSpec> >, ConcatVirtual<TSpec> > 
			(*concat.set, length(*concat.set), 0);
	}

//////////////////////////////////////////////////////////////////////////////
// Concatenator metafunction

	template < typename TString, typename TSpec >
	struct Concatenator< StringSet<TString, ConcatVirtual<TSpec> > > {
		typedef ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > > Type;
	};

	template < typename TString, typename TSpec >
	struct Concatenator< StringSet<TString, ConcatVirtual<TSpec> > const > {
		typedef ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > const > Type;
	};

	template < typename TString, typename TSpec >
	struct Concatenator< StringSet<TString, ConcatDirect<TSpec> > > {
		typedef TString Type;
	};

//////////////////////////////////////////////////////////////////////////////
// concat

	template <typename TString>
	inline typename Concatenator<TString>::Type & 
	concat(TString &string) {
		return string;
	}
//____________________________________________________________________________

	template <typename TString, typename TSpec>
	inline typename Concatenator< StringSet<TString, TSpec> >::Type
	concat(StringSet<TString, TSpec> &set) {
		return typename Concatenator< StringSet<TString, TSpec> >::Type (set);
	}

	template <typename TString, typename TSpec>
	inline typename Concatenator< StringSet<TString, TSpec> const>::Type
	concat(StringSet<TString, TSpec> const &set) {
		return typename Concatenator< StringSet<TString, TSpec> const>::Type (set);
	}
//____________________________________________________________________________

	template <typename TString, typename TSpec>
	inline typename Concatenator< StringSet<TString, ConcatDirect<TSpec> > >::Type & 
	concat(StringSet<TString, ConcatDirect<TSpec> > &set) {
		return set.concat;
	}

	template <typename TString, typename TSpec>
	inline typename Concatenator< StringSet<TString, ConcatDirect<TSpec> > const>::Type & 
	concat(StringSet<TString, ConcatDirect<TSpec> > const &set) {
		return set.concat;
	}


	//////////////////////////////////////////////////////////////////////////////
    // StringSet Iterator
    //////////////////////////////////////////////////////////////////////////////

    template < typename TStringSet, typename TSpec >
    class Iter< TStringSet, ConcatVirtual<TSpec> > : public ::std::iterator <
                                    ::std::bidirectional_iterator_tag,
                                    typename TStringSet::value_type,
                                    typename TStringSet::size_type,
                                    typename TStringSet::pointer,
                                    typename TStringSet::reference > {
	public:

        typedef Iter iterator;
		typedef typename TStringSet::obj_iterator obj_iterator;

        //////////////////////////////////////////////////////////////////////////////
		// STL compatible public iterator interface

		typedef ::std::bidirectional_iterator_tag		iterator_category;
		typedef typename TStringSet::const_reference	const_reference;
		typedef typename TStringSet::value_type		    value_type;
		typedef typename TStringSet::difference_type	difference_type;
		typedef typename TStringSet::pointer			pointer;
		typedef typename TStringSet::reference		    reference;
		
        TStringSet *host;
        unsigned objNo;
        obj_iterator _begin, _cur, _end;

		Iter() {}

        Iter(TStringSet &_host):
            host(&_host)
        {
            objNo = 0;
            _begin = _cur = begin(_host[objNo]);
            _end = end(_host[objNo]);
            testEnd();
        }

        Iter(TStringSet &_host, unsigned _objNo, difference_type _offset):
            host(&_host)
        {
            objNo = _objNo;
            if (objNo < length(_host)) {
				_begin = _cur = begin(_host[objNo]);
				_end = end(_host[objNo]);
                goFurther(_cur, _offset);
                testEnd();
            } else {
                _begin = _cur = obj_iterator();
                _end = obj_iterator();
            }
        }

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic
   		difference_type operator- (const iterator &I) const {
            return tell() - I.tell();
        }

        iterator operator- (difference_type delta) const {
            Pair< unsigned, difference_type > pos;
            host->convert(pos, tell() - delta);
            return iterator(*host, pos.i1, pos.i2);
        }

        iterator operator+ (difference_type delta) const {
            Pair< unsigned, difference_type > pos;
            host->convert(pos, tell() + delta);
            return iterator(*host, pos.i1, pos.i2);
        }

		inline iterator& operator++ () {
            ++_cur;
            testEnd();
			return *this;
        }

		inline iterator operator++ (int) {
			iterator before = *this;
            ++*this;
			return before;
		}

		inline iterator& operator-- () {
            testBegin();
            --_cur;
			return *this;
        }

		inline iterator operator-- (int) {
			iterator before = *this;
            --*this;
			return before;
		}

        inline reference operator* () {
            return *_cur;
        }

        inline const_reference operator* () const {
            return *_cur;
        }

		bool operator== (const iterator &I) const {
			SEQAN_ASSERT(host == I.host);
			return objNo == I.objNo && _cur == I._cur;
		}

		bool operator!= (const iterator &I) const {
			SEQAN_ASSERT(host == I.host);
			return objNo != I.objNo || _cur != I._cur;
		}

		bool operator< (const iterator &I) const {
			SEQAN_ASSERT(host == I.host);
			return objNo < I.objNo || (objNo == I.objNo && _cur < I._cur);
		}

        inline operator obj_iterator() {
            return _cur;
        }

    protected:

        inline void testBegin() {
            while (_cur == _begin && objNo != 0) {
                --objNo;
				_begin = _cur = host->_begin(objNo);
				_end = host->_end(objNo);
            }
        }

        inline void testEnd() {
            if (_cur == _end && objNo != length(*host))
                do {
                    if (++objNo == length(*host)) {
                        _begin = _cur = obj_iterator();
                        _end = obj_iterator();
                        break;
                    }
					_begin = _cur = begin((*host)[objNo]);
					_end = end((*host)[objNo]);
                } while (_cur == _end);
        }

        inline difference_type tell() {
			return host->limits[objNo] + distance(_begin, _cur);
        }
    };


//////////////////////////////////////////////////////////////////////////////
// Sequence - StringSet
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec = void>
struct GenerousStorage;

template<typename TSpec = void>
struct TightStorage;

// Default id holder string set
template<typename TSpec = GenerousStorage<> >
struct IdHolder;



//////////////////////////////////////////////////////////////////////////////

/**
.Class.StringSet:
..cat:Sequences
..summary:A container class for a set of strings.
..signature:StringSet<TString, TSpec>
..param.TString:The string type.
...type:Class.String
..param.TSpec:The specializing type for the StringSet.
...metafunction:Metafunction.Spec
...remarks:Possible values are IdHolder<TightStorage<> > or IdHolder<GenerousStorage<> >
...remarks:TightStorage is very space efficient whereas GenerousStorage provides fast access to the strings in the container via ids.
...default:$GenerousStorage$.
..include:sequence.h
*/

template<typename TString, typename TSpec>
class StringSet<TString, IdHolder<TightStorage<TSpec> > >
{
	public:
		typedef typename Id<StringSet>::Type TIdType;
        typedef String<TString*> TStrings;
		typedef String<TIdType> TIds;

		TStrings strings;
		TIds ids;

		template <typename TPos>
		inline typename Reference<StringSet>::Type
		operator [] (TPos pos)
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<StringSet const>::Type 
		operator [] (TPos pos) const
		{
			return value(*this, pos);
		}
};

//////////////////////////////////////////////////////////////////////////////


template<typename TString, typename TSpec>
class StringSet<TString, IdHolder<GenerousStorage<TSpec> > >
{
	public:
        typedef String<TString*> TStrings;
		typedef typename Size<StringSet>::Type TSize;
		TStrings strings;
		TSize counter;

		StringSet() : counter(0)
		{
			SEQAN_CHECKPOINT
		}

		template <typename TPos>
		inline typename Reference<StringSet>::Type
		operator [] (TPos pos)
		{
			SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<StringSet const>::Type 
		operator [] (TPos pos) const
		{
			return value(*this, pos);
		}
};


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

///.Function.appendValue.param.target.type:Class.StringSet

template<typename TString, typename TSpec>
inline void 
appendValue(StringSet<TString, IdHolder<GenerousStorage<TSpec> > >& me, 
			TString& obj) 
{
	SEQAN_CHECKPOINT
	appendValue(me.strings, &obj);
	++me.counter;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec>
inline void 
appendValue(StringSet<TString, IdHolder<TightStorage<TSpec> > >& me, 
			TString& obj) 
{
	SEQAN_CHECKPOINT
	appendValue(me.strings, &obj);
	appendValue(me.ids, length(me.strings) - 1);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addString:
..cat:Sequences
..summary:Adds a new string to the StringSet and retrieves an id.
..signature:addString(dest, str [,id])
..signature:addString(dest, source, id)
..param.dest:A StringSet.
...type:Class.StringSet
..param.source:A StringSet.
...type:Class.StringSet
..param.str:A new string.
...type:Metafunction.Value
..param.id:An associated id.
...type:Metafunction.Id
..returns:A new id
...type:Metafunction.Id
..see:Function.appendValue
..see:Function.getString
*/

template<typename TString, typename TSpec>
inline typename Id<StringSet<TString, IdHolder<TSpec> > >::Type 
addString(StringSet<TString, IdHolder<TSpec> >& me, 
		  TString& obj) 
{
	SEQAN_CHECKPOINT
	appendValue(me, obj);
	return length(me.strings) - 1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId>
inline typename Id<StringSet<TString, IdHolder<GenerousStorage<TSpec> > > >::Type 
addString(StringSet<TString, IdHolder<GenerousStorage<TSpec> > >& me, 
		  TString& obj,
		  TId id) 
{
	SEQAN_CHECKPOINT
	if (id >= (TId) length(me.strings)) fill(me.strings, id+1, (TString*) 0);
	if ((TString*) me.strings[id] == (TString*) 0) ++me.counter;
	me.strings[id] = &obj;
	return id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId>
inline typename Id<StringSet<TString, IdHolder<TightStorage<TSpec> > > >::Type 
addString(StringSet<TString, IdHolder<TightStorage<TSpec> > >& me, 
		  TString& obj,
		  TId id) 
{
	SEQAN_CHECKPOINT
	for(unsigned int i=0;i<length(me.ids);++i) {
		if ((TId) me.ids[i] == (TId) id) {
			me.strings[i]=&obj;
			return id;
		}
	}
	appendValue(me.strings, &obj);
	appendValue(me.ids, id);
	return id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId>
inline typename Id<StringSet<TString, IdHolder<GenerousStorage<TSpec> > > >::Type 
addString(StringSet<TString, IdHolder<GenerousStorage<TSpec> > >& dest, 
		  StringSet<TString, IdHolder<GenerousStorage<TSpec> > >& source,
		  TId id) 
{
	SEQAN_CHECKPOINT
	if (id >= length(dest.strings)) fill(dest.strings, id+1, (TString*) 0);
	if (dest.strings[id] == (TString*) 0) ++dest.counter;
	dest.strings[id] = source.strings[id];
	return id;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId>
inline typename Id<StringSet<TString, IdHolder<TightStorage<TSpec> > > >::Type 
addString(StringSet<TString, IdHolder<TightStorage<TSpec> > >& dest, 
		  StringSet<TString, IdHolder<TightStorage<TSpec> > >& source,
		  TId id) 
{
	SEQAN_CHECKPOINT
	for(unsigned int i=0;i<length(source.ids);++i) {
		if (source.ids[i] == id) {
			appendValue(dest.strings, source.strings[i]);
			appendValue(dest.ids, source.ids[i]);
			return id;
		}
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getString:
..cat:Sequences
..summary:Retrieves a string from the StringSet given an id.
..signature:getString(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:A reference to a string.
..see:Function.addString
*/

template<typename TString, typename TSpec, typename TId>
inline typename Reference<StringSet<TString, IdHolder<GenerousStorage<TSpec> > > >::Type
getString(StringSet<TString, IdHolder<GenerousStorage<TSpec> > >& me, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	return value(me, id);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId>
inline typename Reference<StringSet<TString, IdHolder<TightStorage<TSpec> > > >::Type
getString(StringSet<TString, IdHolder<TightStorage<TSpec> > >&me, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	for(unsigned int i=0;i<length(me.strings);++i) {
		if ((TId) me.ids[i] == (TId) id) {
			return value(me, i);
		}
	}
	static TString tmp = "";
	return tmp;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.removeString:
..cat:Sequences
..summary:Removes a string from the StringSet given an id.
..signature:removeString(me, id)
..param.me:A StringSet.
...type:Class.StringSet
..param.id:An id.
...type:Metafunction.Id
..returns:void
..see:Function.addString
*/

template<typename TString, typename TSpec, typename TId>
inline void
removeString(StringSet<TString, IdHolder<GenerousStorage<TSpec> > >& me, 
			 TId const id) 
{
	SEQAN_CHECKPOINT
	if (me.strings[id] != (TString*) 0) --me.counter;
	me.strings[id] = 0;
	while((!empty(me.strings)) &&
		(me.strings[length(me.strings) - 1] == 0))
	{
		resize(me.strings, length(me.strings) - 1);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId>
inline void
removeString(StringSet<TString, IdHolder<TightStorage<TSpec> > >& me, 
			 TId const id) 
{
	SEQAN_CHECKPOINT
	for(unsigned int i=0;i<length(me.strings);++i) {
		if (me.ids[i] == id) {
			for(unsigned int j=i+1;j<length(me.strings);++j) {
				me.strings[j-1] = me.strings[j];
				me.ids[j-1] = me.ids[j];
			}
			resize(me.strings, length(me.strings)-1);
			resize(me.ids, length(me.ids)-1);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

///.Function.clear.param.object.type:Class.StringSet

template<typename TString, typename TSpec>
inline void
clear(StringSet< TString, IdHolder<GenerousStorage<TSpec> > > & me) 
{
	SEQAN_CHECKPOINT
	clear(me.strings);
	me.counter = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec>
inline void
clear(StringSet<TString, IdHolder<TightStorage<TSpec> > >& me) 
{
	SEQAN_CHECKPOINT
	clear(me.strings);
	clear(me.ids);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.StringSet

template<typename TString, typename TSpec>
inline typename Size<StringSet<TString, IdHolder<TightStorage<TSpec> > > >::Type 
length(StringSet<TString, IdHolder<TightStorage<TSpec> > > const &me) 
{
	return length(me.strings);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TString, typename TSpec>
inline typename Size<StringSet<TString, IdHolder<TightStorage<TSpec> > > >::Type 
length(StringSet<TString, IdHolder<TightStorage<TSpec> > > &me) 
{
	SEQAN_CHECKPOINT
	return length(me.strings);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec>
inline typename Size<StringSet<TString, IdHolder<GenerousStorage<TSpec> > > >::Type 
length(StringSet<TString, IdHolder<GenerousStorage<TSpec> > > const &me) 
{
	return me.counter;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec>
inline typename Size<StringSet<TString, IdHolder<GenerousStorage<TSpec> > > >::Type 
length(StringSet<TString, IdHolder<GenerousStorage<TSpec> > > &me) 
{
	SEQAN_CHECKPOINT
	return me.counter;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.value.param.object.type:Class.StringSet

template<typename TString, typename TSpec, typename TPos>
inline typename Reference<StringSet< TString, IdHolder<TSpec> > >::Type
value(StringSet< TString, IdHolder<TSpec> > & me, 
	  TPos pos)
{
	SEQAN_CHECKPOINT
	if (me.strings[pos]!=0) return *me.strings[pos];
	else {
		static TString tmp = "";
		return tmp;
	}
}

//////////////////////////////////////////////////////////////////////////////

template < typename TString, typename TSpec, typename TPos >
inline typename Reference< StringSet< TString, IdHolder<TSpec> > const >::Type
value(StringSet< TString, IdHolder<TSpec> > const & me, 
	  TPos pos)
{
	if (me.strings[pos]!=0) return *me.strings[pos];
	else {
		static TString tmp = "";
		return tmp;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.subset:
..cat:Sequences
..summary:Creates a subset of a given StringSet.
..signature:subset(source, dest, id_array [, len])
..param.source:In-parameter:The source StringSet.
...type:Class.StringSet
..param.dest:Out-parameter:The destination StringSet (the subset).
...type:Class.StringSet
..param.id_array:In-parameter:An array of ids. Each id corresponds to a sequence that is supposed to be in the subset.
..param.len:In-parameter:Optional length of the id array.
...remarks:If len is not defined the length function must be valid for this array or string type.
..returns:void
*/

template <typename TString, typename TSpec, typename TIds, typename TLength>
inline void
subset(StringSet<TString, IdHolder<GenerousStorage<TSpec> > >& source,
	   StringSet<TString, IdHolder<GenerousStorage<TSpec> > >& dest,
	   TIds ids,
	   TLength len)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, IdHolder<GenerousStorage<TSpec> > > TStringSet;
	typedef typename Id<TStringSet>::Type TId;

	clear(dest);
	dest.counter = len;
	fill(dest.strings, length(source.strings), (TString*) 0);
	for(unsigned int i=0;i<len;++i) dest.strings[ids[i]] = source.strings[ids[i]];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TString, typename TSpec, typename TIds, typename TLength>
inline void
subset(StringSet<TString, IdHolder<TightStorage<TSpec> > >& source,
	   StringSet<TString, IdHolder<TightStorage<TSpec> > >& dest,
	   TIds ids,
	   TLength len)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, IdHolder<TightStorage<TSpec> > > TStringSet;
	typedef typename Id<TStringSet>::Type TId;

	clear(dest);
	TLength upperBound = length(source.ids);
	for(unsigned int i=0;i<len;++i) {
		TId id = ids[i];
		if ((upperBound > id) &&
			(source.ids[id] == id)) {
				appendValue(dest.strings, source.strings[id]);
				appendValue(dest.ids, id);
		} else {
			typedef String<TId> TIdString;
			typedef typename Iterator<TIdString>::Type TIter;
			TIter it = begin(source.ids);
			for(;!atEnd(it);goNext(it)) {
				if (*it == id) {
					appendValue(dest.strings, source.strings[position(it)]);
					appendValue(dest.ids, id);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TString, typename TSpec, typename TIds>
inline void
subset(StringSet<TString, IdHolder<TSpec> >& source,
		StringSet<TString, IdHolder<TSpec> >& dest,
		TIds ids)
{
	SEQAN_CHECKPOINT
	subset(source,dest,ids,length(ids));
}

//////////////////////////////////////////////////////////////////////////////


}

#endif
