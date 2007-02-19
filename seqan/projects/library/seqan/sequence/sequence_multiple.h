/*
 *  string_set.h
 *  genindex
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

	template < typename TString, typename TSpec = ConcatVirtual<> >
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
	inline typename StringSetLimits<TStringSet>::Type & 
	stringSetLimits(TStringSet &stringSet) {
		return Nothing();
	}

    template <typename TStringSet>
	inline typename StringSetLimits<TStringSet const>::Type & 
	stringSetLimits(TStringSet const &stringSet) {
		return Nothing();
	}

    template <typename TString, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSpec> >::Type & 
	stringSetLimits(StringSet<TString, TSpec> &stringSet) {
		return stringSet.limits;
	}

    template <typename TString, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSpec> const>::Type & 
	stringSetLimits(StringSet<TString, TSpec> const &stringSet) {
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
		typedef typename Iterator<TLimitsString const>::Type TIter;
        TIter begin = limits.begin();
        TIter upper = upper_bound(begin, limits.end(), pos) - 1;
        return distance(begin, upper);
	}

	template <typename TPos, typename TLimitsString>
	inline TPos getSeqOffset(TPos const &pos, TLimitsString const &limits) {
		typedef typename Iterator<TLimitsString const>::Type TIter;
        TIter begin = limits.begin();
        TIter upper = upper_bound(begin, limits.end(), pos) - 1;
        return pos - *upper;
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
		typedef typename Iterator<String<TSize> const>::Type TIter;
        TIter begin = limits.begin();
        TIter upper = upper_bound(begin, limits.end(), pos) - 1;
        result.i1 = distance(begin, upper);
        result.i2 = pos - *upper;
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

        typedef String<TString>								TStrings;
		typedef typename StringSetLimits<StringSet>::Type	TLimits;
		typedef typename Iterator<TLimits const>::Type		TLimitIterator;

		typedef typename Value<TString>::Type			Type;
        typedef typename Size<TString>::Type			SizeType;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator types

        friend class StringSetIterator<StringSet>;

		typedef StringSetIterator<StringSet>	StringSetConstIterator;
		typedef StringSetIterator<StringSet>	StringSetIterator;

		typedef StringSetIterator		iterator;
		typedef StringSetConstIterator	const_iterator;
        typedef Type					value_type;
		typedef Type&					reference;
		typedef Type const &			const_reference;
		typedef Type*					pointer;
		typedef Type const *			const_pointer;
		typedef SizeType				size_type;
		typedef SizeType				difference_type;

		// used by StringSetIterator
        typedef typename Iterator<TString>::Type		obj_iterator;
        typedef typename Iterator<TString const>::Type	const_obj_iterator;

		TStrings	strings;
        TLimits		limits;

		typename Concatenator<StringSet>::Type	concat;

		StringSet(): concat(*this) {};

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
		typedef typename StringSetLimits<StringSet>::Type	TLimits;
        typedef typename Iterator<TLimits const>::Type		TLimitIterator;

		typedef typename Value<TString>::Type			Type;
        typedef typename Size<TString>::Type			SizeType;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator types

        friend class StringSetIterator<StringSet>;

		typedef StringSetIterator<StringSet>	StringSetConstIterator;
		typedef StringSetIterator<StringSet>	StringSetIterator;

		typedef StringSetIterator		iterator;
		typedef StringSetConstIterator	const_iterator;
        typedef Type					value_type;
		typedef Type&					reference;
		typedef Type const &			const_reference;
		typedef Type*					pointer;
		typedef Type const *			const_pointer;
		typedef SizeType				size_type;
		typedef SizeType				difference_type;

		// used by StringSetIterator
        typedef typename Iterator<TString>::Type		obj_iterator;
        typedef typename Iterator<TString const>::Type	const_obj_iterator;

        TLimits		limits;

		typename Concatenator<StringSet>::Type	concat;

		StringSet() {
			resize(limits, 1);
			limits[0] = 0;
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

	template < typename TString, typename TSpec >
    inline typename Size<TString>::Type lengthSum(StringSet< TString, TSpec > const &me) {
        return back(me.limits);
    }

//____________________________________________________________________________

	template < typename TString, typename TSpec >
    inline void appendValue(StringSet< TString, ConcatVirtual<TSpec> > &me, TString &obj) {
        appendValue(me.strings, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }


	template < typename TString >
    inline void appendValue(StringSet< TString, ConcatDirect<void> > &me, TString &obj) {
        append(me.concat, obj);
        appendValue(me.limits, lengthSum(me) + length(obj));
    }

    template < typename TString, typename TLimiter >
    inline void appendValue(StringSet< TString, ConcatDirect<TLimiter> > &me, TString &obj) {
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

/*
        inline void convert(Pair< unsigned, difference_type > &dst, difference_type src) const
        {
            TLimitIterator begin = limits.begin();
            TLimitIterator upper = upper_bound(begin, limits.end(), src) - 1;
            dst.i1 = distance(begin, upper);
            dst.i2 = src - *upper;
        }

        inline iterator begin() {
            return iterator(*this);
        }

        inline const_iterator begin() const {
            return const_iterator(*this);
        }

		inline obj_iterator begin(unsigned i) {
            return ::seqan::begin(TString());
		}

		//inline const_obj_iterator begin(unsigned i) const {
  //          return ::seqan::begin(TString());
		//}

		inline iterator end() {
            return iterator(*this, objs.size(), 0);
        }

        inline const_iterator end() const {
            return const_iterator(*this, objs.size(), 0);
        }

		inline obj_iterator end(unsigned i) {
			return ::seqan::end(TString());
		}

		//inline const_obj_iterator end(unsigned i) const {
		//	return ::seqan::end(TString());
		//}
        
    };

*/

    //////////////////////////////////////////////////////////////////////////////
    // StringSet to String converter
    //////////////////////////////////////////////////////////////////////////////

    template <typename TStringSet>
	struct ConcatenatorNto1 {
		TStringSet &set;
		ConcatenatorNto1 (TStringSet &_set): set(_set) {}
	};

	template < typename TString, typename TSpec >
	struct Concatenator< StringSet<TString, ConcatVirtual<TSpec> > > {
		typedef ConcatenatorNto1< StringSet<TString, ConcatVirtual<TSpec> > > Type;
	};

	template < typename TString, typename TSpec >
	struct Concatenator< StringSet<TString, ConcatDirect<TSpec> > > {
		typedef TString Type;
	};


	template <typename TString>
	inline typename Concatenator<TString>::Type & concat(TString &string) {
		return string;
	}

	template <typename TString, typename TSpec>
	inline typename Concatenator< StringSet<TString, TSpec> >::Type & concat(StringSet<TString, TSpec> &set) {
		return set.concat;
	}

	template <typename TString, typename TSpec>
	inline typename Concatenator< StringSet<TString, TSpec> const>::Type & concat(StringSet<TString, TSpec> const &set) {
		return set.concat;
	}


//____________________________________________________________________________

	template < typename TStringSet >
	inline Iter<TStringSet, void> begin(ConcatenatorNto1<TStringSet> &Nto1) {
		return Iter<TStringSet, void> (Nto1.set);
	}

	template < typename TStringSet >
	inline Iter<TStringSet, void> end(ConcatenatorNto1<TStringSet> &Nto1) {
		return Iter<TStringSet, void> (Nto1.set, length(Nto1.set), 0);
	}

//____________________________________________________________________________

	template < typename TString, typename TSSetSpec, typename TSpec >
	inline Iter< StringSet<TString, TSSetSpec>, ConcatVirtual<TSpec> >
	begin(StringSet<TString, TSSetSpec> &set) {
		return Iter< StringSet<TString, TSSetSpec>, ConcatVirtual<TSpec> > (value(set.holder));
	}

	template < typename TString, typename TSSetSpec, typename TSpec >
	inline Iter< StringSet<TString, TSSetSpec>, ConcatVirtual<TSpec> >
	end(StringSet<TString, TSSetSpec> &set) {
		return Iter< StringSet<TString, TSSetSpec>, ConcatVirtual<TSpec> > (value(set.holder), length(value(set.holder)), 0);
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
		// public iterator interface

		// in fact this is also a random iterator, but you better use
		// the VectorIterator class for *real* random access
        typedef ::std::bidirectional_iterator_tag		iterator_category;
		typedef typename TStringSet::const_reference	const_reference;
		typedef typename TStringSet::value_type		    value_type;
		typedef typename TStringSet::difference_type	difference_type;
		typedef typename TStringSet::pointer			pointer;
		typedef typename TStringSet::reference		    reference;
		
        TStringSet *host;
        unsigned objNo;
        obj_iterator begin, cur, end;

        Iter(TStringSet &_host):
            host(&_host)
        {
            objNo = 0;
            begin = cur = _host->begin(objNo);
            end = _host->end(objNo);
            testEnd();
        }

        Iter(TStringSet &_host, unsigned _objNo, difference_type _offset):
            host(&_host)
        {
            objNo = _objNo;
            if (objNo < _host->element_count()) {
				begin = cur = _host->begin(objNo);
				end = _host->end(objNo);
                goFurther(cur, _offset);
                testEnd();
            } else {
                begin = cur = NULL;
                end = NULL;
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
            ++cur;
            testEnd();
        }

		inline iterator operator++ (int) {
			iterator before = *this;
            ++*this;
			return before;
		}

		inline iterator& operator-- () {
            testBegin();
            --cur;
        }

		inline iterator operator-- (int) {
			iterator before = *this;
            --*this;
			return before;
		}

        //inline reference operator* () {
        //    return *cur;
        //}

        inline const_reference operator* () const {
            return *cur;
        }

		bool operator== (const iterator &I) const {
			SEQAN_ASSERT(host == I.host);
			return objNo == I.objNo && cur == I.cur;
		}

		bool operator!= (const iterator &I) const {
			SEQAN_ASSERT(host == I.host);
			return objNo != I.objNo || cur != I.cur;
		}

		bool operator< (const iterator &I) const {
			SEQAN_ASSERT(host == I.host);
			return objNo < I.objNo || (objNo == I.objNo && cur < I.cur);
		}

        inline operator obj_iterator() {
            return cur;
        }

    protected:

        inline void testBegin() {
            while (cur == begin && objNo != 0) {
                --objNo;
				begin = cur = host->begin(objNo);
				end = host->end(objNo);
            }
        }

        inline void testEnd() {
            if (cur == end && objNo != host->element_count())
                do {
                    if (++objNo == host->element_count()) {
                        begin = cur = NULL;
                        end = NULL;
                        break;
                    }
					begin = cur = host->begin(objNo);
					end = host->end(objNo);
                } while (cur == end);
        }

        inline difference_type tell() {
			return host->limits[objNo] + distance(begin, cur);
        }
    };

/*
    //////////////////////////////////////////////////////////////////////////////
    // Global Interface
    //////////////////////////////////////////////////////////////////////////////

    template < typename TString >
    struct Value< StringSet< TString > > {
        typedef typename Value< TString >::Type Type;
    };

    template < typename TString >
    struct Size< StringSet< TString > > {
        typedef typename Size< TString >::Type Type;
    };

    template < typename TString >
    struct Iterator< StringSet< TString > > {
        typedef StringSetIterator< StringSet< TString> > Type;
    };


    template < typename TString >
    inline typename Size< StringSet<TString> >::Type
    size(StringSet<TString> &me)
    {
        return me.size();
    }

    template < typename TString >
    inline typename Size< StringSet<TString> >::Type
    length(StringSet<TString> &me)
    {
        return me.size();
    }

    template < typename TString >
	inline void add_element(StringSet<TString> &me, TString &obj) {
		return me.add(obj);
	}

    template < typename TString >
	inline TString & element_at(StringSet<TString> &me, unsigned i) {
		return me.element_at(i);
	}

    template < typename TString >
	inline TString const & element_at(StringSet<TString> const &me, unsigned i) {
		return me.element_at(i);
	}

    template < typename TString >
    inline unsigned
	element_count(StringSet<TString> &me)
    {
        return me.element_count();
    }

    template < typename TString >
    inline typename Iterator< StringSet<TString> >::Type
    begin(StringSet<TString> &me) {
        return me.begin();
    }

    template < typename TString >
    inline typename Iterator< StringSet<TString> const >::Type
    begin(StringSet<TString> const &me) {
        return me.begin();
    }

    template < typename TString >
    inline typename Iterator< StringSet<TString> >::Type
    end(StringSet<TString> &me) {
        return me.end();
    }

    template < typename TString >
    inline typename Iterator< StringSet<TString> const >::Type
    end(StringSet<TString> const &me) {
        return me.end();
    }

    template < typename TString, typename TPos >
    inline typename Reference< StringSet<TString> >::Type 
    at(StringSet<TString> &me, TPos pos)
    {
	    return me[pos];
    }

    template < typename TString >
	std::ostream& operator<<(std::ostream &out, StringSet<TString> &p) {
        typename Iterator< StringSet<TString> const >::Type _cur = begin(p), _end = end(p);
        while (_cur != _end) {
		    out << *_cur<< "\n";
            ++_cur;
        }
		return out;
	}


    //////////////////////////////////////////////////////////////////////////////
    // StringSet Buffer Handler
    //////////////////////////////////////////////////////////////////////////////

	template < typename TValue,
               typename TConfig,
			   typename TSpec >
	struct BufferHandler< Pipe< StringSet< String<TValue, External<TConfig> >, TSpec>, Source<ContainerSpec> > >
    {
        typedef TValue																	Type;
        typedef typename Size< StringSet< String<TValue, External<TConfig> > > >::Type	SizeType;
        typedef SimpleBuffer<TValue, SizeType>											SimpleBuffer;

        typedef Pipe< StringSet< String<TValue, External<TConfig> >, TSpec>, Source<ContainerSpec> >	    Pipe;
        typedef typename StringSet< String<TValue, External<TConfig> >, TSpec>::StringSetConstIterator	Iterator;

		Pipe			&pipe;
        int             pageNo;
		unsigned		bufferSize;
        SimpleBuffer	buffer;
		Iterator		iter;

		BufferHandler(Pipe &_pipe, unsigned requestedSize):
			pipe(_pipe),
			bufferSize(requestedSize) {}

        inline SimpleBuffer& begin() {
			allocPage(buffer, Min(bufferSize, length(pipe.in)));
			iter = begin(pipe.in);
			for(unsigned i = 0, len = length(buffer); i < len; ++i) {
				buffer[i] = *iter;
				++iter;
			}
			return buffer;
        }

        inline SimpleBuffer& next() {
			unsigned len = Min(length(buffer), end(pipe.in) - iter);
			for(unsigned i = 0; i < len; ++i) {
				buffer[i] = *iter;
				++iter;
			}
			resize(buffer, len);
			return buffer;
        }

        inline void process() {}

        inline void end() {
			cancel();
		}

        inline void cancel() {
			freePage(buffer);
		}
    };

    template < typename TValue, typename TSpec, typename TConfig >
    struct Value< BufferHandler< String<TValue, External<TConfig> >, TSpec > > {
        typedef SimpleBuffer< TValue > Type;
    };
/*
    // sequence -> external string
    template < typename TValue,
               typename TConfig,
               typename TSource >
    inline void assign(String<TValue, External<TConfig> > &target, TSource const &source) {

        typedef typename Iterator<TSource const>::Type                          ISource;
        typedef typename String<TValue, External<TConfig> >::VectorFwdIterator  ITarget;

        resize(target, length(source));

        ISource it_source       = begin(source);
        ITarget it_target       = begin(target);
		ITarget it_target_end   = end(target);
		while (it_target != it_target_end)
		{
			*it_target = *it_source;
			++it_target;
			++it_source;
		}
    }*/

}

#endif
