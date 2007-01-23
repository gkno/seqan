#ifndef SEQAN_HEADER_BASIC_AGGREGATES_H
#define SEQAN_HEADER_BASIC_AGGREGATES_H

namespace SEQAN_NAMESPACE_MAIN
{

//____________________________________________________________________________

    struct _Compressed;
	typedef Tag<_Compressed> Compressed;

/**
.Class.Pair:
..cat:Aggregates
..summary:Stores two arbitrary objects.
..signature:Pair<T1, T2[, Compression]>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..param.Compression:If $Compressed$, the pair is stored in a more space efficient way (useful for external storage).
...note:When compression is enabled, referring to members is not allowed.
...default:$void$, no compression (faster access).
.Memfunc.Pair#Pair:
..class:Class.Pair
..summary:Constructor
..signature:Pair<T1, T2> ()
..signature:Pair<T1, T2> (pair)
..signature:Pair<T1, T2> (i1, i2)
..param.pair:Other Pair object. (copy constructor)
..param.i1:T1 object.
..param.i2:T2 object.
.Memvar.Pair#i1:
..class:Class.Pair
..summary:T1 object
.Memvar.Pair#i2:
..class:Class.Pair
..summary:T2 object
*/

	// standard storage 
	template <typename _T1, typename _T2 = _T1, typename _Compression = void>
    struct Pair {
        typedef _T1 T1;
        typedef _T2 T2;
	    _T1 i1;
	    _T2 i2;
		Pair() {}
		Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		Pair(_T1 const &_i1, _T2 const &_i2): i1(_i1), i2(_i2) {}
    };

#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
    template <typename _T1, typename _T2>
    struct Pair<_T1, _T2, Compressed> {
        typedef _T1 T1;
        typedef _T2 T2;
	    _T1 i1;
	    _T2 i2;
		Pair() {}
		Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		Pair(_T1 const &_i1, _T2 const &_i2): i1(_i1), i2(_i2) {}
    };
    #pragma pack(pop)
#else
    template <typename _T1, typename _T2>
    struct Pair<_T1, _T2, Compressed> {
        typedef _T1 T1;
        typedef _T2 T2;
        _T1 i1;
        _T2 i2;
		Pair() {}
		Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		Pair(_T1 const &_i1, _T2 const &_i2): i1(_i1), i2(_i2) {}
    }__attribute__((packed));
#endif

    // do not pack two identical types
    template <typename _T>
    struct Pair<_T, _T, Compressed> {
        typedef _T T1;
        typedef _T T2;
        _T i1;
        _T i2;
		Pair() {}
		Pair(Pair const &_p): i1(_p.i1), i2(_p.i2) {}
		Pair(_T const &_i1): i1(_i1), i2() {}
		Pair(_T const &_i1, _T const &_i2): i1(_i1), i2(_i2) {}
    };

    template <typename _T1, typename _T2, typename _Compression>
	std::ostream& operator<<(std::ostream &out, Pair<_T1,_T2,_Compression> const &p) {
		out << "< " << p.i1 << " , " << p.i2 << " >";
		return out;
	}


//____________________________________________________________________________

/**
.Class.Triple:
..cat:Aggregates
..summary:Stores three arbitrary objects.
..signature:Triple<T1, T2, T3[, Compression]>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..param.T3:The type of the third object.
..param.Compression:If $Compressed$, the triple is stored in a more space efficient way (useful for external storage).
...note:When compression is enabled, referring to members is not allowed.
...default:$void$, no compression (faster access).
.Memfunc.Triple#Triple:
..class:Class.Triple
..summary:Constructor
..signature:Triple<T1, T2, T3> ()
..signature:Triple<T1, T2, T3> (triple)
..signature:Triple<T1, T2, T3> (i1, i2, i3)
..param.triple:Other Triple object. (copy constructor)
..param.i1:T1 object.
..param.i2:T2 object.
..param.i3:T3 object.
.Memvar.Triple#i1:
..class:Class.Triple
..summary:T1 object
.Memvar.Triple#i2:
..class:Class.Triple
..summary:T2 object
.Memvar.Triple#i3:
..class:Class.Triple
..summary:T3 object
*/
    template <typename _T1, typename _T2 = _T1, typename _T3 = _T1, typename _Compression = void>
    struct Triple {
        typedef _T1 T1;
        typedef _T2 T2;
        typedef _T3 T3;
        _T1 i1;
        _T2 i2;
        _T3 i3;
		Triple() {}
		Triple(Triple const &_p): i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
		Triple(_T1 const &_i1, _T2 const &_i2, _T3 const &_i3): i1(_i1), i2(_i2), i3(_i3) {}
    };


#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
    template <typename _T1, typename _T2, typename _T3>
    struct Triple<_T1, _T2, _T3, Compressed> {
        typedef _T1 T1;
        typedef _T2 T2;
        typedef _T3 T3;
        _T1 i1;
        _T2 i2;
        _T3 i3;
		Triple() {}
		Triple(Triple const &_p): i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
		Triple(_T1 const &_i1, _T2 const &_i2, _T3 const &_i3): i1(_i1), i2(_i2), i3(_i3) {}
    };
    #pragma pack(pop)
#else
    template <typename _T1, typename _T2, typename _T3>
    struct Triple<_T1, _T2, _T3, Compressed> {
        typedef _T1 T1;
        typedef _T2 T2;
        typedef _T3 T3;
        _T1 i1;
        _T2 i2;
        _T3 i3;
		Triple() {}
		Triple(Triple const &_p): i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
		Triple(_T1 const &_i1): i1(_i1) {}
		Triple(_T1 const &_i1, _T2 const &_i2): i1(_i1), i2(_i2) {}
		Triple(_T1 const &_i1, _T2 const &_i2, _T3 const &_i3): i1(_i1), i2(_i2), i3(_i3) {}
    }__attribute__((packed));
#endif

    // do not pack three identical types
    template <typename _T>
    struct Triple<_T, _T, _T, Compressed> {
        typedef _T T1;
        typedef _T T2;
        typedef _T T3;
        _T i1;
        _T i2;
        _T i3;
		Triple() {}
		Triple(Triple const &_p): i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
		Triple(_T const &_i1, _T const &_i2, _T const &_i3): i1(_i1), i2(_i2), i3(_i3) {}
    };

    template <typename _T1, typename _T2, typename _T3, typename _Compression>
	std::ostream& operator<<(std::ostream &out, Triple<_T1,_T2,_T3,_Compression> const &t) {
		out << "< " << t.i1 << " , " << t.i2 << " , " << t.i3 << " >";
		return out;
	}


//____________________________________________________________________________

/**
.Class.Tuple:
..cat:Aggregates
..summary:A plain fixed-length string.
..signature:Tuple<T, size[, compress]>
..param.T:The value type, that is the type of characters stored in the tuple.
..param.size:The size/length of the tuple.
...remarks:In contrast to @Class.String@ the length of Tuple is fixed.
..param.compress:Enable/Disable compression.
..param.compress:If $void$, no compression is used.
..param.compress:If $Compressed$, the characters are stored as a bit sequence in an ordinal type (char, ..., __int64)
...remarks:Only useful for small alphabets and small tuple sizes (|Sigma|^size <= 2^64) as for DNA or protein m-grams)
...default:void.
..see:Spec.Sampler
*/
    template <typename _T, int _size, typename _Compression = void>
    struct Tuple {
        typedef _T T;
        enum { size = _size };
        _T i[_size];
//		friend std::ostream& operator<<(std::ostream&, const Tuple&);

        _T& operator[](const int k) {
            SEQAN_ASSERT(k >= 0 && k < size);
            return i[k];
        }
        const _T& operator[](const int k) const {
            SEQAN_ASSERT(k >= 0 && k < size);
            return i[k];
        }
		_T* operator&() { return i; }
		const _T* operator&() const { return i; }

		// has to be inline because elements (like this tuple) of packed structs can't be arguments
		template <typename _S>
		inline _S const assignAt(int k, _S const source) {
			return i[k] = source;
		}
    };


    template < int _size >
    struct BitVector {
        typedef typename BitVector<_size + 1>::Type Type;
    };

    template <> struct BitVector<8> { typedef unsigned char Type; };
    template <> struct BitVector<16> { typedef unsigned short Type; };
    template <> struct BitVector<32> { typedef unsigned int Type; };
    template <> struct BitVector<64> { typedef __int64 Type; };

    template <typename _T, int _size>
    struct Tuple<_T, _size, Compressed> {
        typedef _T T;
        enum { size = _size };
//        enum { bitSize = Log2<ValueSize<_T>::VALUE>::VALUE };
        enum { bitSize = BitsPerValue<_T>::VALUE };
        enum { bitMask = (1 << bitSize) - 1 };
        enum { mask = (1 << (size * bitSize)) - 1 };
        typedef typename BitVector< bitSize * size >::Type CT;
        
        CT i;

        const _T operator[](const int k) const {
            SEQAN_ASSERT(k >= 0 && k < size);
            return (i >> (size - 1 - k) * bitSize) & bitMask;
        }
        inline CT operator<<=(int shift) {
            return i = (i << (shift * bitSize)) & mask;
        }
        inline CT operator<<(int shift) const {
            return (i << (shift * bitSize)) & mask;
        }
        inline CT operator>>(int shift) const {
            return i >> (shift * bitSize);
        }
        template <typename T>
        inline void operator|=(T const &t) {
            i |= t;
        }
        template <typename T, typename TSpec>
        inline void operator|=(SimpleType<T, TSpec> const &t) {
            i |= t.value;
        }
		CT* operator&() { return &i; }
		const CT* operator&() const { return &i; }

		// has to be inline because elements (like this tuple) of packed structs can't be arguments
		template <typename _S>
		inline _S const assignAt(int k, _S const source) {
			typedef Tuple<_T, _size, Compressed> Tup;
			typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * bitSize);
			i = (i & ~mask) | ((CT)source << ((_size - 1 - k) * bitSize));
			return source;
		}
    };

    template <typename _T, int _size, typename _Compression>
	inline int length(Tuple<_T, _size, _Compression> const &me) { return _size; }

    template <typename _T, int _size, typename _S>
    inline _S const assignAt(Tuple<_T, _size> &me, int k, _S const source) {
        return me.i[k] = source;
    }

    template <typename _T, int _size, typename _S>
    inline _S const assignAt(Tuple<_T, _size, Compressed> &me, int k, _S const source) {
        typedef Tuple<_T, _size, Compressed> Tup;
        typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * me.bitSize);
        me.i = (me.i & ~mask) | source << ((_size - 1 - k) * me.bitSize);
        return source;
    }

    template <typename _T, typename _S, typename _Spec, int _size>
    inline SimpleType<_S, _Spec> const & assignAt(Tuple<_T, _size, Compressed> &me, int k, SimpleType<_S, _Spec> const &source) {
        typedef Tuple<_T, _size, Compressed> Tup;
        typename Tup::CT mask = Tup::bitMask << ((_size - 1 - k) * me.bitSize);
        me.i = (me.i & ~mask) | source.value << ((_size - 1 - k) * me.bitSize);
        return source;
    }

    template <typename _T, int _size, typename _Compression>
	std::ostream& operator<<(std::ostream& out, Tuple<_T,_size,_Compression> const &a) {
		out << "[";
		if (a.size > 0)
			out << a[0];
		for(int j = 1; j < a.size; ++j)
			out << " " << a[j];
		out << "]";
		return out;
	}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
