/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007-2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ============================================================================
  Implementation of the LongWord class.  It provides an almost
  transparent interface to arbitrarily long words.  In contrast to
  std::bitset, it provides a specialization that allows to set the
  number of bits in the virtual word at runtime.
 ==========================================================================*/

// TODO(holtgrew): This should probably not be called LongWord, maybe BitSet/BitVector/BitString is more appropriate?
// TODO(holtgrew): Optimized implementation of operator>>(TWord, 1) and operator<<(TWord, 1).
// TODO(holtgrew): Optimized implementation of value(TWord, index), operator[] always seems to call the proxy version!

#ifndef SEQAN_MISC_MISC_LONG_WORD_H_
#define SEQAN_MISC_MISC_LONG_WORD_H_

#include <bitset>

namespace seqan {

// TODO(holtgrew): Document this.
template <typename TSpec>
struct LongWord;


template <typename TWord, typename TSpec>
struct LongWordBitProxy;


struct _NativeWidth;
typedef Tag<_NativeWidth> NativeWidth;


struct _NativeWideWidth;
typedef Tag<_NativeWideWidth> NativeWideWidth;


struct _DynamicWidth;
typedef Tag<_DynamicWidth> DynamicWidth;


template <unsigned LENGTH>
struct StaticWidth;


template <typename TWord>
struct LongWordBitProxy<TWord, NativeWidth> {
    TWord &_word;
    unsigned _bitIndex;

    LongWordBitProxy(TWord & word, unsigned bitIndex)
            : _word(word), _bitIndex(bitIndex) {
        SEQAN_CHECKPOINT;
    }

    LongWordBitProxy & operator=(unsigned x) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(x, sizeof(unsigned) * 8);
        _word._value = (_word._value ^ (1 << _bitIndex)) | (x << _bitIndex);
        return *this;
    }

    operator unsigned() const {
        return (_word & (1u << _bitIndex)) >> _bitIndex;
    }
};


template <>
struct LongWord<NativeWidth> {
    typedef LongWord<NativeWidth> TWord;
    unsigned _value;

    LongWord() {
        SEQAN_CHECKPOINT;
        _value = 0;
    }

    LongWord(unsigned const & value) : _value(value) {
        SEQAN_CHECKPOINT;
    }

    // Conversion operators must be member functions.
    operator unsigned() const {
        SEQAN_CHECKPOINT;
        return _value;
    }

    unsigned operator[](unsigned index) const {
        SEQAN_ASSERT_LT(index, sizeof(unsigned) * 8);
        return (_value & (1u << index)) >> index;
    }

    LongWordBitProxy<TWord, NativeWidth> operator[](unsigned index) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(index, sizeof(unsigned) * 8);
        return LongWordBitProxy<TWord, NativeWidth>(*this, index);
    }
};


inline size_t length(LongWord<NativeWidth> const & x) {
    SEQAN_CHECKPOINT;
    return sizeof(unsigned) * 8;
}


inline LongWord<NativeWidth> & operator>>=(LongWord<NativeWidth> & x, unsigned const shift) {
    SEQAN_CHECKPOINT;
    x._value >>= shift;
    return x;
}


inline LongWord<NativeWidth> & operator<<=(LongWord<NativeWidth> & x, unsigned const shift) {
    SEQAN_CHECKPOINT;
    x._value <<= shift;
    return x;
}


inline LongWord<NativeWidth> & operator&=(LongWord<NativeWidth> & x, LongWord<NativeWidth> const & mask) {
    SEQAN_CHECKPOINT;
    x._value &= mask;
    return x;
}


inline LongWord<NativeWidth> & operator|=(LongWord<NativeWidth> & x, LongWord<NativeWidth> const & mask) {
    SEQAN_CHECKPOINT;
    x._value |= mask;
    return x;
}


inline LongWord<NativeWidth> & operator^=(LongWord<NativeWidth> & x, LongWord<NativeWidth> const & mask) {
    SEQAN_CHECKPOINT;
    x._value ^= mask;
    return x;
}


template <>
struct LongWord<NativeWideWidth> {
};


template <typename TWord, unsigned LENGTH>
struct LongWordBitProxy<TWord, StaticWidth<LENGTH> > {
    TWord &_word;
    unsigned _bitIndex;

    LongWordBitProxy(TWord & word, unsigned bitIndex)
            : _word(word), _bitIndex(bitIndex) {
        SEQAN_CHECKPOINT;
    }

    LongWordBitProxy & operator=(unsigned x) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(x, LENGTH);
        unsigned & block = _word._data[_bitIndex / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = _bitIndex % BitsPerValue<unsigned>::VALUE;
        block = (block ^ (1 << bitIndex)) | (x << bitIndex);
        return *this;
    }

    operator unsigned() const {
        SEQAN_CHECKPOINT;
        unsigned const & block = _word._data[_bitIndex / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = _bitIndex % BitsPerValue<unsigned>::VALUE;
//         std::cerr << "LongWordBitProxy<... StaticWidth<" << LENGTH << "> >" << std::endl;
//         std::cerr << "  _word      = " << _word << std::endl;
//         std::cerr << "  _bitIndex  = " << _bitIndex << std::endl;
//         std::cerr << "  block      = " << block << std::endl;
//         std::cerr << "  bitIndex   = " << bitIndex << std::endl;
//         std::cerr << "  (block & (1 << bitIndex)) >> bitIndex = " << ((block & (1 << bitIndex)) >> bitIndex) << std::endl;
//         std::cerr << std::endl;
        return (block & (1 << bitIndex)) >> bitIndex;
    }
};


template <unsigned LENGTH>
struct LongWord<StaticWidth<LENGTH> > {
    typedef LongWord<StaticWidth<LENGTH> > TWord;
    enum TDummy { UNSIGNED_COUNT = (LENGTH / BitsPerValue<unsigned>::VALUE + ((LENGTH % BitsPerValue<unsigned>::VALUE > 0) ? 1 : 0)) };
    unsigned _data[UNSIGNED_COUNT];

    LongWord() {
        SEQAN_CHECKPOINT;
        // Clear _data.
        for (size_t i = 0; i < UNSIGNED_COUNT; ++i)
            _data[i] = 0u;
//         std::cout << "bits per unsigned = " << BitsPerValue<unsigned>::VALUE << std::endl;
//         std::cout << "# of unsigneds = " << ((LENGTH / BitsPerValue<unsigned>::VALUE) + ((LENGTH % BitsPerValue<unsigned>::VALUE > 0) ? 1 : 0)) << std::endl;
//         std::cout << "UNSIGNED_COUNT = " << UNSIGNED_COUNT << std::endl;
//         std::cout << "sizeof(_data) = " << sizeof(_data) << std::endl;
    }

    LongWord(LongWord const & other) {
        SEQAN_CHECKPOINT;
        for (size_t i = 0; i < UNSIGNED_COUNT; ++i)
            _data[i] = other._data[i];
    }

    LongWord & operator=(LongWord const & other) {
        SEQAN_CHECKPOINT;
        for (size_t i = 0; i < UNSIGNED_COUNT; ++i)
            _data[i] = other._data[i];
    }

    unsigned operator[](unsigned index) const {
        SEQAN_CHECKPOINT;
        unsigned const & block = _data[index / BitsPerValue<unsigned>::VALUE];
        unsigned bitIndex = index % BitsPerValue<unsigned>::VALUE;
        return (block & (1 << bitIndex)) >> bitIndex;
    }

    LongWordBitProxy<TWord, StaticWidth<LENGTH> > operator[](unsigned index) {
        SEQAN_CHECKPOINT;
        SEQAN_ASSERT_LT(index, LENGTH);
        return LongWordBitProxy<TWord, StaticWidth<LENGTH> >(*this, index);
    }
};


template <unsigned LENGTH>
size_t length(LongWord<StaticWidth<LENGTH> > const &) {
    SEQAN_CHECKPOINT;
    return LENGTH;
}


template <unsigned LENGTH>
bool operator==(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
        if (a._data[i] != b._data[i]) {
            return false;
        }
    }
    return true;    
}


template <unsigned LENGTH>
bool operator!=(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
//     for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
//         if (a._data[i] == b._data[i]) {
//             return false;
//         }
//     }
//     return true;
    return !(a == b);
}


template <unsigned LENGTH>
bool operator<=(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
        if (a._data[i] > b._data[i])
            return false;
    }
    return true;    
}


template <unsigned LENGTH>
bool operator>=(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
        if (a._data[i] <= b._data[i])
            return false;
    }
    return true;    
}


template <unsigned LENGTH>
bool operator<(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    return !(a >= b);
//     typedef LongWord<StaticWidth<LENGTH> > TLongWord;
//     // TODO(holtgrew): Roll out loop?
//     for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
//         if (a._data[i] >= b._data[i])
//             return false;
//     }
//     return true;    
}


template <unsigned LENGTH>
bool operator>(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    return !(a <= b);
//     typedef LongWord<StaticWidth<LENGTH> > TLongWord;
//     // TODO(holtgrew): Roll out loop?
//     for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i) {
//         if (a._data[i] <= b._data[i])
//             return false;
//     }
//     return true;    
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator|(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        result._data[i] = a._data[i] | b._data[i];
    return result;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator|=(LongWord<StaticWidth<LENGTH> > & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        a._data[i] |= b._data[i];
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator&(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        result._data[i] = a._data[i] & b._data[i];
    return result;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator&=(LongWord<StaticWidth<LENGTH> > & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        a._data[i] &= b._data[i];
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator^(LongWord<StaticWidth<LENGTH> > const & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        result._data[i] = a._data[i] ^ b._data[i];
    return result;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator^=(LongWord<StaticWidth<LENGTH> > & a, LongWord<StaticWidth<LENGTH> > const & b) {
    SEQAN_CHECKPOINT;
    typedef LongWord<StaticWidth<LENGTH> > TLongWord;
    // TODO(holtgrew): Roll out loop?
    TLongWord result;
    for (unsigned i = 0; i < TLongWord::UNSIGNED_COUNT; ++i)
        a._data[i] ^= b._data[i];
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator>>(LongWord<StaticWidth<LENGTH> > const & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator>>=(LongWord<StaticWidth<LENGTH> > & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > operator<<(LongWord<StaticWidth<LENGTH> > const & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return a;
}


template<unsigned LENGTH>
LongWord<StaticWidth<LENGTH> > & operator<<=(LongWord<StaticWidth<LENGTH> > & a, unsigned shift) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_FAIL("Implement me!");
    return a;
}


template <>
struct LongWord<DynamicWidth> {
};

}  // namespace seqan

namespace std {
template <typename TStream, unsigned LENGTH>
TStream & operator<<(TStream & stream, seqan::LongWord<seqan::StaticWidth<LENGTH> > const & a) {
    SEQAN_CHECKPOINT;
    for (size_t i = 0; i < LENGTH; ++i) {
        unsigned x = a[LENGTH - 1 - i];
        stream << x;
    }
    return stream;
}
}  // namespace std

#endif  // SEQAN_MISC_MISC_LONG_WORD_H_
