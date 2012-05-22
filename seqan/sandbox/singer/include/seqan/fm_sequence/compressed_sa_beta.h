// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef COMPRESSED_SA_BETA_H_
#define COMPRESSED_SA_BETA_H_

namespace seqan {

template <typename TSparseString, typename TLfTable, typename TSpec>
class CompressedSA;

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>
{
    typedef Iter<CompressedSA<TSparseString, TLfTable, TSpec> const, PositionIterator> Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>
{
    typedef Iter<CompressedSA<TSparseString, TLfTable, TSpec>, PositionIterator> Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Rooted>:
    Iterator<CompressedSA<TSparseString, TLfTable, TSpec>, Standard>{};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Rooted>:
    Iterator<CompressedSA<TSparseString, TLfTable, TSpec> const, Standard>{};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Value<CompressedSA<TSparseString, TLfTable, TSpec> >
{
    typedef typename Value<TSparseString>::Type Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Value<CompressedSA<TSparseString, TLfTable, TSpec> const>
{
    typedef typename Value<TSparseString>::Type const Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Reference<CompressedSA<TSparseString, TLfTable, TSpec> >
{
    typedef typename Value<CompressedSA<TSparseString, TLfTable, TSpec> >::Type Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Reference<const CompressedSA<TSparseString, TLfTable, TSpec> >
{
    typedef typename Value<CompressedSA<TSparseString, TLfTable, TSpec> >::Type const Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
class CompressedSA
{
public:
    TSparseString   sparseSA;
    TLfTable * lfTable;

    CompressedSA() :
        sparseSA(),
        lfTable()
    {}

    CompressedSA(unsigned blockSize) :
        sparseSA(blockSize),
        lfTable()
    {}

    CompressedSA(unsigned blockSize, TLfTable & lfTable) :
        sparseSA(blockSize),
        lfTable(&lfTable)
    {}

    inline CompressedSA & operator=(CompressedSA const & other)
    {
        sparseSA = other.sparseSA;
        lfTable = other.lfTable;
        return *this;
    }

    typedef typename Value<typename Fibre<TSparseString, FibreValueString>::Type>::Type TCompressedSaValue;
    typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;

    template <typename TPos>
    inline TCompressedSaValue const operator[](TPos pos)
    {
        TIndicatorString const & indicatorString = getFibre(sparseSA, FibreIndicatorString());
        TPos counter = 0;

        while (!getBit(indicatorString, pos))
        {
            pos = lfMapping(*lfTable, pos);
            ++counter;
        }
        return getValue(sparseSA, getRank(indicatorString, pos) - 1) + counter;
    }

    template <typename TPos>
    inline TCompressedSaValue operator[](TPos pos) const
    {
        TIndicatorString const & indicatorString = getFibre(sparseSA, FibreIndicatorString());
        TPos counter = 0;
        while (!getBit(indicatorString, pos))
        {
            pos = lfMapping(*lfTable, pos);
            ++counter;
        }
        return getValue(sparseSA, getRank(indicatorString, pos) - 1) + counter;
    }

    inline bool operator==(const CompressedSA & other) const
    {
        return sparseSA == other.sparseSA &&
               *lfTable == *(other.lfTable);
    }

};

//template <typename TSpecPairI1, typename TSpecPairI2, typename TSpecPairSpec, typename TStringSpec, typename TSparseStringSpec, typename TLfTable, typename TSpec>
//struct CompressedSA<SparseString<String<Pair<TSpecPairI1, TSpecPairI2, TSpecPairSpec>, TStringSpec>, TSparseStringSpec>, TLfTable, TSpec>
//{
//    typedef SparseString<String<Pair<TSpecPairI1, TSpecPairI2, TSpecPairSpec>, TStringSpec>, TSparseStringSpec>   TSparseString;
//    typedef typename Value<typename Fibre<SparseString<TSparseString, TStringSpec>, FibreSparseString>::Type>::Type TCompressedSaValue;
//    typedef typename Fibre<SparseString<TSparseString, TStringSpec>, FibreFibreIndicatorString>::Type TFibreIndicatorString;
//
//    TSparseString   sparseSA;
//    TLfTable *        lfTable;
//
//    CompressedSA() :
//      sparseSA(),
//      lfTable()
//    {}
//
//    template <typename TPos>
//    inline TCompressedSaValue const operator[](TPos pos)
//    {
//        TFibreIndicatorString const & FibreIndicatorString = getFibre(sparseSA, FibreFibreIndicatorString());
//        TPos counter = 0;
//        while (!getBit(FibreIndicatorString, pos))
//        {
//            pos = lfMapping(*lfTable, pos);
//            ++counter;
//        }
//        TCompressedSaValue temp = getValue(sparseSA, getRank(FibreIndicatorString, pos) - 1);
//        temp.i2 += counter;
//        return temp;
//    }
//
//    template <typename TPos>
//    inline TCompressedSaValue operator[](TPos pos) const
//    {
//        TFibreIndicatorString const & FibreIndicatorString = getFibre(sparseSA, FibreFibreIndicatorString());
//        TPos counter = 0;
//        while (!getBit(FibreIndicatorString, pos))
//        {
//            pos = lfMapping(*lfTable, pos);
//            ++counter;
//        }
//        TCompressedSaValue temp = getValue(sparseSA, getRank(FibreIndicatorString, pos) - 1);
//        temp.i2 += counter;
//        return temp;
//    }
//
//    inline bool operator==(const CompressedSA & b) const
//    {
//        return sparseSA == b.sparseSA &&
//               *lfTable == *(b.lfTable);
//    }
//
//};

struct FibreSparseString_;
typedef Tag<FibreSparseString_> const FibreSparseString;

template <typename TSparseString, typename TLfTable, typename TSpec>
struct Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>
{
    typedef TSparseString Type;
};

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type const &
getFibre(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, FibreSparseString)
{
    return compressedSA.sparseSA;
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type &
getFibre(CompressedSA<TSparseString, TLfTable, TSpec>&compressedSA, FibreSparseString)
{
    return compressedSA.sparseSA;
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline bool getNextPos(CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA, TPos & pos)
{
    typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;
    TIndicatorString const & indicatorString = compressedSA.sparseSA.indicatorString;

    if (getBit(indicatorString, pos))
    {
        return true;
    }
    pos = lfMapping(*compressedSA.lfTable, pos);
    return false;
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TValue>
void assignBlockSize(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TValue value)
{
    assignBlockSize(getFibre(compressedSA, FibreSparseString()), value);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
void assignLfTable(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TLfTable & lfTable)
{
    compressedSA.lfTable = &lfTable;
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos, typename TValue>
void assignValue(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos, TValue value)
{
    assignValue(getFibre(compressedSA, FibreSparseString()), pos, value);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Size<typename Fibre<CompressedSA<TSparseString, TLfTable, TSpec>, FibreSparseString>::Type>::Type
length(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return length(getFibre(compressedSA, FibreSparseString()));
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline void clear(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    clear(getFibre(compressedSA, FibreSparseString()));
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline typename Size<typename Fibre<TSparseString, FibreValueString>::Type>::Type
getBlockSize(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return getBlockSize(getFibre(compressedSA, FibreSparseString()));
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
inline void resize(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TSize size)
{
    resize(getFibre(compressedSA, FibreSparseString()), size);
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
inline void resize(CompressedSA<String<TSparseString>, TLfTable, TSpec> & compressedSA, TSize size)
{
    resize(getFibre(compressedSA, FibreSparseString()), size);
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSize>
inline void reserve(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
                    TSize size)
{
    resize(getFibre(compressedSA, FibreSparseString()), size);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >
begin(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >(compressedSA, 0);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >
end(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA)
{
    return Iterator<CompressedSA<TSparseString, TLfTable, TSpec> >(compressedSA, length(compressedSA.compressedSA));
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TSA>
void fill(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TSA const & completeSA)
{
    typedef CompressedSA<TSparseString, TLfTable, TSpec> TCompressedSA;
    typedef typename GetValue<TSA>::Type                            TSAValue;
    typedef typename Size<TSA>::Type                                TSize;
    typedef typename Fibre<TCompressedSA, FibreSA>::Type            TSparseSA;
    typedef typename Fibre<TSparseSA, FibreIndicatorString>::Type   TIndicatorString;

    TSparseSA & sparseSA = getFibre(compressedSA, FibreSA());
    TIndicatorString & indicatorString = getFibre(sparseSA, FibreIndicatorString());

    TSize n = length(completeSA);
    resize(compressedSA, n);

    for (TSize i = 0; i < n; i++)
    {
        TSAValue sa = getValue(completeSA, i);
        if (sa % getBlockSize(compressedSA) == 0)
        {
            setBit(indicatorString, i + 1, 1);
        }
    }

    completeRankSupportBitString(indicatorString);

    TSize counter = 0;
    for (TSize i = 0; i < n; i++)
    {
        if (getBit(indicatorString, i + 1))
        {
            assignValue(compressedSA, counter, completeSA[i]);
            ++counter;
        }
    }
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline typename Value<typename Fibre<TSparseString, FibreSparseString>::Type>::Type
value(CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos)
{
    return compressedSA[pos];
}

template <typename TSparseString, typename TLfTable, typename TSpec, typename TPos>
inline typename Value<typename Fibre<TSparseString, FibreSparseString>::Type>::Type const
value(const CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA, TPos pos)
{
    return compressedSA[pos];
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool open(
    CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".sstring");
    if (!open(getFibre(compressedSA, FibreSparseString()), toCString(name), openMode))
    {
        return false;
    }
    return true;

}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool open(
    CompressedSA<TSparseString, TLfTable, TSpec> & compressedSA,
    const char * fileName)
{
    return open(compressedSA, fileName, DefaultOpenMode<CompressedSA<TSparseString, TLfTable, TSpec> >::VALUE);
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool save(
    CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".sstring");
    if (!save(getFibre(compressedSA, FibreSparseString()), toCString(name), openMode))
    {
        return false;
    }
    return true;
}

template <typename TSparseString, typename TLfTable, typename TSpec>
inline bool save(
    CompressedSA<TSparseString, TLfTable, TSpec> const & compressedSA,
    const char * fileName)
{
    return save(compressedSA, fileName, DefaultOpenMode<CompressedSA<TSparseString, TLfTable, TSpec> >::VALUE);
}

}


#endif // COMPRESSED_SA_BETA_H_