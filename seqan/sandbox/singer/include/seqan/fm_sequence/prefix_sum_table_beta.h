// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <your.email@example.net>
// ==========================================================================

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_PREFIX_SUM_TABLE_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_PREFIX_SUM_TABLE_H_

#include <seqan/fm_index/fm_index_beta.h>

namespace seqan {

template <typename TChar, typename TSpec>
class PrefixSumTable;

struct FibreEntries_;
typedef Tag<FibreEntries_> const FibreEntries;

template <typename TChar, typename TSpec>
struct GetValue<PrefixSumTable<TChar, TSpec> >
{
    typedef typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type TEntries;
    typedef typename Value<TEntries>::Type Type;
};

template <typename TChar, typename TSpec>
struct GetValue<PrefixSumTable<TChar, TSpec> const>
{
    typedef typename Fibre<PrefixSumTable<TChar, TSpec> const, FibreEntries>::Type TEntries;
    typedef typename Value<TEntries>::Type Type;
};

template <typename TChar, typename TSpec>
struct Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>
{
    typedef String<unsigned> Type;
};

template <typename TChar, typename TSpec>
struct Fibre<PrefixSumTable<TChar, TSpec> const, FibreEntries>
{
    typedef String<unsigned> const Type;
};

template <typename TChar, typename TSpec>
struct Size<PrefixSumTable<TChar, TSpec> >
{
    typedef typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type TEntries;
    typedef typename Size<TEntries>::Type Type;
};

template <typename TChar, typename TSpec>
struct Infix<PrefixSumTable<TChar, TSpec> >
{
    typedef typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type TEntries;
    typedef typename Infix<TEntries>::Type Type;
};

template <typename TChar, typename TSpec>
struct Infix<PrefixSumTable<TChar, TSpec> const>
{
    typedef typename Fibre<PrefixSumTable<TChar, TSpec> const, FibreEntries>::Type TEntries;
    typedef typename Infix<TEntries>::Type Type;
};

template <typename TSpec>
struct CharacterValue;

template <typename TChar, typename TSpec>
struct CharacterValue<PrefixSumTable<TChar, TSpec> >
{
    typedef TChar Type;
};

template <typename TChar, typename TSpec>
struct CharacterValue<PrefixSumTable<TChar, TSpec> const>
{
    typedef TChar const Type;
};

template <typename TChar, typename TSpec>
struct Value<PrefixSumTable<TChar, TSpec> >
{
    typedef typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type TEntries;
    typedef typename Value<TEntries>::Type Type;
};

template <typename TChar, typename TSpec>
struct Value<PrefixSumTable<TChar, TSpec> const>
{
    typedef typename Fibre<PrefixSumTable<TChar, TSpec> const, FibreEntries>::Type TEntries;
    typedef typename Value<TEntries>::Type const Type;
};

template <typename TChar, typename TSpec>
class PrefixSumTable
{
    typedef typename Fibre<PrefixSumTable, FibreEntries>::Type TEntries;
    typedef typename Value<TEntries>::Type TEntry;

public:
    TEntries entries;

    PrefixSumTable() :
        entries()
    {}

    PrefixSumTable(String<TChar> const & text) :
        entries()
    {
        createPrefixSumTable(*this, text);
    }

    template <typename TPos>
    inline TEntry & operator[](TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline TEntry const & operator[](TPos pos) const
    {
        return value(*this, pos);
    }

};

template <typename TChar, typename TSpec, typename TText>
inline void createPrefixSumTable(PrefixSumTable<TChar, TSpec> & prefixSumTable, TText & text)
{
    typedef PrefixSumTable<TChar, TSpec> TPrefixSumTable;
    typedef typename Value<typename Fibre<TPrefixSumTable, FibreEntries>::Type>::Type TPrefixSumValue;

    TPrefixSumTable freq;
    getFrequencies(freq, text);

    unsigned alpSize = length(freq);
    resize(prefixSumTable, alpSize + 1, 0);


    TPrefixSumValue temp = 0;
    TPrefixSumValue sum = 0; // = getNumSequences(text);  // $ is smaller than all other chararcters
    for (TPrefixSumValue i = 0; i < alpSize; ++i)
    {
        //std::cerr << i << " " << alpSize << std::endl;
        temp = getPrefixSum(freq, i);
        setPrefixSum(prefixSumTable, sum, i);
        sum += temp;
    }

    //std::cerr << length(freq) << " " << sum << " " << alpSize << std::endl;

    setPrefixSum(prefixSumTable, sum, alpSize);
}

template <typename TChar, typename TSpec>
unsigned getAlphabetSize(PrefixSumTable<TChar, TSpec> const & pst)
{
    return length(pst.entries) - 1;
}

template <typename TDummy, typename TChar>
inline unsigned getCharacterPosition(TDummy const & /*tag*/, TChar character)
{
    return ordValue(character);
}

template <typename TChar, typename TSpec, typename TChar2>
inline unsigned getCharacterPosition(PrefixSumTable<TChar, TSpec> const & /*tag*/, TChar2 character)
{
    //static_cast<Nothing>(character);
    return ordValue(character);
}

template <typename TChar, typename TSpec>
inline int getCharacterPosition(PrefixSumTable<TChar, TSpec> const & /*tag*/, char character)
{
    return static_cast<int>(character) + 128;
}

template <typename TSpec>
inline int getCharacterPosition(PrefixSumTable<char, TSpec> const & /*tag*/, char character)
{
    return static_cast<int>(character) + 128;
}

template <typename TSpec, typename TPos>
inline char getCharacter(PrefixSumTable<char, TSpec> const & /*tag*/, TPos const pos)
{
    return static_cast<char>(pos - 128);
}

template <typename TChar, typename TSpec, typename TPos>
inline TChar getCharacter(PrefixSumTable<TChar, TSpec> const & /*tag*/, TPos const pos)
{
    return static_cast<TChar>(pos);
}

template <typename TText>
unsigned getNumSequences(TText const & /*tag*/)
{
    return 1;
}

template <typename TSequence>
unsigned getNumSequences(StringSet<TSequence> const & stringSet)
{
    return length(stringSet);
}

template <typename TChar, typename TSpec, typename TBeginPos, typename TEndPos>
unsigned getPivotPosition(PrefixSumTable<TChar, TSpec> const & pst, TBeginPos beginPos, TEndPos endPos)
{
    TBeginPos realBeginPos = beginPos + 1;
    TEndPos realEndPos = endPos + 1;
    unsigned lengthRange = realEndPos - realBeginPos + 1;
    unsigned pivotPos = realBeginPos + lengthRange / 2;

    unsigned tooSmallValues = pst[beginPos];
    unsigned pivotValue = (pst[realEndPos] - tooSmallValues) / 2;

    int direction;
    (pst[pivotPos] - tooSmallValues > pivotValue) ? direction = -1 : direction = 1;
    if (direction == -1)
    {
        //while(pivotPos > 0 && pst[pivotPos] - tooSmallValues > pivotValue)
        while (pivotPos > realBeginPos && pst[pivotPos] - tooSmallValues > pivotValue)
        {
            --pivotPos;
        }
        if (pivotPos == beginPos)
            ++pivotPos;
    }
    else
    {
        while (pst[pivotPos] - tooSmallValues <= pivotValue)
        {
            ++pivotPos;
        }
        --pivotPos;
    }
    return pivotPos;
}

template <typename TChar, typename TSpec, typename TPos>
typename Value<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type &
prefixSum(PrefixSumTable<TChar, TSpec>&pst, TPos const pos)
{
    return value(pst, pos);
}

template <typename TChar, typename TSpec, typename TPos>
typename Value<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type const &
prefixSum(PrefixSumTable<TChar, TSpec> const & pst, TPos const pos)
{
    return value(pst, pos);
}

template <typename TChar, typename TSpec, typename TPos>
typename Value<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type
getPrefixSum(PrefixSumTable<TChar, TSpec> const & pst, TPos const pos)
{
    return getValue(pst, pos);
}

template <typename TChar, typename TSpec, typename TPosBegin, typename TPosEnd>
typename Infix<PrefixSumTable<TChar, TSpec> >::Type
getPrefixSumTableRange(PrefixSumTable<TChar, TSpec> const & pst, TPosBegin posBegin, TPosEnd posEnd)
{
    return infix(pst, posBegin + 1, posEnd + 1);
}

template <typename TChar, typename TSpec, typename TPosBegin, typename TPosEnd>
typename Infix<PrefixSumTable<TChar, TSpec> >::Type
getPrefixSumTableRange(PrefixSumTable<TChar, TSpec> & pst, TPosBegin posBegin, TPosEnd posEnd)
{
    return infix(pst, posBegin + 1, posEnd + 1);
}

template <typename TChar, typename TSpec, typename TPos>
inline typename GetValue<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type
getValue(PrefixSumTable<TChar, TSpec> & pst, TPos const pos)
{
    return pst.entries[pos];
}

template <typename TChar, typename TSpec, typename TPos>
inline typename GetValue<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type const
getValue(PrefixSumTable<TChar, TSpec> const & pst, TPos const pos)
{
    return pst.entries[pos];
}

template <typename TChar, typename TSpec, typename TPosBegin, typename TPosEnd>
typename Infix<PrefixSumTable<TChar, TSpec> >::Type
infix(PrefixSumTable<TChar, TSpec> & pst, TPosBegin start, TPosEnd end)
{
    return infix(pst.entries, start, end);
}

template <typename TChar, typename TSpec, typename TPosBegin, typename TPosEnd>
typename Infix<PrefixSumTable<TChar, TSpec> >::Type
infix(PrefixSumTable<TChar, TSpec> const & pst, TPosBegin start, TPosEnd end)
{
    return infix(pst.entries, start, end);
}

template <typename TChar, typename TSpec, typename TNumDollar>
void insertDollar(PrefixSumTable<TChar, TSpec> & pst, TNumDollar const numDollar)
{
    for (unsigned i = 0; i < length(pst); ++i)
        prefixSum(pst, i) = getPrefixSum(pst, i) + numDollar;
}

template <typename TChar, typename TSpec>
inline typename Size<PrefixSumTable<TChar, TSpec> >::Type
length(PrefixSumTable<TChar, TSpec> const & pst)
{
    return length(pst.entries);
}

template <typename TChar, typename TSpec, typename TSize>
inline void resize(PrefixSumTable<TChar, TSpec> const & pst, TSize size)
{
    resize(pst.entries, size);
}

template <typename TChar, typename TSpec, typename TSize, typename TValue>
inline void resize(PrefixSumTable<TChar, TSpec> & pst, TSize size, TValue value)
{
    resize(pst.entries, size, value);
}

template <typename TChar, typename TSpec, typename TPos>
inline void setCharacter(PrefixSumTable<TChar, TSpec> & /*pst*/, TChar /*tag*/, TPos /*tag*/)
{}

template <typename TChar, typename TSpec, typename TValue, typename TPos>
inline void setPrefixSum(PrefixSumTable<TChar, TSpec> & pst, TValue value, TPos const pos)
{
    pst.entries[pos] = value;
}

template <typename TChar, typename TSpec, typename TPos>
inline typename Value<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type &
value(PrefixSumTable<TChar, TSpec>&pst, TPos const pos)
{
    //return value(pst.entries, static_cast<unsigned>(pos));
    return pst.entries[pos];
}

template <typename TChar, typename TSpec, typename TPos>
inline typename Value<typename Fibre<PrefixSumTable<TChar, TSpec>, FibreEntries>::Type>::Type const &
value(PrefixSumTable<TChar, TSpec> const & pst, TPos const pos)
{
    return pst.entries[pos];
}

/*
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<PrefixSumTable<char, TSpec>, FibreEntries>::Type>::Type &
value(PrefixSumTable<char, TSpec> & pst, TPos const pos)
{
    return pst.entries[pos];
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<PrefixSumTable<char, TSpec>, FibreEntries>::Type>::Type const &
value(PrefixSumTable<char, TSpec> const & pst, TPos const pos)
{
    return pst.entries[pos];
}
*/

}


#endif // SANDBOX_MY_SANDBOX_APPS_FMINDEX_PREFIX_SUM_TABLE_H_