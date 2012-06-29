// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_BITSTRING_BETA_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_BITSTRING_BETA_H_

namespace seqan {

//Rank Support Bit String
template <typename TSpec = void>
struct RankSupportBitString;

// FM index fibres

/**
.Tag.Rank Support Bit String Index Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.RankSupportBitString@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a rank support bit string. 
..cat:Index

..tag.FibreBitString:The bit string. 
..tag.FibreBucketString:The bucket string.
..tag.FibreSuperBucketString:The super bucket string.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

struct FibreBitString_;
struct FibreBucketString_;
struct FibreSuperBucketString_;

typedef Tag<FibreBitString_> const          FibreBitString;         typedef FibreBitString          RankSupportBitStringBitString;
typedef Tag<FibreBucketString_> const       FibreBucketString;      typedef FibreBucketString       RankSupportBitStringBucketString;
typedef Tag<FibreSuperBucketString_> const  FibreSuperBucketString; typedef FibreSuperBucketString  RankSupportBitStringSuperBucketString;

// ==========================================================================
// Metafunctions
// ==========================================================================

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreBitString>
{
    typedef String<unsigned long> Type;
};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreBucketString>
{
    typedef String<unsigned short> Type;
};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>
{
    typedef String<unsigned int> Type;
};

//The limiting factor of the size is the underlying data type of the super bucket string
template <typename TSpec>
struct Size<RankSupportBitString<TSpec> >
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type>::Type Type;
};

// ==========================================================================
// Classes
// ==========================================================================

/**
.Class.RankSupportBitString:
..summary:A bit string supporting rank queries in constant time.
..cat:Index
..signature:RankSupportBitString<TSpec>
..param.TSpec:Specialisation tag. 
...default:void
..remarks:The constant rank query time is achieved by evaluating precomputed subsolutions. In order to do so, the bit string is divided into buckets of length l. A super bucket string stores for each block of l buckets the number of bits set from the beginning. In addition a bucket string stores the number of bits set in each bucket from the start of the last super bucket block. Therefore it is possible to compute the result of a rank query in constant time by adding information from the bit, bucket and super bucket string.
..include:seqan/index.h
*/
template <typename TSpec>
struct RankSupportBitString
{
    typedef typename Fibre<RankSupportBitString, FibreBitString>::Type      TBitString;
    typedef typename Fibre<RankSupportBitString, FibreBucketString>::Type     TBucketString;
    typedef typename Fibre<RankSupportBitString, FibreSuperBucketString>::Type    TSuperBucketString;

    TBitString                  bString;
    TBucketString               buString;
    TSuperBucketString          sBuString;
    typename Size<RankSupportBitString>::Type   length_;

    RankSupportBitString() :
        bString(),
        buString(),
        sBuString(),
        length_(0)
    {}

    template <typename TString>
    RankSupportBitString(TString const & input) :
        bString(),
        buString(),
        sBuString(),
        length_(length(input))
    {
        typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type     TFibreSuperBucketString;
        typedef typename Value<TFibreSuperBucketString>::Type                 TFibreSuperBucketStringValue;

        resize(*this, length_);
        for (TFibreSuperBucketStringValue i = 0; i < length(*this); ++i)
        {
            setBit(*this, i, input[i]);
        }
        completeRankSupportBitString(*this);
    }

    inline RankSupportBitString & operator=(RankSupportBitString const & other)
    {
        bString = other.bString;
        buString = other.buString;
        sBuString = other.sBuString;
        length_ = other.length_;
        return *this;
    }

    inline bool operator==(const RankSupportBitString & other) const
    {
        return length_ == other.length_ &&
               bString == other.bString &&
               buString == other.buString &&
               sBuString == other.sBuString;
    }

};

/**
.Function.appendValue:
..param.target:
...type:Class.RankSupportBitString
*/
template <typename TSpec, typename TBit>
inline void appendValue(RankSupportBitString<TSpec> & bitString, TBit const bit)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBitString>::Type   TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBucketString>::Type        TFibreBucketString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type       TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type               TFibreBitStringValue;
    typedef typename Value<TFibreBucketString>::Type              TFibreBucketStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type             TFibreSuperBucketStringValue;

    TFibreBucketStringValue const bitsPerValue_ = BitsPerValue<TFibreBitStringValue>::VALUE;
    TFibreSuperBucketStringValue const length_ = length(bitString);

    //check if the new size exceeds the reserved memory
    if ((length(bitString.bString) * bitsPerValue_) <= length_)
    {
        reserve(bitString, 2 * (length_ + 1));
    }

    //initialize current (super)bucket with preceding one
    takeBuValueAlong_(bitString, length_, bitsPerValue_);
    takeSBuValueAlong_(bitString, length_, bitsPerValue_);

    if (bit)
    {
        setBit(bitString, length_, bit);

        //determine whether
        TFibreSuperBucketStringValue buPos = getBuPos_(bitString, length_);
        if (((buPos + 1) % (bitsPerValue_)))
        {
            ++getFibre(bitString, FibreBucketString())[buPos + 1];
        }
        TFibreSuperBucketStringValue sBuPos = getSBuPos_(bitString, length_);
        ++getFibre(bitString, FibreSuperBucketString())[sBuPos + 1];
    }
    ++bitString.length_;
}

/**
.Function.clear
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec>
inline void clear(RankSupportBitString<TSpec> & bitString)
{
    clear(bitString.bString);
    clear(bitString.buString);
    clear(bitString.sBuString);
}

// TODO (singer): Why do I need this forward?
template <typename TValue>
inline unsigned getRankInBucket_(TValue const value);

/**
.Function.completeRankSupportBitString
..summary:Adds the bucket and super bucket information to the bit string.
..signature:completeRankSupportBitString(bitString)
..param.bitString:The bit string to be completed.
...type:Class.RankSupportBitString
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));

for (unsigned i = 0; i < length(genome); ++i)
    if(genome[i] < Dna5('c'))
        setBit(bitString, 1);

completeRankSupportBitString(bitString);
*/
template <typename TSpec>
inline void completeRankSupportBitString(RankSupportBitString<TSpec> & bitString)
{
    if (length(bitString))
    {

        typedef RankSupportBitString<TSpec>                                         TRankSupportBitString;
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreBitString>::Type   TFibreBitString;
        typedef typename Fibre<TRankSupportBitString, FibreBucketString>::Type      TFibreBucketString;
        typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type TFibreSuperBucketString;
        typedef typename Value<TFibreBitString>::Type                               TFibreBitStringValue;
        typedef typename Value<TFibreBucketString>::Type                            TFibreBucketStringValue;
        typedef typename Value<TFibreSuperBucketString>::Type                       TFibreSuperBucketStringValue;

        TFibreBucketStringValue buSum_ = 0;
        TFibreSuperBucketStringValue sBuSum_ = 0;
        TFibreSuperBucketStringValue superBucketCounter = 1;
        TFibreBucketStringValue const bitsPerValue = BitsPerValue<TFibreBitStringValue>::VALUE;

        for (TFibreSuperBucketStringValue i = 0; i < length(bitString.bString) - 1; ++i)
        {
            buSum_ += getRankInBucket_(bitString.bString[i]);
            bitString.buString[i + 1] = buSum_;
            if (!((i + 1) % bitsPerValue))
            {
                sBuSum_ += buSum_;
                bitString.sBuString[superBucketCounter] = sBuSum_;
                buSum_ = 0;
                bitString.buString[i + 1] = buSum_;
                ++superBucketCounter;
            }
        }
    }
}

/**
.Function.empty
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec>
inline bool empty(RankSupportBitString<TSpec> & bitString)
{
    return empty(getFibre(bitString, FibreBitString()))
        && empty(getFibre(bitString, FibreBucketString()))
        && empty(getFibre(bitString, FibreSuperBucketString()));
}

/**
.Function.getBit
..summary:Returns whether a specified bit is set or not.
..signature:getBit(bitString, pos)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of the bit.
..returns:Returns whether a specified bit is set or not.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));
...
mark all 'a's
...

for (unsigned i = 0; i < length(bitString); ++i)
    if(getBit(bitString, i))
        std::cout << "a found at: " << i << std::endl;
*/
template <typename TSpec, typename TPos>
inline bool getBit(RankSupportBitString<TSpec> & bitString, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBitString>::Type         TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBucketString>::Type        TFibreBucketString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type       TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type               TFibreBitStringValue;
    typedef typename Value<TFibreBucketString>::Type              TFibreBucketStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type             TFibreSuperBucketStringValue;

    TFibreBucketStringValue const bitsPerValue = BitsPerValue<TFibreBitStringValue>::VALUE;
    TFibreBitStringValue const one = 1;
    TFibreBitStringValue const shiftValue = pos % bitsPerValue;
    TFibreBitStringValue const buPos = getBuPos_(bitString, pos);
    return (bitString.bString[buPos] >> shiftValue) & one;
}

template <typename TSpec, typename TPos>
inline bool getBit(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBitString>::Type   TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBucketString>::Type        TFibreBucketString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type       TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type               TFibreBitStringValue;
    typedef typename Value<TFibreBucketString>::Type              TFibreBucketStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type             TFibreSuperBucketStringValue;

    TFibreBucketStringValue const bitsPerValue = BitsPerValue<TFibreBitStringValue>::VALUE;
    TFibreBitStringValue const one = 1;
    TFibreBitStringValue const shiftValue = pos % bitsPerValue;
    TFibreBitStringValue const buPos = getBuPos_(bitString, pos);
    return (bitString.bString[buPos] >> shiftValue) & one;
}

// This function returns the position in the bucket string of the corresponding bucket.
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type>::Type
getBuPos_(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                             TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBitString>::Type     TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type   TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type           TFibreBitStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type         TFibreSuperBucketStringValue;

    TFibreSuperBucketStringValue const bitsPerValue_ = BitsPerValue<TFibreBitStringValue>::VALUE;
    return pos / bitsPerValue_;
}

/**
.Function.getFibre
..param.container:
...type:Class.RankSupportBitString
..param.fibreTag:
...type:Tag.Rank Support Bit String Fibres
*/
template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBitString>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreBitString)
{
    return string.bString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBitString>::Type const &
getFibre(RankSupportBitString<TSpec> const & string, const FibreBitString)
{
    return string.bString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBucketString>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreBucketString)
{
    return string.buString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBucketString>::Type const &
getFibre(RankSupportBitString<TSpec> const & string, const FibreBucketString)
{
    return string.buString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreSuperBucketString)
{
    return string.sBuString;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type const &
getFibre(const RankSupportBitString<TSpec>&string, const FibreSuperBucketString)
{
    return string.sBuString;
}

// This function returns the position of a specified bit within a bucket.
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreBucketString>::Type>::Type
getPosInBu_(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBitString>::Type     TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type   TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type           TFibreBitStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type         TFibreSuperBucketStringValue;

    TFibreSuperBucketStringValue bitsPerValue = BitsPerValue<TFibreBitStringValue>::VALUE;
    return pos % bitsPerValue;
}


/**
.Function.getRank
..summary:Returns the rank (the number of bits set from the start of the bit string) of a specified position.
..signature:getRank(bitString, pos)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of a bit.
..returns:Value type of the super bucket fibre (default unsigned long).
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));
    
    ...
    mark all 'a's
    ...

for (unsigned i = 0; i < length(bitString); ++i)
    if(getBit(bitString, i))
        std::cout << "found the << getRank(bitString, i) << " a at: " << i << std::endl;
*/
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type>::Type
getRank(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBucketString>::Type        TFibreBucketString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type       TFibreSuperBucketString;
    typedef typename Value<TFibreBucketString>::Type              TFibreBucketStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type             TFibreSuperBucketStringValue;

    TFibreSuperBucketStringValue const buPos = getBuPos_(bitString, pos);
    TFibreSuperBucketStringValue const sBuPos = getSBuPos_(bitString, pos);

    return getRankInBucket_(bitString, pos)
           + bitString.buString[buPos]
           + bitString.sBuString[sBuPos];
}

// This function returns the number of bits set in a bucket until a specified position
template <typename TValue>
inline unsigned getRankInBucket_(TValue const value)
{
    return getRankInBucket_(value, typename Eval<(BitsPerValue<TValue>::VALUE > 32)>::Type());
}

template <typename TValue>
inline unsigned getRankInBucket_(TValue const value, False)
{
    return __builtin_popcountl(static_cast<int32_t>(value));
}

template <typename TValue>
inline unsigned getRankInBucket_(TValue const value, True)
{
    return __builtin_popcountll(static_cast<int64_t>(value));
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type>::Type
getRankInBucket_(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBitString>::Type   TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBucketString>::Type        TFibreBucketString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type       TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type               TFibreBitStringValue;
    typedef typename Value<TFibreBucketString>::Type              TFibreBucketStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type             TFibreSuperBucketStringValue;

    TFibreBucketStringValue const posInBu = getPosInBu_(bitString, pos);
    TFibreSuperBucketStringValue const buPos = getBuPos_(bitString, pos);
    TFibreBucketStringValue const bitsPerValue = BitsPerValue<TFibreBitStringValue>::VALUE;
    TFibreBitStringValue const one = -1;
    TFibreBitStringValue const mask = one >> (bitsPerValue - posInBu - 1);
    //std::cerr << "one: " << one << " mask1: " << mask << " " << (getRankInBucket_(bitString.bString[buPos] & mask)) <<  std::endl;
    return getRankInBucket_(bitString.bString[buPos] & mask);
}

// This function returns the number of bits set in a bucket until a specified position
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type>::Type
getSBuPos_(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                             TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBitString>::Type     TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type   TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type           TFibreBitStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type         TFibreSuperBucketStringValue;

    TFibreSuperBucketStringValue const bitsPerValue_ = BitsPerValue<TFibreBitStringValue>::VALUE;

    return pos / (bitsPerValue_ * bitsPerValue_);
}

/**
.Function.length
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type>::Type
length(RankSupportBitString<TSpec> const & bitString)
{
    return bitString.length_;
}

template <typename TSpec>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBucketString>::Type>::Type
length(RankSupportBitString<TSpec> & bitString)
{
    return bitString.length_;
}


/**
.Function.reserve
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec, typename TSize, typename TValue>
inline void reserve(RankSupportBitString<TSpec> & bitString, TSize const size, TValue const /*dummy*/)
{
    reserve(bitString, size);
}

template <typename TSpec, typename TSize>
inline void reserve(RankSupportBitString<TSpec> & bitString, TSize const size)
{
    typedef RankSupportBitString<TSpec>                             TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBitString>::Type     TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBucketString>::Type    TFibreBucketString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type   TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type           TFibreBitStringValue;
    typedef typename Value<TFibreBucketString>::Type          TFibreBucketStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type         TFibreSuperBucketStringValue;

    if (size)
    {
        TFibreBucketStringValue bitsPerBucket_ = BitsPerValue<TFibreBitStringValue>::VALUE;
        TFibreSuperBucketStringValue numberOfBuckets_ = size / bitsPerBucket_;

        resize(bitString.bString, numberOfBuckets_ + 1, 0);
        resize(bitString.buString, numberOfBuckets_ + 2, 0);
        resize(bitString.sBuString, numberOfBuckets_ / bitsPerBucket_ + 2, 0);
    }
    else
    {
        resize(bitString.bString, 0, 0);
        resize(bitString.buString, 0, 0);
        resize(bitString.sBuString, 0, 0);
    }
}

/**
.Function.resize
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec, typename TSize>
inline void resize(RankSupportBitString<TSpec> & bitString, TSize const size)
{
    reserve(bitString, size);
    bitString.length_ = size;
}

template <typename TSpec, typename TSize, typename TValue>
inline void resize(RankSupportBitString<TSpec> & bitString, TSize const size, TValue const value)
{
    reserve(bitString, size, value);
    bitString.length_ = size;
}

/**
.Function.setBit
..summary:Set a specified bit to true or false.
..signature:setBit(bitString, pos, bit)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of the bit.
..param.bit:The value of the bit.
...remarks:Note that values different from 0 are interpreted as 1.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));

for (unsigned i = 0; i < length(genome); ++i)
    if(genome[i] < Dna5('c'))
        setBit(bitString, 1);

completeRankSupportBitString(bitString);
*/
template <typename TSpec, typename TPos, typename TBit>
inline void setBit(RankSupportBitString<TSpec> & bitString, TPos const pos, TBit const setBit)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBitString>::Type   TFibreBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBucketString>::Type        TFibreBucketString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type       TFibreSuperBucketString;
    typedef typename Value<TFibreBitString>::Type               TFibreBitStringValue;
    typedef typename Value<TFibreBucketString>::Type              TFibreBucketStringValue;
    typedef typename Value<TFibreSuperBucketString>::Type             TFibreSuperBucketStringValue;

    TFibreSuperBucketStringValue const buPos = getBuPos_(bitString, pos);
    TFibreBitStringValue const posInBucket = getPosInBu_(bitString, pos);
    TFibreBitStringValue const one = 1;
    TFibreBitStringValue const shiftValue = one << posInBucket;
    if (!setBit)
    {
        bitString.bString[buPos] &= ~(shiftValue);
        return;
    }
    bitString.bString[buPos] |= shiftValue;
}


// This function checks if the specified position corresponds to the first
// position in a new bucket in the current superBucket. If this is the case the former
// bucket values have to be copied.
template <typename TSpec, typename TPos, typename TBitsPerBucket>
inline void takeBuValueAlong_(RankSupportBitString<TSpec> & bitString, TPos const pos, TBitsPerBucket const bpb)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type       TFibreSuperBucketString;
    typedef typename Value<TFibreSuperBucketString>::Type             TFibreSuperBucketStringValue;

    TFibreSuperBucketStringValue buPos = getBuPos_(bitString, pos);
    if (buPos)
    {
        if (!(pos % bpb) && ((pos + bpb) % (bpb * bpb)))
        {
            getFibre(bitString, FibreBucketString())[buPos + 1] = getFibre(bitString, FibreBucketString())[buPos];
        }
    }
}

// This function checks if the specified position corresponds to the first
// position in a new superBucket. If this is the case the former
// superBucket values have to be copied.
template <typename TSpec, typename TPos, typename TBitsPerBucket>
inline void takeSBuValueAlong_(RankSupportBitString<TSpec> & bitString, TPos const pos, TBitsPerBucket const bpb)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBucketString>::Type       TFibreSuperBucketString;
    typedef typename Value<TFibreSuperBucketString>::Type             TFibreSuperBucketStringValue;

    TFibreSuperBucketStringValue sBuPos = getSBuPos_(bitString, pos);
    if (sBuPos)
    {
        if (!(pos % (bpb * bpb)))
        {
            getFibre(bitString, FibreSuperBucketString())[sBuPos + 1] = getFibre(bitString, FibreSuperBucketString())[sBuPos];
        }
    }
}

/*template <typename TSpec, typename TBit>
inline void appendBitOnly(RankSupportBitString<TSpec> & rankSupportBitString, TBit bit)
{
    setBit(rankSupportBitString, length(rankSupportBitString), bit);
    ++rankSupportBitString.length;
}

template <typename TValue>
inline unsigned getRankInBucket(const TValue value)
{
    return getRankInBucket(value, typename Eval < (BitsPerValue<TValue>::VALUE > 32) > ::Type());
}

template <typename TValue>
inline unsigned getRankInBucket(const TValue value, False)
{
    return __builtin_popcountl(static_cast<int32_t>(value));
}

template <typename TValue>
inline unsigned getRankInBucket(const TValue value, True)
{
    return __builtin_popcountll(static_cast<int64_t>(value));
}

template <typename TValue, typename TPos>
inline TValue getRankInBucket(const String<TValue> & bString, const TPos pos)
{

    unsigned short const bitsPerValue = BitsPerValue<TValue>::VALUE;
    TValue const one = -1;
   // TValue const mask = one >> (bitsPerValue - (pos % bitsPerValue) - 1);
    TValue const mask = one >> (bitsPerValue - (pos & (bitsPerValue - 1)) - 1);
    return getRankInBucket(bString[pos / bitsPerValue] & mask);
}

template <typename TSpec>
inline void printBits(const TSpec entrie, const int blocks);

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type
getRank(const RankSupportBitString<TSpec> & rankSupportBitString, const TPos pos)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type        TSuperBucketString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;

    typename Value<TSuperBucketString>::Type sum = 0;
    unsigned bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;
    unsigned long bucketPos = pos / bitsPerBucket;
    if (bucketPos)
    {
        sum += rankSupportBitString.bucketString[bucketPos];
        unsigned superBucketPos = (bucketPos - 1) / bitsPerBucket;
        if (superBucketPos)
        {
            --superBucketPos;
            sum += rankSupportBitString.sBuString[superBucketPos];
        }
    }

    sum += getRankInBucket(rankSupportBitString.bString, pos);
    return sum;
}

template <typename TBitString>
inline void completeRankSupportBitString(TBitString & bString){}

template <typename TSpec>
//, typename TRankSupportBitString>
inline void completeRankSupportBitString(RankSupportBitString<TSpec> & bString)
{
    if (length(bString))
    {
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBucketString>::Type         TBucketString;
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type        TSuperBucketString;

        typedef typename Value<TBucketString>::Type TBucketValue;
        typedef typename Value<TSuperBucketString>::Type TSuperBucketValue;

        resize(bString, length(bString));

        unsigned bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;
        TSuperBucketValue superBucketCounter = 0;

        TBucketValue tempSum = 0;
        TBucketValue bucketSum = 0;
        TSuperBucketValue superBucketSum = 0;

        TBitString & bString_ = bString.bString;
        TBucketString & bucketString = bString.bucketString;
        TSuperBucketString & sBuString = bString.sBuString;

        typedef typename Value<TSuperBucketString>::Type TSize;
        for (TSize i = 0; i < length(bString_) - 1; i++)
        {
            tempSum = getRankInBucket(bString_[i]);
            bucketSum += tempSum;
            bucketString[i + 1] = bucketSum;
            if (!((i + 1) & (bitsPerBucket - 1)))
            {
                superBucketSum += bucketSum;
                sBuString[superBucketCounter] = superBucketSum;
                bucketSum = 0;
                ++superBucketCounter;
            }
        }
    }
}

//Manuel Forwards
template <typename TText, typename TSpec>
struct WaveletTree;

template <typename TChar, typename TPointer, typename TSpec>
struct WaveletTreeStructure;

struct FibreSplitValues_;
typedef Tag<FibreSplitValues_> const FibreSplitValues;
//template <typename TText>
//struct WaveletTreeStructure;

struct FibreBitStrings_;
typedef Tag<FibreBitStrings_> const FibreBitStrings;

//fills ONE bit string with 0 for characters that will appear in the left branch and
//1 for characters that will appear in the right branch
//template < typename TBitString, typename TCharacterValue, typename TPosInSubTree, typename TText, typename TWaveletTreeSpec >//, typename TRankSupportBitString>
template <typename TBitString, typename TText>
//, typename TRankSupportBitString>
inline void fillBitString(
    const unsigned lowestValue,
    const unsigned splitValue,
    const unsigned highestValue,
    TBitString & bString,        //TRankSupportBitString &counterBitString,
    const TText & text)
{
    unsigned short character;
    unsigned long long pos = 0;
    if (length(bString) == 0)
    {
        return;
    }
    for (unsigned i = 0; i < length(text); ++i)
    {
        character = text[i];
        if (character >= lowestValue && character <= highestValue)
        {
            if (character >= splitValue)
            {
                setBit(bString, pos, 1);
            }
            ++pos;
        }
    }
    completeRankSupportBitString(bString);
}

template <typename TCharacterValue, typename TWaveletTreeSpec, typename TText>
//, typename TRankSupportBitString>
inline void fillBitString(
    WaveletTree<TText, TWaveletTreeSpec> & waveletTree,
    typename Iterator<typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type>::Type & iter,
    const TText & text,
    const TCharacterValue lowerBound,
    const TCharacterValue upperBound)
{
    TCharacterValue lowestValue = lowerBound;
    TCharacterValue splitValue = iter.waveletTreeStructure->treeNodes[iter.position].i1;
    TCharacterValue highestValue = upperBound;
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreBitStrings>::Type TBitStrings;
    typedef typename Value<TBitStrings>::Type TBitString;
    TBitString & bString = waveletTree.bStrings[iter.position];

    TCharacterValue character;
    unsigned long long pos = 0;
    //	TBitString &bString = tree.bStrings[treePos];
    if (length(bString) == 0)
    {
        return;
    }
    for (unsigned i = 0; i < length(text); ++i)
    {
        character = text[i];
        if ((character >= lowestValue) && (character <= highestValue))
        {
            if (character >= splitValue)
            {
                setBit(bString, pos, 1);
            }
            ++pos;
        }
    }
    completeRankSupportBitString(bString);
}

template <typename TText>
inline unsigned short nearestPowOfTwo(TText & text)
{
    unsigned l = length(text);
    --l;
    for (unsigned i = 1; i <= 64; i *= 2)
    {
        l |= (l >> i);
    }
    ++l;
    return l;
}

template <typename TSpec>
inline bool open(
    RankSupportBitString<TSpec> & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1);

    String<char> name;
    name = fileName;    append(name, ".bit");
    if (!open(getFibre(string, FibreRankSupportBitString()), toCString(name), openMode))
        return false;

    name = fileName;    append(name, ".bucket");    open(getFibre(string, FibreRankSupportBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".sbucket");   open(getFibre(string, FibreRankSupportSuperBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".length");    open(lengthString, toCString(name), openMode);
    string.length = lengthString[0];
    return true;
}

template <typename TSpec>
inline bool open(
    RankSupportBitString<TSpec> & string,
    const char * fileName)
{
    return open(string, fileName, OPEN_RDONLY);
}

// ATTENTION:
// This implementation of open doesn't work with external memory StringSets (External<>, MMap<>)
// If you need a persistent external StringSet you have to use a Owner<ConcatDirect<> > StringSet.
template <typename TSpec>
inline bool open(StringSet<RankSupportBitString<TSpec> > & multi, const char * fileName, int openMode)
{
    SEQAN_CHECKPOINT

    typedef typename Size<RankSupportBitString<TSpec> >::Type TSize;
    typedef String<TSize> TSizeString;
    TSizeString sizeString;
    resize(sizeString, length(multi));
    CharString name = fileName;
    name = fileName;    append(name, ".ssize"); open(sizeString, toCString(name), openMode);

    char id[12];     // 2^32 has 10 decimal digits + 1 (0x00)
    unsigned i = 0;
    clear(multi);
    while (true)
    {
        sprintf(id, ".%u", i);
        name = fileName;
        append(name, id);
        {
            resize(multi, i + 1);
            if (!open(multi[i], toCString(name), (openMode & ~OPEN_CREATE) | OPEN_QUIET))
            {
                resize(multi, i);
                break;
            }
        }
        ++i;
    }

    for (TSize i = 0; i < length(multi); ++i)
    {
        multi[i].length = sizeString[i];
    }

    return i > 1;
}

template <typename TSpec>
inline bool open(
    StringSet<RankSupportBitString<TSpec> > & strings,
    const char * fileName)
{
    return open(strings, fileName, OPEN_RDONLY);
}

template <typename TSpec>
inline bool save(
    RankSupportBitString<TSpec> const & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1);
    lengthString[0] = string.length;
    String<char> name;
    name = fileName;    append(name, ".length");    save(lengthString, toCString(name), openMode);
    name = fileName;    append(name, ".bit");       save(getFibre(string, FibreRankSupportBitString()), toCString(name), openMode);
    name = fileName;    append(name, ".bucket");    save(getFibre(string, FibreRankSupportBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".sbucket");   save(getFibre(string, FibreRankSupportSuperBucketString()), toCString(name), openMode);
    return true;
}

template <typename TSpec>
inline bool save(
    RankSupportBitString<TSpec> const & string,
    const char * fileName)
{
    return save(string, fileName, DefaultOpenMode<RankSupportBitString<TSpec> >::VALUE);
}

template <typename TSpec>
inline bool save(StringSet<RankSupportBitString<TSpec> > const & multi, const char * fileName, int openMode)
{
    SEQAN_CHECKPOINT

    typedef typename Size<RankSupportBitString<TSpec> >::Type TSize;
    typedef String<TSize> TSizeString;
    TSizeString sizeString;
    resize(sizeString, length(multi));
    for (TSize i = 0; i < length(multi); ++i)
    {
        sizeString[i] = length(multi[i]);
    }

    CharString name;
    name = fileName;    append(name, ".ssize"); save(sizeString, toCString(name), openMode);
    if (length(multi) == 0) return true;

    char id[12];     // 2^32 has 10 decimal digits + 2 ('.' and 0x00)
    for (unsigned i = 0; i < length(multi); ++i)
    {
        sprintf(id, ".%u", i);
        name = fileName;
        append(name, &(id[0]));
        if (!save(multi[i], toCString(name), openMode))
            return false;
    }
    return true;
}

template <typename TSpec>
inline bool save(
    StringSet<RankSupportBitString<TSpec> > const & strings,
    const char * fileName)
{
    return save(strings, fileName, OPEN_RDONLY);
}
*/

/*    template <typename TPos>
    inline typename Value<TSuperBucketString>::Type const operator[](TPos pos)
    {
        unsigned const bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;
        return bString[pos / bitsPerBucket];
    }
*/


/*    template <typename TText>
    RankSupportBitString(TText & text) :
        bString(),
        butString(),
        sButString(),
        length_(length(text))
    {
        unsigned const bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;

        TSuperBucketString textLength = length(text);
        (length(text) % bitsPerBucket == 0) ? resize(bString, textLength / bitsPerBucket) : resize(bString, textLength / bitsPerBucket + 1);

        unsigned const bitsPerBucketStringEntrie = BitsPerValue<TBucketString>::VALUE;
        resize(bucketString, length(bString) / bitsPerBucketStringEntrie);

        unsigned const bitsPerSuperBucketStringEntrie = bitsPerBucketStringEntrie * bitsPerBucketStringEntrie;
        resize(sBuString, length(bucketString) / bitsPerSuperBucketStringEntrie);
    }
*/


/*template <typename TValue>
inline void printBits(TValue entrie)
{
    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
    TValue one = 1;
    std::cerr << "entrie: " << entrie << std::endl;
    std::cerr << bitsPerValue << std::endl;
    for (TValue i = 0; i < bitsPerValue; ++i)
    {
        std::cerr << ((entrie >> i) & one);
    }
    std::cerr << std::endl;
}
*/

/*
template <typename TValue, typename TSize>
inline std::ostream & printBits(std::ostream & stream, TValue entrie, TSize blockSize)
{
    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
    bool temp;
    for (int i = bitsPerValue - 1; i >= 0; --i)
    {
        temp = (entrie >> i) & 1;
        stream << temp;
        if ((bitsPerValue - i) % blockSize == 0)
            stream << " ";
    }
    return stream;
}
*/

/*
template <typename TSpec>
inline std::ostream & operator<<(std::ostream & stream, const RankSupportBitString<TSpec> & rankSupportBitString)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type                TBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBucketString>::Type             TBucketString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportSuperBucketString>::Type        TSuperBucketString;

    typedef typename Value<TBitString>::Type          TBitStringValue;
    typedef typename Value<TBucketString>::Type       TBucketStringValue;
    typedef typename Value<TSuperBucketString>::Type  TSuperBucketStringValue;

    unsigned bitsPerBucket = BitsPerValue<typename Value<TBitString>::Type>::VALUE;

    TBitString const & bString = rankSupportBitString.bString;
    TBucketString const & bucketString = rankSupportBitString.bucketString;
    TSuperBucketString const & sBuString = rankSupportBitString.sBuString;

    stream << "  ";
    for (TBitStringValue i = 0; i < length(bString); i++)
    {
        printBits(stream, bString[i], bitsPerBucket);
    }
    stream << std::endl;

    for (TBucketStringValue i = 0; i < length(bucketString); i++)
    {
        stream << bucketString[i] << " ";
    }
    stream << std::endl;

    for (TSuperBucketStringValue i = 0; i < length(sBuString); i++)
    {
        stream << sBuString[i] << " ";
    }
    return stream;
    //	return print(stream, rankSupportBitString);
}
*/
/*template <typename TSpec, typename TSize>
inline void resize(RankSupportBitString<TSpec> & bitString, TSize size)
{
    reserve(rankSupportBitString, size);
    rankSupportBitString.length = size;
}

template <typename TSpec, typename TSize, typename TValue>
inline void resize(RankSupportBitString<TSpec> & rankSupportBitString, TSize size, TValue value)
{
    rankSupportBitString.length = size;

    typedef typename Fibre<RankSupportBitString<TSpec>, FibreRankSupportBitString>::Type    TBitString;
    typedef typename Value<TBitString>::Type TBitStringValue;
    unsigned bitsPerBucket = BitsPerValue<TBitStringValue>::VALUE;
    unsigned long long numberOfBuckets;
    (size) ? numberOfBuckets = size / bitsPerBucket + 1 : numberOfBuckets = 0;

    resize(rankSupportBitString.bString, numberOfBuckets, value);
    resize(rankSupportBitString.bucketString, numberOfBuckets, value);
    resize(rankSupportBitString.sBuString, numberOfBuckets / bitsPerBucket + 1, value);
}
*/



}


#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_BITSTRING_BETA_H_
