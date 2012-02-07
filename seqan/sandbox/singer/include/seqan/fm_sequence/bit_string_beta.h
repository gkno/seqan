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
struct Rsbs;

struct FibreRsbsBS_;
struct FibreRsbsBuS_;
struct FibreRsbsSBuS_;

typedef Tag<FibreRsbsBS_> const RsbsBS;
typedef Tag<FibreRsbsBuS_> const RsbsBuS;
typedef Tag<FibreRsbsSBuS_> const RsbsSBuS;

template <typename TSpec>
struct Fibre<Rsbs<TSpec>, RsbsBS>
{
    typedef String<unsigned long> Type;
};

template <typename TSpec>
struct Fibre<Rsbs<TSpec>, RsbsBuS>
{
    typedef String<unsigned short> Type;
};

template <typename TSpec>
struct Fibre<Rsbs<TSpec>, RsbsSBuS>
{
    typedef String<unsigned int> Type;
};

//The limiting factor of the size is the underlying data type of the super bucket string
template <typename TSpec>
struct Size<Rsbs<TSpec> >
{
    typedef typename Value<typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type>::Type Type;
};

template <typename TSpec>
struct Rsbs
{
    typedef typename Fibre<Rsbs, RsbsBS>::Type      TBitString;
    typedef typename Fibre<Rsbs, RsbsBuS>::Type     TBucketString;
    typedef typename Fibre<Rsbs, RsbsSBuS>::Type    TSuperBucketString;

    TBitString                  bString;
    TBucketString               buString;
    TSuperBucketString          sBuString;
    typename Size<Rsbs>::Type   length_;

    Rsbs() :
        bString(),
        buString(),
        sBuString(),
        length_(0)
    {}

    template <typename TString>
    Rsbs(TString const & input) :
        bString(),
        buString(),
        sBuString(),
        length_(length(input))
    {
        typedef Rsbs<TSpec>                                     TRsbs;
        typedef typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type     TRsbsSBuS;
        typedef typename Value<TRsbsSBuS>::Type                 TRsbsSBuSValue;

        resize(*this, length_);
        for (TRsbsSBuSValue i = 0; i < length(*this); ++i)
        {
            setBit(*this, i, input[i]);
        }
        completeRsbs(*this);
    }

    inline Rsbs & operator=(Rsbs const & other)
    {
        bString = other.bString;
        buString = other.buString;
        sBuString = other.sBuString;
        length_ = other.length_;
        return *this;
    }

    inline bool operator==(const Rsbs & other) const
    {
        return length_ == other.length_ &&
               bString == other.bString &&
               buString == other.buString &&
               sBuString == other.sBuString;
    }

};

template <typename TSpec>
inline typename Fibre<Rsbs<TSpec>, RsbsBS>::Type &
getFibre(Rsbs<TSpec>&string, const RsbsBS)
{
    return string.bString;
}

template <typename TSpec>
inline typename Fibre<Rsbs<TSpec>, RsbsBS>::Type const &
getFibre(Rsbs<TSpec> const & string, const RsbsBS)
{
    return string.bString;
}

template <typename TSpec>
inline typename Fibre<Rsbs<TSpec>, RsbsBuS>::Type &
getFibre(Rsbs<TSpec>&string, const RsbsBuS)
{
    return string.buString;
}

template <typename TSpec>
inline typename Fibre<Rsbs<TSpec>, RsbsBuS>::Type const &
getFibre(Rsbs<TSpec> const & string, const RsbsBuS)
{
    return string.buString;
}

template <typename TSpec>
inline typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type &
getFibre(Rsbs<TSpec>&string, const RsbsSBuS)
{
    return string.sBuString;
}

template <typename TSpec>
inline typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type const &
getFibre(const Rsbs<TSpec>&string, const RsbsSBuS)
{
    return string.sBuString;
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<Rsbs<TSpec>, RsbsBuS>::Type>::Type
getPosInBu_(Rsbs<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<TRsbs, RsbsBS>::Type     TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type   TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type           TRsbsBSValue;
    typedef typename Value<TRsbsSBuS>::Type         TRsbsSBuSValue;

    TRsbsSBuSValue bitsPerValue = BitsPerValue<TRsbsBSValue>::VALUE;
    return pos % bitsPerValue;
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type>::Type
getBuPos_(Rsbs<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef Rsbs<TSpec>                             TRsbs;
    typedef typename Fibre<TRsbs, RsbsBS>::Type     TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type   TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type           TRsbsBSValue;
    typedef typename Value<TRsbsSBuS>::Type         TRsbsSBuSValue;

    TRsbsSBuSValue const bitsPerValue_ = BitsPerValue<TRsbsBSValue>::VALUE;
    return pos / bitsPerValue_;
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type>::Type
getSBuPos_(Rsbs<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef Rsbs<TSpec>                             TRsbs;
    typedef typename Fibre<TRsbs, RsbsBS>::Type     TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type   TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type           TRsbsBSValue;
    typedef typename Value<TRsbsSBuS>::Type         TRsbsSBuSValue;

    TRsbsSBuSValue const bitsPerValue_ = BitsPerValue<TRsbsBSValue>::VALUE;

    return pos / (bitsPerValue_ * bitsPerValue_);
}

template <typename TSpec, typename TSize>
inline void reserve(Rsbs<TSpec> & bitString, TSize const size)
{
    typedef Rsbs<TSpec>                             TRsbs;
    typedef typename Fibre<TRsbs, RsbsBS>::Type     TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsBuS>::Type    TRsbsBuS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type   TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type           TRsbsBSValue;
    typedef typename Value<TRsbsBuS>::Type          TRsbsBuSValue;
    typedef typename Value<TRsbsSBuS>::Type         TRsbsSBuSValue;

    if (size)
    {
        TRsbsBuSValue bitsPerBucket_ = BitsPerValue<TRsbsBSValue>::VALUE;
        TRsbsSBuSValue numberOfBuckets_ = size / bitsPerBucket_;

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

template <typename TSpec, typename TSize>
inline void resize(Rsbs<TSpec> & bitString, TSize const size)
{
    reserve(bitString, size);
    bitString.length_ = size;
}

template <typename TSpec>
inline typename Value<typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type>::Type
length(Rsbs<TSpec> const & bitString)
{
    return bitString.length_;
}

template <typename TSpec>
inline typename Value<typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type>::Type
length(Rsbs<TSpec> & bitString)
{
    return bitString.length_;
}

template <typename TSpec>
inline void clear(Rsbs<TSpec> & bitString)
{
    clear(bitString.bString);
    clear(bitString.buString);
    clear(bitString.sBuString);
}

template <typename TSpec, typename TPos>
inline bool getBit(Rsbs<TSpec> & bitString, TPos const pos)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<TRsbs, RsbsBS>::Type         TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsBuS>::Type        TRsbsBuS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type               TRsbsBSValue;
    typedef typename Value<TRsbsBuS>::Type              TRsbsBuSValue;
    typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

    TRsbsBuSValue const bitsPerValue = BitsPerValue<TRsbsBSValue>::VALUE;
    TRsbsBSValue const one = 1;
    TRsbsBSValue const shiftValue = pos % bitsPerValue;
    TRsbsBSValue const buPos = getBuPos_(bitString, pos);
    return (bitString.bString[buPos] >> shiftValue) & one;
}

template <typename TSpec, typename TPos>
inline bool getBit(Rsbs<TSpec> const & bitString, TPos const pos)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<Rsbs<TSpec>, RsbsBS>::Type   TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsBuS>::Type        TRsbsBuS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type               TRsbsBSValue;
    typedef typename Value<TRsbsBuS>::Type              TRsbsBuSValue;
    typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

    TRsbsBuSValue const bitsPerValue = BitsPerValue<TRsbsBSValue>::VALUE;
    TRsbsBSValue const one = 1;
    TRsbsBSValue const shiftValue = pos % bitsPerValue;
    TRsbsBSValue const buPos = getBuPos_(bitString, pos);
    return (bitString.bString[buPos] >> shiftValue) & one;
}

template <typename TSpec, typename TPos, typename TBit>
inline void setBit(Rsbs<TSpec> & bitString, TPos const pos, TBit const setBit)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<Rsbs<TSpec>, RsbsBS>::Type   TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsBuS>::Type        TRsbsBuS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type               TRsbsBSValue;
    typedef typename Value<TRsbsBuS>::Type              TRsbsBuSValue;
    typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

    TRsbsSBuSValue const buPos = getBuPos_(bitString, pos);
    TRsbsBSValue const posInBucket = getPosInBu_(bitString, pos);
    TRsbsBSValue const one = 1;
    TRsbsBSValue const shiftValue = one << posInBucket;
    if (!setBit)
    {
        bitString.bString[buPos] &= ~(shiftValue);
        return;
    }
    bitString.bString[buPos] |= shiftValue;
}

/*
 * This function checks if the specified position corresponds to the first
 * position in a new bucket in the current superBucket. If this is the case the former
 * bucket values have to be copied.
 */
template <typename TSpec, typename TPos, typename TBitsPerBucket>
inline void takeBuValueAlong_(Rsbs<TSpec> & bitString, TPos const pos, TBitsPerBucket const bpb)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
    typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

    TRsbsSBuSValue buPos = getBuPos_(bitString, pos);
    if (buPos)
    {
        if (!(pos % bpb) && ((pos + bpb) % (bpb * bpb)))
        {
            getFibre(bitString, RsbsBuS())[buPos + 1] = getFibre(bitString, RsbsBuS())[buPos];
        }
    }
}

/*
 * This function checks if the specified position corresponds to the first
 * position in a new superBucket. If this is the case the former
 * superBucket values have to be copied.
 */
template <typename TSpec, typename TPos, typename TBitsPerBucket>
inline void takeSBuValueAlong_(Rsbs<TSpec> & bitString, TPos const pos, TBitsPerBucket const bpb)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
    typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

    TRsbsSBuSValue sBuPos = getSBuPos_(bitString, pos);
    if (sBuPos)
    {
        if (!(pos % (bpb * bpb)))
        {
            getFibre(bitString, RsbsSBuS())[sBuPos + 1] = getFibre(bitString, RsbsSBuS())[sBuPos];
        }
    }
}

template <typename TSpec, typename TBit>
inline void append(Rsbs<TSpec> & bitString, TBit const bit)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<Rsbs<TSpec>, RsbsBS>::Type   TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsBuS>::Type        TRsbsBuS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type               TRsbsBSValue;
    typedef typename Value<TRsbsBuS>::Type              TRsbsBuSValue;
    typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

    TRsbsBuSValue const bitsPerValue_ = BitsPerValue<TRsbsBSValue>::VALUE;
    TRsbsSBuSValue const length_ = length(bitString);

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
        TRsbsSBuSValue buPos = getBuPos_(bitString, length_);
        if (((buPos + 1) % (bitsPerValue_)))
        {
            ++getFibre(bitString, RsbsBuS())[buPos + 1];
        }
        TRsbsSBuSValue sBuPos = getSBuPos_(bitString, length_);
        ++getFibre(bitString, RsbsSBuS())[sBuPos + 1];
    }
    ++bitString.length_;
}

template <typename TValue>
inline unsigned getRankInBucket_(const TValue value, False)
{
    return __builtin_popcountl(static_cast<int32_t>(value));
}

template <typename TValue>
inline unsigned getRankInBucket_(const TValue value, True)
{
    return __builtin_popcountll(static_cast<int64_t>(value));
}

template <typename TValue>
inline unsigned getRankInBucket_(const TValue value)
{
    return getRankInBucket_(value, typename Eval<(BitsPerValue<TValue>::VALUE > 32)>::Type());
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type>::Type
getRankInBucket_(Rsbs<TSpec> & bitString, TPos const pos)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<Rsbs<TSpec>, RsbsBS>::Type   TRsbsBS;
    typedef typename Fibre<TRsbs, RsbsBuS>::Type        TRsbsBuS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
    typedef typename Value<TRsbsBS>::Type               TRsbsBSValue;
    typedef typename Value<TRsbsBuS>::Type              TRsbsBuSValue;
    typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

    TRsbsBuSValue const posInBu = getPosInBu_(bitString, pos);
    TRsbsSBuSValue const buPos = getBuPos_(bitString, pos);
    TRsbsBuSValue const bitsPerValue = BitsPerValue<TRsbsBSValue>::VALUE;
    TRsbsBSValue const one = -1;
    TRsbsBSValue const mask = one >> (bitsPerValue - posInBu - 1);
    //std::cerr << "one: " << one << " mask1: " << mask << " " << (getRankInBucket_(bitString.bString[buPos] & mask)) <<  std::endl;
    return getRankInBucket_(bitString.bString[buPos] & mask);
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<Rsbs<TSpec>, RsbsSBuS>::Type>::Type
getRank(Rsbs<TSpec> & bitString, TPos const pos)
{
    typedef Rsbs<TSpec>                                 TRsbs;
    typedef typename Fibre<TRsbs, RsbsBuS>::Type        TRsbsBuS;
    typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
    typedef typename Value<TRsbsBuS>::Type              TRsbsBuSValue;
    typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

    TRsbsSBuSValue const buPos = getBuPos_(bitString, pos);
    TRsbsSBuSValue const sBuPos = getSBuPos_(bitString, pos);

    //std::cerr << "buPos: " <<  buPos << " " << sBuPos << std::endl;
//	std::cerr << "getRankInBucket_(bitString, pos): " <<  getRankInBucket_(bitString, pos) << " " <<  bitString.buString[buPos] << " " << bitString.sBuString[sBuPos] << std::endl;


    return getRankInBucket_(bitString, pos)
           + bitString.buString[buPos]
           + bitString.sBuString[sBuPos];
}

template <typename TSpec>
inline void completeRsbs(Rsbs<TSpec> & bitString)
{
    if (length(bitString))
    {

        typedef Rsbs<TSpec>                                 TRsbs;
        typedef typename Fibre<Rsbs<TSpec>, RsbsBS>::Type   TRsbsBS;
        typedef typename Fibre<TRsbs, RsbsBuS>::Type        TRsbsBuS;
        typedef typename Fibre<TRsbs, RsbsSBuS>::Type       TRsbsSBuS;
        typedef typename Value<TRsbsBS>::Type               TRsbsBSValue;
        typedef typename Value<TRsbsBuS>::Type              TRsbsBuSValue;
        typedef typename Value<TRsbsSBuS>::Type             TRsbsSBuSValue;

        TRsbsBuSValue buSum_ = 0;
        TRsbsSBuSValue sBuSum_ = 0;
        TRsbsSBuSValue superBucketCounter = 1;
        TRsbsBuSValue const bitsPerValue = BitsPerValue<TRsbsBSValue>::VALUE;

        for (TRsbsSBuSValue i = 0; i < length(bitString.bString) - 1; ++i)
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

/*template <typename TSpec, typename TBit>
inline void appendBitOnly(Rsbs<TSpec> & rankSupportBitString, TBit bit)
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
inline typename Value<typename Fibre<Rsbs<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type
getRank(const Rsbs<TSpec> & rankSupportBitString, const TPos pos)
{
    typedef typename Fibre<Rsbs<TSpec>, FibreRankSupportSuperBucketString>::Type        TSuperBucketString;
    typedef typename Fibre<Rsbs<TSpec>, FibreRsbs>::Type                TBitString;

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
inline void completeRsbs(TBitString & bString){}

template <typename TSpec>
//, typename TRsbs>
inline void completeRsbs(Rsbs<TSpec> & bString)
{
    if (length(bString))
    {
        typedef typename Fibre<Rsbs<TSpec>, FibreRsbs>::Type                TBitString;
        typedef typename Fibre<Rsbs<TSpec>, FibreRankSupportBucketString>::Type         TBucketString;
        typedef typename Fibre<Rsbs<TSpec>, FibreRankSupportSuperBucketString>::Type        TSuperBucketString;

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
//template < typename TBitString, typename TCharacterValue, typename TPosInSubTree, typename TText, typename TWaveletTreeSpec >//, typename TRsbs>
template <typename TBitString, typename TText>
//, typename TRsbs>
inline void fillBitString(
    const unsigned lowestValue,
    const unsigned splitValue,
    const unsigned highestValue,
    TBitString & bString,        //TRsbs &counterBitString,
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
    completeRsbs(bString);
}

template <typename TCharacterValue, typename TWaveletTreeSpec, typename TText>
//, typename TRsbs>
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
    completeRsbs(bString);
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
    Rsbs<TSpec> & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<Rsbs<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1);

    String<char> name;
    name = fileName;    append(name, ".bit");
    if (!open(getFibre(string, FibreRsbs()), toCString(name), openMode))
        return false;

    name = fileName;    append(name, ".bucket");    open(getFibre(string, FibreRankSupportBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".sbucket");   open(getFibre(string, FibreRankSupportSuperBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".length");    open(lengthString, toCString(name), openMode);
    string.length = lengthString[0];
    return true;
}

template <typename TSpec>
inline bool open(
    Rsbs<TSpec> & string,
    const char * fileName)
{
    return open(string, fileName, OPEN_RDONLY);
}

// ATTENTION:
// This implementation of open doesn't work with external memory StringSets (External<>, MMap<>)
// If you need a persistent external StringSet you have to use a Owner<ConcatDirect<> > StringSet.
template <typename TSpec>
inline bool open(StringSet<Rsbs<TSpec> > & multi, const char * fileName, int openMode)
{
    SEQAN_CHECKPOINT

    typedef typename Size<Rsbs<TSpec> >::Type TSize;
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
    StringSet<Rsbs<TSpec> > & strings,
    const char * fileName)
{
    return open(strings, fileName, OPEN_RDONLY);
}

template <typename TSpec>
inline bool save(
    Rsbs<TSpec> const & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<Rsbs<TSpec>, FibreRankSupportSuperBucketString>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1);
    lengthString[0] = string.length;
    String<char> name;
    name = fileName;    append(name, ".length");    save(lengthString, toCString(name), openMode);
    name = fileName;    append(name, ".bit");       save(getFibre(string, FibreRsbs()), toCString(name), openMode);
    name = fileName;    append(name, ".bucket");    save(getFibre(string, FibreRankSupportBucketString()), toCString(name), openMode);
    name = fileName;    append(name, ".sbucket");   save(getFibre(string, FibreRankSupportSuperBucketString()), toCString(name), openMode);
    return true;
}

template <typename TSpec>
inline bool save(
    Rsbs<TSpec> const & string,
    const char * fileName)
{
    return save(string, fileName, DefaultOpenMode<Rsbs<TSpec> >::VALUE);
}

template <typename TSpec>
inline bool save(StringSet<Rsbs<TSpec> > const & multi, const char * fileName, int openMode)
{
    SEQAN_CHECKPOINT

    typedef typename Size<Rsbs<TSpec> >::Type TSize;
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
    StringSet<Rsbs<TSpec> > const & strings,
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
    Rsbs(TText & text) :
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
inline std::ostream & operator<<(std::ostream & stream, const Rsbs<TSpec> & rankSupportBitString)
{
    typedef typename Fibre<Rsbs<TSpec>, FibreRsbs>::Type                TBitString;
    typedef typename Fibre<Rsbs<TSpec>, FibreRankSupportBucketString>::Type             TBucketString;
    typedef typename Fibre<Rsbs<TSpec>, FibreRankSupportSuperBucketString>::Type        TSuperBucketString;

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
inline void resize(Rsbs<TSpec> & bitString, TSize size)
{
    reserve(rankSupportBitString, size);
    rankSupportBitString.length = size;
}

template <typename TSpec, typename TSize, typename TValue>
inline void resize(Rsbs<TSpec> & rankSupportBitString, TSize size, TValue value)
{
    rankSupportBitString.length = size;

    typedef typename Fibre<Rsbs<TSpec>, FibreRsbs>::Type    TBitString;
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
