/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  Author: David Weese <david.weese@fu-berlin.de>
 ============================================================================
  Some useful functors for modified strings.
 ==========================================================================*/

#ifndef SEQAN_MODIFIER_MODIFIER_FUNCTORS_H_
#define SEQAN_MODIFIER_MODIFIER_FUNCTORS_H_

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// text transformation
//////////////////////////////////////////////////////////////////////////////

/**
.Class.FunctorUpcase:
..cat:Modifier
..summary:Functor that returns the upper case character to a given character.
..signature:FunctorUpcase<TValue>
..param.TValue:The input value type.
..remarks:This Functor is a derivation of the STL unary function.
..include:seqan/modifier.h
*/
template <typename InType, typename Result = InType>
struct FunctorUpcase : public ::std::unary_function<InType,Result> 
{
    inline Result operator()(InType x) const {
        SEQAN_CHECKPOINT;
        if (('a' <= x) && (x <= 'z')) return (x + ('A' - 'a'));
        return x; 
    }
};


/**
.Class.FunctorLowcase:
..cat:Modifier
..summary:Functor that returns the lower case character to a given character.
..signature:FunctorLowcase<TValue>
..param.TValue:The input value type.
..remarks:This Functor is a derivation of the STL unary function.
..include:seqan/modifier.h
*/
template <typename InType, typename Result = InType>
struct FunctorLowcase : public ::std::unary_function<InType,Result> 
{
    inline Result operator()(InType x) const {
        SEQAN_CHECKPOINT;
        if (('A' <= x) && (x <= 'Z')) return (x + ('a' - 'A'));
        return x; 
    }
};


//////////////////////////////////////////////////////////////////////////////
// alphabet transformation
//////////////////////////////////////////////////////////////////////////////

/**
.Class.FunctorConvert:
..cat:Modifier
..summary:Functor that converts a $TInValue$ type to a $TOutValue$ type character.
..signature:FunctorConvert<TInValue, TOutValue>
..param.TInValue:The input value type.
..param.TOutValue:The output value type.
..remarks:This Functor is a derivation of the STL unary function.
..include:seqan/modifier.h
*/
template <typename InType, typename OutType>
struct FunctorConvert : public ::std::unary_function<InType,OutType> 
{
    inline OutType operator()(InType x) const {
        SEQAN_CHECKPOINT;
        return x; 
    }
};


//////////////////////////////////////////////////////////////////////////////
// DNA/RNA complement
//////////////////////////////////////////////////////////////////////////////

/**
.Class.FunctorComplement:
..cat:Modifier
..summary:Functor that returns the complement nucleotide to a given nucleotide.
..signature:FunctorComplement<TValue>
..param.TValue:The input value type.
...type:Spec.Dna
...type:Spec.Dna5
...type:Spec.Rna
...type:Spec.Rna5
..remarks:This Functor is a derivation of the STL unary function.
..include:seqan/modifier.h
*/
template <typename TValue>
struct FunctorComplement;


template <typename T = void>
struct TranslateTableDna5ToDna5Complement_
{
    static char const VALUE[5];
};


template <typename T>
char const TranslateTableDna5ToDna5Complement_<T>::VALUE[5] = {'T', 'G', 'C', 'A', 'N'};


template <typename T = void>
struct TranslateTableRna5ToRna5Complement_
{
    static char const VALUE[5];
};


template <typename T>
char const TranslateTableRna5ToRna5Complement_<T>::VALUE[5] = {'U', 'G', 'C', 'A', 'N'};


template <>
struct FunctorComplement<char> : public ::std::unary_function<Dna5,Dna5> 
{
    inline Dna5 operator()(Dna5 x) const {
        SEQAN_CHECKPOINT;
        return TranslateTableDna5ToDna5Complement_<>::VALUE[x.value]; 
    }
};

template <>
struct FunctorComplement<Dna> : public ::std::unary_function<Dna,Dna> 
{
    inline Dna operator()(Dna x) const {
        SEQAN_CHECKPOINT;
        return TranslateTableDna5ToDna5Complement_<>::VALUE[x.value]; 
    }
};


template <>
struct FunctorComplement<Dna5> : public ::std::unary_function<Dna5,Dna5> 
{
    inline Dna5 operator()(Dna5 x) const {
        SEQAN_CHECKPOINT;
        return TranslateTableDna5ToDna5Complement_<>::VALUE[x.value]; 
    }
};


template <>
struct FunctorComplement<Rna> : public ::std::unary_function<Rna,Rna> 
{
    inline Rna operator()(Rna x) const {
        SEQAN_CHECKPOINT;
        return TranslateTableRna5ToRna5Complement_<>::VALUE[x.value]; 
    }
};


template <>
struct FunctorComplement<Rna5> : public ::std::unary_function<Rna5,Rna5> 
{
    inline Dna5 operator()(Rna5 x) const {
        SEQAN_CHECKPOINT;
        return TranslateTableRna5ToRna5Complement_<>::VALUE[x.value]; 
    }
};


template <>
struct FunctorComplement<DnaQ> : public ::std::unary_function<DnaQ,DnaQ> 
{
    inline DnaQ operator()(DnaQ x) const {
        SEQAN_CHECKPOINT;
        int qual = getQualityValue(x);
        x = TranslateTableDna5ToDna5Complement_<>::VALUE[ordValue((Dna)x)];
        assignQualityValue(x, qual);
        return x;
    }
};


template <>
struct FunctorComplement<Dna5Q> : public ::std::unary_function<Dna5Q,Dna5Q> 
{
    inline Dna5Q operator()(Dna5Q x) const {
        SEQAN_CHECKPOINT;
        int qual = getQualityValue(x);
        x = TranslateTableDna5ToDna5Complement_<>::VALUE[ordValue((Dna5)x)];
        assignQualityValue(x, qual);
        return x;
    }
};

}  // namespace seqan

#endif  // SEQAN_MODIFIER_MODIFIER_FUNCTORS_H_
