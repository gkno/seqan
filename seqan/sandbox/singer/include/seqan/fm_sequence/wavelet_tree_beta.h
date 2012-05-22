// ==========================================================================
//                                  wavelet tree
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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

template <typename TSpec = void>
struct FmiDollarSubstituted;

//////////////////////////////////////////////////////////////////////////////
// WaveletTree fibres

/**
.Tag.WaveletTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.WaveletTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a WaveletTree.
..cat:WaveletTree

..tag.FibreBitStrings:The string set containing a bit string for each node.

..tag.FibreWaveletTreeStructure:The wavelet tree structure of the wavelet tree. 

..tag.FibreDollarPositions:The bit string encoding the position of the dollar sign.
...remarks:This fibre is only available if the wavelet tree is used as the 
occurrence table data structure of a FM index.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

///.Metafunction.Fibre.param.TSpec.type:Tag.WaveletTree Fibres
struct FibreBitStrings_;
struct FibreWaveletTreeStructure_;
struct FibreDollarPositions_;
struct FibreTreeNodes_;

template <typename TText, typename TSpec>
class WaveletTree;

// ==========================================================================
//Tags, Classes, Enums
// ==========================================================================

typedef Tag<FibreWaveletTreeStructure_> const FibreWaveletTreeStructure;
typedef Tag<FibreBitStrings_> const FibreBitStrings;
typedef Tag<FibreDollarPositions_> const FibreDollarPositions;
typedef Tag<FibreTreeNodes_> const FibreTreeNodes;

// ==========================================================================
//Metafunctions
// ==========================================================================

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>
{
    typedef StringSet<RankSupportBitString<void> > Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>
{
    typedef typename Value<TText>::Type TChar;
    typedef WaveletTreeStructure<TChar, void> Type;
};

template <typename TStringSpec, typename TSpec>
struct Fibre<WaveletTree<String<unsigned char, TStringSpec>, TSpec>, FibreWaveletTreeStructure>
{
    typedef WaveletTreeStructure<unsigned short, void> Type;
};

template <typename TText, typename TSpec>
struct Value<WaveletTree<TText, TSpec> >
{
    typedef typename Value<TText>::Type Type;
};

template <typename TText, typename TSpec>
struct Value<WaveletTree<TText, TSpec> const>
{
    typedef typename Value<TText>::Type const Type;
};

// ==========================================================================
//Tags, Classes, Enums
// ==========================================================================

/**
.Class.WaveletTree:
..cat:Graph
..summary:A wavelet tree is a tree like binary encoding of a text. 
..signature:WaveletTree<TText, TSpec>
..param.TText:The value type of the text.
..param.TSpec:The wavelet tree specialisation.
...tag:FmiDollarSubstituted
...default:void.
..include:seqan/index.h
*/
template <typename TText, typename TSpec = void>
class WaveletTree
{
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type TBitStrings;
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;

public:
    TBitStrings bitStrings;
    TWaveletTreeStructure waveletTreeStructure;

    WaveletTree() :
        bitStrings(),
        waveletTreeStructure()
    {}

    WaveletTree(TText const & text) :
        bitStrings(),
        waveletTreeStructure()
    {
        createWaveletTree(*this, text);
    }

    template <typename TFreqTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable) :
        bitStrings(),
        waveletTreeStructure()
    {
        createWaveletTree(*this, text, freqTable);
    }

    template <typename TFreqTable, typename TPrefixSumTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable, TPrefixSumTable const & prefixSumTable) :
        bitStrings(),
        waveletTreeStructure()
    {
        createWaveletTree(*this, text, freqTable, prefixSumTable);
    }

    inline WaveletTree & operator=(WaveletTree const & other)
    {
        bitStrings = other.bitStrings;
        waveletTreeStructure = other.waveletTreeStructure;
        return *this;
    }

    inline bool operator==(const WaveletTree & b) const
    {
        typedef typename Size<TText>::Type                            TSize;
        bool test = true;
        if (length(bitStrings) == length(b.bitStrings))
        {
            for (TSize i = 0; i < length(bitStrings); ++i)
            {
                if (!(bitStrings[i] == b.bitStrings[i]))
                {
                    test = false;
                }
            }
        }
        else
        {
            test = false;
        }

        return test &&
               waveletTreeStructure == b.waveletTreeStructure;
    }

};

template <typename TText, typename TSpec>
class WaveletTree<TText, FmiDollarSubstituted<TSpec> >:
    public WaveletTree<TText, void>
{
    typedef WaveletTree<TText, void> TBase;
    typedef typename Size<TText>::Type TSize;
    typedef typename Value<TText>::Type TChar;

public:
    TSize dollarPosition;
    TChar dollarSubstitute;

    WaveletTree() :
        TBase(),
        dollarPosition(),
        dollarSubstitute()
    {}

    WaveletTree(TText const & text) :
        TBase(text),
        dollarPosition(),
        dollarSubstitute()
    {}

    inline bool operator==(const WaveletTree & b) const
    {
        return static_cast<typename WaveletTree::TBase>(*this) == static_cast<typename WaveletTree::TBase>(b) &&
               dollarPosition == b.dollarPosition &&
               dollarSubstitute == b.dollarSubstitute;
    }

};


// ==========================================================================
//Functions
// ==========================================================================

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type &
getFibre(WaveletTree<TText, TSpec>&tree, const FibreBitStrings)
{
    return tree.bitStrings;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type const &
getFibre(const WaveletTree<TText, TSpec>&tree, const FibreBitStrings)
{
    return tree.bitStrings;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>::Type &
getFibre(WaveletTree<TText, TSpec>&tree, FibreWaveletTreeStructure)
{
    return tree.waveletTreeStructure;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>::Type const &
getFibre(WaveletTree<TText, TSpec> const & tree, const FibreWaveletTreeStructure)
{
    return tree.waveletTreeStructure;
}

template <typename TText, typename TSpec>
inline void clear(WaveletTree<TText, TSpec> & tree)
{
    for (unsigned i = 0; i < length(tree.bitStrings); ++i)
    {
        clear(tree.bitStrings[i]);
    }
    resize(tree.bitStrings, 0);
    clear(tree.waveletTreeStructure);
}

template <typename TText, typename TSpec>
inline unsigned getNumNodes(WaveletTree<TText, TSpec> & tree)
{
    return length(tree.waveletTreeStructure.treeNodes);
}

template <
    typename TText,
    typename TWaveletTreeSpec,
    typename TCharIn,
    typename TPos>
inline unsigned getOccImpl(const WaveletTree<TText, TWaveletTreeSpec> & tree,
                           const TCharIn character,
                           const TPos pos)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeNodes>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
    typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;

    TPos sum = pos + 1;
    TPointer treePos = 0;


    typename Iterator<const TWaveletTreeStructure>::Type it(tree.waveletTreeStructure, treePos);
    TChar charInTree = tree.waveletTreeStructure.minCharValue;
    do
    {
        TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
        if (character < getFibre(tree, FibreWaveletTreeStructure()).treeNodes[treePos].i1)
        {
            sum -= addValue;
            if (!goLeftChild(it))
                break;
        }
        else
        {
            charInTree = getCharacter(it);
            sum = addValue;
            if (!goRightChild(it))
                break;
        }
        treePos = getPosition(it);
    }
    while (sum);

    if (character == charInTree)
        return sum;

    return 0;
}

template <typename TText, typename TChar, typename TPos, typename TWaveletTreeSpec>
inline unsigned getOcc(const WaveletTree<TText, TWaveletTreeSpec> & tree,
                       const TChar character,
                       const TPos pos)
{
    return getOccImpl(tree, character, pos);
}

template <typename TText, typename TChar, typename TPos>
inline unsigned getOcc(const WaveletTree<TText, FmiDollarSubstituted<> > & tree,
                       const TChar character,
                       const TPos pos)
{
    unsigned occ = getOccImpl(tree, character, pos);
    if (character == tree.dollarSubstitute && pos >= tree.dollarPosition)
    {
        return occ - 1;
    }
    return occ;
}

template <typename TText, typename TWaveletTreeSpec, typename TPos>
inline typename Value<TText>::Type
getCharacter(const WaveletTree<TText, TWaveletTreeSpec> & tree,
             const TPos pos)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreWaveletTreeStructure>::Type TSplitValues;
    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
    typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;

    TPos sum = pos + 1;
    TPointer treePos = 0;

    typename Iterator<const TSplitValues>::Type iter(tree.waveletTreeStructure, treePos);
    bool direction;
    TChar character = tree.waveletTreeStructure.minCharValue;
    do
    {
        direction = getBit(tree.bitStrings[treePos], sum - 1);
        TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
        if (direction)
        {
            character = getCharacter(iter); // + 1;
            sum = addValue;
            if (!goRightChild(iter))
                break;
        }
        else
        {
            sum -= addValue;
            if (!goLeftChild(iter))
                break;
        }
        treePos = getPosition(iter);
    }
    while (true);

    return character;
}

template <typename TText, typename TWaveletTreeSpec>
inline typename Value<TText>::Type
getDollarSub(const WaveletTree<TText, TWaveletTreeSpec> & tree)
{
    return tree.dollarSub;
}

template <typename TText, typename TWaveletTreeSpec>
inline TText getAlphabet(WaveletTree<TText, TWaveletTreeSpec> & tree)
{
    getAlphabet(tree.splitValues);
}

inline unsigned getTreeLevel(unsigned treePosition)
{
    return floor(log(treePosition + 1) / log(2));
}

template <typename TText, typename TChar>
inline void setDollarSubstitute(WaveletTree<TText, FmiDollarSubstituted<> > & tree,
                                TChar dollarSubstitute)
{
    tree.dollarSubstitute = dollarSubstitute;
}

template <typename TText, typename TPos>
inline void setDollarPosition(WaveletTree<TText, FmiDollarSubstituted<> > & tree,
                              TPos position)
{
    tree.dollarPosition = position;
}

template <typename TText, typename TWaveletTreeSpec>
inline void fillWaveletTree(WaveletTree<TText, TWaveletTreeSpec> & tree,
                            const TText & text)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreBitStrings>::Type             TFibreRankSupportBitStrings;
    typedef typename Value<TFibreRankSupportBitStrings>::Type                                       TFibreRankSupportBitString;
    typedef typename Fibre<TFibreRankSupportBitString, FibreBitString>::Type                        TFibreBitString;
    typedef typename Size<TFibreBitString>::Type                                                    TSize;
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreWaveletTreeStructure>::Type   TSplitValues;

    resize(tree.bitStrings, length(tree.waveletTreeStructure));

    for (TSize i = 0; i < length(text); ++i)
    {
        typename Iterator<TSplitValues>::Type it(tree.waveletTreeStructure, 0);
        bool bit;

        do
        {
            if (value(text, i) < getCharacter(it))
            {
                bit = 0;
                append(getFibre(tree, FibreBitStrings())[getPosition(it)], bit);
                if (!goLeftChild(it))
                    break;
            }
            else
            {
                bit = 1;
                append(getFibre(tree, FibreBitStrings())[getPosition(it)], bit);
                if (!goRightChild(it))
                    break;
            }
        }
        while (true);
    }

    TFibreRankSupportBitStrings & bitStrings = getFibre(tree, FibreBitStrings());
    for (TSize i = 0; i < length(bitStrings); ++i)
    {
        TFibreRankSupportBitString & temp = bitStrings[i];
        completeRankSupportBitString(temp);
    }
}

template <typename TText, typename TWaveletTreeSpec, typename TChar>
inline unsigned getNodePosition(WaveletTree<TText, TWaveletTreeSpec> & tree, TChar character)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreWaveletTreeStructure>::Type TSplitValues;
    typename Iterator<TSplitValues>::Type iter(tree.splitValues, 0);
    return getNodePosition(iter, character);
}

template <typename TText, typename TWaveletTreeSpec>
inline void createWaveletTree(WaveletTree<TText, TWaveletTreeSpec> & waveletTree,
                              TText const & text)
{
    computeWaveletTreeStructure(getFibre(waveletTree, FibreWaveletTreeStructure()), text);
    fillWaveletTree(waveletTree, text);
}

template <typename TText, typename TWaveletTreeSpec, typename TPrefixSumTable, typename TDollarSub, typename TDollarPos>
inline void createWaveletTree(LFTable<WaveletTree<TText, TWaveletTreeSpec>, TPrefixSumTable> & lfTable,
                              TText const & text,
                              TDollarSub const dollarSub,
                              TDollarPos const dollarPos)
{
    setDollarSubstitute(lfTable.occTable, dollarSub);
    setDollarPosition(lfTable.occTable, dollarPos);

    computeWaveletTreeStructure(lfTable);
    fillWaveletTree(getFibre(lfTable, FibreOccTable()), text);
}

}
#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_
