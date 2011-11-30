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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREESTRUCTURE_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREESTRUCTURE_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

template <typename TValue, typename TSpec = void>
struct WaveletTreeStructure;
struct FibreTreeNodes_;
//struct FibreTreeNodesPlusDollar_;

// ==========================================================================
//Tags, Classes, Enums
// ==========================================================================
typedef Tag<FibreTreeNodes_> const FibreTreeNodes;
//typedef Tag<FibreTreeNodesPlusDollar_> const FibreTreeNodesPlusDollar;

template <typename TValue, typename TSpec>
struct Fibre<WaveletTreeStructure<TValue, TSpec>, FibreTreeNodes>
{
	//TODO
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
	//typedef unsigned TValue;
    typedef String<Pair<TValue, TValue> > Type;
    typedef TValue Value;
};

//template <typename TText>
//struct Fibre<WaveletTreeStructure<TValue, TSpec>, FibreTreeNodesPlusDollar>
//{
//    typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
//	//typedef unsigned TValue;
//    typedef String<Pair<TValue, TValue> > Type;
//};

template <typename TValue, typename TSpec>
struct WaveletTreeStructure
{
    typename Fibre<WaveletTreeStructure, FibreTreeNodes>::Type treeNodes;

    WaveletTreeStructure() :
    	treeNodes()
    {};

    template <typename TFreqString>
    WaveletTreeStructure(TFreqString & freqString) :
    	treeNodes()
    {
        //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
        TValue sigmaSize = length(freqString);

        resize(treeNodes, sigmaSize - 1);
    	TValue numChildNodes, lowest, highest, posInTree;
        initStartTreeStructureValue(numChildNodes, lowest, highest, posInTree, sigmaSize);

        computeTreeEntries(freqString,
                           *this,
                           lowest,
                           highest,
                           posInTree,
                           numChildNodes);
    }

    bool operator==(const WaveletTreeStructure & b) const
    {
        return treeNodes == b.treeNodes;
    }
};

template <typename TValue>
void initStartTreeStructureValue(TValue & numChildNodes,
								 TValue & lowest,
								 TValue & highest,
								 TValue & posInTree,
								 TValue const sigmaSize)
{
	numChildNodes = 0;
	lowest = 0;
	highest = sigmaSize - 1;
	posInTree = 0;
}


// ==========================================================================
// Functions
// ==========================================================================

template <typename TValue, typename TSpec>
typename Fibre<WaveletTreeStructure<TValue, TSpec>, FibreTreeNodes>::Type &
getFibre(WaveletTreeStructure<TValue, TSpec> & treeStructure, FibreTreeNodes)
{
    return treeStructure.treeNodes;
}

template <typename TValue, typename TSpec>
typename Fibre<WaveletTreeStructure<TValue, TSpec>, FibreTreeNodes>::Type const &
getFibre(WaveletTreeStructure<TValue, TSpec> const & treeStructure, FibreTreeNodes)
{
    return treeStructure.treeNodes;
}

/*template< typename TBitString, typename TText, typename TSpec >
typename Fibre< WaveletTree< TBitString, TSplitValue, TPosInSubTree, TSpec >, FibreSplitValues >::Type const &
getFibre(WaveletTree< TBitString, TSplitValue, TPosInSubTree, TSpec > const &tree, const FibreSplitValues)
{
    return tree.splitValues;
}*/

template <typename TValue, typename TSpec>
unsigned length(WaveletTreeStructure<TValue, TSpec> & tree)
{
    return length(tree.treeNodes);
}

template <typename TValue, typename TSpec, typename TSize>
void resize(WaveletTreeStructure<TValue, TSpec> & treeStructure, TSize size)
{
    resize(treeStructure.treeNodes, size);
}

template <typename TValue, typename TSpec, typename TSize>
void resize(WaveletTreeStructure<TValue, TSpec> & treeStructure, TSize size,
			typename Value<typename Fibre<WaveletTreeStructure<TValue, TSpec>,FibreTreeNodes>::Type>::Type value)
{
    resize(treeStructure.treeNodes, size, value);
}

template <typename TValue, typename TSpec>
void clear(WaveletTreeStructure<TValue, TSpec> & treeStructure)
{
    clear(treeStructure.treeNodes);
}

template <typename TValue, typename TSpec>
struct Iter<WaveletTreeStructure<TValue, TSpec> const, Standard>
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    TValue position;
    const WaveletTreeStructure<TValue, TSpec> * waveletTreeStructure;

    template <typename TPos>
    Iter(const WaveletTreeStructure<TValue, TSpec> & treeStructure, const TPos pos) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};



template <typename TValue, typename TSpec>
struct Iter<WaveletTreeStructure<TValue, TSpec>, Standard>
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    TValue position;
    WaveletTreeStructure<TValue, TSpec> * waveletTreeStructure;

    template <typename TPos>
    Iter(WaveletTreeStructure<TValue, TSpec> & treeStructure, TPos pos) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};


//error: no matching function for call to Ô
//
//Iter<const WaveletTreeStructure<unsigned char, void>, const seqan::Tag<seqan::Standard_> >::
//
//Iter(const seqan::WaveletTreeStructure<unsigned char, void>&, seqan::getOccImpl(const seqan::WaveletTree<TText, TSpec>&, TTreeSplitValue, TPos)
//[with TText = seqan::String<char, seqan::Alloc<void> >,
// TWaveletTreeSpec = seqan::Tag<seqan::DollarSubstituted_>, TTreeSplitValue = seqan::String<long unsigned int, seqan::Alloc<void> >, TPos = unsigned int]::TValue&)Õ
///home/fenn/jsinger/Master/Development/seqan-trunk/sandbox/singer/include/seqan/fm_sequence/wavelet_tree_structure.h:177:
//note: candidates are: seqan::Iter<const seqan::WaveletTreeStructure<TValue, TSpec>, const seqan::Tag<seqan::Standard_> >::Iter(const seqan::WaveletTreeStructure<TValue, TSpec>&, TValue) [with TValue = unsigned char, TSpec = void]

template <typename TValue, typename TSpec>
struct Iterator<WaveletTreeStructure<TValue, TSpec>, Standard>
{
    typedef Iter<WaveletTreeStructure<TValue, TSpec>, Standard> Type;
};

template <typename TValue, typename TSpec>
struct Iterator<WaveletTreeStructure<TValue, TSpec> const, Standard>
{
    typedef Iter<WaveletTreeStructure<TValue, TSpec> const, Standard> Type;
};

template <typename TValue, typename TSpec>
struct Iterator<WaveletTreeStructure<TValue, TSpec>, Rooted>:
    Iterator<WaveletTreeStructure<TValue, TSpec>, Standard>{};

template <typename TValue, typename TSpec>
struct Iterator<WaveletTreeStructure<TValue, TSpec> const, Rooted>:
    Iterator<WaveletTreeStructure<TValue, TSpec> const, Standard>{};


template <typename TValue, typename TSpec>
struct Value<WaveletTreeStructure<TValue, TSpec> >
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    typedef Pair<TValue, TValue> Type;
};

template <typename TValue, typename TSpec>
struct Value<WaveletTreeStructure<TValue, TSpec> const>
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    typedef Pair<TValue, TValue> const Type;
};

template <typename TValue, typename TSpec>
struct Reference<WaveletTreeStructure<TValue, TSpec> >
{
    typedef typename Value<WaveletTreeStructure<TValue, TSpec> >::Type Type;
};

template <typename TValue, typename TSpec>
struct Reference<const WaveletTreeStructure<TValue, TSpec> >
{
    typedef typename Value<WaveletTreeStructure<TValue, TSpec> >::Type const Type;
};


template <typename TValue, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TValue, TSpec> const>::Type
begin(WaveletTreeStructure<TValue, TSpec> const & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type(waveletTreeStructure, 0);
}

template <typename TValue, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TValue, TSpec> const>::Type
end(WaveletTreeStructure<TValue, TSpec> const & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type(waveletTreeStructure, length(waveletTreeStructure.treeNodes));
}

template <typename TValue, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type
begin(WaveletTreeStructure<TValue, TSpec> & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type(waveletTreeStructure.treeNodes, 0);
}

template <typename TValue, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type
end(WaveletTreeStructure<TValue, TSpec> & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type(waveletTreeStructure.treeNodes, length(waveletTreeStructure.treeNodes));
}

template <typename TValue, typename TSpec, typename TPos>
void setPosition(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter, TPos pos)
{
    iter.position = pos;
}

template <typename TValue, typename TSpec>
TValue getPosition(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    return iter.position;
}

template <typename TValue, typename TSpec>
TValue getPosition(Iter<const WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    return iter.position;
}

template <typename TValue, typename TSpec>
void setCharacter(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter,
                  TValue character)
{
    iter.waveletTreeStructure->treeNodes[iter.position].i1 = character;
}

template <typename TValue, typename TSpec>
TValue getCharacter(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    return iter.waveletTreeStructure->treeNodes[iter.position].i1;
}

template <typename TValue, typename TSpec>
TValue getCharacter(Iter<const WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    return iter.waveletTreeStructure->treeNodes[iter.position].i1;
}

template <typename TValue, typename TSpec>
void setLeftChildPos(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 0)
    {
        iter.waveletTreeStructure->treeNodes[iter.position].i2 = 2;
        return;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
    {
        return;
    }
    std::cerr << "ERROR: The right child has just been deleted!" << std::endl;
}

template <typename TValue, typename TSpec>
TValue getLeftChildPos(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 > 1)
    {
        return iter.position + 1;
    }
    return 0;
}

template <typename TValue, typename TSpec>
TValue getLeftChildPos(Iter<const WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 > 1)
    {
        return iter.position + 1;
    }
    return 0;
}

template <typename TValue, typename TSpec>
void setRightChildPosOnly(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    iter.waveletTreeStructure->treeNodes[iter.position].i2 = 1;
}

template <typename TValue, typename TSpec, typename TPos>
void setRightChildPos(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter, TPos rightChildPosition)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 0)
    {
        iter.waveletTreeStructure->treeNodes[iter.position].i2 = 1;
        return;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
    {
        iter.waveletTreeStructure->treeNodes[iter.position].i2 = rightChildPosition + 2;
        return ;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
    {
        return;
    }
    iter.waveletTreeStructure->treeNodes[iter.position].i2 = rightChildPosition + 2;
}

template <typename TValue, typename TSpec>
TValue getRightChildPos(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 > 2)
    {
        return iter.waveletTreeStructure->treeNodes[iter.position].i2 - 2;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
    {
        return iter.position + 1;
    }
    return 0;
}

template <typename TValue, typename TSpec>
TValue getRightChildPos(Iter<const WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 > 2)
    {
        return iter.waveletTreeStructure->treeNodes[iter.position].i2 - 2;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
    {
        return iter.position + 1;
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename TValue, typename TSpec>
void setNodeToLeaf(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 != 0)
    {
        std::cerr << "You just deleted ";
        if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
        {
            std::cerr << "the right sub tree!" << std::endl;
        }
        if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
        {
            std::cerr << "the left sub tree!" << std::endl;
        }
        else
        {
            std::cerr << "both sub trees!" << std::endl;
        }
    }
    iter.waveletTreeStructure->treeNodes[iter.position].i2 = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename TValue, typename TSpec>
void goLeft(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    iter.position = getLeftChildPos(iter);
}

template <typename TValue, typename TSpec>
void goLeft(Iter<const WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    iter.position = getLeftChildPos(iter);
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename TValue, typename TSpec>
void goRight(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    iter.position = getRightChildPos(iter);
}

template <typename TValue, typename TSpec>
void goRight(Iter<const WaveletTreeStructure<TValue, TSpec>, Standard> & iter)
{
    iter.position = getRightChildPos(iter);
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename TValue, typename TSpec, typename TPosInSubTree>
void printTreeLevel(WaveletTreeStructure<TValue, TSpec> & tree, TPosInSubTree posInTree)
{
    TPosInSubTree leftChild = tree.treeNodes[posInTree].i2;
    TPosInSubTree rightChild = tree.treeNodes[posInTree].i3;
    std::cout << "(" << posInTree << ": " << tree.treeNodes[posInTree].i1 << ", " << leftChild << ", " << rightChild << ") | ";
    if (leftChild)
    {
        printTreeLevel(tree, leftChild);
    }
    if (rightChild)
    {
        printTreeLevel(tree, rightChild);
    }
}

template <typename TValue, typename TSpec>
void printTree(WaveletTreeStructure<TValue, TSpec> & tree)
{
    printTreeLevel(tree, 0);
}

template <typename TValue, typename TSpec, typename TString>
void writeGraphImpl(Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter, TString name)
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type iter2 = iter;
    std::ofstream stream(toCString(name), std::ios::app);
    TValue pos = getLeftChildPos(iter);
    if (pos)
    {
        stream << ordValue(iter.waveletTreeStructure->treeNodes[iter.position].i1) << " -> " << ordValue(iter.waveletTreeStructure->treeNodes[pos].i1) << ";" << std::endl;
        goLeft(iter);
        writeGraphImpl(iter, name);
    }
    else
    {
        stream << ordValue(iter.waveletTreeStructure->treeNodes[iter.position].i1) << " -> " << "leave1" << ordValue(iter.position) << ";" << std::endl;
    }

    pos = getRightChildPos(iter2);
    if (pos)
    {
        stream << ordValue(iter2.waveletTreeStructure->treeNodes[iter2.position].i1) << " -> " << ordValue(iter2.waveletTreeStructure->treeNodes[pos].i1) << ";" << std::endl;
        goRight(iter2);
        writeGraphImpl(iter2, name);
    }
    else
    {
        stream << ordValue(iter2.waveletTreeStructure->treeNodes[iter2.position].i1) << " -> " << "leave2" << ordValue(iter2.position) << ";" << std::endl;
    }
    stream.close();
}

template <typename TValue, typename TSpec>
void writeGraph(WaveletTreeStructure<TValue, TSpec> & tree)
{

    typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type iter(tree, 0);

    String<char> name = "testfile.dot";
    std::ofstream stream(toCString(name), std::ios::out);
    stream << "digraph G {" << std::endl;
    stream.close();
    writeGraphImpl(iter, name);

    stream.open(toCString(name), std::ios::app);
    stream << "}" << std::endl;
    stream.close();
}

template <typename TValue, typename TSpec, typename TString>
void writeGraph(WaveletTreeStructure<TValue, TSpec> & tree, TString name)
{

    typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type iter(tree, 0);
    std::ofstream stream(toCString(name), std::ios::out);
    stream << "digraph G {" << std::endl;
    stream.close();
    writeGraphImpl(iter, name);

    stream.open(toCString(name), std::ios::app);
    stream << "}" << std::endl;
    stream.close();
}

template <typename TStringLengthString, typename TPrefixSumTable, typename TIter>//, typename TTreeSplitValue>
void computeSingleStringLengthFromTree(TStringLengthString & lengthString,
                                       TPrefixSumTable & prefixSumTable,
                                       TIter & iter,
                                       unsigned lowerBound,
                                       unsigned upperBound)
//										TreeSplitValue lowerBound,
//                                       TTreeSplitValue upperBound)
{

    lengthString[iter.position] = prefixSumTable[upperBound + 1] - prefixSumTable[lowerBound];
    TIter iter2 = iter;

    std::cerr << "iter.position: " << (int)iter.position << " upperBound: " << ordValue(upperBound) << " lowerBound: " << ordValue(lowerBound) << " " <<prefixSumTable[upperBound + 1] << " " << prefixSumTable[lowerBound] << std::endl;

    //TTreeSplitValue newSplit = getCharacter(iter);
    unsigned newSplit = getCharacter(iter);
    goLeft(iter);
    if (getPosition(iter))
    {
        computeSingleStringLengthFromTree(lengthString,
                                          prefixSumTable,
                                          iter,
                                          lowerBound,
                                          newSplit);
    }

    newSplit = getCharacter(iter2) + 1;
    goRight(iter2);
    if (getPosition(iter2))
    {
        computeSingleStringLengthFromTree(lengthString,
                                          prefixSumTable,
                                          iter2,
                                          newSplit,
                                          upperBound);
    }
}

template <typename TStringLengthString, typename TPrefixSumTable, typename TValue, typename TSpec>
void computeStringLengthFromTree(TStringLengthString & lengthString,
                                 TPrefixSumTable & prefixSumTable,
                                 WaveletTreeStructure<TValue, TSpec> & tree)
{
    resize(lengthString, length(tree.treeNodes), 0);
    typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type iter(tree, 0);
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    //typedef unsigned TValue;

    computeSingleStringLengthFromTree(lengthString,
                                      prefixSumTable,
                                      iter,
                                      (TValue)0,
                                      (TValue)(length(prefixSumTable) - 2));
}

/*
    template< typename TFreqString,
        typename TText,
        typename TTreeNodeWidthString,
        typename TValue,
        typename TPosInTree,
        typename TNumOfChildNodes >
    void computeTreeEntries_(
        TFreqString &freq,
        String< Pair<  > > &treeNodes,
        TTreeNodeWidthString &treeNodeWidth,
        TValue smallestValue,
        TValue biggestValue,
        TPosInTree posInTree,
        TNumOfChildNodes &numChildNodes)
    {
        TValue oldSmallestValue = smallestValue;
        TValue oldBiggestValue = biggestValue;
        typedef typename Value< TFreqString >::Type TSize;
        TSize leftSize = freq[smallestValue];
        TSize rightSize = freq[biggestValue];

        if(smallestValue == biggestValue)
        {
            return;
        }

        if(smallestValue == biggestValue - 1)
        {
            treeNodes[posInTree].i1 = smallestValue;
            treeNodes[posInTree].i2 = 0;
            treeNodes[posInTree].i3 = 0;
            treeNodeWidth[posInTree] = freq[smallestValue] + freq[biggestValue];
            ++numChildNodes;
            return;
        }

        TValue leftCounter = 0;
        TValue rightCounter = 0;
        while(biggestValue - 1 > smallestValue)
        {
            if(leftSize < rightSize)
            {
                leftSize += freq[smallestValue + 1];
                ++smallestValue;
                ++leftCounter;
            }
            else if(leftSize > rightSize)
            {
                rightSize += freq[biggestValue - 1];
                --biggestValue;
                ++rightCounter;
            }
            else
            {
                if(leftCounter < rightCounter)
                {
                    leftSize += freq[smallestValue + 1];
                    ++smallestValue;
                    ++leftCounter;
                }
                else
                {
                    rightSize += freq[biggestValue - 1];
                    --biggestValue;
                    ++rightCounter;
                }
            }
        }

        treeNodeWidth[posInTree] = leftSize + rightSize;
        treeNodes[posInTree].i1 = smallestValue;

        computeTreeEntries_(freq, treeNodes, treeNodeWidth, oldSmallestValue, smallestValue, posInTree + 1, numChildNodes);
    //	treeNodes[posInTree].i2 = posInTree + 1;
        treeNodes[posInTree].i2 = posInTree + numChildNodes + 1;
        ++numChildNodes;

        TNumOfChildNodes numChildNodes2 = 0;
        computeTreeEntries_(freq, treeNodes, treeNodeWidth, biggestValue, oldBiggestValue, treeNodes[posInTree].i3, numChildNodes2);
        numChildNodes += numChildNodes2;
        return;
    }
*/

template <typename TFreqString,
          typename TValue,
          typename TSpec,
          typename TNumOfChildNodes>
void computeTreeEntries(
    TFreqString & freq,
    Iter<WaveletTreeStructure<TValue, TSpec>, Standard> & iter,
    TValue smallestValue,
    TValue biggestValue,
    TNumOfChildNodes & numChildNodes)
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TTreeValue;

    TValue oldSmallestValue = smallestValue;
    TValue oldBiggestValue = biggestValue;
    typedef typename Value<TFreqString>::Type TSize;
    TSize leftSize = freq[smallestValue];
    TSize rightSize = freq[biggestValue];


   // std::cerr << "pos: " << (int)getPosition(iter) << "smallestValue: " << (int) smallestValue << " biggestValue: " << (int)biggestValue << std::endl;

    if (smallestValue == biggestValue - 1)
    {
    	//std::cerr << "smallestValue: " << (int) smallestValue << " biggestValue: " << (int)biggestValue << std::endl;
        setCharacter(iter, (TValue)smallestValue);
    	//std::cerr << "computeTreeEntries01" << std::endl;
        setNodeToLeaf(iter);
    	//std::cerr << "computeTreeEntries02" << std::endl;
        ++numChildNodes;
        return;
    }
    //std::cerr << "computeTreeEntries0" << std::endl;
    TValue leftCounter = 0;
    TValue rightCounter = 0;
    while (biggestValue - 1 > smallestValue)
    {
        if (leftSize < rightSize)
        {
            ++smallestValue;
            leftSize += freq[smallestValue];
            ++leftCounter;
        }
        else if (leftSize > rightSize)
        {
            --biggestValue;
            rightSize += freq[biggestValue];
            ++rightCounter;
        }
        else
        {
            if (leftCounter < rightCounter)
            {
                ++smallestValue;
                leftSize += freq[smallestValue];
                ++leftCounter;
            }
            else
            {
                --biggestValue;
                rightSize += freq[biggestValue];
                ++rightCounter;
            }
        }
    }

    //std::cerr << "computeTreeEntries1" << std::endl;
    setCharacter(iter, (TValue)smallestValue);
    //if the elements of the node do not appear in the text make a leaf here
    if (!leftSize && !rightSize)
    {
        setNodeToLeaf(iter);
        ++numChildNodes;
        return;
    }

    //std::cerr << "computeTreeEntries2" << std::endl;
    //there must be either a left or right subtree or both
    typename Iterator<WaveletTreeStructure<TValue, TSpec> >::Type iter2 = iter;

    //is there only one element on the left side?
    if (oldSmallestValue == smallestValue)
    {
//        setRightChildPos(iter, iter.position + numChildNodes + 1);
//        ++numChildNodes;
    		setRightChildPos(iter2, iter2.position + numChildNodes + 1);
    	        TNumOfChildNodes numChildNodes2 = 0;
    	        goRight(iter2);
    	        computeTreeEntries(freq, iter2, biggestValue, oldBiggestValue, numChildNodes2);
    	        numChildNodes += numChildNodes2 + 1;
    	        return;
    }
    else
    {
        setLeftChildPos(iter);
        goLeft(iter);
        computeTreeEntries(freq, iter, oldSmallestValue, smallestValue, numChildNodes);
    }

    //std::cerr << "computeTreeEntries3" << std::endl;
    //is there only one element on the right side?
    if (biggestValue == oldBiggestValue)
    {
        setLeftChildPos(iter2);
        ++numChildNodes;
    }
    else
    {
        setRightChildPos(iter2, iter2.position + numChildNodes + 1);
        TNumOfChildNodes numChildNodes2 = 0;
        goRight(iter2);
        computeTreeEntries(freq, iter2, biggestValue, oldBiggestValue, numChildNodes2);
        numChildNodes += numChildNodes2 + 1;
        return;
    }
    //std::cerr << "computeTreeEntries4" << std::endl;
}

template <typename TValue, typename TSpec>
inline bool open(
    WaveletTreeStructure<TValue, TSpec> & structure,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".treestruct");    open(getFibre(structure, FibreTreeNodes()), toCString(name), openMode);
    return true;
}

template <typename TValue, typename TSpec>
inline bool open(WaveletTreeStructure<TValue, TSpec> & structure,
				 const char * fileName)
{
    return open(structure, fileName, DefaultOpenMode<WaveletTreeStructure<TValue, TSpec> >::VALUE);
}

template <typename TValue, typename TSpec>
inline bool save(
    WaveletTreeStructure<TValue, TSpec> const & structure,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".treestruct");    save(getFibre(structure, FibreTreeNodes()), toCString(name), openMode);
    return true;
}

template <typename TValue, typename TSpec>
inline bool save(
    WaveletTreeStructure<TValue, TSpec> const & structure,
    const char * fileName)
{
    return save(structure, fileName, DefaultOpenMode<WaveletTreeStructure<TValue, TSpec> >::VALUE);
}

}

#endif
