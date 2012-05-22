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
//       from this software withoFIut specific prior written permission.
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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELET_TREE_STRUCTURE_ITERATOR_BETA_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELET_TREE_STRUCTURE_ITERATOR_BETA_H_

namespace seqan {

// ==========================================================================
// Metafunctions
// ==========================================================================
template <typename TChar, typename TSpec, typename TIterSpec>
struct Container<Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<TIterSpec> > >
{
    typedef WaveletTreeStructure<TChar, TSpec> Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Container<Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<TIterSpec> > >
{
    typedef WaveletTreeStructure<TChar, TSpec> const Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Container<Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > >:
    Container<Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<> > >
{};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Container<Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<ParentLinks<TIterSpec> > > >:
    Container<Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<> > >
{};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<WaveletTreeStructure<TChar, TSpec>, TIterSpec>
{
    typedef Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<> > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<WaveletTreeStructure<TChar, TSpec> const, TIterSpec>
{
    typedef Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<> > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<WaveletTreeStructure<TChar, TSpec>, TopDown<TIterSpec> >
{
    typedef Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<> > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<WaveletTreeStructure<TChar, TSpec> const, TopDown<TIterSpec> >
{
    typedef Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<> > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > >
{
    typedef Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<> > > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<WaveletTreeStructure<TChar, TSpec> const, TopDown<ParentLinks<TIterSpec> > >
{
    typedef Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<ParentLinks<> > > Type;
};

// ==========================================================================
//Tags, Classes, Enums
// ==========================================================================

template <typename TChar, typename TSpec>
class Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<> >
{
    typedef typename Fibre<WaveletTreeStructure<TChar, TSpec>, FibreTreeNodes>::Type TWaveletTreeNodes;
    typedef typename Value<TWaveletTreeNodes>::Type TWaveletTreeNode;
    typedef typename Value<TWaveletTreeNode, 2>::Type TPos;

public:
    TPos position;
    WaveletTreeStructure<TChar, TSpec> const * waveletTreeStructure;

    template <typename TPos>
    Iter(WaveletTreeStructure<TChar, TSpec> const & treeStructure, const TPos pos) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};

template <typename TChar, typename TSpec>
class Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<> >
{
    typedef typename Fibre<WaveletTreeStructure<TChar, TSpec>, FibreTreeNodes>::Type TWaveletTreeNodes;
    typedef typename Value<TWaveletTreeNodes>::Type TWaveletTreeNode;
    typedef typename Value<TWaveletTreeNode, 2>::Type TPos;

public:
    TPos position;
    WaveletTreeStructure<TChar, TSpec> * waveletTreeStructure;

    template <typename TPos>
    Iter(WaveletTreeStructure<TChar, TSpec> & treeStructure, TPos pos) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};

template <typename TChar, typename TSpec, typename TIterSpec>
class Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > >
{
    typedef typename Fibre<WaveletTreeStructure<TChar, TSpec>, FibreTreeNodes>::Type TWaveletTreeNodes;
    typedef typename Value<TWaveletTreeNodes>::Type TWaveletTreeNode;
    typedef typename Value<TWaveletTreeNode, 2>::Type TPos;

public:
    String<TPos> position;
    WaveletTreeStructure<TChar, TSpec> * waveletTreeStructure;

    template <typename TPos>
    Iter(WaveletTreeStructure<TChar, TSpec> & treeStructure, TPos pos) :
        position(),
        waveletTreeStructure(&treeStructure)
    {
        appendValue(position, pos);
    }

};

template <typename TChar, typename TSpec, typename TIterSpec>
class Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<ParentLinks<TIterSpec> > >
{
    typedef typename Fibre<WaveletTreeStructure<TChar, TSpec>, FibreTreeNodes>::Type TWaveletTreeNodes;
    typedef typename Value<TWaveletTreeNodes>::Type TWaveletTreeNode;
    typedef typename Value<TWaveletTreeNode, 2>::Type TPos;

public:
    String<TPos> position;
    WaveletTreeStructure<TChar, TSpec> const * waveletTreeStructure;

    template <typename TPos>
    Iter(const WaveletTreeStructure<TChar, TSpec> & treeStructure, const TPos pos) :
        position(),
        waveletTreeStructure(&treeStructure)
    {
        appendValue(position, pos);
    }

};

template <typename TChar, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TChar, TSpec> const>::Type
begin(WaveletTreeStructure<TChar, TSpec> const & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TChar, TSpec> >::Type(waveletTreeStructure, 0);
}

template <typename TChar, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TChar, TSpec> >::Type
begin(WaveletTreeStructure<TChar, TSpec> & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TChar, TSpec> >::Type(waveletTreeStructure.treeNodes, 0);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline WaveletTreeStructure<TChar, TSpec> &
container(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & it)
{
    return *it.waveletTreeStructure;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline WaveletTreeStructure<TChar, TSpec> const &
container(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > const & it)
{
    return *it.waveletTreeStructure;
}

template <typename TChar, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TChar, TSpec> >::Type
end(WaveletTreeStructure<TChar, TSpec> & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TChar, TSpec> >::Type(waveletTreeStructure.treeNodes, length(waveletTreeStructure.treeNodes));
}

template <typename TChar, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TChar, TSpec> const>::Type
end(WaveletTreeStructure<TChar, TSpec> const & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TChar, TSpec> >::Type(waveletTreeStructure, length(waveletTreeStructure.treeNodes));
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline TChar getCharacter(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    return iter.waveletTreeStructure->treeNodes[getPosition(iter)].i1;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline TChar getCharacter(Iter<const WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    return iter.waveletTreeStructure->treeNodes[getPosition(iter)].i1;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getLeftChildPos(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 > 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getLeftChildPos(Iter<const WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 > 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

template <typename TText, typename TIter>
inline void getNexNode(TText & alphabet, TIter & it)
{
    TIter copyIt = it;
    if (getLeftChildPos(it))
    {
        goLeftChild(it);
        getNexNode(alphabet, it);
    }
    appendValue(alphabet, getCharacter(copyIt));
    if (getRightChildPos(it))
    {
        goRightChild(copyIt);
        getNexNode(alphabet, copyIt);
    }
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getNodePosition(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter, TChar character)
{
    while (!isLeaf(iter))
    {
        if (character < getCharacter(iter))
            goLeftChild(iter);
        else
            goRightChild(iter);
    }
    return getPosition(iter);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getNodePosition(Iter<WaveletTreeStructure<TChar, TSpec> const, TIterSpec> & iter, TChar character)
{
    while (!isLeaf(iter))
    {
        if (character < getCharacter(iter))
            goLeftChild(iter);
        else
            goRightChild(iter);
    }
    return getPosition(iter);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned getNumChildNodes(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & it)
{
    unsigned counter = 0;
    unsigned pos = getPosition(it);

    if (goDown(it))
    {
        do
        {
            ++counter;

            if (goDown(it))
            {
                continue;
            }
            else
            {
                if (goRight(it))
                {
                    continue;
                }
                else
                    while (goUp(it) && !goRight(it) && pos < getPosition(it))
                        ;
            }
        }
        while (pos < getPosition(it));
        return counter;
    }

    return 0;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getPosition(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    return iter.position;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getPosition(Iter<const WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    return iter.position;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getPosition(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & iter)
{
    return iter.position[length(iter.position) - 1];
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getPosition(Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<ParentLinks<TIterSpec> > > & iter)
{
    return iter.position[length(iter.position) - 1];
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getRightChildPos(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 > 2)
    {
        return iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 - 2;
    }
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getRightChildPos(Iter<const WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 > 2)
    {
        return iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 - 2;
    }
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

/**
.Function.WaveletTreeStructure#goDown
..signature:bool goDown(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.TopDown Iterator
..param.char:$iterator$ goes down the edge beginning with $char$.
...type:Class.WaveletTreeStructure
..remarks:$goDown(iterator)$ goes down the left edge if it exist, the right edge otherwise. 
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goDown(it); // go to left child of root node
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goDown(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    if (goLeftChild(iter))
        return true;

    if (goRightChild(iter))
        return true;

    return false;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goDown(Iter<WaveletTreeStructure<TChar, TSpec> const, TIterSpec> & iter)
{
    if (goLeftChild(iter))
        return true;

    if (goRightChild(iter))
        return true;

    return false;
}

/**
.Function.goLeftChild
..signature:bool goLeftChild(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.TopDown Iterator
...type:Class.WaveletTreeStructure
..remarks:$goLeftChild(iterator)$ goes down the left edge if it exist. 
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goLeftChild(it); // go to left child of root node
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goLeftChild(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    unsigned leftChildPos = getLeftChildPos(iter);
    if (leftChildPos == 0)
        return false;

    setPosition(iter, leftChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goLeftChild(Iter<WaveletTreeStructure<TChar, TSpec> const, TIterSpec> & iter)
{
    unsigned leftChildPos = getLeftChildPos(iter);
    if (leftChildPos == 0)
        return false;

    setPosition(iter, leftChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goLeftChild(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & iter)
{
    for (unsigned i = 0; i < length(iter.position); ++i)
        std::cerr << iter.position[i] << "\t";
    std::cerr << std::endl;

    unsigned leftChildPos = getLeftChildPos(iter);
    if (leftChildPos == 0)
        return false;

    appendValue(iter.position, leftChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goLeftChild(Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<ParentLinks<TIterSpec> > > & iter)
{
    unsigned leftChildPos = getLeftChildPos(iter);
    if (leftChildPos == 0)
        return false;

    appendValue(iter.position, leftChildPos);
    return true;
}

/**
.Function.goRight
..param.iterator:
...type:Class.WaveletTreeStructure
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goDown(it); // go to left child of root node
goRight(it); // go to right child of root node
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRight(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    unsigned pos = getPosition(iter);
    if (goUp(iter))
    {
        if (goRightChild(iter))
        {
            if (pos != getPosition(iter))
                return true;
        }
        else
        {
            goToPosition(iter, pos);
        }
    }

    return false;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRight(Iter<WaveletTreeStructure<TChar, TSpec> const, TIterSpec> & iter)
{
    unsigned pos = getPosition(iter);
    if (goUp(iter))
    {
        if (goRightChild(iter))
        {
            if (pos != getPosition(iter))
                return true;
        }
        else
        {
            goToPosition(iter, pos);
        }

    }

    return false;
}

/**
.Function.goRightChild
..signature:bool goRightChild(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.TopDown Iterator
...type:Class.WaveletTreeStructure
..remarks:$goRightChild(iterator)$ goes down the right edge if it exist. 
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRightChild(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    unsigned rightChildPos = getRightChildPos(iter);
    if (rightChildPos == 0)
        return false;

    setPosition(iter, rightChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRightChild(Iter<const WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    unsigned rightChildPos = getRightChildPos(iter);
    if (rightChildPos == 0)
        return false;

    setPosition(iter, rightChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRightChild(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & iter)
{
    unsigned rightChildPos = getRightChildPos(iter);
    if (rightChildPos == 0)
        return false;

    appendValue(iter.position, rightChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRightChild(Iter<const WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & iter)
{
    unsigned rightChildPos = getRightChildPos(iter);
    if (rightChildPos == 0)
        return false;

    appendValue(iter.position, rightChildPos);
    return true;
}


/**
.Function.goToPosition
..signature:bool goToPosition(iterator, pos)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.TopDown Iterator
..param.pos:A position.
..remarks:$goToPosition(iterator)$ goes to position pos if it exist. 
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goToPosition(it, 2); // go to right child of root node
*/

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline bool goToPosition(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter, TPos pos)
{
    if (pos >= length(container(iter).treeNodes))
        return false;

    iter.position = pos;
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline bool goToPosition(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & iter, TPos pos)
{
    if (pos >= length(container(iter).treeNodes))
        return false;

    appendValue(iter.position, pos);
    return true;
}

/**
.Function.WaveletTreeStructure#goUp
..signature:goUp(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.TopDownHistory Iterator
..remarks:$goUp(iterator)$ goes to the parent node.
..returns:$true$ if the current node is not the root node. 
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
goUp(it); // go to root node
*/

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goUp(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & it)
{
    unsigned treeLevel = length(it.position);

    if (isRoot(it))
        return false;

    resize(it.position, treeLevel - 1);
    setPosition(it, it.position[treeLevel - 2]);

    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goUp(Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<ParentLinks<TIterSpec> > > & it)
{
    unsigned treeLevel = length(it.position);

    if (isRoot(it))
        return false;

    resize(it.position, treeLevel - 1);
    setPosition(it, it.position[treeLevel - 2]);

    return true;
}

// This function implements the functionality of go up and 
// resizes the borderString of the structure construction.
template <typename TChar, typename TSpec, typename TIterSpec, typename TBorderString>
inline bool goUpStructureConstruction_(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & it, TBorderString & borderString)
{
    unsigned treeLevel = length(it.position);

    if (isRoot(it))
        return false;

    resize(borderString, treeLevel - 1);
    resize(it.position, treeLevel - 1);
    setPosition(it, it.position[treeLevel - 2]);

    return true;
}

/**
.Function.isLeaf
..param.iterator:An iterator of a wavelet tree structure.
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
goUp(it); // go to root node
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool isLeaf(Iter<WaveletTreeStructure<TChar, TSpec> const, TIterSpec> & iter)
{
    return !((*iter.waveletTreeStructure).treeNodes[getPosition(iter)].i2);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool isLeaf(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    return !((*iter.waveletTreeStructure).treeNodes[getPosition(iter)].i2);
}

// This function creates the right sibling of the current node
// and goes to that one.
template <typename TChar, typename TSpec, typename TIterSpec, typename TBorderString, typename TPrefixSumTable>
inline bool setAndGoRight(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & it, TBorderString & borderString, TPrefixSumTable & pst)
{
    if (isRoot(it) || (borderString[length(borderString) - 1].i2 == borderString[length(borderString) - 2].i2))
        return false;

    unsigned pos = getPosition(it);
    unsigned numChildNodes = getNumChildNodes(it);
    goUp(it);

    TChar pivot = getCharacter(it);
    setRightChildPos(it, pos + numChildNodes + 1);
    goRightChild(it);
    borderString[length(borderString) - 1].i1 = getCharacterPosition(pst, pivot);
    borderString[length(borderString) - 1].i2 = borderString[length(borderString) - 2].i2;

    return true;
}

/**
.Function.setCharacter
..signature:bool setCharacter(iterator, character)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.TopDown Iterator
..param.character:The character to be assigned to a node.
...type:TChar
..remarks:$setCharacter(iterator, character)$ sets the character of the node the iterator points to to character.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
setCharacter(it,'T'); // sets the character of the root's
                      // right child to 'T'
*/
template <typename TChar, typename TSpec, typename TIterSpec, typename TChar2>
inline void setCharacter(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter,
                         TChar2 character)
{
    iter.waveletTreeStructure->treeNodes[getPosition(iter)].i1 = character;
}

// This function sets the left child of the current node, or the right if there is no left child.
template <typename TChar, typename TSpec, typename TIterSpec, typename TBorderString, typename TCharPST, typename TSpecPST>
void setChildNodes_(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & it,
                   TBorderString & borderString,
                   PrefixSumTable<TCharPST, TSpecPST> & pst)
{
    typedef typename Value<TBorderString>::Type TBorderStringValue;
    unsigned leftBorder = borderString[length(borderString) - 1].i1;
    unsigned rightBorder = borderString[length(borderString) - 1].i2;
    unsigned pivotPosition = getPivotPosition(pst, leftBorder, rightBorder);

    setCharacter(it, getCharacter(pst, pivotPosition));

    if (leftBorder == pivotPosition - 1)
    {
        setRightChildPosOnly(it);
        appendValue(borderString, TBorderStringValue(pivotPosition, borderString[length(borderString) - 1].i2));
        return;
    }

    setLeftChildPos(it);
    appendValue(borderString, TBorderStringValue(borderString[length(borderString) - 1].i1, pivotPosition - 1));
}

// This functions sets the pointer to the left child.
template <typename TChar, typename TSpec, typename TIterSpec>
inline void setLeftChildPos(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 0)
    {
        iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 = 2;
        return;
    }
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 2)
    {
        return;
    }
    std::cerr << "ERROR: The right child has just been deleted!" << std::endl;
}

// This function sets the current node to be a node.
template <typename TChar, typename TSpec, typename TIterSpec>
inline void setNodeToLeaf(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 != 0)
    {
        std::cerr << "You just deleted ";
        if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 1)
        {
            std::cerr << "the right sub tree!" << std::endl;
        }
        if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 2)
        {
            std::cerr << "the left sub tree!" << std::endl;
        }
        else
        {
            std::cerr << "both sub trees!" << std::endl;
        }
    }
    iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 = 0;
}

// This function sets the position of iter to pos, 
template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline void setPosition(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter, TPos pos)
{
    iter.position = pos;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline void setPosition(Iter<WaveletTreeStructure<TChar, TSpec> const, TIterSpec> & iter, TPos pos)
{
    iter.position = pos;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline void setPosition(Iter<WaveletTreeStructure<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > > & iter, TPos pos)
{
    iter.position[length(iter.position) - 1] = pos;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline void setPosition(Iter<WaveletTreeStructure<TChar, TSpec> const, TopDown<ParentLinks<TIterSpec> > > & iter, TPos pos)
{
    iter.position[length(iter.position) - 1] = pos;
}

// This functions sets the pointer to the left child.
template <typename TChar, typename TSpec, typename TPos, typename TIterSpec>
inline void setRightChildPos(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter, TPos rightChildPosition)
{
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 0)
    {
        iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 = 1;
        return;
    }
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 2)
    {
        iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 = rightChildPosition + 2;
        return;
    }
    if (iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 == 1)
    {
        return;
    }
    iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 = rightChildPosition + 2;
}

// This function sets the pointer to the right child such that it is cleat that there is no left child.
template <typename TChar, typename TSpec, typename TIterSpec>
inline void setRightChildPosOnly(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & iter)
{
    iter.waveletTreeStructure->treeNodes[getPosition(iter)].i2 = 1;
}

/**
.Function.isRoot
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.TopDown Iterator
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTreeStructure<Dna5> waveletTreeStructure(genome);

Iterator<WaveletTreeStructure<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

isRoot(it) // returns true
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool isRoot(Iter<WaveletTreeStructure<TChar, TSpec>, TIterSpec> & it)
{
    return !getPosition(it);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool isRoot(Iter<WaveletTreeStructure<TChar, TSpec> const, TIterSpec> & it)
{
    return !getPosition(it);
}

}
#endif // SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELET_TREE_STRUCTURE_ITERATOR_BETA_H_
