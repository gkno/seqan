// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Code for read/write access to BAM tag dicts.
// ==========================================================================

// TODO(holtgrew): Test me!
// TODO(holtgrew): assignValue() is missing.

#ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
#define EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.BamTagsDict
..cat:Fragment Store
..signature:BamTagsDict
..summary:Indexes start positions of BAM tags in a @Class.CharString@ and provides a dict-like API.
..example.code:
CharString str = "AA:value1\tAB:value2";
BamTagsDict tags(str);
std::cerr << length(tags) << std::endl;  // #=> "2"
for (unsigned i = 0; i < length(tags); ++i)
    std::cerr << getKey(tags, i) << " -> " << getValue(tags, i) << std::endl;
// #=> "AA -> value1"
// #=> "AB -> value2"
..include:seqan/bam_io.h
*/

class BamTagsDict
{
public:
    Holder<CharString> _host;
    String<unsigned> _positions;

    BamTagsDict() {}

    BamTagsDict(CharString & tags) : _host(tags) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <>
struct Host<BamTagsDict>
{
    typedef CharString Type;
};

template <>
struct Host<BamTagsDict const>
{
    typedef CharString const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

inline Host<BamTagsDict>::Type &
host(BamTagsDict & bamTags)
{
    return value(bamTags._host);
}

inline Host<BamTagsDict const>::Type &
host(BamTagsDict const & bamTags)
{
    return value(bamTags._host);
}

// ----------------------------------------------------------------------------
// Function hasIndex()
// ----------------------------------------------------------------------------

/**
.Function.hasIndex
..cat:Fragment Store
..summary:Return $true$ if @Class.BamTagsDict@ has an index.
..signature:hasIndex(bamTags)
..param.bamTags:SAM Tags to query
...type:Class.BamTagsDict
..returns:$bool$
..include:<seqan/store_ex.h>
*/

inline bool
hasIndex(BamTagsDict const & bamTags)
{
    return length(bamTags._positions) != 0u;
}

inline bool
hasIndex(BamTagsDict & bamTags)
{
    return hasIndex(const_cast<BamTagsDict const &>(bamTags));
}

// ----------------------------------------------------------------------------
// Function getTypeSize()
// ----------------------------------------------------------------------------

// Return sizeof() of the type identified with the given char.  Returns -1 if not
// valid, -2 if of variable length.

inline int
getTypeSize(char c)
{
    switch (c)
    {
        case 'A':
            return 1;
        case 'f':
            return 4;
        case 'Z':
        case 'H':
        case 'B':
            return -1;
        case 'c':
        case 'C':
            return 1;
        case 's':
        case 'S':
            return 2;
        case 'i':
        case 'I':
            return 4;
    }
    return -2;
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

/**
.Function.buildIndex
..cat:Fragment Store
..summary:Build index for a @Class.BamTagsDict@ object.
..signature:buildIndex(bamTags)
..param.bamTags:SAM Tags to build index for.
...type:Class.BamTagsDict
..returns:$void$
..include:<seqan/store_ex.h>
*/

inline void
buildIndex(BamTagsDict & bamTags)
{
    typedef Host<BamTagsDict>::Type TCharString;
    typedef Iterator<TCharString, Rooted>::Type TCharStringIter;

    clear(bamTags._positions);
    appendValue(bamTags._positions, 0);
    for (TCharStringIter it = begin(host(bamTags)); !atEnd(it);)
    {
        it += 2;
        char c = *it;
        if (c == 'H' || c == 'Z')
        {
            while (!atEnd(it) && *it != '\0')
                ++it;
            ++it;
        }
        else if (c == 'B')
        {
            ++it;
            c = *it;
            ++it;
            __uint32 len = 0;
            memcpy(&len, &*it, 4);
            it += 4;
            it += len * getTypeSize(c);
        }
        else
        {
            ++it;
            it += getTypeSize(c);
        }

        appendValue(bamTags._positions, position(it));
    }
    if (length(host(bamTags)) != 0u)
        appendValue(bamTags._positions, length(host(bamTags)) + 1);  // +1 since there is not tab at the end
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

inline void setHost(BamTagsDict & bamTags, CharString & newHost)
{
    set(bamTags._host, newHost);
    buildIndex(bamTags);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

inline unsigned
length(BamTagsDict const & tags)
{
    if (!hasIndex(tags))
        buildIndex(const_cast<BamTagsDict &>(tags));
    return length(tags._positions) - 1;
}

// ----------------------------------------------------------------------------
// Function getType()
// ----------------------------------------------------------------------------

template <typename TPos>
inline char
getType(BamTagsDict & tags, TPos idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);
    return host(tags)[tags._positions[idx] + 2];
}

// ----------------------------------------------------------------------------
// Function getKey()
// ----------------------------------------------------------------------------

template <typename TPos>
inline Infix<Host<BamTagsDict>::Type>::Type
getKey(BamTagsDict & tags, TPos idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);
    return infix(host(tags), tags._positions[idx], tags._positions[idx] + 2);
}

template <typename TPos>
inline Infix<Host<BamTagsDict const>::Type>::Type
getKey(BamTagsDict const & tags, TPos idx)
{
    return getKey(const_cast<BamTagsDict &>(tags), idx);
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

// Note that you will get <type char> + payload.

template <typename TIdx>
inline CharString
getValue(BamTagsDict & tags, TIdx idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);
        
    typedef typename Position<CharString>::Type TPos;
    TPos beginPos = tags._positions[idx] + 2;
    TPos endPos = beginPos + 1;
    
    char theType = getType(tags, idx);
    if (theType == 'Z' || theType == 'H')
    {
        typedef typename Iterator<CharString, Rooted>::Type TIterator;
        TIterator it = begin(host(tags), Rooted()) + beginPos + 1;
        for (; !atEnd(it) && *it != '\0'; goNext(it))
            endPos += 1;
        endPos += 1;
    }
    else if (theType == 'B')
    {
        __uint32 len = 0;
        memcpy(&len, &host(tags)[tags._positions[idx]] + 4, 4);
        char c = host(tags)[tags._positions[idx] + 3];
        int typeSize = getTypeSize(c);
        SEQAN_ASSERT_GT(typeSize, 0);
        endPos += 5 + len * typeSize;
    }
    else
    {
        endPos += getTypeSize(theType);
    }
    
    return infix(host(tags), beginPos, endPos);
}

template <typename TPos>
inline CharString //Infix<Host<BamTagsDict const>::Type>::Type
getValue(BamTagsDict const & tags, TPos idx)
{
    return getValue(const_cast<BamTagsDict &>(tags), idx);
}

// ----------------------------------------------------------------------------
// Function extractValue()
// ----------------------------------------------------------------------------

template <typename TDest, typename TIdx>
inline bool
extractValue(TDest & dest, BamTagsDict & tags, TIdx idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    char * ptr = reinterpret_cast<char *>(&dest);
    int typeSize = getTypeSize(host(tags)[tags._positions[idx] + 2]);
    if (typeSize < 0)
        return false;
    memcpy(ptr, &host(tags)[tags._positions[idx] + 3], typeSize);
    return true;
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
