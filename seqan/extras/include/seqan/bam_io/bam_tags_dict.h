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

// TODO(holtgrew): assignValue() is missing.
// TODO(holtgrew): eraseValue() is missing.

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
..summary:Indexes start positions of BAM tags in a $Shortcut.CharString@ and provides a dict-like API.
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
// Function getBamTypeSize()
// ----------------------------------------------------------------------------

// Return sizeof() of the type identified with the given char.  Returns -2 if not
// valid, -1 if of variable length.

/**
.Function.getBamTypeSize
..signature:getBamTypeSize(c)
..summary:Return size of the type identified by $c$.
..param.c:The BAM type identifier
..returns:$int$ with the $sizeof()$ of the type, -1 for variable sized types, -2 for invalid parameters.
..include:seqan/bam_io.h
*/

inline int
getBamTypeSize(char c)
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
            it += len * getBamTypeSize(c);
        }
        else
        {
            ++it;
            it += getBamTypeSize(c);
        }

        appendValue(bamTags._positions, position(it));
    }
    if (length(host(bamTags)) != 0u)
        appendValue(bamTags._positions, length(host(bamTags)) + 1);  // +1 since there is not tab at the end
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

///.Function.setHost.param.object.type:Class.BamTagsDict

inline void setHost(BamTagsDict & bamTags, CharString & newHost)
{
    set(bamTags._host, newHost);
    buildIndex(bamTags);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

///.Function.length.param.object.type:Class.BamTagsDict

inline unsigned
length(BamTagsDict const & tags)
{
    if (!hasIndex(tags))
        buildIndex(const_cast<BamTagsDict &>(tags));
    return length(tags._positions) - 1;
}

// ----------------------------------------------------------------------------
// Function getTagType()
// ----------------------------------------------------------------------------

/**
.Function.getTagType
..cat:BAM I/O
..signature:getTagType(tagsDict, idx)
..summary:Get key of a tag by index.
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
..param.idx:Index of the tag whose key to retrieve.
..returns:$char$, the SAM/BAM identifier of the type.
..include:seqan/bam_io.h
*/

template <typename TPos>
inline char
getTagType(BamTagsDict & tags, TPos idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);
    return host(tags)[tags._positions[idx] + 2];
}

// ----------------------------------------------------------------------------
// Function getTagKey()
// ----------------------------------------------------------------------------

/**
.Function.getTagKey
..cat:BAM I/O
..signature:getTagKey(tagsDict, idx)
..summary:Return key of a tag by index.
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
..param.idx:Index of the tag whose key to retrieve.
..returns:Infix of the underlying string.
..include:seqan/bam_io.h
*/

template <typename TPos>
inline Infix<Host<BamTagsDict>::Type>::Type
getTagKey(BamTagsDict & tags, TPos idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);
    return infix(host(tags), tags._positions[idx], tags._positions[idx] + 2);
}

template <typename TPos>
inline Infix<Host<BamTagsDict const>::Type>::Type
getTagKey(BamTagsDict const & tags, TPos idx)
{
    return getTagKey(const_cast<BamTagsDict &>(tags), idx);
}

// ----------------------------------------------------------------------------
// Function findTagKey()
// ----------------------------------------------------------------------------

/**
.Function.findTagKey
..cat:BAM I/O
..signature:findTagKey(idx, tagsDict, name)
..summary:Return key of a tag by index.
..param.idx:Index of the tag with the given key.
...type:nolink:$unsigned$
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
..param.name:Name of the key to find.
...type:Shortcut.CharStrin
..returns:$bool$, indicating whether such a key could be found.
..include:seqan/bam_io.h
*/

inline bool
findTagKey(unsigned & idx, BamTagsDict & tags, CharString const & name)
{
    for (idx = 0; idx < length(tags); ++idx)
        if (getTagKey(tags, idx) == name)
            return true;
    return false;
}

inline bool
findTagKey(unsigned & idx, BamTagsDict const & tags, CharString const & name)
{
    return findTagKey(idx, const_cast<BamTagsDict &>(tags), name);
}

// ----------------------------------------------------------------------------
// Function getTagValue()
// ----------------------------------------------------------------------------

/**
.Function.getTagValue
..cat:BAM I/O
..signature:getTagValue(tagsDict, idx)
..param.tagsDict:The @Class.BamTagsDict@ to retrieve data from.
..param.idx:Index of the tag whose value to retrieve.
..returns:@Shortcut.CharString@ with the raw tags data.
..remarks:Note that you will get $<type char> + payload$ in case of @Class.BamTagsDict@.
..include:seqan/bam_io.h
*/

template <typename TIdx>
inline CharString
getTagValue(BamTagsDict & tags, TIdx idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);
        
    typedef typename Position<CharString>::Type TPos;
    TPos beginPos = tags._positions[idx] + 2;
    TPos endPos = beginPos + 1;
    
    char theType = getTagType(tags, idx);
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
        int typeSize = getBamTypeSize(c);
        SEQAN_ASSERT_GT(typeSize, 0);
        endPos += 5 + len * typeSize;
    }
    else
    {
        endPos += getBamTypeSize(theType);
    }
    
    return infix(host(tags), beginPos, endPos);
}

template <typename TPos>
inline CharString //Infix<Host<BamTagsDict const>::Type>::Type
getTagValue(BamTagsDict const & tags, TPos idx)
{
    return getValue(const_cast<BamTagsDict &>(tags), idx);
}

// ----------------------------------------------------------------------------
// Function extractValue()
// ----------------------------------------------------------------------------

/**
.Function.extractValue
..cat:BAM I/O
..signature:extractValue(dest, tags, idx)
..summary:Extract and cast "atomic" value from tags string with index $idx$.
..param.dest:The variable to write the value to.
...remarks:The value is first copied in a variable of the type indicated in the BAM file. Then it is cast into the type of $dest$.
..param.tags:Raw tags string as in BAM.
...type:Shortcut.CharString
..params.idx:Index of the tag in the tag list.
..returns:$bool$, indicating the success.
..remarks:The function only works for atomic types such as $int$, not for $char*$ or arrays.
..see:Function.getTagValue
..include:seqan/bam_io.h
*/

template <typename TDest, typename TIdx>
inline bool
extractValue(TDest & dest, BamTagsDict & tags, TIdx idx)
{
    if (!hasIndex(tags))
        buildIndex(tags);

    char typeC = host(tags)[tags._positions[idx] + 2];
    if (typeC == 'A')
    {
        char x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 1);
        dest = x;
    }
    else if (typeC == 'c')
    {
        __int8 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 1);
        dest = x;
    }
    else if (typeC == 'C')
    {
        __uint8 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 1);
        dest = x;
    }
    else if (typeC == 's')
    {
        __int16 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 2);
        dest = x;
    }
    else if (typeC == 'S')
    {
        __uint16 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 2);
        dest = x;
    }
    else if (typeC == 'i')
    {
        __int32 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 4);
        dest = x;
    }
    else if (typeC == 'I')
    {
        __uint32 x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 4);
        dest = x;
    }
    else if (typeC == 'f')
    {
        float x = 0;
        char * ptr = reinterpret_cast<char *>(&x);
        memcpy(ptr, &host(tags)[tags._positions[idx] + 3], 4);
        dest = x;
    }
    else // variable sized type or invald
    {
        return false;
    }
    return true;
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_TAGS_DICT_H_
