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
// Wrapper for zlib streams.  We use a wrapper because gzFile is an alias
// to void* which could be used in other places, too.
// ==========================================================================

#include <zlib.h>

#ifndef SEQAN_STREAM_STREAM_GZ_FILE_H_
#define SEQAN_STREAM_STREAM_GZ_FILE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.GZ File Stream
..cat:Input / Output
..signature:Stream<GZFile>
..general:Class.Stream
..summary:Adaption from $gzFile$ of $<zlib.h>$ to streams.
..remarks:This is only available if @Macro.SEQAN_HAS_ZLIB@ is set to 1.
..remarks:Not default and copy constructable, not assignable.
..include:seqan/stream.h
 */

template <>
class Stream<GZFile>
{
public:
    gzFile & _gzFile;

    Stream(gzFile & gzFile) : _gzFile(gzFile) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference<Stream<GZFile> >
{
    typedef z_off_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<Stream<GZFile> >
{
    typedef z_off_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<Stream<GZFile> >
{
    typedef z_off_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value<Stream<GZFile> >
{
    typedef char Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<Stream<GZFile>, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<GZFile>, Tell>
{
    typedef True Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------

inline int
streamPeek(char & c, Stream<GZFile> & stream)
{
    int x = gzgetc(stream._gzFile);
    if (x < 0)
        return x;
    c = x;
    gzungetc(c, stream._gzFile);
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, Stream<GZFile> & stream)
{
    int x = gzgetc(stream._gzFile);
    if (x < 0)
        return x;
    c = x;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(Stream<GZFile> & stream)
{
    // std::cerr << "gzeof(stream._gzFile) == " << gzeof(stream._gzFile) << std::endl;
    return gzeof(stream._gzFile) != 0;
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(Stream<GZFile> & stream)
{
    int result = 0;
    gzerror(stream._gzFile, &result);
    // Anything >= is OK/stream end etc.
    return result < 0;
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, Stream<GZFile> & stream, size_t maxLen)
{
    return gzread(stream._gzFile, target, maxLen);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(Stream<GZFile> & stream, char const c)
{
    return gzwrite(stream._gzFile, &c, 1) != 1;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(Stream<GZFile> & stream, char const * source, size_t count)
{
    return gzwrite(stream._gzFile, const_cast<char const *>(source), count);
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

///.Function.streamFlush.remarks:If the stream is of type @Spec.GZ File Stream@ then this function calls $gzflush()$ with $Z_SYNC_FLUSH$. Note that many flush calls to such compressed streams reduce the compression rate.

inline int
streamFlush(Stream<GZFile> & stream)
{
    return gzflush(stream._gzFile, Z_SYNC_FLUSH);
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
streamSeek(Stream<GZFile> & stream, long int delta, int origin)
{
    SEQAN_ASSERT_NEQ(origin, SEEK_END);  // Not supported.
    return gzseek(stream._gzFile, delta, origin) < 0;
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position<Stream<GZFile> >::Type
streamTell(Stream<GZFile> & stream)
{
    return gztell(stream._gzFile);
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_STREAM_GZ_FILE_H_
