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
// Wrapper to BZFILE * that fulfills the Stream concept.  Note that, different
// from zlib, bzlib does not support seek, tell or position on BZFILE*.  We do
// not emulate support for this either.  The most case in SeqAn is sequential
// reading and writing and we have full support for this.  We do not support
// peek either, however.  There is no flush(), it is a null implementation but
// should never be required for the use cases of bzlib in SeqAn.
// ==========================================================================

#include <bzlib.h>

#ifndef SEQAN_STREAM_BZ2_FILE_WRAPPER_H_
#define SEQAN_STREAM_BZ2_FILE_WRAPPER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.BZ2 File Stream
..cat:Input / Output
..signature:Stream<BZ2Stream>
..general:Class.Stream
..summary:Wrapper for $BZFILE *$ streams from bzlib.
..remarks:This is only available if @Macro.SEQAN_HAS_ZLIB@ is set to 1.
..remarks:Not default and copy constructable.
..include:seqan/stream.h
 */

template <>
class Stream<BZ2File>
{
public:
    BZFILE * _file;
    int _error;

    Stream(BZFILE * file) : _file(file), _error(0)
    {}

private:
    // Disable default, copy construction and assignment.
    Stream() {}
    Stream(Stream const & /*other*/) {}
    Stream & operator=(Stream const & /*other*/) { return *this; }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference<Stream<BZ2File> >
{
    // bzlib streams rely on FILE * streams
    typedef Difference<Stream<FILE *> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<Stream<BZ2File> >
{
    // bzlib streams rely on FILE * streams
    typedef Position<Stream<FILE *> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<Stream<BZ2File> >
{
    // bzlib streams rely on FILE * streams
    typedef Size<Stream<FILE *> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value<Stream<BZ2File> >
{
    // bzlib streams rely on FILE * streams
    typedef Value<Stream<FILE *> >::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, HasPeek>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<Stream<BZ2File>, Seek<TSpec> >
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<BZ2File>, Tell>
{
    typedef False Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(Stream<BZ2File> & stream)
{
    // std::cerr << "stream._error == " << stream._error << std::endl;
    // std::cerr << "BZ_FINISH == " << BZ_STREAM_END << std::endl;
    return stream._error == BZ_STREAM_END;
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, Stream<BZ2File> & stream)
{
    if (streamEof(stream))
        return 1;
    return BZ2_bzRead(&stream._error, stream._file, &c, 1) != 1;
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(Stream<BZ2File> & stream)
{
    // Anything >= means OK.
    if (stream._error < 0)
        return stream._error;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, Stream<BZ2File> & stream, size_t maxLen)
{
    return BZ2_bzRead(&stream._error, stream._file, target, maxLen);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(Stream<BZ2File> & stream, char const c)
{
    BZ2_bzWrite(&stream._error, stream._file, const_cast<char *>(&c), 1);
    return stream._error;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(Stream<BZ2File> & stream, char const * source, size_t count)
{
    BZ2_bzWrite(&stream._error, stream._file, const_cast<char *>(source), count);
    if (stream._error)
        return 0;
    else
        return count;
}


// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

inline int
streamPut(Stream<BZ2File> & stream, char const c)
{
    return streamWriteChar(stream, c);
}

inline int
streamPut(Stream<BZ2File> & stream, char const * source)
{
    return (streamWriteBlock(stream, source, strlen(source))
                == strlen(source) )  ?   0 : 1;
}

template <typename TSpec>
inline int
streamPut(Stream<BZ2File> & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source))
                == length(source))  ?   0 : 1;
}

template <typename TSource>
inline int
streamPut(Stream<BZ2File> & stream, TSource const & source)
{
    char buffer[1024] = "";
    ::std::stringstream s;

    s << source;
    if (s.fail())
        return s.fail();

    s >> buffer;
    if (s.fail())
        return s.fail();

    buffer[1023] = 0;

    return (streamWriteBlock(stream, buffer, strlen(buffer))
                == strlen(buffer) )  ?   0 : 1;
    //TODO(h4nn3s): might not be the fastest way, eh?
}


// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

///.Function.streamFlush.remarks:If the stream is of type @Spec.BZ2 File Stream@ then this function does nothing.

inline int
streamFlush(Stream<BZ2File> & /*stream*/)
{
    // Null implementation, function is not so important in the use cases for
    // SeqAn anyway.
    return 0;
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_BZ2_FILE_WRAPPER_H_
