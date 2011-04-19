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
// Adaption for <cstdio> streams: std::FILE * to the stream concept.
// ==========================================================================

#include <cstdio>

#ifndef SEQAN_STREAM_ADAPT_CSTIO_H_
#define SEQAN_STREAM_ADAPT_CSTIO_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Adaption."FILE *"
..summary:Adaption from $FILE *$ of $<cstdio>$ to streams.
 */

// ============================================================================
// Metafunctions
// ============================================================================

/* // Clashes with definition of these metafunctions in file module.
// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference<FILE *>
{
    typedef long int Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<FILE *>
{
    typedef long int Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<FILE *>
{
    typedef long int Type;
};
*/

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<FILE *, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<FILE *, Tell>
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
streamPeek(char & c, FILE * stream)
{
    c = fgetc(stream);
    ungetc(c, stream);
    if (c == EOF)
        return EOF;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, FILE * stream)
{
    c = fgetc(stream);
    if (c == EOF)
        return EOF;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(FILE * stream)
{
    return ::std::feof(stream) != 0;
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(FILE * stream)
{
    return ::std::ferror(stream);
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, FILE * stream, size_t maxLen)
{
    return ::std::fread(target, sizeof(char), maxLen, stream);
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(FILE * stream, char const c)
{
    int x = ::std::fputc(c, stream);
    if (x == EOF)
        return EOF;
    return c != x;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(FILE * stream, char const * source, size_t count)
{
    return ::std::fwrite(source, sizeof(char), count, stream);
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

inline int
streamPut(FILE * stream, char const c)
{
    return streamWriteChar(stream, c);
}

inline char const *
_streamPutChar(char const*/**/)
{
    return "%s";
}

inline char const *
_streamPutChar(int const/**/)
{
    return "%d";
}

inline char const *
_streamPutChar(unsigned int const/**/)
{
    return "%u";
}

inline char const *
_streamPutChar(long const/**/)
{
    return "%D";
}

inline char const *
_streamPutChar(unsigned long const/**/)
{
    return "%U";
}

inline char const *
_streamPutChar(float const/**/)
{
    return "%.2f"; 
}

inline char const *
_streamPutChar(double const/**/)
{
    return "%.2lf";
}

// TODO(h4nn3s) according to man fprintf's point character is locale dependent,
// maybe overload for doubles and floats to avoid that?
template <typename TSource>
inline int
streamPut(FILE * stream, TSource const & source)
{
    int result = fprintf(stream, _streamPutChar(source), source);
    if (result == -1)
        return errno;
    return 0;
}


// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

inline int
streamFlush(FILE * stream)
{
    return ::std::fflush(stream);
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
streamSeek(FILE * stream, long int delta, int origin)
{
    return ::std::fseek(stream, delta, origin);
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position<FILE *>::Type
streamTell(FILE * stream)
{
    return ::std::ftell(stream);
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_ADAPT_CSTIO_H_
