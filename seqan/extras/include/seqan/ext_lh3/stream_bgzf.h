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
// Wrapper for BGZF streams.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_EXT_LH3_STREAM_BGZF_H_
#define EXTRAS_INCLUDE_SEQAN_EXT_LH3_STREAM_BGZF_H_

// TODO(holtgrew): Check file for truncation via function.

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Bgzf_;
typedef Tag<Bgzf_> Bgzf;
template <> class Stream<Bgzf>;
inline void close(Stream<Bgzf> & stream);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.BGZF Stream
..cat:Input / Output
..signature:Stream<Bgzf>
..general:Class.Stream
..summary:Adaption from $BZGf *$ of $bgzf.h$ to streams.
..remarks:Not copy constructable, not assignable.
..include:seqan/stream.h
..example.code:
Stream<Bgzf> stream(bgzf);
if (!open("myfile.bam", "r"))
    return 1;  // error
// ... work
 */

struct Bgzf_;
typedef Tag<Bgzf_> Bgzf;

template <>
class Stream<Bgzf>
{
public:
    bool _bgzfOwned;
    BGZF * _bgzf;
    int _error;

    Stream() : _bgzfOwned(false), _bgzf(0), _error(0) {}
    
    Stream(BGZF * bgzf) : _bgzfOwned(false), _bgzf(bgzf), _error(0) {}

    ~Stream()
    {
        if (this->_bgzfOwned)
            close(*this);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================


// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <>
struct Difference<Stream<Bgzf> >
{
    typedef __int64 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <>
struct Position<Stream<Bgzf> >
{
    typedef __int64 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <>
struct Size<Stream<Bgzf> >
{
    typedef __int64 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <>
struct Value<Stream<Bgzf> >
{
    typedef char Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsInput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, IsInput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, IsOutput>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, IsOutput>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasPeek>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, HasPeek>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, HasFilename>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, HasFilename>
{
    typedef False Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Seek<TSpec> >
// ----------------------------------------------------------------------------

template <typename TSpec>
struct HasStreamFeature<Stream<Bgzf>, Seek<TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStreamFeature<, Tell>
// ----------------------------------------------------------------------------

template <>
struct HasStreamFeature<Stream<Bgzf>, Tell>
{
    typedef True Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.open
..signature:open(stream, fileName, mode)
..param.stream:Stream to open.
...type:Class.Stream
..param.mode:Mode string for opening the file, e.g. $"r"$, $"w"$.
...remarks:Do not specify for @Spec.Char Array Stream@.
...type:nolink:$char const *$
 */

inline bool
open(Stream<Bgzf> & stream, char const * filename, char const * mode)
{
    if (stream._bgzfOwned)
        close(stream);
    stream._bgzf = bgzf_open(filename, mode);
    stream._bgzfOwned = true;
    return stream._bgzf != 0;
}

// ----------------------------------------------------------------------------
// Function close()
// ----------------------------------------------------------------------------

/**
.Function.close
..signature:close(stream)
..param.stream:Stream to close.
...type:Class.Stream
 */

inline void
close(Stream<Bgzf> & stream)
{
    if (stream._bgzf == 0)
        return;
    int res = bgzf_close(stream._bgzf);
    (void) res;
    stream._bgzf = 0;
    return;
}

// ----------------------------------------------------------------------------
// Function streamPeek()
// ----------------------------------------------------------------------------

inline int
streamPeek(char & c, Stream<Bgzf> & stream)
{
    int x = bgzf_peek(stream._bgzf);
    if (x < 0)
    {
        stream._error = x;
        return x;
    }
    c = x;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamReadChar()
// ----------------------------------------------------------------------------

inline int
streamReadChar(char & c, Stream<Bgzf> & stream)
{
    int x = bgzf_getc(stream._bgzf);
    if (x < 0)
    {
        stream._error = x;
        return x;
    }
    c = x;
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamEof()
// ----------------------------------------------------------------------------

inline bool
streamEof(Stream<Bgzf> & stream)
{
    // Peek once to set EOF state.
    char c;
    stream._error = streamPeek(c, stream);
    (void)c;
    return stream._error == -1;
}

// ----------------------------------------------------------------------------
// Function streamError
// ----------------------------------------------------------------------------

inline int
streamError(Stream<Bgzf> & stream)
{
    if (stream._error == -1)
        return 0;
    return stream._error;
}

// ----------------------------------------------------------------------------
// Function streamReadBlock()
// ----------------------------------------------------------------------------

inline size_t
streamReadBlock(char * target, Stream<Bgzf> & stream, size_t maxLen)
{
    int x = bgzf_read(stream._bgzf, target, maxLen);
    if (x < 0)
    {
        stream._error = x;
        return 0;
    }
    return x;
}

// ----------------------------------------------------------------------------
// Function streamWriteChar
// ----------------------------------------------------------------------------

inline int
streamWriteChar(Stream<Bgzf> & stream, char const c)
{
    int x = bgzf_write(stream._bgzf, &c, 1);
    if (x < 0)
    {
        stream._error = x;
        return x;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function streamWriteBlock()
// ----------------------------------------------------------------------------

inline size_t
streamWriteBlock(Stream<Bgzf> & stream, char const * source, size_t count)
{
    int x = bgzf_write(stream._bgzf, source, count);
    if (x < 0)
    {
        stream._error = x;
        return 0;
    }
    return count;
}

// ----------------------------------------------------------------------------
// Function streamFlush()
// ----------------------------------------------------------------------------

///.Function.streamFlush.remarks:If the stream is of type @Spec.GZ File Stream@ then this function calls $gzflush()$ with $Z_SYNC_FLUSH$. Note that many flush calls to such compressed streams reduce the compression rate.

inline int
streamFlush(Stream<Bgzf> & stream)
{
    return bgzf_flush(stream._bgzf);
}

// ----------------------------------------------------------------------------
// Function streamSeek()
// ----------------------------------------------------------------------------

inline int
streamSeek(Stream<Bgzf> & stream, __uint64 delta, int origin)
{
    SEQAN_ASSERT_EQ(origin, SEEK_SET);  // Only SET is supported.
    return bgzf_seek(stream._bgzf, delta, origin) < 0;
}

// ----------------------------------------------------------------------------
// Function streamTell()
// ----------------------------------------------------------------------------

inline Position<Stream<Bgzf> >::Type
streamTell(Stream<Bgzf> & stream)
{
    return bgzf_tell(stream._bgzf);
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_EXT_LH3_STREAM_BGZF_H_
