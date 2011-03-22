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
// The Single-Pass Record Reader using the Memory Mapped String
// specialization.  It uses that data does not have to be read explicitely but
// is automatically loaded by the OS on pagefaults.
// ==========================================================================

#ifndef SEQAN_STREAM_RECORD_READER_SINGLE_MMAP_H_
#define SEQAN_STREAM_RECORD_READER_SINGLE_MMAP_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Also use mmadvise to mark used memory as such?

// TODO(holtgrew): This could easily be adjusted to work for any string specialization by adding another layer, signature would then be RecordReader<TFile, SinglePass<StringReader<TStringSpec> > >.

/**
.Spec.Single-Pass MMap RecordReader
..cat:Input / Output
..general:Spec.Single-Pass RecordReader
..summary:Record reader specialization for single-pass reading from memory mapped files.
..signature:RecordReader<TStream, SinglePass<Mapped> >
..param.TStream:The @Concept.Stream@ type to work on.
..remarks:This record reader does not have any buffers but uses the memory mapped string directly.
..remarks:Is not default or copy constructable.
..remarks:The buffer size is the granularity in which @Function.mmapAdvise@ will be called.
..see:Spec.MMap String
..include:seqan/stream.h
 */

template <typename TMMapString>
class RecordReader<TMMapString, SinglePass<Mapped> >
{
public:
    typedef TMMapString TString;
    typedef typename Iterator<TString, Standard>::Type TIter;
    typedef typename Size<TString>::Type TSize;

    TString & _string;
    TIter _current, _end;
    TSize _bufferSize;

    enum {
        OK = 0,
        INVALID_FORMAT
    };

    RecordReader(TString & string)
            : _string(string), _current(begin(string)), _end(end(string)),
              _bufferSize(BUFSIZ)
    {}

    RecordReader(TString & string, unsigned bufferSize)
            : _string(string), _current(begin(string)), _end(end(string)),
              _bufferSize(bufferSize)
    {}

private:
    // No default or copy constructor.
    RecordReader() {}
    RecordReader(RecordReader const &) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function hasMore()
// ----------------------------------------------------------------------------

template <typename TMMapString>
inline bool
hasMore(RecordReader<TMMapString, SinglePass<Mapped> > & recordReader)
{
    return recordReader._current != recordReader._end;
}

// ----------------------------------------------------------------------------
// Function resultCode()
// ----------------------------------------------------------------------------

template <typename TMMapString>
inline int
resultCode(RecordReader<TMMapString, SinglePass<Mapped> > & /*recordReader*/)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TMMapString>
inline bool
goNext(RecordReader<TMMapString, SinglePass<Mapped> > & recordReader)
{
    SEQAN_ASSERT(recordReader._current != recordReader._end);
    recordReader._current += 1;
    return true;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TMMapString>
inline char
value(RecordReader<TMMapString, SinglePass<Mapped> > & recordReader)
{
    SEQAN_ASSERT(recordReader._current != recordReader._end);
    return *recordReader._current;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_RECORD_READER_SINGLE_MMAP_H_
