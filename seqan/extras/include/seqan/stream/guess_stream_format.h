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
// Author: Manuel Holtgrewe <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// File Format detection with Stream Class
// ==========================================================================


#ifndef SEQAN_STREAM_GUESSSTREAMFORMAT_H_
#define SEQAN_STREAM_GUESSSTREAMFORMAT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(h4nn3s): do we want this here?
// typedef
//     TagList<Fastq,
//     TagList<Fasta,
//     TagList<QSeq,
//     TagList<Raw> > > >                      SeqFormats; //IOREV
// typedef TagSelector<SeqFormats>             AutoSeqFormat; //IOREV _doc_
// 


template <typename TStream, typename TPass>
class ReduceToOneBuffer_
{
public:
    RecordReader<TStream, TPass> & _recordreader;
    typename RecordReader<TStream, TPass>::TIter _currentBeforeSuspend;
    // the following is only needed for MMapStrings
    typename RecordReader<TStream, TPass>::TIter _endBeforeSuspend;

    ReduceToOneBuffer_(RecordReader<TStream, TPass> & reader)
            : _recordreader(reader),
              _currentBeforeSuspend(reader._current),
              _endBeforeSuspend(reader._end)
    {
        _suspendRefill(TPass());
    }

    ~ReduceToOneBuffer_()
    {
        _resumeRefillAndReset(TPass());
    }
private:
    inline void
    _suspendRefill(SinglePass<void> const & /* tag */)
    {
        _recordreader._stayInOneBuffer = true;
    }

    inline void
    _suspendRefill(SinglePass<Mapped> const & /* tag */)
    {
        if (_recordreader._end - _recordreader._current > BUFSIZ)
            _recordreader._end = _recordreader._current + BUFSIZ;
    }

    template <typename TSpec>
    inline void
    _suspendRefill(DoublePass<TSpec> const & /* tag */)
    {
        _suspendRefill(SinglePass<TSpec>());
    }


    inline void
    _resumeRefillAndReset(SinglePass<void> const & /* tag */)
    {
        _recordreader._stayInOneBuffer = false;
        if (_currentBeforeSuspend == 0) // there had been no Buffer at start
            _recordreader._current = begin(_recordreader._buffer, Standard());
        else
            _recordreader._current = _currentBeforeSuspend;
    }

    inline void
    _resumeRefillAndReset(DoublePass<void> const & /* tag */)
    {
        _recordreader._stayInOneBuffer = false;
        if (_currentBeforeSuspend == 0) // there had been no Buffer at start
            _recordreader._current = begin(*_recordreader._currentBuffer,
                                           Standard());
        else
            _recordreader._current = _currentBeforeSuspend;
    }

    inline void
    _resumeRefillAndReset(SinglePass<Mapped> const & /* tag */)
    {
        _recordreader._end = _endBeforeSuspend;
        _recordreader._current = _currentBeforeSuspend;
    }

    template <typename TSpec>
    inline void
    _resumeRefillAndReset(DoublePass<TSpec> const & /* tag */)
    {
        _resumeRefillAndReset(SinglePass<TSpec>());
    }
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


// ----------------------------------------------------------------------------
// Function guessStreamFormats()
// ----------------------------------------------------------------------------


//TODO(h4nn3s): test + dddoc
template < typename TRecordReader >
inline bool
checkStreamFormats(TRecordReader & reader, TagSelector<> & format)
{
    format.tagId = 0;
    return false;
}

//TODO(h4nn3s): test + dddoc
template <typename TRecordReader, typename TFileSeq, typename TTagList >
inline bool
checkStreamFormats(TRecordReader & reader, TagSelector<TTagList> & format)
{
    if (checkStreamFormat(reader, typename TTagList::Type()))
    {
        if (format.tagId == 0)
        {
            // if tagId == 0 then store detected format
            format.tagId = LENGTH<TTagList>::VALUE;
            return true;
        } else
            // if tagId != 0 then compare detected format with tagId
            return format.tagId == LENGTH<TTagList>::VALUE;
    }
    return checkStreamFormats(reader,
                              static_cast<
                                typename TagSelector<TTagList>::Base &
                                         > (format) );
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_FSTREAM_H_
