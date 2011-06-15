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

/**
.ShortCut.SeqStreamFormats
..cat:Input/Output
..summary:A Tag list of the currently implemented Sequence-Formats (in RecordReader/Stream-IO)
..signature:SeqStreamFormats
..shortcutfor:Tag.TagList
...signature:TagList<Fastq, TagList<Fasta > >
..include:seqan/stream.h
*/
typedef TagList<Fastq, TagList<Fasta > >    SeqStreamFormats;

/**
.ShortCut.AutoSeqStreamFormat
..cat:Input/Output
..summary:A TagSelector for @ShortCut.SeqStreamFormats@. T list of the currently implemented Sequence-Formats (in RecordReader/Stream-IO)
..signature:AutoSeqStreamFormat
..shortcutfor:Class.TagSelector
..shortcutfor:Tag.TagList
...signature:TagSelector<TagList<Fastq, TagList<Fasta > > >
..remarks:can be passed to @Function.checkStreamFormat@ and will offer the index of the detected FileFormat in its member tagId
..see:ShortCut.SeqStreamFormats
..see:Function.checkStreamFormat
..include:seqan/stream.h
*/
typedef TagSelector<SeqStreamFormats>       AutoSeqStreamFormat;


/**
.Class.LimitRecordReaderInScope
..cat:Input/Output
..summary:manipulates a @Class.RecordReader@ -Object so that it operates only on one buffer
..signature:LimitRecordReaderInScope<TStream, TSpec>
..param.TStream:The @Concept.Stream@ of the @Class.RecordReader@.
...type:Concept.Stream
..param.TSpec:The specialization of the @Class.RecordReader@.
...type:Class.RecordReader
..see:Class.RecordReader
..see:Function.checkStreamFormat
..include:seqan/stream.h
..remarks:This class is intended for situations, where you do not wish the RecordReader to rebuffer and where you wish to return to the original reading position after reading, e.g. when detecting the file format of the stream.
..remarks:It is used by passing the RecordReader-object on construction (this already does the necessary changes in the RecordReader). Upon deconstruction of this object, the RecordReader is reset to its original state, including all iterators.
..remarks:This works on all RecordReader-objects, independent of the underlying stream-object. It also works, if the underlying stream does not support seeking.
..include:seqan/stream.h
*/

template <typename TStream, typename TPass>
class LimitRecordReaderInScope
{
public:
    RecordReader<TStream, TPass> & _recordreader;
    typename RecordReader<TStream, TPass>::TIter _currentBeforeSuspend;
    // the following is only needed for MMapStrings
    typename RecordReader<TStream, TPass>::TIter _endBeforeSuspend;

    LimitRecordReaderInScope(RecordReader<TStream, TPass> & reader)
            : _recordreader(reader),
              _currentBeforeSuspend(reader._current),
              _endBeforeSuspend(reader._end)
    {
        _suspendRefill(TPass());
    }

    ~LimitRecordReaderInScope()
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
// Function checkStreamFormat()
// ----------------------------------------------------------------------------

/**
.Funtion.checkStreamFormat
..cat:Input/Output
..summary:check whether the data provided by reader is (one of) the specified format(s).
..signature:checkStreamFormat(TRecordReader & reader, TTag const &)
...param.reader:The @Class.RecordReader@ to read from
...param.TTag:The tag to check against.
..signature:checkStreamFormat(TRecordReader & reader, TagSelector<TTagList> & formats)
...param.reader:The @Class.RecordReader@ to read from
...param.formats:A @Tag.TagSelector@ object that contains the list of tags to check and provides a tagId member with index of the detected tag.
..returns: $True$ if (one of) the specified Tag(s) tested positive and $False$ otherwise
...type:nolink:$bool$
..remarks:With the help of @Class.LimitRecordReaderInScope@ these functions do
not (permanently) alter the position in the stream.
..remarks:The tagId-member of the TagSelector holds the index in inside-to-outside order and begins counting at one. E.g. The Index of FASTQ in TagList<Fastq, TagList<Fasta > > would be 2
..include:seqan/stream.h
*/

template < typename TRecordReader >
inline bool
checkStreamFormat(TRecordReader & reader, TagSelector<> & formats)
{
    (void)reader;
    formats.tagId = 0;
    return false;
}

template <typename TRecordReader, typename TTagList >
inline bool
checkStreamFormat(TRecordReader & reader, TagSelector<TTagList> & formats)
{
    if (checkStreamFormat(reader, typename TTagList::Type()))
    {
        if (formats.tagId == 0)
        {
            // if tagId == 0 then store detected format
            formats.tagId = LENGTH<TTagList>::VALUE;
            return true;
        } else
            // if tagId != 0 then compare detected format with tagId
            return formats.tagId == LENGTH<TTagList>::VALUE;
    }
    return checkStreamFormat(reader,
                              static_cast<
                                typename TagSelector<TTagList>::Base &
                                         > (formats) );
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_ADAPT_FSTREAM_H_