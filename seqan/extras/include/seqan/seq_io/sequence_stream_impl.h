// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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

// TODO(holtgrew): Add #if guards for bz2 and zlib.

#ifndef EXTRAS_INCLUDE_SEQAN_SEQ_IO_SEQ_IO_IMPL_H_
#define EXTRAS_INCLUDE_SEQAN_SEQ_IO_SEQ_IO_IMPL_H_

#include <memory>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ---------------------------------------------------------------------------
// Class SeqIOFileType_
// ---------------------------------------------------------------------------

struct SeqIOFileType_
{
    // Used for storing the file type guess.
    enum Type
    {
        FILE_TYPE_ERROR,
        FILE_TYPE_TEXT,
        FILE_TYPE_GZ,
        FILE_TYPE_BZ2
    };
};

// ---------------------------------------------------------------------------
// Class SeqIOFileFormat_
// ---------------------------------------------------------------------------

struct SeqIOFileFormat_
{
    enum Type
    {
        FILE_FORMAT_ERROR,
        FILE_FORMAT_FASTA,
        FILE_FORMAT_FASTQ,
        FILE_FORMAT_EMBL
    };
};

// ----------------------------------------------------------------------------
// Class SequenceStreamImpl_
// ----------------------------------------------------------------------------

// Implementation of EasySeqIO class.
//
// This class does the multiplexing from settings known at runtime (file format) to types at compile-time.

class SequenceStreamImpl_
{
public:
#if SEQAN_HAS_ZLIB
    std::auto_ptr<Stream<GZFile> > _gzStream;
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
    std::auto_ptr<Stream<BZ2File> > _bz2Stream;
#endif  // #if SEQAN_HAS_BZIP2
    std::auto_ptr<String<char, MMap<> > > _mmapString;

    // TODO(holtgrew): We could get rid of some of these with type erasure on streams and record readers. Would this be enough?

#if SEQAN_HAS_ZLIB
    std::auto_ptr<RecordReader<Stream<GZFile>, SinglePass<> > > _gzReader;
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
    std::auto_ptr<RecordReader<Stream<BZ2File>, SinglePass<> > > _bz2Reader;
#endif  // #if SEQAN_HAS_BZIP2
    std::auto_ptr<RecordReader<String<char, MMap<> >, SinglePass<Mapped> > > _mmapReaderSinglePass;
    std::auto_ptr<RecordReader<String<char, MMap<> >, DoublePass<Mapped> > > _mmapReaderDoublePass;

    SeqIOFileType_::Type _fileType;

    // Whether or not we are at the end of the underlying file.
    bool _atEnd;
    // Whether or not there was no error.
    bool _isGood;
    // Whether or not we open the file for reading only.
    bool _readMode;
    // Whether or not to use double pass record reader.
    bool _doublePass;

    SequenceStreamImpl_(CharString const & filename, SeqIOFileType_::Type fileType, bool readMode, bool doublePass = false) :
        _fileType(fileType), _atEnd(false), _isGood(true), _readMode(readMode), _doublePass(doublePass)
    {
        if (!readMode)
        {
            // Open for writing.
            switch (fileType)
            {
#if SEQAN_HAS_ZLIB
            case SeqIOFileType_::FILE_TYPE_GZ:
            {
                _gzStream.reset(new Stream<GZFile>());
                if (!open(*_gzStream, toCString(filename), "w"))
                    _isGood = false;
            }
            break;

#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
            case SeqIOFileType_::FILE_TYPE_BZ2:
            {
                _bz2Stream.reset(new Stream<BZ2File>());
                if (!open(*_bz2Stream, toCString(filename), "w"))
                    _isGood = false;
            }
            break;

#endif  // #if SEQAN_HAS_BZIP2
            case SeqIOFileType_::FILE_TYPE_TEXT:
            default:      // Error cannot appear here.
            {
                _mmapString.reset(new String<char, MMap<> >());
                if (!open(*_mmapString, toCString(filename), OPEN_CREATE | OPEN_RDWR))
                    _isGood = false;
            }
            break;
            }
        }
        else
        {
            // Open for reading.

            // Open stream or memory mapped file and create MMap reader.
            switch (fileType)
            {
            case SeqIOFileType_::FILE_TYPE_TEXT:
            {
                _mmapString.reset(new String<char, MMap<> >());
                if (!open(*_mmapString, toCString(filename), OPEN_RDONLY))
                {
                    _isGood = false;
                }
                else
                {
                    _doublePass = doublePass;          // Only enable for MMap if set so.
                    if (!doublePass)
                        _mmapReaderSinglePass.reset(new RecordReader<String<char, MMap<> >, SinglePass<Mapped> >(*_mmapString));
                    else
                        _mmapReaderDoublePass.reset(new RecordReader<String<char, MMap<> >, DoublePass<Mapped> >(*_mmapString));
                }
            }
            break;

#if SEQAN_HAS_ZLIB
            case SeqIOFileType_::FILE_TYPE_GZ:
            {
                _gzStream.reset(new Stream<GZFile>());
                if (!open(*_gzStream, toCString(filename), "r"))
                    _isGood = false;
                else
                    _gzReader.reset(new RecordReader<Stream<GZFile>, SinglePass<> >(*_gzStream));
            }
            break;

#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
            case SeqIOFileType_::FILE_TYPE_BZ2:
            {
                _bz2Stream.reset(new Stream<BZ2File>());
                if (!open(*_bz2Stream, toCString(filename), "r"))
                    _isGood = false;
                else
                    _bz2Reader.reset(new RecordReader<Stream<BZ2File>, SinglePass<> >(*_bz2Stream));
            }
            break;

#endif  // #if SEQAN_HAS_BZIP2
            default:
            {
                _isGood = false;
            }
            }
        }
    }

    bool atEnd()
    {
        return _atEnd;
    }

    bool isGood()
    {
        return _isGood;
    }

    // -----------------------------------------------------------------------
    // Function readRecord()
    // -----------------------------------------------------------------------

    // These two template functions implement the record-wise reading of sequence files.

    template <typename TId, typename TSequence, typename TQualities, typename TFormatTag>
    int readRecord(TId & id, TSequence & seq, TQualities & qual, TFormatTag const & tag)
    {
        int res = 0;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            if (!_doublePass)
            {
                res = seqan::readRecord(id, seq, qual, *_mmapReaderSinglePass, tag);
                _atEnd = seqan::atEnd(*_mmapReaderSinglePass);
            }
            else
            {
                res = seqan::readRecord(id, seq, qual, *_mmapReaderDoublePass, tag);
                _atEnd = seqan::atEnd(*_mmapReaderDoublePass);
            }
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
            res = seqan::readRecord(id, seq, qual, *_gzReader, tag);
            _atEnd = seqan::atEnd(*_gzReader);
            break;      // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
            res = seqan::readRecord(id, seq, qual, *_bz2Reader, tag);
            _atEnd = seqan::atEnd(*_bz2Reader);
            break;      // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    template <typename TId, typename TSequence, typename TFormatTag>
    int readRecord(TId & id, TSequence & seq, TFormatTag const & tag)
    {
        int res = 0;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            if (!_doublePass)
            {
                res = seqan::readRecord(id, seq, *_mmapReaderSinglePass, tag);
                _atEnd = seqan::atEnd(*_mmapReaderSinglePass);
            }
            else
            {
                res = seqan::readRecord(id, seq, *_mmapReaderDoublePass, tag);
                _atEnd = seqan::atEnd(*_mmapReaderDoublePass);
            }
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
            res = seqan::readRecord(id, seq, *_gzReader, tag);
            _atEnd = seqan::atEnd(*_gzReader);
            break;      // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
            res = seqan::readRecord(id, seq, *_bz2Reader, tag);
            _atEnd = seqan::atEnd(*_bz2Reader);
            break;      // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    // -----------------------------------------------------------------------
    // Function readBatch()
    // -----------------------------------------------------------------------

    // These two template functions implement the batch-wise reading of sequence files.

    template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TQualities,
              typename TQualSpec, typename TFormatTag>
    int readBatch(StringSet<TId, TIdSpec> & ids,
                  StringSet<TSequence, TSeqSpec> & seqs,
                  StringSet<TQualities, TQualSpec> & quals,
                  unsigned num,
                  TFormatTag const & tag)
    {
        clear(ids);
        clear(seqs);
        clear(quals);

        int res = 0;
        TId id;
        TSequence seq;
        TQualities qual;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            if (!_doublePass)
            {
                for (unsigned i = 0; (res == 0) && (i < num) && !seqan::atEnd(*_mmapReaderSinglePass); ++i)
                {
                    res = seqan::readRecord(id, seq, qual, *_mmapReaderSinglePass, tag);
                    appendValue(ids, id);
                    appendValue(seqs, seq);
                    appendValue(quals, qual);
                }
                _atEnd = seqan::atEnd(*_mmapReaderSinglePass);
            }
            else
            {
                for (unsigned i = 0; (res == 0) && (i < num) && !seqan::atEnd(*_mmapReaderDoublePass); ++i)
                {
                    res = seqan::readRecord(id, seq, qual, *_mmapReaderDoublePass, tag);
                    appendValue(ids, id);
                    appendValue(seqs, seq);
                    appendValue(quals, qual);
                }
                _atEnd = seqan::atEnd(*_mmapReaderDoublePass);
            }
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
        {
            for (unsigned i = 0; (res == 0) && (i < num) && !seqan::atEnd(*_gzReader); ++i)
            {
                res = seqan::readRecord(id, seq, qual, *_gzReader, tag);
                appendValue(ids, id);
                appendValue(seqs, seq);
                appendValue(quals, qual);
            }
            _atEnd = seqan::atEnd(*_gzReader);
        }
        break;          // end of case

#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
        {
            for (unsigned i = 0; (res == 0) && (i < num) && !seqan::atEnd(*_bz2Reader); ++i)
            {
                res = seqan::readRecord(id, seq, qual, *_bz2Reader, tag);
                appendValue(ids, id);
                appendValue(seqs, seq);
                appendValue(quals, qual);
            }
            _atEnd = seqan::atEnd(*_bz2Reader);
        }
        break;          // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TFormatTag>
    int readBatch(StringSet<TId, TIdSpec> & ids,
                  StringSet<TSequence, TSeqSpec> & seqs,
                  unsigned num,
                  TFormatTag const & tag)
    {
        clear(ids);
        clear(seqs);

        int res = 0;
        TId id;
        TSequence seq;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            if (!_doublePass)
            {
                for (unsigned i = 0; (res == 0) && (i < num) && !seqan::atEnd(*_mmapReaderSinglePass); ++i)
                {
                    res = seqan::readRecord(id, seq, *_mmapReaderSinglePass, tag);
                    appendValue(ids, id);
                    appendValue(seqs, seq);
                }
                _atEnd = seqan::atEnd(*_mmapReaderSinglePass);
            }
            else
            {
                for (unsigned i = 0; (res == 0) && (i < num) && !seqan::atEnd(*_mmapReaderDoublePass); ++i)
                {
                    res = seqan::readRecord(id, seq, *_mmapReaderDoublePass, tag);
                    appendValue(ids, id);
                    appendValue(seqs, seq);
                }
                _atEnd = seqan::atEnd(*_mmapReaderDoublePass);
            }
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
        {
            for (unsigned i = 0; (res == 0) && (i < num) && !seqan::atEnd(*_gzReader); ++i)
            {
                res = seqan::readRecord(id, seq, *_gzReader, tag);
                appendValue(ids, id);
                appendValue(seqs, seq);
            }
            _atEnd = seqan::atEnd(*_gzReader);
        }
        break;          // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
        {
            for (unsigned i = 0; (res == 0) && (i < num) && !seqan::atEnd(*_bz2Reader); ++i)
            {
                res = seqan::readRecord(id, seq, *_bz2Reader, tag);
                appendValue(ids, id);
                appendValue(seqs, seq);
            }
            _atEnd = seqan::atEnd(*_bz2Reader);
        }
        break;          // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    // -----------------------------------------------------------------------
    // Function readAll()
    // -----------------------------------------------------------------------

    // These two template functions implement the reading of whole sequence files.

    template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TQualities,
              typename TQualSpec, typename TFormatTag>
    int readAll(StringSet<TId, TIdSpec> & ids,
                StringSet<TSequence, TSeqSpec> & seqs,
                StringSet<TQualities, TQualSpec> & quals,
                TFormatTag const & tag)
    {
        clear(ids);
        clear(seqs);
        clear(quals);

        int res = 0;
        TId id;
        TSequence seq;
        TQualities qual;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            if (!_doublePass)
            {
                while (!seqan::atEnd(*_mmapReaderSinglePass))
                {
                    res = seqan::readRecord(id, seq, qual, *_mmapReaderSinglePass, tag);
                    appendValue(ids, id);
                    appendValue(seqs, seq);
                    appendValue(quals, qual);
                }
                _atEnd = seqan::atEnd(*_mmapReaderSinglePass);
            }
            else
            {
                while (!seqan::atEnd(*_mmapReaderDoublePass))
                {
                    res = seqan::readRecord(id, seq, qual, *_mmapReaderDoublePass, tag);
                    appendValue(ids, id);
                    appendValue(seqs, seq);
                    appendValue(quals, qual);
                }
                _atEnd = seqan::atEnd(*_mmapReaderDoublePass);
            }
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
        {
            while (!seqan::atEnd(*_gzReader))
            {
                res = seqan::readRecord(id, seq, qual, *_gzReader, tag);
                appendValue(ids, id);
                appendValue(seqs, seq);
                appendValue(quals, qual);
            }
            _atEnd = seqan::atEnd(*_gzReader);
        }
        break;          // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
        {
            while (!seqan::atEnd(*_bz2Reader))
            {
                res = seqan::readRecord(id, seq, qual, *_bz2Reader, tag);
                appendValue(ids, id);
                appendValue(seqs, seq);
                appendValue(quals, qual);
            }
            _atEnd = seqan::atEnd(*_bz2Reader);
        }
        break;          // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TFormatTag>
    int readAll(StringSet<TId, TIdSpec> & ids,
                StringSet<TSequence, TSeqSpec> & seqs,
                TFormatTag const & tag)
    {
        clear(ids);
        clear(seqs);

        int res = 0;
        TId id;
        TSequence seq;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            if (!_doublePass)
            {
                while (!seqan::atEnd(*_mmapReaderSinglePass))
                {
                    res = seqan::readRecord(id, seq, *_mmapReaderSinglePass, tag);
                    appendValue(ids, id);
                    appendValue(seqs, seq);
                }
                _atEnd = seqan::atEnd(*_mmapReaderSinglePass);
            }
            else
            {
                while (!seqan::atEnd(*_mmapReaderDoublePass))
                {
                    res = seqan::readRecord(id, seq, *_mmapReaderDoublePass, tag);
                    appendValue(ids, id);
                    appendValue(seqs, seq);
                }
                _atEnd = seqan::atEnd(*_mmapReaderDoublePass);
            }
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
        {
            while (!seqan::atEnd(*_gzReader))
            {
                res = seqan::readRecord(id, seq, *_gzReader, tag);
                appendValue(ids, id);
                appendValue(seqs, seq);
            }
            _atEnd = seqan::atEnd(*_gzReader);
        }
        break;          // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
        {
            while (!seqan::atEnd(*_bz2Reader))
            {
                res = seqan::readRecord(id, seq, *_bz2Reader, tag);
                appendValue(ids, id);
                appendValue(seqs, seq);
            }
            _atEnd = seqan::atEnd(*_bz2Reader);
        }
        break;          // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    // -----------------------------------------------------------------------
    // Function writeRecord()
    // -----------------------------------------------------------------------

    // These two template functions implement the writing of one record-by-record writing.

    template <typename TId, typename TSequence, typename TQualities, typename TFormatTag>
    int writeRecord(TId const & id, TSequence const & seq, TQualities const & qual, TFormatTag const & tag)
    {
        int res = 0;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            res = seqan::writeRecord(*_mmapString, id, seq, qual, tag);
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
            res = seqan::writeRecord(*_gzStream, id, seq, qual, tag);
            break;      // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
            res = seqan::writeRecord(*_bz2Stream, id, seq, qual, tag);
            break;      // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    template <typename TId, typename TSequence, typename TFormatTag>
    int writeRecord(TId const & id, TSequence const & seq, TFormatTag const & tag)
    {
        int res = 0;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            res = seqan::writeRecord(*_mmapString, id, seq, tag);
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
            res = seqan::writeRecord(*_gzStream, id, seq, tag);
            break;      // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
            res = seqan::writeRecord(*_bz2Stream, id, seq, tag);
            break;      // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    // -----------------------------------------------------------------------
    // Function writeAll()
    // -----------------------------------------------------------------------

    // These two template functions implement the writing of multiple strings to a sequence file.

    template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TQualities,
              typename TQualSpec, typename TFormatTag>
    int writeAll(StringSet<TId, TIdSpec> const & ids,
                 StringSet<TSequence, TSeqSpec> const & seqs,
                 StringSet<TQualities, TQualSpec> const & quals,
                 TFormatTag const & tag)
    {
        int res = 0;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            res = seqan::write2(*_mmapString, ids, seqs, quals, tag);
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
            res = seqan::write2(*_gzStream, ids, seqs, quals, tag);
            break;      // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
            res = seqan::write2(*_bz2Stream, ids, seqs, quals, tag);
            break;      // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    template <typename TId, typename TIdSpec, typename TSequence, typename TSeqSpec, typename TFormatTag>
    int writeAll(StringSet<TId, TIdSpec> const & ids,
                 StringSet<TSequence, TSeqSpec> const & seqs,
                 TFormatTag const & tag)
    {
        int res = 0;
        TId id;
        TSequence seq;

        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            res = seqan::write2(*_mmapString, ids, seqs, tag);
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
            res = seqan::write2(*_gzStream, ids, seqs, tag);
            break;      // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
            res = seqan::write2(*_bz2Stream, ids, seqs, tag);
            break;      // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            return 1;
        }

        _isGood = _isGood && (res == 0);
        return res;
    }

    // -----------------------------------------------------------------------
    // Functions for closing/flushing.
    // -----------------------------------------------------------------------

    void flush()
    {
        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            seqan::flush(*_mmapString);
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
            seqan::streamFlush(*_gzStream);
            break;      // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
            seqan::streamFlush(*_bz2Stream);
            break;      // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            break;
        }
    }

    void close()
    {
        switch (_fileType)
        {
        case SeqIOFileType_::FILE_TYPE_TEXT:
            seqan::close(*_mmapString);
            break;      // end of case

#if SEQAN_HAS_ZLIB
        case SeqIOFileType_::FILE_TYPE_GZ:
            seqan::close(*_gzStream);
            break;      // end of case

#endif // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
        case SeqIOFileType_::FILE_TYPE_BZ2:
            seqan::close(*_bz2Stream);
            break;      // end of case

#endif  // #if SEQAN_HAS_BZIP2
        default:
            break;
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_SEQ_IO_SEQ_IO_IMPL_H_
