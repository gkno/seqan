// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Definition of the weighted interval type with supporting functions,
// i.e. for I/O.
// ==========================================================================

// TODO(holtgrew): Should maybe better be called LabeledInterval?

#ifndef WIT_BUILDER_INTERVALS_H_
#define WIT_BUILDER_INTERVALS_H_

#include <iostream>

#include <seqan/basic.h>
#include <seqan/misc/misc_parsing.h>
#include <seqan/sequence.h>

using namespace seqan;

// The columns in a WIT file.
#define WIT_COLUMN_NAMES "QNAME\tSCORE\tRNAME\tRDIR\tSPOS\tEPOS"

struct Gsi_;
typedef seqan::Tag<Gsi_> Gsi;

// Dummy only at the moment.
class WitHeader
{
public:
    WitHeader() {}
};


// A simple record class that stores the information of a weighted
// interval for WIT files.
struct WitRecord {
    // Name of the read in question.
    CharString readName;

    enum {
      FLAG_PAIRED = 0x01,
      FLAG_FIRST_MATE = 0x40,
      FLAG_SECOND_MATE = 0x80
    };

    // Flags, 0x01 - paired, 0x40 - first read in pair, 0x80 - second read in
    // pair.
    int flags;

    // Id of the read, possibly not set.
    size_t readId;

    // Original distance of the interval before lowering.
    int originalDistance;

    // Distance associated with the interval.
    int distance;

    // Name of the contig the interval is defined on.
    CharString contigName;

    // Id of the contig, possibly not set.
    size_t contigId;

    // true iff the interval is on the forward strand.
    bool isForward;

    // First position on the interval.
    size_t firstPos;

    // Last position on the interval.
    size_t lastPos;

    // Default constructor.
    WitRecord() : flags(0), readId(0), originalDistance(0), distance(0), contigId(0), isForward(0), firstPos(0), lastPos(0)
    {}

    // Complete constructor for all properties.
    WitRecord(CharString const & _readName, int const & _flags,
              int const & _distance,
              CharString const & _contigName, bool const & _isForward,
              size_t const & _firstPos, size_t const & _lastPos)
            : readName(_readName), flags(_flags), readId(0), originalDistance(_distance), distance(_distance),
              contigName(_contigName), contigId(0), isForward(_isForward),
              firstPos(_firstPos), lastPos(_lastPos) {}

    // Lexicographic comparison.
    bool operator<(WitRecord const & other) const {
        if (readId < other.readId)
            return true;
        if (readId == other.readId && distance < other.distance)
            return true;
        if (readId == other.readId && distance == other.distance &&
            contigId < other.contigId)
            return true;
        if (readId == other.readId && distance == other.distance &&
            contigId == other.contigId && firstPos < other.firstPos)
            return true;
        if (readId == other.readId && distance == other.distance &&
            contigId == other.contigId && firstPos == other.firstPos &&
            lastPos < other.lastPos)
            return true;
        return false;
    }
};


struct WitRecord_Lt_ContigIdReadIdLastPos {
    bool operator()(WitRecord const & a, WitRecord const & b) {
        if (a.contigId < b.contigId)
            return true;
        if (a.contigId == b.contigId && a.readId < b.readId)
            return true;
        if (a.contigId == b.contigId && a.readId == b.readId && a.lastPos < b.lastPos)
            return true;
        return false;
    }
};

void clear(WitRecord & record)
{
    clear(record.readName);
    record.flags = 0;
    record.readId = 0;
    record.originalDistance = 0;
    record.distance = 0;
    clear(record.contigName);
    record.contigId = 0;
    record.isForward = true;
    record.firstPos = 0;
    record.lastPos = 0;
}


// Output-stream operator for WitRecord objects.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.
template <typename TStream>
TStream & operator<<(TStream & stream, WitRecord const & record) {
    stream << record.readName;
    if (record.flags & WitRecord::FLAG_PAIRED) {
      if (record.flags & WitRecord::FLAG_FIRST_MATE) {
        stream << "/0";
      } else if (record.flags & WitRecord::FLAG_SECOND_MATE) {
        stream << "/1";
      } else {
        stream << "/?";
      }
    }
    stream << "\t"
           << record.distance << "\t"
           << record.contigName << "\t"
           << (record.isForward ? "F" : "R") << "\t"
           << record.firstPos << "\t"
           << record.lastPos;
    return stream;
}


// Write WIT header line ("@WIT\tVN:1.0") to stream.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.
template <typename TStream>
TStream &writeWitHeader(TStream & stream) {
//IOREV why does this return a stream reference? is wit more widely used?
    stream << "@WIT\tVN:1.0" << std::endl;
    stream << "@MATES\tSEP:/\tTYPE:01" << std::endl;
    return stream;
}


// Write a WIT comment line ("# %s") to stream.
//
// TStream -- an output stream.
// TString -- the type of the comment string.
//
// stream -- stream to write to.
// str    -- string to write to stream.
template <typename TStream, typename TString>
TStream &writeWitComment(TStream &stream, TString const & str) {
//IOREV why does this return a stream reference? is wit more widely used?
    stream << "# " << str << std::endl;
    return stream;
}


// Write a WIT data line to stream.
//
// TStream -- an output stream.
//
// stream -- stream to write to.
// record -- the WitRecord to write out.
template <typename TStream, typename TString, typename TScore, typename TPos>
TStream &writeWitRecord(TStream & stream, WitRecord const & record) {
//IOREV why does this return a stream reference? is wit more widely used?
    return stream << record << std::endl;
}


// Read WIT header line ("@WIT\tVN:1.0") from record reader.

template <typename TStream, typename TSpec>
int readRecord(WitHeader & header, RecordReader<TStream, TSpec> & reader, Gsi const & /*tag*/)
{
    (void) header;
    
    CharString tmp;
    // Read "@WIT".
    if (readUntilTabOrLineBreak(tmp, reader) != 0)
        return 1;  // Could not read header.
    if (tmp != "@WIT")
        std::cerr << "WARNING: File did not begin with \"@WIT\", was: \"" << tmp << "\"" << std::endl;
    // Skip "\t".
    if (skipChar(reader, '\t') != 0)
        return 1;  // Next char was not TAB.
    // Read "VN:1.0".
    clear(tmp);
    if (readUntilTabOrLineBreak(tmp, reader) != 0)
        return 1;  // Could not read version.
    if (tmp != "VN:1.0")
        std::cerr << "WARNING: Version is not \"VN:1.0\", was: \"" << tmp << "\"" << std::endl;
    // Skip to and after end of line.
    skipLine(reader);
    // Maybe read/skip additional header lines.
    while (!atEnd(reader) && value(reader) == '@')
        if (skipLine(reader) != 0)
            return 1;

    // Skipy any trailing comment lines.
    while (!atEnd(reader) && value(reader) == '#')
        if (skipLine(reader) != 0)
            return 1;
    return 0;
}


// Read WIT record, skipping comments.
//
// TStream -- type of the stream.
// TChar   -- lookahead for the parser.
//
// stream -- Stream to read from.
// record -- the Witrecord to read.
// c -- lookahead for the parser.
//
// Note that the WIT file contains 1-based entries but we subtract 1
// at this location.
//
// Returns true iff the record could be successfully read from the file.

template <typename TStream, typename TSpec>
int readRecord(WitRecord & record, RecordReader<TStream, TSpec> & reader, Gsi const & /*tag*/)
{
    CharString buffer;

    // No more records in file.
    if (atEnd(reader))
        return 1;

    // Read read name.
    clear(record);
    if (readUntilTabOrLineBreak(record.readName, reader) != 0)
        return 1;
    if (value(reader) != '\t')
        return 1;
    skipChar(reader, '\t');

    // Interpret trailing characters for mate-pair identifier in read name.
    if (length(record.readName) >= 2u && record.readName[length(record.readName) - 2] == '/')
    {
        char c = back(record.readName);
        if (c == '0')
            record.flags = WitRecord::FLAG_PAIRED | WitRecord::FLAG_FIRST_MATE;
        else if (c == '1')
            record.flags = WitRecord::FLAG_PAIRED | WitRecord::FLAG_SECOND_MATE;
        else
            return 1;  // Could not interpret trailing mate indicator.
        resize(record.readName, length(record.readName) - 2);
    }

    // Read distance.
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
        return 1;
    if (!lexicalCast2(record.distance, buffer))
        return 1;  // Could not convert distance.
    record.originalDistance = record.distance;
    if (value(reader) != '\t')
        return 1;
    skipChar(reader, '\t');

    // Read contig name.
    if (readUntilTabOrLineBreak(record.contigName, reader) != 0)
        return 1;
    if (value(reader) != '\t')
        return 1;
    skipChar(reader, '\t');

    // Read 'F'/'R'.
    clear(buffer);
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
        return 1;
    if (buffer != "F" && buffer != "R")
        return 1;
    record.isForward = (buffer[0] == 'F');
    if (value(reader) != '\t')
        return 1;
    skipChar(reader, '\t');

    // Read first pos.
    clear(buffer);
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
        return 1;
    if (!lexicalCast2(record.firstPos, buffer))
        return 1;
    if (value(reader) != '\t')
        return 1;
    skipChar(reader, '\t');
    
    // Read last pos.
    clear(buffer);
    if (readUntilTabOrLineBreak(buffer, reader) != 0)
        return 1;
    if (!lexicalCast2(record.lastPos, buffer))
        return 1;
    if (value(reader) != '\t' && value(reader) != '\r' && value(reader) != '\n')  // TODO(holtgrew): \t only here because output buggy.
        return 1;
    if (skipLine(reader) != 0)
        return 1;  // Skip line.

    // Skipy any trailing comment lines.
    while (!atEnd(reader) && value(reader) == '#')
        if (skipLine(reader) != 0)
            return 1;

    return 0;
}

#endif  // WIT_BUILDER_INTERVALS_H_

