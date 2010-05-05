/*
  Definition of the weighted interval type with supporting functions,
  i.e. for I/O.
 */

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


// A simple record class that stores the information of a weighted
// interval for WIT files.
struct WitRecord {
    // Name of the read in question.
    CharString readName;

    // Id of the read, possibly not set.
    size_t readId;

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
    WitRecord() {}

    // Complete constructor for all properties.
    WitRecord(CharString const & _readName, int const & _distance,
              CharString const & _contigName, bool const & _isForward,
              size_t const & _firstPos, size_t const & _lastPos)
            : readName(_readName), readId(0), distance(_distance),
              contigName(_contigName), contigId(0), isForward(_isForward),
              firstPos(_firstPos), lastPos(_lastPos) {}

    // Lexicographic comparison.
    bool operator<(WitRecord const & other) const {
        if (readId < other.readId)
            return true;
        if (readId == other.readId and distance < other.distance)
            return true;
        if (readId == other.readId and distance == other.distance and
            contigId < other.contigId)
            return true;
        if (readId == other.readId and distance == other.distance and
            contigId == other.contigId and firstPos < other.firstPos)
            return true;
        if (readId == other.readId and distance == other.distance and
            contigId == other.contigId and firstPos == other.firstPos and
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


// Output-stream operator for WitRecord objects.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.
template <typename TStream>
TStream & operator<<(TStream & stream, WitRecord const & record) {
    stream << record.readName << "\t"
           << record.distance << "\t"
           << record.contigName << "\t"
           << (record.isForward ? "F" : "R") << "\t"
           << record.firstPos << "\t"
           << record.lastPos << "\t";
    return stream;
}


// Write WIT header line ("@HD\tVN:1.0") to stream.
//
// TStream -- type of the stream.
//
// stream -- Stream to write to.
template <typename TStream>
TStream &writeWitHeader(TStream & stream) {
    stream << "@HD\tVN:1.0" << std::endl;
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
    return stream << record << std::endl;
}


// Read WIT header line ("@HD\tVN:1.0") from stream.
//
// TStream -- type of the stream.
// TChar   -- type of lookahead character.
//
// stream -- Stream to read from.
// c -- Lookahead character.
template <typename TStream, typename TChar>
void readWitHeader(TStream &stream, TChar &c) {
    CharString tmp;
    // Read "@HD".
    c = _streamGet(stream);
    tmp = _parse_readWordUntilWhitespace(stream, c);
    if (tmp != CharString("@HD"))
        std::cerr << "WARNING: File did not begin with \"@HD\", was: \"" << tmp << "\"" << std::endl;
    // Skip "\t".
    _parse_skipWhitespace(stream, c);
    // Read "VN:1.0".
    tmp = _parse_readWordUntilWhitespace(stream, c);
    if (tmp != CharString("VN:1.0"))
        std::cerr << "WARNING: Version is not \"VN:1.0\"" << std::endl;
    // Skip to and after end of line.
    _parse_skipLine(stream, c);
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
template <typename TStream, typename TChar>
bool readWitRecord(TStream & stream, WitRecord & record, TChar & c) {
    String<char> tmp;

    // Maybe skip comments.
    while (not _streamEOF(stream) and c == '#')
        _parse_skipLine(stream, c);
    if (_streamEOF(stream))
        return false;

    // Read line.
    _parse_readIdentifier(stream, record.readName, c);
    _parse_skipWhitespace(stream, c);
    record.distance = _parse_readNumber(stream, c);
    _parse_skipWhitespace(stream, c);
    _parse_readIdentifier(stream, record.contigName, c);
    _parse_skipWhitespace(stream, c);
    record.isForward = (_parse_readChar(stream, c) == 'F');
    _parse_skipWhitespace(stream, c);
    record.firstPos = _parse_readNumber(stream, c);
    _parse_skipWhitespace(stream, c);
    record.lastPos = _parse_readNumber(stream, c);
    
    // Skip to and after end of line.
    _parse_skipLine(stream, c);
    return true;
}

#endif  // WIT_BUILDER_INTERVALS_H_
