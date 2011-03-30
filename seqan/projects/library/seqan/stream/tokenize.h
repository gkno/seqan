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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Functions for tokenizing streams
// ==========================================================================

#ifndef SEQAN_STREAM_TOKENIZE_H
#define SEQAN_STREAM_TOKENIZE_H

// NOTHING HERE IS TESTED YET


#define IO_TOKENIZE_EOF_BEFORE_SUCCESS 1024
#define IO_TOKONIZE_NO_SUCCESS 1025

namespace seqan {

    
// ----------------------- Helper structs ------------------------------------
struct Whitespace_;
struct Blank_;
struct Char_;
struct Digit_;
struct Alpha_;
struct AlphaNum_;
struct UnixEOL_;
struct BackslashR_;
struct Graph_;


// template <typename TSpec>
// struct CharCompare_;


template <>
struct CharCompare_<Whitespace_>
{
    static inline int _compare(const int c)
    {
        return isspace(c);
    }
}

template <>
struct CharCompare_<Blank_>
{
    static inline int _compare(const int c)
    {
        return isblank(c);
    }
}

template <>
struct CharCompare_<Alpha_>
{
    static inline int _compare(const int c)
    {
        return isalpha(c);
    }
}

template <>
struct CharCompare_<AlphaNum_>
{
    static inline int _compare(const int c)
    {
        return isalnum(c);
    }
}


template <>
struct CharCompare_<Digit_>
{
    static inline int _compare(const int c)
    {
        return isdigit(c);
    }
}

template <>
struct CharCompare_<Graph_>
{
    static inline int _compare(const int c)
    {
        return isgraph(c);
    }
}

template <>
struct CharCompare_<UnixEOL_>
{
    static inline int _compare(const int c)
    {
        return (c == '\n');
    }
}

template <>
struct CharCompare_<BackslashR_>
{
    static inline int _compare(const int c)
    {
        return (c == '\r');
    }
}

template <>
struct CharCompare_<Dna_>
{
    static inline int _compare(const int c)
    {
        switch (c)
        {
            case 'a':
            case 'c':
            case 'g':
            case 't':
            case 'A':
            case 'C':
            case 'G':
            case 'T': return true;
        }
        return false;
    }
}

template <>
struct CharCompare_<Dna5_>
{
    static inline int _compare(const int c)
    {
        switch (c)
        {
            case 'a':
            case 'c':
            case 'g':
            case 't':
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'n':
            case 'N':
                return true;
        }
        return false;
    }
}

// template <>
// struct CharCompare_<Char_, int c2>
// {
//     static inline int compare(const int c)
//     {
//         return c == c2;
//     }
// }

// ----------------------- Helper functions-----------------------------------

// read chars from record reader depending on condition
template <typename TCompare, // specialization of character comparison
          typename TRecordReader, // record reader
          typename TBuffer> // usually charstring, but maybe DnaString or so
inline int _readHelper(TRecordReader & reader,
                       TBuffer & buffer,
                       bool desiredOutcomeOfComparison = true) 
/*   desired behaviour of loop -> "readUntil()" or "readwhile()"  */
{
    clear(buffer);
    typedef Value<TRecordReader::_buffer>::Type TChar;
    for (TChar c = value(reader); hasMore(reader); goNext(reader))
    {
        if (resultCode(reader) != 0)
            return resultCode(reader);
        if (TCompare::_compare(c) == desiredOutcomeOfComparison)
            return 0;
        append(buffer, c, Generous()); // TODO Generous() is right?
    }
    returnIO_TOKENIZE_EOF_BEFORE_SUCCESS;
}

// same as above, just don't save characters read
template <typename TCompare, typename TRecordReader>
inline int _skipHelper(TRecordReader & reader,
                       bool desiredOutcomeOfComparison = true)
{
    typedef Value<RecordReader<TStream, TPass>::_buffer>::Type TChar;
    for (TChar c = value(reader); hasMore(reader); goNext(reader))
    {
        if (resultCode(reader) != 0)
            return resultCode(reader);
        if (TCompare::_compare(c) == desiredOutcomeOfComparison)
            return 0;
    }
    returnIO_TOKENIZE_EOF_BEFORE_SUCCESS;
}


// same as above but allows for a second character type to be ignored
template <typename TCompare, 
          typename TCompare2, // specialization of character class to be ignored
          typename TRecordReader, 
          typename TBuffer> 
inline int _readHelper(TRecordReader & reader,
                       TBuffer & buffer,
                       bool desiredOutcomeOfComparison = true) 

{
    clear(buffer);
    typedef Value<TRecordReader::_buffer>::Type TChar;
    for (TChar c = value(reader); hasMore(reader); goNext(reader))
    {
        if (resultCode(reader) != 0)
            return resultCode(reader);
        if (TCompare2::_compare(c))
            continue;
        if (TCompare::_compare(c) == desiredOutcomeOfComparison)
            return 0;
        append(buffer, c, Generous()); // TODO Generous() is right?
    }
    returnIO_TOKENIZE_EOF_BEFORE_SUCCESS;
}


template <typename TCompare, typename TRecordReader>
inline int _readAndCompareWithStr(TRecordReader & reader,
                                  const TString & str)
{
    bool win = false;
    for (int i = 0; i < length(str); ++i)
    {
        if (value(reader) != str[i]) //mismatch
            break;
        if (i == length(str)) // win
            win = true
        else if (!hasMore(reader))
            returnIO_TOKENIZE_EOF_BEFORE_SUCCESS;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return win ? 0 : IO_TOKONIZE_NO_SUCCESS;
}



// ----------------------- the real deal ------------------------------------

//TODO document
template <typename TStream, typename TPass>
inline int readUntilWhitespace(RecordReader<TStream, TPass> & reader,
                               TBuffer & buffer)
{
    return _readHelper< CharCompare_<Whitespace_>,
                        RecordReader<TStream, TPass>,
                        TBuffer >
                        (reader, buffer);
}

/* OLD FUNCTION
template<typename TFile, typename TChar>
inline String<char>
_parseReadWordUntilWhitespace(TFile& file, TChar& c)
{
    String<char> str(c);
    if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
        c = _streamGet(file);
        return str;
    }
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
        append(str, c);
    }
    return str;
}*/


//TODO document
template <typename TStream, typename TPass>
inline int readUntilBlank(RecordReader<TStream, TPass> & reader,
                          TBuffer & buffer)
{
    return _readHelper< CharCompare_<Blank_>,
                        RecordReader<TStream, TPass>,
                        TBuffer >
                        (reader, buffer);
}

//TODO document
template <typename TStream, typename TPass, typename TChar>
inline int readUntilX(RecordReader<TStream, TPass> & reader,
                      TBuffer & buffer,
                      TChar x)
{
    clear(buffer);
    typedef Value<TRecordReader::_buffer>::Type TChar;
    for (TChar c = value(reader); hasMore(reader); goNext(reader))
    {
        if (resultCode(reader) != 0)
            return resultCode(reader);
        if (c == x)
            return 0;
        append(buffer, c, Generous()); // TODO Generous() is right?
    }
    returnIO_TOKENIZE_EOF_BEFORE_SUCCESS;
}

//TODO document
template <typename TStream, typename TPass, typename TChar>
inline int readNChars(RecordReader<TStream, TPass> & reader,
                      TBuffer & buffer,
                      uint n)
{
    clear(buffer);
    resize(buffer, n);
    int i = 0;
    typedef Value<TRecordReader::_buffer>::Type TChar;
    for (int i = 0; i < n; ++i)
    {
        assignValue(buffer, i, value(reader));
        if (!hasMore(reader))
            returnIO_TOKENIZE_EOF_BEFORE_SUCCESS;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return 0;
}

// ---

//TODO document
template <typename TStream, typename TPass>
inline int skipUntilWhitespace(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper< CharCompare_<Whitespace_>,
                        RecordReader<TStream, TPass> >
                        (reader);
}

//TODO document
template <typename TStream, typename TPass>
inline int skipUntilBlank(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper< CharCompare_<Blank_>,
                        RecordReader<TStream, TPass> >
                        (reader);
}

//TODO document
template <typename TStream, typename TPass>
inline int skipUntilGraph(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper< CharCompare_<Graph_>,
                        RecordReader<TStream, TPass> >
                        (reader);
}
/* OLD FUNCTION, note that old function has wrong name -> see isgraph()
template<typename TFile, typename TChar>
inline void 
_parseSkipWhitespace(TFile& file, TChar& c)
{
    if ((unsigned) c > 32) return;
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if ((unsigned) c > 32) break;
    }
}*/


//TODO document
template <typename TStream, typename TPass>
inline int skipUntilX(RecordReader<TStream, TPass> & reader, const char x)
{
    typedef Value<RecordReader<TStream, TPass>::_buffer>::Type TChar;
    for (TChar c = value(reader); hasMore(reader); goNext(reader))
    {
        if (resultCode(reader) != 0)
            return resultCode(reader);
        if (c == x)
            return 0;
    }
    returnIO_TOKENIZE_EOF_BEFORE_SUCCESS;
}

/* OLD FUNCTION
template<typename TFile, typename TChar>
inline void 
_parseSkipUntilChar(TFile& file, const TChar &x, TChar& c)
{
    if (c == x) return;
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if (c == x) break;
    }
}

2nd OLD FUNCTION

template<typename TFile, typename TChar>
inline bool
_parseUntil(TFile & file, TChar& c, TChar x)
{
SEQAN_CHECKPOINT
    typename Position<TFile>::Type pos = _streamTellG(file);
    TChar c_before = c;
    while (!_streamEOF(file) && c != x){
        c = _streamGet(file);
    }
    if(!_streamEOF(file)) return true;
    _streamSeekG(file,pos);
    c = c_before;
    return false;
}

*/


//TODO document
// I think this function is stupid, but it is implemented in misc_parsing.h
// so I offer replacement; with proper buffering shiftOr would be faster
template <typename TStream, typename TPass, typename TString>
inline int skipUntilString(RecordReader<TStream, TPass> & reader,
                           const TString & str)
{
    if (length(str) < 1)
        return -1; //TODO some better error code
    while(skipUntilX(reader, str[0])==0)
    {
        switch(int r = _readAndCompareWithStr(reader, str))
        {
            case 0: return 0;
            case IO_TOKONIZE_NO_SUCCESS: break;
            default: return r;
        }
    }
    returnIO_TOKENIZE_EOF_BEFORE_SUCCESS;
    /* NOTE I am not 100% sure this will get every pattern with repitions in it.
     * it should be much better than the original, since it escapes on mismatch
     * directly and doesnt read (and possibly discard) N charecters in a row,
     * which definitely leads to misses. */
}
/* OLD FUNCTION
template<typename TFile, typename TChar, typename TSize>
inline bool
_parseUntil(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
//IOREV _nodoc_ _hasCRef_ _requiresSeek_ does "reset" if word not found
SEQAN_CHECKPOINT
    typename Position<TFile>::Type pos = _streamTellG(file);
    TChar c_before = c;
    while (!_streamEOF(file)){
        if(c == word[0])
            if(word == _parseReadWord(file,c,len))
                break;
        c = _streamGet(file);
    }
    if(!_streamEOF(file)) return true;
    _streamSeekG(file,pos);
    c = c_before;
    return false;
}*/


// ---

//TODO document
template <typename TStream, typename TPass, typename TBuffer>
inline int readLetters(RecordReader<TStream, TPass> & reader,
                      TBuffer & buffer)
{
    return _readHelper< CharCompare_<Alpha_>,
                        RecordReader<TStream, TPass>,
                        TBuffer >
                        (reader, buffer, false);
}
/* OLD FUNCTION, note that it has a wrong name (word is usually alphanum)
template<typename TFile, typename TChar>
inline String<char>
_parseReadWord(TFile & file, TChar& c)
{
    // Read word
    String<char> str(c);
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if (!_parseIsLetter(c)) break;
        append(str, c);
    }
    return str;
}*/


//TODO document
template <typename TStream, typename TPass, typename TBuffer>
inline int readAlphaNums(RecordReader<TStream, TPass> & reader,
                         TBuffer & buffer)
{
    return _readHelper< CharCompare_<AlphaNum_>,
                        RecordReader<TStream, TPass>,
                        TBuffer >
                        (reader, buffer, false);
}
/* OLD FUNCTION 
template<typename TFile, typename TString, typename TChar>
inline void
_parseReadIdentifier(TFile & file, TString& str, TChar& c)
{
    // Read identifier
    append(str, c, Generous());
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if (!_parseIsAlphanumericChar(c)) break;
        append(str, c, Generous());
    }
}*/


//TODO document
template <typename TStream, typename TPass>
inline int skipWhitespaces(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper< CharCompare_<Whitespace_>,
                        RecordReader<TStream, TPass> >
                        (reader, false);
}

//TODO document
template <typename TStream, typename TPass>
inline int skipBlanks(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper< CharCompare_<Blank_>,
                        RecordReader<TStream, TPass> >
                        (reader, false);
}

/* OLD FUNCTION, note that it has a wrong name
template<typename TFile, typename TChar>
inline void 
_parseSkipSpace(TFile& file, TChar& c)
{
    if (c != '\t' && c != ' ') return;
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if (c != '\t' && c != ' ') break;
    }
}

2nd OLD FUNCTION

template<typename TFile, typename TChar>
inline void 
_parseSkipBlanks(TFile& file, TChar& c)
{
    if ((c != ' ') && (c != '\t')) return;
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if ((c != ' ') && (c != '\t')) break;
    }
}

*/



//TODO document
template <typename TStream, typename TPass, typename TBuffer>
inline int readLine(RecordReader<TStream, TPass> & reader,
                          TBuffer & buffer)
{
    int r = _readHelper< CharCompare_<UnixEOL_>,
                         CharCompare_<BackslashR_>, // just ignore carriage return
                         RecordReader<TStream, TPass>,
                         TBuffer >
                         (reader, buffer);
    if (r != 0)
        return r;

//      replaced by ignoring BackslashR directly
//     // if file is in windows format we have \r in our buffer that has to go
//     if (value(buffer, length(buffer)-1) == '\r')
//         resize(buffer, length(buffer)-1); 
//         // TODO no realloc happens here, right?

    if (hasMore(reader))
        goNext(reader); // go to beginning of next line

    return resultCode(reader);
}

//TODO document
template <typename TStream, typename TPass, typename TBuffer>
inline int readLineStripTrailingBlanks(RecordReader<TStream, TPass> & reader,
                                       TBuffer & buffer)
{
    int r = readLine(reader, buffer);

    if (r != 0)
        return r;

    LENGTH<buffer> pos = length(buffer) -1;

    while (isblank(value(buffer, pos)) == 0) --pos;

    if (pos + 1 != length(buffer))
        resize(buffer, pos+1);
    return 0;
}
/* OLD FUNCTION
template<typename TFile, typename TChar>
inline String<char>
_parseReadFilepath(TFile& file, TChar& c)
{
//IOREV _nodoc_ _hasCRef_ _requiresSeek_ rename to something more generic _parseLineStripTrailingWhitespace; simplify code
    String<char> str(c);
    if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
        c = _streamGet(file);
        return str;
    }
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
        append(str, c);
    }
    typename Iterator<String<char>,Rooted >::Type str_it = end(str);    
    while(str_it != begin(str)) {
        --str_it;
        if(*str_it != ' ' && *str_it != '\t'){
        ++str_it;
        break;
        }
    }
    resize(str,position(str_it));
    return str;
}*/



//TODO document
template <typename TStream, typename TPass>
inline int skipLine(RecordReader<TStream, TPass> & reader)
{
    int r = _skipHelper< CharCompare_<UnixEOL_>,
                         RecordReader<TStream, TPass> >
                         (reader);
    if (r != 0)
        return r;

    if (hasMore(reader))
        goNext(reader); // go to beginning of next line

    return resultCode(reader);
}

/* OLD FUNCTION
template<typename TFile, typename TChar>
inline void 
_parseSkipLine(TFile& file, TChar& c)
{
//IOREV _nodoc_ _hasCRef_ treats \n as EOL 
    if (c == '\n') {
        c = _streamGet(file);
        return;
    }
    while (!_streamEOF(file)) {
        c = _streamGet(file);
        if (c == '\n') break;
    }
    c = _streamGet(file);
}*/



//TODO document
template <typename TStream, typename TPass, typename TBuffer>
inline int readDna5IgnoreWhitespaces(RecordReader<TStream, TPass> & reader,
                                     TBuffer & buffer)
{
    return _readHelper< CharCompare_<Dna5_>,
                        CharCompare_<Whitespace_>,
                        RecordReader<TStream, TPass>,
                        TBuffer >
                        (reader, buffer);
} 
// this would read a fasta or fastq sequence, since meta and qualities begin 
// with special chars


//TODO document
// also skips non-graph characters at linestart, c must not be with non-graph
// NOTE that position is on c, like in old function
// NOTE also, that we DO NOT reset in case of search failure as that would
// require seek
template <typename TStream, typename TPass>
inline int skipUntilLineBeginsWithC(RecordReader<TStream, TPass> & reader,
                                    const TChar & c )
)
{
    int r = 0;
    while ((r = skipLine(reader)) == 0 )
    {
        r = skipUntilGraph(reader);
        if (r != 0)
            return r;
        if (value(reader) == c)
            return 0;
    }
    return r;
}

/* OLD FUNCTION
/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with character x (skip whitespaces)
// zeigt am ende darauf!!!
template<typename TFile, typename TChar>
inline bool
_parseUntilBeginLine(TFile & file, TChar& c, TChar x)
{
//IOREV _nodoc_ _hasCRef_ _requiresSeek_ rename to something more understandable; does a "reset", if seek fails; unify behaviour and doc with similar functions
SEQAN_CHECKPOINT
    _parseSkipWhitespace(file,c);
    typename Position<TFile>::Type pos = _streamTellG(file);
    TChar c_before = c;
    while (!_streamEOF(file) && c != x){
        _parseSkipLine(file, c);
        _parseSkipWhitespace(file,c);
    }
    if(!_streamEOF(file)) return true;
    _streamSeekG(file,pos);
    c = c_before;
    return false;
}*/

//TODO document
// also skips non-graph characters at linestart, str must not begin with non-graph
// NOTE that position is behind str, like in old function, if a char behind
// string exists. NOTE also, that we DO NOT reset in case of search failure
// as that would require seek
template <typename TStream, typename TPass, typename TString>
inline int skipUntilLineBeginsWithStr(RecordReader<TStream, TPass> & reader,
                                      const TString & str )
)
{
    int r = 0;
    while ((r = skipLine(reader)) == 0 )
    {
        r = skipUntilGraph(reader);
        if (r != 0)
            return r;
        switch(r = _readAndCompareWithStr(reader, str))
        {
            case 0: return 0;
            case IO_TOKONIZE_NO_SUCCESS: break;
            default: return r;
        }
    }
    return r;
}
/* OLD FUNCTION 
/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with word
//zeigt am ende dahinter!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parseUntilBeginLine(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
//IOREV _nodoc_ _hasCRef_ _requiresSeek_ rename to something more understandable; does a "reset", if seek fails; unify behaviour and doc with similar functions
SEQAN_CHECKPOINT
    _parseSkipWhitespace(file,c);
    typename Position<TFile>::Type pos = _streamTellG(file);
    TChar c_before = c;
    while (!_streamEOF(file)){
        if(c == word[0])
            if(word == _parseReadWord(file,c,len))
                break;
        _parseSkipLine(file, c);
        _parseSkipWhitespace(file,c);
    }
    if(!_streamEOF(file)) return true;
    _streamSeekG(file,pos);
    c = c_before;
    return false;
}*/


// TODO missing functions from misc_parsing.h:
// _parseUntilBeginLine with max num lines
//_parseUntilBeginLineOneOf(
// all functions that require lexical_cast<>() which isnt there yet

// TODO stuff from other files (shouldnt be that much and much is duplicate)


}

#endif // def SEQAN_STREAM_RECORD_READER_SINGLE_H_
