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



namespace seqan {

// ----------------------- enum ------------------------------------

/**
.Enum.Tokenize::Result
..cat:Input / Output
..summary:Enum with return values for Tokenizing operations.
..value.SUCCESS:Reading the specified data succeeded.
..value.EOF_BEFORE_SUCCESS:End of file was reached before the pattern was found.
..value.NO_SUCCESS:The pattern was not found.
..see:Function.mmapAdvise
..include:seqan/stream.h
 */
struct Tokenize
{
    enum Result {
        SUCCESS = 0,
        EOF_BEFORE_SUCCESS = 1024,
        NO_SUCCESS = 1025
    };
};



    
// ----------------------- Helper structs ------------------------------------
struct Whitespace__;
struct Blank__;
struct Char__;
struct Digit__;
struct Alpha__;
struct AlphaNum__;
struct UnixEOL__;
struct BackslashR__;
struct Graph__;

typedef Tag<Whitespace__> Whitespace_;
typedef Tag<Blank__> Blank_;
typedef Tag<Char__> Char_;
typedef Tag<Digit__> Digit_;
typedef Tag<Alpha__> Alpha_;
typedef Tag<AlphaNum__> AlphaNum_;
typedef Tag<UnixEOL__> UnixEOL_;
typedef Tag<BackslashR__> BackslashR_;
typedef Tag<Graph__> Graph_;


// template <typename TSpec>
// struct CharCompare_;


inline int
_charCompare(int const c, Whitespace_ const & /* tag*/)
{
    return isspace(c);
}


inline int
_charCompare(int const c, Blank_ const & /* tag*/)
{
    return isblank(c);
}

inline int
_charCompare(int const c, Alpha_ const & /* tag*/)
{
    return isalpha(c);
}

inline int
_charCompare(int const c, AlphaNum_ const & /* tag*/)
{
    return isalnum(c);
}

inline int
_charCompare(int const c, Digit_ const & /* tag*/)
{
    return isdigit(c);
}

inline int
_charCompare(int const c, Graph_ const & /* tag*/)
{
    return isgraph(c);
}

inline int
_charCompare(int const c, UnixEOL_ const & /* tag*/)
{
    return (c == '\n');
}

inline int
_charCompare(int const c, BackslashR_ const & /* tag*/)
{
    return (c == '\r');
}

inline int
_charCompare(int const c, Dna const & /* tag*/)
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

inline int
_charCompare(int const c, Dna5 const & /* tag*/)
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

// template <>
// struct CharCompare_<Char_, int c2>
// {
//     static inline int compare(int const c)
//     {
//         return c == c2;
//     }
// }

// ----------------------- Helper functions-----------------------------------

// read chars from record reader depending on condition
template <typename TTagSpec, // specialization of character comparison
          typename TRecordReader, // record reader
          typename TBuffer> // usually charstring, but maybe DnaString or so
inline int
_readHelper(TBuffer & buffer,
            TRecordReader & reader,
            Tag<TTagSpec> const & tag,
            bool const desiredOutcomeOfComparison) 
/*   desired behaviour of loop -> "readUntil()" or "readwhile()"  */
{
    clear(buffer);
    typedef Value<TRecordReader::_buffer>::Type TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (_charCompare(c, tag) == desiredOutcomeOfComparison)
            return 0;
        append(buffer, c, Generous()); // TODO Generous() is right?
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return Tokenize::EOF_BEFORE_SUCCESS;
}

template <typename TSpec, // specialization of character comparison
          typename TRecordReader, // record reader
          typename TBuffer> // usually charstring, but maybe DnaString or so
inline int
_readHelper(TBuffer & buffer,
            TRecordReader & reader,
            Tag<TSpec> const & tag)
/*   default behaviour of loop is "readUntil()" */
{
    _readHelper(buffer, reader, tag, true);
}

// same as above, just don't save characters read
template <typename TTagSpec, typename TRecordReader>
inline int
_skipHelper(TRecordReader & reader,
            Tag<TSpec> const & tag,
            bool const desiredOutcomeOfComparison)
{
    typedef Value<TRecordReader::_buffer>::Type TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (_charCompare(c, tag) == desiredOutcomeOfComparison)
            return 0;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return Tokenize::EOF_BEFORE_SUCCESS;
}

template <typename TTagSpec, typename TRecordReader>
inline int
_skipHelper(TRecordReader & reader,
            Tag<TSpec> const & tag)
{
    _skipHelper(reader, tag, true)
}

// same as above but allows for a second character type to be ignored
template <typename TTagSpec, 
          typename TTagSpec2, // specialization of character class to be ignored
          typename TRecordReader, 
          typename TBuffer> 
inline int
_readHelper(TBuffer & buffer,
            TRecordReader & reader,
            Tag<TTagSpec> const & compTag,
            Tag<TTagSpec2> const & skipTag,
            bool const desiredOutcomeOfComparison)
{
    clear(buffer);
    typedef Value<TRecordReader::_buffer>::Type TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (!_charCompare(c, skipTag))
        {
            if (_charCompare(c, compTag) == desiredOutcomeOfComparison)
                return 0;
            append(buffer, c, Generous()); // TODO Generous() is right?
        }
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return Tokenize::EOF_BEFORE_SUCCESS;
}


template <typename TTagSpec,
          typename TTagSpec2, // specialization of character class to be ignored
          typename TRecordReader,
          typename TBuffer>
inline int
_readHelper(TBuffer & buffer,
            TRecordReader & reader,
            Tag<TTagSpec> const & compTag,
            Tag<TTagSpec> const & skipTag)

{
    _readHelper(buffer, reader, compTag, skipTag, true);
}

template <typename TCompare, typename TRecordReader, typename TString>
inline int
_readAndCompareWithStr(TRecordReader & reader,
                       TString const & str)
{
    bool win = false;
    for (int i = 0; i < length(str); ++i)
    {
        if (atEnd(reader))
            return Tokenize::EOF_BEFORE_SUCCESS;
        if (value(reader) != str[i]) //mismatch
            break;
        if (i == length(str)) // win
            win = true;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return win ? 0 : Tokenize::NO_SUCCESS;
}



// ----------------------- the real deal ------------------------------------

/**
.Function.readUntilWhitespace
..cat:Input/Output
..summary:Read characters from stream into buffer until Whitespace is encountered
..signature:readUntilWhitespace(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString, Shortcut.DnaString or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:Whitespace is more than '' and '\t', see @Function.isspace@
..remarks:This function stops *on* the whitespace character. The whitespace is not written to buffer.
..include:seqan/stream.h
..see:Function.isspace
..see:Function.skipUntilWhitespace
..see:Enum.Tokenize::Result
 */

template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilWhitespace(TBuffer & buffer,
                    RecordReader<TStream, TPass> & reader)
{
    return _readHelper(buffer,
                       reader,
                       Whitespace_());
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


/**
.Function.readUntilBlank
..cat:Input/Output
..summary:Read characters from stream into buffer until Blank is encountered
..signature:readUntilBlank(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString, Shortcut.DnaString or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:Blank is '' and '\t', see @Function.isblank@
..remarks:This function stops *on* the blank character. The blank is not written to buffer.
..include:seqan/stream.h
..see:Function.isblank
..see:Function.skipUntilBlank
..see:Enum.Tokenize::Result
 */
template <typename TBuffer, typename TStream, typename TPass>
inline int
readUntilBlank(TBuffer & buffer.
               RecordReader<TStream, TPass> & reader)
{
    return _readHelper(buffer,
                       reader,
                       Blank_());
}

/**
.Function.readUntilChar
..cat:Input/Output
..summary:Read characters from stream into buffer until Char is encountered
..signature:readUntilChar(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader, TChar const & x)
..param.buffer:The buffer to write to
...type:Shortcut.CharString, Shortcut.DnaString or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.x:The character to stop on
...type:$char$ or similar
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops *on* the character x. It is not written to buffer.
..include:seqan/stream.h
..see:Enum.Tokenize::Result
..see:skipUntilChar
 */
template <typename TBuffer, typename TStream, typename TPass, typename TChar>
inline int
readUntilChar(TBuffer & buffer,
              RecordReader<TStream, TPass> & reader,
              TChar const & x)
{
    typedef Value<RecordReader<TStream, TPass>::_buffer>::Type TChar;

    clear(buffer);

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (c == x)
            return 0;
        append(buffer, c, Generous()); // TODO Generous() is right?
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return Tokenize::EOF_BEFORE_SUCCESS;
}

/**
.Function.readNChars
..cat:Input/Output
..summary:Read exactly n characters from stream into buffer
..signature:readNChars(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader, uint const n)
..param.buffer:The buffer to write to
...type:Shortcut.CharString, Shortcut.DnaString or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.n:The number of characters to read
...type:$uint$
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..include:seqan/stream.h
..see:Function.isblank
..see:Enum.Tokenize::Result
 */
template <typename TBuffer, typename TStream, typename TPass, typename TChar>
inline int
readNChars(TBuffer & buffer,
           RecordReader<TStream, TPass> & reader,
           uint const n)
{
    typedef Value<RecordReader<TStream, TPass>::_buffer>::Type TChar;

    clear(buffer);
    resize(buffer, n);

    for (int i = 0; i < n; ++i)
    {
        if (atEnd(reader))
            return Tokenize::EOF_BEFORE_SUCCESS;
        assignValue(buffer, i, value(reader));

        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return 0;
}

// ---
/**
.Function.skipUntilWhitespace
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until Whitespace is encountered
..signature:skipUntilWhitespace(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error skiping
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:Whitespace is more than '' and '\t', see @Function.isspace@
..remarks:This function stops *on* the whitespace character. The whitespace is not skipped.
..include:seqan/stream.h
..see:Function.isspace
..see:Function.skipWhitespaces
..see:Function.readUntilWhitespace
..see:Enum.Tokenize::Result
 */
template <typename TStream, typename TPass>
inline int
skipUntilWhitespace(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper(reader, Whitespace_());
}

/**
.Function.skipUntilBlank
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until Blank is encountered
..signature:skipUntilBlank(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error skiping
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops *on* the blank character. The blank is not skipped.
..include:seqan/stream.h
..see:Function.isblank
..see:Enum.Tokenize::Result
..see:Function.readUntilBlank
 */
template <typename TStream, typename TPass>
inline int
skipUntilBlank(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper(reader, Blank_());
}
/* OLD FUNCTION
template <typename TIterator>
inline bool
_seekWhiteSpace(TIterator &it, TIterator itEnd)
{
    while (!_isWhiteSpace(*it))
        if (++it == itEnd) return false;
    return true;
}*/


/**
.Function.skipUntilGraph
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until printable, non-' ' character is encountered
..signature:skipUntilGraph(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error skiping
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops *on* the graph character. The graph is not skipped.
..remarks:See @Function.isgraph@ for details on the "graph"-group of characters
..include:seqan/stream.h
..see:Function.isgraph
..see:Enum.Tokenize::Result
 */
template <typename TStream, typename TPass>
inline int
skipUntilGraph(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper(reader, Graph_());
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


/**
.Function.skipUntilChar
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until Char is encountered
..signature:skipUntilChar(RecordReader<TStream, TPass> & recordReader, TChar const & x)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.x:The character to stop on
...type:$char$ or similar
..returns:$int$, 0 if there was no error skiping
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops *on* the character x. x is not skipped.
..include:seqan/stream.h
..see:Enum.Tokenize::Result
..see:Function.readUntilChar
 */
template <typename TStream, typename TPass, typename TChar>
inline int
skipUntilChar(RecordReader<TStream, TPass> & reader,
              TChar const & x)
{
    typedef Value<RecordReader<TStream, TPass>::_buffer>::Type TChar;

    while (!atEnd(reader))
    {
        TChar c = value(reader);
        if (c == x)
            return 0;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return Tokenize::EOF_BEFORE_SUCCESS;
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

can also replace things like:
template <typename TIterator>
inline bool
_seekTab(TIterator& it, TIterator itEnd)
{
    for (; it != itEnd; ++it)
        if (*it == '\t') return true;
    return false;
}

*/


/**
.Function.skipUntilString
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until String is encountered
..signature:skipUntilString(RecordReader<TStream, TPass> & recordReader, TString const & str)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..param.str:The string to stop on
...type:$char$ or similar
..returns:$int$, 0 if there was no error skiping
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops *behind* the character string
..include:seqan/stream.h
..see:Enum.Tokenize::Result
..see:Function.skipUntilChar
 */
// I think this function is stupid, but it is implemented in misc_parsing.h
// so I offer replacement; with proper buffering shiftOr would be faster
template <typename TStream, typename TPass, typename TString>
inline int
skipUntilString(RecordReader<TStream, TPass> & reader,
                           TString const & str)
{
    if (length(str) < 1)
        return -1; //TODO some better error code
    while(skipUntilChar(reader, str[0])==0)
    {
        switch(int r = _readAndCompareWithStr(reader, str))
        {
            case 0: return 0;
            case Tokenize::NO_SUCCESS: break;
            default: return r;
        }
    }
    return Tokenize::EOF_BEFORE_SUCCESS;
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

/**
.Function.readLetters
..cat:Input/Output
..summary:Read characters from stream as long as characters are letters
..signature:skipLetters(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString, Shortcut.DnaString or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops *behind* the last letter read.
..include:seqan/stream.h
..see:Enum.Tokenize::Result
..see:Function.isalpha
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readLetters(TBuffer & buffer,
                       RecordReader<TStream, TPass> & reader)
{
    return _readHelper(buffer,
                       reader,
                       Alpha_(),
                       false);
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


/**
.Function.readAlphaNums
..cat:Input/Output
..summary:Read characters from stream as long as characters are letters
..signature:skipAlphaNums(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops *behind* the last letter read.
..include:seqan/stream.h
..see:Enum.Tokenize::Result
..see:Function.isalnum
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readAlphaNums(TBuffer & buffer,
                         RecordReader<TStream, TPass> & reader)
{
    return _readHelper(buffer,
                       reader,
                       AlphaNum_(),
                       false);
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


/**
.Function.skipWhitespaces
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until non-Whitespace is encountered
..signature:skipWhitespaces(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error skiping
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:Whitespace is more than '' and '\t', see @Function.isspace@
..remarks:This function stops *behind* the last whitespace character.
..include:seqan/stream.h
..see:Function.isspace
..see:Function.readUntilWhitespace
..see:Function.skipUntilWhitespace
..see:Enum.Tokenize::Result
 */
template <typename TStream, typename TPass>
inline int
skipWhitespaces(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper(reader, Whitespace_(), false);
}

/**
.Function.skipBlanks
..cat:Input/Output
..summary:Skip (i.e. read without saving) characters from stream until non-Blank is encountered
..signature:skipBlanks(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error skiping
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:Blank is more than '' and '\t', see @Function.isspace@
..remarks:This function stops *behind* the last whitespace character.
..include:seqan/stream.h
..see:Function.isspace
..see:Function.readUntilBlank
..see:Function.skipUntilBlank
..see:Enum.Tokenize::Result
 */
template <typename TStream, typename TPass>
inline int
skipBlanks(RecordReader<TStream, TPass> & reader)
{
    return _skipHelper(reader, Blank_(), false);
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


/**
.Function.readLine
..cat:Input/Output
..summary:Read a line from stream and save it to buffer
..signature:readLine(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops on the beginning of the next line, if there is a next line
..remarks:End-line characters are not written to buffer.
..remarks:Works on ANSI EOL and on Unix EOL.
..include:seqan/stream.h
..see:Enum.Tokenize::Result
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readLine(TBuffer & buffer, RecordReader<TStream, TPass> & reader)
{
    int r = _readHelper(buffer,
                        reader,
                        UnixEOL_(), // abort on Newline
                        BackslashR_()); // just ignore carriage return
    if (r != 0)
        return r;

    if (!atEnd(reader))
        goNext(reader); // go to beginning of next line

    return resultCode(reader);
}

/**
.Function.readLineStripTrailingBlanks
..cat:Input/Output
..summary:Read a line from stream and save it to buffer, remove trailing blanks
..signature:readLineStripTrailingBlanks(TBuffer & buffer, RecordReader<TStream, TPass> & recordReader)
..param.buffer:The buffer to write to
...type:Shortcut.CharString or similar
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops on the beginning of the next line, if there is a next line
..remarks:End-line characters and all trailing blanks are not written to buffer.
..remarks:Works on ANSI EOL and on Unix EOL.
..include:seqan/stream.h
..see:Enum.Tokenize::Result
..see:Function.isblank
 */
template <typename TStream, typename TPass, typename TBuffer>
inline int
readLineStripTrailingBlanks(TBuffer & buffer,
                            RecordReader<TStream, TPass> & reader)
{
    int r = readLine(reader, buffer);

    if (r != 0)
        return r;   // 1234567890

    LENGTH<buffer> pos = length(buffer) -1;

    if (pos < 0)
        return 0;

    while ((isblank(value(buffer, pos)) == 0) && (pos >= 0)) --pos;

    if (pos + 1 != length(buffer))
        resize(buffer, pos+1, Exact());
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


/**
.Function.skipLine
..cat:Input/Output
..summary:Skip a line in stream and go to beginning of next
..signature:readLine(RecordReader<TStream, TPass> & recordReader)
..param.recordReader:The @Class.RecordReader@ to read from.
...type:Class.RecordReader
..returns:$int$, 0 if there was no error reading
..returns:non-zero value on errors, especially Tokenize::EOF_BEFORE_SUCCESS
..remarks:This function stops on the beginning of the next line, if there is a next line
..remarks:Works on ANSI EOL and on Unix EOL.
..include:seqan/stream.h
..see:Enum.Tokenize::Result
..see:Function.readLine
 */
template <typename TStream, typename TPass>
inline int
skipLine(RecordReader<TStream, TPass> & reader)
{
    int r = _skipHelper(reader, UnixEOL_());

    if (r != 0)
        return r;

    if (!atEnd(reader))
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
inline int
readDna5IgnoringWhitespaces(TBuffer & buffer,
                          RecordReader<TStream, TPass> & reader)
{
    return _readHelper(buffer, reader, Dna5_(), Whitespace_(), false);
} 
// this would read a fasta or fastq sequence, since meta and qualities begin 
// with special chars


//TODO document
// also skips non-graph characters at linestart, c must not be non-graph
// NOTE that position is on c, like in old function
// NOTE also, that we DO NOT reset in case of search failure as that would
// require seek
template <typename TStream, typename TPass, typename TChar>
inline int
skipUntilLineBeginsWithChar(RecordReader<TStream, TPass> & reader,
                            TChar const & c )
{
    if (!isgraph(c))
        return 
    int r = 0;
    while (!atEnd(reader) && (r = skipLine(reader)) == 0 )
    {
        r = skipUntilGraph(reader);
        if (r != 0)
            return r;
        if (value(reader) == c)
            return 0;
    }
    if (atEnd(reader))
        return Tokenize::EOF_BEFORE_SUCCESS;
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
inline int
skipUntilLineBeginsWithStr(RecordReader<TStream, TPass> & reader,
                           const TString & str )
{
    int r = 0;
    while (!atEnd(reader) && (r = skipLine(reader)) == 0 )
    {
        r = skipUntilGraph(reader);
        if (r != 0)
            return r;
        switch(r = _readAndCompareWithStr(reader, str))
        {
            case 0: return 0;
            case Tokenize::NO_SUCCESS: break;
            default: return r;
        }
    }
    if (atEnd(reader))
        return Tokenize::EOF_BEFORE_SUCCESS;
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

//TODO document
// also skips non-graph characters at linestart, c must not be non-graph
// NOTE that position is on c, like in old function
// NOTE also, that we DO NOT reset in case of search failure as that would
// require seek
template <typename TStream, typename TPass, typename TString>
inline int
skipUntilLineBeginsWithOneCharOfStr(RecordReader<TStream, TPass> & reader,
                                    TString const & str )
{
    if (!isgraph(c))
        return
    int r = 0;
    while (!atEnd(reader) && (r = skipLine(reader)) == 0 )
    {
        r = skipUntilGraph(reader);
        if (r != 0)
            return r;
        for (int i = 0; i < length(str); ++i)
        {
            if (value(reader) == value(str,i))
                return 0;
        }

    }
    if (atEnd(reader))
        return Tokenize::EOF_BEFORE_SUCCESS;
    return r;
}
/* OLD FUNCTION:
/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with one of the characters in string x (skip whitespaces)
//zeigt am ende darauf!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parseUntilBeginLineOneOf(TFile & file, TChar& c, String<TChar> & x, TSize len)
{
//IOREV _nodoc_ _hasCRef_ _requiresSeek_ rename to something more understandable; does a "reset", if seek fails; unify behaviour and doc with similar functions
SEQAN_CHECKPOINT
    _parseSkipWhitespace(file,c);
    typename Position<TFile>::Type pos = _streamTellG(file);
    TChar c_before = c;
    bool found = false;
    while (!_streamEOF(file)){
        for(int i = 0; i < len; ++i)
            if(c == x[i])
            {
                found = true;
                break;
            }
        if(found) break;
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

#endif // def SEQAN_STREAM_RECORD_READER_SINGLE_H_
