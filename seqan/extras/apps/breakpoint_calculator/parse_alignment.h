// ==========================================================================
//                           Breakpoint Calculator
// ==========================================================================
// Copyright (C) 2012 by Birte Kehr
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
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_
#define SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_

#include <seqan/stream.h>
#include <seqan/align.h>

#include "parse_alignment.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tags Xmfa, XmfaSwap, Maf
// ----------------------------------------------------------------------------

struct Xmfa_;
typedef Tag<Xmfa_> Xmfa;

struct XmfaSwap_;
typedef Tag<XmfaSwap_> XmfaSwap;

struct Maf_;
typedef Tag<Maf_> Maf;

// ----------------------------------------------------------------------------
// Class AlignmentBlockRow
// ----------------------------------------------------------------------------

template <typename TPosition, typename TSize>
struct AlignmentBlockRow
{
    TSize rowNum;             // row number
    CharString chromosomeId;  // chromosome identifier

    TPosition startPos;      // in source seq
    TPosition endPos;        // in source seq

    bool orientation;
    // startPos is always smaller or equal endPos also if orientation is false

    AlignmentBlockRow() {}

    AlignmentBlockRow(TSize id, CharString chr, TPosition start, TPosition end, bool ori)
    {
        rowNum = id;
        chromosomeId = chr;
        startPos = start;
        endPos = end;
        orientation = ori;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function skipHeader()
// ----------------------------------------------------------------------------

template <typename TStream, typename TPassSpec, typename TTag>
int
skipHeader(RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           TTag const &)
{
    CharString buffer;
    int res = 0;

    // Skip comment lines
    while (value(recordReader) == '#')
    {
        res = skipLine(recordReader);
        if (res)
            return res;

        if (atEnd(recordReader))
            return 1;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                                 [Xmfa]
// ----------------------------------------------------------------------------

// Reads one collinear alignment block from a file in XMFA format

template <typename TSequence, typename TSize, typename TStream, typename TPassSpec>
int
readRecord(Align<TSequence, ArrayGaps> & align,
           std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > & idToRowsMap,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           bool swapPos,
           Xmfa const &)
{
    typedef typename Position<TSequence>::Type TPosition;
    typedef AlignmentBlockRow<TSize, TSize> TRow;

    CharString buffer;
    int res = 0;

    TPosition startPos, endPos;
    CharString id, chromosomeId = "";
    bool orientation;

    while (!skipChar(recordReader, '>'))
    {
        typename Row<Align<TSequence, ArrayGaps> >::Type gapseq;

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read seq id
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        id = buffer;

        res = skipChar(recordReader, ':');
        if (res)
            return res;

        // read start position
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        startPos = lexicalCast<TPosition>(buffer);

        res = skipChar(recordReader, '-');
        if (res)
            return res;

        // read end position
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        endPos = lexicalCast<TPosition>(buffer);

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read orientation
        clear(buffer);
        res = readNChars(buffer, recordReader, 1);
        if (res)
            return res;

        if (buffer  == "+")
            orientation = true;
        else if (buffer == "-")
            orientation = false;
        else
            return 1;

        if (swapPos && orientation == false)
        {
            TPosition help = startPos;
            startPos = endPos;
            endPos = help;
        }
        if (startPos != 0)
            --startPos;

        // skip rest of line
        res = skipLine(recordReader);
        if (res)
            return res;

        // skip gapped sequence
        clear(buffer);
        while (value(recordReader) != '>' && value(recordReader) != '=')
        {
            res = skipLine(recordReader);
            if (res)
                return res;
        }

        if (endPos != startPos)
        {
            appendValue(idToRowsMap[id], TRow(length(rows(align)), chromosomeId, startPos, endPos, orientation));
            appendValue(rows(align), gapseq);
        }
    }

    if (value(recordReader) != '=')
        return 1;

    skipUntilLineBeginsWithChar(recordReader, '>');

    return 0;
}

template <typename TSequence, typename TSize, typename TStream, typename TPassSpec>
int
readRecord(Align<TSequence, ArrayGaps> & align,
           std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > & idToRowsMap,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           Xmfa const &)
{
    return readRecord(align, idToRowsMap, recordReader, false, Xmfa());
}

template <typename TSequence, typename TSize, typename TStream, typename TPassSpec>
int
readRecord(Align<TSequence, ArrayGaps> & align,
           std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > & idToRowsMap,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           XmfaSwap const &)
{
    return readRecord(align, idToRowsMap, recordReader, true, Xmfa());
}

// ----------------------------------------------------------------------------
// Function readRecord()                                                  [Maf]
// ----------------------------------------------------------------------------

// Reads one collinear alignment block from a file in MAF format

template <typename TSequence, typename TSize, typename TStream, typename TPassSpec>
int
readRecord(Align<TSequence, ArrayGaps> & align,
           std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > & idToRowsMap,
           RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
           Maf const &)
{
    typedef typename Position<TSequence>::Type TPosition;
    typedef AlignmentBlockRow<TSize, TSize> TRow;

    CharString buffer;
    int res = 0;

    TPosition startPos, len, seqLen;
    CharString seqId, chromosomeId = "";
    bool orientation;

    skipWhitespaces(recordReader);
    res = skipChar(recordReader, 'a');
    if (res)
        return res;

    res = skipLine(recordReader);
    if (res)
        return res;

    skipWhitespaces(recordReader);
    while (!skipChar(recordReader, 's'))
    {
        typename Row<Align<TSequence, ArrayGaps> >::Type gapseq;

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read seq id
        clear(buffer);
        res = readUntilOneOf(buffer, recordReader, '.', ' ', '\t');
        if (res)
            return res;

        seqId = buffer;

        res = skipChar(recordReader, '.');
        if (!res)
        {
            // read chromosome id
            clear(buffer);
            res = readUntilWhitespace(buffer, recordReader);
            if (res)
                return res;

            chromosomeId = buffer;
        }

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read start position
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        startPos = lexicalCast<TPosition>(buffer);

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // read length
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        len = lexicalCast<TPosition>(buffer);

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // skip orientation
        clear(buffer);
        res = readNChars(buffer, recordReader, 1);
        if (res)
            return res;

        if (buffer  == '+')
            orientation = true;
        else if (buffer == '-')
            orientation = false;
        else
            return 1;

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // skip seq length
        clear(buffer);
        res = readDigits(buffer, recordReader);
        if (res)
            return res;

        seqLen = lexicalCast<TPosition>(buffer);

        res = skipWhitespaces(recordReader);
        if (res)
            return res;

        // skip gapped seq
        skipLine(recordReader);

        if (!orientation)
            startPos = seqLen - (startPos + len);

        if (len > 0)
        {
            appendValue(idToRowsMap[seqId],
                        TRow(length(rows(align)), chromosomeId, startPos, startPos + len, orientation));
            appendValue(rows(align), gapseq);
        }

        skipWhitespaces(recordReader);
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function parseAlignment()
// ----------------------------------------------------------------------------

template <typename TSequence, typename TStringSet, typename TFile, typename TTag>
int
parseAlignment(String<Align<TSequence, ArrayGaps> > & aligns,
               String<std::map<CharString, String<AlignmentBlockRow<typename Size<TStringSet>::Type, typename Size<TStringSet>::Type> > > > & idToRowsMaps,
               TStringSet & /*seqs*/,
               TFile & file,
               bool verbose,
               TTag const tag)
{
    typedef typename Size<TStringSet>::Type TSize;
    typedef AlignmentBlockRow<TSize, TSize> TRow;
    typedef std::map<CharString, String<TRow> > TMap;

    typedef Align<TSequence, ArrayGaps> TAlign;

    RecordReader<TFile, SinglePass<> > recordReader(file);
    int res = 0;

    res = skipHeader(recordReader, tag);
    if (res)
        return res;

    while (!atEnd(recordReader))
    {
        TAlign align;
        TMap idToRows;
        res = readRecord(align, idToRows, recordReader, tag);
        if (res)
            return res;

        appendValue(aligns, align);
        appendValue(idToRowsMaps, idToRows);

        skipWhitespaces(recordReader);
        while (!atEnd(recordReader) && value(recordReader) == '#')
            res = skipLine(recordReader);
    }
    if (!atEnd(recordReader))
        return 1;

    if (verbose)
    {
        std::cout << length(aligns) << " alignment block";
        if (length(aligns) != 1)
            std::cout << "s";
        std::cout << " loaded." << std::endl;
    }

    return 0;
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_
