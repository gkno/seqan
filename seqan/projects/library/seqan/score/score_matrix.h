/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
 ============================================================================
  Code for score matrices with data from files or built-in data.
 ==========================================================================*/

#ifndef SEQAN_SCORE_SCORE_MATRIX_H_
#define SEQAN_SCORE_SCORE_MATRIX_H_

#include <seqan/file.h>

// TODO(holtgrew): If the complex type conversions are necessary, a static_cast<> is more C++ and explicit.


namespace SEQAN_NAMESPACE_MAIN {

template <typename TValue, typename TSequenceValue, typename TSpec>
struct _ScoringMatrixData;


template <typename TSequenceValue = AminoAcid, typename TSpec = Default>
struct ScoreMatrix;


/**
.Tag.File Format.tag.ScoreMatrixFile:Score matrix file.
..include:seqan/score.h
*/
struct TagScoreMatrixFile_;
typedef Tag<TagScoreMatrixFile_> const ScoreMatrixFile;


/**
.Spec.Score Matrix:
..cat:Scoring
..summary:A general scoring matrix.
..general:Class.Score
..signature:Score<TValue, ScoreMatrix<TSequenceValue, TSpec> >
..param.TValue:Type of the score values.
...default:$int$
..param.TSequenceValue:Type of alphabet underlying the matrix.
...default:$AminoAcid$
..include:seqan/score.h
 */
template <typename TValue, typename TSequenceValue, typename TSpec>
class Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > {
public:
    // Static computation of the required array size.
    enum {
        VALUE_SIZE = ValueSize<TSequenceValue>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    // The data table.
    TValue data_tab[TAB_SIZE];

    // The gap extension score.
    TValue data_gap_extend;

    // The gap open score.
    TValue data_gap_open;

    /**
.Memfunc.Score Matrix#Score
..cat:Scoring
..summary:Constructor.
..class:Spec.Score Matrix
..signature:Score(gapExtend)
..param.gapExtend:The gap extension penalty.
...remark:TValue
     */
    explicit Score(TValue _gap_extend = -1)
        : data_gap_extend(_gap_extend),
          data_gap_open(_gap_extend) {
        SEQAN_CHECKPOINT;
        setDefaultScoreMatrix(*this, TSpec());
    }

    /**
.Memfunc.Score Matrix#Score
..signature:Score(gapExtend, gapOpen)
..param.gapOpen:The gap open penalty.
...remark:TValue
     */
    Score(TValue _gap_extend, TValue _gap_open)
        : data_gap_extend(_gap_extend), data_gap_open(_gap_open) {
        SEQAN_CHECKPOINT;
        setDefaultScoreMatrix(*this, TSpec());
    }

    /**
.Memfunc.Score Matrix#Score
..signature:Score(filename, gapExtend)
..param.filename:The path to the file to load.
...type:Class.String
..see:Function.loadScoreMatrix
     */
    template <typename TString>
    Score(TString const & filename, TValue _gap_extend = -1)
        : data_gap_extend(_gap_extend), data_gap_open(_gap_extend) {
        SEQAN_CHECKPOINT;
        loadScoreMatrix(*this, filename);
    }

    /**
.Memfunc.Score Matrix#Score
..signature:Score(filename, gapExtend, gapOpen)
     */
    template <typename TString>
    Score(TString const & filename, TValue _gap_extend, TValue _gap_open)
        : data_gap_extend(_gap_extend), data_gap_open(_gap_open) {
        SEQAN_CHECKPOINT;
        loadScoreMatrix(*this, filename);
    }
};


// TODO(holtgrew): Does it make sense to document each Score specialization?  Should dddoc show a list of all specializations of a class?
template <typename TValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2>
inline TValue
score(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc, TVal1 val1, TVal2 val2) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    // TODO(holtgrew): Why not implicit cast?
    unsigned int i = (TSequenceValue) val1;  // conversion TVal1 => TSequenceValue => integral
    unsigned int j = (TSequenceValue) val2;  // conversion TVal2 => TSequenceValue => integral
    return sc.data_tab[i * TScore::VALUE_SIZE + j];
}


/**
.Function.setScore:
..cat:Scoring
..summary:Set the substitution score between two values.
..signature:setScore(scoreMatrix, val1, val2, score)
..param.scoreMatrix:
...type:Spec.Score Matrix
..param.val1:First value.
..param.val2:Second value.
..param.score:The value to set the score to.
..include:seqan/score.h
 */
template <typename TValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2, typename T>
inline void
setScore(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TVal1 val1, TVal2 val2, T score) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    // TODO(holtgrew): Why not implicit cast?
    unsigned int i = (TSequenceValue) val1;  // conversion TVal1 => TSequenceValue => integral
    unsigned int j = (TSequenceValue) val2;  // conversion TVal2 => TSequenceValue => integral
    sc.data_tab[i * TScore::VALUE_SIZE + j] = score;
}


/**
.Function.setDefaultScoreMatrix:
..cat:Scoring
..summary:Set the value of the given matrix to the default value.
..signature:setDefaultScoreMatrix(scoreMatrix, tag)
..param.scoreMatrix:The @Spec.Score Matrix@ to set.
...type:Spec.Score Matrix
..param.tag:The tag to specify the matrix.
...type:Shortcut.Blosum30
...type:Shortcut.Blosum62
...type:Shortcut.Blosum80
..include:seqan/score.h
 */
template <typename TValue, typename TSequenceValue, typename TSpec, typename TTag>
inline void
setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TTag) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    TValue const * tab = _ScoringMatrixData<TValue, TSequenceValue, TTag>::getData();
    arrayCopy(tab, tab + TScore::TAB_SIZE, sc.data_tab);
}


/**
.Function.setDefaultScoreMatrix
..param.tag:
...type:Tag.Default
...remark:If @Tag.Default@, then the matrix will be filled with default constructed $TValue$ values.
..include:seqan/score.h
 */
template <typename TValue, typename TSequenceValue, typename TSpec>
inline void
setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, Default) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    arrayFill(sc.data_tab, sc.data_tab + TScore::TAB_SIZE, TValue());
}


/*
.Function._sscanfValue:
..cat:Input/Output
..summary:Use sscanf to parse a value from a $char *$ buffer.
..signature:_sscanfValue(buffer, value)
..param.buffer:Buffer to parse into.
...type:const char *
..param.value:Variable to parse the value from $buffer$ to.
...type:unsigned int
..include:seqan/score.h
 */
inline void
_sscanfValue(const char * buf, unsigned int & val) {
    SEQAN_CHECKPOINT;
    std::sscanf(buf, "%u", & val);
}


/*
.Function._sscanfValue.param.value.type:int
..include:seqan/score.h
 */
inline void
_sscanfValue(const char * buf, int & val) {
    SEQAN_CHECKPOINT;
    std::sscanf(buf, "%i", & val);
}


/*
.Function._sscanfValue.param.value.type:float
..include:seqan/score.h
 */
inline void
_sscanfValue(const char * buf, float & val) {
    SEQAN_CHECKPOINT;
    std::sscanf(buf, "%f", & val);
}


/*
.Function._sscanfValue.param.value.type:double
..include:seqan/score.h
 */
inline void
_sscanfValue(const char * buf, double & val) {
    SEQAN_CHECKPOINT;
    std::sscanf(buf, "%lf", & val);
}


template <typename TFile, typename TMeta>
void
readMeta(TFile & fl, TMeta & meta, ScoreMatrixFile) {
    SEQAN_CHECKPOINT;
    clear(meta);
    if (_streamEOF(fl)) return;

    typedef typename Value<TMeta>::Type TValue;
    TValue c = _streamGet(fl);

    while (!_streamEOF(fl) && (c == '#')) {
        c = _streamGet(fl);
        _stream_appendLine(fl, meta, c);
        appendValue(meta, '\n');
    }
    _streamUnget(fl);
}


template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec>
void
read(TFile & fl, Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, ScoreMatrixFile) {
    // TODO(holtgrew): The following is not very stable, does not interpret lines as whitespace separated numbers but infers column widths from the labels.  Should be fixed.
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    typedef typename Value<TFile>::Type TFileValue;

    // Clear the matrix.
    arrayFill(sc.data_tab, sc.data_tab + TScore::TAB_SIZE, TValue());

    // Start reading.
    if (_streamEOF(fl)) return;

    TFileValue c = _streamGet(fl);
    String<TFileValue> s;

    // Search for alphabet line, the first line that does not start
    // with '#', e.g.
    // " A R N D C Q E G H I L K M F P S T W Y V B Z X *"
    do {
        clear(s);
        _stream_appendLine(fl, s, c);
    } while (!_streamEOF(fl) && (empty(s) || (s[0] == '#')));

    if (_streamEOF(fl)) return;

    // Build table to map the alphabet line to the TSequenceValue values.
    typedef typename Iterator<String<TFileValue>, Standard>::Type TIterator;
    TIterator it = begin(s);
    TIterator it_end = end(s);
    String<unsigned int> mapping_;
    String<unsigned int> column_;
    for (unsigned int i = 0; it != it_end; ++i) {
        if ((*it) != ' ') {
        // TODO(holtgrew): This kind of type conversion really necessary?
            unsigned int pos = (TSequenceValue) *it;  // Conversion TFileValue => TSequenceValue => integral
            appendValue(mapping_, pos);
            appendValue(column_, i);  // This marks the end of the column
        }
        ++it;
    }

    // Read the matrix itself.
    while (!_streamEOF(fl)) {
        clear(s);
        _stream_appendLine(fl, s, c);

        if (empty(s) || (s[0] == '#')) continue;  // Skip empty lines and comments.

        // Read first character = alphabet column.
        // TODO(holtgrew): This kind of type conversion really necessary?
        unsigned int row = (TSequenceValue) s[0];  // Conversion TFileValue => TSequenceValue => integral
        unsigned int offset = row * TScore::VALUE_SIZE;

        // Read rest of the line.
        unsigned int right;
        unsigned int left = 0;

        TFileValue buf[100];  // 100 is enough, believe me!
        buf[99] = 0;

        for (unsigned int ii = 0; ii < length(column_); ++ii) {
            // Read column ii.
            //
            // Scan cell into buffer.
            right = column_[ii];
            TFileValue * it;
            for (it = buf + 99; it >= buf; --right) {
                if (right <= left) break;
                if (s[right] == ' ') break;
                --it;
                *it = s[right];
            }

            // Parse buffer.
            TValue val;
            _sscanfValue(it, val);

            sc.data_tab[offset + mapping_[ii]] = val;

            left = column_[ii];
        }
    }
}


template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec>
inline void
read(TFile & fl, Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc) {
    SEQAN_CHECKPOINT;
    read(fl, sc, ScoreMatrixFile());
}


/**
.Function.loadScoreMatrix
..cat:Input/Output
..summary:Load a score matrix from a file.
..signature:loadScoreMatrix(score, filename)
..remark:TODO, comment this better/at all.
..includes:seqan/score.h
..include:seqan/score.h
**/
template <typename TValue, typename TSequenceValue, typename TSpec, typename TString>
inline void
loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TString & filename) {
    SEQAN_CHECKPOINT;
    FILE * fl;
    _streamOpen(fl, filename);
    read(fl, sc);
    _streamClose(fl);
}


/**
.Function.loadScoreMatrix
..signature:loadScoreMatrix(score, filename, meta)
..remark:TODO, comment this better/at all.
..includes:seqan/score.h
..include:seqan/score.h
**/
template <typename TValue, typename TSequenceValue, typename TSpec, typename TString, typename TMeta>
inline void
loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TString & filename, TMeta & meta) {
    SEQAN_CHECKPOINT;
    FILE * fl;
    _streamOpen(fl, filename);
    readMeta(fl, meta, ScoreMatrixFile());
    read(fl, sc, ScoreMatrixFile());
    _streamClose(fl);
}


/*
.Function._sprintfValue:
..cat:Input/Output
..summary:Use sprintf to print a value into a $const char *$ buffer.
..signature:_sprintf(buffer, value)
..param.buffer:Buffer to write to.
...type:char *
...remark:Must be of sufficient size.
..param.value:Variable for which to write the string value to $buffer$.
...type:unsigned int
..includes:seqan/score.h
..include:seqan/score.h
 */
inline void
_sprintfValue(char * buf, unsigned int val) {
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): sprintf is unsafe, the C++ idiom is to use a std::string, we should probably use String.
    std::sprintf(buf, "%u", val);
}


/*
.Function._sprintfValue.param.value.type:int
..include:seqan/score.h
 */
inline void
_sprintfValue(char * buf, int val) {
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): sprintf is unsafe, the C++ idiom is to use a std::string, we should probably use String.
    std::sprintf(buf, "%d", val);
}


/*
.Function._sprintfValue.param.value.type:float
..include:seqan/score.h
 */
inline void
_sprintfValue(char * buf, float val) {
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): sprintf is unsafe, the C++ idiom is to use a std::string, we should probably use String.
    double d = val;
    std::sprintf(buf, "%G", d);
}


/*
.Function._sprintfValue.param.value.type:float
..include:seqan/score.h
 */
inline void
_sprintfValue(char * buf, double val) {
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): sprintf is unsafe, the C++ idiom is to use a std::string, we should probably use String.
    std::sprintf(buf, "%G", val);
}


/*
.Function._writeScoringMatrix:
..cat:Input/Output
..summary:Write the data of a scoring matrix to a file.
..signature:_writeScoringMatrix(file, table, meta)
..remark:TODO(holtgrew):More documentation.
..include:seqan/score.h
*/
template <typename TSequenceValue, typename TFile, typename TValue, typename TMeta>
void
_writeScoringMatrix(TFile & fl, TValue * tab, TMeta & meta) {
    SEQAN_CHECKPOINT;
    typedef typename Value<TFile>::Type TFileValue;

    enum {
        VALUE_SIZE = ValueSize<TSequenceValue>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    typedef typename Value<TFile>::Type TFileValue;

    // Write meta data.
    if (!empty(meta)) {
        bool line_begin = true;
        for (unsigned int i = 0; i < length(meta); ++i) {
            if (line_begin) {
                // Escape each line with a starting '#'.
                _streamPut(fl, '#');
                line_begin = false;
            }
            if (meta[i] == '\r') continue;
            if (meta[i] == '\n') line_begin = true;
            _streamPut(fl, meta[i]);
        }
        if (!line_begin) _streamPut(fl, '\n');
    }

    // Determine column width.
    unsigned int col_width = 1;
    char buf[100];  // 100 is enough, believe me!
    for (unsigned int i = 0; i < TAB_SIZE; ++i) {
        _sprintfValue(buf, tab[i]);
        unsigned int cell_width = std::strlen(buf);
        if (cell_width > col_width)
            col_width = cell_width;  // Compute maximum.
    }

    ++col_width;  // Increase col_width for an additional blank.

    // Write alphabet line.
    _streamPut(fl, ' ');  // A blank for alphabet column.
    for (unsigned int j = 0; j < VALUE_SIZE; ++j) {
        TFileValue val = (TSequenceValue) j;  // Conversion integral => TSequenceValue => TFileValue.
        // Leading blanks for column j.
        for (unsigned int k = 1; k < col_width; ++k) _streamPut(fl, ' ');
        _streamPut(fl, val);
    }
    _streamPut(fl, '\n');

    // Write rest of matrix.
    for (unsigned int i = 0; i < VALUE_SIZE; ++i) {
        // Write alphabet column cell.
        TFileValue val = (TSequenceValue) i;  // Conversion integral => TSequenceValue => TFileValue
        _streamPut(fl, val);

        // Write rest of line i.
        unsigned int offset = i * VALUE_SIZE;
        for (unsigned int j = 0; j < VALUE_SIZE; ++j) {
            _sprintfValue(buf, tab[offset + j]);
            unsigned int len = strlen(buf);

            // Leading blanks.
            for (unsigned int k = 0; k < col_width - len; ++k) _streamPut(fl, ' ');

            // Write cell.
            for (unsigned int k = 0; k < len; ++k) _streamPut(fl, buf[k]);
        }
        _streamPut(fl, '\n');
    }
}


/**
.Function.write:
..cat:Input/Output
..signature:write(file, scoreMatrix, meta)
..remark:TODO, comment better/at all.
..includes:seqan/score.h
..include:seqan/score.h
 */
template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec, typename TMeta>
inline void
write(TFile & fl, Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc, TMeta & meta) {
    SEQAN_CHECKPOINT;
    _writeScoringMatrix<TSequenceValue>(fl, sc.data_tab, meta);
}


/**
.Function.write:
..cat:Input/Output
..signature:write(file, scoreMatrix)
..remark:TODO, comment better/at all.
..includes:seqan/score.h
..include:seqan/score.h
 */
template <typename TFile, typename TValue, typename TSequenceValue, typename TSpec>
inline void
write(TFile & fl, Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc) {
    SEQAN_CHECKPOINT;
    write(fl, sc, "");
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SCORE_SCORE_MATRIX_H_
