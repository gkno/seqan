/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 207-010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  Author: Stefan Aiche <aiche@fu-berlin.de>
 ============================================================================
  Wildcard pattern matching using a modification of the Shift-And Algorithm.
 ==========================================================================*/

#ifndef SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_
#define SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_

namespace seqan {

struct _WildShiftAnd;
typedef Tag<_WildShiftAnd> WildShiftAnd;


template <typename _TNeedle>
struct Pattern<_TNeedle, WildShiftAnd> : _FindState {
    typedef _TNeedle TNeedle;

    // The pattern's state.
    TState _state;

    // The needle we store.
    Holder<TNeedle> _host;

    Pattern() : _state(STATE_EMPTY) {
        SEQAN_CHECKPOINT;
    }

    explicit
    Pattern(TNeedle & ndl)
        : _state(STATE_INITIAL),
          _host(ndl) {
        SEQAN_CHECKPOINT;
    }
};


template <typename TNeedle>
struct Needle<Pattern<TNeedle, WildShiftAnd> > {
    typedef typename Value<TNeedle>::Type Value;
};


template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


// TODO(holtgrew): Workaround for default host implementation for non-const type.
template <typename TNeedle>
TNeedle const & host(Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return host(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
TNeedle const & needle(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    return value(pattern._host);
}


template <typename TNeedle>
typename Position<TNeedle>::Type length(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return infix(needle(pattern), 0u, length(needle(pattern)));
}


// TODO(holtgrew): Workaround for default infix implementation for non-const type.
template <typename TNeedle>
Segment<TNeedle const, InfixSegment> infix(Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return infix(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, WildShiftAnd> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return begin(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default begin implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type begin(Pattern<TNeedle, WildShiftAnd> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return begin(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, WildShiftAnd> const & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return end(needle(pattern), spec);
}


// TODO(holtgrew): Workaround for default end implementation for non-const type.
template <typename TNeedle, typename TTag>
typename Iterator<TNeedle const, Tag<TTag> const>::Type end(Pattern<TNeedle, WildShiftAnd> & pattern, Tag<TTag> const & spec) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return end(const_cast<TPattern const &>(pattern), spec);
}


template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return 0u;
}


// TODO(holtgrew): Workaround for default beginPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type beginPosition(Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return beginPosition(const_cast<TPattern const &>(pattern));
}


template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, WildShiftAnd> const & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    SEQAN_ASSERT_EQ(TPattern::STATE_BEGIN_FOUND, pattern._state);
    return length(needle(pattern));
}


// TODO(holtgrew): Workaround for default endPosition implementation for non-const type.
template <typename TNeedle>
typename Position<TNeedle>::Type endPosition(Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return endPosition(const_cast<TPattern const &>(pattern));
}


template <typename THaystack, typename TNeedle>
bool find(Finder<THaystack, Default> & finder,
          Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    typedef typename Position<TNeedle>::Type TPosition;

    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);

    // Do not continue if the state is "not found".
    if (finder._state == TPattern::STATE_NOTFOUND)
        return false;
    // Initialize finder if state is "initial".  Otherwise advance at
    // least by one (if not set to of haystack with setEndPosition()).
    if (finder._state == TPattern::STATE_INITIAL) {
        finder._beginPosition = 0u;
        finder._endPosition = length(needle(pattern));
    } else if (finder._state == TPattern::STATE_NO_HIT) {
        // Only advance if not at end if set manually to a "no hit" position.
        if (finder._endPosition == length(haystack(finder)))
            return false;
        finder._beginPosition += 1;
    } else {
        finder._beginPosition += 1;
    }

    // Search the needle in the haystack naively.
    for (TPosition i = 0u; i < length(needle(pattern));) {
        // Break out of loop if no more match is possible.
        if (finder._beginPosition >= length(haystack(finder)) - length(needle(pattern))) {
            finder._state = TFinder::STATE_NOTFOUND;
            pattern._state = TPattern::STATE_NOTFOUND;
            return false;
        }
        // Otherwise, go on searching.
        if (needle(pattern)[i] != haystack(finder)[finder._beginPosition + i]) {
            finder._beginPosition += 1;
            i = 0u;
            continue;
        }
        i += 1;
    }
    finder._endPosition = finder._beginPosition + length(needle(pattern));
    finder._state = TFinder::STATE_BEGIN_FOUND;
    pattern._state = TPattern::STATE_BEGIN_FOUND;
    return true;
}


template <typename THaystack, typename TNeedle>
bool findBegin(Finder<THaystack, Default> & finder,
               Pattern<TNeedle, WildShiftAnd> & pattern) {
    SEQAN_CHECKPOINT;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    return finder._state == TPattern::STATE_BEGIN_FOUND;
}


template <typename THaystack, typename TNeedle, typename TPosition>
bool setEndPosition(Finder<THaystack, Default> & finder,
                    Pattern<TNeedle, WildShiftAnd> & pattern,
                    TPosition const & pos) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    // State of finder and pattern should be in sync.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // End position must not be right of the end of the haystack.
    SEQAN_ASSERT_LEQ(static_cast<typename _MakeUnsigned<TPosition>::Type>(pos), length(haystack(finder)));
    // Begin position must not be left of the beginning of the haystack.
    SEQAN_ASSERT_GEQ(static_cast<typename _MakeUnsigned<TPosition>::Type>(pos), length(needle(pattern)));

    // Set the end position.
    finder._endPosition = pos;
    finder._beginPosition = pos - length(needle(pattern));

    // Check whether there is a hit at this position and update the
    // state accordingly.
    typedef typename Position<THaystack>::Type THaystackPos;
    for (THaystackPos i = 0u; i < length(needle(pattern)); ++i) {
        if (needle(pattern)[i] != haystack(finder)[finder._beginPosition + i]) {
            finder._state = TPattern::STATE_NO_HIT;
            pattern._state = TFinder::STATE_NO_HIT;
            return false;
        }
    }
    finder._state = TPattern::STATE_BEGIN_FOUND;
    pattern._state = TPattern::STATE_BEGIN_FOUND;
    return true;
}


// Build the alignment resulting from the search result as specified by the
// finder and the pattern.  If the state is not "begin found" then no alignment
// is built and false is returned.
template <typename THaystack, typename TNeedle, typename TAlignSeq, typename TAlignSpec>
bool buildAlignment(Finder<THaystack, Default> &finder,
                    Pattern<TNeedle, WildShiftAnd> &pattern,
                    Align<TAlignSeq, TAlignSpec> &outAlignment) {
    SEQAN_CHECKPOINT;
    typedef Finder<THaystack, Default> TFinder;
    typedef Pattern<TNeedle, WildShiftAnd> TPattern;
    typedef Align<TAlignSeq, TAlignSpec> TAlign;
    typedef typename Row<TAlign>::Type TRow;

    // Both finder and pattern must be in the "found begin position" state.
    SEQAN_ASSERT_EQ(finder._state, pattern._state);
    // Can only build alignment if the state is "begin found".
    if (finder._state != TFinder::STATE_BEGIN_FOUND)
        return false;

    // Initialize alignment with the two sequences.
    resize(rows(outAlignment), 2);
    assignSource(row(outAlignment, 0), haystack(finder));
    assignSource(row(outAlignment, 1), needle(pattern));
    insertGaps(row(outAlignment, 1), 0, beginPosition(finder));

    return true;
}


// TODO(holtgrew): Should probably to into some utility header.
inline bool _isUnsigned(CharString const & number) {
    SEQAN_CHECKPOINT;
    if (length(number) == 0u)
        return false;

    typedef Position<CharString>::Type TPosition;
	for (TPosition i = 0; i < length(number); ++i) {
        if (number[i] > '9' || number[i] < '0')
            return false;
	}
	return true;
}


// The states used by the finite state machine that validates the regular
// expressions of the WildShiftAnd algorithm.
//
// These regular expressions consist of sequences of characters, character
// classes or the placeholder ".".  Each can be followed by a quantifier,
// where quantifiers are characters from the set {*, +, ?} or explicit numbers
// as in a{i, j}.
enum _Find_WildShiftAnd_ParserStates {
    STATE_NO_CHAR,                    // Expecting character, character class or ".".
    STATE_CHAR,                       // Read character, character class or ".".
    STATE_ESCAPED,                    // Just read a backslash, next character is escaped.
    STATE_CHAR_CLASS,                 // In a character class, read at least one character.
    STATE_CHAR_CLASS_ESCAPED,         // Escape character in character class, a character must follow.
    STATE_CHAR_CLASS_DASH,            // Just read a dash in a character class.
    STATE_CHAR_CLASS_NO_CHAR,         // In a character class, expecting a char, no "-".
    STATE_QUANTIFIER_BEGIN,           // Just read "{" from a quantifier.
    STATE_QUANTIFIER_NUM1,            // Read at least one character from the first number.
    STATE_QUANTIFIER_COMMA,           // Read comma, expecting first character of second number in quantifier.
    STATE_QUANTIFIER_NUM2,            // Read at least one character from second number.
};


// Check whether the pattern is valid.
//
// We use a finite state machine for the validation with the states of
// the enum _Find_WildShiftAnd_ParserStates.
inline bool _find_WildShiftAnd_isValid(CharString const & needle) {
    SEQAN_CHECKPOINT;

    if (length(needle) == 0)
        return false;

    typedef Iterator<CharString, Standard>::Type TIterator;

    _Find_WildShiftAnd_ParserStates state = STATE_NO_CHAR;

    // TODO(holtgrew): Check i <= j in a{i, j}?

    for (TIterator it = begin(needle, Standard()); it != end(needle, Standard()); ++it) {
        switch (state) {
            case STATE_NO_CHAR:
                if (*it == '\\') {
                    state = STATE_ESCAPED;
                } else if (*it == '[') {
                    state = STATE_CHAR_CLASS_NO_CHAR;
                } else if (*it != '+' && *it != '*' && *it != '?' && *it != '{' && *it != '}' && *it != '-') {
                    // No character, char class or "." here, must not be quantifier.
                    state = STATE_CHAR;
                } else {
                    return false;
                }
                break;
            case STATE_CHAR:
                if (*it == '\\') {
                    state = STATE_ESCAPED;
                } else if (*it == '{') {
                    state = STATE_QUANTIFIER_BEGIN;
                } else if (*it == '[') {
                    state = STATE_CHAR_CLASS_NO_CHAR;
                } else if (*it == '+' || *it == '*' || *it == '?') {
                    state = STATE_NO_CHAR;
                } else if (*it == ']' || *it == '}' || *it == '-') {
                    // Other special character have to be escaped.
                    return false;
                }
                break;
            case STATE_ESCAPED:
                state = STATE_CHAR;
                break;
            case STATE_CHAR_CLASS:
                if (*it == ']') {
                    state = STATE_CHAR;
                } else if (*it == '\\') {
                    state = STATE_CHAR_CLASS_ESCAPED;
                } else if (*it == '-') {
                    state = STATE_CHAR_CLASS_DASH;
                } else if (*it == '+' || *it == '*' || *it == '?' || *it == '{' || *it == '}' || *it == '.' || *it == '[') {
                    // Other special character have to be escaped.
                    return false;
                }
                break;
            case STATE_CHAR_CLASS_ESCAPED:
                state = STATE_CHAR_CLASS;
                break;
            case STATE_CHAR_CLASS_DASH:
                if (*it == '+' || *it == '*' || *it == '?' || *it == '{' ||
                    *it == '}' || *it == '.' || *it == '[' || *it == ']' ||
                    *it == '-') {
                    return false;
                } else if (*it == '\\') {
                    state = STATE_CHAR_CLASS_ESCAPED;
                } else {
                    state = STATE_CHAR_CLASS_NO_CHAR;
                }
                break;
            case STATE_CHAR_CLASS_NO_CHAR:
                if (*it == ']') {
                    state = STATE_CHAR;
                } else if (*it == '\\') {
                    state = STATE_CHAR_CLASS_ESCAPED;
                } else if (*it == '+' || *it == '*' || *it == '?' || *it == '{' || *it == '}' || *it == '.' || *it == '-' || *it == '[') {
                    return false;
                } else {
                    state = STATE_CHAR_CLASS;
                }
                break;
            case STATE_QUANTIFIER_BEGIN:
                if (*it >= '0' && *it <= '9') {
                    state = STATE_QUANTIFIER_NUM1;
                } else {
                    return false;
                }
                break;
            case STATE_QUANTIFIER_NUM1:
                if (*it == ',') {
                    state = STATE_QUANTIFIER_COMMA;
                } else if (*it < '0' || *it > '9') {
                    return false;
                }
                break;
            case STATE_QUANTIFIER_COMMA:
                if (*it >= '0' && *it <= '9') {
                    state = STATE_QUANTIFIER_NUM2;
                } else {
                    return false;
                }
                break;
            case STATE_QUANTIFIER_NUM2:
                if (*it == '}') {
                    state = STATE_NO_CHAR;
                } else if (*it < '0' || *it > '9') {
                    return false;
                }
                break;
        }
    }

    if (state == STATE_CHAR || state == STATE_NO_CHAR)
        return true;
    return false;
}


// Determine the pattern length without wildcard characters.  needle
// must be a valid pattern.
inline Position<CharString>::Type _find_WildShiftAnd_lengthWithoutWildcards(CharString const & needle) {
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_TRUE(_find_WildShiftAnd_isValid(needle));

    typedef Position<CharString>::Type TPosition;
    TPosition result = 0u;

    for (TPosition i = 0u; i < length(needle); ++i) {
        char c = needle[i];

        // Used in the switch statement below, must be defined before.
        CharString number;
        unsigned n, m;

        switch (c) {
            case '+':
            case '*':
            case '?':
                // Skip quantifier characters.
                continue;
            case '[':
                // Skip over character set, means one character.
                while (needle[i] != ']')
                    ++i;
                result += 1;
                break;
            case '{':
                // Skip over explicit quantifier, extract number of quantified
                // characters.
                ++i;
                // Read first number.
                while (needle[i] != ',') {
                    appendValue(number, needle[i]);
                    i += 1;
                }
                n = atoi(toCString(number));
                i += 1;  // Skip comma.
                // Read second number.
                clear(number);
                while (needle[i] != '}') {
                    appendValue(number, needle[i]);
                    i += 1;
                }
                m = atoi(toCString(number));

                SEQAN_ASSERT_GEQ(m, n);
                result += m - 1;
                break;
            case '\\':
                // Escape-sequences count as one character.
                result += 1;
                i += 1;
            default:
                result += 1;
                break;
        }
    }

    return result;
}


// Build a string with all characters in the range [begin, end) in host into
// result.
//
// TODO(holtgrew): We could simply use segments in the caller and a template argument for host instead of begin/end here.
inline void _find_WildShiftAnd_getCharacterClass(
        CharString & result, CharString const & host,
        Position<CharString>::Type begin, Position<CharString>::Type end) {
    SEQAN_CHECKPOINT;
    typedef Position<CharString>::Type TPosition;
    clear(result);

    for (TPosition pos = begin; pos < end; ++pos) {
        if (host[pos] == '\\') {
            SEQAN_ASSERT_LT(pos + 1, end);
            pos += 1;
            appendValue(result, host[pos]);
        } else if (host[pos] != '-') {
            appendValue(result, host[pos]);
        } else {
            // Read end of range character.
            SEQAN_ASSERT_LT(pos + 1, end);
            char last;
            if (host[pos + 1] == '\\') {
                SEQAN_ASSERT_LT(pos + 2, end);
                pos += 2;
                last = host[pos];
            } else {
                pos += 1;
                last = host[pos];
            }
            SEQAN_ASSERT_LEQ(back(result), last);
            for (char c = back(result) + 1; c <= last; ++c)
                appendValue(result, c);
        }
    }
}

}  // namespace seqan
          
#endif  // SEQAN_FIND2_FIND_PATTERN_WILD_SHIFTAND_H_
