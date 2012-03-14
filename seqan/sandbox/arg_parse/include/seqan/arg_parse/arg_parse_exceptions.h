// ==========================================================================
//                           arg_parse_exceptions.h
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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_EXCEPTIONS_H_
#define SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_EXCEPTIONS_H_

#include <string>
#include <sstream>
#include <exception>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

class ParseException: public std::exception
{
protected:
    std::string _what;

public:
    ParseException(std::string const & what)
        : _what(what)
    {
    }

    // we need to define this one to avoid looser throw specifier error
    virtual ~ParseException() throw()
    {}

    virtual const char* what() const throw()
    {
        return _what.c_str();
    }
};


class InvalidOptionException : public ParseException
{

public:
    InvalidOptionException(std::string const & option)
        : ParseException("")
    {
        std::stringstream what;
        what << "illegal option -- " << option;
        _what = what.str();
    }

    // we need to define this one to avoid looser throw specifier error
    virtual ~InvalidOptionException() throw()
    {}
};

class MissingArgumentException : public ParseException
{
public:
    MissingArgumentException(std::string const & option)
        : ParseException("")
    {
        std::stringstream what;
        what << "option requires an argument -- " << option;
        _what = what.str();
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_EXCEPTIONS_H_
