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

#ifndef SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_
#define SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_

#include <seqan/sequence.h>

#include <seqan/arg_parse/arg_parse_type_support.h>

namespace seqan
{

// ----------------------------------------------------------------------------
// Class CommandLineArgument
// ----------------------------------------------------------------------------

/**
.Class.CommandLineArgument:
..cat:Miscellaneous
..summary:Stores information for a specific command line argument. It can be either an argument of
a CommandLineOption or directly an Argument on the command line.
..signature:CommandLineArgument
..remarks: .
..include:seqan/misc/misc_cmdparser.h
*/

/**
.Memfunc.CommandLineArgument#CommandLineArgument:
..class:Class.CommandLineArgument
..summary:Constructor
..signature:CommandLineOption (argumentType [, isListArgument, argumentLabel, numberOfArguments, _default])
..param.argumentType:A CommandLineArgument.ArgumentType value defining the type (e.g., String) of the
CommandLineArgument.
...tableheader:Flag|Description
...table:$CommandLineArgument::STRING$|Argument is a string
...table:$CommandLineArgument::INTEGER$|Argument is an integer
...table:$CommandLineArgument::DOUBLE$|A float
...table:$CommandLineArgument::INPUTFILE$|An input file
...table:$CommandLineArgument::OUTPUTFILE$|An output file
..param.isListArgument:Defines if the argument can be given multiple times.
...default:false.
..param.argumentLabel:Defines a user defined argument label for the help output. If this option is
not set, CommandLineArgument will automatically define a label based on the ArgumentType.
..param.numberOfArguments: Defines if the argument consists of defined number of elements (e.g., if
you want to provide an interval you would set this option to 2, so the parser knows that he needs
to search for exactly 2 values).
...default:1.
..param:_default:Sets the default value for this argument.
*/

class CommandLineArgument
{
public:
    enum ArgumentType
    {
                    // argument is
        STRING,     // ..  a string
        INTEGER,    // .. an integer
        DOUBLE,     // .. a float
        INPUTFILE,  // .. an inputfile (implicitly also a string, since paths/filenames are strings
        OUTPUTFILE  // .. an outputfile (implicitly also a string, since paths/filenames are strings)
    };


    // ----------------------------------------------------------------------------
    // Members to store type information
    // ----------------------------------------------------------------------------
    ArgumentType _argumentType;
    int          _numberOfArguments;
    CharString   _argumentLabel;
    bool         _isListArgument;

    // ----------------------------------------------------------------------------
    // Members to store the values
    // ----------------------------------------------------------------------------
    String<CharString>  defaultValue;
    String<CharString>  value;

    // ----------------------------------------------------------------------------
    // Members for restrictions
    // ----------------------------------------------------------------------------
    CharString            minValue;
    CharString            maxValue;
    StringSet<CharString> validValues;

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------
    CommandLineArgument(ArgumentType argumentType,
                        bool isListArgument = false,
                        CharString const & argumentLabel = "",
                        int numberOfArguments = 1) :
                        _argumentType(argumentType),
                        _numberOfArguments(numberOfArguments),
                        _argumentLabel(argumentLabel),
                        _isListArgument(isListArgument),
                        minValue(""),
                        maxValue("")
    {}


    template <typename TValue>
    CommandLineArgument(ArgumentType argumentType,
                        bool isListArgument,
                        CharString const & argumentLabel,
                        int numberOfArguments,
                        TValue const & _default) :
                        _argumentType(argumentType),
                        _numberOfArguments(numberOfArguments),
                        _argumentLabel(argumentLabel),
                        _isListArgument(isListArgument),
                        minValue(""),
                        maxValue("")
    {
        std::stringstream strm;
        strm << _default;
        appendValue(defaultValue, strm.str());
    }
};

// ----------------------------------------------------------------------------
// Function isListArgument()
// ----------------------------------------------------------------------------
/**
.Function.isListArgument
..summary:Returns whether the argument can be given multiple times.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..returns:$true$ if the argument argument can be given multiple times.
..see:Memfunc.CommandLineArgument#CommandLineArgument.param.isListArgument
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isListArgument(CommandLineArgument const & me)
{
    return me._isListArgument;
}

// ----------------------------------------------------------------------------
// Function isStringArgument()
// ----------------------------------------------------------------------------
/**
.Function.isStringArgument
..summary:Returns whether the argument is a string.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..returns:$true$ if the argument argument is a string argument.
..see:Memfunc.CommandLineArgument#CommandLineArgument.param.argumentType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isStringArgument(CommandLineArgument const & me)
{
    return (me._argumentType == CommandLineArgument::STRING) ||
           (me._argumentType == CommandLineArgument::INPUTFILE) ||
           (me._argumentType == CommandLineArgument::OUTPUTFILE);
}

// ----------------------------------------------------------------------------
// Function isIntegerArgument()
// ----------------------------------------------------------------------------
/**
.Function.isIntegerArgument
..summary:Returns whether the argument is an integer.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..returns:$true$ if the argument argument is an integer argument.
..see:Memfunc.CommandLineArgument#CommandLineArgument.param.argumentType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isIntegerArgument(CommandLineArgument const & me)
{
    return (me._argumentType == CommandLineArgument::INTEGER);
}

// ----------------------------------------------------------------------------
// Function isDoubleArgument()
// ----------------------------------------------------------------------------
/**
.Function.isDoubleArgument
..summary:Returns whether the argument is a double.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..returns:$true$ if the argument argument is a double argument.
..see:Memfunc.CommandLineArgument#CommandLineArgument.param.argumentType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isDoubleArgument(CommandLineArgument const & me)
{
    return (me._argumentType == CommandLineArgument::DOUBLE);
}

// ----------------------------------------------------------------------------
// Function isInputFileArgument()
// ----------------------------------------------------------------------------
/**
.Function.isInputFileArgument
..summary:Returns whether the argument is an input file.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..returns:$true$ if the argument argument is an input file argument.
..see:Memfunc.CommandLineArgument#CommandLineArgument.param.argumentType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isInputFileArgument(CommandLineArgument const & me)
{
    return (me._argumentType == CommandLineArgument::INPUTFILE);
}

// ----------------------------------------------------------------------------
// Function isOutputFileArgument()
// ----------------------------------------------------------------------------
/**
.Function.isOutputFileArgument
..summary:Returns whether the argument is an output file.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..returns:$true$ if the argument argument is an output file argument.
..see:Memfunc.CommandLineArgument#CommandLineArgument.param.argumentType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isOutputFileArgument(CommandLineArgument const & me)
{
    return (me._argumentType == CommandLineArgument::OUTPUTFILE);
}

// ----------------------------------------------------------------------------
// Function getArgumentLabel()
// ----------------------------------------------------------------------------
/**
.Function.getArgumentLabel
..summary:Returns the label for the given @Class.CommandLineArgument@. Either the user defined label
is returned or a default label (based on the ArgumentType is used).
..cat:Miscellaneous
..signature:getArgumentLabel(argument)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..returns:A $ShortCut.CharString$ containing the label.
..include:seqan/misc/misc_cmdparser.h
*/

inline CharString const
getArgumentLabel(CommandLineArgument const & me)
{
    if(me._argumentLabel != "")
    {
        return me._argumentLabel;
    }
    else
    {
        // infer from argument type
        CharString baseLabel = "";
        if(isInputFileArgument(me) || isOutputFileArgument(me)) baseLabel = "FILE";
        else if(isStringArgument(me)) baseLabel = "STR";
        else if(isIntegerArgument(me) || isDoubleArgument(me)) baseLabel = "NUM";

        CharString finalLabel;

        if(me._numberOfArguments != 1)
        {
            for(int i = 0; i < me._numberOfArguments; ++i)
            {
                if(i != 0) append(finalLabel, " ");
                append(finalLabel, baseLabel);
            }
        }
        else if(isListArgument(me)) finalLabel = baseLabel; // maybe we want to customize list labels
        else finalLabel = baseLabel;

        return finalLabel;
    }
}

// ----------------------------------------------------------------------------
// Helper Function _intervalAssert()
// ----------------------------------------------------------------------------

// this methods ensures that the given arguments define a non emtpy value interval
// otherwise it will trigger a SEQAN_CHECK failure
template<typename TIntervalBorder>
inline void
_intervalAssert(const CharString minValueAsString, const CharString maxValueAsString)
{
    if(minValueAsString != "" && maxValueAsString != "")
        SEQAN_CHECK(_cast<TIntervalBorder>(minValueAsString) < _cast<TIntervalBorder>(maxValueAsString),
                    "The interval [%s:%s] is empty. Please specify a valid, non-empty interval.",
                    toCString(minValueAsString),
                    toCString(maxValueAsString));
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/**
.Function.setMinValue
..summary:Sets the minimum value of a @Class.CommandLineArgument@ object.
..cat:Miscellaneous
..signature:setMinValue(argument,minValue)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..param.minValue:A @Shortcut.CharString@ containing a string representation of the minimum value
of the @Class.CommandLineArgument@.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setMinValue(CommandLineArgument & me, const CharString _minValue)
{
    if(isDoubleArgument(me))
    {
        SEQAN_CHECK(_isCastable<double>(_minValue), "The maximal value for a double argument must be double.");
        _intervalAssert<double>(_minValue, me.maxValue);
        me.minValue = _minValue;
    }
    else if(isIntegerArgument(me))
    {
        SEQAN_CHECK(_isCastable<int>(_minValue), "The maximal value for an integer argument must be an integer");
        _intervalAssert<int>(_minValue, me.maxValue);
        me.minValue = _minValue;
    }
    else
        SEQAN_FAIL("min/max values are not applicable to non numeric arguments");
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/**
.Function.setMaxValue
..summary:Sets the maximum value of a @Class.CommandLineArgument@ object.
..cat:Miscellaneous
..signature:setMaxValue(argument,maxValue)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..param.maxValue:A @Shortcut.CharString@ containing a string representation of the maximum value
of the @Class.CommandLineArgument@.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setMaxValue(CommandLineArgument & me, const CharString _maxValue)
{
    if(isDoubleArgument(me))
    {
        SEQAN_CHECK(_isCastable<double>(_maxValue), "The maximal value for a double argument must be double.");
        _intervalAssert<double>(me.minValue, _maxValue);
        me.maxValue = _maxValue;
    }
    else if(isIntegerArgument(me))
    {
        SEQAN_CHECK(_isCastable<int>(_maxValue), "The maximal value for an integer argument must be an integer");
        _intervalAssert<int>(me.minValue, _maxValue);
        me.maxValue = _maxValue;
    }
    else
        SEQAN_FAIL("min/max values are not applicable to non numeric arguments");
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/**
.Function.setValidValues
..summary:Sets the set of allowed values of a @Class.CommandLineArgument@ object.
..cat:Miscellaneous
..signature:setValidValues(argument,values)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..param.values:A $String<CharString>$ containing all valid entries for the option.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setValidValues(CommandLineArgument & me, StringSet<CharString> const & _values)
{
    if(isDoubleArgument(me) || isIntegerArgument(me))
        SEQAN_FAIL("CommandLineArgument does not support setting valid values for numeric arguments.");

    me.validValues = _values;
}

// ----------------------------------------------------------------------------
// Helper Function _isInInterval()
// ----------------------------------------------------------------------------

// check if the given value is in the provided interval
template<typename TTarget, typename TString>
inline bool
_isInInterval(TString value, TString lowerIntervalBound, TString upperIntervalBound)
{
    bool isInInterval = true;

    if(lowerIntervalBound != "")
        isInInterval &= (_cast<TTarget>(lowerIntervalBound) <= _cast<TTarget>(value));
    if(upperIntervalBound != "")
        isInInterval &= (_cast<TTarget>(value) <= _cast<TTarget>(upperIntervalBound));

    return isInInterval;
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

/**
.Function.assignValue
..summary:Assigns the given value (if applicable) to the @Class.CommandLineArgument@ object. If
the @Class.CommandLineArgument@ is a list or can hold multiple values
(@Memfunc.CommandLineArgument#CommandLineArgument.param.numberOfArguments@) the value will be appended.
Otherwise the value will be overwritten.
If the value can be assigned the method returns $true$, $false$ otherwise.
..cat:Miscellaneous
..signature:assignValue(argument,value)
..param.option:The @Class.CommandLineArgument@ object.
...type:Class.CommandLineArgument
..param.value:A @ShortCut.CharString>@ containing the value that should be assigned.
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
assignValue(CommandLineArgument & me, CharString const & value)
{
    typedef Iterator<StringSet<CharString>, Rooted>::Type TStringSetIterator;

    // type checks
    if(isIntegerArgument(me))
    {
        if(!_isCastable<int>(value))
            return false;
        if(!_isInInterval<int>(value, me.minValue, me.maxValue))
            return false;
    }

    if(isDoubleArgument(me))
    {
        if(!_isCastable<double>(value))
            return false;
        if(!_isInInterval<double>(value, me.minValue, me.maxValue))
            return false;
    }

    // check valid values
    if(isStringArgument(me))
    {
        if(!empty(me.validValues))
        {
            bool isContained = false;
            for(TStringSetIterator validValue = begin(me.validValues);
                validValue != end(me.validValues);
                goNext(validValue))
            {
                if(isInputFileArgument(me) || isOutputFileArgument(me))
                {
                    if (length(*validValue) > length(value))
                        continue;
                    else
                        isContained |= (suffix(value, length(value) - length(*validValue)) == *validValue);
                }
                else
                {
                    isContained |= (*validValue == value);
                }
                if(isContained) break;
            }
            if(!isContained) return false;
        }
    }

    // assignment
    if(isListArgument(me)) // just append
        appendValue(me.value, value, Exact());
    else if(me._numberOfArguments != 1)
    {

    }

}



} // namespace seqan

#endif // SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_

