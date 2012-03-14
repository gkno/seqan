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

#ifndef SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_
#define SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_

#include <string>
#include <vector>
#include <seqan/arg_parse/arg_parse_argument.h>

namespace seqan {

// ----------------------------------------------------------------------------
// Class ArgParseOption
// ----------------------------------------------------------------------------

/**
.Class.ArgParseOption:
..cat:Miscellaneous
..summary:Stores information for a specific command line option.
..signature:ArgParseOption
..remarks:A @Class.ArgParseOption@ object can be added to a @Class.CommandLineParser@ via @Function.addOption@.
..include:seqan/misc/misc_cmdparser.h
*/

/**
.Memfunc.ArgParseOption#ArgParseOption:
..class:Class.ArgParseOption
..summary:Constructor
..signature:ArgParseOption (shortName, longName, helpText [, argument])
..param.shortName:A std::string containing the short-name option identifier (e.g. $"h"$ for the $-h/--help$ option).
Although not suggested the short-name can contain more than 1 character.
...remarks:Note that the leading "-" is not passed.
..param.longName:A std::string containing the long-name option identifier (e.g. $"help"$ for the $-h/--help$ option).
...remarks:Note that the leading "--" is not passed.
..param.helpText:A std::string containing the help text associated with this option.
..param.argument:A ArgParseArgument for the option (e.g., an integer argument).
...type:Class.ArgParseArgument
*/

class ArgParseOption
{
public:
    // ----------------------------------------------------------------------------
    // Members to specify the names of the ArgParseOption
    // ----------------------------------------------------------------------------
    std::string         shortName;     // short option name
    std::string         longName;      // long option name

    // ----------------------------------------------------------------------------
    // Members representing type, content and restrictions of the ArgParseOption
    // ----------------------------------------------------------------------------
    ArgParseArgument    _argument;    // the argument of the option (if aplicable)
    bool                _isFlag;      // true if this a bool option, that has no
                                      // argument we will internally represent it as a
                                      // string option set to either "true" or "false"
    bool                _isRequired; // true if this ArgParseOption must be set
    bool                _isVisible;    // true if this ArgParseOption should not be
                                      // shown on the command line

    // ----------------------------------------------------------------------------
    // Members to help text
    // ----------------------------------------------------------------------------
    std::string         _helpText;    // The help text shown on the command
                                      // line

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------
    ArgParseOption(std::string const & _shortName,
                   std::string const & _longName,
                   std::string const & _help,
                   ArgParseArgument const & argument) :
                   shortName(_shortName),
                   longName(_longName),
                   _argument(argument),
                   _isFlag(false),
                   _isRequired(false),
                   _isVisible(false),
                   _helpText(_help)
    {}

    ArgParseOption(std::string const & _shortName,
                   std::string const & _longName,
                   std::string const & _help) :
                   shortName(_shortName),
                   longName(_longName),
                   _argument(ArgParseArgument::STRING, false, "", 1, "true"),
                   _isFlag(true),
                   _isRequired(false),
                   _isVisible(false),
                   _helpText(_help)
    {
        setValidValues(_argument, "true false");
    }
};

// ----------------------------------------------------------------------------
// Function isStringOption()
// ----------------------------------------------------------------------------

/**
.Function.isStringOption
..summary:Returns whether option argument can be a string.
..cat:Miscellaneous
..signature:isStringOption(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option argument can be a string.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isStringOption(ArgParseOption const & me)
{
    return isStringArgument(me._argument) && !me._isFlag;
}

// ----------------------------------------------------------------------------
// Function isBooleanOption()
// ----------------------------------------------------------------------------

/**
.Function.isBooleanOption
..summary:Returns whether option is a switch.
..cat:Miscellaneous
..signature:isBooleanOption(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is a switch.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isBooleanOption(ArgParseOption const & me)
{
    return me._isFlag;
}

// ----------------------------------------------------------------------------
// Function isDoubleOption()
// ----------------------------------------------------------------------------

/**
.Function.isDoubleOption
..summary:Returns whether option argument can be a double.
..cat:Miscellaneous
..signature:isDoubleOption(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option argument can be a double.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isDoubleOption(ArgParseOption const & me)
{
    return isDoubleArgument(me._argument);
}

// ----------------------------------------------------------------------------
// Function isIntOption()
// ----------------------------------------------------------------------------

/**
.Function.isIntOption
..summary:Returns whether option argument can be an integer.
..cat:Miscellaneous
..signature:isIntOption(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option argument can be an integer.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isIntOption(ArgParseOption const & me)
{
    return isIntegerArgument(me._argument);
}

// ----------------------------------------------------------------------------
// Function isVisible()
// ----------------------------------------------------------------------------

/**
.Function.isVisible
..summary:Returns whether option is hidden on the help screen.
..cat:Miscellaneous
..signature:isHiddenOption(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is hidden on the help screen.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isVisible(ArgParseOption const & me)
{
    return me._isVisible;
}

// ----------------------------------------------------------------------------
// Function setVisibility()
// ----------------------------------------------------------------------------

/**
.Function.setVisibility
..summary:Sets the visibility of the ArgParseOption.
..cat:Miscellaneous
..signature:setVisibility(option, visibility)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.visibility:The new visibility of the option.
...type:Bool
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setVisibility(ArgParseOption & me, bool visibility)
{
    me._isVisible = visibility;
}

// ----------------------------------------------------------------------------
// Function isRequired()
// ----------------------------------------------------------------------------

/**
.Function.isRequired
..summary:Returns whether the option is mandatory.
..cat:Miscellaneous
..signature:isRequired(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is mandatory.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isRequired(ArgParseOption const & me)
{
    return me._isRequired;
}

// ----------------------------------------------------------------------------
// Function setRequired()
// ----------------------------------------------------------------------------

/**
.Function.setRequired
..summary:Sets whether or not the option is mandatory.
..cat:Miscellaneous
..signature:setRequired(option, required)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.required:The new required value of the option.
...type:Bool
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setRequired(ArgParseOption & me, bool required)
{
    me._isRequired = required;
}

// ----------------------------------------------------------------------------
// Function isOptionList()
// ----------------------------------------------------------------------------

/**
.Function.isOptionList
..summary:Returns whether the option can be given multiple times.
..cat:Miscellaneous
..signature:isOptionList(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option can be given multiple times on command line.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isListOption(ArgParseOption const & me)
{
    return isListArgument(me._argument);
}

// ----------------------------------------------------------------------------
// Function isInputFile()
// ----------------------------------------------------------------------------

/**
.Function.isInputFile
..summary:Returns whether the argument of the given option is an input file.
..cat:Miscellaneous
..signature:isInputFile(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the argument of the option is an input file.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isInputFile(ArgParseOption const & me)
{
    return isInputFileArgument(me._argument);
}

// ----------------------------------------------------------------------------
// Function isOutputFile()
// ----------------------------------------------------------------------------

/**
.Function.isOutputFile
..summary:Returns whether the argument of the given option is an output file.
..cat:Miscellaneous
..signature:isOutputFile(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the argument of the option is an output file.
..see:Memfunc.ArgParseOption#ArgParseOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isOutputFile(ArgParseOption const & me)
{
    return isOutputFileArgument(me._argument);
}

// ----------------------------------------------------------------------------
// Helper Function _writeOptName()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
_writeOptName(TStream & target, ArgParseOption const & me)
{
    //IOREV _notio_ irrelevant for iorev
    _streamWrite(target, empty(me.shortName) ? "" : "-");
    _streamWrite(target, me.shortName);
    _streamWrite(target, (empty(me.shortName) || empty(me.longName)) ? "" : ", ");
    if (!empty(me.longName))
    {
        _streamWrite(target, "--");
        _streamWrite(target, me.longName);
    }
}

// ----------------------------------------------------------------------------
// Function write()                                           ArgParseOption
// ----------------------------------------------------------------------------

/**
.Function.write
..summary:Writes the basic information about the @Class.ArgParseOption@ to the provided stream.
..cat:Miscellaneous
..signature:write(stream,option)
..param.stream:The target stream.
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TStream>
inline void
write(TStream & target, ArgParseOption const & me)
{
    //IOREV _nodoc_ this specialization is not documented
    _streamPut(target, '\t');
    _writeOptName(target, me);
    _streamPut(target, '\t');
    _streamPut(target, '\t');
    _streamWrite(target, me._helpText);
}

// ----------------------------------------------------------------------------
// operator<<()                                               ArgParseOption
// ----------------------------------------------------------------------------

template <typename TStream>
inline TStream &
operator<<(TStream & target, ArgParseOption const & source)
{
    //IOREV _nodoc_ this specialization is not documented
    write(target, source);
    return target;
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/**
.Function.setMinValue
..summary:Sets the minimum value of a @Class.ArgParseOption@ object.
..cat:Miscellaneous
..signature:setMinValue(option,minValue)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.minValue:A std::string containing a string representation of the minimum value of the @Class.ArgParseOption@.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setMinValue(ArgParseOption & me, const std::string _minValue)
{
    setMinValue(me._argument, _minValue);
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/**
.Function.setMaxValue
..summary:Sets the maximum value of a @Class.ArgParseOption@ object.
..cat:Miscellaneous
..signature:setMaxValue(option,maxValue)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.maxValue:A std::string containing a string representation of the maximum value of the @Class.ArgParseOption@.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setMaxValue(ArgParseOption & me, const std::string _maxValue)
{
    setMaxValue(me._argument, _maxValue);
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/**
.Function.setValidValues
..summary:Sets the set of allowed values of a @Class.ArgParseOption@ object.
..cat:Miscellaneous
..signature:setValidValues(option,values)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.values:A std::vector<std::string> containing all valid entries for the option or a
std::string with valid values separated by spaces.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setValidValues(ArgParseOption & me, std::vector<std::string> const & _values)
{
    setValidValues(me._argument, _values);
}

inline void
setValidValues(ArgParseOption & me, std::string const & _values)
{
    setValidValues(me._argument, _values);
}

// ----------------------------------------------------------------------------
// Function getArgumentLabel()
// ----------------------------------------------------------------------------
/**
.Function.getArgumentLabel
..summary:Returns the label for the given @Class.ArgParseArgument@. Either the user defined label
is returned or a default label (based on the ArgumentType is used).
..cat:Miscellaneous
..signature:getArgumentLabel(argument)
..param.option:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:A $ShortCut.std::string$ containing the label.
..include:seqan/arg_parse.h
*/

inline std::string const
getArgumentLabel(ArgParseOption const & me)
{
    if(isBooleanOption(me))
        return "";
    else
        return getArgumentLabel(me._argument);
}

// ----------------------------------------------------------------------------
// Function isSet()
// ----------------------------------------------------------------------------
/**
.Function.isSet
..summary:Returns whether the option was set on the command line or not.
..cat:Miscellaneous
..signature:isSet(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseArgument
..returns:True if the option was used on the command line, otherwise false.
..include:seqan/arg_parse.h
*/

inline bool
isSet(ArgParseOption const & me)
{
    return isSet(me._argument);
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

/**
.Function.assignValue
..summary:Assigns the given value (if applicable) to the @Class.ArgParseOption@ object. If
the @Class.ArgParseOption@ is a list or can hold multiple values
(@Memfunc.ArgParseArgument#ArgParseArgument.param.numberOfArguments@) the value will be appended.
Otherwise the value will be overwritten.
..cat:Miscellaneous
..signature:assignArgumentValue(option,value [, argNo])
..param.argument:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.value:A std::string containing the value that should be assigned.
..include:seqan/arg_parse.h
*/

inline void
assignArgumentValue(ArgParseOption & me, std::string const & value) throw (ParseException)
{
    assignArgumentValue(me._argument, value);
}

// ----------------------------------------------------------------------------
// Function getArgumentValue()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValue
..summary:Returns the value of the @Class.ArgParseOption@ object. If
the @Class.ArgParseOption@ is a list or can hold multiple values
(@Memfunc.ArgParseArgument#ArgParseArgument.param.numberOfArguments@) you can specify which value
you want to get. If not set the first value will be returned.
..cat:Miscellaneous
..signature:getArgumentValue(option [, position])
..param.argument:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.position:A unsigned int defining the which value should be returned.
..returns:The value set at position $position$.
..include:seqan/arg_parse.h
*/

inline std::string const &
getArgumentValue(ArgParseOption const & me, unsigned position)
{
    return getArgumentValue(me._argument, position);
}

inline std::string const &
getArgumentValue(ArgParseOption const & me)
{
    return getArgumentValue(me, 0);
}

// ----------------------------------------------------------------------------
// Function getArgumentValues()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValues
..summary:Returns all values of the @Class.ArgParseOption@ object as const std::vector.
..cat:Miscellaneous
..signature:getArgumentValues(option)
..param.argument:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$std::vector<std::string>$ containing the values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const &
getArgumentValues(ArgParseOption const & me)
{
    return getArgumentValues(me._argument);
}

// ----------------------------------------------------------------------------
// Function numberOfArguments
// ----------------------------------------------------------------------------

/**
.Function.numberOfArguments
..summary:Returns the number of arguments for this @Class.ArgParseOption@.
..cat:Miscellaneous
..signature:numberOfArguments(argument)
..param.argument:The @Class.ArgParseOption@ object.
...type:Class.ArgParseArgument
..returns:The number of allowed arguments for this @Class.ArgParseOption@.
..include:seqan/arg_parse.h
*/

inline unsigned
numberOfArguments(ArgParseOption const & me)
{
    return numberOfArguments(me._argument);
}



} // namespace seqan

#endif // SANDBOX_ARG_PARSE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_
