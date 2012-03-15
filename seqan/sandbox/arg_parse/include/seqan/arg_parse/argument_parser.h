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

#ifndef SANDBOX_ARG_PARSE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_
#define SANDBOX_ARG_PARSE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_

#include <sstream>
#include <seqan/map.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <seqan/arg_parse/arg_parse_argument.h>
#include <seqan/arg_parse/arg_parse_option.h>
#include <seqan/arg_parse/arg_parse_type_support.h>

#include <seqan/misc/misc_terminal.h>
#include <seqan/misc/tool_doc.h>

#include <iostream>
#include <fstream>

namespace seqan
{

/**
.Class.ArgumentParser
..cat:Miscellaneous
..summary:Stores multiple @Class.ArgParseOption@ objects and parses the command line _arguments for these options.
..signature:ArgumentParser
..include:seqan/arg_parse.h
..remarks:
See the documentation of @Class.ToolDoc@ on how to format text.
Where possible, formatting is added automatically for you.
You have to use formatting in the following places: (1) usage lines, (2) option help texts, (3) description and additional text sections.
..example.text:
The following gives a simple example of how to use the @Class.ArgumentParser@.
..example.code:
ArgumentParser parser("alf");
setShortDescription(parser, "Alignment free sequence comparison");
setVersion(parser, "1.0");
setDate(parser, "Jan 2010");


addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN\\fP \\fB-o\\fP \\fIOUT\\fP");

addDescription(parser,
               "ALF can be used to calculate the pairwise similarity of sequences "
               "using alignment-free methods.  All methods which are implemented are "
               "based on k-mer counts.");

ArgParseOption optionInputFile("i", "inputFile", "Name of the multi-FASTA input.",
                                  OptionType::String | OptionType::Mandatory);
addOption(parser, ArgParseOption("i", "inputFile", "Name of the multi-FASTA input.", ));

ArgParseOption optionInputFile("o", "outputFile", "Name of the output file.",
                                  OptionType::String | OptionType::Mandatory);
optionInputFile = addArgumentText(optionInputFile, "OUT");
addOption(parser, optionInputFile);

addTextSection(parser, "See Also");
addText(parser, "http://www.seqan.de/projects/alf");
..see:Class.ToolDoc

.Memfunc.ArgumentParser#ArgumentParser
..class:Class.ArgumentParser
..summary:Constructor
..signature:ArgumentParser ()
..signature:ArgumentParser (applicationName)
..param.applicationName:A std::string containing the name of the application.
..remarks:If the name of the application is not passed to the constructor it will be extracted from the command line.
*/

class ArgumentParser
{
public:

    // ----------------------------------------------------------------------------
    // Enum ParseResult
    // ----------------------------------------------------------------------------
    // will be used as return value of parse(..) to indicate whether parsing worked
    enum ParseResult
    {
        OK,
        ERROR,
        HELP,
        VERSION,
        WRITE_CTD,
        EXPORT_HELP
    };

    // ----------------------------------------------------------------------------
    // Class Typedefs
    // ----------------------------------------------------------------------------
    typedef std::vector<ArgParseOption>   TOptionMap;
    typedef std::vector<ArgParseArgument> TArgumentMap;
    typedef Size<TOptionMap>::Type        TOptionMapSize;
    typedef Size<TArgumentMap>::Type      TArgumentMapSize;

    typedef std::map<std::string, TOptionMapSize> TStringMap;
    typedef std::vector<std::string>              TValueMap;

    // ----------------------------------------------------------------------------
    // Mapping of option names to options
    // ----------------------------------------------------------------------------
    TStringMap   shortNameMap;
    TStringMap   longNameMap;
    TOptionMap   optionMap;
    TArgumentMap argumentList;

    // ----------------------------------------------------------------------------
    // Members
    // ----------------------------------------------------------------------------
    std::vector<std::string> _description;
    std::string              _appName;
    std::vector<std::string> _titleText;
    std::vector<std::string> _usageText;
    std::vector<std::string> _versionText;
    
    ToolDoc _toolDoc;

    // ----------------------------------------------------------------------------
    // return values for unset parameters
    // ----------------------------------------------------------------------------
    const std::string              _null;
    const std::vector<std::string> _nullSet;

    // friend declaration to make addOption() available in init function
    friend inline void addOption(ArgumentParser & me, ArgParseOption const & opt);

    // ----------------------------------------------------------------------------
    // Function init()
    // ----------------------------------------------------------------------------
    void init()
    {
        addOption(*this, ArgParseOption("h", "help", "Displays this help message."));
        addOption(*this, ArgParseOption("", "write-ctd", "Exports the app's interface description to a .ctd file.", ArgParseArgument(ArgParseArgument::OUTPUTFILE)));
        addOption(*this, ArgParseOption("", "export-help", "Export help to a format. One of {'html', 'man', 'txt'}.", ArgParseArgument(ArgParseArgument::STRING, false, "FORMAT")));
    }

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    ArgumentParser()
    {
        init();
    }

    ArgumentParser(std::string _appName) : _appName(_appName)
    {
        setName(_toolDoc, _appName);
        init();
    }
};

// ----------------------------------------------------------------------------
// Function hasOption()
// ----------------------------------------------------------------------------

/**
.Function.hasOption:
..summary:Returns whether a certain option is registered in the parser.
..cat:Miscellaneous
..signature:hasOption(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the option.
..returns:$true$ if the option is registered.
..include:seqan/arg_parse.h
*/

inline bool
hasOption(ArgumentParser const & me, std::string const & _name)
{
    return (hasKey(me.shortNameMap, _name) || hasKey(me.longNameMap, _name));
}

// ----------------------------------------------------------------------------
// Function addOption()
// ----------------------------------------------------------------------------

/**
.Function.addOption
..summary:Adds a @Class.ArgParseOption@ object to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addOption(parser, option)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The new @Class.ArgParseOption@ object that should be added.
...type:Class.ArgParseOption
..include:seqan/arg_parse.h
*/

inline void
addOption(ArgumentParser & me, ArgParseOption const & opt)
{
    // check if an option with the same identifiers was already registered
    SEQAN_CHECK(!hasOption(me, opt.shortName), "There already is an option with the name %s!", toCString(opt.shortName));
    SEQAN_CHECK(!hasOption(me, opt.longName), "There already is an option with the name %s!", toCString(opt.longName));

    // finally append the option
    appendValue(me.optionMap, opt);

    if (!empty(opt.shortName))
        me.shortNameMap.insert(std::make_pair<std::string, ArgumentParser::TOptionMapSize>(opt.shortName, length(me.optionMap) - 1));
    if (!empty(opt.longName))
        me.longNameMap.insert(std::make_pair<std::string, ArgumentParser::TOptionMapSize>(opt.longName, length(me.optionMap) - 1));
}

/**
.Function.addArgument
..summary:Adds a @Class.ArgParseArgument@ object to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addArgument(parser, argument)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The new @Class.ArgParseArgument@ object that should be added.
...type:Class.ArgParseArgument
..include:seqan/arg_parse.h
*/

inline void
addArgument(ArgumentParser & me, ArgParseArgument const & arg)
{
    // check previous arguments
    //  .. lists can only be last argument
    if(!me.argumentList.empty())
    {
        SEQAN_CHECK(!isListArgument(me.argumentList[me.argumentList.size() - 1]),
                    "You cannot add an additional argument after a list argument.");
    }

    // check current argument
    //  .. arguments should not have default values
    SEQAN_CHECK(arg.defaultValue.empty(), "Arguments cannot have default values.");
    SEQAN_CHECK(arg._numberOfArguments == 1, "n-Tuple of arguments are not supported.");

    me.argumentList.push_back(arg);
}

// ----------------------------------------------------------------------------
// Function addLine()
// ----------------------------------------------------------------------------

/**
.Function.addLine:
..summary:Adds a line of text to the help output of the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addLine(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A line of text that will be added to the help output.
...type:Shortcut.CharString
..include:seqan/arg_parse.h
*/

template <typename TString>
inline void
addLine(ArgumentParser & me, TString const & line)
{
    addOption(me, ArgParseOption("", "", line));
}

// ----------------------------------------------------------------------------
// Function addHelpLine()
// ----------------------------------------------------------------------------

/**
.Function.addHelpLine:
..summary:Adds an extra line of text below the help text of an option.
..cat:Miscellaneous
..signature:addHelpLine(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A line of text that will be added below the help text of an option.
...type:Shortcut.CharString
..include:seqan/arg_parse.h
*/

template <typename TString>
inline void
addHelpLine(ArgumentParser & me, TString const & line)
{
    addOption(me, ArgParseOption("", "", line));
}

// ----------------------------------------------------------------------------
// Function addSection()
// ----------------------------------------------------------------------------

/**
.Function.addSection:
..summary:Adds a new section the help output of the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addSection(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A section header that will be added to the help output.
...type:Shortcut.CharString
..include:seqan/arg_parse.h
*/

template <typename TString>
inline void
addSection(ArgumentParser & me, TString const & line)
{
    addLine(me, "");
    addLine(me, line);
}

// ----------------------------------------------------------------------------
// Function addTitleLine()
// ----------------------------------------------------------------------------

/**
.Function.addTitleLine:
..summary:Adds a line of text to the title output of the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addTitleLine(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A text line that will be added to the title output.
...type:Shortcut.CharString
..include:seqan/arg_parse.h
*/

template <typename TString>
inline void
addTitleLine(ArgumentParser & me, TString const & line)
{
    me._titleText.push_back(line);
}

// ----------------------------------------------------------------------------
// Function addVersionLine()
// ----------------------------------------------------------------------------

/**
.Function.addVersionLine:
..summary:Adds a line of text to the version output of the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addVersionLine(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A text line that will be added to the version output.
...type:Shortcut.CharString
..include:seqan/arg_parse.h
*/

template <typename TString>
inline void
addVersionLine(ArgumentParser & me, TString const & line)
{
    if (empty(me._versionText))
        addOption(me, ArgParseOption("V", "version", "Print version information."));
    me._versionText.push_back(line);
}

// ----------------------------------------------------------------------------
// Function addUsageLine()
// ----------------------------------------------------------------------------

/**
.Function.addUsageLine:
..summary:Adds a line of text to the usage output of the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addUsageLine(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:A text line that will be added to the usage output.
..include:seqan/arg_parse.h
*/

inline void
addUsageLine(ArgumentParser & me, std::string const & line)
{
    me._usageText.push_back(line);
}

// ----------------------------------------------------------------------------
// Function _getOptionIndex()
// ----------------------------------------------------------------------------
// note that it is assumed that the option exists if this method is called

inline Size< String<ArgParseOption> >::Type
_getOptionIndex(ArgumentParser const & me, std::string const & _name)
{
    typedef Size< String<ArgParseOption> >::Type TOptionPosition;
    TOptionPosition option_index;
    if (me.shortNameMap.find(_name) != me.shortNameMap.end())
    {
        option_index = me.shortNameMap.find(_name)->second;
    }
    else
    {
        option_index = me.longNameMap.find(_name)->second;
    }
    return option_index;
}

// ----------------------------------------------------------------------------
// Function getOption()
// ----------------------------------------------------------------------------

inline ArgParseOption &
getOption(ArgumentParser & me, std::string const & _name)
{
    SEQAN_CHECK(hasOption(me, _name), "Unknown option: %s", toCString(_name));
    return me.optionMap[_getOptionIndex(me,_name)];
}

inline ArgParseOption const &
getOption(ArgumentParser const & me, std::string const & _name)
{
    SEQAN_CHECK(hasOption(me, _name), "Unknown option: %s", toCString(_name));
    return me.optionMap[_getOptionIndex(me,_name)];
}

// ----------------------------------------------------------------------------
// Function setRequired()
// ----------------------------------------------------------------------------

/**
.Function.setRequired
..summary:Sets whether or not the option defined by the parameter $name$ (which can be
 either the short or the long name) is mandatory.
..cat:Miscellaneous
..signature:setRequired(parser, optionName, required)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.required:The new required value of the option.
...type:Bool
..include:seqan/arg_parse.h
*/

inline void
setRequired(ArgumentParser & me, std::string const & _name, bool required)
{
    SEQAN_CHECK(hasOption(me, _name), "Unknown option: %s", toCString(_name));
    return setRequired(getOption(me, _name), required);
}

// ----------------------------------------------------------------------------
// Function getArgument()
// ----------------------------------------------------------------------------

inline ArgParseArgument &
getArgument(ArgumentParser & me, unsigned position)
{
    SEQAN_CHECK(position < me.argumentList.size(),
                "ArgumentParser: Only %d arguments available", me.argumentList.size());
    return me.argumentList[position];
}

inline ArgParseArgument const &
getArgument(ArgumentParser const & me, unsigned position)
{
    SEQAN_CHECK(position < me.argumentList.size(),
                "ArgumentParser: Only %d arguments available", me.argumentList.size());
    return me.argumentList[position];
}


// ----------------------------------------------------------------------------
// Function _printStringSet()
// ----------------------------------------------------------------------------

template <typename TStringSet, typename TStream>
inline void
_printStringSet(TStringSet const & set, TStream & target)
{
    for (unsigned r = 0; r < length(set); ++r)
    {
        _streamWrite(target, set[r]);
        _streamPut(target, '\n');
    }
}

template <typename TStream>
inline void
_printStringSet(std::vector<std::string> const & set, TStream & target)
{
    for (unsigned r = 0; r < length(set); ++r)
    {
        target << set[r] << "\n";
    }
}


// ----------------------------------------------------------------------------
// Function _printUsage()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
_printUsage(ArgumentParser const & me, TStream & target)
{
    _streamWrite(target, "Usage: ");
    if (empty(me._usageText))
    {
        target << me._appName;
        _streamWrite(target, " [OPTION]... ");

        for (unsigned r = 0; r < me.argumentList.size(); ++r)
        {
            target << getArgumentLabel(getArgument(me, r));
        }
        _streamPut(target, '\n');
    }
    else
    {
        for (unsigned r = 0; r < length(me._usageText); ++r)
        {
            if (r)
                _streamWrite(target, "       ");
            target << me._appName;
            _streamPut(target, ' ');
            target << me._usageText[r];
            _streamPut(target, '\n');
        }
    }
}

// ----------------------------------------------------------------------------
// Function _printTitle()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
_printTitle(ArgumentParser const & me, TStream & target)
{
    _printStringSet(me._titleText, target);
}

// ----------------------------------------------------------------------------
// Function printShortHelp()
// ----------------------------------------------------------------------------

/**
.Function.printShortHelp
..summary:Prints a short help message for the parser to a stream
..cat:Miscellaneous
..signature:printShortHelp(parser[, stream])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.stream:Target stream (e.g. $std::cerr$).
..include:seqan/arg_parse.h
*/

template <typename TStream>
inline void
printShortHelp(ArgumentParser const & me, TStream & target)
{
    _printTitle(me, target);
    _printUsage(me, target);
    _streamWrite(target, "Try '");
    target << me._appName;
    _streamWrite(target, " --help' for more information.\n");
}

// ----------------------------------------------------------------------------
// Function printHelp()
// ----------------------------------------------------------------------------

/**
.Function.printHelp
..summary:Prints the complete help message for the parser to a stream.
..cat:Miscellaneous
..signature:printHelp(parser[, stream][, format])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
..param.format:Format to print, one of "html", "man", "txt".
..include:seqan/arg_parse.h
*/

template <typename TStream>
inline void
printHelp(ArgumentParser const & me, TStream & target, CharString const & format)
{
    ToolDoc toolDoc(me._toolDoc);
    clearEntries(toolDoc);  // We will append me._toolDoc later.

    // Build synopsis section.
    addSection(toolDoc, "Synopsis");
    for (unsigned i = 0; i < length(me._usageText); ++i)
    {
        std::string text = "\\fB";
        append(text, me._appName);
        append(text, "\\fP ");
        append(text, me._usageText[i]);
        addText(toolDoc, text, false);
    }

    // Build Description section, include options.
    addSection(toolDoc, "Description");

    // Add description to tool documentation.
    for (unsigned i = 0; i < length(me._description); ++i)
        addText(toolDoc, me._description[i]);

    // Add options to description section.
    for (unsigned i = 0; i < length(me.optionMap); ++i)
    {
        ArgParseOption const & opt = me.optionMap[i];
        if (empty(opt.shortName) && empty(opt.longName))  // this is not an option but a text line
        {
            if (empty(opt._helpText))  // TODO(holtgrew): Should go away in future.
                continue;  // Skip empty lines.

            // Is command line parser section, maps to ToolDoc subsection.
            std::string title = opt._helpText;
            append(title, ":");
            addSubSection(toolDoc, title);
        }
        else
        {
            // Build list item term.
            std::string term;
            if (!empty(opt.shortName))
            {
                term = "\\fB-";
                append(term, opt.shortName);
                append(term, "\\fP");
            }
            if (!empty(opt.shortName) && !empty(opt.longName))
                append(term, ", ");
            if (!empty(opt.longName))
            {
                append(term, "\\fB--");
                append(term, opt.longName);
                append(term, "\\fP");
            }
            // Get arguments, autogenerate if necessary.
            std::string arguments = getArgumentLabel(opt);

            // Write arguments to term line -> only exception, boolean flags
            if (!empty(arguments))
            {
                // Tokenize argument names.
                std::istringstream iss(toCString(arguments));
                std::vector<std::string> tokens;
                std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                          std::back_inserter<std::vector<std::string> >(tokens));
                // Append them, formatted in italic.
                for (unsigned i = 0; i < length(tokens); ++i)
                {
                    append(term, " \\fI");
                    append(term, tokens[i]);
                    append(term, "\\fP");
                }
            }

            // Add list item.
            addListItem(toolDoc, term, opt._helpText);
        }
    }

    append(toolDoc, me._toolDoc);
    print(target, toolDoc, format);
}

template <typename TStream>
inline void
printHelp(ArgumentParser const & me, TStream & target)
{
    printHelp(me, target, "txt");
}

// ----------------------------------------------------------------------------
// Function printVersion()
// ----------------------------------------------------------------------------

/**
.Function.printVersion
..summary:Prints a version text to a stream.
..cat:Miscellaneous
..signature:version(parser[, stream])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
..include:seqan/arg_parse.h
*/

template <typename TStream>
inline void
printVersion(ArgumentParser const & me, TStream & target)
{
    _printStringSet(me._versionText, target);
}

// ----------------------------------------------------------------------------
// Function isSet()
// ----------------------------------------------------------------------------

/**
.Function.isSet
..summary:Returns whether an option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSet(parser,optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that identifies the option (either short or long name).
..returns:$true$ if the option was set.
..include:seqan/arg_parse.h
*/

inline bool
isSet(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return isSet(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function _allRequiredSet()
// ----------------------------------------------------------------------------

inline bool
_allRequiredSet(ArgumentParser const & me)
{
    for (unsigned o = 0; o < length(me.optionMap); ++o)
        if (!isSet(me.optionMap[o]) && isRequired(me.optionMap[o]))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function _allArgumentsSet()
// ----------------------------------------------------------------------------

inline bool
_allArgumentsSet(ArgumentParser const & me)
{
    for(unsigned a = 0; a < me.argumentList.size(); ++a)
    {
        if(!isSet(me.argumentList[a]))
            return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Function _parseAppName()
// ----------------------------------------------------------------------------

inline std::string
_parseAppName(std::string const & candidate)
{
    //IOREV _notio_ irrelevant for io-revision
    int i = length(candidate) - 1;

    for (; i >= 0; --i)
        if (candidate[i] == '\\' || candidate[i] == '/')
            break;
    return candidate.substr(i+1);
}

// ----------------------------------------------------------------------------
// Function getOptionValue()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValue:
..summary:Retrieves the value of an option given either the short or long name.
..cat:Miscellaneous
..signature:getOptionValue(value, parser, optionIdentifier[, argNo])
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A @Shortcut.std::string@ that is either the short or long name of the option.
..param.argNo:If the option is list, the $argNo$-th list element is returned.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline bool
getOptionValue(TValue & val, ArgumentParser const & me, std::string const & name,
               unsigned argNo)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return _convertOptionValue(val, getOption(me, name), getArgumentValue(getOption(me, name), argNo));
}

template <typename TValue>
inline bool
getOptionValue(TValue & val, ArgumentParser const & me, std::string const & name)
{
    return getOptionValue(val, me, name, 0);
}

// ----------------------------------------------------------------------------
// Function getArgumentValue()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValue:
..summary:Retrieves the value of an argument given by its position.
..cat:Miscellaneous
..signature:getArgumentValue(value, parser, argumentPosition[, argNo])
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..param.argNo:If the argument is a list, the $argNo$-th list element is returned.
..returns: $true$ if the requested argument is set and has the requested type, $false$ otherwise.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline bool
getArgumentValue(TValue & value, ArgumentParser & me, unsigned argumentPosition, unsigned argNo)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    return _convertArgumentValue(value, getArgument(me, argumentPosition), getArgumentValue(getArgument(me, argumentPosition), argNo));
}

template <typename TValue>
inline bool
getArgumentValue(TValue & value, ArgumentParser & me, unsigned argumentPosition)
{
    return getArgumentValue(value, me, argumentPosition, 0);
}

// ----------------------------------------------------------------------------
// Function getOptionValues()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValues
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getOptionValues(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A @Shortcut.std::string@ that is either the short or long name of the option.
..returns: A $String<std::string>$ of option values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const &
getOptionValues(ArgumentParser & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return getArgumentValues(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function getArgumentValues()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValues
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getArgumentValues(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: A $String<std::string>$ of argument values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const &
getOptionValues(ArgumentParser & me, unsigned argumentPosition)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition, "Argument Parser has only %d arguments.", me.argumentList.size());
    return getArgumentValues(getArgument(me, argumentPosition));
}


// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/**
.Function.setMinValue
..summary:Sets the minimum value of a @Class.ArgParseOption@ object identified by .
..cat:Miscellaneous
..signature:setMinValue(parser,optionName,minValue)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The identifier of the command line option.
..param.minValue:A @Shortcut.std::string@ containing a string representation of the minimum value of the @Class.ArgParseOption@.
..include:seqan/arg_parse.h
*/
inline void
setMinValue(ArgumentParser & me, std::string const & name, std::string const & _minValue)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMinValue(getOption(me, name), _minValue);
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/**
.Function.setMaxValue
..summary:Sets the maximum value of a @Class.ArgParseOption@ object.
..cat:Miscellaneous
..signature:setMaxValue(parser,optionName,maxValue)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The identifier of the command line option.
..param.maxValue:A @Shortcut.std::string@ containing a string representation of the maximum value of the @Class.ArgParseOption@.
..include:seqan/arg_parse.h
*/

inline void
setMaxValue(ArgumentParser & me, std::string const & name,
            std::string const & _maxValue)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMaxValue(getOption(me, name), _maxValue);
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/**
.Function.setValidValues
..summary:Sets the set of allowed values of a @Class.ArgParseOption@ object.
..cat:Miscellaneous
..signature:setValidValues(parser,optionName,values)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The identifier of the command line option.
..param.values:A $String<std::string>$ containing all valid entries for the option.
..include:seqan/arg_parse.h
*/

inline void
setValidValues(ArgumentParser & me, std::string const & name,
               std::vector<std::string> const & _values)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), _values);
}

inline void
setValidValues(ArgumentParser & me, std::string const & name,
               std::string const & _values)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), _values);
}

// ----------------------------------------------------------------------------
// Function addDescription()
// ----------------------------------------------------------------------------

/**
.Function.addDescription
..summary:Appends a description paragraph to the @Class.ArgumentParser@ documentation.
..cat:Miscellaneous
..signature:addDescription(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:The description paragraph.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void addDescription(ArgumentParser & me, std::string const & description)
{
    appendValue(me._description, description);
}

// ----------------------------------------------------------------------------
// Function setAppName()
// ----------------------------------------------------------------------------

/**
.Function.setAppName
..summary:Sets application name of @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:setAppName(parser, appName)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.appName:The name of the application.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void setAppName(ArgumentParser & me, std::string const & name)
{
    setName(me._toolDoc, name);
}

// ----------------------------------------------------------------------------
// Function setShortDescription()
// ----------------------------------------------------------------------------

/**
.Function.setShortDescription
..summary:Sets short description test of @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:setShortDescription(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:The short description text.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void setShortDescription(ArgumentParser & me, std::string const & description)
{
    setShortDescription(me._toolDoc, description);
}

// ----------------------------------------------------------------------------
// Function setVersion()
// ----------------------------------------------------------------------------

/**
.Function.setVersion
..summary:Sets version string of @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:setVersion(parser, versionString)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.versionString:The version string to set.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void setVersion(ArgumentParser & me, std::string const & versionString)
{
    setVersion(me._toolDoc, versionString);
}

// ----------------------------------------------------------------------------
// Function setDate()
// ----------------------------------------------------------------------------

/**
.Function.setDate
..summary:Sets date string of @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:setDate(parser, date)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.date:The date string.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void setDate(ArgumentParser & me, std::string const & date)
{
    setDate(me._toolDoc, date);
}

// ----------------------------------------------------------------------------
// Function addTextSection()
// ----------------------------------------------------------------------------

/**
.Function.addTextSection
..summary:Adds a text section to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addTextSection(parser, title)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.title:The section title.
..returns:$void$
..remarks:This will result in an additional section heading to be printed.
..include:seqan/arg_parse.h
*/

inline void addTextSection(ArgumentParser & me, std::string const & title)
{
    addSection(me._toolDoc, title);
}

// ----------------------------------------------------------------------------
// Function addTextSubSection()
// ----------------------------------------------------------------------------

/**
.Function.addTextSubSection
..summary:Adds a text subsection to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addTextSubSection(parser, title)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.title:The subsection title.
..returns:$void$
..remarks:This will result in an additional subsection heading to be printed.
..include:seqan/arg_parse.h
*/

inline void addTextSubSection(ArgumentParser & me, std::string const & title)
{
    addSubSection(me._toolDoc, title);
}

// ----------------------------------------------------------------------------
// Function addText()
// ----------------------------------------------------------------------------

/**
.Function.addText
..summary:Appends a text paragraph to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addText(parser, text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.text:The content of the text.
..returns:$void$
..include:seqan/arg_parse.h
*/

inline void addText(ArgumentParser & me, std::string const & text)
{
    addText(me._toolDoc, text);
}

}  // namespace seqan

#endif // SANDBOX_ARG_PARSE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_
