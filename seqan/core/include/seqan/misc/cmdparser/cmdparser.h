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

#ifndef CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_H_
#define CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_H_

#include <sstream>
#include <seqan/map.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/misc/cmdparser/cmdoption.h>
#include <seqan/misc/cmdparser/cmdparser_type_support.h>
#include <seqan/misc/cmdparser/cmdparser_doc.h>

#include <seqan/misc/misc_terminal.h>

#include <iostream>
#include <fstream>

namespace seqan {

/**
 * TODO: arguments also need a type if they are used as input/output file
 * TODO: correct parameter ordering of nearly all cmdparser function to be seqan conform f(out,in)
 */

// ---------------------------------------------------------------------------
// Class EnvVarDescription
// ---------------------------------------------------------------------------

/**
.Class.EnvVarDescription
..cat:Miscellaneous
..summary:Description of environment variable.
..signature:EnvVarDescription
..remarks:Don't use directly, use @Function.describeEnvVar@ instead.
..see:Function.describeEnvVar
..include:seqan/misc/misc_cmdparser.h

.Memfunc.EnvVarDescription#EnvVarDescription
..class:Class.CommandLineParser
..summary:Constructor
..signature:CommandLineParser(name, text)
..param.name:Name of the environment variable.
...type:Shortcut.CharString
..param.text:Description of the environment variable.
...type:Shortcut.CharString
*/

class EnvVarDescription
{
public:
    CharString _name;
    CharString _text;

    EnvVarDescription(CharString const & name, CharString const & doc) :
            _name(name), _text(doc)
    {}
};

/**
.Class.CommandLineParser
..cat:Miscellaneous
..summary:Stores multiple @Class.CommandLineOption@ objects and parses the command line _arguments for these options.
..signature:CommandLineParser
..include:seqan/misc/misc_cmdparser.h
*/

/**
.Memfunc.CommandLineParser#CommandLineParser
..class:Class.CommandLineParser
..summary:Constructor
..signature:CommandLineParser ()
..signature:CommandLineParser (applicationName)
..param.applicationName:A @Shortcut.CharString@ containing the name of the application.
..remarks:If the name of the application is not passed to the constructor it will be extracted from the command line.
*/

class CommandLineParser
{
public:
    // ----------------------------------------------------------------------------
    // Class Typedefs
    // ----------------------------------------------------------------------------
    typedef String<CommandLineOption>   TOptionMap;
    typedef Size<TOptionMap>::Type      TSize;

    typedef std::map<CharString, TSize> TStringMap;
    typedef String<CharString>          TValueMap;

    // ----------------------------------------------------------------------------
    // Mapping of option names to options
    // ----------------------------------------------------------------------------
    TStringMap shortNameMap;
    TStringMap longNameMap;
    TOptionMap optionMap;
    
    String<EnvVarDescription> envVars;

    // ----------------------------------------------------------------------------
    // Members
    // ----------------------------------------------------------------------------
    unsigned           _requiredArguments;
    String<CharString> _arguments;
    CharString         _appName;
    String<CharString> _titleText;
    String<CharString> _usageText;
    String<CharString> _versionText;

    // ----------------------------------------------------------------------------
    // Command line formating members
    // ----------------------------------------------------------------------------

    // The following is an example for a command line options list:
    //
    //    -i,   --insert-foo   This is the description of the insert-foo parameters.
    //                         Its description is pretty long.
    //          --long-only    This is a parameter with a long variant only.
    //    -vvl, --very-very-long-param
    //                         This is a very very long parameter.
    //     `-'  `----------'
    //      A        B
    //    `----------------'
    // `-'       C          `-'
    //  D                    E
    // `----------------------------------------------------------------------------'
    //                                      F
    //
    // We have the following measures:
    //
    // A: Longest short option's length, determined from options given by user.
    // B: Longest long options' length, determined from options given by user.
    // C: Left column width, determined by A and B, there is an upper limit.
    // D: Padding on left-hand side, configuration.
    // E: Center padding, configuration.
    // F: Maximal screen width, given by terminal window width, configuration.
    //
    // Additionally, the configuration contains a minimal screen width, and a
    // default screen width.

    unsigned _paddingLeft;
    unsigned _paddingCenter;
    unsigned _paddingRight;
    unsigned _maximalLeftColumnWidth;
    unsigned _minimalScreenWidth;
    unsigned _maximalScreenWidth;
    unsigned _defaultScreenWidth;

    // ----------------------------------------------------------------------------
    // return values for unset parameters
    // ----------------------------------------------------------------------------
    const CharString         _null;
    const String<CharString> _nullSet;

    // friend declaration to make addOption() available in init function
    friend inline void addOption(CommandLineParser & me, CommandLineOption const & opt);

    // ----------------------------------------------------------------------------
    // Function init()
    // ----------------------------------------------------------------------------
    void init()
    {
        // Set defaults for the layout.
        _paddingLeft = 2;
        _paddingRight = 1;
        _paddingCenter = 4;
        _maximalLeftColumnWidth = 25;
        _minimalScreenWidth = 40;
        _maximalScreenWidth = 200;
        _defaultScreenWidth = 80;

        _requiredArguments = 0;
        addOption(*this, CommandLineOption("h", "help", "Displays this help message.", OptionType::Boolean));
        addOption(*this, CommandLineOption("", "write-ctd", "Exports the app's interface description to a .ctd file.", OptionType::OUTPUTFILE));
    }

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    CommandLineParser()       
    {
        init();
    }

    CommandLineParser(CharString _appName) : _appName(_appName)
    {
        init();
    }
};

// ----------------------------------------------------------------------------
// Function hasOptionLong()
// ----------------------------------------------------------------------------

/**
.Function.hasOptionLong:
..summary:Returns whether a certain long-name option is registered in the parser.
..cat:Miscellaneous
..signature:hasOptionLong(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the long-name option.
..returns:$true$ if the option is registered.
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
hasOptionLong(CommandLineParser const & me, CharString const & _long)
{
    return hasKey(me.longNameMap, _long);
}

// ----------------------------------------------------------------------------
// Function hasOptionShort()
// ----------------------------------------------------------------------------

/**
.Function.hasOptionShort:
..summary:Returns whether a certain short-name option is registered in the parser.
..cat:Miscellaneous
..signature:hasOptionShort(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the short-name option.
..returns:$true$ if the option is registered.
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
hasOptionShort(CommandLineParser const & me, CharString const & _short)
{
    return hasKey(me.shortNameMap, _short);
}

// ----------------------------------------------------------------------------
// Function hasOption()
// ----------------------------------------------------------------------------

/**
.Function.hasOption:
..summary:Returns whether a certain option is registered in the parser.
..cat:Miscellaneous
..signature:hasOption(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the option.
..returns:$true$ if the option is registered.
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
hasOption(CommandLineParser const & me, CharString const & _name)
{
    return (hasKey(me.shortNameMap, _name) || hasKey(me.longNameMap, _name));
}

// ----------------------------------------------------------------------------
// Function describeEnvVar()
// ----------------------------------------------------------------------------

/**
.Function.describeEnvVar
..summary:Adds a @Class.EnvVarDescription@ object to the @Class.CommandLineParser@.
..cat:Miscellaneous
..signature:describeEnvVar(parser, name, description)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.name:The name of the environment variable.
...type:Shortcut.CharString
..param.name:The description for the environment variable.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

void describeEnvVar(CommandLineParser & parser, CharString const & name, CharString const & txt)
{
    appendValue(parser.envVars, EnvVarDescription(name, txt));
}

// ----------------------------------------------------------------------------
// Function addOption()
// ----------------------------------------------------------------------------

/**
.Function.addOption
..summary:Adds a @Class.CommandLineOption@ object to the @Class.CommandLineParser@.
..cat:Miscellaneous
..signature:addOption(parser, option)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.option:The new @Class.CommandLineOption@ object that should be added.
...type:Class.CommandLineOption
..include:seqan/misc/misc_cmdparser.h
*/

inline void
addOption(CommandLineParser & me, CommandLineOption const & opt)
{
    // check if an option with the same identifiers was already registered
    SEQAN_CHECK(!hasOption(me, opt.shortName), "There already is an option with the name %s!", toCString(opt.shortName));
    SEQAN_CHECK(!hasOption(me, opt.longName), "There already is an option with the name %s!", toCString(opt.longName));

    // finally append the option
    appendValue(me.optionMap, opt);

    if (!empty(opt.shortName))
        insert(me.shortNameMap, opt.shortName, length(me.optionMap) - 1);
    if (!empty(opt.longName))
        insert(me.longNameMap, opt.longName, length(me.optionMap) - 1);
}

// ----------------------------------------------------------------------------
// Function addLine()
// ----------------------------------------------------------------------------

/**
.Function.addLine:
..summary:Adds a line of text to the help output of the @Class.CommandLineParser@.
..cat:Miscellaneous
..signature:addLine(parser, text)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.text:A line of text that will be added to the help output.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TString>
inline void
addLine(CommandLineParser & me, TString const & line)
{
    addOption(me, CommandLineOption("", "", line, 0));
}

// ----------------------------------------------------------------------------
// Function addHelpLine()
// ----------------------------------------------------------------------------

/**
.Function.addHelpLine:
..summary:Adds an extra line of text below the help text of an option.
..cat:Miscellaneous
..signature:addHelpLine(parser, text)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.text:A line of text that will be added below the help text of an option.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TString>
inline void
addHelpLine(CommandLineParser & me, TString const & line)
{
    addOption(me, CommandLineOption("", "", line, 1));
}

// ----------------------------------------------------------------------------
// Function addSection()
// ----------------------------------------------------------------------------

/**
.Function.addSection:
..summary:Adds a new section the help output of the @Class.CommandLineParser@.
..cat:Miscellaneous
..signature:addSection(parser, text)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.text:A section header that will be added to the help output.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TString>
inline void
addSection(CommandLineParser & me, TString const & line)
{
    addLine(me, "");
    addLine(me, line);
}

// ----------------------------------------------------------------------------
// Function addTitleLine()
// ----------------------------------------------------------------------------

/**
.Function.addTitleLine:
..summary:Adds a line of text to the title output of the @Class.CommandLineParser@.
..cat:Miscellaneous
..signature:addTitleLine(parser, text)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.text:A text line that will be added to the title output.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TString>
inline void
addTitleLine(CommandLineParser & me, TString const & line)
{
    appendValue(me._titleText, line);
}

// ----------------------------------------------------------------------------
// Function addVersionLine()
// ----------------------------------------------------------------------------

/**
.Function.addVersionLine:
..summary:Adds a line of text to the version output of the @Class.CommandLineParser@.
..cat:Miscellaneous
..signature:addVersionLine(parser, text)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.text:A text line that will be added to the version output.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TString>
inline void
addVersionLine(CommandLineParser & me, TString const & line)
{
    if (empty(me._versionText))
        addOption(me, CommandLineOption("V", "version", "Print version information.", OptionType::Boolean));
    appendValue(me._versionText, line);
}

// ----------------------------------------------------------------------------
// Function addUsageLine()
// ----------------------------------------------------------------------------

/**
.Function.addUsageLine:
..summary:Adds a line of text to the usage output of the @Class.CommandLineParser@.
..cat:Miscellaneous
..signature:addUsageLine(parser, text)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.text:A text line that will be added to the usage output.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

inline void
addUsageLine(CommandLineParser & me, CharString const & line)
{
    appendValue(me._usageText, line);
}

// ----------------------------------------------------------------------------
// Function _getOptionIndex()
// ----------------------------------------------------------------------------

inline Size< String<CommandLineOption> >::Type
_getOptionIndex(CommandLineParser const & me, CharString const & _name)
{
    typedef Size< String<CommandLineOption> >::Type TOptionPosition;
    TOptionPosition option_index;
    if (hasKey(me.shortNameMap, _name))
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

inline CommandLineOption &
getOption(CommandLineParser & me, CharString const & _name)
{
    SEQAN_ASSERT_MSG(hasOption(me, _name), "Unknown option: %s", toCString(_name));
    return me.optionMap[_getOptionIndex(me,_name)];
}

inline CommandLineOption const &
getOption(CommandLineParser const & me, CharString const & _name)
{
    SEQAN_ASSERT_MSG(hasOption(me, _name), "Unknown option: %s", toCString(_name));
    return me.optionMap[_getOptionIndex(me,_name)];
}

// ----------------------------------------------------------------------------
// Function setRequiredArguments()
// ----------------------------------------------------------------------------

/**
.Function.setRequiredArguments
..summary:Sets the number of _arguments (non-parameterized options) are required by the program.
..cat:Miscellaneous
..signature:requiredArguments(parser, count)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.count:A $unsigned int$ defining the amount of non-parameterized options requried by the program.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setRequiredArguments(CommandLineParser & me, unsigned count)
{
    me._requiredArguments = count;
}

/**
.Function.requiredArguments
..summary:Sets the number of _arguments (non-parameterized options) are required by the program.
..cat:Miscellaneous
..signature:requiredArguments(parser, count)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.count:A $unsigned int$ defining the amount of non-parameterized options requried by the program.
..status:deprecated, use $Function.setRequiredArguments$
..see:Function.setRequiredArguments
..include:seqan/misc/misc_cmdparser.h
*/

inline void
requiredArguments(CommandLineParser & me, unsigned count)
{
    setRequiredArguments(me, count);
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

// ----------------------------------------------------------------------------
// Function _printUsage()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
_printUsage(CommandLineParser const & me, TStream & target)
{
    _streamWrite(target, "Usage: ");
    if (empty(me._usageText))
    {
        _streamWrite(target, me._appName);
        _streamWrite(target, " [OPTION]... ");
        for (unsigned r = 0; r < me._requiredArguments; ++r)
        {
            _streamWrite(target, "<ARG");
            _streamPutInt(target, r + 1);
            _streamWrite(target, "> ");
        }
        _streamPut(target, '\n');
    }
    else
    {
        for (unsigned r = 0; r < length(me._usageText); ++r)
        {
            if (r)
                _streamWrite(target, "       ");
            _streamWrite(target, me._appName);
            _streamPut(target, ' ');
            _streamWrite(target, me._usageText[r]);
            _streamPut(target, '\n');
        }
    }
}

// ----------------------------------------------------------------------------
// Function _printTitle()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
_printTitle(CommandLineParser const & me, TStream & target)
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
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TStream>
inline void
printShortHelp(CommandLineParser const & me, TStream & target)
{
    _printTitle(me, target);
    _printUsage(me, target);
    _streamWrite(target, "Try '");
    _streamWrite(target, me._appName);
    _streamWrite(target, " --help' for more information.\n");
}


/**
.Function.shortHelp
..summary:Prints a short help message for the parser to a stream
..cat:Miscellaneous
..signature:shortHelp(parser[, stream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
..status:deprecated, use $Function.printShortHelp$
..see:Function.printShortHelp
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TStream>
inline void
shortHelp(CommandLineParser const & me, TStream & target)
{
    printShortHelp(me, target);
}

template <typename TStream>
inline void
shortHelp(CommandLineParser const & me)
{
    printShortHelp(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function printHelp()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Reorder?
template <typename TStream>
void printEnvVarDocumentation(TStream & stream, CommandLineParser const & me);

/**
.Function.printHelp
..summary:Prints the complete help message for the parser to a stream.
..cat:Miscellaneous
..signature:printHelp(parser[, stream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TStream>
inline void
printHelp(CommandLineParser const & me, TStream & target)
{
    _printTitle(me, target);
    _streamPut(target, '\n');
    _printUsage(me, target);
    _streamPut(target, '\n');

    // Compute length of longest short option.
    unsigned longestShort = 0;
    for (unsigned i = 0; i < length(me.optionMap); ++i)
    {
        CommandLineOption const & opt = me.optionMap[i];
        if (length(opt.shortName) > longestShort)
            longestShort = length(opt.shortName) + 1;
    }

    // Now, get longest long option that does not violate the maximal length of the left column.
    const unsigned COMMA_PADDING = 2;  // == length(", ")
    unsigned longestLong = 0;
    for (unsigned i = 0; i < length(me.optionMap); ++i)
    {
        CommandLineOption const & opt = me.optionMap[i];
        unsigned len = length(opt.longName) + 2;
        if (me._paddingLeft + longestShort + COMMA_PADDING + len > me._maximalLeftColumnWidth)
            continue;  // Ignore, will put description one line below.
        if (len > longestLong)
            longestLong = len;
    }

    // Compute virtual tabs to shift to.
    unsigned tabShort = me._paddingLeft;
    unsigned tabLong  = tabShort + COMMA_PADDING + longestShort;
    unsigned tabDescription = tabLong + me._paddingCenter + longestLong;

    // Get screen width from terminal size, default and minimal screen width.
    unsigned screenWidth = 0, screenHeight = 0;
    bool success = getTerminalSize(screenWidth, screenHeight);
    if (!success)
        screenWidth = me._defaultScreenWidth;
    if (screenWidth < me._minimalScreenWidth)
        screenWidth = me._minimalScreenWidth;
    if (me._maximalScreenWidth > 0u && screenWidth > me._maximalScreenWidth)
        screenWidth = me._maximalScreenWidth;
    screenWidth -= me._paddingRight;
    // std::cerr << "screen-width  == " << screenWidth << std::endl;
    // std::cerr << "tab-short     == " << tabShort << std::endl;
    // std::cerr << "tab-long      == " << tabLong << std::endl;
    // std::cerr << "tab-descript  == " << tabDescription << std::endl;
    // std::cerr << "longest-short == " << longestShort << std::endl;
    // std::cerr << "longest-long  == " << longestLong << std::endl;

    // Print the options.
    for (unsigned o = 0; o < length(me.optionMap); ++o)
    {
        unsigned s = 0;  // Column of current caret.
        CommandLineOption const & opt = me.optionMap[o];
        if (isHiddenOption(opt))
            continue;  // do not print hidden options
        if (opt.optionType == 0)
            continue;  // TODO(holtgrew): Can this happen, isn't this an error?

        // Print short and long option.
        for (; s < tabShort; ++s)
            _streamPut(target, ' ');
        if (!empty(opt.shortName))
        {
            _streamPut(target, '-');
            _streamWrite(target, opt.shortName);
            s += length(opt.shortName) + 1;
            if (!empty(opt.longName))
            {
                _streamWrite(target, ", ");
                s += 2;
            }
        }
        for (; s < tabLong; ++s)
            _streamPut(target, ' ');
        if (!empty(opt.longName))
        {
            _streamWrite(target, "--");
            _streamWrite(target, opt.longName);
            s += length(opt.longName) + 2;
        }

        // Start description on next line if necessary.
        if (s > me._maximalLeftColumnWidth)
        {
            _streamPut(target, '\n');
            s = 0;
        }
        for (; s < tabDescription; ++s)
            _streamPut(target, ' ');

        // Print description.

        // First, tokenize text into words.
        StringSet<CharString> words;
        CharString word;
        for (unsigned i = 0; i < length(opt.helpText); ++i)
        {
            char c = opt.helpText[i];
            if (isspace(c))
            {
                if (!empty(word))
                    appendValue(words, word);
                clear(word);
            }
            else
            {
                appendValue(word, c);
            }
        }
        if (!empty(word))
            appendValue(words, word);

        // Now, print all the words.
        bool firstWord = true;
        for (unsigned i = 0; i < length(words); ++i)
        {
            if (s + 1 + length(words[i]) <= screenWidth)
            {
                // Fits on this line.  First, print space if not first word.
                if (!firstWord)
                {
                    _streamPut(target, ' ');
                    s += 1;
                }
                firstWord = false;

                // Print word itself.
                _streamWrite(target, words[i]);
                s += length(words[i]);
            }
            else
            {
                // Does not fit on this line.  Print on next.
                _streamPut(target, '\n');
                for (s = 0; s < tabDescription; ++s)
                    _streamPut(target, ' ');
                _streamWrite(target, words[i]);
                s += length(words[i]);
            }
        }
        _streamPut(target, '\n');
    }

    if (!empty(me.envVars))
    {
        _streamWrite(target, "\nEnvironment Variables\n\n");
        printEnvVarDocumentation(target, me);
    }
}

/**
.Function.help
..summary:Prints the complete help message for the parser to a stream.
..cat:Miscellaneous
..signature:help(parser[, stream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
..status:deprecated, use $Function.printHelp$
..see:Function.printHelp
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TStream>
inline void
help(CommandLineParser const & me, TStream & target)
{
    printHelp(me, target);
}

inline void
help(CommandLineParser const & me)
{
    printHelp(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function printVersion()
// ----------------------------------------------------------------------------

/**
.Function.printVersion
..summary:Prints a version text to a stream.
..cat:Miscellaneous
..signature:version(parser[, stream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TStream>
inline void
printVersion(CommandLineParser const & me, TStream & target)
{
    _printStringSet(me._versionText, target);
}

/**
.Function.version
..summary:Prints a version text to a stream.
..cat:Miscellaneous
..signature:version(parser[, stream])
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.stream:Target stream (e.g. $std::cerr$).
...default: $std::cerr$
..status:deprecated, use $Function.printVersion$
..see:Function.printVersion
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TStream>
inline void
version(CommandLineParser const & me, TStream & target)
{
    printVersion(me, target);
}

inline void
version(CommandLineParser const & me)
{
    printVersion(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function isSet()
// ----------------------------------------------------------------------------

/**
.Function.isSet
..summary:Returns whether an option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSet(parser,optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the option (either short or long name).
..returns:$true$ if the option was set.
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isSet(CommandLineParser const & me, CharString const & name)
{
    SEQAN_ASSERT_MSG(hasOption(me, name), "Unknown option: %s", toCString(name));
    return !empty(getOption(me, name).value);
}

/**
.Function.isSetShort
..summary:Returns whether a short-name option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSetShort(parser,optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the short-name option.
..returns:$true$ if the option was set.
..status:deprecated, use $Function.isSet$
..see:Function.isSet
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isSetShort(CommandLineParser & me, CharString const & shortName)
{
    return isSet(me, shortName);
}

/**
.Function.isSetLong
..summary:Returns whether a long-name option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSetLong(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the long-name option.
..returns:$true$ if the option was set.
..status:deprecated, use $Function.isSet$
..see:Function.isSet
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isSetLong(CommandLineParser & me, CharString const & longName)
{
    return isSet(me, longName);
}

// ----------------------------------------------------------------------------
// Function _allMandatorySet()
// ----------------------------------------------------------------------------

inline bool
_allMandatorySet(CommandLineParser const & me)
{
    for (unsigned o = 0; o < length(me.optionMap); ++o)
        if (empty(me.optionMap[o].value) && isOptionMandatory(me.optionMap[o]))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function _parseAppName()
// ----------------------------------------------------------------------------

inline CharString
_parseAppName(CharString const & candidate)
{
    //IOREV _notio_ irrelevant for io-revision
    int i = length(candidate) - 1;

    for (; i >= 0; --i)
        if (candidate[i] == '\\' || candidate[i] == '/')
            break;

    return suffix(candidate, i + 1);
}

// ----------------------------------------------------------------------------
// Function _reportInvalidType()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
inline void
_reportInvalidType(CommandLineParser const & me, CommandLineOption const & opt,
                   CharString const & val, TErrorStream & estream)
{
    _streamWrite(estream, me._appName);
    _streamWrite(estream, ": \"");
    _streamWrite(estream, val);
    _streamWrite(estream, "\" is not a valid ");

    // there should be no other situation then those two
    if (isIntOption(opt))
        _streamWrite(estream, "integer");
    else if (isDoubleOption(opt))
        _streamWrite(estream, "double");

    _streamWrite(estream, " value for '");
    _writeOptName(estream, opt);
    _streamWrite(estream, "'\n");
}

// ----------------------------------------------------------------------------
// Function _reportMissingArguments()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
inline void
_reportMissingArgument(CommandLineParser const & me, 
                       CommandLineOption const & opt, TErrorStream & estream)
{
    _streamWrite(estream, me._appName);
    _streamWrite(estream, ": \'");
    _writeOptName(estream, opt);
    _streamWrite(estream, "\' requires ");
    _streamPutInt(estream, opt.argumentsPerOption);
    _streamWrite(estream, " value(s)\n");
}

// ----------------------------------------------------------------------------
// Function _reportInvalidOption()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
inline void
_reportInvalidOption(CommandLineParser const & me, CharString const & option, 
                     TErrorStream & estream)
{
    _streamWrite(estream, me._appName);
    _streamWrite(estream, ": invalid option '");
    _streamWrite(estream, option);
    _streamWrite(estream, "\'\n");
}

// ----------------------------------------------------------------------------
// Function _reportValueNotInRange()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
inline void
_reportValueNotInRange(CommandLineOption const & opt, CharString const & val, 
                       TErrorStream & estream)
{
    _writeOptName(estream, opt);
    _streamWrite(estream, ": given argument \"");
    _streamWrite(estream, val);
    _streamWrite(estream, "\" is not in the required range [");
    _streamWrite(estream, (opt.minValue != "" ? opt.minValue : "-inf"));
    _streamWrite(estream, ":");
    _streamWrite(estream, (opt.maxValue != "" ? opt.maxValue : "+inf"));
    _streamWrite(estream, "]\n");
}

// ----------------------------------------------------------------------------
// Function _reportInvalidValue()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
inline void
_reportInvalidValue(CommandLineOption const & opt, CharString const & val, 
                    TErrorStream & estream)
{
    typedef Iterator<StringSet<CharString> const, Rooted>::Type TStringSetIter;

    _writeOptName(estream, opt);
    _streamWrite(estream, ": given argument \"");
    _streamWrite(estream, val);
    _streamWrite(estream, "\" is not a valid value [");
    for (TStringSetIter valid = begin(opt.validValues);; )
    {
        _streamWrite(estream, *valid);

        goNext(valid);
        if (valid == end(opt.validValues))
            break;
        else
            _streamWrite(estream, ", ");
    }
    _streamWrite(estream, "]\n");
}

// ----------------------------------------------------------------------------
// Function _reportInvalidFileType()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
inline void
_reportInvalidFileType(CommandLineOption const & opt, CharString const & val, 
                       TErrorStream & estream)
{
    _writeOptName(estream, opt);
    _streamWrite(estream, ": given argument \"");
    _streamWrite(estream, val);
    _streamWrite(estream, "\" is not a valid file type [");
    for (Iterator<StringSet<CharString> const, Rooted>::Type valid = begin(opt.validValues);; )
    {
        _streamWrite(estream, *valid);

        goNext(valid);
        if (valid == end(opt.validValues))
            break;
        else
            _streamWrite(estream, ", ");
    }
    _streamWrite(estream, "]\n");
}

// ----------------------------------------------------------------------------
// Function _checkMinMaxValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TErrorStream>
bool _checkMinMaxValue(CommandLineOption const & opt, CharString const & val, 
                       TErrorStream & estream)
{
    TValue d_value = 0;
    if (!_convertOptionValue(opt, d_value, val)) 
        SEQAN_FAIL("Conversion should work");

    if (opt.minValue != "")
    {
        // check min and max
        TValue minVal = 0;
        _convertOptionValue(opt, minVal, opt.minValue);

        if (d_value < minVal)
        {
            _reportValueNotInRange(opt, val, estream);
            return false;
        }
    }

    if (opt.maxValue != "")
    {
        TValue maxVal = 0;
        _convertOptionValue(opt, maxVal, opt.maxValue);

        if (d_value > maxVal)
        {
            _reportValueNotInRange(opt, val, estream);
            return false;
        }
    }

    return true;
    
}

// ----------------------------------------------------------------------------
// Function _checkValidValue()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
bool _checkValidValues(CommandLineOption const & opt, CharString const & val, 
                       TErrorStream & estream)
{
    typedef Iterator<StringSet<CharString> const, Rooted>::Type TStringSetIter;

    if (length(opt.validValues) == 0) return true; // no restrictions

    if (isInputFile(opt) || isOutputFile(opt))
    {
        // check if our val is in the valid strings
        for (TStringSetIter valid = begin(opt.validValues); valid != end(opt.validValues); goNext(valid))
        {
            if (length(*valid) > length(val))
                continue;
            else if (suffix(val, length(val) - length(*valid)) == *valid)
                return true;
        }
        // couldn't find our target
        _reportInvalidFileType(opt, val, estream);
        return false;
    }
    else
    {
        // check if our val is in the valid strings
        for (TStringSetIter valid = begin(opt.validValues); valid != end(opt.validValues); goNext(valid))
        {
            if (*valid == val)
                return true;
        }
        // couldn't find our target
        _reportInvalidValue(opt, val, estream);
        return false;
    }
}

// ----------------------------------------------------------------------------
// Function _checkRestrictions()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
bool _checkRestrictions(CommandLineOption const & opt, CharString const & val,
                        TErrorStream & estream)
{
    // we already know that the value can be converted to double/int, so we only need to check if it is
    // also in the required range
    if (isDoubleOption(opt))
    {
        return _checkMinMaxValue<double, TErrorStream>(opt, val, estream) && _checkValidValues(opt, val, estream);
    }
    if (isIntOption(opt))
    {
        return _checkMinMaxValue<int, TErrorStream>(opt, val, estream) && _checkValidValues(opt, val, estream);
    }
    if (isStringOption(opt))
    {
        return _checkValidValues(opt, val, estream);
    }

    // no restrictions to check
    return true;
}

// ----------------------------------------------------------------------------
// Function _assignOptionValue()
// ----------------------------------------------------------------------------

template <typename TErrorStream>
bool _assignOptionValue(CommandLineParser & me, CommandLineOption & opt, 
                        CharString const & val, unsigned argNo, TErrorStream & estream)
{
    if (isDoubleOption(opt) && !_isDouble(val))
    {
        _reportInvalidType(me, opt, val, estream);
        return false;
    }
    else if (isIntOption(opt) && !_isInt(val))
    {
        _reportInvalidType(me, opt, val, estream);
        return false;
    }

    if (!_checkRestrictions(opt, val, estream)) 
        return false;
    
    if (isOptionList(opt))
    {
        appendValue(opt.value, val, Generous());
    }
    else
    {
        if (argNo == 0) clear(opt.value);
        appendValue(opt.value, val, Exact());
    }
    return true;
}

template <typename TErrorStream>
bool
_assignOptionValue(CommandLineParser & me, CommandLineOption & opt, 
                   CharString const & val, TErrorStream & estream)
{
    return _assignOptionValue(me, opt, val, 0, estream);
}

template <typename TErrorStream>
bool
_assignOptionValue(CommandLineParser & me, unsigned option_index, 
                   CharString const & val, unsigned argNo, TErrorStream & estream)
{
    // get the option object
    CommandLineOption & opt = me.optionMap[option_index];
    return _assignOptionValue(me, opt, val, argNo, estream);
}


template <typename TErrorStream>
inline bool
_assignOptionValue(CommandLineParser & me, unsigned option_index, 
                   CharString const & val, TErrorStream & estream)
{
    return _assignOptionValue(me, option_index, val, 0, estream);
}

// ----------------------------------------------------------------------------
// Function _getOptionValues()
// ----------------------------------------------------------------------------

inline String<CharString> const &
_getOptionValues(CommandLineParser const &, CommandLineOption const & opt)
{
    if (empty(opt.value))
        return opt.defaultValue;
    else
        return opt.value;
}

inline String<CharString> const &
_getOptionValues(CommandLineParser const & me, unsigned option_index)
{
    return _getOptionValues(me, me.optionMap[option_index]);
}

// ----------------------------------------------------------------------------
// Function _getOptionValue()
// ----------------------------------------------------------------------------

inline CharString const &
_getOptionValue(CommandLineParser const & me, 
                CommandLineOption const & opt, unsigned argNo)
{
    if (argNo < length(opt.value))
        return opt.value[argNo];
    if (argNo < length(opt.defaultValue))
        return opt.defaultValue[argNo];
    
    return me._null; 
}

inline CharString const &
_getOptionValue(CommandLineParser const & me, CommandLineOption const & opt)
{
    return _getOptionValue(me, opt, 0);
}

inline CharString const &
_getOptionValue(CommandLineParser const & me, unsigned option_index, unsigned argNo)
{
    CommandLineOption const & opt = me.optionMap[option_index];
    return _getOptionValue(me, opt, argNo);
}

inline CharString const &
_getOptionValue(CommandLineParser const & me, unsigned option_index)
{
    return _getOptionValue(me, option_index, 0);
}

// ----------------------------------------------------------------------------
// Function _convertOptionValue()
// ----------------------------------------------------------------------------

inline bool
_convertOptionValue(CommandLineOption const & opt, bool & dst, CharString const & src)
{
    if (!isBooleanOption(opt))
        return false;

    dst = !empty(src);
    return true;
}

inline bool
_convertOptionValue(CommandLineOption const & opt, int & dst, CharString const & src)
{
    if (!isIntOption(opt))
        return false;

    std::istringstream stream(toCString(src));
    return !(stream >> dst).fail();
}

inline bool
_convertOptionValue(CommandLineOption const & opt, unsigned int & dst, CharString const & src)
{
    if (!isIntOption(opt))
        return false;

    std::istringstream stream(toCString(src));
    return !(stream >> dst).fail();
}

inline bool
_convertOptionValue(CommandLineOption const & opt, __int64 & dst, CharString const & src)
{
    if (!isIntOption(opt))
        return false;

    std::istringstream stream(toCString(src));
    return !(stream >> dst).fail();
}

inline bool
_convertOptionValue(CommandLineOption const & opt, __uint64 & dst, CharString const & src)
{
    if (!isIntOption(opt))
        return false;

    std::istringstream stream(toCString(src));
    return !(stream >> dst).fail();
}

inline bool
_convertOptionValue(CommandLineOption const & opt, float & dst, CharString const & src)
{
    if (!isDoubleOption(opt))
        return false;

    std::istringstream stream(toCString(src));
    return !(stream >> dst).fail();
}

inline bool
_convertOptionValue(CommandLineOption const & opt, double & dst, CharString const & src)
{
    if (!isDoubleOption(opt))
        return false;

    std::istringstream stream(toCString(src));
    return !(stream >> dst).fail();
}

template <typename TObject>
inline bool
_convertOptionValue(CommandLineOption const & opt, TObject & dst, CharString const & src)
{
    if (!isStringOption(opt))
        return false;

    assign(dst, src);
    return true;
}

// ----------------------------------------------------------------------------
// Function getOptionValue()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValue:
..summary:Retrieves the value of an option given either the short or long name.
..cat:Miscellaneous
..signature:getOptionValue(parser, optionIdentifier[, argNo], value)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that is either the short or long name of the option.
..param.argNo:If the option is list, the $argNo$-th list element is returned.
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TValue>
inline bool
getOptionValue(CommandLineParser const & me, CharString const & name, 
               unsigned argNo, TValue & val)
{
    SEQAN_ASSERT_MSG(hasOption(me, name), "Unknown option: %s", toCString(name));
    return _convertOptionValue(getOption(me, name), val, _getOptionValue(me, getOption(me, name), argNo));
}

template <typename TValue>
inline bool
getOptionValue(CommandLineParser const & me, CharString const & name,
               TValue & val)
{
    return getOptionValue(me, name, 0, val);
}

// ----------------------------------------------------------------------------
// Function getOptionValues()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValues
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getOptionValues(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that is either the short or long name of the option.
..returns: A $String<CharString>$ of option values.
..include:seqan/misc/misc_cmdparser.h
*/

inline String<CharString> const &
getOptionValues(CommandLineParser & me, CharString const & name)
{
    SEQAN_ASSERT_MSG(hasOption(me, name), "Unknown option: %s", toCString(name));
    return _getOptionValues(me, getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function getOptionValueShort()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValueShort
..summary:Retrieves the value of a short-name option given on the command line.
..cat:Miscellaneous
..signature:getOptionValueShort(parser, optionIdentifier[, argNo], value)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the short-name of the option.
..param.argNo:If the option is list, the $argNo$-th list element is returned.
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
..status:deprecated, use $Function.getOptionValue$
..see:Function.getOptionValue
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TValue>
inline bool
getOptionValueShort(CommandLineParser const & me, CharString
                    const & shortName, unsigned argNo, TValue & val)
{
    return getOptionValue(me, shortName, argNo, val);
}

template <typename TValue>
inline bool
getOptionValueShort(CommandLineParser const & me, 
                    CharString const & shortName, TValue & val)
{
    return getOptionValueShort(me, shortName, 0, val);
}

// ----------------------------------------------------------------------------
// Function getOptionValuesShort()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValuesShort
..summary:Returns all values of a short-name option given on the command line.
..cat:Miscellaneous
..signature:getOptionValuesShort(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the short-name of the option.
..returns: A $String<CharString>$ of option values.
..status:deprecated, use $Function.getOptionValues$
..see:Function.getOptionValues
..include:seqan/misc/misc_cmdparser.h
*/

inline String<CharString> const &
getOptionValuesShort(CommandLineParser & me, CharString const & shortName)
{
    return getOptionValues(me, shortName);
}

// ----------------------------------------------------------------------------
// Function getOptionValueLong()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValueLong
..summary:Retrieves the value of a long-name option given on the command line.
..cat:Miscellaneous
..signature:getOptionValueLong(parser, optionIdentifier[, argNo], value)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the long-name of the option.
..param.argNo:If the option is list, the $argNo$-th list element is returned.
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
..status:deprecated, use $Function.getOptionValue$
..see:Function.getOptionValue
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TValue>
inline bool
getOptionValueLong(CommandLineParser const & me, 
                   CharString const & longName, unsigned argNo, TValue & val)
{
    return getOptionValue(me, longName, argNo, val);
}

template <typename TValue>
inline bool
getOptionValueLong(CommandLineParser const & me, 
                   CharString const & longName, TValue & val)
{
    return getOptionValueLong(me, longName, 0, val);
}

// ----------------------------------------------------------------------------
// Function getOptionValuesLong()
// ----------------------------------------------------------------------------

/**
.Function.getOptionValuesLong
..summary:Returns all values of a long-name option given on the command line.
..cat:Miscellaneous
..signature:getOptionValuesLong(parser, optionIdentifier)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the long-name of the option.
..returns: A $String<CharString>$ of option values.
..include:seqan/misc/misc_cmdparser.h
..status:deprecated, use $Function.getOptionValues$
..see:Function.getOptionValues
*/

inline String<CharString> const &
getOptionValuesLong(CommandLineParser & me, CharString const & longName)
{
    return getOptionValues(me, longName);
}

// ----------------------------------------------------------------------------
// Function getArgumentValue()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValue
..summary:Returns an argument set on the command line.
..cat:Miscellaneous
..signature:getArgumentValue(parser, position)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.position:A zero based $int$ indicating which argument you want to get.
..returns:The command line argument or an empty string if it doesn't exist.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

inline CharString const &
getArgumentValue(CommandLineParser const & me, unsigned position)
{
    if (position < length(me._arguments))
        return me._arguments[position];
    else
        return me._null;
}

// ----------------------------------------------------------------------------
// Function getArgumentValues()
// ----------------------------------------------------------------------------

/**
.Function.getArgumentValues
..summary:Returns all _arguments set on the command line.
..cat:Miscellaneous
..signature:getArgumentValues(parser)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..returns:All command line _arguments as a $String<CharString>$.
..see:Function.getArgumentValue
..include:seqan/misc/misc_cmdparser.h
*/

inline String<CharString> const &
getArgumentValues(CommandLineParser const & me)
{
    return me._arguments;
}

// ----------------------------------------------------------------------------
// Function argumentCount()
// ----------------------------------------------------------------------------

/**
.Function.argumentCount
..summary:Returns the count of passed _arguments.
..cat:Miscellaneous
..signature:argumentCount(parser)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..include:seqan/misc/misc_cmdparser.h
*/

inline Size<String<CharString> >::Type
argumentCount(CommandLineParser const & me)
{
    return length(me._arguments);
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/**
.Function.setMinValue
..summary:Sets the minimum value of a @Class.CommandLineOption@ object identified by .
..cat:Miscellaneous
..signature:setMinValue(parser,optionName,minValue)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.option:The identifier of the command line option.
...type:Shortcut.CharString
..param.minValue:A @Shortcut.CharString@ containing a string representation of the minimum value of the @Class.CommandLineOption@.
..include:seqan/misc/misc_cmdparser.h
*/
inline void
setMinValue(CommandLineParser & me, CharString const & name, CharString const & _minValue)
{
    SEQAN_ASSERT_MSG(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMinValue(getOption(me, name), _minValue);
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/**
.Function.setMaxValue
..summary:Sets the maximum value of a @Class.CommandLineOption@ object.
..cat:Miscellaneous
..signature:setMaxValue(parser,optionName,maxValue)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.option:The identifier of the command line option.
...type:Shortcut.CharString
..param.maxValue:A @Shortcut.CharString@ containing a string representation of the maximum value of the @Class.CommandLineOption@.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setMaxValue(CommandLineParser & me, CharString const & name, 
            CharString const & _maxValue)
{
    SEQAN_ASSERT_MSG(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMaxValue(getOption(me, name), _maxValue);
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/**
.Function.setValidValues
..summary:Sets the set of allowed values of a @Class.CommandLineOption@ object.
..cat:Miscellaneous
..signature:setValidValues(parser,optionName,values)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..param.option:The identifier of the command line option.
...type:Shortcut.CharString
..param.values:A $String<CharString>$ containing all valid entries for the option.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setValidValues(CommandLineParser & me, CharString const & name, 
               StringSet<CharString> const & _values)
{
    SEQAN_ASSERT_MSG(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), _values);
}

inline void
setValidValues(CommandLineParser & me, CharString const & name, 
               CharString const & _values)
{
    // convert array to String<CharString>
    StringSet<CharString> values;
    CharString current_argument;

    for (Iterator<CharString const, Rooted>::Type ch  = begin(_values); ch != end(_values); goNext(ch))
    {
        if (*ch == ' ')
        {
            appendValue(values, current_argument);
            current_argument = "";
        }
        else
        {
            append(current_argument, *ch);
        }
    }
    if (current_argument != "")
        appendValue(values, current_argument);

    // add as restriction
    setValidValues(me, name, values);
}

// ----------------------------------------------------------------------------
// Function printEnvVarDocumentation()
// ----------------------------------------------------------------------------

template <typename TStream>
void printEnvVarDocumentation(TStream & stream, CommandLineParser const & me)
{
    // Get screen width from terminal size, default and minimal screen width.
    unsigned screenWidth = 0, screenHeight = 0;
    bool success = getTerminalSize(screenWidth, screenHeight);
    if (!success)
        screenWidth = me._defaultScreenWidth;
    if (screenWidth < me._minimalScreenWidth)
        screenWidth = me._minimalScreenWidth;
    if (me._maximalScreenWidth > 0u && screenWidth > me._maximalScreenWidth)
        screenWidth = me._maximalScreenWidth;
    screenWidth -= me._paddingRight;

    // Compute length of longest env var name.
    unsigned longestName = 0;
    for (unsigned i = 0; i < length(me.envVars); ++i)
    {
        EnvVarDescription const & desc = me.envVars[i];
        unsigned len = length(desc._name) + 1;
        if (me._paddingLeft + len + 1 > me._maximalLeftColumnWidth)
            continue;  // Ignore, will put description one line below.
        if (len > longestName)
            longestName  = len;
    }

    for (unsigned i = 0; i < length(me.envVars); ++i)
    {
        EnvVarDescription const & desc = me.envVars[i];

        // Tokenize text into words.
        StringSet<CharString> words;
        CharString word;
        for (unsigned i = 0; i < length(desc._text); ++i)
        {
            char c = desc._text[i];
            if (isspace(c))
            {
                if (!empty(word))
                    appendValue(words, word);
                clear(word);
            }
            else
            {
                appendValue(word, c);
            }
        }
        if (!empty(word))
            appendValue(words, word);

        // Compute tabs.
        unsigned tabVar = me._paddingLeft;
        unsigned tabDesc = tabVar + longestName + me._paddingCenter;
        
        // Print variable and description.
        unsigned s = 0;
        for (; s < tabVar; ++s)
            _streamPut(stream, ' ');
        _streamPut(stream, '$');
        _streamWrite(stream, desc._name);
        s += length(desc._name) + 1;
        for (; s < tabDesc; ++s)
            _streamPut(stream, ' ');

        // Now, print all the words.
        bool firstWord = true;
        for (unsigned i = 0; i < length(words); ++i)
        {
            if (s + 1 + length(words[i]) <= screenWidth)
            {
                // Fits on this line.  First, print space if not first word.
                if (!firstWord)
                {
                    _streamPut(stream, ' ');
                    s += 1;
                }
                firstWord = false;

                // Print word itself.
                _streamWrite(stream, words[i]);
                s += length(words[i]);
            }
            else
            {
                // Does not fit on this line.  Print on next.
                _streamPut(stream, '\n');
                for (s = 0; s < tabDesc; ++s)
                    _streamPut(stream, ' ');
                _streamWrite(stream, words[i]);
                s += length(words[i]);
            }
        }
        _streamPut(stream, '\n');
    }
}

}  // namespace seqan

#endif // CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_H_
