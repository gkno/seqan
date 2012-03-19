// ==========================================================================
//                                 arg_parse
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

#ifndef SANDBOX_ARG_PARSE_TESTS_ARG_PARSE_TEST_ARG_PARSE_H_
#define SANDBOX_ARG_PARSE_TESTS_ARG_PARSE_TEST_ARG_PARSE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>


#include <seqan/arg_parse.h> // Header under test

#include <seqan/basic.h>
#include <iostream>

using namespace std;

const char * A_INT_0 = "test";
const char * A_INT_1 = "-i";
const char * A_INT_2 = "--integer";
const char * A_INT_3 = "1";
const char * A_INT_5 = "not-an-int";

const char * A_DOUBLE_0 = "test";
const char * A_DOUBLE_1 = "-d";
const char * A_DOUBLE_2 = "--double";
const char * A_DOUBLE_3 = "1.56";
const char * A_DOUBLE_5 = "not-a-double";
const char * A_DOUBLE_6 = "6.0221418e23";

const char * A_STRING_0 = "test";
const char * A_STRING_1 = "-s";
const char * A_STRING_2 = "--string";
const char * A_STRING_3 = "this-is-a-string-value";

const char * A_IN_FILE_0 = "test";
const char * A_IN_FILE_1 = "-i";
const char * A_IN_FILE_2 = "--in";
const char * A_IN_FILE_3 = "input.fasta";

const char * A_OUT_FILE_0 = "test";
const char * A_OUT_FILE_1 = "-o";
const char * A_OUT_FILE_2 = "--out";
const char * A_OUT_FILE_3 = "output.fasta";

const char * A_ARGUMENT_0 = "test";
const char * A_ARGUMENT_1 = "argument1";
const char * A_ARGUMENT_2 = "argument2";
const char * A_ARGUMENT_3 = "argument3";
const char * A_ARGUMENT_INT_4 = "-10";
const char * A_ARGUMENT_DOUBLE_5 = "6.0221418e23";

namespace seqan {

// moved initialization of cmd parser out of the test functions
// to have single place to change in case of interface changes
// or test extensions
void testInitDoubleParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("d", "double","set a double option", ArgParseArgument(ArgParseArgument::DOUBLE)));
}

void testInitIntegerParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("i", "integer","set an integer option", ArgParseArgument(ArgParseArgument::INTEGER)));
}

void testInitStringParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("s", "string", "set a string option", ArgParseArgument(ArgParseArgument::STRING, true)));
}

void testInFileTypeParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("i", "in", "set an input file", ArgParseArgument(ArgParseArgument::INPUTFILE)));
}

void testOutFileTypeParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("o", "out", "set an output file", ArgParseArgument(ArgParseArgument::OUTPUTFILE)));
}

SEQAN_DEFINE_TEST(test_int_short_argument)
{

    ArgumentParser parser;
    testInitIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_1, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    int integerValue = 0;
    SEQAN_ASSERT(getOptionValue(integerValue, parser, "integer"));
    SEQAN_ASSERT_EQ(integerValue, 1);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_int_long_argument)
{

    ArgumentParser parser;
    testInitIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    int integerValue = 0;
    SEQAN_ASSERT(getOptionValue(integerValue, parser, "integer"));
    SEQAN_ASSERT_EQ(integerValue, 1);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_non_int_argument)
{

    ArgumentParser parser;
    testInitIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_1, A_INT_5};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'not-an-int' cannot be casted to integer\n");
}

SEQAN_DEFINE_TEST(test_double_short_argument)
{

    ArgumentParser parser;
    testInitDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    double doubleValue = 0.0;
    SEQAN_ASSERT(getOptionValue(doubleValue, parser, "double"));
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
    SEQAN_ASSERT_EQ(error_stream.str(), "");

}

SEQAN_DEFINE_TEST(test_double_long_argument)
{

    ArgumentParser parser;
    testInitDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    double doubleValue = 0.0;
    SEQAN_ASSERT(getOptionValue(doubleValue, parser, "double"));
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_non_double_argument)
{

    ArgumentParser parser;
    testInitDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_5};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'not-a-double' cannot be casted to double\n");
}

SEQAN_DEFINE_TEST(test_double_scientific_notation)
{

    ArgumentParser parser;
    testInitDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_6};

    std::stringstream error_stream;

    double doubleValue = 0.0;
    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);
    SEQAN_ASSERT(getOptionValue(doubleValue, parser, "double"));
    SEQAN_ASSERT_EQ(doubleValue, 6.0221418e23);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_string_short_argument)
{

    ArgumentParser parser;
    testInitStringParser(parser);

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_1, A_STRING_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "string"));
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_string_long_argument)
{

    ArgumentParser parser;
    testInitStringParser(parser);

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_2, A_STRING_3};

    std::stringstream error_stream;
    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "string"));
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_string_missing_argument)
{

    ArgumentParser parser;
    testInitStringParser(parser);

    int argc = 2;
    const char * argv[2] = {A_STRING_0, A_STRING_2};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: option requires an argument -- string\n");
}

SEQAN_DEFINE_TEST(test_string_list)
{
    ArgumentParser parser;
    testInitStringParser(parser);

    int argc = 7;
    const char * argv[7] = {A_STRING_0, A_STRING_1, A_STRING_3, A_STRING_2, A_STRING_3, A_STRING_1, A_STRING_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    std::vector<std::string> const & values = getOptionValues(parser, "string");

    SEQAN_ASSERT_EQ(length(values), 3u);

    for (unsigned i = 0; i < length(values); ++i)
    {
        SEQAN_ASSERT_EQ(value(values, i), "this-is-a-string-value");
    }

    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_double_values_in_range)
{
    ArgumentParser parser;
    testInitDoubleParser(parser);

    setMinValue(parser, "double", "1.0");
    setMaxValue(parser, "double", "2.0");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    double doubleValue = 0.0;
    SEQAN_ASSERT(getOptionValue(doubleValue, parser, "double"));
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_double_values_to_small)
{
    ArgumentParser parser;
    testInitDoubleParser(parser);

    setMinValue(parser, "double", "1.6");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value '1.56' is not in the interval [1.6:+inf]\n");
}

SEQAN_DEFINE_TEST(test_min_max_double_values_to_big)
{
    ArgumentParser parser;
    testInitDoubleParser(parser);

    setMaxValue(parser, "double", "1.5");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value '1.56' is not in the interval [-inf:1.5]\n");
}

SEQAN_DEFINE_TEST(test_min_max_int_values_in_range)
{
    ArgumentParser parser;
    testInitIntegerParser(parser);

    setMinValue(parser, "integer", "-10");
    setMaxValue(parser, "integer", "2");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    int integerValue = 0;
    SEQAN_ASSERT(getOptionValue(integerValue, parser, "integer"));
    SEQAN_ASSERT_EQ(integerValue, 1);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_int_values_to_small)
{
    ArgumentParser parser;
    testInitIntegerParser(parser);

    setMinValue(parser, "integer", "3");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value '1' is not in the interval [3:+inf]\n");
}

SEQAN_DEFINE_TEST(test_min_max_int_values_to_big)
{
    ArgumentParser parser;
    testInitIntegerParser(parser);

    setMaxValue(parser, "integer", "-3");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value '1' is not in the interval [-inf:-3]\n");
}

SEQAN_DEFINE_TEST(test_allowed_values_contained)
{
    ArgumentParser parser;
    testInitStringParser(parser);

    setValidValues(parser, "string", "a b c this-is-a-string-value");

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_2, A_STRING_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "string"));
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_allowed_values_not_contained)
{
    ArgumentParser parser;
    testInitStringParser(parser);

    setValidValues(parser, "string", "a b c");

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_2, A_STRING_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'this-is-a-string-value' is not in the list of allowed values [a, b, c]\n");
}

SEQAN_DEFINE_TEST(test_input_file_short)
{
    ArgumentParser parser;
    testInFileTypeParser(parser);

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_1, A_IN_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "in"));
    SEQAN_ASSERT_EQ(value, "input.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_input_file_long)
{
    ArgumentParser parser;
    testInFileTypeParser(parser);

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_2, A_IN_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "in"));
    SEQAN_ASSERT_EQ(value, "input.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_input_file_missing)
{
    ArgumentParser parser;
    testInFileTypeParser(parser);

    int argc = 2;
    const char * argv[2] = {A_IN_FILE_0, A_IN_FILE_1};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: option requires an argument -- i\n");
}

SEQAN_DEFINE_TEST(test_input_file_invalid_type)
{
    ArgumentParser parser;
    testInFileTypeParser(parser);

    setValidValues(parser, "in", "FASTA fa");

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_2, A_IN_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'input.fasta' is not in the list of allowed file extensions [*.FASTA, *.fa]\n");
}

SEQAN_DEFINE_TEST(test_input_file_valid_type)
{
    ArgumentParser parser;
    testInFileTypeParser(parser);

    setValidValues(parser, "in", "fasta FASTA fa");

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_2, A_IN_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "in"));
    SEQAN_ASSERT_EQ(value, "input.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_short)
{
    ArgumentParser parser;
    testOutFileTypeParser(parser);

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_1, A_OUT_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "out"));
    SEQAN_ASSERT_EQ(value, "output.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_long)
{
    ArgumentParser parser;
    testOutFileTypeParser(parser);

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_2, A_OUT_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "out"));
    SEQAN_ASSERT_EQ(value, "output.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_missing)
{
    ArgumentParser parser;
    testOutFileTypeParser(parser);

    int argc = 2;
    const char * argv[2] = {A_OUT_FILE_0, A_OUT_FILE_1};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: option requires an argument -- o\n");
}

SEQAN_DEFINE_TEST(test_output_file_invalid_type)
{
    ArgumentParser parser;
    testOutFileTypeParser(parser);

    setValidValues(parser, "out", "FASTA fa");

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_2, A_OUT_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'output.fasta' is not in the list of allowed file extensions [*.FASTA, *.fa]\n");
}

SEQAN_DEFINE_TEST(test_output_file_valid_type)
{
    ArgumentParser parser;
    testOutFileTypeParser(parser);

    setValidValues(parser, "out", "fasta FASTA fa");

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_2, A_OUT_FILE_3};

    std::stringstream error_stream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "out"));
    SEQAN_ASSERT_EQ(value, "output.fasta");
    SEQAN_ASSERT_EQ(error_stream.str(), "");
}

SEQAN_DEFINE_TEST(test_argument_string)
{
    ArgumentParser parser;
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));

    int argc = 3;
    const char * argv[3] = {A_ARGUMENT_0, A_ARGUMENT_1, A_ARGUMENT_2};

    std::stringstream error_stream;
    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);

}

SEQAN_DEFINE_TEST(test_argument_not_all_set)
{
    ArgumentParser parser;
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));

    int argc = 2;
    const char * argv[2] = {A_ARGUMENT_0, A_ARGUMENT_1};

    std::stringstream error_stream;
    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
}

SEQAN_DEFINE_TEST(test_argument_double)
{
    ArgumentParser parser;
    addArgument(parser, ArgParseArgument(ArgParseArgument::DOUBLE));

    int argc = 2;
    const char * argv[2] = {A_ARGUMENT_0, A_ARGUMENT_DOUBLE_5};

    std::stringstream error_stream;
    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::OK);
    double doubleValue = 0.0;
    SEQAN_ASSERT(getArgumentValue(doubleValue, parser, 0));
    SEQAN_ASSERT_EQ(doubleValue, 6.0221418e23);
}

SEQAN_DEFINE_TEST(test_argument_not_a_double)
{
    ArgumentParser parser;
    addArgument(parser, ArgParseArgument(ArgParseArgument::DOUBLE));

    int argc = 2;
    const char * argv[2] = {A_ARGUMENT_0, A_ARGUMENT_1};

    std::stringstream error_stream;
    SEQAN_ASSERT_EQ(parse(parser, argc, argv, error_stream), ArgumentParser::ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'argument1' cannot be casted to double\n");
}

SEQAN_DEFINE_TEST(test_isDouble)
{
    CharString a = "this is not a double";
    CharString b = "2.5";
    CharString c = "-45.5245";
    CharString d = "6.0221418e23";
    CharString e = "-45.5245aeeeb";

    SEQAN_ASSERT(!_isDouble(a));
    SEQAN_ASSERT(_isDouble(b));
    SEQAN_ASSERT(_isDouble(c));
    SEQAN_ASSERT(_isDouble(d));
    SEQAN_ASSERT(!_isDouble(e));
}

SEQAN_DEFINE_TEST(test_isInt)
{
    CharString a = "this is not an int";
    CharString b = "2";
    CharString c = "-4253252";
    CharString d = "6aaefgeag";

    SEQAN_ASSERT(!_isInt(a));
    SEQAN_ASSERT(_isInt(b));
    SEQAN_ASSERT(_isInt(c));
    SEQAN_ASSERT(!_isInt(d));
}

} // namespace seqan

#endif  // SANDBOX_ARG_PARSE_TESTS_ARG_PARSE_TEST_ARG_PARSE_H_