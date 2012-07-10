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

#ifndef SEQAN_CORE_TESTS_ARG_PARSE_TEST_ARG_PARSE_H_
#define SEQAN_CORE_TESTS_ARG_PARSE_TEST_ARG_PARSE_H_

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

const char * A_BOOL = "test-boolean-flags";
const char * A_BOOL_1 = "-b";
const char * A_BOOL_2 = "-c";
const char * A_BOOL_3 = "-bc";
const char * A_BOOL_4 = "-cb";

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

const char * A_TUPLE_LIST = "test_tuple_list";
const char * A_TUPLE_LIST_L = "-l";
const char * A_TUPLE_LIST_L_1 = "10";
const char * A_TUPLE_LIST_L_2 =  "20";
const char * A_TUPLE_LIST_L_3 =  "30";
const char * A_TUPLE_LIST_L_4 =  "40";

const char * A_TUPLE_LIST_DL = "-k";
const char * A_TUPLE_LIST_DL_1 = "5.1";
const char * A_TUPLE_LIST_DL_2 = "6.2";
const char * A_TUPLE_LIST_DL_3 = "7.3";
const char * A_TUPLE_LIST_DL_4 = "5.5";
const char * A_TUPLE_LIST_DL_5 = "6.6";
const char * A_TUPLE_LIST_DL_6 = "7.7";

namespace seqan
{

// moved initialization of cmd parser out of the test functions
// to have single place to change in case of interface changes
// or test extensions
void setupDoubleParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("d", "double","set a double option", ArgParseArgument::DOUBLE));
}

void setupIntegerParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("i", "integer","set an integer option", ArgParseArgument::INTEGER));
}

void setupStringParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("s", "string", "set a string option", ArgParseArgument::STRING, true));
}

void setupInputFileParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("i", "in", "set an input file", ArgParseArgument::INPUTFILE));
}

void setupOutputFileParser(ArgumentParser & parser)
{
    addOption(parser, ArgParseOption("o", "out", "set an output file", ArgParseArgument::OUTPUTFILE));
}

SEQAN_DEFINE_TEST(test_unset_value)
{
  ArgumentParser parser;
  setupIntegerParser(parser);

  int argc = 1;
  const char * argv[1] = {A_INT_0};

  std::stringstream error_stream;
  std::stringstream outputStream;

  SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
  SEQAN_ASSERT_EQ(error_stream.str(), "");
  SEQAN_ASSERT_EQ(outputStream.str(), "");

  int integerValue = 584864836;
  SEQAN_ASSERT(!getOptionValue(integerValue, parser, "integer"));
  SEQAN_ASSERT_EQ(integerValue, 584864836);
}

SEQAN_DEFINE_TEST(test_unset_values)
{
  ArgumentParser parser;
  setupIntegerParser(parser);

  int argc = 1;
  const char * argv[1] = {A_INT_0};

  std::stringstream error_stream;
  std::stringstream outputStream;

  SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
  SEQAN_ASSERT_EQ(error_stream.str(), "");
  SEQAN_ASSERT_EQ(outputStream.str(), "");

  std::vector<std::string> values = getOptionValues(parser, "integer");
  SEQAN_ASSERT_EQ(values.size(), 0u);
}


SEQAN_DEFINE_TEST(test_int_short_argument)
{

    ArgumentParser parser;
    setupIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_1, A_INT_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    int integerValue = 0;
    SEQAN_ASSERT(getOptionValue(integerValue, parser, "integer"));
    SEQAN_ASSERT_EQ(integerValue, 1);
}

SEQAN_DEFINE_TEST(test_int_long_argument)
{

    ArgumentParser parser;
    setupIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    int integerValue = 0;
    SEQAN_ASSERT(getOptionValue(integerValue, parser, "integer"));
    SEQAN_ASSERT_EQ(integerValue, 1);
}

SEQAN_DEFINE_TEST(test_non_int_argument)
{

    ArgumentParser parser;
    setupIntegerParser(parser);

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_1, A_INT_5};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'not-an-int' cannot be casted to integer\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_double_short_argument)
{

    ArgumentParser parser;
    setupDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    double doubleValue = 0.0;
    SEQAN_ASSERT(getOptionValue(doubleValue, parser, "double"));
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
}

SEQAN_DEFINE_TEST(test_double_long_argument)
{

    ArgumentParser parser;
    setupDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    double doubleValue = 0.0;
    SEQAN_ASSERT(getOptionValue(doubleValue, parser, "double"));
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
}

SEQAN_DEFINE_TEST(test_non_double_argument)
{

    ArgumentParser parser;
    setupDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_5};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'not-a-double' cannot be casted to double\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_double_scientific_notation)
{

    ArgumentParser parser;
    setupDoubleParser(parser);

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_1, A_DOUBLE_6};

    std::stringstream error_stream;
    std::stringstream outputStream;

    double doubleValue = 0.0;
    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    SEQAN_ASSERT(getOptionValue(doubleValue, parser, "double"));
    SEQAN_ASSERT_EQ(doubleValue, 6.0221418e23);
}

SEQAN_DEFINE_TEST(test_string_short_argument)
{

    ArgumentParser parser;
    setupStringParser(parser);

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_1, A_STRING_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "string"));
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
}

SEQAN_DEFINE_TEST(test_string_long_argument)
{

    ArgumentParser parser;
    setupStringParser(parser);

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_2, A_STRING_3};

    std::stringstream error_stream;
    std::stringstream outputStream;
    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "string"));
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
}

SEQAN_DEFINE_TEST(test_string_missing_argument)
{

    ArgumentParser parser;
    setupStringParser(parser);

    int argc = 2;
    const char * argv[2] = {A_STRING_0, A_STRING_2};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: option requires an argument -- string\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_string_list)
{
    ArgumentParser parser;
    setupStringParser(parser);

    int argc = 7;
    const char * argv[7] = {A_STRING_0, A_STRING_1, A_STRING_3, A_STRING_2, A_STRING_3, A_STRING_1, A_STRING_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    std::vector<std::string> const & values = getOptionValues(parser, "string");

    SEQAN_ASSERT_EQ(length(values), 3u);

    for (unsigned i = 0; i < length(values); ++i)
    {
        SEQAN_ASSERT_EQ(value(values, i), "this-is-a-string-value");
    }
}

SEQAN_DEFINE_TEST(test_min_max_double_values_in_range)
{
    ArgumentParser parser;
    setupDoubleParser(parser);

    setMinValue(parser, "double", "1.0");
    setMaxValue(parser, "double", "2.0");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    double doubleValue = 0.0;
    SEQAN_ASSERT(getOptionValue(doubleValue, parser, "double"));
    SEQAN_ASSERT_EQ(doubleValue, 1.56);
}

SEQAN_DEFINE_TEST(test_min_max_double_values_to_small)
{
    ArgumentParser parser;
    setupDoubleParser(parser);

    setMinValue(parser, "double", "1.6");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value '1.56' is not in the interval [1.6:+inf]\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_double_values_to_big)
{
    ArgumentParser parser;
    setupDoubleParser(parser);

    setMaxValue(parser, "double", "1.5");

    int argc = 3;
    const char * argv[3] = {A_DOUBLE_0, A_DOUBLE_2, A_DOUBLE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value '1.56' is not in the interval [-inf:1.5]\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_int_values_in_range)
{
    ArgumentParser parser;
    setupIntegerParser(parser);

    setMinValue(parser, "integer", "-10");
    setMaxValue(parser, "integer", "2");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    int integerValue = 0;
    SEQAN_ASSERT(getOptionValue(integerValue, parser, "integer"));
    SEQAN_ASSERT_EQ(integerValue, 1);
}

SEQAN_DEFINE_TEST(test_min_max_int_values_to_small)
{
    ArgumentParser parser;
    setupIntegerParser(parser);

    setMinValue(parser, "integer", "3");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value '1' is not in the interval [3:+inf]\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_min_max_int_values_to_big)
{
    ArgumentParser parser;
    setupIntegerParser(parser);

    setMaxValue(parser, "integer", "-3");

    int argc = 3;
    const char * argv[3] = {A_INT_0, A_INT_2, A_INT_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value '1' is not in the interval [-inf:-3]\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_allowed_values_contained)
{
    ArgumentParser parser;
    setupStringParser(parser);

    setValidValues(parser, "string", "a b c this-is-a-string-value");

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_2, A_STRING_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "string"));
    SEQAN_ASSERT_EQ(value, "this-is-a-string-value");
}

SEQAN_DEFINE_TEST(test_allowed_values_not_contained)
{
    ArgumentParser parser;
    setupStringParser(parser);

    setValidValues(parser, "string", "a b c");

    int argc = 3;
    const char * argv[3] = {A_STRING_0, A_STRING_2, A_STRING_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'this-is-a-string-value' is not in the list of allowed values [a, b, c]\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_input_file_short)
{
    ArgumentParser parser;
    setupInputFileParser(parser);

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_1, A_IN_FILE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "in"));
    SEQAN_ASSERT_EQ(value, "input.fasta");
}

SEQAN_DEFINE_TEST(test_input_file_long)
{
    ArgumentParser parser;
    setupInputFileParser(parser);

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_2, A_IN_FILE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "in"));
    SEQAN_ASSERT_EQ(value, "input.fasta");
}

SEQAN_DEFINE_TEST(test_input_file_missing)
{
    ArgumentParser parser;
    setupInputFileParser(parser);

    int argc = 2;
    const char * argv[2] = {A_IN_FILE_0, A_IN_FILE_1};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: option requires an argument -- i\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_input_file_invalid_type)
{
    ArgumentParser parser;
    setupInputFileParser(parser);

    setValidValues(parser, "in", "FASTA fa");

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_2, A_IN_FILE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'input.fasta' is not in the list of allowed file extensions [*.FASTA, *.fa]\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_input_file_valid_type)
{
    ArgumentParser parser;
    setupInputFileParser(parser);

    setValidValues(parser, "in", "fasta FASTA fa");

    int argc = 3;
    const char * argv[3] = {A_IN_FILE_0, A_IN_FILE_2, A_IN_FILE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "in"));
    SEQAN_ASSERT_EQ(value, "input.fasta");
}

SEQAN_DEFINE_TEST(test_output_file_short)
{
    ArgumentParser parser;
    setupOutputFileParser(parser);

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_1, A_OUT_FILE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "out"));
    SEQAN_ASSERT_EQ(value, "output.fasta");
}

SEQAN_DEFINE_TEST(test_output_file_long)
{
    ArgumentParser parser;
    setupOutputFileParser(parser);

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_2, A_OUT_FILE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "out"));
    SEQAN_ASSERT_EQ(value, "output.fasta");
}

SEQAN_DEFINE_TEST(test_output_file_missing)
{
    ArgumentParser parser;
    setupOutputFileParser(parser);

    int argc = 2;
    const char * argv[2] = {A_OUT_FILE_0, A_OUT_FILE_1};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: option requires an argument -- o\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_invalid_type)
{
    ArgumentParser parser;
    setupOutputFileParser(parser);

    setValidValues(parser, "out", "FASTA fa");

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_2, A_OUT_FILE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'output.fasta' is not in the list of allowed file extensions [*.FASTA, *.fa]\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_output_file_valid_type)
{
    ArgumentParser parser;
    setupOutputFileParser(parser);

    setValidValues(parser, "out", "fasta FASTA fa");

    int argc = 3;
    const char * argv[3] = {A_OUT_FILE_0, A_OUT_FILE_2, A_OUT_FILE_3};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    CharString value;
    SEQAN_ASSERT(getOptionValue(value, parser, "out"));
    SEQAN_ASSERT_EQ(value, "output.fasta");
}

SEQAN_DEFINE_TEST(test_argument_string)
{
    ArgumentParser parser;
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));

    int argc = 3;
    const char * argv[3] = {A_ARGUMENT_0, A_ARGUMENT_1, A_ARGUMENT_2};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_argument_not_all_set)
{
    ArgumentParser parser;
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING));

    int argc = 2;
    const char * argv[2] = {A_ARGUMENT_0, A_ARGUMENT_1};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: Not all required arguments or options were set\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

SEQAN_DEFINE_TEST(test_argument_double)
{
    ArgumentParser parser;
    addArgument(parser, ArgParseArgument(ArgParseArgument::DOUBLE));

    int argc = 2;
    const char * argv[2] = {A_ARGUMENT_0, A_ARGUMENT_DOUBLE_5};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

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
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test: the given value 'argument1' cannot be casted to double\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
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

// Testing lists of tuples
SEQAN_DEFINE_TEST(test_int_list_option)
{
    ArgumentParser parser;
    addOption(parser, ArgParseOption("l", "list", "this is a list option", ArgParseArgument::INTEGER, true, "",2));

    int argc = 7;
    const char* argv[7] = {A_TUPLE_LIST, A_TUPLE_LIST_L, A_TUPLE_LIST_L_1, A_TUPLE_LIST_L_2, A_TUPLE_LIST_L, A_TUPLE_LIST_L_3, A_TUPLE_LIST_L_4};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    SEQAN_ASSERT_EQ(getOptionValueCount(parser, "list"), 4u);

    int value = 0;
    SEQAN_ASSERT(getOptionValue(value, parser, "list", 0));
    SEQAN_ASSERT_EQ(value, 10);
    SEQAN_ASSERT(getOptionValue(value, parser, "list", 1));
    SEQAN_ASSERT_EQ(value, 20);

    SEQAN_ASSERT(getOptionValue(value, parser, "list", 2));
    SEQAN_ASSERT_EQ(value, 30);
    SEQAN_ASSERT(getOptionValue(value, parser, "list", 3));
    SEQAN_ASSERT_EQ(value, 40);
}

SEQAN_DEFINE_TEST(test_double_list_option)
{
    ArgumentParser parser;
    addOption(parser, ArgParseOption("k", "double-list", "this is a list option", ArgParseArgument::DOUBLE, true, "", 3));

    int argc = 9;
    const char* argv[9] = {A_TUPLE_LIST, A_TUPLE_LIST_DL, A_TUPLE_LIST_DL_1, A_TUPLE_LIST_DL_2, A_TUPLE_LIST_DL_3, A_TUPLE_LIST_DL, A_TUPLE_LIST_DL_4, A_TUPLE_LIST_DL_5, A_TUPLE_LIST_DL_6};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    SEQAN_ASSERT_EQ(getOptionValueCount(parser, "double-list"), 6u);

    double value = 0;
    SEQAN_ASSERT(getOptionValue(value, parser, "double-list", 0));
    SEQAN_ASSERT_EQ(value, 5.1);
    SEQAN_ASSERT(getOptionValue(value, parser, "double-list", 1));
    SEQAN_ASSERT_EQ(value, 6.2);
    SEQAN_ASSERT(getOptionValue(value, parser, "double-list", 2));
    SEQAN_ASSERT_EQ(value, 7.3);

    SEQAN_ASSERT(getOptionValue(value, parser, "double-list", 3));
    SEQAN_ASSERT_EQ(value, 5.5);
    SEQAN_ASSERT(getOptionValue(value, parser, "double-list", 4));
    SEQAN_ASSERT_EQ(value, 6.6);
    SEQAN_ASSERT(getOptionValue(value, parser, "double-list", 5));
    SEQAN_ASSERT_EQ(value, 7.7);

}

SEQAN_DEFINE_TEST(test_double_list_option_not_enough_arguments)
{
    ArgumentParser parser;
    addOption(parser, ArgParseOption("k", "double-list", "this is a list option", ArgParseArgument::DOUBLE, true, "", 3));

    int argc = 8;
    const char* argv[8] = {A_TUPLE_LIST, A_TUPLE_LIST_DL, A_TUPLE_LIST_DL_1, A_TUPLE_LIST_DL_2, A_TUPLE_LIST_DL_3, A_TUPLE_LIST_DL, A_TUPLE_LIST_DL_4, A_TUPLE_LIST_DL_5};

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_ERROR);
    SEQAN_ASSERT_EQ(error_stream.str(), "test_tuple_list: option requires an argument -- k\n");
    SEQAN_ASSERT_EQ(outputStream.str(), "");
}

void setUpBoolParser(ArgumentParser& parser)
{
    addOption(parser, ArgParseOption("b", "", "This is a boolean flag"));
    addOption(parser, ArgParseOption("c", "", "This is a boolean flag"));
    addOption(parser, ArgParseOption("d", "", "This is a boolean flag that we will not set."));
}

SEQAN_DEFINE_TEST(test_boolean_flags)
{
    ArgumentParser parser;
    setUpBoolParser(parser);

    int argc = 3;
    const char* argv[3] = { A_BOOL, A_BOOL_1, A_BOOL_2 };

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    bool isSet = false;
    SEQAN_ASSERT(getOptionValue(isSet, parser, "b"));
    SEQAN_ASSERT(isSet);

    isSet = false;
    SEQAN_ASSERT(getOptionValue(isSet, parser, "c"));
    SEQAN_ASSERT(isSet);

    isSet = true;
    SEQAN_ASSERT(getOptionValue(isSet, parser, "d"));
    SEQAN_ASSERT(!isSet);
}

SEQAN_DEFINE_TEST(test_combined_boolean_flags)
{
   ArgumentParser parser;
   setUpBoolParser(parser);

   int argc = 2;
   const char* argv[2] = { A_BOOL, A_BOOL_3 };

   std::stringstream error_stream;
   std::stringstream outputStream;

   SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
   SEQAN_ASSERT_EQ(error_stream.str(), "");
   SEQAN_ASSERT_EQ(outputStream.str(), "");

   bool isSet = false;
   SEQAN_ASSERT(getOptionValue(isSet, parser, "b"));
   SEQAN_ASSERT(isSet);

   isSet = false;
   SEQAN_ASSERT(getOptionValue(isSet, parser, "c"));
   SEQAN_ASSERT(isSet);

   isSet = true;
   SEQAN_ASSERT(getOptionValue(isSet, parser, "d"));
   SEQAN_ASSERT(!isSet);
}

SEQAN_DEFINE_TEST(test_long_short_flag_name)
{
    ArgumentParser parser;
    addOption(parser, ArgParseOption("bc", "", "This is a boolean flag"));
    addOption(parser, ArgParseOption("d", "", "This is a boolean flag"));

    int argc = 2;
    const char* argv[2] = { A_BOOL, A_BOOL_3 };

    std::stringstream error_stream;
    std::stringstream outputStream;

    SEQAN_ASSERT_EQ(parse(parser, argc, argv, outputStream, error_stream), ArgumentParser::PARSE_OK);
    SEQAN_ASSERT_EQ(error_stream.str(), "");
    SEQAN_ASSERT_EQ(outputStream.str(), "");

    bool isSet = false;
    SEQAN_ASSERT(getOptionValue(isSet, parser, "bc"));
    SEQAN_ASSERT(isSet);

    isSet = true;
    SEQAN_ASSERT(getOptionValue(isSet, parser, "d"));
    SEQAN_ASSERT(!isSet);
}

} // namespace seqan

#endif  // SEQAN_CORE_TESTS_ARG_PARSE_TEST_ARG_PARSE_H_
