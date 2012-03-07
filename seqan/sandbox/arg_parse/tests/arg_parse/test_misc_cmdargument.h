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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================
// Tests for misc/cmdparser/cmdargument.h.
// ==========================================================================

#ifndef SEQAN_HEADER_TEST_MISC_CMDARGUMENT_H_
#define SEQAN_HEADER_TEST_MISC_CMDARGUMENT_H_

#include <seqan/basic.h>

#include <seqan/arg_parse/cmdargument.h>

namespace seqan
{

SEQAN_DEFINE_TEST(test_argument_string_label)
{
    CommandLineArgument arg1(CommandLineArgument::STRING);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "STR");

    arg1._numberOfArguments = 2;
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "STR STR");

    CommandLineArgument arg2(CommandLineArgument::STRING);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg2), "STR");
}

SEQAN_DEFINE_TEST(test_argument_int_label)
{
    CommandLineArgument arg1(CommandLineArgument::INTEGER);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "NUM");

    arg1._numberOfArguments = 2;
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "NUM NUM");

    CommandLineArgument arg2(CommandLineArgument::INTEGER);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg2), "NUM");
}

SEQAN_DEFINE_TEST(test_argument_double_label)
{
    CommandLineArgument arg1(CommandLineArgument::DOUBLE);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "NUM");

    arg1._numberOfArguments = 2;
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "NUM NUM");

    CommandLineArgument arg2(CommandLineArgument::DOUBLE);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg2), "NUM");
}

SEQAN_DEFINE_TEST(test_argument_inputfile_label)
{
    CommandLineArgument arg1(CommandLineArgument::INPUTFILE);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "FILE");

    arg1._numberOfArguments = 2;
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "FILE FILE");

    CommandLineArgument arg2(CommandLineArgument::INPUTFILE);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg2), "FILE");
}

SEQAN_DEFINE_TEST(test_argument_outputfile_label)
{
    CommandLineArgument arg1(CommandLineArgument::OUTPUTFILE);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "FILE");

    arg1._numberOfArguments = 2;
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "FILE FILE");

    CommandLineArgument arg2(CommandLineArgument::OUTPUTFILE);
    SEQAN_ASSERT_EQ(getArgumentLabel(arg2), "FILE");
}

SEQAN_DEFINE_TEST(test_argument_user_defined_label)
{
    CommandLineArgument arg1(CommandLineArgument::STRING, false, "my_label");
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "my_label");

    arg1._numberOfArguments = 2;
    SEQAN_ASSERT_EQ(getArgumentLabel(arg1), "my_label");
}

SEQAN_DEFINE_TEST(test_argument_min_max_boundaries)
{
    CommandLineArgument arg(CommandLineArgument::INTEGER);
    setMinValue(arg, "1");
    setMaxValue(arg, "10");


}

}  // namespace seqan

#endif // SEQAN_HEADER_TEST_MISC_CMDARGUMENT_H_
