/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ===========================================================================
  Copyright (C) 2007-2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  The SeqAn testing infrastructure.  Based on ideas from the OpenMS
  "ClassTest.h".
  ===========================================================================*/
#ifndef SEQAN_BASIC_BASIC_TESTING_H_
#define SEQAN_BASIC_BASIC_TESTING_H_

#include <iostream>  // stdout, stderr

namespace seqan {

// Namespace for the testing infrastructure.
//
// This namespace contains the variables and functions that are used
// in the macros below to perform the tests.
namespace ClassTest {
    // Number of tests that were run.
    int testCount;

    // Number of errors that occured.
    int errorCount;

    // Number of skipped tests.
    int skippedCount;

    // Flag whether there was an error in this test.
    bool thisTestOk;

    // Flag whether this test was skipped.
    bool thisTestSkipped;

    // Name of the current test.
    const char *currentTestName;

    // Initialize the testing infrastructure.
    //
    // Used through SEQAN_BEGIN_TESTSUITE(test_name)
    void beginTestSuite(const char *testName) {
        testCount = 0;
        skippedCount = 0;
        errorCount = 0;
    }

    // Run test suite finalization.
    //
    // Used through SEQAN_END_TESTSUITE
    //
    // Prints a bottom banner with the error count and returns the
    // program's return code.
    int endTestSuite() {
        std::cout << "**************************************" << std::endl;
        std::cout << " Total Tests: " << testCount << std::endl;
        std::cout << " Skipped:     " << skippedCount << std::endl;
        std::cout << " Errors:      " << errorCount << std::endl;
        std::cout << "**************************************" << std::endl;
        return errorCount > 0;
    }

    // Run test initialization.
    void beginTest(const char *testName) {
        currentTestName = testName;
        thisTestOk = true;
        thisTestSkipped = false;
        testCount += 1;
    }

    // Run test finalization.
    void endTest() {
        if (thisTestSkipped) {
            std::cout << currentTestName << " SKIPPED" << std::endl;
        } else if (thisTestOk) {
            std::cout << currentTestName << " OK" << std::endl;
        } else {
            std::cerr << currentTestName << " FAILED" << std::endl;
        }
    }

    // Marks the current test as "skipped".
    void skipCurrentTest() {
        thisTestSkipped = true;
        skippedCount += 1;
    }

    // Called by the macro SEQAN_ASSERT_EQ.
    //
    // Tests that the given two value are equal.  Returns true iff the
    // two values are equal.
    template <typename T1, typename T2>
    bool testEqual(const char *file, int line,
                   const T1 &value1, const char *expression1,
                   const T2 &value2, const char *expression2,
                   const char *comment = 0) {
        if (not (value1 == value2)) {
            // Increase global error count.
            thisTestOk = false;
            errorCount += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " == " << expression2 << " was: " << value1
                      << " != " << value2;
            if (comment)
                std::cerr << " (" << comment << ")";
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Called by the macro SEQAN_ASSERT_NEQ.
    //
    // Tests that the given two value are not equal.  Returns true iff
    // the two values are equal.
    template <typename T1, typename T2>
    bool testNotEqual(const char *file, int line,
                   const T1 &value1, const char *expression1,
                   const T2 &value2, const char *expression2,
                   const char *comment = 0) {
        if (not (value1 != value2)) {
            // Increase global error count.
            thisTestOk = false;
            errorCount += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression1 << " != " << expression2 << " was: " << value1
                      << " == " << value2;
            if (comment)
                std::cerr << " (" << comment << ")";
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Called by the macro SEQAN_ASSERT.
    //
    // Test that the given argument evaluates to true.
    template <typename T>
    bool testTrue(const char *file, int line,
                  const T &value_, const char *expression_,
                  const char *comment = 0) {
        if (not (value_)) {
            // Increase global error count.
            thisTestOk = false;
            errorCount += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression_ << " should be true but was " << (value_);
            if (comment)
                std::cerr << " (" << comment << ")";
            std::cerr << std::endl;
            return false;
        }
        return true;
    }


    // Called by the macro SEQAN_ASSERT.
    //
    // Test that the given argument evaluates to false.
    template <typename T>
    bool testFalse(const char *file, int line,
                   const T &value_, const char *expression_,
                   const char *comment = 0) {
        if (value_) {
            // Increase global error count.
            thisTestOk = false;
            errorCount += 1;
            // Print assertion failure text, with comment if any is given.
            std::cerr << file << ":" << line << " Assertion failed : "
                      << expression_ << " should be false but was " << (value_);
            if (comment)
                std::cerr << " (" << comment << ")";
            std::cerr << std::endl;
            return false;
        }
        return true;
    }
}  // namespace ClassTest


// This macro expands to startup code for a test file.
#define SEQAN_BEGIN_TESTSUITE(suite_name)       \
    int main(int argc, char **argv) {           \
    ::ClassTest::beginTestSuite(#suite_name);


// This macro expands to shutdown code for a test file.
#define SEQAN_END_TESTSUITE                     \
    return ::ClassTest::endTestSuite();         \
    }


// This macro expands to function header for one test.
#define SEQAN_DEFINE_TEST(test_name)            \
    void SEQAN_TEST_ ## test_name ()            \


// This macro expands to code to call a given test.
#define SEQAN_CALL_TEST(test_name)              \
    do {                                        \
        ::ClassTest::beginTest(#test_name);     \
        SEQAN_TEST_ ## test_name();             \
        ::ClassTest::endTest();                 \
    } while (false)


// This macro returns from the current function and logs a "skipped"
// event for the current test.
#define SEQAN_SKIP_TEST                         \
    do {                                        \
        ::ClassTest::skipCurrentTest();         \
        return;                                 \
    } while (false)


// Equality assertion with an optional comment.
//
// Usage:  SEQAN_ASSERT_EQ(4, 4);
// Usage:  SEQAN_ASSERT_EQ(4, 5, "Wheee...");
#define SEQAN_ASSERT_EQ(_arg1, _arg2, ...)                              \
    do {                                                                \
        if (not ::ClassTest::testEqual(__FILE__, __LINE__,              \
                                       (_arg1), #_arg1,                 \
                                       (_arg2), #_arg2,                 \
                                       ## __VA_ARGS__)) {               \
        }                                                               \
    } while (false)


// Inequality assertion with an optional comment.
//
// Usage:  SEQAN_ASSERT_NEQ(4, 5);
// Usage:  SEQAN_ASSERT_NEQ(4, 4, "Wheee...");
#define SEQAN_ASSERT_NEQ(_arg1, _arg2, ...)                             \
    do {                                                                \
        if (not ::ClassTest::testNotEqual(__FILE__, __LINE__,           \
                                          (_arg1), #_arg1,              \
                                          (_arg2), #_arg2,              \
                                          ## __VA_ARGS__)) {            \
        }                                                               \
    } while (false)


// TODO(holtgrew): Rename to SEQAN_TASSERT once that name is free.
// Trueness assertion with an optional comment.
//
// Usage:  SEQAN_ASSERT_TRUE(true, "Yay!");
// Usage:  SEQAN_ASSERT_TRUE(false);
#define SEQAN_ASSERT_TRUE(_arg1, ...)                         \
    do {                                                      \
        if (not ::ClassTest::testTrue(__FILE__, __LINE__,     \
                                      (_arg1), #_arg1,        \
                                      ##__VA_ARGS__)) {       \
        }                                                     \
    } while (false)


// Falseness assertion with an optional comment.
//
// Usage:  SEQAN_ASSERT_NOT(true, "Yay!");
// Usage:  SEQAN_ASSERT_NOT(false);
#define SEQAN_ASSERT_NOT(_arg1, ...)                           \
    do {                                                       \
        if (not ::ClassTest::testFalse(__FILE__, __LINE__,     \
                                       (_arg1), #_arg1,        \
                                       ##__VA_ARGS__)) {       \
        }                                                      \
    } while (false)

}  // namespace seqan

#endif  // SEQAN_BASIC_BASIC_TESTING_H_
