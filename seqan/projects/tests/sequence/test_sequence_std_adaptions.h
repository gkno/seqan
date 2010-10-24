/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ============================================================================
  Tests for the STL adaptions of seqan.
  ==========================================================================*/

// TODO(holtgrew): Split into a file for STL strings and one for STL lists?  Writing templatized tests do not make too much sense, I guess, becauses lists and strings are so dissimilar.

#ifndef TEST_SEQUENCE_TEST_SEQUENCE_STD_ADAPTIONS_H_
#define TEST_SEQUENCE_TEST_SEQUENCE_STD_ADAPTIONS_H_

// Tests the return types and existence of the metafunctions for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_metafunctions_std_string)
{
    using namespace seqan;
    
    typedef int TElement;
    typedef std::basic_string<TElement> TString;
    typedef TString const TConstString;

    // Test IsContiguous<>::VALUE
    {
        bool b = IsContiguous<TString>::VALUE;
        SEQAN_ASSERT_NOT(b);
        b = IsContiguous<TConstString>::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    // Test Value<>::VALUE
    {
        typedef Value<TString>::Type TValue;
        bool b = TYPECMP<TValue, TElement>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Value<TConstString>::Type TConstValue;
        b = TYPECMP<TConstValue, TElement>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test GetValue<>::VALUE
    {
        typedef GetValue<TString>::Type TGetValue;
        bool b = TYPECMP<TGetValue, TElement &>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef GetValue<TConstString>::Type TConstGetValue;
        b = TYPECMP<TConstGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test GetReference<>::VALUE
    {
        typedef Reference<TString>::Type TReference;
        bool b = TYPECMP<TReference, TElement &>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Reference<TConstString>::Type TConstReference;
        b = TYPECMP<TConstReference, TElement const &>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test Iterator<, Rooted>::VALUE
    {
        typedef Iterator<TString, Rooted>::Type TIterator;
        typedef Iter<TString, AdaptorIterator<Iter<TString, StdIteratorAdaptor> > > TExpected;
        bool b = TYPECMP<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Iterator<TConstString, Rooted>::Type TConstIterator;
        typedef Iter<TConstString, AdaptorIterator<Iter<TConstString, StdIteratorAdaptor> > > TExpectedConst;
        b = TYPECMP<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test Iterator<, Standard>::VALUE
    {
        typedef Iterator<TString, Standard>::Type TIterator;
        typedef Iter<TString, StdIteratorAdaptor> TExpected;
        bool b = TYPECMP<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Iterator<TConstString, Standard>::Type TConstIterator;
        typedef Iter<TConstString, StdIteratorAdaptor> TExpectedConst;
        b = TYPECMP<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test Position<>::VALUE
    {
        typedef Position<TString>::Type TPosition;
        bool b = TYPECMP<TPosition, TString::size_type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Position<TConstString>::Type TConstPosition;
        b = TYPECMP<TConstPosition, TString::size_type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test Size<>::VALUE
    {
        typedef Size<TString>::Type TPosition;
        bool b = TYPECMP<TPosition, TString::size_type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Size<TConstString>::Type TConstPosition;
        b = TYPECMP<TConstPosition, TString::size_type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
}


// Test iterators for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_iterators_std_string)
{
    using namespace seqan;
    
    // Test const iterator.
    {
        std::string const str = "Unimportant contents.";
        typedef Iterator<std::string const>::Type TIterator;

        std::string strCopy;
        for (TIterator it = begin(str, Standard()); it != end(str, Standard()); ++it)
            appendValue(strCopy, value(it));

        SEQAN_ASSERT_EQ(str, strCopy);
    }

    // Test non-const iterator.
    {
        std::string str = "Unimportant contents.";
        typedef Iterator<std::string>::Type TIterator;

        std::string strCopy;
        for (TIterator it = begin(str, Standard()); it != end(str, Standard()); ++it)
            appendValue(strCopy, value(it));

        SEQAN_ASSERT_EQ(str, strCopy);
    }
}


// Tests for the basic sequence functions for STL strings,
// e.g. value(), front(), back().
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_interface_std_string)
{
    using namespace seqan;

    std::string str = "Hello World!";

    // value(str, i), getValue(str, i)
    SEQAN_ASSERT_EQ(value(str, 0), 'H');
    SEQAN_ASSERT_EQ(value(str, 4), 'o');
    SEQAN_ASSERT_EQ(getValue(str, 0), 'H');
    SEQAN_ASSERT_EQ(getValue(str, 4), 'o');

    // front(), back()
    SEQAN_ASSERT_EQ(front(str), 'H');
    SEQAN_ASSERT_EQ(back(str), '!');

    // length()
    SEQAN_ASSERT_EQ(length(str), 12u);

    // TODO(holtgrew): Anything else missing? Probably...
}


// Tests for the memory allocation and reservation related functions
// for STL strings.
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_memory_std_string)
{
    using namespace seqan;

    // Test resize function -- resize down.
    {
        std::string str = "Hello world!";
        resize(str, 5);
        SEQAN_ASSERT_EQ(str, "Hello");
    }

    // Test resize function -- resize up.
    {
        std::string str = "12345";
        resize(str, 6);
        // The following gives an assertion in positional setValue() if not resized properly.
        str[5] = '6';
        SEQAN_ASSERT_EQ(str, "123456");
    }

    // Tests reserve function.
    {
        std::string str;
        reserve(str, 10);
        SEQAN_ASSERT_GEQ(capacity(str), 10u);
    }
    {
        std::string str;
        reserve(str, 10, Generous());
        SEQAN_ASSERT_GEQ(capacity(str), 10u);
    }
    {
        std::string str;
        reserve(str, 10, Exact());
        SEQAN_ASSERT_EQ(capacity(str), 10u);
    }
}


// Tests the return types and existence of the metafunctions for STL lists.
SEQAN_DEFINE_TEST(test_sequence_adaptions_metafunctions_std_list)
{
    using namespace seqan;
    
    typedef int TElement;
    typedef std::list<TElement> TList;
    typedef TList const TConstList;

    // Test IsContiguous<>::VALUE
    {
        bool b = IsContiguous<TList>::VALUE;
        SEQAN_ASSERT_NOT(b);
        b = IsContiguous<TConstList>::VALUE;
        SEQAN_ASSERT_NOT(b);
    }
    // Test Value<>::VALUE
    {
        typedef Value<TList>::Type TValue;
        bool b = TYPECMP<TValue, TElement>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Value<TConstList>::Type TConstValue;
        b = TYPECMP<TConstValue, TElement>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test GetValue<>::VALUE
    {
        typedef GetValue<TList>::Type TGetValue;
        bool b = TYPECMP<TGetValue, TElement &>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef GetValue<TConstList>::Type TConstGetValue;
        b = TYPECMP<TConstGetValue, TElement const &>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test GetReference<>::VALUE
    {
        typedef Reference<TList>::Type TReference;
        bool b = TYPECMP<TReference, TElement &>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Reference<TConstList>::Type TConstReference;
        b = TYPECMP<TConstReference, TElement const &>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test Iterator<, Rooted>::VALUE
    {
        typedef Iterator<TList, Rooted>::Type TIterator;
        typedef Iter<TList, AdaptorIterator<Iter<TList, StdIteratorAdaptor> > > TExpected;
        bool b = TYPECMP<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Iterator<TConstList, Rooted>::Type TConstIterator;
        typedef Iter<TConstList, AdaptorIterator<Iter<TConstList, StdIteratorAdaptor> > > TExpectedConst;
        b = TYPECMP<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test Iterator<, Standard>::VALUE
    {
        typedef Iterator<TList, Standard>::Type TIterator;
        typedef Iter<TList, StdIteratorAdaptor> TExpected;
        bool b = TYPECMP<TIterator, TExpected>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Iterator<TConstList, Standard>::Type TConstIterator;
        typedef Iter<TConstList, StdIteratorAdaptor> TExpectedConst;
        b = TYPECMP<TConstIterator, TExpectedConst>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test Position<>::VALUE
    {
        typedef Position<TList>::Type TPosition;
        bool b = TYPECMP<TPosition, TList::size_type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Position<TConstList>::Type TConstPosition;
        b = TYPECMP<TConstPosition, TList::size_type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
    // Test Size<>::VALUE
    {
        typedef Size<TList>::Type TPosition;
        bool b = TYPECMP<TPosition, TList::size_type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
        typedef Size<TConstList>::Type TConstPosition;
        b = TYPECMP<TConstPosition, TList::size_type>::VALUE;
        SEQAN_ASSERT_TRUE(b);
    }
}


// Test iterators for STL lists.
SEQAN_DEFINE_TEST(test_sequence_adaptions_iterators_std_list)
{
    using namespace seqan;
    
    typedef int TElement;
    
    // Test Standard, non-const iterators.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList, Standard>::Type TIterator;

        TList list;
        appendValue(list, 1);
        appendValue(list, 2);
        appendValue(list, 3);

        // The command sequence in the following is a bit arbitrary
        // but should robustly test that the iterators work correctly.
        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        --it;
        SEQAN_ASSERT_EQ(1, *it);
        TIterator itEnd = seqan::end(list);
        SEQAN_ASSERT_NOT(it == itEnd);
        ++it;
        ++it;
        ++it;
        SEQAN_ASSERT_TRUE(it == itEnd);

        // The following does not apply to const iterators.
        it = seqan::begin(list);
        *seqan::begin(list) = 4;
        SEQAN_ASSERT_EQ(4, *it);
    }
    // Test Standard, const iterators.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList const, Standard>::Type TIterator;

        TList mutableList;
        appendValue(mutableList, 1);
        appendValue(mutableList, 2);
        appendValue(mutableList, 3);

        TList const & list = mutableList;

        // The command sequence in the following is a bit arbitrary
        // but should robustly test that the iterators work correctly.
        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        --it;
        SEQAN_ASSERT_EQ(1, *it);
        TIterator itEnd = seqan::end(list);
        SEQAN_ASSERT_NOT(it == itEnd);
        ++it;
        ++it;
        ++it;
        SEQAN_ASSERT_TRUE(it == itEnd);
    }
    // Test Rooted, non-const iterators.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList, Rooted>::Type TIterator;

        TList list;
        appendValue(list, 1);
        appendValue(list, 2);
        appendValue(list, 3);

        // The command sequence in the following is a bit arbitrary
        // but should robustly test that the iterators work correctly.
        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        --it;
        SEQAN_ASSERT_EQ(1, *it);
        TIterator itEnd = seqan::end(list);
        SEQAN_ASSERT_NOT(it == itEnd);
        ++it;
        ++it;
        ++it;
        SEQAN_ASSERT_TRUE(it == itEnd);

        // The following does not apply to const iterators.
        it = seqan::begin(list);
        *seqan::begin(list) = 4;
        SEQAN_ASSERT_EQ(4, *it);
    }
    // Test Rooted, const iterators.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList const, Rooted>::Type TIterator;

        TList mutableList;
        appendValue(mutableList, 1);
        appendValue(mutableList, 2);
        appendValue(mutableList, 3);

        TList const & list = mutableList;

        // The command sequence in the following is a bit arbitrary
        // but should robustly test that the iterators work correctly.
        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        --it;
        SEQAN_ASSERT_EQ(1, *it);
        TIterator itEnd = seqan::end(list);
        SEQAN_ASSERT_NOT(it == itEnd);
        ++it;
        ++it;
        ++it;
        SEQAN_ASSERT_TRUE(it == itEnd);
    }
}


// Test the basic sequence interface implemented for STL list, e.g. front() and back().
SEQAN_DEFINE_TEST(test_sequence_adaptions_sequence_interface_std_list)
{
    using namespace seqan;
    
    typedef int TElement;
    
    // Test with non-const container.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList>::Type TIterator;

        // Prepare list...
        TList list;
        appendValue(list, 1);
        appendValue(list, 2);
        appendValue(list, 3);

        // Test reading front and back.
        SEQAN_ASSERT_EQ(1, front(list));
        SEQAN_ASSERT_EQ(3, back(list));

        // Test assigning to front and back.
        front(list) = -1;
        back(list) = -3;

        TIterator it = seqan::begin(list);
        SEQAN_ASSERT_EQ(-1, *it);
        ++it;
        SEQAN_ASSERT_EQ(2, *it);
        ++it;
        SEQAN_ASSERT_EQ(-3, *it);

        // Test appending and prepending values.
        prependValue(list, 42);
        appendValue(list, 43);
        SEQAN_ASSERT_EQ(42, front(list));
        SEQAN_ASSERT_EQ(43, back(list));

        // Test length().
        SEQAN_ASSERT_EQ(5u, length(list));

        // Test clear().
        {
            TList listCopy(list);
            SEQAN_ASSERT_EQ(5u, length(listCopy));
            clear(listCopy);
            SEQAN_ASSERT_EQ(0u, length(listCopy));
        }
    }
    // Test with const container.
    {
        typedef std::list<TElement> TList;
        typedef Iterator<TList>::Type TIterator;

        // Prepare list...
        TList mutableList;        
        appendValue(mutableList, 1);
        appendValue(mutableList, 2);
        appendValue(mutableList, 3);

        TList const & list = mutableList;

        // Test reading front and back.
        SEQAN_ASSERT_EQ(1, front(list));
        SEQAN_ASSERT_EQ(3, back(list));

        // Test length().
        SEQAN_ASSERT_EQ(3u, length(list));
    }
}

#endif  // TEST_SEQUENCE_TEST_SEQUENCE_STD_ADAPTIONS_H_
