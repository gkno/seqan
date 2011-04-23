// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef TESTS_BASIC_TEST_BASIC_ITERATOR_H_
#define TESTS_BASIC_TEST_BASIC_ITERATOR_H_

// --------------------------------------------------------------------------
// Tests for Pointer Adaption to Iterator Concept
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_metafunctions)
{
    using namespace seqan;
    
    // Pointers.
    {
        typedef int * TIterator;
        
        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        
        b = IsSameType<typename Value<TIterator>::Type, int>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int &>::VALUE;
        SEQAN_ASSERT(b);
    }
    // Const-Pointers.
    {
        typedef int const * TIterator;
        
        bool b = IsSameType<typename Difference<TIterator>::Type, ptrdiff_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Position<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Size<TIterator>::Type, size_t>::VALUE;
        SEQAN_ASSERT(b);
        
        // TODO(holtgrew): This is inconsistent
        b = IsSameType<typename Value<TIterator>::Type, int const>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename GetValue<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
        b = IsSameType<typename Reference<TIterator>::Type, int const &>::VALUE;
        SEQAN_ASSERT(b);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_transport)
{
    using namespace seqan;
    
    // assign()
    {
        int x = 1, y = 2;
        int * ptr = &x;
        assign(ptr, &y);
        SEQAN_ASSERT_EQ(ptr, &y);
    }
    // move()
    {
        int x = 1, y = 2;
        int * ptr = &x;
        move(ptr, &y);
        SEQAN_ASSERT_EQ(ptr, &y);
    }
    // set()
    {
        int x = 1, y = 2;
        int * ptr = &x;
        set(ptr, &y);
        SEQAN_ASSERT_EQ(ptr, &y);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_transport_value)
{
    using namespace seqan;
    
    // assignValue()
    {
        CDStruct cs1, cs2;
        resetCDStructStatics();

        CDStruct * ptr = &cs1;
        assignValue(ptr, cs2);

        SEQAN_ASSERT_EQ(ptr->copiedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->movedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->setFrom, -1);
        SEQAN_ASSERT_EQ(ptr->assignedFrom, cs2.id);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &cs2);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 1);
    }
    // moveValue()
    {
        CDStruct cs1, cs2;
        resetCDStructStatics();

        CDStruct * ptr = &cs1;
        moveValue(ptr, cs2);

        SEQAN_ASSERT_EQ(ptr->copiedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->movedFrom, cs2.id);
        SEQAN_ASSERT_EQ(ptr->setFrom, -1);
        SEQAN_ASSERT_EQ(ptr->assignedFrom, -1);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, &cs2);
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 1);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
    // setValue()
    {
        CDStruct cs1, cs2;
        resetCDStructStatics();

        CDStruct * ptr = &cs1;
        setValue(ptr, cs2);

        SEQAN_ASSERT_EQ(ptr->copiedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->movedFrom, -1);
        SEQAN_ASSERT_EQ(ptr->setFrom, -1);
        SEQAN_ASSERT_EQ(ptr->assignedFrom, -1);
        
        SEQAN_ASSERT_EQ(ptr, &cs2);

        SEQAN_ASSERT_EQ(CDStruct::lastOther, static_cast<CDStruct *>(0));
        SEQAN_ASSERT_EQ(CDStruct::defaultConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::copyConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moveConstructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::moves, 0);
        SEQAN_ASSERT_EQ(CDStruct::destructions, 0);
        SEQAN_ASSERT_EQ(CDStruct::assignments, 0);
        SEQAN_ASSERT_EQ(CDStruct::sets, 0);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_movement)
{
    using namespace seqan;
    
    // goNext/operator++
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        int * ptr = &arr[4];
        goNext(ptr);
        SEQAN_ASSERT_EQ(ptr, &arr[5]);
        
        ptr = &arr[4];
        ptr++;
        SEQAN_ASSERT_EQ(ptr, &arr[5]);

        ptr = &arr[4];
        ++ptr;
        SEQAN_ASSERT_EQ(ptr, &arr[5]);
    }
    // goPrevious/operator--
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        int * ptr = &arr[4];
        goPrevious(ptr);
        SEQAN_ASSERT_EQ(ptr, &arr[3]);
        
        ptr = &arr[4];
        ptr--;
        SEQAN_ASSERT_EQ(ptr, &arr[3]);

        ptr = &arr[4];
        --ptr;
        SEQAN_ASSERT_EQ(ptr, &arr[3]);
    }
    // goFurther/operator+=/operator-=
    {
        int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        int * ptr = &arr[4];
        goFurther(ptr, 2);
        SEQAN_ASSERT_EQ(ptr, &arr[6]);

        ptr = &arr[4];
        goFurther(ptr, -2);
        SEQAN_ASSERT_EQ(ptr, &arr[2]);
        
        ptr = &arr[4];
        ptr += 2;
        SEQAN_ASSERT_EQ(ptr, &arr[6]);

        ptr = &arr[4];
        ptr -= 2;
        SEQAN_ASSERT_EQ(ptr, &arr[2]);
    }
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_pointer_arithmetics)
{
    using namespace seqan;
    
    int arr[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    int * ptr = &arr[4];

    int * ptr2 = ptr + 2;
    SEQAN_ASSERT_EQ(ptr2, &arr[6]);

    ptr2 = ptr - 2;
    SEQAN_ASSERT_EQ(ptr2, &arr[2]);
    
    ptr2 = ptr + 2;
    SEQAN_ASSERT_EQ(ptr2 - ptr, 2);
}

// --------------------------------------------------------------------------
// Tests for STL Iterator Adaption to Iterator Concept
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_metafunctions)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_constructors)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_transport)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_transport_value)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_movement)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adapt_std_iterator_arithmetics)
{
    using namespace seqan;
}

// --------------------------------------------------------------------------
// Tests for Adaptor Iterator
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_metafunctions)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_constructors)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_transport)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_transport_value)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_movement)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_arithmetics)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_rooted_metafunctions)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_adaptor_rooted_functions)
{
    using namespace seqan;
}

// --------------------------------------------------------------------------
// Tests for Positional Iterator
// --------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_basic_iterator_position_metafunctions)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_constructors)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_transport)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_transport_value)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_movement)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_arithmetics)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_rooted_metafunctions)
{
    using namespace seqan;
}

SEQAN_DEFINE_TEST(test_basic_iterator_position_rooted_functions)
{
    using namespace seqan;
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_ITERATOR_H_
