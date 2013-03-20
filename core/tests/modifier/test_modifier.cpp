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

#include <seqan/file.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>

SEQAN_DEFINE_TEST(test_modifier_modified_string_metafunctions)
{
    typedef seqan::Dna5String                           TString;
    typedef seqan::ModifiedString<TString>              TInnerModifiedString;
    typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

    SEQAN_ASSERT(+(seqan::IsSameType<typename seqan::Pointer_<TString>::Type, TString *>::VALUE));
    SEQAN_ASSERT(+(seqan::IsSameType<typename seqan::Pointer_<TInnerModifiedString>::Type, TInnerModifiedString>::VALUE));
    SEQAN_ASSERT(+(seqan::IsSameType<typename seqan::Pointer_<TOuterModifiedString>::Type, TOuterModifiedString>::VALUE));
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_construct)
{
    typedef seqan::Dna5String                           TString;
    typedef seqan::ModifiedString<TString>              TInnerModifiedString;
    typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

    // Default Construction with one level.
    {
        TInnerModifiedString modified;
    }

    // Construct with underlying string and one level.
    {
        TString original;
        TInnerModifiedString inner(original);
    }

    // Construct with underlying string and two levels.
    {
        TString original;
        TInnerModifiedString inner(original);
        TOuterModifiedString outer(inner);
    }

    // Construct with underlying string and two levels, direct construction.
    {
        TString original;
        TOuterModifiedString outer((TInnerModifiedString(original)));
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_assignment)
{
    typedef seqan::Dna5String                           TString;
    typedef seqan::ModifiedString<TString>              TInnerModifiedString;
    typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

    // Copy with one level.
    {
        TString original = "CGAT";
        TInnerModifiedString  inner(original);

        TInnerModifiedString inner2;
        inner2 = inner;
    }

    // Copy with two levels.
    {
        TString original = "CGAT";
        TInnerModifiedString  inner(original);

        TInnerModifiedString inner2;
        inner2 = inner;
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_length)
{
    typedef seqan::Dna5String                           TString;
    typedef seqan::ModifiedString<TString>              TInnerModifiedString;
    typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;

    // Default Construction with one level.
    {
        TInnerModifiedString modified;
    }

    // Construct with underlying string and one level.
    {
        TString original = "CGAT";
        TInnerModifiedString inner(original);

        SEQAN_ASSERT_EQ(length(original), 4u);
        SEQAN_ASSERT_EQ(length(inner), 4u);
    }

    // Construct with underlying string and two levels.
    {
        TString original = "CGAT";
        TInnerModifiedString inner(original);
        TOuterModifiedString outer(inner);

        SEQAN_ASSERT_EQ(length(original), 4u);
        SEQAN_ASSERT_EQ(length(inner), 4u);
        SEQAN_ASSERT_EQ(length(outer), 4u);
    }

    // Construct with underlying string and two levels, direct construction.
    {
        TString original = "CGAT";
        TOuterModifiedString outer((TInnerModifiedString(original)));

        SEQAN_ASSERT_EQ(length(original), 4u);
        SEQAN_ASSERT_EQ(length(host(outer)), 4u);
        SEQAN_ASSERT_EQ(length(outer), 4u);
    }
}

// Construct modified string cascade from innermost host.
SEQAN_DEFINE_TEST(test_modifier_modified_string_cascade)
{
    // Host is a string.
    {
        typedef seqan::Dna5String                           TString;
        typedef seqan::ModifiedString<TString>              TInnerModifiedString;
        typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;
        
        TString original = "CGAT";
        TInnerModifiedString inner(original);
        TOuterModifiedString outer(original);
    }

    // Host is a segment.
    {
        typedef seqan::Dna5String                           TString;
        typedef typename seqan::Infix<TString>::Type        TInfix;
        typedef seqan::ModifiedString<TInfix>               TInnerModifiedString;
        typedef seqan::ModifiedString<TInnerModifiedString> TOuterModifiedString;
        
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TInnerModifiedString inner(origInfix);
        TOuterModifiedString outer(origInfix);
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_iterator_construct)
{
    typedef seqan::Dna5String                                 TString;
    typedef typename seqan::Iterator<seqan::Dna5String>::Type TIter;
    typedef seqan::ModifiedIterator<TIter>                    TModifiedIter;

    // Default construction.
    {
        TModifiedIter itM;
    }

    // Construction from host iterator.
    {
        TString seq = "CGAT";
        TIter it = begin(seq);
        TModifiedIter itM(it);
    }
}

struct LowerFunctor : std::unary_function<char, char>
{
    char operator()(char c) const
    {
        return tolower(c);
    }
};

struct CaesarFunctor : std::unary_function<char, char>
{
    char operator()(char c) const
    {
        int x = c;
        x += 3;
        return x % 256;
    }
};

SEQAN_DEFINE_TEST(test_modifier_modified_string_mod_view)
{
    typedef seqan::CharString TString;
    typedef seqan::ModView<LowerFunctor> TModViewLower;
    typedef seqan::ModifiedString<TString, TModViewLower> TInnerModifiedString;
    typedef seqan::ModView<CaesarFunctor> TModViewCaesar;
    typedef seqan::ModifiedString<TInnerModifiedString, TModViewCaesar> TOuterModifiedString;

    // One level only.
    {
        TString original = "CGAT";
        TInnerModifiedString modified(original);

        SEQAN_ASSERT_EQ(modified, "cgat");
    }
    // Two levels.
    {
        TString original = "CGAT";
        TInnerModifiedString modified(original);
        TOuterModifiedString outer(modified);

        SEQAN_ASSERT_EQ(outer, "fjdw");
    }
    // Two levels directly constructed with original.
    {
        TString original = "CGAT";
        TOuterModifiedString outer(original);

        SEQAN_ASSERT_EQ(outer, "fjdw");
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_mod_view_segment)
{
    typedef seqan::CharString TString;
    typedef seqan::Infix<TString>::Type TInfix;
    typedef seqan::ModView<LowerFunctor> TModViewLower;
    typedef seqan::ModifiedString<TInfix, TModViewLower> TInnerModifiedString;
    typedef seqan::ModView<CaesarFunctor> TModViewCaesar;
    typedef seqan::ModifiedString<TInnerModifiedString, TModViewCaesar> TOuterModifiedString;

    // One level only.
    {
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TInnerModifiedString modified(origInfix);

        SEQAN_ASSERT_EQ(modified, "ga");
    }
    // Two levels.
    {
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TInnerModifiedString modified(origInfix);
        TOuterModifiedString outer(modified);

        SEQAN_ASSERT_EQ(outer, "jd");
    }
    // Two levels directly constructed with original.
    {
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TOuterModifiedString outer(origInfix);

        SEQAN_ASSERT_EQ(outer, "jd");
    }
}

SEQAN_DEFINE_TEST(test_modifier_modified_string_reverse_segment)
{
    // Inner is lower, outer is reverse.
    {
        typedef seqan::CharString TString;
        typedef seqan::Infix<TString>::Type TInfix;
        
        typedef seqan::ModView<LowerFunctor> TModViewLower;
        typedef seqan::ModifiedString<TInfix, TModViewLower> TInnerModifiedString;
        
        typedef seqan::ModifiedString<TInnerModifiedString, seqan::ModReverse> TOuterModifiedString;
        
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TOuterModifiedString outer(origInfix);
        SEQAN_ASSERT_EQ(outer, "ag");
    }
    // Inner is reverse, outer is lower.
    {
        typedef seqan::CharString TString;
        typedef seqan::Infix<TString>::Type TInfix;
        
        typedef seqan::ModifiedString<TInfix, seqan::ModReverse> TInnerModifiedString;

        typedef seqan::ModView<LowerFunctor> TModViewLower;
        typedef seqan::ModifiedString<TInnerModifiedString, TModViewLower> TOuterModifiedString;
        
        TString original = "CGAT";
        TInfix origInfix = infix(original, 1, 3);
        TOuterModifiedString outer(origInfix);
        SEQAN_ASSERT_EQ(outer, "ag");
    }
}

SEQAN_DEFINE_TEST(test_modifier_minimal)
{
    typedef seqan::CharString TString;
    typedef seqan::ModifiedString<TString, seqan::ModReverse> TInnerModifiedString;

    TString original = "The QUICK brown fox.";
    TInnerModifiedString inner(original);
    // static_cast<typename seqan::Parameter_<TString>::Type *>(seqan::Nothing());  // #=> non-const reference
    // static_cast<typename TInnerModifiedString::THostPointer_ *>(seqan::Nothing());  // #=> const pointer
    // static_cast<typename seqan::Pointer_<TString>::Type *>(seqan::Nothing());
}

SEQAN_BEGIN_TESTSUITE(test_modifier) 
{
    SEQAN_CALL_TEST(test_modifier_modified_string_metafunctions);
    SEQAN_CALL_TEST(test_modifier_modified_string_construct);
    SEQAN_CALL_TEST(test_modifier_modified_string_assignment);
    SEQAN_CALL_TEST(test_modifier_modified_string_length);
    SEQAN_CALL_TEST(test_modifier_modified_string_cascade);

    // TODO(holtgrew): SEQAN_CALL_TEST(test_modifier_modified_iterator_metafunctions);
    SEQAN_CALL_TEST(test_modifier_modified_iterator_construct);

    SEQAN_CALL_TEST(test_modifier_modified_string_mod_view);
    SEQAN_CALL_TEST(test_modifier_modified_string_mod_view_segment);

    SEQAN_CALL_TEST(test_modifier_modified_string_reverse_segment);

    SEQAN_CALL_TEST(test_modifier_minimal);
}
SEQAN_END_TESTSUITE
