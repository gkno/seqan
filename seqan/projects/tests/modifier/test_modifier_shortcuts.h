/*
  Tests for modifier/modifier_shortcuts.h
*/

#ifndef TEST_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_
#define TEST_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_reverse) {
    DnaString str = "CGAT";
    DnaString kExpectedStr = "TAGC";

    DnaStringReverseComplement modifiedString(str);
    DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_reverse) {
    Dna5String str = "CGATN";
    Dna5String kExpectedStr = "NTAGC";

    Dna5StringReverseComplement modifiedString(str);
    Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_reverse) {
    RnaString str = "CGAU";
    RnaString kExpectedStr = "UAGC";

    RnaStringReverseComplement modifiedString(str);
    RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_reverse) {
    Rna5String str = "CGAUN";
    Rna5String kExpectedStr = "NUAGC";

    Rna5StringReverseComplement modifiedString(str);
    Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_complement) {
    DnaString str = "CGAT";
    DnaString kExpectedStr = "GCTA";

    DnaStringReverseComplement modifiedString(str);
    DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_complement) {
    Dna5String str = "CGATN";
    Dna5String kExpectedStr = "GCTAN";

    Dna5StringReverseComplement modifiedString(str);
    Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_complement) {
    RnaString str = "CGAU";
    RnaString kExpectedStr = "GCUA";

    RnaStringReverseComplement modifiedString(str);
    RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_complement) {
    Rna5String str = "CGAUN";
    Rna5String kExpectedStr = "GCUAN";

    Rna5StringReverseComplement modifiedString(str);
    Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna_string_reverse_complement) {
    DnaString str = "CGAT";
    DnaString kExpectedStr = "ATCG";

    DnaStringReverseComplement modifiedString(str);
    DnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_dna5_string_reverse_complement) {
    Dna5String str = "CGATN";
    Dna5String kExpectedStr = "NATCG";

    Dna5StringReverseComplement modifiedString(str);
    Dna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna_string_reverse_complement) {
    RnaString str = "CGAU";
    RnaString kExpectedStr = "AUCG";

    RnaStringReverseComplement modifiedString(str);
    RnaString modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_rna5_string_reverse_complement) {
    Rna5String str = "CGAUN";
    Rna5String kExpectedStr = "NAUCG";

    Rna5StringReverseComplement modifiedString(str);
    Rna5String modifiedStringCopy = modifiedString;
    SEQAN_ASSERT_EQ(kExpectedStr, modifiedStringCopy);
}


// TODO(holtgrew): The following could be made non-redundant with a helper function.


SEQAN_DEFINE_TEST(test_modifer_shortcuts_complement_in_place_string) {
    Dna5String str = "CGATN";
    Dna5String const kStr = str;
    Dna5String const kExpectedResult = "GCTAN";

    // Test non-const version.
    complementInPlace(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    complementInPlace(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_complement_in_place_string_set) {
    Dna5String str1 = "CCGGTTAANN";
    Dna5String str2 = "CGTANCGTAN";
    Dna5String const kExpectedStr1 = "GGCCAATTNN";
    Dna5String const kExpectedStr2 = "GCATNGCATN";

    // Test non-const version.
    StringSet<Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    complementInPlace(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, strSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, strSet[1]);

    // Test const version.
    StringSet<Dna5String> const kStrSet = strSet;

    complementInPlace(kStrSet);
    SEQAN_ASSERT_EQ(2u, length(kStrSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, kStrSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, kStrSet[1]);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_in_place_string) {
    Dna5String str = "CGATN";
    Dna5String const kStr = str;
    Dna5String const kExpectedResult = "NATCG";

    // Test non-const version.
    reverseComplementInPlace(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    reverseComplementInPlace(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_complement_in_place_string_set) {
    Dna5String str1 = "CCGGTTAANN";
    Dna5String str2 = "CGTANCGTAN";
    Dna5String const kExpectedStr1 = "NNTTAACCGG";
    Dna5String const kExpectedStr2 = "NTACGNTACG";

    // Test non-const version.
    StringSet<Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    reverseComplementInPlace(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, strSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, strSet[1]);

    // Test const version.
    StringSet<Dna5String> const kStrSet = strSet;

    reverseInPlace(kStrSet);
    SEQAN_ASSERT_EQ(2u, length(kStrSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, kStrSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, kStrSet[1]);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_in_place_string) {
    Dna5String str = "CGATN";
    Dna5String const kStr = str;
    Dna5String const kExpectedResult = "NTAGC";

    // Test non-const version.
    reverseInPlace(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    reverseInPlace(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_reverse_in_place_string_set) {
    Dna5String str1 = "CCGGTTAANN";
    Dna5String str2 = "CGTANCGTAN";
    Dna5String const kExpectedStr1 = "NNAATTGGCC";
    Dna5String const kExpectedStr2 = "NATGCNATGC";

    // Test non-const version.
    StringSet<Dna5String> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    reverseInPlace(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, strSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, strSet[1]);

    // Test const version.
    StringSet<Dna5String> const kStrSet = strSet;

    reverseInPlace(kStrSet);
    SEQAN_ASSERT_EQ(2u, length(kStrSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, kStrSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, kStrSet[1]);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_lower_in_place_string) {
    CharString str = "This is a test!";
    CharString const kStr = str;
    CharString const kExpectedResult = "this is a test!";

    // Test non-const version.
    toLowerInPlace(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    toLowerInPlace(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_lower_in_place_string_set) {
    CharString str1 = "This is a test!";
    CharString str2 = "This is also a test!";
    CharString const kExpectedStr1 = "this is a test!";
    CharString const kExpectedStr2 = "this is also a test!";

    // Test non-const version.
    StringSet<CharString> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    toLowerInPlace(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, strSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, strSet[1]);

    // Test const version.
    StringSet<CharString> const kStrSet = strSet;

    toLowerInPlace(kStrSet);
    SEQAN_ASSERT_EQ(2u, length(kStrSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, kStrSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, kStrSet[1]);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_upper_in_place_string) {
    CharString str = "This is a test!";
    CharString const kStr = str;
    CharString const kExpectedResult = "THIS IS A TEST!";

    // Test non-const version.
    toUpperInPlace(str);
    SEQAN_ASSERT_EQ(kExpectedResult, str);

    // Test const version.
    // TODO(holtgrew): This should not be possible!
    toUpperInPlace(kStr);
    SEQAN_ASSERT_EQ(kExpectedResult, kStr);
}


SEQAN_DEFINE_TEST(test_modifer_shortcuts_to_upper_in_place_string_set) {
    CharString str1 = "This is a test!";
    CharString str2 = "This is also a test!";
    CharString const kExpectedStr1 = "THIS IS A TEST!";
    CharString const kExpectedStr2 = "THIS IS ALSO A TEST!";

    // Test non-const version.
    StringSet<CharString> strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);

    toUpperInPlace(strSet);
    SEQAN_ASSERT_EQ(2u, length(strSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, strSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, strSet[1]);

    // Test const version.
    StringSet<CharString> const kStrSet = strSet;

    toUpperInPlace(kStrSet);
    SEQAN_ASSERT_EQ(2u, length(kStrSet));
    SEQAN_ASSERT_EQ(kExpectedStr1, kStrSet[0]);
    SEQAN_ASSERT_EQ(kExpectedStr2, kStrSet[1]);
}



#endif  // TEST_MODIFIER_TEST_MODIFIER_SHORTCUTS_H_
