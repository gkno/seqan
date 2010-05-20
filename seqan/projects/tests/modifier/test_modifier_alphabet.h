#ifndef TESTS_MODIFIER_TEST_MODIFIER_ALPHABET_H_
#define TESTS_MODIFIER_TEST_MODIFIER_ALPHABET_H_

#include <seqan/basic.h>
#include <seqan/modifier.h>

using namespace seqan;

// Test the size metafunctions.
SEQAN_DEFINE_TEST(test_modifier_alphabet_size_metafunctions) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;
    
    // Add special gap symbol to Dna.
    SEQAN_ASSERT_EQ((int)BitsPerValue<TDnaWithGap>::VALUE, 3);
    SEQAN_ASSERT_EQ((int)ValueSize<TDnaWithGap>::VALUE, 5);

    // Add an arbitrary symbol to Dna.
    SEQAN_ASSERT_EQ((int)BitsPerValue<TDnaWithX>::VALUE, 3);
    SEQAN_ASSERT_EQ((int)ValueSize<TDnaWithX>::VALUE, 5);
    
}

SEQAN_DEFINE_TEST(test_modifier_DnaQ) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna5Q, ModExpand<'-'> > TDna5QWithGap;
    typedef ModifiedAlphabet<DnaQ, ModExpand<'X'> > TDnaQWithX;

    // Call convertImpl() directly.
    // TODO(holtgew): Call this?

    // Call convertImpl() through convert().
    // TDnaWithGap;
    SEQAN_ASSERT_EQ(Dna5('A'), convert<Dna5>(TDna5QWithGap('A')));
    SEQAN_ASSERT_EQ(Dna5('C'), convert<Dna5>(TDna5QWithGap('C')));
    SEQAN_ASSERT_EQ(Dna5('G'), convert<Dna5>(TDna5QWithGap('G')));
    SEQAN_ASSERT_EQ(Dna5('T'), convert<Dna5>(TDna5QWithGap('T')));
    SEQAN_ASSERT_EQ(Dna5('N'), convert<Dna5>(Dna5Q('N')));
    SEQAN_ASSERT_EQ(Dna5('N'), convert<Dna5>(TDna5QWithGap('N')));
    SEQAN_ASSERT_EQ(Dna('T'), convert<Dna>(TDna5QWithGap('T')));
    SEQAN_ASSERT_EQ('-', convert<char>(TDna5QWithGap('-')));

    SEQAN_ASSERT_EQ(Dna('A'), convert<Dna>(TDnaQWithX('A')));
    SEQAN_ASSERT_EQ(Dna('C'), convert<Dna>(TDnaQWithX('C')));
    SEQAN_ASSERT_EQ(Dna('G'), convert<Dna>(TDnaQWithX('G')));
    SEQAN_ASSERT_EQ(Dna('T'), convert<Dna>(TDnaQWithX('T')));
    SEQAN_ASSERT_EQ(Dna5('T'), convert<Dna5>(TDnaQWithX('T')));
    SEQAN_ASSERT_EQ('X', convert<char>(TDnaQWithX('X')));

	SEQAN_ASSERT_EQ('A', convert<char>(TDna5QWithGap(0)));
	SEQAN_ASSERT_EQ('C', convert<char>(TDna5QWithGap(1)));
	SEQAN_ASSERT_EQ('G', convert<char>(TDna5QWithGap(2)));
	SEQAN_ASSERT_EQ('T', convert<char>(TDna5QWithGap(3)));
	SEQAN_ASSERT_EQ('N', convert<char>(TDna5QWithGap(4)));
	SEQAN_ASSERT_EQ('-', convert<char>(TDna5QWithGap(5)));
}

SEQAN_DEFINE_TEST(test_modifier_alphabet_enumerate) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna5, ModExpand<'-'> > TDna5WithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;

    // TDnaWithGap;
	SEQAN_ASSERT_EQ('A', convert<char>(TDnaWithGap(0)));
	SEQAN_ASSERT_EQ('C', convert<char>(TDnaWithGap(1)));
	SEQAN_ASSERT_EQ('G', convert<char>(TDnaWithGap(2)));
	SEQAN_ASSERT_EQ('T', convert<char>(TDnaWithGap(3)));
	SEQAN_ASSERT_EQ('-', convert<char>(TDnaWithGap(4)));

    // TDnaWithX;
	SEQAN_ASSERT_EQ('A', convert<char>(TDnaWithX(0)));
	SEQAN_ASSERT_EQ('C', convert<char>(TDnaWithX(1)));
	SEQAN_ASSERT_EQ('G', convert<char>(TDnaWithX(2)));
	SEQAN_ASSERT_EQ('T', convert<char>(TDnaWithX(3)));
	SEQAN_ASSERT_EQ('X', convert<char>(TDnaWithX(4)));

    // TDna5WithGap;
	SEQAN_ASSERT_EQ('A', convert<char>(TDna5WithGap(0)));
	SEQAN_ASSERT_EQ('C', convert<char>(TDna5WithGap(1)));
	SEQAN_ASSERT_EQ('G', convert<char>(TDna5WithGap(2)));
	SEQAN_ASSERT_EQ('T', convert<char>(TDna5WithGap(3)));
	SEQAN_ASSERT_EQ('N', convert<char>(TDna5WithGap(4)));
	SEQAN_ASSERT_EQ('-', convert<char>(TDna5WithGap(5)));
}

SEQAN_DEFINE_TEST(test_modifier_alphabet_convert) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna5, ModExpand<'-'> > TDna5WithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;

    // Call convertImpl() directly.
    // TODO(holtgew): Call this?

    // Call convertImpl() through convert().
    // TDnaWithGap;
    SEQAN_ASSERT_EQ(Dna('A'), convert<Dna>(TDnaWithGap('A')));
    SEQAN_ASSERT_EQ(Dna('C'), convert<Dna>(TDnaWithGap('C')));
    SEQAN_ASSERT_EQ(Dna('G'), convert<Dna>(TDnaWithGap('G')));
    SEQAN_ASSERT_EQ(Dna('T'), convert<Dna>(TDnaWithGap('T')));
    SEQAN_ASSERT_EQ(Dna('A'), convert<Dna>(TDnaWithGap('-')));
    // TDnaWithX
    SEQAN_ASSERT_EQ(Dna('A'), convert<Dna>(TDnaWithX('A')));
    SEQAN_ASSERT_EQ(Dna('C'), convert<Dna>(TDnaWithX('C')));
    SEQAN_ASSERT_EQ(Dna('G'), convert<Dna>(TDnaWithX('G')));
    SEQAN_ASSERT_EQ(Dna('T'), convert<Dna>(TDnaWithX('T')));
    SEQAN_ASSERT_EQ(Dna('A'), convert<Dna>(TDnaWithX('X')));
	
	
	typedef String<TDnaWithGap> TString;

	TString str = "aCgT-AcGt";
//	std::cout << str << std::endl;

	SEQAN_ASSERT_EQ(str, "aCgT-AcGt");
	SEQAN_ASSERT_EQ(str, "AcGt-aCgT");
}


SEQAN_DEFINE_TEST(test_modifier_alphabet_ord_value) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;

    // Call ordValue().
    // TDnaWithGap;
    SEQAN_ASSERT_EQ(0u, ordValue(TDnaWithGap('A')));
    SEQAN_ASSERT_EQ(1u, ordValue(TDnaWithGap('C')));
    SEQAN_ASSERT_EQ(2u, ordValue(TDnaWithGap('G')));
    SEQAN_ASSERT_EQ(3u, ordValue(TDnaWithGap('T')));
    SEQAN_ASSERT_EQ(4u, ordValue(TDnaWithGap('-')));
    // TDnaWithX
    SEQAN_ASSERT_EQ(0u, ordValue(TDnaWithX('A')));
    SEQAN_ASSERT_EQ(1u, ordValue(TDnaWithX('C')));
    SEQAN_ASSERT_EQ(2u, ordValue(TDnaWithX('G')));
    SEQAN_ASSERT_EQ(3u, ordValue(TDnaWithX('T')));
    SEQAN_ASSERT_EQ(4u, ordValue(TDnaWithX('X')));
}


// Test the operator== implementations.
SEQAN_DEFINE_TEST(test_modifier_alphabet_operator_eq) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;
    typedef ModifiedAlphabet<TDnaWithX, ModExpand<'Y'> > TDnaWithXY;

    // TODO(holtgrew): Should that many alphabets be comparable?
    // TODO(holtgrew): Does not compile with Dna instead of Dna5 here.

    // Tests with character existing in both alphabets -- true.
//     SEQAN_ASSERT_TRUE(TDnaWithGap('A') == Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') == TDnaWithGap('A'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('A') == Dna5('A'));
    SEQAN_ASSERT_TRUE(Dna5('A') == TDnaWithGap('A'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('A') == Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') == TDnaWithX('A'));
    SEQAN_ASSERT_TRUE(TDnaWithX('A') == Dna5('A'));
    SEQAN_ASSERT_TRUE(Dna5('A') == TDnaWithX('A'));
    // Tests with character existing in both alphabets -- false
//     SEQAN_ASSERT_NOT(TDnaWithGap('A') == Dna('C'));
//     SEQAN_ASSERT_NOT(Dna('C') == TDnaWithGap('A'));
    SEQAN_ASSERT_NOT(TDnaWithGap('A') == Dna5('C'));
    SEQAN_ASSERT_NOT(Dna5('C') == TDnaWithGap('A'));
//     SEQAN_ASSERT_NOT(TDnaWithX('A') == Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('A') == TDnaWithX('A'));
    SEQAN_ASSERT_NOT(TDnaWithX('A') == Dna5('C'));
    SEQAN_ASSERT_NOT(Dna5('C') == TDnaWithX('A'));

    // Tests with new character -- true.
//     SEQAN_ASSERT_TRUE(TDnaWithGap('-') == Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') == TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('-') == Dna5('N'));
    SEQAN_ASSERT_TRUE(Dna5('N') == TDnaWithGap('-'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('X') == Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') == TDnaWithX('X'));
    SEQAN_ASSERT_TRUE(TDnaWithX('X') == Dna5('N'));
    SEQAN_ASSERT_TRUE(Dna5('N') == TDnaWithX('X'));
    // Tests with new character -- false.
//     SEQAN_ASSERT_NOT(TDnaWithGap('-') == Dna('G'));
//     SEQAN_ASSERT_NOT(Dna('G') == TDnaWithGap('-'));
    SEQAN_ASSERT_NOT(TDnaWithGap('-') == Dna5('G'));
    SEQAN_ASSERT_NOT(Dna5('G') == TDnaWithGap('-'));
//     SEQAN_ASSERT_NOT(TDnaWithX('X') == Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('A') == TDnaWithX('X'));
    SEQAN_ASSERT_NOT(TDnaWithX('X') == Dna5('G'));
    SEQAN_ASSERT_NOT(Dna5('G') == TDnaWithX('X'));

    // Test with new character only -- true
    SEQAN_ASSERT_FAIL("The following does not compile.");
    /*
    SEQAN_ASSERT_TRUE(TDnaWithX('X') == TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('-') == TDnaWithX('X'));
    */
    // Test with new character only -- false
    /*
    SEQAN_ASSERT_NOT(TDnaWithXY('Y') == TDnaWithGap('-'));
    SEQAN_ASSERT_NOT(TDnaWithGap('-') == TDnaWithXY('Y'));
    */
}


// Test the operator!= implementations.
SEQAN_DEFINE_TEST(test_modifier_alphabet_operator_neq) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;
    typedef ModifiedAlphabet<TDnaWithX, ModExpand<'Y'> > TDnaWithXY;

    // TODO(holtgrew): Should that many alphabets be comparable?
    // TODO(holtgrew): Does not compile with Dna instead of Dna5 here.

    // Tests with character existing in both alphabets -- true.
//     SEQAN_ASSERT_NOT(TDnaWithGap('A') != Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('A') != TDnaWithGap('A'));
    SEQAN_ASSERT_NOT(TDnaWithGap('A') != Dna5('A'));
    SEQAN_ASSERT_NOT(Dna5('A') != TDnaWithGap('A'));
//     SEQAN_ASSERT_NOT(TDnaWithX('A') != Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('A') != TDnaWithX('A'));
    SEQAN_ASSERT_NOT(TDnaWithX('A') != Dna5('A'));
    SEQAN_ASSERT_NOT(Dna5('A') != TDnaWithX('A'));
    // Tests with character existing in both alphabets -- false
//     SEQAN_ASSERT_TRUE(TDnaWithGap('A') != Dna('C'));
//     SEQAN_ASSERT_TRUE(Dna('C') != TDnaWithGap('A'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('A') != Dna5('C'));
    SEQAN_ASSERT_TRUE(Dna5('C') != TDnaWithGap('A'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('A') != Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') != TDnaWithX('A'));
    SEQAN_ASSERT_TRUE(TDnaWithX('A') != Dna5('C'));
    SEQAN_ASSERT_TRUE(Dna5('C') != TDnaWithX('A'));

    // Tests with new character -- true.
//     SEQAN_ASSERT_NOT(TDnaWithGap('-') != Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('A') != TDnaWithGap('-'));
    SEQAN_ASSERT_NOT(TDnaWithGap('-') != Dna5('N'));
    SEQAN_ASSERT_NOT(Dna5('N') != TDnaWithGap('-'));
//     SEQAN_ASSERT_NOT(TDnaWithX('X') != Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('A') != TDnaWithX('X'));
    SEQAN_ASSERT_NOT(TDnaWithX('X') != Dna5('N'));
    SEQAN_ASSERT_NOT(Dna5('N') != TDnaWithX('X'));
    // Tests with new character -- false.
//     SEQAN_ASSERT_TRUE(TDnaWithGap('-') != Dna('G'));
//     SEQAN_ASSERT_TRUE(Dna('G') != TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('-') != Dna5('G'));
    SEQAN_ASSERT_TRUE(Dna5('G') != TDnaWithGap('-'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('X') != Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') != TDnaWithX('X'));
    SEQAN_ASSERT_TRUE(TDnaWithX('X') != Dna5('G'));
    SEQAN_ASSERT_TRUE(Dna5('G') != TDnaWithX('X'));

    // Test with new character only -- true
    SEQAN_ASSERT_FAIL("The following does not compile.");
    /*
    SEQAN_ASSERT_NOT(TDnaWithX('X') != TDnaWithGap('-'));
    SEQAN_ASSERT_NOT(TDnaWithGap('-') != TDnaWithX('X'));
    */
    // Test with new character only -- false
    /*
    SEQAN_ASSERT_TRUE(TDnaWithXY('Y') != TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('-') != TDnaWithXY('Y'));
    */
}


// Test the operator< implementations.
SEQAN_DEFINE_TEST(test_modifier_alphabet_operator_lt) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;
    typedef ModifiedAlphabet<TDnaWithX, ModExpand<'Y'> > TDnaWithXY;

    // TODO(holtgrew): Should that many alphabets be comparable?
    // TODO(holtgrew): Does not compile with Dna instead of Dna5 here.

    // Tests with character existing in both alphabets -- true.
//     SEQAN_ASSERT_TRUE(TDnaWithGap('A') < Dna('C'));
//     SEQAN_ASSERT_TRUE(Dna('A') < TDnaWithGap('C'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('A') < Dna5('C'));
    SEQAN_ASSERT_TRUE(Dna5('A') < TDnaWithGap('C'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('A') < Dna('C'));
//     SEQAN_ASSERT_TRUE(Dna('A') < TDnaWithX('C'));
    SEQAN_ASSERT_TRUE(TDnaWithX('A') < Dna5('C'));
    SEQAN_ASSERT_TRUE(Dna5('A') < TDnaWithX('C'));
    // Tests with character existing in both alphabets -- false
//     SEQAN_ASSERT_NOT(TDnaWithGap('C') < Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('C') < TDnaWithGap('A'));
    SEQAN_ASSERT_NOT(TDnaWithGap('C') < Dna5('A'));
    SEQAN_ASSERT_NOT(Dna5('C') < TDnaWithGap('A'));
//     SEQAN_ASSERT_NOT(TDnaWithX('C') < Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('C') < TDnaWithX('A'));
    SEQAN_ASSERT_NOT(TDnaWithX('C') < Dna5('A'));
    SEQAN_ASSERT_NOT(Dna5('C') < TDnaWithX('A'));

    // Tests with new character -- true.
//     SEQAN_ASSERT_TRUE(Dna('T') < TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(Dna5('T') < TDnaWithGap('-'));
//     SEQAN_ASSERT_TRUE(Dna('T') < TDnaWithX('X'));
    SEQAN_ASSERT_TRUE(Dna5('T') < TDnaWithX('X'));
    // Tests with new character -- false.
//     SEQAN_ASSERT_NOT(Dna('T') < TDnaWithGap('-'));
    SEQAN_ASSERT_NOT(Dna5('T') < TDnaWithGap('-'));
//     SEQAN_ASSERT_NOT(Dna('T') < TDnaWithX('X'));
    SEQAN_ASSERT_NOT(Dna5('T') < TDnaWithX('X'));

    // Test with new character only -- false.
    SEQAN_ASSERT_FAIL("The following does not compile.");
    /*
    SEQAN_ASSERT_TRUE(TDnaWithX('X') < TDnaWithXY('Y'));
    */
    // Test with new character only -- false
    /*
    SEQAN_ASSERT_NOT(TDnaWithXY('Y') < TDnaWithX('X'));
    */
}


// Test the operator> implementations.
SEQAN_DEFINE_TEST(test_modifier_alphabet_operator_gt) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;
    typedef ModifiedAlphabet<TDnaWithX, ModExpand<'Y'> > TDnaWithXY;

    // TODO(holtgrew): Should that many alphabets be comparable?
    // TODO(holtgrew): Does not compile with Dna instead of Dna5 here.

    // Tests with character existing in both alphabets -- true.
//     SEQAN_ASSERT_TRUE(TDnaWithGap('C') > Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('C') > TDnaWithGap('A'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('C') > Dna5('A'));
    SEQAN_ASSERT_TRUE(Dna5('C') > TDnaWithGap('A'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('C') > Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('C') > TDnaWithX('A'));
    SEQAN_ASSERT_TRUE(TDnaWithX('C') > Dna5('A'));
    SEQAN_ASSERT_TRUE(Dna5('C') > TDnaWithX('A'));
    // Tests with character existing in both alphabets -- false.
//     SEQAN_ASSERT_NOT(TDnaWithGap('A') > Dna('C'));
//     SEQAN_ASSERT_NOT(Dna('A') > TDnaWithGap('C'));
    SEQAN_ASSERT_NOT(TDnaWithGap('A') > Dna5('C'));
    SEQAN_ASSERT_NOT(Dna5('A') > TDnaWithGap('C'));
//     SEQAN_ASSERT_NOT(TDnaWithX('A') > Dna('C'));
//     SEQAN_ASSERT_NOT(Dna('A') > TDnaWithX('C'));
    SEQAN_ASSERT_NOT(TDnaWithX('A') > Dna5('C'));
    SEQAN_ASSERT_NOT(Dna5('A') > TDnaWithX('C'));

    // Test with new character only -- true.
    SEQAN_ASSERT_FAIL("The following does not compile.");
    /*
    SEQAN_ASSERT_NOT(TDnaWithX('X') > TDnaWithXY('Y'));
    */
    // Test with new character only -- false
    /*
    SEQAN_ASSERT_TRUE(TDnaWithXY('Y') > TDnaWithX('X'));
    */

    // Tests with new character only -- false.
//     SEQAN_ASSERT_NOT(Dna('T') > TDnaWithGap('-'));
    SEQAN_ASSERT_NOT(Dna5('T') > TDnaWithGap('-'));
//     SEQAN_ASSERT_NOT(Dna('T') > TDnaWithX('X'));
    SEQAN_ASSERT_NOT(Dna5('T') > TDnaWithX('X'));
    // Tests with new character -- false.
//     SEQAN_ASSERT_TRUE(Dna('T') > TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(Dna5('T') > TDnaWithGap('-'));
//     SEQAN_ASSERT_TRUE(Dna('T') > TDnaWithX('X'));
    SEQAN_ASSERT_TRUE(Dna5('T') > TDnaWithX('X'));
}


// Test the operator<= implementations.
SEQAN_DEFINE_TEST(test_modifier_alphabet_operator_leq) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;
    typedef ModifiedAlphabet<TDnaWithX, ModExpand<'Y'> > TDnaWithXY;

    // TODO(holtgrew): Should that many alphabets be comparable?
    // TODO(holtgrew): Does not compile with Dna instead of Dna5 here.

    // Tests with character existing in both alphabets -- true.
//     SEQAN_ASSERT_TRUE(TDnaWithGap('A') <= Dna('A'));
//     SEQAN_ASSERT_TRUE(TDnaWithGap('A') <= Dna('C'));
//     SEQAN_ASSERT_TRUE(Dna('A') <= TDnaWithGap('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') <= TDnaWithGap('C'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('A') <= Dna5('A'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('A') <= Dna5('C'));
    SEQAN_ASSERT_TRUE(Dna5('A') <= TDnaWithGap('A'));
    SEQAN_ASSERT_TRUE(Dna5('A') <= TDnaWithGap('C'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('A') <= Dna('A'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('A') <= Dna('C'));
//     SEQAN_ASSERT_TRUE(Dna('A') <= TDnaWithX('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') <= TDnaWithX('C'));
    SEQAN_ASSERT_TRUE(TDnaWithX('A') <= Dna5('A'));
    SEQAN_ASSERT_TRUE(TDnaWithX('A') < Dna5('C'));
    SEQAN_ASSERT_TRUE(Dna5('A') <= TDnaWithX('A'));
    SEQAN_ASSERT_TRUE(Dna5('A') < TDnaWithX('C'));
    // Tests with character existing in both alphabets -- false
//     SEQAN_ASSERT_NOT(TDnaWithGap('C') <= Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('C') <= TDnaWithGap('A'));
    SEQAN_ASSERT_NOT(TDnaWithGap('C') <= Dna5('A'));
    SEQAN_ASSERT_NOT(Dna5('C') <= TDnaWithGap('A'));
//     SEQAN_ASSERT_NOT(TDnaWithX('C') <= Dna('A'));
//     SEQAN_ASSERT_NOT(Dna('C') <= TDnaWithX('A'));
    SEQAN_ASSERT_NOT(TDnaWithX('C') <= Dna5('A'));
    SEQAN_ASSERT_NOT(Dna5('C') <= TDnaWithX('A'));

    // Tests with new character -- true.
    SEQAN_ASSERT_TRUE(TDnaWithGap('-') <= TDnaWithGap('-'));
//     SEQAN_ASSERT_TRUE(Dna('T') <= TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(Dna5('N') <= TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(Dna5('T') <= TDnaWithGap('-'));
//     SEQAN_ASSERT_TRUE(Dna('T') <= TDnaWithX('X'));
//     SEQAN_ASSERT_TRUE(Dna('T') <= TDnaWithX('X'));
    SEQAN_ASSERT_TRUE(Dna5('N') <= TDnaWithX('X'));
    SEQAN_ASSERT_TRUE(Dna5('T') <= TDnaWithX('X'));
    // Tests with new character -- false.
//     SEQAN_ASSERT_NOT(Dna('T') <= TDnaWithGap('-'));
    SEQAN_ASSERT_NOT(Dna5('T') <= TDnaWithGap('-'));
//     SEQAN_ASSERT_NOT(Dna('T') <= TDnaWithX('X'));
    SEQAN_ASSERT_NOT(Dna5('T') <= TDnaWithX('X'));

    // Test with new character only -- false.
    SEQAN_ASSERT_FAIL("The following does not compile.");
    /*
    SEQAN_ASSERT_TRUE(TDnaWithX('X') <= TDnaWithXY('Y'));
    */
    // Test with new character only -- false
    /*
    SEQAN_ASSERT_NOT(TDnaWithXY('Y') <= TDnaWithX('X'));
    */
}

// Test the operator>= implementations.
SEQAN_DEFINE_TEST(test_modifier_alphabet_operator_geq) {
    // Add special gap and an arbitrary symbol to Dna.
    typedef ModifiedAlphabet<Dna, ModExpand<'-'> > TDnaWithGap;
    typedef ModifiedAlphabet<Dna, ModExpand<'X'> > TDnaWithX;
    typedef ModifiedAlphabet<TDnaWithX, ModExpand<'Y'> > TDnaWithXY;

    // TODO(holtgrew): Should that many alphabets be comparable?
    // TODO(holtgrew): Does not compile with Dna instead of Dna5 here.

    // Tests with character existing in both alphabets -- true.
//     SEQAN_ASSERT_TRUE(TDnaWithGap('A') >= Dna('A'));
//     SEQAN_ASSERT_TRUE(TDnaWithGap('C') >= Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') >= TDnaWithGap('A'));
//     SEQAN_ASSERT_TRUE(Dna('C') >= TDnaWithGap('A'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('A') >= Dna5('A'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('C') >= Dna5('A'));
    SEQAN_ASSERT_TRUE(Dna5('A') >= TDnaWithGap('A'));
    SEQAN_ASSERT_TRUE(Dna5('C') >= TDnaWithGap('A'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('A') >= Dna('A'));
//     SEQAN_ASSERT_TRUE(TDnaWithX('C') >= Dna('A'));
//     SEQAN_ASSERT_TRUE(Dna('A') >= TDnaWithX('A'));
//     SEQAN_ASSERT_TRUE(Dna('C') >= TDnaWithX('A'));
    SEQAN_ASSERT_TRUE(TDnaWithX('A') >= Dna5('A'));
    SEQAN_ASSERT_TRUE(TDnaWithX('C') >= Dna5('A'));
    SEQAN_ASSERT_TRUE(Dna5('A') >= TDnaWithX('A'));
    SEQAN_ASSERT_TRUE(Dna5('C') > TDnaWithX('A'));
    // Tests with character existing in both alphabets -- false.
//     SEQAN_ASSERT_NOT(TDnaWithGap('A') >= Dna('C'));
//     SEQAN_ASSERT_NOT(Dna('A') >= TDnaWithGap('C'));
    SEQAN_ASSERT_NOT(TDnaWithGap('A') >= Dna5('C'));
    SEQAN_ASSERT_NOT(Dna5('A') >= TDnaWithGap('C'));
//     SEQAN_ASSERT_NOT(TDnaWithX('A') >= Dna('C'));
//     SEQAN_ASSERT_NOT(Dna('A') >= TDnaWithX('C'));
    SEQAN_ASSERT_NOT(TDnaWithX('A') >= Dna5('C'));
    SEQAN_ASSERT_NOT(Dna5('A') >= TDnaWithX('C'));

    // Test with new character only -- true.
    SEQAN_ASSERT_FAIL("The following does not compile.");
    /*
    SEQAN_ASSERT_NOT(TDnaWithX('X') >= TDnaWithXY('Y'));
    */
    // Test with new character only -- false
    /*
    SEQAN_ASSERT_TRUE(TDnaWithXY('Y') >= TDnaWithX('X'));
    */

    // Tests with new character -- true.
//     SEQAN_ASSERT_TRUE(Dna('T') >= TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(TDnaWithGap('-') >= TDnaWithGap('-'));
    SEQAN_ASSERT_TRUE(Dna5('T') > TDnaWithGap('-'));
//     SEQAN_ASSERT_TRUE(Dna('T') >= TDnaWithX('X'));
    SEQAN_ASSERT_TRUE(TDnaWithX('X') > TDnaWithX('X'));
    SEQAN_ASSERT_TRUE(Dna5('T') > TDnaWithX('X'));
    // Tests with new character only -- false.
//     SEQAN_ASSERT_NOT(Dna('T') > TDnaWithGap('-'));
    SEQAN_ASSERT_NOT(Dna5('T') > TDnaWithGap('-'));
//     SEQAN_ASSERT_NOT(Dna('T') > TDnaWithX('X'));
    SEQAN_ASSERT_NOT(Dna5('T') > TDnaWithX('X'));
}

#endif  // TESTS_MODIFIER_TEST_MODIFIER_ALPHABET_H_
