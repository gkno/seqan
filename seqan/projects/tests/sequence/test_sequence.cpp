#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

#include "test_string.h"
#include "test_stringset.h"
#include "test_segment.h"
#include "test_sequence_std_adaptions.h"

SEQAN_BEGIN_TESTSUITE(Sequence tests)
{
    // -----------------------------------------------------------------------
    // Tests for STL adaptions.
    // -----------------------------------------------------------------------
    //
    // Test adaptions for std::string.
    SEQAN_CALL_TEST(test_sequence_adaptions_metafunctions_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_iterators_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_interface_std_string);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_memory_std_string);

    // Test adaptions for std::list.
    SEQAN_CALL_TEST(test_sequence_adaptions_metafunctions_std_list);
    SEQAN_CALL_TEST(test_sequence_adaptions_iterators_std_list);
    SEQAN_CALL_TEST(test_sequence_adaptions_sequence_interface_std_list);
    
    // Use the constant EMPTY_STRING once to get rid of linker error in Visual Studio.
    (void)seqan::String<char const ,struct seqan::CStyle>::EMPTY_STRING;

	SEQAN_CALL_TEST(Sequence_Interface);
	SEQAN_CALL_TEST(String_Base);
	SEQAN_CALL_TEST(String_Alloc);
	SEQAN_CALL_TEST(String_Array);
	SEQAN_CALL_TEST(String_Stack);
	SEQAN_CALL_TEST(String_Pointer);
	SEQAN_CALL_TEST(String_CStyle);
	SEQAN_CALL_TEST(String_Packed);
	SEQAN_CALL_TEST(Std_String);

	SEQAN_CALL_TEST(Lexical);
	SEQAN_CALL_TEST(Combinatoric);

	SEQAN_CALL_TEST(Segment);

    SEQAN_CALL_TEST(StringSet_Owner_Default);
    SEQAN_CALL_TEST(StringSet_Concat_Owner_Default);
    SEQAN_CALL_TEST(StringSet_Concat_Owner_ConcatDirect);
    SEQAN_CALL_TEST(StringSet_Id_Dependent_Tight);
    SEQAN_CALL_TEST(StringSet_Id_Dependent_Generous);
    SEQAN_CALL_TEST(StringSetIdHolder_Char_Dependent_Tight);
    SEQAN_CALL_TEST(StringSetIdHolder_Char_Dependent_Generous);

//	debug::verifyCheckpoints("projects/library/seqan/sequence/sequence_multiple.h");

    SEQAN_CALL_TEST(Infix);
    SEQAN_CALL_TEST(Suffix);
    SEQAN_CALL_TEST(Ticket317);

    // -----------------------------------------------------------------------
    // Checkpoint Verification
    // -----------------------------------------------------------------------
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence/std_string.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/sequence/adapt_std_list.h");
    // TODO(holtgrew): Add more checkpoints.
}
SEQAN_END_TESTSUITE
