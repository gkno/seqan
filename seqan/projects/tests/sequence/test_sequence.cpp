#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

#include "test_string.h"
#include "test_stringset.h"
#include "test_segment.h"

SEQAN_BEGIN_TESTSUITE(Sequence tests)
{
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

//	debug::verifyCheckpoints("projects/library/seqan/sequence.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/sequence_interface.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/string_base.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/string_alloc.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/string_pointer.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/string_array.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/string_cstyle.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/lexical.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/std_string.h");

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

//	debug::verifyCheckpoints("projects/library/seqan/sequence/segment_infix.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/segment_suffix.h");
//	debug::verifyCheckpoints("projects/library/seqan/sequence/segment_prefix.h");
	//debug::verifyCheckpoints("projects/library/seqan/sequence/segment_gram.h");

#define ROOT "projects/tests/sequence/"
    SEQAN_VERIFY_CHECKPOINTS(ROOT"test_string.h");
    SEQAN_VERIFY_CHECKPOINTS(ROOT"test_stringset.h");
    SEQAN_VERIFY_CHECKPOINTS(ROOT"test_segment.h");
#undef ROOT
} SEQAN_END_TESTSUITE
