#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <ctime>
#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>

#include <seqan/map.h>

#include <seqan/misc/edit_environment.h>
#include <seqan/misc/misc_base.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/misc/misc_dequeue.h>
#include <seqan/misc/misc_long_word.h>
#include <seqan/misc/misc_map.h>
#include <seqan/misc/misc_random.h>
#include <seqan/misc/misc_set.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

#include "test_misc_long_word.h"

using namespace std;
using namespace seqan;


SEQAN_DEFINE_TEST(test_misc_random) {
	mtRandInit();

	for (unsigned int i=0; i<100; ++i)
	{
		cout << mtRand() << ", ";
	}
	cout << "\n\n";

	for (unsigned int i=0; i<100; ++i)
	{
		cout << geomRand<int>() << ", ";
	}

}


SEQAN_BEGIN_TESTSUITE(test_misc) {
//     SEQAN_CALL_TEST(test_misc_random);

    SEQAN_CALL_TEST(test_misc_long_word_native_interface);
    SEQAN_CALL_TEST(test_misc_long_word_static_interface);

//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/edit_environment.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_base.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_cmdparser.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_dequeue.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_map.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_long_word.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_parsing.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_random.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/misc_set.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/priority_type_base.h");
//     SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/misc/priority_type_heap.h");
}
SEQAN_END_TESTSUITE

