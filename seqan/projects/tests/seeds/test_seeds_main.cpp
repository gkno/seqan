#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/seeds.h>

void Main_BandedAlign();
void Main_GlobalSeedChain();
void Main_MemoryManager();
void Main_Seeds();
void Main_SeedSet();

SEQAN_DEFINE_TEST(test_seed_banded_align) {
    Main_BandedAlign();
}


SEQAN_DEFINE_TEST(test_seed_global_seed_chain) {
    Main_GlobalSeedChain();
}


SEQAN_DEFINE_TEST(test_seed_memory_manager) {
    Main_MemoryManager();
}


SEQAN_DEFINE_TEST(test_seed_seeds) {
    Main_Seeds();
}


SEQAN_DEFINE_TEST(test_seed_seed_set) {
    Main_SeedSet();
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_BEGIN_TESTSUITE(test_seed) {
    SEQAN_CALL_TEST(test_seed_banded_align);
    SEQAN_CALL_TEST(test_seed_global_seed_chain);
    SEQAN_CALL_TEST(test_seed_memory_manager);
    SEQAN_CALL_TEST(test_seed_seeds);
    SEQAN_CALL_TEST(test_seed_seed_set);
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/banded_align.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/banded_chain_align.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/banded_chain_align_affine.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/seed_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/seed_multi.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/global_seed_chain.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/memoryManager_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/memoryManager_int.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/seedSet_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/seeds/seedSet_score.h");
}
SEQAN_END_TESTSUITE