#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

#include <seqan/consensus.h>
#include "test_consensus.h"


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_Consensus();		// Test Graph Consensus

	SEQAN_TREPORT("TEST END")

	return 0;
}
