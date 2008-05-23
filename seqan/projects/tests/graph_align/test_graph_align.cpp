#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

#include <seqan/graph_align.h>
#include "test_graph_align.h"

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_GraphAlignment();		// Test Graph Alignment

	SEQAN_TREPORT("TEST END")

	return 0;
}
