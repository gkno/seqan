#ifndef SEQAN_HEADER_TEST_GRAPH_CONSENSUS_H
#define SEQAN_HEADER_TEST_GRAPH_CONSENSUS_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{



//////////////////////////////////////////////////////////////////////////////


void Test_GraphConsensus() {
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_consensus_base.h");
	debug::verifyCheckpoints("projects/library/seqan/graph/graph_consensus_library.h");
}


}

#endif

