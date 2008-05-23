#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

//STL
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>
#include <string>
#include <ctime>

//SeqAn
#include <seqan/graph_utils.h>
#include "test_graph_utils.h"

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_GraphUtils();		// Test Utilities

	SEQAN_TREPORT("TEST END")

	return 0;
}
