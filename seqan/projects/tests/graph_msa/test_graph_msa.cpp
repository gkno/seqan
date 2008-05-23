#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

#define TEST_PATH "projects/tests/graph/"
#define LIB_PATH "projects/library/seqan/graph/"

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
#include <seqan/graph_msa.h>
#include "test_graph_tcoffee.h"

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_GraphTCoffee();		// Test T-Coffee

	SEQAN_TREPORT("TEST END")

	return 0;
}
