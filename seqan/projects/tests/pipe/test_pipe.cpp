#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/pipe.h>
#include "test_pipe.h"

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////


int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	testContainers();

	debug::verifyCheckpoints("projects/library/seqan/pipe/pool_base.h");

	SEQAN_TREPORT("TEST END")

	return 0;
}
