#define SEQAN_DEBUG
#define SEQAN_TEST
//#define SEQAN_VVERBOSE

#define TEST_PATH "projects/tests/blast/"
#define LIB_PATH "projects/library/seqan/blast/"


#include <seqan/file.h>
#include <seqan/blast.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>
#include <string>
#include <ctime>


//#include "test_blast_calling.h"
#include "test_blast_parsing.h"

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	//Test_BlastCalling();
	//Test_BlastHsp<int>(); 
	//Test_BlastSclGCoffee<int>();
	Test_BlastStoreReport<int>();
	//Test_BlastParsing<int>();
	//debug::verifyCheckpoints("projects/library/seqan/blast/blast_parsing.h");
	//debug::verifyCheckpoints("projects/tests/blast/test_blast_parsing.h");


	SEQAN_TREPORT("TEST END")

	return 0;
}
