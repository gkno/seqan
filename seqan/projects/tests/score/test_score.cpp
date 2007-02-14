#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <time.h>
#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/score.h>

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")
//____________________________________________________________________________

	SEQAN_TREPORT("TEST END")

	return 0;
}
