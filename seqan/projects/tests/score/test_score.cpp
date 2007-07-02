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
//compare to http://www.bioinformatics.nl/tools/pam.html

void testScorePAM()
{
	PamDayhoff pam;
	//Score<int, Pam<> > pam(250, -1, 0);
	//Score<int, Pam<AminoAcid, Pam_Data_Dayhoff_MDM78> > pam;

	SEQAN_TASSERT(getDist(pam) == 250);
	SEQAN_TASSERT(scoreGapExtend(pam) == -1);
	SEQAN_TASSERT(scoreGapOpen(pam) == 0);

	//print out the matrix
	for (unsigned int i = 0; i < 24; ++i)
	{
		AminoAcid a = i;
		cout << a << " ";
		for (unsigned int j = 0; j < 24; ++j)
		{
			AminoAcid b = j;
			printf("%2i ", score(pam, a, b));
		}
		cout << "\n";
	}
}


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	testScorePAM();

//____________________________________________________________________________

	SEQAN_TREPORT("TEST END")

	return 0;
}
