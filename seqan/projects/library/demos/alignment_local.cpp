#include <iostream>
#include <seqan/align.h>
#include <seqan/graph.h>

using namespace std;
using namespace seqan;

int main()
{
///This program computes the best local alignment between two given sequences.
	Align< String<char> > ali;
	appendValue(rows(ali), "aphilologicaltheorem");
	appendValue(rows(ali), "bizarreamphibology");
    cout << "Score = " << localAlignment(ali, Score<int>(3,-3,-2), SmithWaterman()) << endl;
	cout << ali;

	return 0;
}
