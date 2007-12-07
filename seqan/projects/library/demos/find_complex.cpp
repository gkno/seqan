#include <iostream>
#include <seqan/find.h>

using namespace seqan;
using namespace std;

///This program uses the algorithm @Spec.WildShiftAnd@ to perform a wildcard search.
int main() 
{
	String<char> haystack = "send more money!";
	String<char> needle = "mo";

	printAllOccs<Horspool>(haystack, needle);
	printAllOccs<BomAlgo> (haystack, needle);
	printAllOccs<BndmAlgo>(haystack, needle);
	printAllOccs<ShiftAnd>(haystack, needle);
	printAllOccs<ShiftOr> (haystack, needle);

	return 0;
}

