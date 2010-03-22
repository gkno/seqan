// FRAGMENT(initialization)
#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main()
{
	CharString haystack = "send more money!";
	CharString needle = "mo";

// FRAGMENT(output)
	Finder<CharString> finder(haystack);
	Pattern<CharString, Horspool> pattern(needle);
	while (find(finder, pattern))
		std::cout << position(finder) << std::endl;

	return 0;
}

