// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main ()
{
	typedef Index<CharString> TIndex;
	TIndex index("How many wood would a woodchuck chuck.");
	Iterator< TIndex, TopDown<> >::Type it(index);
// FRAGMENT(iteration)

	CharString pattern = "wooud";
	while (repLength(it) < length(pattern))
	{
		// go down edge starting with the next pattern character
		if (!goDown(it, pattern[repLength(it)])) break;
		unsigned endPos = _min(repLength(it), length(pattern));
		// compare remaining edge characters with pattern
		std::cout << representative(it) << std::endl;
		if (infix(representative(it), parentRepLength(it) + 1, endPos) !=
			infix(pattern, parentRepLength(it) + 1, endPos)) break;
	}

// FRAGMENT(output)
	// output match positions if pattern was found
	if (repLength(it) >= length(pattern))
		for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
			std::cout << getOccurrences(it)[i] << std::endl;

	return 0;
}
