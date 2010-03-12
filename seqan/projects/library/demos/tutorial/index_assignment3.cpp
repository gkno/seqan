// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main ()
{
	typedef Index<CharString> TIndex;
	TIndex index("How many wood would a woodchuck chuck.");
	Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
// FRAGMENT(iteration)

	do {
		std::cout << '"' << representative(it) << '"' << std::endl;
		if (!goDown(it) && !goRight(it))
			while (goUp(it) && !goRight(it)) ;
	} while (!isRoot(it));

	return 0;
}
