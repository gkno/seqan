// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main ()
{
	StringSet<CharString> myStringSet;
	appendValue(myStringSet, "CDFGHC");
	appendValue(myStringSet, "CDEFGAHC");

	typedef Index< StringSet<CharString> > TMyIndex;
	TMyIndex myIndex(myStringSet);

// FRAGMENT(iteration)
	Iterator< TMyIndex, MUMs >::Type myIterator(myIndex);

	while (!atEnd(myIterator))
	{
		std::cout << representative(myIterator) << std::endl;
		++myIterator;
	}

	return 0;
}
