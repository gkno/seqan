// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main ()
{
	StringSet<CharString> myStringSet;
	appendValue(myStringSet, "tobeornottobe");
	appendValue(myStringSet, "thebeeonthecomb");
	appendValue(myStringSet, "beingjohnmalkovich");

	typedef Index< StringSet<CharString> > TMyIndex;
	TMyIndex myIndex(myStringSet);

// FRAGMENT(iteration1)
	Iterator< TMyIndex, TopDown< ParentLinks<Postorder> > >::Type myIterator(myIndex);

	// top-down iterators start in the root node which is not the first node
	// of a postorder DFS thus we have to manually go the DFS start with goBegin
	goBegin(myIterator);
	while (!atEnd(myIterator))
	{
		std::cout << representative(myIterator) << std::endl;
		++myIterator;
	}

// FRAGMENT(iteration2)
	Iterator< TMyIndex, BottomUp<> >::Type myIterator2(myIndex);

	while (!atEnd(myIterator2))
	{
		std::cout << representative(myIterator2) << std::endl;
		++myIterator2;
	}

	return 0;
}
