// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
template < typename TIndexSpec >
void constrainedDFS ()
{
	typedef Index<CharString, TIndexSpec> TIndex;
	TIndex index("How many wood would a woodchuck chuck.");
	typename Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
// FRAGMENT(iteration)

	do {
		std::cout << '"' << representative(it) << '"' << std::endl;
		do {
			if (!(goDown(it) && repLength(it) <= 4) && !goRight(it))
				while (goUp(it) && !goRight(it)) ;
		} while (repLength(it) > 4); 
	} while (!isRoot(it));
}


int main ()
{
	constrainedDFS< Index_ESA<> > ();
	constrainedDFS< Index_Wotd<> > ();
	return 0;
}
