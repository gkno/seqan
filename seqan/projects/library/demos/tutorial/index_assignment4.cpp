// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>
#include <typeinfo>

using namespace seqan;

// FRAGMENT(initialization)
template < typename TIndexSpec >
void constrainedDFS ()
{
	typedef Index<CharString, TIndexSpec> TIndex;
	TIndex index("tobeornottobe");
	typename Iterator< TIndex, TopDown<ParentLinks<> > >::Type it(index);
// FRAGMENT(iteration)

	do {
		std::cout << representative(it) << std::endl;
		do {
			if (!(goDown(it) && repLength(it) <= 3) && !goRight(it))
				while (goUp(it) && !goRight(it)) ;
		} while (repLength(it) > 3); 
	} while (!isRoot(it));
	std::cout << std::endl;
}


int main ()
{
	constrainedDFS< Index_ESA<> > ();
	constrainedDFS< Index_Wotd<> > ();
	return 0;
}
