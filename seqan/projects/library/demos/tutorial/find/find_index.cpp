// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(index_esa)
int main()
{
	Index<CharString> index_esa("tobeornottobe");
	Finder< Index<CharString> > finder_esa(index_esa);

	while (find(finder_esa, "be"))
		std::cout << position(finder_esa) << std::endl;

// FRAGMENT(index_qgram)
	typedef Index<CharString, Index_QGram< UngappedShape<2> > > TQGramIndex;
	TQGramIndex index_2gram("tobeornottobe");
	Finder< TQGramIndex > finder_2gram(index_2gram);

	while (find(finder_2gram, "be"))
		std::cout << position(finder_2gram) << std::endl;

	return 0;
}
