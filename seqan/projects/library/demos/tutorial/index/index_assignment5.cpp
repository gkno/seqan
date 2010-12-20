// FRAGMENT(includes)
#include <iostream>
#include <seqan/index.h>

using namespace seqan;

// FRAGMENT(initialization)
int main ()
{
	Index<DnaString, Index_QGram<OneGappedShape> > index("CATGATTACATA");
//	Index<DnaString, Index_QGram<GenericShape> > index("CATGATTACATA");
	stringToShape(indexShape(index), "1101");
//	Index<DnaString, Index_QGram<GappedShape<HardwiredShape<1,2> > > > index("CATGATTACATA");

// FRAGMENT(output)
	hash(indexShape(index), "ATCA");
	for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
		std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;

	return 0;
}
