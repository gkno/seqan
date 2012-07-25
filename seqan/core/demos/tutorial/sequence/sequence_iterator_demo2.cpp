// The comment lines containing FRAGMENT(fragment-line) are there for the
// documentation system.  You can ignore them when reading this file.i
// This is the draft for the new iterator tutorial
// FRAGMENT(includes)
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main() {
// FRAGMENT(metafunctions)
	Dna5String genome = "TATANNNGCGCG";
	Iterator<Dna5String >::Type it = begin(genome);
	Iterator<Dna5String >::Type itEnd = end(genome);
// FRAGMENT(iterators)
	while (it != itEnd) {
		std::cout << *it;
		++it;
	}
	std::cout << std::endl;
// FRAGMENT(rooted-iterators)
	Iterator<Dna5String, Rooted >::Type it2 = begin(genome);
	for (goBegin(it2); !atEnd(it2); goNext(it2)) {
		if (getValue(it2) == 'N')
		    value(it2) = 'A';
	}
// FRAGMENT(iterator-reverse)
	goEnd(it2);
        while (!atBegin(it2)) {
		goPrevious(it2);
		std::cout << getValue(it2);
	}
	std::cout << std::endl;
// FRAGMENT(assign-value)
	assignValue(begin(genome), 'N');
	std::cout << genome << std::endl;
	
	return 0;
}
