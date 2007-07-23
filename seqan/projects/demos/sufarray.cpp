#include <iostream>
#include <seqan/sequence.h>
#include <seqan/index.h>
using namespace seqan;

int main ()
{
	String<char> text = "hello world!";
	String<char> pattern = "l";
	String<unsigned> sa;

	resize(sa, length(text));
	createSuffixArray(sa, text, Skew7());

	Pair< Iterator<String<unsigned> const>::Type > hits;
	hits = equalRangeSA(text, sa, pattern);

	for(; hits.i1 != hits.i2; ++hits.i1)
		std::cout << *(hits.i1) << " ";
	std::cout << ::std::endl;
 
	return 0;
}
