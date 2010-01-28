// part 1
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/basic.h>

using namespace seqan;
// part 2
template<typename TAlphabet>
void
showAllLetterOfMyAlphabet(TAlphabet const&) {
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	for(TSize i = 0; i < alphSize; ++i) {
		::std::cout << TAlphabet(i) << ',';
	}
	::std::cout << ::std::endl;
}
//part 3
int main()
{
	showAllLetterOfMyAlphabet(AminoAcid());
	showAllLetterOfMyAlphabet(Dna());
	showAllLetterOfMyAlphabet(Dna5());
	return 0;
}
