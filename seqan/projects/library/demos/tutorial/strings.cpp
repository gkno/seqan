#include <iostream>
#include <seqan/basic.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
	typedef String<AminoAcid> TAminoAcidString;
	TAminoAcidString sourceSeq = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
// part 2
	typedef Iterator<TAminoAcidString>::Type TIter;

	TIter it = begin(sourceSeq);
	TIter itEnd = end(sourceSeq);
	for(;it != itEnd; goNext(it)) {
		if (value(it) == 'R') value(it) = 'A';
		::std::cout << value(it) << ',';
	}
	::std::cout << ::std::endl;
// part 3
	typedef Size<TAminoAcidString>::Type TSize;
	typedef String<TSize> TCounterString;
	TCounterString counter;
	TSize alphSize = ValueSize<AminoAcid>::VALUE;
	fill(counter, alphSize, 0);
	it = begin(sourceSeq);
	for(;it != itEnd; goNext(it)) {
		value(counter, ordValue(value(it))) += 1;
	}
// part 4
	typedef Iterator<TCounterString>::Type TCounterIter;
	TCounterIter countIt = begin(counter);
	TCounterIter countItEnd = end(counter);
	TSize pos = 0;
	for(;countIt != countItEnd; ++countIt, ++pos) {
		::std::cout << AminoAcid(pos) << ':' << value(countIt) << ::std::endl;
	}


	return 0;
}

