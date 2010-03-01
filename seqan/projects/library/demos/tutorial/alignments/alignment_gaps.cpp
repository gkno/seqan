// FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
// FRAGMENT(typedefs)
	typedef String<char> TSequence;
	typedef Gaps<String<char>,ArrayGaps> TGaps;


// FRAGMENT(init)
	TSequence seq1 = "CDFGHC";
	TSequence seq2 = "CDEFGAHC";

	TGaps gaps1;
	TGaps gaps2;
	setSource(gaps1,seq1);
	setSource(gaps2,seq2);

// FRAGMENT(manipulation)
	insertGap(gaps1,2);
	insertGap(gaps1,5);

	::std::cout << (TSequence)gaps1 << ::std::endl;
	::std::cout << (TSequence)gaps2 << ::std::endl;

// FRAGMENT(printingViewPos)
	::std::cout << ::std::endl << "ViewToSource1: ";
	for(unsigned i = 0; i < length(gaps1); ++i)
		::std::cout << toSourcePosition(gaps1, i) << ",";
	
	::std::cout << ::std::endl << "ViewToSource2: ";
	for(unsigned i = 0; i < length(gaps2); ++i)
		::std::cout << toSourcePosition(gaps2, i) << ",";
	::std::cout << ::std::endl;

// FRAGMENT(printingSourcePos)
	::std::cout << ::std::endl << "SourceToView1: ";
	for(unsigned i = 0; i < length(source(gaps1)); ++i)
		::std::cout << toViewPosition(gaps1, i) << ",";

	::std::cout << ::std::endl << "SourceToView2: ";
	for(unsigned i = 0; i < length(source(gaps2)); ++i)
		::std::cout << toViewPosition(gaps2, i) << ",";
	::std::cout << ::std::endl;



	return 0;
}
