#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

//#define SEQAN_DEBUG
//#define SEQAN_DEBUG_INDEX

//#define SEQAN_TEST
//#define SEQAN_TEST_SKEW3
//#define SEQAN_TEST_SKEW7

#include <seqan/index.h>
#include "test_index_creation.h"


using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////


bool testIndexCreation()
{
		typedef String<char> TText;
		typedef String<unsigned> TArray;

		TText	text;
		TArray	sa;
		TArray	lcp;
		TArray	child, childExt;
		TText	bwt;

		const int runs = 10;					// conduct 10 test runs 
		const int maxSize = 20 * 1024 * 1024;	// max text size is 20 megabyte

		for(int i = 0; i < runs; ++i) {

			cout << "*** RUN " << i << " ***" << endl;
			
			int size = rand() % maxSize;

//___randomize_text___________________________________________________________

			resize(text,size);
			if (i < runs/2)
				randomize(text);
			else
				textRandomize(text);
/*			String<char,External<> > errorText;	// read in text causing an error
			open(errorText,"error.txt");
			text = errorText;
*/
			size = length(text);
			cout << "text created (n=" << size << ")" << endl;

//___create_suffix_array______________________________________________________

			resize(sa, size);
			createSuffixArray(sa, text, Skew3());
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (internal Skew3) failed" << endl;
				return false;
			}
			cout << ".";

			blank(sa);
			createSuffixArray(sa, text, Skew7());
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (internal Skew7) failed" << endl;
				return false;
			}
			cout << ".";

			blank(sa);
			createSuffixArrayExt(sa, text, Skew3());
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (external Skew3) failed" << endl;
				return false;
			}
			cout << ".";

			blank(sa);
			createSuffixArrayExt(sa, text, Skew7());
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (external Skew7) failed" << endl;
				return false;
			}
			cout << ".";

//___create_lcp_table_________________________________________________________

			resize(lcp, size);
			createLCPTable(lcp, text, sa, Kasai());
			if (!isLCPTable(lcp, sa, text)) {
				cout << "suffix array creation (internal Kasai) failed" << endl;
				return false;
			}
			cout << ".";

			blank(lcp);
			createLCPTable(lcp, text, sa, KasaiInPlace());
			if (!isLCPTable(lcp, sa, text)) {
				cout << "suffix array creation (internal in-place Kasai) failed" << endl;
				return false;
			}
			cout << ".";

			blank(lcp);
			createLCPTableExt(lcp, text, sa, Kasai());
			if (!isLCPTable(lcp, sa, text)) {
				cout << "suffix array creation (external Kasai) failed" << endl;
				return false;
			}
			cout << ".";

//___create_child_table_______________________________________________________

			resize(child, size);
			for(int i=0; i<size; ++i)
				child[i] = supremumValue<unsigned>();
			createChildTable(child, lcp);
			cout << ".";

			unsigned undefs=0;
			for(int i=0; i<size; ++i)
				if (child[i] == supremumValue<unsigned>()) ++undefs;
			if (undefs) ::std::cout << undefs << " undefined values";

			resize(childExt, size);
			createChildTableExt(childExt, lcp);
			cout << ".";

			if (!isEqual(child, childExt)) {
				cout << "child table creation failed" << endl;
				return false;
			}

			cout << " OK!" << endl;

		}
		return true;
}

/*
int main() {
	return (testIndexCreation())? 0: 1;
}
*/

