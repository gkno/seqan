#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/index.h>
#include "test_index_creation.h"


using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////



void testBuild()
{
		typedef String<char> TText;
		typedef StringSet< TText, ConcatDirect<> > TMulti;

		String<char> gen1, gen2, gen3;
		std::cout << open(gen1, "corpus/NC_000117.txt");
		std::cout << open(gen2, "corpus/NC_002620.txt");
		std::cout << open(gen3, "corpus/NC_007429.txt");

        Index<TMulti> esa;
		appendValue(indexText(esa), gen1);
		appendValue(indexText(esa), gen2);
		appendValue(indexText(esa), gen3);



		indexRequire(esa, ESA_SA());
		indexRequire(esa, ESA_LCP());
		indexRequire(esa, ESA_BWT());
		indexRequire(esa, ESA_ChildTab());

        save(esa, "corpus/chlamydia");
}
