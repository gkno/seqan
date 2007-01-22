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


void testSTreeIterators()
{
		typedef Index<String<char>, Index_ESA<> > TIndex;

		String<char> text("acaaacatatz");
		TIndex index(text);
		//Iter<TIndex, VSTree< TopDown< ParentLinks<> > > > it(index);
		Iter<TIndex, VSTree< BottomUp<> > > it(index);

		cout << "SA       = ";
		for(unsigned i=0; i<length(indexSA(index)); ++i)
			cout << saAt(i, index) << "  ";
		cout << endl;

		cout << "LCP      = ";
		for(unsigned i=0; i<length(indexLCP(index)); ++i)
			cout << lcpAt(i, index) << "  ";
		cout << endl;

		cout << "ChildTab = ";
		for(unsigned i=0; i<length(indexChildTab(index)); ++i)
			cout << childAt(i, index) << "  ";
		cout << endl;

//		while (goDown(it));
		while (!atEnd(it)) {
			std::cout << value(it) << " = " << _repLength(it) << "  " << representative(it) << std::endl;
			goNext(it,Postorder());
		}
}


void testSuperMaxRepeats()
{

		typedef String<char, External<> > TText;

        Index<TText> esa;
        open(esa, "corpus/chlamydia");

        Iterator< Index<TText>, SuperMaxRepeats >::Type it(esa, 20);
        unsigned counter = 0;
        while (!atEnd(it)) {
                ++it;
                ++counter;
        }

        std::cout << "supermaximal repeats: " << counter << std::endl;

}


void testMUMs()
{
/*		typedef StringSet< String<char, External<> > > TText;

        Index<TText> esa;
        open(esa, "corpus/chr1");

        Iterator< Index<TText>, MUMs >::Type it(esa, 20);
        unsigned counter = 0;
        while (!atEnd(it)) {
                ++it;
                ++counter;
        }

        std::cout << "supermaximal repeats: " << counter << std::endl;
*/
}


template <typename TAlgorithmSpec>
void testFind()
{
		String<unsigned int> pos;

	//____________________________________________________________________________
	// Test1 - small needle

		String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
		Index<String<char> > index(haystack);

		Finder<Index< String<char> >, TAlgorithmSpec> finder(index);

		String<char> needle1("ist");
		Pattern<String<char>,void > pattern(needle1);

		while (find(finder, pattern))
			append(pos,position(finder));

		SEQAN_TASSERT(pos[0] == 5);
		SEQAN_TASSERT(pos[1] == 31);

	//____________________________________________________________________________
	// Test2 - large needle

		haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
		clear(index);
		clear(finder);

		needle1 = "abcdefghijklmnopqrstuvwxyzabcdefg";
		setNeedle(pattern, needle1);

		clear(pos);
		while (find(finder, pattern))
			append(pos,position(finder));

		SEQAN_TASSERT(pos[1] == 0);
		SEQAN_TASSERT(pos[0] == 26);
}


void Main_TestQGram();


int main()
{
	SEQAN_TREPORT("TEST BEGIN")

		Main_TestQGram();
//		testFind<ESA_MLR>();
//		testBuild();
		testSTreeIterators();
/*		testSuperMaxRepeats();
		testMUMs();	
*/
	SEQAN_TREPORT("TEST END")
		return 0;
}
