#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

//#define SEQAN_DEBUG
//#define SEQAN_TEST

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


void testMultiIndex()
{
		typedef String<Dna5> TText;
		typedef StringSet< TText, ConcatDirect<> > TMulti;

		String<Dna5> t[6];
		//t[0] = "caterpillar";
		//t[1] = "catwoman";
		//t[2] = "pillow";
		//t[3] = "willow";
		//t[4] = "ill";
		//t[5] = "wow";

		t[0] = "caggctcgcgt";
		t[1] = "caggaacg";
		t[2] = "tcgttg";
		t[3] = "tggtcg";
		t[4] = "agg";
		t[5] = "ctg";

        Index<TMulti> esa;
		for(unsigned i=0; i<3; ++i)
			appendValue(indexText(esa), t[i]);

		Iter<Index<TMulti>, VSTree< BottomUp<> > > it(esa);
		while (!atEnd(it)) {
			cout << countOccurences(it) << "\t";
			cout << representative(it) << "\t";
//			cout << "parentEdgeLabel:" << parentEdgeLabel(it);
			cout << endl;
			goNext(it);
		}

/*
		indexRequire(esa, ESA_SA());
		indexRequire(esa, ESA_BWT());
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << saAt(i,esa) << " " << bwtAt(i,esa) << "    " << suffix(t[getValueI1(saAt(i,esa))], getValueI2(saAt(i,esa))) << endl;

//		resize(indexLCP(esa), length(indexRawText(esa)));
//		createLCPTableExt(indexLCP(esa), indexText(esa), indexSA(esa), Kasai());
		indexRequire(esa, ESA_LCP());
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << lcpAt(i,esa) << "    " << suffix(t[getValueI1(saAt(i,esa))], getValueI2(saAt(i,esa))) << endl;

		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << saAt(i,esa) << " = " << indexRawSA(esa)[i] << "    " << endl;
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << bwtAt(i,esa) << " = " << indexBWT(esa).tab[i] << "    " << endl;
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << lcpAt(i,esa) << " = " << indexLCP(esa)[i] << "    " << endl;
*/

/*

		resize(sa, length(indexRawText(esa)));
		createSuffixArrayExt(sa, indexText(esa), Skew7());

		for(int i=0; i<length(indexRawText(esa)); ++i)
			cout << indexRawText(esa)[i] << "    ";

		String<unsigned> lcp;
		resize(lcp, length(indexRawText(esa)));
		createLCPTableExt(lcp, indexText(esa), sa, Kasai());
*/
}



void testSTreeIterators()
{
		typedef Index<String<char>, Index_ESA<> > TIndex;

		String<char> text("acaaacatatz");
//		String<char> text("AAAAAGGGGG");
		TIndex index(text);
		Iter<TIndex, VSTree< TopDown< ParentLinks<Preorder> > > > it(index);
		Iter<TIndex, VSTree< TopDown<> > > itNoLinks(it);	// test conversion
		//Iter<TIndex, VSTree< BottomUp<> > > it(index);

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
			cout << countOccurences(it) << "\t";
			cout << representative(it) << "\t";
			cout << "parentEdgeLabel: " << parentEdgeLabel(it);
			cout << endl;
			goNext(it);
		}
}


void testMaxRepeats()
{
//		typedef String<char, External<> > TText;
		typedef String<char> TText;

        Index<TText> esa;
//        open(esa, "corpus/NC_000117.txt");
		indexText(esa) = "HALLOBALLOHALLEBALLO";

		FILE* dotFile = fopen("stree.dot","w");
		write(dotFile, esa, DotDrawing());
		fclose(dotFile);

        Iterator< Index<TText>, MaxRepeats >::Type it(esa, 3);
		typedef MaxRepeat< Index<TText> > TRepeat;
        while (!atEnd(it)) {
			cout << representative(it) << ":";
			Iterator<TRepeat>::Type mit(it);
			while (!atEnd(mit)) {
				cout << "\t" << *mit;
				++mit;
			}
			cout << endl;
            ++it;
        }
}


template <typename TIteratorSpec>
void testSuperMaxRepeats()
{
//		typedef String<char, External<> > TText;
		typedef String<char> TText;

        Index<TText> esa;
//        open(esa, "corpus/NC_000117.txt");
		indexText(esa) = "HALLOBALLOHALLEBALLO";

        typename Iterator< Index<TText>, TIteratorSpec >::Type it(esa);
        while (!atEnd(it)) {
			cout << representative(it) << ":";
			for(typename Size<Index<TText> >::Type i = 0; i < countOccurences(it); ++i)
				cout << "\t" << getOccurences(it)[i];
			cout << endl;
            ++it;
        }
}


void testMUMs()
{
		typedef String<char> TText;
		typedef StringSet< TText, ConcatDirect<> > TMulti;
		typedef Index<TMulti, Index_ESA<> > TIndex;

		String<char> t[3];

		t[0] = "fefhalloballo";
		t[1] = "halloballefser";
		t[2] = "grballoballo";


        TIndex esa;
		for(int i = 0; i < 3; ++i)
			appendValue(indexText(esa), t[i]);			// add sequences to multiple index

		Iterator<TIndex, MUMs>::Type it(esa, 3);		// set minimum MUM length to 3
		String< SAValue<TIndex>::Type > occs;			// temp. string storing the hit positions

		cout << resetiosflags(ios::left);
		while (!atEnd(it)) 
		{
			cout << representative(it) << ":";
			occs = getOccurences(it);					// gives hit positions (seqNo,seqOfs)
			orderOccurences(occs);						// order them by seqNo
			
			for(unsigned i = 0; i < length(occs); ++i)
				cout << "\t" << getValueI2(occs[i]);
			cout << endl;
			cout << alignment(it) << endl;

			++it;
		}
}


void simpleMUMmer()
{
		typedef String<char> TText;
		typedef StringSet< TText, ConcatDirect<> > TMulti;
		typedef Index<TMulti, Index_ESA<> > TIndex;

		String<char> t[3];
/*		if (!open(t[0], "corpus/NC_000117.txt") ||
		    !open(t[1], "corpus/NC_002620.txt") ||
			!open(t[2], "corpus/NC_007429.txt")) 
		{
			cout << "Could not read in all datasets" << endl;
			return;
		}
*/

		t[0] = "fefhalloballo";
		t[1] = "halloballefser";
		t[2] = "grballoballo";


        TIndex esa;
		for(int i = 0; i < 3; ++i)
			appendValue(indexText(esa), t[i]);			// add sequences to multiple index

		Iterator<TIndex, MUMs>::Type it(esa, 20);		// set minimum MUM length to 20
		String< SAValue<TIndex>::Type > occs;			// temp. string storing the hit positions

		cout << resetiosflags(ios::left);
		while (!atEnd(it)) 
		{
			occs = getOccurences(it);					// gives hit positions (seqNo,seqOfs)
			orderOccurences(occs);						// order them by seqNo
			
			for(unsigned i = 0; i < length(occs); ++i)
			   	cout << setw(8)
				     << getValueI2(occs[i])	+ 1			// output them in MUMmer's output format
					 << "  ";

			cout << setw(8) 
				 << repLength(it)
				 << endl;

			cout << "== " << alignment(it) << endl;

			++it;
		}
		cout << setiosflags(ios::left) << endl;
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


bool testIndexCreation();
void Main_TestQGram();

int main()
{
	SEQAN_TREPORT("TEST BEGIN")

	std::cout << "===================================" << ::std::endl;
	std::cout << "----Basic Suffix Tree iterators----" << ::std::endl;
	testSTreeIterators();
//	testIndexCreation();
//	Main_TestQGram();
//	testBuild();

	std::cout << "===================================" << ::std::endl;
	std::cout << "----Suffix Array based Finder------" << ::std::endl;
	testFind<ESA_FIND_MLR>();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Multiple Sequence Index--------" << ::std::endl;
	testMultiIndex();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----MUMs---------------------------" << ::std::endl;
	testMUMs();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Maximal Repeats----------------" << ::std::endl;
	testMaxRepeats();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Supermaximal Repeats-----------" << ::std::endl;
	testSuperMaxRepeats<SuperMaxRepeats>();
	std::cout << "===================================" << ::std::endl;
	std::cout << "----Supermaximal Repeats---fast----" << ::std::endl;
	testSuperMaxRepeats<SuperMaxRepeatsFast>();
	std::cout << "===================================" << ::std::endl;

//		simpleMUMmer();	
	SEQAN_TREPORT("TEST END")
		return 0;
}
