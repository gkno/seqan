/*==========================================================================
 SeqAn - The Library for Sequence Analysis
 http://www.seqan.de 
 ============================================================================
 Copyright (C) 2007
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.
 
 ============================================================================
 Author: David Weese <david.weese@fu-berlin.de>
 ==========================================================================*/

#ifndef TESTS_INDEX_TEST_STREE_ITERATORS_H
#define TESTS_INDEX_TEST_STREE_ITERATORS_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

SEQAN_DEFINE_TEST(testBuild)
{
		typedef String<char> TText;
		typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;

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

template <typename TIter>
inline void _printNode(TIter const it) 
{
		std::cout << countOccurrences(it) << "\t";
		std::cout << representative(it) << "\t";
//		std::cout << "parentEdgeLabel:" << parentEdgeLabel(it);
		std::cout << std::endl;
}

SEQAN_DEFINE_TEST(testMultiIndex)
{
		typedef String<Dna5> TText;
		typedef StringSet< TText, Owner<> > TMulti;

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

		t[0] = "ac";
		t[1] = "ac";
//		t[2] = "aatt";

        Index<TMulti> esa;
		for(unsigned i=0; i<2; ++i)
			appendValue(indexText(esa), t[i]);

		// efficient dfs iterator (hiding edges with empty labels)
		{
			std::cout << "BottomUp without empty edges" << std::endl;
			Iter<Index<TMulti>, VSTree< BottomUp<> > > it(esa);
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// efficient dfs iterator
		{
			std::cout << std::endl << "BottomUp with empty edges" << std::endl;
			Iter<Index<TMulti>, VSTree< BottomUp<PostorderEmptyEdges> > > it(esa);
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown dfs iterator (hiding edges with empty labels)
		{
			std::cout << std::endl << "TopDown postorder without empty edges" << std::endl;
			Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<Postorder> > > > it(esa);
			while (goDown(it)) ;
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown dfs iterator
		{
			std::cout << std::endl << "TopDown postorder with empty edges" << std::endl;
			Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<PostorderEmptyEdges> > > > it(esa);
			while (goDown(it)) ;
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown dfs iterator (hiding edges with empty labels)
		{
			std::cout << std::endl << "TopDown preorder without empty edges" << std::endl;
			Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<Preorder> > > > it(esa);
			goDown(it,'c');
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown dfs iterator
		{
			std::cout << std::endl << "TopDown preorder with empty edges" << std::endl;
			Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<PreorderEmptyEdges> > > > it(esa);
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown iterator w/o parent links (hiding edges with empty labels)
		{
			std::cout << std::endl << "TopDown with empty edges" << std::endl;
			Iter<Index<TMulti>, VSTree< TopDown<HideEmptyEdges> > > it(esa);
			_printNode(it);
			while (goDown(it))
				_printNode(it);
		}

		// topdown iterator w/o parent links
		{
			std::cout << std::endl << "TopDown with empty edges" << std::endl;
			Iter<Index<TMulti>, VSTree< TopDown<EmptyEdges> > > it(esa);
			_printNode(it);
			while (goDown(it))
				_printNode(it);
		}
/*
		indexRequire(esa, ESA_SA());
		indexRequire(esa, ESA_BWT());
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			std::cout << saAt(i,esa) << " " << bwtAt(i,esa) << "    " << suffix(t[getValueI1(saAt(i,esa))], getValueI2(saAt(i,esa))) << std::endl;

//		resize(indexLCP(esa), length(indexRawText(esa)));
//		createLCPTableExt(indexLCP(esa), indexText(esa), indexSA(esa), Kasai());
		indexRequire(esa, ESA_LCP());
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			std::cout << lcpAt(i,esa) << "    " << suffix(t[getValueI1(saAt(i,esa))], getValueI2(saAt(i,esa))) << std::endl;

		for(int i=0; i<length(indexRawSA(esa)); ++i)
			std::cout << saAt(i,esa) << " = " << indexRawSA(esa)[i] << "    " << std::endl;
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			std::cout << bwtAt(i,esa) << " = " << indexBWT(esa).tab[i] << "    " << std::endl;
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			std::cout << lcpAt(i,esa) << " = " << indexLCP(esa)[i] << "    " << std::endl;
*/

/*

		resize(sa, length(indexRawText(esa)));
		createSuffixArrayExt(sa, indexText(esa), Skew7());

		for(int i=0; i<length(indexRawText(esa)); ++i)
			std::cout << indexRawText(esa)[i] << "    ";

		String<unsigned> lcp;
		resize(lcp, length(indexRawText(esa)));
		createLCPTableExt(lcp, indexText(esa), sa, Kasai());
*/
}

template <typename TIndex1, typename TIndex2>
void compareTreeIterators(TIndex1 &index1, TIndex2 &index2)
{
	Iter<TIndex1, VSTree< TopDown< ParentLinks<Preorder> > > > it1(index1);
	Iter<TIndex2, VSTree< TopDown< ParentLinks<Preorder> > > > it2(index2);

	while (!atEnd(it1) && !atEnd(it2)) 
	{
		SEQAN_ASSERT_EQ(representative(it1), representative(it2));
		SEQAN_ASSERT_EQ(parentEdgeLabel(it1), parentEdgeLabel(it2));
		goNext(it1);
		goNext(it2);
	}

	SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it1));
}

template <typename TIndexSpec1, typename TIndexSpec2>
void compareIndices()
{
	{
		CharString text("mississippi");
		Index<CharString, TIndexSpec1> index1(text);
		Index<CharString, TIndexSpec2> index2(text);
		compareTreeIterators(index1, index2);
	}
	{
		DnaString text("acaaacatat");
		Index<DnaString, TIndexSpec1> index1(text);
		Index<DnaString, TIndexSpec2> index2(text);
		compareTreeIterators(index1, index2);
	}
	{
		StringSet<CharString> t;
		resize(t, 6);
		t[0] = "caterpillar";
		t[1] = "catwoman";
		t[2] = "pillow";
		t[3] = "willow";
		t[4] = "ill";
		t[5] = "wow";
		Index<StringSet<CharString>, TIndexSpec1> index1(t);
		Index<StringSet<CharString>, TIndexSpec2> index2(t);
		compareTreeIterators(index1, index2);
	}
	{
		StringSet<DnaString> t;
		resize(t, 6);
		t[0] = "caggctcgcgt";
		t[1] = "caggaacg";
		t[2] = "tcgttg";
		t[3] = "tggtcg";
		t[4] = "agg";
		t[5] = "ctg";
		Index<StringSet<DnaString>, TIndexSpec1> index1(t);
		Index<StringSet<DnaString>, TIndexSpec2 > index2(t);
		compareTreeIterators(index1, index2);
	}
}

SEQAN_DEFINE_TEST(testCompareIndices_Esa_Wotd)
{
	compareIndices<Index_ESA<>, Index_Wotd<> >();
}


template <typename TIndexSpec>
void testSTreeIterators()
{
		typedef Index<String<char>, TIndexSpec> TIndex;

		String<char> text("acaaacatatz");
//		String<char> text("AAAAAGGGGG");
		TIndex index(text);
		Iter<TIndex, VSTree< TopDown< ParentLinks<Preorder> > > > it(index);
		Iter<TIndex, VSTree< TopDown<> > > itNoLinks(it);	// test conversion
		//Iter<TIndex, VSTree< BottomUp<> > > it(index);

//		while (goDown(it));
		while (!atEnd(it)) {
//			std::cout << countOccurrences(it) << "\t";
			std::cout << representative(it) << "\t";
			std::cout << "parentEdgeLabel: " << parentEdgeLabel(it); // << " " << value(it).node << "  " << value(it).range;
			std::cout << std::endl;
			goNext(it);
		}
			std::cout << std::endl;
		_dump(index);
/*		goBegin(it);
		while (!atEnd(it)) {
			std::cout << countOccurrences(it) << "\t";
			std::cout << representative(it) << "\t";
			std::cout << "parentEdgeLabel: " << parentEdgeLabel(it) << " " << value(it).node << "  " << value(it).range;
			std::cout << std::endl;
			goNext(it);
		}
		_dump(index);
*/}

SEQAN_DEFINE_TEST(testSTreeIterators_Wotd)
{
	testSTreeIterators<Index_Wotd<> >();
}

SEQAN_DEFINE_TEST(testSTreeIterators_WotdOriginal)
{
	testSTreeIterators<Index_Wotd<WotdOriginal> >();
}

SEQAN_DEFINE_TEST(testSTreeIterators_Esa)
{
	testSTreeIterators<Index_ESA<> >();
}


SEQAN_DEFINE_TEST(testMaxRepeats)
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
//			std::cout << representative(it) << ":";
			Iterator<TRepeat, MaxRepeatOccurrences>::Type mit(it);
			while (!atEnd(mit)) {
//				std::cout << "\t" << *mit;
				++mit;
			}
//			std::cout << std::endl;
            ++it;
        }
}


SEQAN_DEFINE_TEST(testMultiMEMs)
{
		typedef String<char> TText;
		typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;

        Index<TMulti> esa;

		String<Dna5> t[6];
		t[0] = "caterpillar";
		t[1] = "catwoman";
		t[2] = "pillow";
		t[3] = "willow";
		t[4] = "ill";
		t[5] = "wow";

		FILE* dotFile = fopen("stree.dot","w");
		write(dotFile, esa, DotDrawing());
		fclose(dotFile);

        Iterator< Index<TMulti>, MultiMEMs >::Type it(esa, 3);
		typedef MultiMEM< Index<TMulti> > TMultiMEM;
        while (!atEnd(it)) {
			std::cout << representative(it) << ":";
			Iterator<TMultiMEM>::Type mit(it);
			while (!atEnd(mit)) {
				std::cout << "\t" << *mit;
				++mit;
			}
			std::cout << std::endl;
            ++it;
        }
}

template <typename TIteratorSpec>
void _testSuperMaxRepeats()
{
//		typedef String<char, External<> > TText;
		typedef String<char> TText;

        Index<TText> esa;
//        open(esa, "corpus/NC_000117.txt");
		indexText(esa) = "HALLOBALLOHALLEBALLO";

        typename Iterator< Index<TText>, TIteratorSpec >::Type it(esa);
        while (!atEnd(it)) {
//			std::cout << representative(it) << ":";
//			for(typename Size<Index<TText> >::Type i = 0; i < countOccurrences(it); ++i)
//				std::cout << "\t" << getOccurrences(it)[i];
//			std::cout << std::endl;
            ++it;
        }
}

SEQAN_DEFINE_TEST(testSuperMaxRepeats)
{
	_testSuperMaxRepeats<SuperMaxRepeats>();
}

SEQAN_DEFINE_TEST(testSuperMaxRepeatsFast)
{
	_testSuperMaxRepeats<SuperMaxRepeatsFast>();
}


SEQAN_DEFINE_TEST(testMUMs)
{
		typedef String<char> TText;
		typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;
		typedef Index<TMulti, Index_ESA<> > TIndex;

		String<char> t[3];

		t[0] = "fefhalloballo";
		t[1] = "halloballefser";
		t[2] = "grballoballo";


        TIndex esa;
		for(int i = 0; i < 3; ++i)
			appendValue(indexText(esa), t[i]);			// add sequences to multiple index

		Iterator<TIndex, MUMs>::Type  it(esa, 3);		// set minimum MUM length to 3
		String< SAValue<TIndex>::Type > occs;			// temp. string storing the hit positions

		std::cout << std::resetiosflags(std::ios::left);
		while (!atEnd(it)) 
		{
//			std::cout << representative(it) << ":";
			occs = getOccurrences(it);					// gives hit positions (seqNo,seqOfs)
			orderOccurrences(occs);						// order them by seqNo
			
//			for(unsigned i = 0; i < length(occs); ++i)
//				std::cout << "\t" << getValueI2(occs[i]);
//			std::cout << std::endl;
//			std::cout << alignment(it) << std::endl;

			++it;
		}
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
		seqan::Pattern<String<char> > pattern(needle1);	

		while (find(finder, pattern))
			append(pos,position(finder));

		SEQAN_ASSERT_TRUE(pos[0] == 5);
		SEQAN_ASSERT_TRUE(pos[1] == 31);

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

		SEQAN_ASSERT_TRUE(pos[1] == 0);
		SEQAN_ASSERT_TRUE(pos[0] == 26);
}

SEQAN_DEFINE_TEST(testFind_Esa_Mlr)
{
	testFind<ESA_FIND_MLR>();
}

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
