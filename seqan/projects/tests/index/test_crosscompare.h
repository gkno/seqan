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

#ifndef TESTS_INDEX_TEST_CROSS_COMPARE_H
#define TESTS_INDEX_TEST_CROSS_COMPARE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


template <typename TIter>
inline void _dumpState(TIter const it) 
{
	std::cout << typeid(it).name() << std::endl;
	std::cout << "  range:            " << value(it).range << std::endl;
	std::cout << "  countOccurrences: " << countOccurrences(it) << std::endl;
	std::cout << "  representative:   " << representative(it)   << std::endl;
	std::cout << "  parentEdgeLabel:  " << parentEdgeLabel(it)  << std::endl;
	std::cout << "  Occurrences:      ";
	for(unsigned i=0;i<countOccurrences(it);++i)
		std::cout << getOccurrences(it)[i] << " ";
	std::cout << std::endl;
}

template <typename TIterSpec, typename TIndex1, typename TIndex2>
void crossBottomUp(TIndex1 &index1, TIndex2 &index2)
{
	typename Iterator<TIndex1, TIterSpec>::Type iter1(index1);
	typename Iterator<TIndex2, TIterSpec>::Type iter2(index2);

	while (!atEnd(iter1) && !atEnd(iter2)) 
	{
		SEQAN_ASSERT_EQ(representative(iter1), representative(iter2));
		SEQAN_ASSERT_EQ(countOccurrences(iter1), countOccurrences(iter2));
		SEQAN_ASSERT_EQ(parentEdgeLength(iter1), parentEdgeLength(iter2));
		goNext(iter1);
		goNext(iter2);
	}
	SEQAN_ASSERT_EQ(atEnd(iter1), atEnd(iter2));
}

template <typename TIndexSpec1, typename TIndexSpec2, typename TText>
void crossIndex(TText &text)
{
	Index<TText, TIndexSpec1> index1(text);
	Index<TText, TIndexSpec2> index2(text);

/*	crossBottomUp< TopDown<ParentLinks<Preorder> > > (index1, index2);
	crossBottomUp< TopDown<ParentLinks<Postorder> > > (index1, index2);
*/	crossBottomUp< TopDown<ParentLinks<PreorderEmptyEdges> > > (index1, index2);
	crossBottomUp< TopDown<ParentLinks<PostorderEmptyEdges> > > (index1, index2);
}

template <typename TIndexSpec1, typename TIndexSpec2>
void crossIndices()
{
/*	{
		CharString text("mississippi");
		crossIndex<TIndexSpec1,TIndexSpec2> (text);
	}
	{
		DnaString text("acaaacatat");
		crossIndex<TIndexSpec1,TIndexSpec2> (text);
	}
*/	{
		StringSet<CharString> t;
		resize(t, 6);
		t[0] = "caterpillar";
		t[1] = "catwoman";
		t[2] = "pillow";
		t[3] = "willow";
		t[4] = "ill";
		t[5] = "wow";
		crossIndex<TIndexSpec1,TIndexSpec2> (t);
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
		crossIndex<TIndexSpec1,TIndexSpec2> (t);
	}
}

SEQAN_DEFINE_TEST(testIndexCrossCompare)
{
	crossIndices<Index_ESA<>, Index_Wotd<> >();
	crossIndices<Index_Wotd<>, Index_Wotd<DFI<> > >();
}

//////////////////////////////////////////////////////////////////////////////


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
