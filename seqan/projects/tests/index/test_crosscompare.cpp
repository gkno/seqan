#include <seqan/platform.h>

#include <iostream>
#include <fstream>

//#define SEQAN_DEBUG
//#define SEQAN_TEST

#include "test_index_creation.h"
#include <seqan/align.h>
#include <seqan/index.h>
#include <typeinfo>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////



template <typename TIter>
inline void _dumpState(TIter const it) 
{
	cout << typeid(it).name() << endl;
	cout << "  range:            " << value(it).range << endl;
	cout << "  countOccurrences: " << countOccurrences(it) << endl;
	cout << "  representative:   " << representative(it)   << endl;
	cout << "  parentEdgeLabel:  " << parentEdgeLabel(it)  << endl;
	cout << "  Occurrences:      ";
	for(unsigned i=0;i<countOccurrences(it);++i)
		cout << getOccurrences(it)[i] << " ";
	cout << endl;
}

template <typename TIterSpec, typename TIndex1, typename TIndex2>
bool crossBottomUp(TIndex1 &index1, TIndex2 &index2)
{
	typename Iterator<TIndex1, TIterSpec>::Type iter1(index1);
	typename Iterator<TIndex2, TIterSpec>::Type iter2(index2);

	while (!atEnd(iter1) && !atEnd(iter2)) 
	{
		if (value(iter1).range.i1 == 9 && value(iter1).range.i2 == 11) {
			cerr << endl;
			_dump(index2);
			cerr << endl;
		}

		if (!( representative(iter1) == representative(iter2) 
			&& countOccurrences(iter1) == countOccurrences(iter2)
			&& parentEdgeLength(iter1) == parentEdgeLength(iter2)))
		{
			cerr << "Iterators differ !!!" << endl;
			_dumpState(iter1);
			cerr << endl << "!=" << endl;
			_dumpState(iter2);
			cerr << endl;
			_dump(index2);
//			return false;
		} else {
			cerr << value(iter1).range << "\t  " << representative(iter1) << endl;
		}
		goNext(iter1);
		goNext(iter2);
	}
	if (!(atEnd(iter1) && atEnd(iter2)))
	{
		cerr << "Only one iterator reached the end !!!" << endl;
		if (!atEnd(iter1))
			_dumpState(iter1);
		if (!atEnd(iter2))
			_dumpState(iter2);
		return false;
	}
	return true;
}

template <typename TIndexSpec1, typename TIndexSpec2, typename TText>
bool crossIndex(TText &text)
{
	Index<TText, TIndexSpec1> index1(text);
	Index<TText, TIndexSpec2> index2(text);

	bool equal = true;
/*	equal &= crossBottomUp< TopDown<ParentLinks<Preorder> > > (index1, index2);
	equal &= crossBottomUp< TopDown<ParentLinks<Postorder> > > (index1, index2);
*/	equal &= crossBottomUp< TopDown<ParentLinks<PreorderEmptyEdges> > > (index1, index2);
	equal &= crossBottomUp< TopDown<ParentLinks<PostorderEmptyEdges> > > (index1, index2);

	return equal;
}

template <typename TIndexSpec1, typename TIndexSpec2>
bool crossIndices()
{
	bool equal = true;
/*	{
		CharString text("mississippi");
		equal &= crossIndex<TIndexSpec1,TIndexSpec2> (text);
	}
	{
		DnaString text("acaaacatat");
		equal &= crossIndex<TIndexSpec1,TIndexSpec2> (text);
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
		equal &= crossIndex<TIndexSpec1,TIndexSpec2> (t);
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
		equal &= crossIndex<TIndexSpec1,TIndexSpec2> (t);
	}
	return equal;
}

bool crossCompare()
{
	bool equal = true;
	equal &= crossIndices<Index_ESA<>, Index_Wotd<> > ();
	equal &= crossIndices<Index_Wotd<>, Index_Wotd<DFI<> > > ();
	return equal;
}

/*
int main() {
	return crossCompare();
}
*/
