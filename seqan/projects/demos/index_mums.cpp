#include <iostream>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

int main ()
{
///We begin with a @Class.StringSet@ that stores multiple strings.
	StringSet< String<char> > mySet;
	resize(mySet, 3);
	mySet[0] = "SeqAn is a library for sequence analysis.";
	mySet[1] = "The String class is the fundamental sequence type in SeqAn.";
	mySet[2] = "Subsequences can be handled with SeqAn's Segment class.";

///Then we create an @Class.Index@ of this @Class.StringSet@.
	typedef Index< StringSet<String<char> > > TMyIndex;
	TMyIndex myIndex(mySet);

///To find maximal unique matches (MUMs),
///we use an @Spec.MUMs Iterator.iterator@ of the index specialized for the search of MUMs.
///We set the minimum MUM length to 3.
	Iterator< TMyIndex, MUMs >::Type myMUMiterator(myIndex, 3);
	String< SAValue<TMyIndex>::Type > occs;

	while (!atEnd(myMUMiterator)) 
	{
///A multiple match can be represented by positions it occurs at in every sequence
///and its length. @Function.getOccurrences@ returns an unordered sequence of pairs
///(seqNo,seqOfs) the match occurs at.
		occs = getOccurrences(myMUMiterator);
///To order them ascending according seqNo we use @Function.orderOccurrences@.
		orderOccurrences(occs);
		
		for(unsigned i = 0; i < length(occs); ++i)
			cout << getValueI2(occs[i]) << ", ";

///@Function.repLength@ returns the length of the match.
		cout << repLength(myMUMiterator) << "   ";

///The match string itself can be determined with @Function.representative@.
		cout << '\"' << representative(myMUMiterator) << '\"' << endl;

		++myMUMiterator;
	}

	return 0;
}

	//output:
	//
	// 0, 53, 45, 5   "SeqAn"
	// 23, 36, 3, 8   "sequence"
