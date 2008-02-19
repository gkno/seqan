#include <iostream>
#include <seqan/index.h>

using namespace std;
using namespace seqan;


	struct TMySpec;

	template <typename TString, typename TSpec>
	bool nodePredicate(Iter<Index<TString, Index_Wotd<TMySpec> >, TSpec> const &it) {
		return (countOccurrences(it) >= 2) && (repLength(it) >= 4);
	}


int main ()
{
///We begin with a @Class.String@ to store our sequence.
	String<char> myString = "How many wood would a woodchuck chuck.";

///Then we create an @Class.Index@ of this @Class.StringSet@
	typedef Index< String<char>, Index_Wotd<TMySpec> > TMyIndex;
	TMyIndex myIndex(myString);

///To find maximal repeats, we use SeqAn's @Spec.MaxRepeats Iterator@
///and set the minimum repeat length to 3.
	typedef Iterator< TMyIndex, TopDown<ParentLinks<> > >::Type TMaxRepeatIterator;
	TMaxRepeatIterator myRepeatIterator(myIndex);

	goBegin(myRepeatIterator);
	while (!atEnd(myRepeatIterator))
	{

///@Function.countOccurrences@ returns the number of hits of the representative.
		cout << countOccurrences(myRepeatIterator) << "x  ";

///The repeat string itself can be determined with @Function.representative@
		cout << "\t\"" << representative(myRepeatIterator) << '\"' << endl;

		goNext(myRepeatIterator);
	}

	return 0;
}
