#include <iostream>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

/// constraint parameters
struct TMyConstraints {
	double p_min;
	int replen_max;
};

/// SeqAn extensions
namespace seqan 
{
	// custom TSpec for our customized wotd-Index
	struct TMyIndex;

	template <typename TText>
	struct Cargo<Index<TText, Index_Wotd<TMyIndex> > > {
		typedef TMyConstraints Type;
	};

	// node predicate
	template <typename TText, typename TSpec>
	bool nodePredicate(Iter<Index<TText, Index_Wotd<TMyIndex> >, TSpec> const &it) 
	{
		TMyConstraints const &cons = cargo(container(it));
		unsigned replen = repLength(it);
		unsigned textLen = length(container(it));

		if (textLen == replen) return false;
		return ((double)countOccurrences(it) / 
				(double)(textLen - replen)) > cons.p_min;
	}

	// monotonic hull
	template <typename TText, typename TSpec>
	bool nodeHullPredicate(Iter<Index<TText, Index_Wotd<TMyIndex> >, TSpec> const &it) 
	{
		Index<TText, Index_Wotd<TMyIndex> > const &index = container(it);
		TMyConstraints const &cons = cargo(index);

		if (repLength(it) > cons.replen_max) return false;
		return ((double)countOccurrences(it) / 
				(double)(length(index) - cons.replen_max)) > cons.p_min;
	}
}

int main ()
{
///We begin with a @Class.String@ to store our sequence.
	String<char> myString = "How many wood would a woodchuck chuck.";

///Then we create our customized index which is a specialization
///of the deferred @Class.Index.wotd-Index@
	typedef Index< String<char>, Index_Wotd<TMyIndex> > TMyIndex;
	TMyIndex myIndex(myString);

	cargo(myIndex).replen_max = 10;
	cargo(myIndex).p_min = 0.05;

///To find all strings that fulfill our constraints,
///we simply do a dfs-traversal via @Function.goBegin@ and @Function.goNext@
	typedef Iterator< TMyIndex, TopDown<ParentLinks<> > >::Type TConstrainedIterator;
	TConstrainedIterator myConstrainedIterator(myIndex);

	goBegin(myConstrainedIterator);
	while (!atEnd(myConstrainedIterator))
	{

///@Function.countOccurrences@ returns the number of hits of the representative.
		cout << countOccurrences(myConstrainedIterator) << "x  ";

///The representative string can be determined with @Function.representative@
		cout << "\t\"" << representative(myConstrainedIterator) << '\"' << endl;

		goNext(myConstrainedIterator);
	}

	return 0;
}
