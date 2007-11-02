#include <iostream>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

int main ()
{
///The following code builds a suffix array and searches in it
    Index< String<char> > index_esa("tobeornottobe");
	Finder< Index< String<char> > > finder_esa(index_esa);

	cout << "hit at ";
	while (find(finder_esa, "be"))
		cout << position(finder_esa) << " ";
	cout << endl;

///Doing the same using a q-gram index (q==2)
	typedef Index< String<char>, Index_QGram< FixedShape<2> > > TQGramIndex;
    TQGramIndex index_2gram("tobeornottobe");
	Finder< TQGramIndex > finder_2gram(index_2gram);

	cout << "hit at ";
	while (find(finder_2gram, "be"))
		cout << position(finder_2gram) << " ";
	cout << endl;

	return 0;
}
