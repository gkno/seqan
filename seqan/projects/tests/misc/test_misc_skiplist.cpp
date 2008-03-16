#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <ctime>
#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_NOSRAN

#include <seqan/basic.h>

#include <seqan/misc/misc_skiplist.h>

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

void Main_TestSkiplist() 
{
	SEQAN_TREPORT("TEST SKIPLIST BEGIN")
//____________________________________________________________________________


	typedef Skiplist<int, int> TSkiplist;
	TSkiplist sl;

	SEQAN_TASSERT(length(sl) == 0)

	insertElement(sl, 3, 1);
	insertElement(sl, 1, 2);
	insertElement(sl, 10, 3);
	insertElement(sl, 2, 4);

	SEQAN_TASSERT(length(sl) == 4)

	int arr1[] = {2, 4, 1, 3};
	int arr2[] = {1, 2, 3, 10};
	typedef Iterator<TSkiplist>::Type TIterator;
	TIterator it = begin(sl);
	for (int i = 0; !atEnd(it); ++i)
	{
		SEQAN_TASSERT(i < 4)
		SEQAN_TASSERT(value(it) == arr1[i])
		SEQAN_TASSERT(key(it) == arr2[i])
		goNext(it);
	}

	sl[8] = 5;
	sl[2] = 6;

	int arr3[] = {2, 6, 1, 5, 3};
	int arr4[] = {1, 2, 3, 8, 10};
	i = 0;
	for (it = begin(sl); it != end(sl); ++it)
	{
		SEQAN_TASSERT(i < 5)
		SEQAN_TASSERT(value(it) == arr3[i])
		SEQAN_TASSERT(key(it) == arr4[i])
		++i;
	}

	it = findElement(sl, 7);
	SEQAN_TASSERT(it)
	SEQAN_TASSERT(key(it) == 8)
	SEQAN_TASSERT(value(it) == 5)

	SEQAN_TASSERT(hasKey(sl, 8));
	SEQAN_TASSERT(length(sl) == 5)
	removeElement(sl, it);
	SEQAN_TASSERT(!hasKey(sl, 8));
	SEQAN_TASSERT(length(sl) == 4)

//____________________________________________________________________________

	SEQAN_TREPORT("TEST SKIPLIST END")
}
