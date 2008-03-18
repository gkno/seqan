#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <ctime>
#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_NOSRAN


#include <seqan/map.h>


#include <seqan/find.h>

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template <typename TMap>
void Test_NoCargo_Single() 
{
	typedef typename Value<TMap>::Type TValue;
	typedef typename Iterator<TMap>::Type TIterator;

	TMap map;

	SEQAN_TASSERT(length(map) == 0)

	insert(map, 2);
	SEQAN_TASSERT(map[2]);
	SEQAN_TASSERT(!map[5]);

	insert(map, 5);
	SEQAN_TASSERT(map[5]);

	insert(map, 1);
	SEQAN_TASSERT(hasKey(map, 1));

	int arr1[] = {1, 2, 5};
	TIterator it = begin(map);
	for (int i = 0; !atEnd(it); ++i)
	{
		SEQAN_TASSERT(i < 3)
		SEQAN_TASSERT(key(it) == arr1[i])
		goNext(it);
	}

	it = find(map, 2);
	erase(map, it);
	SEQAN_TASSERT(!map[2]);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TMap>
void Test_Cargo_Single() 
{
	typedef typename Value<TMap>::Type TValue;
	typedef typename Iterator<TMap>::Type TIterator;

	TMap map;

	SEQAN_TASSERT(length(map) == 0)

	insert(map, TValue(3, 1) );
	insert(map, TValue(1, 2) );
	insert(map, 10, 3);
	insert(map, 2, 4);

	SEQAN_TASSERT(length(map) == 4)

	int arr1[] = {2, 4, 1, 3};
	int arr2[] = {1, 2, 3, 10};
	TIterator it = begin(map);
	int i;
	for (i = 0; !atEnd(it); ++i)
	{
		SEQAN_TASSERT(i < 4)
		SEQAN_TASSERT(cargo(it) == arr1[i])
		SEQAN_TASSERT(key(it) == arr2[i])
		goNext(it);
	}

	map[8] = 5;
	map[2] = 6;

	SEQAN_TASSERT(mapValue(map, 8) == 5);

	int arr3[] = {2, 6, 1, 5, 3};
	int arr4[] = {1, 2, 3, 8, 10};
	i = 0;
	for (it = begin(map); it != end(map); ++it)
	{
		SEQAN_TASSERT(i < 5)
		SEQAN_TASSERT(cargo(it) == arr3[i])
		SEQAN_TASSERT(key(it) == arr4[i])
		++i;
	}

	it = find(map, 7);
	SEQAN_TASSERT(it)
	SEQAN_TASSERT(key(it) == 8)
	SEQAN_TASSERT(cargo(it) == 5)

	SEQAN_TASSERT(hasKey(map, 8));
	SEQAN_TASSERT(length(map) == 5)
	erase(map, it);
	SEQAN_TASSERT(!hasKey(map, 8));
	SEQAN_TASSERT(length(map) == 4)

	erase(map, 2);
	SEQAN_TASSERT(!hasKey(map, 2));
	SEQAN_TASSERT(length(map) == 3)

	erase(map, 8);
	SEQAN_TASSERT(length(map) == 3)

	clear(map);
	SEQAN_TASSERT(length(map) == 0)
}

//////////////////////////////////////////////////////////////////////////////

template <typename TMap>
void Test_Cargo_Multiple() 
{
	typedef typename Value<TMap>::Type TValue;
	typedef typename Iterator<TMap>::Type TIterator;

	TMap map;

	insert(map, 2, 4);
	SEQAN_TASSERT(length(map) == 1)

	insert(map, 2, 5);
	SEQAN_TASSERT(length(map) == 1)

	add(map, 2, 6);
	SEQAN_TASSERT(length(map) == 2)

	add(map, 3, 7);
	SEQAN_TASSERT(length(map) == 3)

	eraseAll(map, 2);
	SEQAN_TASSERT(length(map) == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Skiplist_Extra()
{
	typedef Pair<char, int> TValue;
	Map<TValue, Skiplist< > > map;

	map[8] = 5;
	SEQAN_TASSERT(value(map, 8) == TValue(8, 5) );

}


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")
//____________________________________________________________________________


	Test_Cargo_Single< Map< Pair<char, int>, Skiplist< > > >();
	Test_Cargo_Single< Map< Pair<char, int>, VectorSet< > > >();

	Test_NoCargo_Single< Map< char, Skiplist< > > >();
	Test_NoCargo_Single< Map< char, VectorSet< > > >();

	Test_Cargo_Multiple< Map< Pair<char, int>, Skiplist< > > >();

	Test_Skiplist_Extra();

//____________________________________________________________________________

	SEQAN_TREPORT("TEST END")

	return 0;
}
