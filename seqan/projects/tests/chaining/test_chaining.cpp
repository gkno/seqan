#include <time.h>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <memory>

//#define SEQAN_SL_DEBUG_PRINTOUT

//#define DEMOZEIT
//#define DEMOLAUF

#define SEQAN_DEBUG
#define SEQAN_TEST

//////////////////////////////////////////////////////////////////////////////

//todo: das muss man anders loesen

typedef int TKey;
typedef char TVal;

//////////////////////////////////////////////////////////////////////////////

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/chaining.h>


//////////////////////////////////////////////////////////////////////////////

#include "skip_list_test.h"

//#include "rt_entry.h"
//#include "rt_test_data.h"
//#include "range_tree_test.h"

//
//#include "range_max_tree_test.h"
//
//#include "rt_entry_test.h"

//////////////////////////////////////////////////////////////////////////////


using namespace seqan;
using namespace std;

//////////////////////////////////////////////////////////////////////////////
//todo: ausgliedern

template <typename TKey, typename TVal, typename TVal2>
inline TVal & 
setValue( std::pair< TKey, TVal > & p, TVal2 const & val_ )
{
	return p.second = val_;
}


//////////////////////////////////////////////////////////////////////////////

template <typename TModus, typename TSpec, typename TStructuring>
void 
testSkiplist()
{
	typedef Pair<int, char> TObject;
	typedef SkipList<TObject, TModus, TSpec, TStructuring> TSkipList;
	typedef typename Iterator<TSkipList>::Type TSkipListIterator;
	typedef String<TObject> TObjects;

	TObjects objs1;
	resize(objs1, 27);
	for (int i = 0; i < 26; ++i)
	{
		assignKey(value(objs1, i), rand());
		value(objs1, i).i2 = 'A' + i;
	}
	assignKey(value(objs1, 26), 10000);
	value(objs1, 26).i2 = '!';

	TSkipList sk1(begin(objs1), end(objs1));
	
	TSkipListIterator it1 = searchElement(sk1, 10000);
	cout << value(it1).i2;
}

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
	SEQAN_TREPORT("TEST BEGIN")

	testSkiplist<Dynamic, Default, Complete>();


	//testInsertsSkipList();
	//testEraseSkipList();
	//testInsertsErasesSkipList();
	//testSkipListInterface();
	//testSkipIterator< Static, Complete >();
	//testSkipIterator< Static, Deferred >();
	//testSkipIterator< Dynamic, Complete >();
	//testSkipIterator< Dynamic, Deferred >();
	//testRangeTree();


	//TODO: weitere Dateien
	debug::verifyCheckpoints("projects/library/seqan/chaining/skip_base_element.h");
	debug::verifyCheckpoints("projects/library/seqan/chaining/skip_element.h");
	debug::verifyCheckpoints("projects/library/seqan/chaining/skip_list_base.h");
	debug::verifyCheckpoints("projects/library/seqan/chaining/skip_list_dynamic.h");
	debug::verifyCheckpoints("projects/library/seqan/chaining/skip_list_impl.h");
	debug::verifyCheckpoints("projects/library/seqan/chaining/skip_list_iterator.h");
	debug::verifyCheckpoints("projects/library/seqan/chaining/skip_list_type.h");

	SEQAN_TREPORT("TEST END")
	

	//srand( time(0) );
	
/*
	//std::cout << "Hallo Welt" << std::endl;
	testDynamicSkipList();

	#ifdef DEMOLAUF
	for( int i = 0; i < 10; ++i )
	{
		std::cout << "Legende: " << std::endl;
		std::cout << "Schlüssel im Turm:       theKey" << std::endl;
		std::cout << "                         theKey" << std::endl;
		std::cout << "                        -------- "<< std::endl;
		std::cout << "Schlüssel des Elementes: theKey" << std::endl;
		std::cout << "Sortierung(x=ja/o=nein):  x " << std::endl;
		std::cout << "Höhe des Turmes:          2 " << std::endl;
		std::cout << "****************************** Dynamisch Deferrred: ***************************" << std::endl;
		testElementsDynamicDeferred();
		std::cout << "****************************** Dynamisch Complete: ***************************" << std::endl;
		testElementsDynamicComplete();
		std::cout << "****************************** Statisch Complete: ***************************" << std::endl;
		testElementsStaticComplete();
		std::cout << "****************************** Statisch Deferrred: ***************************" << std::endl;
		testElementsStaticDeferred();
		
	}
	#endif
	#ifdef DEMOZEIT
		timetest();
	#endif

*/
	return 0;
}

