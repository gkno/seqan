#include <time.h>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>

//#define SEQAN_SL_DEBUG_PRINTOUT

//#define DEMOZEIT
//#define DEMOLAUF

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/chaining.h>


typedef int TKey;
typedef char TVal;

using namespace seqan;
using namespace std;

/*
// Zeittest für dynamische SkipListe
template<typename TKey, typename TVal>
void 
checkSL_def_dyn(	std::vector< std::pair<TKey, TVal> > & pairs,
					TKey * search_keys,
					TKey searches )
{
	SkipList< std::pair<TKey, TVal>, Dynamic, Default, Deferred > SL( begin( pairs ), end( pairs ) );
	volatile TKey vol;
	for( TKey i = 0; i < searches; ++i){
		vol = key( searchElement( SL, search_keys[i] ) );
	}
}



// Zeittest für statische SkipListe

template<typename TKey, typename TVal>
void 
checkSL_def_stat(	std::vector< std::pair<TKey, TVal> > & pairs,
					TKey * search_keys,
					TKey searches )
{
	SkipList< std::pair<TKey, TVal>, Static, Default, Deferred > SL( begin( pairs ), end( pairs ) );
	volatile TKey vol;
	for( TKey i = 0; i < searches; ++i ){
		SkipBaseElement< std::pair<TKey, TVal>, Static, Default, Deferred > * elem = _pointsTo( searchElement( SL, search_keys[i] ) );
		vol = key( *elem, SL );
	}
}

// Zeittest für dynamische SkipListe
template<typename TKey, typename TVal>
void 
checkSL_compl_dyn(	std::vector< std::pair<TKey, TVal> > & pairs,
					TKey * search_keys,
					TKey searches )
{
	SkipList< std::pair<TKey, TVal>, Dynamic > SL( begin( pairs ), end( pairs ) );
	volatile TKey vol;
	for( TKey i = 0; i < searches; ++i ){
		vol = key( searchElement( SL, search_keys[i] ) );
	}
}



// Zeittest für statische SkipListe

template<typename TKey, typename TVal>
void 
checkSL_compl_stat(	std::vector< std::pair<TKey, TVal> > & pairs,
					TKey * search_keys,
					TKey searches )
{
	SkipList< std::pair<TKey, TVal>, Static > SL( begin( pairs ), end( pairs ) );
	volatile TKey vol;
	for( TKey i = 0; i < searches; ++i ){
		SkipBaseElement< std::pair<TKey, TVal>, Static > * elem = _pointsTo( searchElement( SL, search_keys[i] ) );
		vol = key( *elem, SL );
		//std::cout<< search_keys[sk] << " " << key( *(elem-1) ) << " " << key(* elem ) << " " << key( *(elem+1) ) << std::endl;
	}
}

// Zeittest für Map

template<typename TKey, typename TVal>
void 
checkMultiMap(	std::vector<std::pair<TKey, TVal> > & pairs_bm, 
			TKey searches_bm,
			TKey * search_keys_bm)
{
	int sk_bm = 0;
	volatile TKey vol;
	std::multimap< TKey, TVal > map_bm(pairs_bm.begin(), pairs_bm.end());
	typename std::multimap< TKey, TVal >::iterator result;
	for( TKey i = 0; i < searches_bm; i++){
		result = map_bm.find( search_keys_bm[i] );
		if( result != map_bm.end() )
			vol = (*result).first;
		sk_bm++;
	}
}

// Zeittest für Map

template<typename TKey, typename TVal>
void 
checkMap(	std::vector<std::pair<TKey, TVal> > & pairs_bm, 
			TKey searches_bm,
			TKey * search_keys_bm)
{
	int sk_bm = 0;
	volatile TKey vol;
	std::map< TKey, TVal > map_bm(pairs_bm.begin(), pairs_bm.end());
	typename std::map< TKey, TVal >::iterator result;
	for( TKey i = 0; i < searches_bm; i++){
		result = map_bm.find( search_keys_bm[i] );
		if( result != map_bm.end() )
			vol = (*result).first;
		sk_bm++;
	}
}

// Zeittest lineare Suche im Array
template<typename TKey, typename TVal>
void checkLinear(	TKey * keys,
					TKey num_keys,
					TKey searches_bm,
					TKey * search_keys_bm)
{
	for( TKey i = 0; i < searches_bm; ++i){
		volatile TKey result;
		for( TKey j = 0; j < num_keys; ++j ){
			if( keys[j] == 	search_keys_bm[i] ){
				result = keys[j];
				break;
			}
		}
	}
}



// Zeittest für Vector

template<typename TKey, typename TVal>
void checkVector(	std::vector<std::pair<TKey, TVal> > & pairs_bm, 
			TKey searches_bm,
			TKey * search_keys_bm)
{
	int sk_bm = 0;
	volatile TKey vol;
	std::vector<TKey> vector_bm;
	vector_bm.reserve( pairs_bm.size());
	for( TKey i = 0; i < pairs_bm.size(); i++ ){
		vector_bm.push_back( pairs_bm[i].first );
	}
	typename std::vector< TKey >::iterator result;
	for( TKey i = 0; i < searches_bm; i++){
		result = std::find( vector_bm.begin(), vector_bm.end(), search_keys_bm[sk_bm] );
		if( result != vector_bm.end() )
		{
			vol = *result;
			break;
		}
		sk_bm++;
	}
}




void testDynamicSkipList()
{	
	int start_elems = 10;
	int ips = 5;
	int i_steps = 10;
	int search_tests = 100;

	std::vector< std::pair<int, char> > elems;
	String< std::pair<int, char> > s_elems;
	elems.reserve( start_elems + i_steps * search_tests );
	reserve( s_elems, start_elems + i_steps * search_tests );

	for( int i = 0; i < start_elems + ips * i_steps; ++i )
	{
		elems.push_back( std::make_pair<int, char>( rand() % ( i_steps * start_elems )  + 1, ('a') + i%26 ) );
	}
	
	std::multimap< int, char > map( elems.begin(), elems.begin() + start_elems );

	//SkipList< std::pair<int, char>, Dynamic, Default, Deferred > DSL( begin( elems ), begin( elems ) + start_elems );
	SkipList< std::pair<int, char>, Dynamic, Default, Complete > SL( begin( s_elems ), begin( s_elems ) + start_elems );

	std::vector< std::pair<int, char> >::iterator it = elems.begin() + start_elems;

	for( int i = 0; i < i_steps; ++i )
	{
		for( int j = 0; j < ips; ++j )
		{
			std::pair<int, char> & obj = *it;
			map.insert( obj );
			insert( SL, obj );
			//dump( SL );
			//insert( DSL, obj );
			++it;
		}
		
		int mm_result = 0;
		int sl_res = 0;
		int dsl_res = 0;
		std::multimap< int, char >::iterator result;
		for( int j = 0; j < search_tests; ++j )
		{
			int search_key= rand() % ( i_steps * start_elems ) + 1;
			result = map.lower_bound( search_key);
			//mm_result = ( *result ).first;
			//dump( SL );
			sl_res = key( searchElement( SL, search_key) );
			////dsl_res = key( searchElement( DSL, search_key) );
			//if( sl_res != dsl_res )
			//	std::cout << std::endl;
		}
	}
}

*/

//
//// Test dynamischer SkipListe mit Einfügen und Löschen
//void testElementsDynamicDeferred()
//{
//	std::vector< std::pair<int, char> > elems;
//	elems.reserve( 20 );
//
//	for( int i = 0; i < 20; i++ )
//	{
//		elems.push_back( std::make_pair<int, char>( rand() % 100 + 1, ('a') + i%26 ) );
//	}
//
//		// Multiples Vorkommen von Schlüsseln ist zugelassen
//	elems[1].first = 25;
//	elems[1].second = '-';
//
//	elems[3].first = 25;
//	elems[3].second = '#';
//
//	elems[5].first = 25;
//	elems[5].second = '.';
//
//	elems[10].first = 25;
//	elems[10].second = '*';
//
//		// Konstruktion der SkipListe mit einer Range
//	SkipList< std::pair<int, char>, Dynamic, Default, Deferred > S( begin( elems ), begin( elems ) + 15 );
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Beginn: " << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		char c;
//		std::cin >> c;
//	#endif
//	
//		// In buffer werden die Ergebnisse gespeichert
//	Iterator< SkipList< std::pair<int, char>, Dynamic, Default, Deferred > >::Type buffer; 
//
//		// Einfügen in unsortierte Liste -> Elemente werden am Ende eingefügt
//	for( int i = 10; i < 12; i++ )
//	{	
//		std::cout << key( elems[i] ) << std::endl;
//		insert( S, elems[i] );
//	}
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach erstem einfügen(in unsortierte Liste): " << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//
//		// Suche nach dem Element mit dem Schlüssel 50
//	buffer = searchElement( S, 50 );
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach Suche nach 50: " << std::endl;
//	std::cout << " Ergebnis: " << std::endl;
//	std::cout << key( buffer ) << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//	
//		// Suche nach 25
//	buffer = searchElement( S, 25 );
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach Suche nach 25: " << std::endl;
//	std::cout << " Ergebnis: " << std::endl;
//	std::cout << key( buffer ) << std::endl;
//	dump( S );
//
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach Einfügen in (teil-)sortierte Liste: " << std::endl;
//	std::cout << " Einfügen von: " << std::endl;
//	/*for( int i = 12; i < 15; i++ ){
//		std::cout << elems[i].first << std::endl;
//		insert( S, elems[i] );
//	}*/
//
//	std::cout << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//		
//	buffer = searchElement( S, 50 );
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach Suche nach 50: " << std::endl;
//	std::cout << " Ergebnis: " << std::endl;
//	std::cout << key( buffer ) << std::endl;
//	std::cout << "Löschen des erhaltenen Elementes:" << std::endl;
//	Key< std::pair< TKey, TVal > >::Type rand_key= key( buffer );
//	//erase( S, buffer );
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//	
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << "Löschen des Elementes, welches das Objekt enthält:" << key( elems[4] ) << std::endl;
//	//erase( S, elems[4] );
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach Suche nach 50: " << std::endl;
//	std::cout << " Ergebnis: " << std::endl;
//	buffer = searchElement( S, 50 );
//	std::cout << key( buffer ) << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//	
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach Suche nach 55: " << std::endl;
//	std::cout << " Ergebnis: " << std::endl;
//	buffer = searchElement( S, 55 );
//	std::cout << key( buffer ) << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//}
//
//
//// Test dynamischer SkipListe mit Einfügen und Löschen
//void testElementsDynamicComplete()
//{
//	std::vector< std::pair<int, char> > elems;
//	elems.reserve( 15 );
//
//	for( int i = 0; i < 15; i++ )
//	{
//		elems.push_back( std::make_pair<int, char>( rand() % 100 + 1, ('a') + i%26 ) );
//	}
//
//		// Multiples Vorkommen von Schlüsseln ist zugelassen
//	elems[1].first = 25;
//	elems[1].second = '-';
//
//	elems[5].first = 25;
//	elems[5].second = '#';
//
//	elems[10].first = 25;
//	elems[10].second = '.';
//
//	elems[14].first = 25;
//	elems[14].second = '*';
//
//		// Konstruktion der SkipListe mit einer Range
//	SkipList< std::pair<int, char>, Dynamic > S( begin( elems ), begin( elems ) + 10 );
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Beginn: " << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		char c;
//		std::cin >> c;
//	#endif
//	
//		// In buffer werden die Ergebnisse gespeichert
//	Iterator< SkipList<std::pair<int, char>, Dynamic > >::Type buffer; 
//
//		// Einfügen in unsortierte Liste -> Elemente werden am Ende eingefügt
//	for( int i = 10; i < 15; i++ ){
//		std::cout << "Inserting " <<  elems[i].first << std::endl;
//		insert( S, elems[i] );
//		dump( S );
//	}
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach erstem einfügen: " << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//
//		// Suche nach dem Element mit dem Schlüssel 50
//	buffer = searchElement( S, 50 );
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach Suche nach 50: " << std::endl;
//	std::cout << " Ergebnis: " << std::endl;
//	std::cout << key( buffer ) << std::endl;
//	dump( S );
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//	
//		// Suche nach 25
//	buffer = searchElement( S, 25 );
//	std::cout << "________________________________________________________________________________" << std::endl;
//	std::cout << " Nach Suche nach 25: " << std::endl;
//	std::cout << " Ergebnis: " << std::endl;
//	std::cout << key( buffer ) << std::endl;
//	dump( S );
//
//		
//	std::cout << "________________________________________________________________________________" << std::endl;
//	for( int i = 0; i < 5; ++i ){
//		std::pair<int, char> * elem = &elems[ rand() % 15 ];
//		std::cout << "Löschen des Elementes, welches das Objekt enthält:" << key( *elem ) << std::endl;
//		erase( S, *elem );
//		dump( S );
//	}
//	#ifdef SEQAN_SL_DEMO
//		std::cin >> c;
//	#endif
//
//	
//}

/*
Siehe "skip_list_type.h"

template <typename T1, typename T2, typename TCompression>
inline T1 &
key(Pair<T1, T2, TCompression> const &pair) 
{
	return pair.i1;
}
template <typename T1, typename T2, typename TCompression>
inline T2 &
value(Pair<T1, T2, TCompression> const &pair) 
{
	return pair.i2;
}
*/
template <typename TKey, typename TVal, typename TVal2> inline
TVal & setValue( std::pair< TKey, TVal > & p, TVal2 const & val_ )
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

