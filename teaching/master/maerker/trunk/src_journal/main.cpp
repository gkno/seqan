#define SEQAN_PROFILE

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "STNode.h"
#include "Tree.h"
#include <seqan.h>
#include "JournalString.h"
#include "JournalIterator.h"
#include "Utility.h"

using namespace std;
using namespace seqan;

class Base{
public:

	Base( int f ) : m( f ) { };
	virtual ~Base() {
		cout << "Base Destroy! " << m << endl;
	};

	int m;
};

class Derived : public Base {
public:
	Derived( int f ) : Base( f ) { };
	virtual ~Derived() {
		cout << "Derived Destroy!" << endl;
	}
};


int main(){

	String< char > wombat = "Wombat";
	String< char, Journal< Alloc<>, Basic > > grombat( wombat );

	std::cout << grombat << std::endl;
	erase( grombat, 1, 4 );
	std::cout << "erasing 1..3:" << std::endl;
	std::cout << grombat << std::endl;
	std::cout << "flattening . . ." << std::endl;
	flatten( grombat );
	std::cout << "journal state:" << std::endl;
	std::cout << grombat << std::endl;
	std::cout << "underlying string:" << std::endl;
	std::cout << wombat << std::endl;

	seqan::util::generate_random_string( 36, wombat, "abcdexx" );
	String< char, Journal< Alloc<>, Basic > > schlombat( wombat );

	std::cout << schlombat << std::endl;
	std::cout << "replacing 'x' with '_':" << std::endl;
	String< char > ins = "_";

	for( unsigned int q = 0; q < length( schlombat ); ++ q ){
		if( value( schlombat, q ) == 'x' ){
			erase( schlombat, q, q+1 );
			insert( schlombat, q, ins );
		}
	}

	std::cout << schlombat << std::endl;
	std::cout << "flattening. . ." << std::endl;
	flatten( schlombat );
	std::cout << "journal state:" << std::endl;
	std::cout << schlombat << std::endl;
	std::cout << "underlying string:" << std::endl;
	std::cout << wombat << std::endl;

	return 23;
	// BEGIN: Tests that are not run anymore at the moment
	cout << "Length\tDeletions\tJournal\tAlloc" << endl;
	SEQAN_PROTIMESTART(_alloc);
	SEQAN_PROTIMESTART(_journal);

	String< char > str_alloc;
	String< unsigned int > positions;
	for( unsigned int i = 100000; i <= 1600000; i *= 2 ){

		for( unsigned int k = 1; k <= 20; ++k ){
			for( unsigned int l = 0; l <= 10; ++l ){
				seqan::util::generate_random_string( i, str_alloc, "acgt" );
				String< char, Journal< Alloc<>, Basic > > str_journal( str_alloc );

				for( unsigned int k = 0; k < 1000; ++k ){
					appendValue( positions, seqan::util::randomNumber( seqan::length( str_alloc ) - 1001 ), Generous() );
				}

				cout << i << "\t" << 20*l << "\t";

				for( unsigned int j = 0; j < 20*l; ++j ){
					erase( str_journal, positions[j], positions[j] + 1 );
				}

				char tmp;
				SEQAN_PROTIMEUPDATE(_journal);
				for( unsigned int k = 0; k < 100; ++k ){
					for( unsigned int j = 0; j < length( positions ); ++j ){
						tmp = value( str_journal, positions[j] );
					}
				}
				cout << SEQAN_PROTIMEUPDATE(_journal) << "\t";

				SEQAN_PROTIMEUPDATE(_alloc);
				for( unsigned int k = 0; k < 100; ++k ){
					for( unsigned int j = 0; j < length( positions ); ++j ){
						tmp = value( str_alloc, positions[j] );
					}
				}
				cout << SEQAN_PROTIMEUPDATE(_alloc) << "\n";

				clear( positions );
				clear( str_alloc );
			}
		}
	}

	cout << endl;
	cout << " Fin " << endl;
}

