#include <iostream>
#include <fstream>
#include <functional>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

    template <typename InType, typename Result = InType>
    struct CaesarChiffre : public unary_function<InType,Result> 
	{
		InType delta;

		CaesarChiffre() {}
		CaesarChiffre(InType _delta): delta(_delta) {}

        Result operator()(InType x) const {
			if (('a' <= x) && (x <= 'z')) return (x - 'a' + delta) % ('z' - 'a' + 1) + 'a';
			if (('A' <= x) && (x <= 'Z')) return (x - 'A' + delta) % ('Z' - 'A' + 1) + 'A';
			return x; 
		}
    };



void testIterators()
{
		String<char> origin = "This is our original string";

	//____________________________________________________________________________
	// Test1 - no modification (default)
	{

		typedef ModifiedIterator< Iterator<String<char> >::Type > TModIterDefault;

		TModIterDefault it, itEnd(end(origin));

		it = begin(origin);
		while (it != itEnd) {
			cout << *it;
			++it;
		}
		cout << endl;

	}
	//____________________________________________________________________________
	// Test2 - Caesar chiffre
	{

		typedef CaesarChiffre<char> TEncode;
		typedef ModifiedIterator< Iterator<String<char> >::Type, ModView<TEncode> > TModIterCaesar;

		TEncode encode(2);
		TModIterCaesar it(encode), itEnd(end(origin));

		it = begin(origin);
		while (it != itEnd) {
			cout << *it;
			it = it + 1;
		}
		cout << endl;

	}
}

void testStrings()
{
		String<char> origin = "Vjku ku qwt qtkikpcn uvtkpi";

	//____________________________________________________________________________
	// Test1 - no modification (default)
	{

		typedef ModifiedString< String<char> > TModStringDefault;

		TModStringDefault modString(origin);

		cout << modString << endl;

	}
	//____________________________________________________________________________
	// Test2 - Caesar chiffre
	{

		typedef CaesarChiffre<char> TDecode;
		typedef ModifiedString< String<char>, ModView<TDecode> > TModStringCaesar;

		TDecode decode(-2);
		TModStringCaesar modString(origin);
		assignModViewFunctor(modString, decode);

		cout << modString << endl;

	}
}

int main()
{
	SEQAN_TREPORT("TEST BEGIN")

		testIterators();
		testStrings();

	SEQAN_TREPORT("TEST END")
		return 0;
}
