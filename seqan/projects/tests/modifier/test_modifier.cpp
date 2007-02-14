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

	// custom functor (Caesar chiffre)
    template <typename InType, typename Result = InType>
    struct CaesarChiffre : public unary_function<InType,Result> 
	{
		InType delta;

		CaesarChiffre() {}
		CaesarChiffre(InType _delta) {
			if (_delta < 0)
				delta = ('z' - 'a' + 1) - (-_delta) % ('z' - 'a' + 1);
			else
				delta = _delta;
		}

        inline Result operator()(InType x) const {
			if (('a' <= x) && (x <= 'z')) return (x - 'a' + delta) % ('z' - 'a' + 1) + 'a';
			if (('A' <= x) && (x <= 'Z')) return (x - 'A' + delta) % ('Z' - 'A' + 1) + 'A';
			return x; 
		}
    };



void testIterators()
{
		String<char> origin = "Vjku ku qwt qtkikpcn uvtkpi";

	//____________________________________________________________________________
	// Test1 - no modification (default)
	{

		typedef ModifiedIterator< Iterator<String<char> >::Type > TModIterDefault;

		TModIterDefault it, itEnd(end(origin));

		cout << "*** Iterator Test: Caesar chiffre ***" << endl;
		cout << "chiffre:  ";
		
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

		TEncode encode(-2);
		TModIterCaesar it(encode), itEnd(end(origin));

		cout << "original: ";

		it = begin(origin);
		while (it != itEnd) {
			cout << *it;
			it = it + 1;
		}
		cout << endl << endl;

	}
}

void testStrings()
{
		String<char> origin = "This is our original string";

	//____________________________________________________________________________
	// Test1 - no modification (default)

		typedef ModifiedString< String<char> > TModStringDefault;

		TModStringDefault nomod(origin);

		cout << "*** Test1/2: Caesar chiffre ***" << endl;
		cout << "origin:  " << nomod << endl;
		
	//____________________________________________________________________________
	// Test2 - Caesar chiffre

		typedef CaesarChiffre<char> TEncode;
		typedef ModifiedString< String<char>, ModView<TEncode> > TModStringCaesar;

		TEncode encode(2);
		TModStringCaesar chiffre(origin);
		assignModViewFunctor(chiffre, encode);

		cout << "chiffre: " << chiffre << endl << endl;

	//____________________________________________________________________________
	// Test3 - upcase/downcase

		typedef ModifiedString< String<char>, ModView< FunctorUpcase<char> > > TModStringUp;
		typedef ModifiedString< String<char>, ModView< FunctorDowncase<char> > > TModStringDown;

		TModStringUp	up(origin);
		TModStringDown	down(origin);

		cout << "*** Test3: upcase/downcase ***" << endl;
		cout << "upcase:   " << up << endl;
		cout << "downcase: " << down << endl << endl;

	//____________________________________________________________________________
	// Test4 - alphabet conversion

		String<char> originDNA = "acgtnACGTN";
		typedef ModifiedString< String<char>, ModView< FunctorConvert<char, Dna5> > > TModStringDNA;

		TModStringDNA	dna(originDNA);

		cout << "*** Test4: alphabet conversion ***" << endl;
		cout << "origin: " << originDNA << endl;
		cout << "Dna5:   " << dna << endl << endl;

	//____________________________________________________________________________
	// Test5 - nested modifiers

		typedef CaesarChiffre<char> TEncode;
		typedef ModifiedString< 
					ModifiedString< 
						ModifiedString< 
							String<char>, 
							ModView<TEncode> 
						>, 
						ModView<TEncode> 
					>, 
					ModView<TEncode>
				> TModStringNested;

		TModStringNested nested(origin);

		TEncode enc2(2), enc3(3), enc_5(-5);

		assignModViewFunctor(nested, enc2);
		assignModViewFunctor(host(nested), enc3);
		assignModViewFunctor(host(host(nested)), enc_5);

		cout << "*** Test5: nested modifiers ***" << endl;
		cout << "nested: " << nested << endl << endl;

}

int main()
{
	SEQAN_TREPORT("TEST BEGIN")

		testStrings();
		testIterators();

	SEQAN_TREPORT("TEST END")
		return 0;
}
