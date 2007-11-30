#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;

int main ()
{
	String<Dna> myString = "attacgg";

	typedef ModifiedString<String<Dna>, ModComplementDna>   TMyComplement;
	typedef ModifiedString<TMyComplement, ModReverse>       TMyReverseComplement;

	TMyReverseComplement myReverseComplement(myString);

	cout << myString << endl;
	cout << myReverseComplement << endl;

	infix(myString, 1, 1) = "cgt";

	cout << myString << endl;
	cout << myReverseComplement << endl;

	cout << DnaStringReverseComplement(myString) << endl;

	return 0;
}
