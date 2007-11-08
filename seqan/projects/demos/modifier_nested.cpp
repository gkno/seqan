#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;

int main ()
{
	String<Dna> myString = "attacgg";

	typedef ModifiedString<DnaString, ModView< FunctorComplement<Dna> > > TMyComplement;
	typedef ModifiedString<TMyComplement, ModReverse> TMyReverseComplement;

	TMyComplement myComplement(myString);
	TMyReverseComplement myReverseComplement(myString);

	cout << myString << endl;
	cout << myReverseComplement << endl;

	infix(myString, 1, 1) = "cgt";

	cout << myString << endl;
	cout << myReverseComplement << endl;

	return 0;
}
