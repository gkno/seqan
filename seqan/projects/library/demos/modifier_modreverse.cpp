#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;


int main ()
{
	String<char> myString = "A man, a plan, a canal-Panama";
	ModifiedString< String<char>, ModReverse > myModifier(myString);

	cout << myString << endl;
	cout << myModifier << endl;

	infix(myString, 9, 9) = "master ";

	cout << myString << endl;
	cout << myModifier << endl;

	return 0;
}
