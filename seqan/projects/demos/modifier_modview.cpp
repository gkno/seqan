#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace std;
using namespace seqan;

struct MyFunctor : public unary_function<char,char> 
{
	inline char operator()(char x) const 
	{
		if (('a' <= x) && (x <= 'z')) return (x + ('A' - 'a'));
		return x; 
	}
};


int main ()
{
	String<char> myString = "A man, a plan, a canal-Panama";
	ModifiedString< String<char>, ModView<MyFunctor> > myModifier(myString);

	cout << myString << endl;
	cout << myModifier << endl;

	infix(myString, 9, 9) = "master ";

	cout << myString << endl;
	cout << myModifier << endl;

	return 0;
}
