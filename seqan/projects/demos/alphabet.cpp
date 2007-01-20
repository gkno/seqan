#include <iostream>
#include <seqan/basic.h>
using namespace seqan;

int main()
{
///The typical alphabet is convertible to $char$.
///Note that a conversion of a $char$ into another alphabet and back can change the value of the $char$.
	Dna a = 'a';
	std::cout << a; //output: 'A'

	Dna5 b = 'f'; //'f' is unknown character
	std::cout << b; //output: 'N'

///Many SeqAn alphabet classes can be converted into each other.
	b = a;
	std::cout << b; //output: 'A'

	Iupac c = b;
	std::cout << c; //output: 'A'

	return 0;
}
