// FRAGMENT(swap-headers)
#include <iostream>
#include <seqan/basic.h>
#include <seqan/file.h>


using namespace seqan;
using namespace std;

// FRAGMENT(swap-declaration)
template <typename T> void swap(T& container, int i, int j)
{
// FRAGMENT(swap-metafunction)
	typename Value<T>::Type help = value(container,i);
// FRAGMENT(swap-work)
	value(container,i) = value(container,j);
	value(container,j) = help;
	return;
}


// FRAGMENT(swap-main)
int main()
{
	typedef String<Dna> TDnaString;
	TDnaString dna = "AATT";
	
	typedef String<int> TIntString;
	typedef Iterator<String<int>, Rooted >::Type TIntIterator;
	
	TIntString numbers;
    appendValue(numbers,1);
	appendValue(numbers,1);
	appendValue(numbers,3);
	appendValue(numbers,3); 
	
// FRAGMENT(swap-apply)
	swap(dna,1,4);
	cout << dna << endl;
	
    swap(numbers,1,2);
	for (TIntIterator it=begin(numbers); !atEnd(it); goNext(it)) {
		::std::cout << *it;
	}
	cout << endl;
	
	return 0;
}

