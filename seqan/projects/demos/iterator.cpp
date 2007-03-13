#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

int main()
{
	seqan::String<char> str = "admn";
///The metafunction @Metafunction.Iterator@ returns the iterator type for a given container type.
	seqan::Iterator<seqan::String<char> >::Type it = begin(str);

	std::cout << *it;                //output: 'a'
///The following lines show a loop through $str$ in the standard library style:
	while (it != end(str))           //output: "admn"
	{
		std::cout << *it;
		++it;
	}
	std::cout << std::endl;
///Seqan offers an alternative style for accessing iterators that avoids operators.
///Note that the functions @Function.goBegin@ and @Function.atEnd@ do net get $str$ as arguments,
/// because $it$ is a @Concept.Rooted Iterator.rooted iterator@.
///The following loop increments each character in $str$:
	for (goBegin(it); !atEnd(it); goNext(it)) 
	{
		++value(it);
	}
///This is a reverse loop through $str$.
///Note that @Function.goPrevious@ is called before the value of $it$ is accessed,
/// because the end position of a container is the position behind the last item in the container:
	goEnd(it);

	while (!atBegin(it))              //output: "oneb"
	{
		goPrevious(it);
		std::cout << getValue(it);
	}
	std::cout << std::endl;
///Another (write only) way to access the value of an iterator is @Function.assignValue@:
	assignValue(begin(str), 'X');

	std::cout << str << std::endl;        //output: "Xeno"
	

	return 0;
}
