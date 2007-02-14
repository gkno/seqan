#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

int main()
{
	seqan::String<char> str = "abcdefg";
///The @Metafunction.Iterator.Iterator metafunction@ returns a @Concept.Rooted Iterator.rooted iterator@ by default. 
///	We can also specify the iterator kind explicitly by passing the @Tag.Iterator Spec.iterator spec Rooted@ as second argument.
	seqan::Iterator<seqan::String<char>, seqan::Rooted>::Type it = begin(str);
///The same iterator spec can be used as a tag for functions that return an iterator, e.g. @Function.begin@ or @Function.end@
	it = end(str, seqan::Rooted());
///A rooted iterator "knows" its container, so it supports the function @Function.container@:
	std::cout << container(it);          //output: "abcdefg"
///Moreover, it is possible to apply functions like @Function.goBegin@ or @Function.position@ to rooted iterators without specify the container.
	std::cout << position(it);           //output: 7

	return 0;
}
