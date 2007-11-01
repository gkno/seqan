// Projekt, mit dem die Demos getestet werden koennen

#define main runAllocator
#include "../../demos/allocator.cpp"
#undef main

#define main runAlphabet
#include "../../demos/alphabet.cpp"
#undef main

/*
//!kompiliert nicht!
#define main runIterator
#include "../../demos/iterator.cpp"
#undef main
*/

#define main runRootedIterator
#include "../../demos/rooted_iterator.cpp"
#undef main

#define main runString
#include "../../demos/string_1.cpp"
#undef main

#define main runSufarray
#include "../../demos/sufarray.cpp"
#undef main

#define main runIndex
#include "../../demos/index.cpp"
#undef main

#define main runGraph
#include "../../demos/graph.cpp"
#undef main

#define main runFind
#include "../../demos/find.cpp"
#undef main

#define main runFileFormat
#include "../../demos/file_format.cpp"
#undef main

int main() 
{
	runAllocator();
	runAlphabet();
//	runIterator();
	runRootedIterator();
	runString();
	runSufarray();
	runIndex();
	runGraph();
	runFind();
	runFileFormat();
}
