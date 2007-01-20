#include <seqan/basic.h>
using namespace seqan;

//An arbitrary class
struct MyClass
{
};

int main()
{
///The following code creates 100 instances of $MyClass$ on the heap.
///The allocator object used is a temporary $Default$.
	MyClass * my_class_arr;
	allocate(Default(), my_class_arr, 100);
	arrayConstruct(my_class_arr, my_class_arr + 100);
///Before the storage is deallocated, the $MyClass$ objects must be destroyed.
	arrayDestruct(my_class_arr, my_class_arr + 100);
	deallocate(Default(), my_class_arr, 100);
///We can use any kind of object as allocator, but dedicated allocators offer more advanced functionality, e.g. @Function.clear@.
	Allocator<SimpleAlloc< > > alloc1;
	allocate(alloc1, my_class_arr, 200);

	char * char_array;
	allocate(alloc1, char_array, 300);

	clear(alloc1); //deallocated all storage at once.
	
	return 0;
}
