// FRAGMENT(definitions)
#define SEQAN_PROFILE // enable time measurements

#include <iostream>
#include <seqan/basic.h>


using namespace seqan;

int main()
{
	Allocator<MultiPool< > > mpa;
	Allocator<SimpleAlloc< > > sa;
	unsigned runs = 100000;
	char *buf;
// FRAGMENT(time-measurements)
	// store blocksizes in an array
	unsigned bs[3];
	bs[0] = 10;
	bs[1] = 100;
	bs[2] = 1000;
	SEQAN_PROTIMESTART(timeAlloc);
	
	// loop through the different block sizes
	for (int i=0; i<3; ++i) {
		for (int j=0; j<runs; ++j) {
			allocate(mpa,buf,bs[i],TagAllocateTemp());
		}
		clear(mpa);
		std::cout << "Allocating and clearing " << runs << " times blocks of size " 
					<< bs[i] << " with MultiPool Allocator took " << SEQAN_PROTIMEDIFF(timeAlloc) << std::endl;	
		// the PROTIMEDIFF macro sets the internal timer clock to zero
		for (int j=0; j<runs; ++j) {
			allocate(sa,buf,bs[i],TagAllocateTemp());
		}
		clear(sa);
		std::cout << "Allocating and clearing " << runs << " times blocks of size " 
		<< bs[i] << " with Standard Allocator took " << SEQAN_PROTIMEDIFF(timeAlloc) << std::endl;	
	}

	return 0;
}

