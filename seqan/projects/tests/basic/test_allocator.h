// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <iostream>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/basic.h"

#include <memory>
#include <vector>
#include <map>

using namespace std;
using namespace seqan;

//____________________________________________________________________________

struct TestAllocator
{
	mutable map<char *, size_t> data_allocated;
	mutable map<char *, size_t> data_deallocated;

	TestAllocator() {}
	~TestAllocator() 
	{
		map<char *, size_t>::iterator it = data_allocated.begin();
		while (it != data_allocated.end())
		{
			SEQAN_ASSERT_TRUE_MSG(data_deallocated.count(it->first), "memory block not deallocated");
            deallocate(int(), it->first, it->second);
			++it;
		}
	}
};

template <typename TValue, typename TSize, typename TUsage>
void allocate(TestAllocator & me, 
			  TValue * & data_, 
			  TSize count, 
			  Tag<TUsage> const)
{
	SEQAN_ASSERT(count)
	allocate(int(), data_, count);
	me.data_allocated[(char *) data_] = count;
}

template <typename TValue, typename TSize, typename TUsage>
void deallocate(TestAllocator & me, 
				TValue * data_, 
				TSize count, 
				Tag<TUsage> const)
{
	SEQAN_ASSERT_TRUE_MSG(me.data_allocated.count((char *) data_), "memory block was not allocated");
	SEQAN_ASSERT_TRUE_MSG(me.data_allocated[(char *) data_] == count, "memory block was allocated with different size");
	SEQAN_ASSERT_TRUE_MSG(!me.data_deallocated.count((char *) data_), "memory block already deallocated");

	me.data_deallocated[(char *) data_] = count;
}

int countAllocs(TestAllocator & me)
{
	return me.data_allocated.size();
}
int countDeallocs(TestAllocator & me)
{
	return me.data_deallocated.size();
} 
 


SEQAN_DEFINE_TEST(testSimpleAllocator) {
	int * dat1;
	int * dat2;

	Allocator<SimpleAlloc<TestAllocator> > allo1;
	allocate(allo1, dat1, 100);
	allocate(allo1, dat2, 105);
	deallocate(allo1, dat1, 100);
	allocate(allo1, dat2, 201);

	SEQAN_ASSERT_TRUE(countAllocs(parentAllocator(allo1)) == 3);
	SEQAN_ASSERT_TRUE(countDeallocs(parentAllocator(allo1)) == 1);

	clear(allo1);

	SEQAN_ASSERT_TRUE(countDeallocs(parentAllocator(allo1)) == 3);
}
//____________________________________________________________________________

SEQAN_DEFINE_TEST(testPoolAllocator) {
	int * dat1;
	int * dat2;

	typedef Allocator<SimpleAlloc<TestAllocator> > TParentAlloc;
	Allocator<SinglePool<20 * sizeof(int), TParentAlloc> > allo1;
	allocate(allo1, dat1, 20);
	allocate(allo1, dat2, 20);
	deallocate(allo1, dat1, 20);
	allocate(allo1, dat2, 20);

	SEQAN_ASSERT_TRUE(dat1 == dat2);

	SEQAN_ASSERT_TRUE(countAllocs(parentAllocator(parentAllocator(allo1))) == 1);
	SEQAN_ASSERT_TRUE(countDeallocs(parentAllocator(parentAllocator(allo1))) == 0);

	allocate(allo1, dat1, 100);
	deallocate(allo1, dat1, 100);

	SEQAN_ASSERT_TRUE(countAllocs(parentAllocator(parentAllocator(allo1))) == 2);
	SEQAN_ASSERT_TRUE(countDeallocs(parentAllocator(parentAllocator(allo1))) == 1);

	clear(allo1);

	SEQAN_ASSERT_TRUE(countDeallocs(parentAllocator(parentAllocator(allo1))) == 2);
}

//____________________________________________________________________________

SEQAN_DEFINE_TEST(testMultiPoolAllocator) {
	int * dat1;
	int * dat2;

	typedef Allocator<SimpleAlloc<TestAllocator> > TParentAlloc;
	Allocator<MultiPool<TParentAlloc> > allo1;
	allocate(allo1, dat1, 20);
	allocate(allo1, dat2, 20);
	deallocate(allo1, dat1, 20);
	allocate(allo1, dat2, 20);

	SEQAN_ASSERT_TRUE(dat1 == dat2);

	SEQAN_ASSERT_TRUE(countAllocs(parentAllocator(parentAllocator(allo1))) == 1);
	SEQAN_ASSERT_TRUE(countDeallocs(parentAllocator(parentAllocator(allo1))) == 0);

	allocate(allo1, dat1, 30);
	deallocate(allo1, dat1, 30);

	SEQAN_ASSERT_TRUE(countAllocs(parentAllocator(parentAllocator(allo1))) == 2);
	SEQAN_ASSERT_TRUE(countDeallocs(parentAllocator(parentAllocator(allo1))) == 0);

	clear(allo1);

	SEQAN_ASSERT_TRUE(countDeallocs(parentAllocator(parentAllocator(allo1))) == 2);
}

//____________________________________________________________________________
