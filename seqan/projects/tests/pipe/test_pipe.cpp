#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>

#include <seqan/file.h>
#include <seqan/pipe.h>
#include <seqan/system.h>

#include "test_pipe.h"

using namespace std;
using namespace seqan;


// Maximum size of data to use for the external tests.
const size_t MAX_SIZE = 1u << 20;


// TODO(holtgrew): The following test* functions should actually be defined with SEQAN_DEFINE_TEST().


template <typename TStringSpec>
void testExternalString(unsigned maxSize = 16*1024*1024) 
{
	typedef String<unsigned, TStringSpec> TExtString;
    typedef typename Iterator<TExtString const, Standard>::Type TIter;

    SimpleBuffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    TExtString extString;
    for(unsigned i = 1; i <= maxSize; i = i << 1) 
	{
        // ::std::cout << i << " "; ::std::cout.flush();
        resize(buf, i);
        randomize(buf);

        Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
        extString << src;

        TIter I = begin(extString);
        for(unsigned *cur = buf.begin; cur != buf.end; ++cur) {
            if (*cur != *I) {
                SEQAN_ASSERT_FAIL("testExternalString failed at position %u", cur - buf.begin);
                // not reached
            }
            ++I;
        }
    }
    freePage(buf, buf);
}



void testPool(unsigned maxSize = 16*1024*1024) {
    SimpleBuffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    Pool<unsigned,PoolSpec<> > pool;
    for(unsigned i = 1; i <= maxSize; i = i << 1) {
        /*
        ::std::cout << i;
        ::std::cout.flush();
        */

        resize(buf, i);
        randomize(buf);

        Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
        pool << src;

        /*
        if (pool.memBuffer.begin)
            ::std::cout << "* ";
        else 
            ::std::cout << " ";
        ::std::cout.flush();
        */

        beginRead(pool);
        for(unsigned *cur = buf.begin; cur != buf.end; cur++) {
            if (*cur != *pool) {
                endRead(pool);
                SEQAN_ASSERT_FAIL("testPool failed at position %u", cur - buf.begin);
                // not reached
            }
            ++pool;
        }
        endRead(pool);
    }
    freePage(buf, buf);
}



void testMapper(unsigned maxSize = 16*1024*1024) {
    SimpleBuffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    Pool<unsigned,MapperSpec<MapperConfig<IdentityMap<unsigned> > > > mapper;
    for(unsigned i = 1; i <= maxSize; i = i << 1) {
        /*
        ::std::cout << i;
        ::std::cout.flush();
        */

        resize(buf, i);
        permute(buf);

        Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
        mapper << src;

        /*(
        if (mapper.memBuffer.begin)
            ::std::cout << "* ";
        else 
            ::std::cout << " ";
        ::std::cout.flush();
        */

        beginRead(mapper);
        for(unsigned j = 0; j < i; ++j) {
            if (*mapper != j) {
                freePage(buf, buf);
                SEQAN_ASSERT_FAIL("testMapper failed at position %u", j);
                // not reached
            }	
            ++mapper;
        }
        endRead(mapper);
    }
    freePage(buf, buf);
}



void testPartiallyFilledMapper(unsigned maxSize = 16*1024*1024) {
    SimpleBuffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    Pool<unsigned,MapperSpec<MapperConfig<IdentityMap<unsigned> > > > mapper;
    for(unsigned i = 1; i <= maxSize; i = i << 1) {
        /*
        ::std::cout << i;
        ::std::cout.flush();
        */

        resize(buf, i);
        permute(buf);

        // partially fill the mapper 
        mapper.undefinedValue = i;	// select i as an undefined value (all defined values are less than i)
        resize(mapper, i);
        resize(buf, i - i/3);
        Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
        beginWrite(mapper) && append(mapper, src) && endWrite(mapper);

        /*
        if (mapper.memBuffer.begin)
            ::std::cout << "* ";
        else 
            ::std::cout << " ";
        ::std::cout.flush();
        */

        unsigned undefCounter = 0, missCounter = 0;
        beginRead(mapper);
        for(unsigned j = 0; j < i; ++j) {
            if (*mapper == i) 
                ++undefCounter;
            else
                if (*mapper != j) {
                    ++missCounter;
                    if (!mapper.memBuffer.begin) { // external mapping -> no misses allowed
                        freePage(buf, buf);
                        SEQAN_ASSERT_FAIL("testPartiallyFilledMapper failed at position %u [ = %u ]", j, *mapper);
                        // not reached
                    }
                }
            ++mapper;
        }
        endRead(mapper);
        if (mapper.memBuffer.begin) {
            if (undefCounter + missCounter > i/3) {
                SEQAN_ASSERT_FAIL("testPartiallyFilledMapper failed [only %u of %u undefind", undefCounter + missCounter, i / 3);
                // not reached
            }
        } else
            if (undefCounter != i/3) {
                SEQAN_ASSERT_FAIL("testPartiallyFilledMapper failed [only %u of %u undefined]", undefCounter + missCounter, i / 3);
                // not reached
            }
    }
    freePage(buf, buf);
}



void testSorter(unsigned maxSize = 16*1024*1024) {
    SimpleBuffer<unsigned> buf;
    allocPage(buf, maxSize, buf);

    Pool<unsigned,SorterSpec<SorterConfig<SimpleCompare<unsigned> > > > sorter;
    for(unsigned i = 1; i <= maxSize; i = i << 1) {
        /*
        ::std::cout << i;
        ::std::cout.flush();
        */

        resize(buf, i);
        permute(buf);

        Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
        sorter << src;

        /*
        if (sorter.memBuffer.begin)
            ::std::cout << "* ";
        else 
            ::std::cout << " ";
        ::std::cout.flush();
        */

        beginRead(sorter);
        unsigned j = *sorter, pos = 0;
        while (!eof(sorter)) {
            if (*sorter < j || /* *sorter < 0 || */ *sorter >= i) {
                freePage(buf, buf);
                SEQAN_ASSERT_FAIL("testSorter failed at position %u", pos);
                // not reached
            }
            j = *sorter;
            ++sorter; ++pos;
        }
        endRead(sorter);
    }
    freePage(buf, buf);
}


SEQAN_DEFINE_TEST(test_pipe_test_external_string) {
    testExternalString<MMap<> >(MAX_SIZE);
    testExternalString<External<> >(MAX_SIZE);
}


SEQAN_DEFINE_TEST(test_pipe_test_simple_pool) {
    testPool(MAX_SIZE);
}


SEQAN_DEFINE_TEST(test_pipe_test_mapper) {
    testMapper(MAX_SIZE);
}


SEQAN_DEFINE_TEST(test_pipe_test_mapper_partially_filled) {
    testPartiallyFilledMapper(MAX_SIZE);
}


SEQAN_DEFINE_TEST(test_pipe_test_sorter) {
    testSorter(MAX_SIZE);
}


SEQAN_BEGIN_TESTSUITE(test_pipe) {
    SEQAN_CALL_TEST(test_pipe_test_external_string);
    SEQAN_CALL_TEST(test_pipe_test_simple_pool);
    SEQAN_CALL_TEST(test_pipe_test_mapper);
    SEQAN_CALL_TEST(test_pipe_test_mapper_partially_filled);
    SEQAN_CALL_TEST(test_pipe_test_sorter);

	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_caster.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_counter.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_echoer.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_edit_environment.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_filter.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_generated_forwards.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_iterator.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_joiner.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_namer.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_sampler.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_shifter.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_source.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pipe_tupler.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pool_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pool_mapper.h");
	SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/pipe/pool_sorter.h");
}
SEQAN_END_TESTSUITE

