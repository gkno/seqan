/*
 *  basic_metaprogramming.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_BASIC_METAPROGRAMMING_H
#define SEQAN_BASIC_METAPROGRAMMING_H

namespace SEQAN_NAMESPACE_MAIN
{

	using ::std::memset;

//////////////////////////////////////////////////////////////////////////////
// generic metaprogramming

///Empty Data Class.
struct Nothing {};
struct Yes { enum { VALUE = true }; };
struct No { enum { VALUE = false }; };

//! \brief IF template metaprogramming statement

//! If \c Flag is \c true then \c IF<>::result is of type Type1
//! otherwise of \c IF<>::result is of type Type2
template <bool Flag,class Type1, class Type2>
struct IF
{
  typedef Type1 Type;
};

template <class Type1, class Type2>
struct IF<false,Type1,Type2>
{
  typedef Type2 Type;
};


// compare two types
template <class Type1, class Type2>
struct TYPECMP
{
    enum { VALUE = false };
};

template <class Type1>
struct TYPECMP<Type1, Type1>
{
    enum { VALUE = true }; 
};

const int DEFAULT = ~(~0u >> 1); // initialize with the smallest int

struct NilCase {};
  
template <int tag_,class Type_,class Next_ = NilCase>
struct CASE
{
  enum { tag = tag_ };
  typedef Type_ Type;
  typedef Next_ Next;
};

template <int tag,class Case>
class SWITCH
{
  typedef typename Case::Next NextCase;
  enum
  {
    caseTag = Case::tag,
    found   = (caseTag == tag || caseTag == DEFAULT)
  };
public:
  typedef typename IF<found,
          typename Case::Type,
          typename SWITCH<tag,NextCase>::Type
          >::Type Type;
  
};

template <int tag>
class SWITCH<tag,NilCase>
{
public:
  typedef NilCase Type;
};


// example of a loop Worker class
struct WorkerNothing
{
    template <typename Arg>
    static inline void body(Arg &arg, int I) {}
};

template <typename Worker, int I>
class LOOP {
public:
    template <typename Arg>
    static inline void run(Arg &arg) {
        LOOP<Worker, I - 1>::run(arg);
        Worker::body(arg, I);
    }
};

template <typename Worker>
class LOOP<Worker, 0> {
public:
    // end of loop
    template <typename Arg>
    static inline void run(Arg &arg) {}
};

template <typename Worker, int I>
class LOOP_REVERSE {
public:
    template <typename Arg>
    static inline void run(Arg &arg) {
        Worker::body(arg, I);
        LOOP_REVERSE<Worker, I - 1>::run(arg);
    }
};

template <typename Worker>
class LOOP_REVERSE<Worker, 0> {
public:
    // end of loop
    template <typename Arg>
    static inline void run(Arg &arg) {}
};


template < __int64 numerus >
struct Log2 {
    enum { VALUE = Log2<(numerus + 1) / 2>::VALUE + 1 };		// ceil(log_2(n))
};

template < __int64 numerus >
struct Log2Floor {
    enum { VALUE = Log2Floor<numerus / 2>::VALUE + 1 };		// floor(log_2(n))
};

template <> struct Log2<1> { enum { VALUE = 0 }; };
template <> struct Log2<0> { enum { VALUE = 0 }; };
template <> struct Log2Floor<1> { enum { VALUE = 0 }; };
template <> struct Log2Floor<0> { enum { VALUE = 0 }; };


// memset with metaprogramming

template <unsigned SIZE, bool direct>
struct MemsetWorker {
	finline static void run(unsigned char* ptr, unsigned char c) { ::std::memset(ptr, c, SIZE); }
};

template <unsigned  SIZE>
struct MemsetWorker<SIZE, true> {
    finline static void run(unsigned char* ptr, unsigned char c) {
		*((unsigned*)ptr) = (unsigned)c << 24 + (unsigned)c << 16 + (unsigned)c << 8 + (unsigned)c;
		MemsetWorker<SIZE - 4, true>::run(ptr + 4, c);
	}
};

template <>
struct MemsetWorker<0, true> {
	finline static void run(unsigned char* ptr, unsigned char c) {}
};

template <>
struct MemsetWorker<1, true> {
	finline static void run(unsigned char* ptr, unsigned char c) { *ptr = c; }
};

template <>
struct MemsetWorker<2, true> {
	finline static void run(unsigned char* ptr, unsigned char c) { *(unsigned short *)ptr = (unsigned short)c << 8 + (unsigned short)c; }
};

template <>
struct MemsetWorker<3, true> {
	finline static void run(unsigned char* ptr, unsigned char c) {
		MemsetWorker<2, true>::run(ptr, c);
		MemsetWorker<1, true>::run(ptr + 2, c);
	}
};

template <unsigned SIZE>
finline void memset(void* ptr, unsigned char c) {
	MemsetWorker<SIZE, SIZE <= 32>::run((unsigned char*)ptr, c);
}


// memset with metaprogramming and constant fill value

template <unsigned SIZE, bool direct, unsigned char c>
struct MemsetConstValueWorker {
	finline static void run(unsigned char* ptr) { ::std::memset(ptr, c, SIZE); }
};

template <unsigned  SIZE, unsigned char c>
struct MemsetConstValueWorker<SIZE, true, c> {
    finline static void run(unsigned char* ptr) {
		*((unsigned*)ptr) = (unsigned)c << 24 + (unsigned)c << 16 + (unsigned)c << 8 + (unsigned)c;
		MemsetConstValueWorker<SIZE - 4, true, c>::run(ptr + 4);
	}
};

template <unsigned char c>
struct MemsetConstValueWorker<0, true, c> {
	finline static void run(unsigned char* ptr) {}
};

template <unsigned char c>
struct MemsetConstValueWorker<1, true, c> {
	finline static void run(unsigned char* ptr) { *ptr = c; }
};

template <unsigned char c>
struct MemsetConstValueWorker<2, true, c> {
	finline static void run(unsigned char* ptr) { *(unsigned short *)ptr = (unsigned short)c << 8 + (unsigned short)c; }
};

template <unsigned char c>
struct MemsetConstValueWorker<3, true, c> {
	finline static void run(unsigned char* ptr) {
		MemsetConstValueWorker<2, true, c>::run(ptr);
		MemsetConstValueWorker<1, true, c>::run(ptr + 2);
	}
};

template <unsigned SIZE, unsigned char c>
finline void memset(void* ptr) {
	MemsetConstValueWorker<SIZE, SIZE <= 32, c>::run((unsigned char*)ptr);
}

}

#endif
