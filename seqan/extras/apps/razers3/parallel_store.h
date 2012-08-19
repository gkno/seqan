#ifndef APPS_RAZERS_PARALLEL_STORE_H
#define APPS_RAZERS_PARALLEL_STORE_H

#ifdef PLATFORM_GCC_DIS
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#include <numeric>
#endif  // #ifdef PLATFORM_GCC_DIS

using namespace seqan;

struct Parallel_;
typedef Tag<Parallel_> Parallel;

#ifdef PLATFORM_GCC_DIS

/*
template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign& alignStore, Tag<TSortSpec> const &, Parallel const &) 
{
  __gnu_parallel::sort(
		begin(alignStore, Standard() ), 
		end(alignStore, Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign const & alignStore, Tag<TSortSpec> const &, Parallel const &) 
{
  __gnu_parallel::sort(
		begin(const_cast<TAlign&>(alignStore), Standard() ), 
		end(const_cast<TAlign&>(alignStore), Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}
*/

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const &less, Parallel const &) 
{
  __gnu_parallel::sort(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const &less, Parallel const &) 
{
  __gnu_parallel::sort(
		begin(const_cast<TAlign&>(alignStore), Standard()), 
		end(const_cast<TAlign&>(alignStore), Standard()), 
		less);
}

template <typename TIntString>
inline void
partialSum(TIntString &intString)
{
    __gnu_parallel::partial_sum(
        begin(intString, Standard()), 
        end(intString, Standard()), 
        begin(intString, Standard()));
}

#else  // #ifdef PLATFORM_GCC

/*
template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign& alignStore, Tag<TSortSpec> const & tag, Parallel const &) 
{
  sortAlignedReads(alignStore, tag);
}

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign const & alignStore, Tag<TSortSpec> const & tag, Parallel const &) 
{
  sortAlignedReads(alignStore, tag);
}
*/

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign & alignStore, TFunctorLess const &less, Parallel const &) 
{
  ::std::sort(begin(alignStore, Standard()), end(alignStore, Standard()), less);
//  sortAlignedReads(alignStore, less);
}

template <typename TAlign, typename TFunctorLess>
inline void
sortAlignedReads(TAlign const & alignStore, TFunctorLess const &less, Parallel const &) 
{
  ::std::sort(begin(alignStore, Standard()), end(alignStore, Standard()), less);
//  sortAlignedReads(alignStore, less);
}

template <typename TIntString>
inline void
partialSum(TIntString &intString)
{
    ::std::partial_sum(
        begin(intString, Standard()), 
        end(intString, Standard()), 
        begin(intString, Standard()));
}

#endif  // #ifdef PLATFORM_GCC

#endif  // ifndef APPS_RAZERS_PARALLEL_STORE_H

