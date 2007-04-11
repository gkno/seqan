#ifndef SEQAN_HEADER_PRIORITY_TYPE_BASE_H
#define SEQAN_HEADER_PRIORITY_TYPE_BASE_H

#include <functional>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct PriorityHeap;

//////////////////////////////////////////////////////////////////////////////


/**
.Class.PriorityType:
..cat:Miscellaneous
..summary:Stores items in such a way that the item with the highest priority is at the top.
..signature:PriorityType<TValue, TLess, TSpec>
..param.TValue:The value type that is stored.
...default:int
..param.TLess:The comparator type that is used for sorting the items stored.
...default:std::less<TValue>
..param.TSpec:The specializing type.
...default:@Spec.PriorityHeap@
*/
template <typename TValue = int, typename TLess = ::std::less<TValue>, typename TSpec = PriorityHeap>
class PriorityType;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TLess, typename TSpec>
struct Value< PriorityType<TValue, TLess, TSpec> >
{
	typedef TValue Type;
};

template <typename TValue, typename TLess, typename TSpec>
struct Size< PriorityType<TValue, TLess, TSpec> >
{
	typedef typename Size<TValue>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
