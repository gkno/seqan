#ifndef SEQAN_HEADER_PRIORITY_TYPE_BASE_H
#define SEQAN_HEADER_PRIORITY_TYPE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct PriorityHeap;

//////////////////////////////////////////////////////////////////////////////


/**
.Class.PriorityType:
..cat:Miscellaneous
..summary:Stores items in such a way that the item with the highest priority is at the top.
..signature:PriorityType<TValue, TComparator, TSpec>
..param.TValue:The value type that is stored.
...default:int
..param.TComparator:The comparator type that is used for sorting the items stored.
...default:std::less<TValue>
..param.TSpec:The specializing type.
...default:@Spec.PriorityHeap@
*/
template <typename TValue = int, typename TComparator = ::std::less<TValue>, typename TSpec = PriorityHeap>
class PriorityType;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TComparator, typename TSpec>
struct Value< PriorityType<TValue, TComparator, TSpec> >
{
	typedef TValue Type;
};

template <typename TValue, typename TComparator, typename TSpec>
struct Size< PriorityType<TValue, TComparator, TSpec> >
{
	typedef typename Size<TValue>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
