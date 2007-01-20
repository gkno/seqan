#ifndef SEQAN_HEADER_PRIORITY_TYPE_TREE_H
#define SEQAN_HEADER_PRIORITY_TYPE_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
/**
.Spec.PriorityHeap
..cat:Miscellaneous
..general:Class.PriorityType
..summary:Stores the priority data on a heap.
..signature:PriorityType<TValue, TComparator, PriorityHeap> >
*/
template < typename TValue, typename TComparator>
class PriorityType<TValue,TComparator,PriorityHeap>
{
public: 
	
	typedef String<TValue> THeap;
	
	TComparator cmp_;
	THeap heap_;


	//--------------------------------------------------------------
	PriorityType ()
		: cmp_(),heap_()
	{
SEQAN_CHECKPOINT
	}

	
//	PriorityType (PriorityType const & other_)
//		: cmp_(other_.cmp_),
//			heap_(other_.heap_)
//	{
//	}
		
//	inline PriorityType const &
//	operator = (PriorityType const & other_)
//	{
//		cmp_ = other_.cmp_;
//		heap = other_.heap_;
//		return *this;
//	}
	
	~PriorityType()
	{
SEQAN_CHECKPOINT
	}

}; // class PriorityType




// Empty the priority queue
///.Function.clear.param.object.type:Class.PriorityType
template <typename TValue, typename TComparator>
void clear (PriorityType<TValue,TComparator, PriorityHeap> & me)
{
	clear(me.heap_); 
}

// true if priority queue is empty 
///.Function.empty.param.object.type:Class.PriorityType
template <typename TValue, typename TComparator>
bool empty(PriorityType<TValue, TComparator, PriorityHeap> & me) 
{
SEQAN_CHECKPOINT
	return empty(me.heap_); 
}

// Number of elements in the priority queue
///.Function.length.param.object.type:Class.PriorityType
template <typename TValue, typename TComparator>
typename Size<PriorityType<TValue, TComparator, PriorityHeap> >::Type
length( PriorityType<TValue, TComparator, PriorityHeap> & me)
{ 
SEQAN_CHECKPOINT
	return length(me.heap_);
}






// Return the `best' element
/**
.Function.PriorityType#top:
..summary:Reference to the item with the highest priority.
..cat:Content Manipulation
..signature:top(object)
..param.object:A priority queue.
...type:Class.PriorityType
..remarks:To delete this item and adjust the priority queue use @Function.PriorityType#pop@.
..see:Function.PriorityType#pop
..see:Function.PriorityType#push
*/
template <typename TValue, typename TComparator>
TValue & 
top(PriorityType<TValue, TComparator, PriorityHeap> & me) 
{
SEQAN_CHECKPOINT
	return value(me.heap_,beginPosition(me.heap_));
}

// Copy heap position i to heap position h.
template <typename TValue, typename TComparator, typename TSize>
void _copyHeapElement ( PriorityType<TValue, TComparator, PriorityHeap> & me, const TSize & i, TSize & h)
{
SEQAN_CHECKPOINT
	_copyHeapElement ( me, me.heap_[i], h );
	h = i;
}

// Copy element to heap position h.
template <typename TValue, typename TComparator, typename TSize>
void _copyHeapElement (PriorityType<TValue, TComparator, PriorityHeap> & me, TValue element, TSize h)
{
SEQAN_CHECKPOINT
	me.heap_[h] = element;
}

/////////////////////////////////////////////////////////////////////////////////
//  lower priority of first element in queue 
/**
.Function.adjustTop
..cat:Miscellaneous
..summary:Adjusts the priority of the first item.
..param.object
...type:Class.PriorityType
*/
template <typename TValue, typename TComparator>
void adjustTop (PriorityType<TValue, TComparator, PriorityHeap> & me)	// so könnte man es dann auch nennen
{
SEQAN_CHECKPOINT
	if(!empty(me.heap_))
		_adjustHeapTowardLeaves(me,me.heap_[0],0,2);
}

////////////
//?? adjustHeapElement...
/////////////

/////////////////////////////////////////////////////////////////////////////////
/// Push a new element
/**
.Function.PriorityType#push:
..summary:Inserts a new item and adjusts the priority queue if necessary.
..cat:Content Manipulation
..signature:push(object, element)
..param.object:A priority queue.
...type:Class.PriorityType
..param.element:The item to be inserted in the priority queue.
...metafunction:Metafunction.Value
..remarks:The result of this operation is stored in $object$.
..see:Function.PriorityType#top
..see:Function.PriorityType#pop
*/
template <typename TValue, typename TComparator>
void push ( PriorityType<TValue, TComparator, PriorityHeap> & me, TValue element )
{
SEQAN_CHECKPOINT
	// root index is zero
	if ( empty(me.heap_) ) {
		resize(me.heap_,1);
		_copyHeapElement (me, element, 0 );
		return;
	}
	typedef typename Size<PriorityType<TValue, TComparator, PriorityHeap> >::Type TSize;
	TSize h = length(me); 
	resize(me.heap_,h+1);
	_adjustHeapTowardRoot(me, element, h); 
}

//////////////////////////////////////////////////////////////////////////////////
/// Priority got better.  Perform a cyclic shift along the tree edges toward root.
template <typename TValue, typename TComparator, typename TSize>
void _adjustHeapTowardRoot(PriorityType<TValue, TComparator, PriorityHeap> & me, TValue element, TSize h )
{
SEQAN_CHECKPOINT
	// root index is zero
	while ( h > 0) {
		const TSize i = (h-1)/2; 
		if ( me.cmp_ ( me.heap_[i], element ) )
		{
			_copyHeapElement ( me, i, h );
		}
		else
		{
			break;
		}
	}
	_copyHeapElement ( me, element, h );
}

//////////////////////////////////////////////////////////////////////////////
/// Pop 'best' element
/**
.Function.PriorityType#pop:
..summary:Deletes item with the highest priority and adjusts the priority queue.
..cat:Content Manipulation
..signature:pop(object)
..param.object:A priority queue.
...type:Class.PriorityType
..remarks:This function only deletes this item, but does not return it. To access the item use @Function.PriorityType#top@.
..see:Function.PriorityType#top
..see:Function.PriorityType#push
*/
template <typename TValue, typename TComparator>
void pop (PriorityType<TValue, TComparator, PriorityHeap> & me)
{
SEQAN_CHECKPOINT
	// root index is zero
	TValue element = getValue(me.heap_,endPosition(me.heap_)-1); 
	typedef typename Size<PriorityType<TValue, TComparator, PriorityHeap> >::Type TSize;
	TSize heap_size = length(me) - 1 ;
	resize (me.heap_,heap_size);
	if ( heap_size > 0 ) 
		_adjustHeapTowardLeaves(me, element, 0, 1 );

}


//////////////////////////////////////////////////////////////////////////////////
/// Priority got worse. Perform a cyclic shift along the tree edges toward leaves.
template <typename TValue, typename TComparator, typename TSize>
void _adjustHeapTowardLeaves(PriorityType<TValue, TComparator, PriorityHeap> & me, TValue element, TSize h, TSize i ) //für mich: h=0, i=1
{
SEQAN_CHECKPOINT
	// root index is zero
	const TSize heap_size = length(me);
	TComparator compare_ = me.cmp_;
	while ( i < heap_size )
	{
		if ( compare_ ( element, me.heap_[i] ) )
		{
			if ( compare_ ( me.heap_[i-1], me.heap_[i] ) )
			{
				_copyHeapElement ( me, i, h );
			}
			else
			{
				_copyHeapElement ( me, i-1, h );
			}
		}
		else
		{
			if ( compare_ ( element, me.heap_[i-1] ) )
			{
				_copyHeapElement ( me, i-1, h );
			}
			else
			{
				break;
			}
		}
		i = 2*h+2;
	}
	if ( i == heap_size && compare_ ( element, me.heap_[i-1] ) )
		_copyHeapElement ( me, i-1, h );
	_copyHeapElement ( me, element, h );
}


	//MetaFunctions

///.Metafunction.Size.param.T.type:Class.PriorityType
template < typename TValue, typename TComparator>
struct Size<PriorityType<TValue, TComparator, PriorityHeap> >
{
	typedef typename Size<typename PriorityType<TValue, TComparator, PriorityHeap>::THeap>::Type Type;
};

///.Metafunction.Value.param.T.type:Class.PriorityType
template < typename TValue, typename TComparator>
struct Value<PriorityType<TValue, TComparator, PriorityHeap> >
{
	typedef TValue Type;
};




////////////////////////
// debug
//template <typename TValue, typename THeap, typename TComparator>
//void check(PriorityType<TValue, THeap, TComparator> & me) { // debug
//	typedef typename Size<PriorityType<TValue, THeap, TComparator> >::Type TSize;
//	bool okay = true;
//	for ( TSize i = 1; i < length(me)-1; ++i )
//		if ( me.cmp_ ( me.heap_[(i-1)/2], me.heap_[i] ) ) {
//			cout << '\n' << (i-1)/2 << " < " << i << " : "<< (me.heap_[(i-1)/2]).value_ << " !< " << (me.heap_[i]).value_;
//			okay = false;
//		}
//	if ( okay )
//		cout << " ... seems okay\n";
//	else
//		cout << "\n... there were errors!\n";
//}

	//////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
