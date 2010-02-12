 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
 Author: Anne-Katrin Emde <emde@fu-berlin.de>
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_TYPES_H
#define SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_TYPES_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Interval Tree Types
//////////////////////////////////////////////////////////////////////////////

///---------------------------------------------------------------///

//////////////////// Interval and ID type ///////////////////
/**
.Class.IntervalAndCargo:
..cat:Miscellaneous
..summary:A simple record type that stores an interval and a cargo value.
..signature:IntervalAndCargo<TValue, TCargo>
..param.TValue:The value type, that is the type of the interval borders.
...default:int.
...metafunction:Metafunction.Value
..param.TCargo:The cargo type.
...default:int.
...metafunction:Metafunction.Cargo
..include:seqan/refinement.h
*/
template<typename TValue = int, typename TCargo = int>
class IntervalAndCargo
{
public:
    /**
.Memvar.IntervalAndCargo#i1:
..class:Class.PointAndCargo
..summary:The first element in the interval of type i1.
     */
	TValue i1;

    /**
.Memvar.IntervalAndCargo#i2:
..class:Class.PointAndCargo
..summary:The last element in the interval of type i2.
     */
	TValue i2;

    /**
.Memvar.IntervalAndCargo#cargo:
..class:Class.PointAndCargo
..summary:The stored cargo of type TCargo.
     */
	TCargo cargo;

    /**
.Memfunc.IntervalAndCargo#IntervalAndCargo:
..signature:IntervalAndCargo()
     */
    IntervalAndCargo()
    {
SEQAN_CHECKPOINT
    }

    /**
.Memfunc.IntervalAndCargo#IntervalAndCargo:
..class:Class.IntervalAndCargo
..summary:Constructor.
..signature:IntervalAndCargo(i1, i2, cargo)
..param.i1:The first element in the interval, of type TValue.
..param.i2:The last element in the interval of type TValue.
..param.cargo:The cargo value of type TCargo.
     */
	IntervalAndCargo(TValue i1, TValue i2, TCargo cargo):
		i1(i1), i2(i2), cargo(cargo)
	{
SEQAN_CHECKPOINT
	}
};


/**
.Function.leftBoundary:
..cat:Miscellaneous
..summary:Access to the right boundary.
..signature:leftBoundary(interval)
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The reference to the left boundary of the interval of type TValue&.
..see:Function.getLeftBoundary
..see:Function.rightBoundary
..see:Function.getRightBoundary
*/
template<typename TValue, typename TCargo>
TValue &
leftBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i1;
}


/**
.Function.rightBoundary:
..cat:Miscellaneous
..summary:Access to the right boundary.
..signature:leftBoundary(interval)
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The reference to the right boundary of the interval of type TValue&.
..see:Function.getRightBoundary
..see:Function.leftBoundary
..see:Function.getLeftBoundary
*/
template<typename TValue, typename TCargo>
TValue &
rightBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i2;
}


/**
.Function.getLeftBoundary:
..cat:Miscellaneous
..summary:Get method for the left boundary.
..signature:leftBoundary(interval)
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The copy of the left boundary of the interval of type TValue.
..see:Function.leftBoundary
..see:Function.getRightBoundary
..see:Function.rightBoundary
*/
template<typename TValue, typename TCargo>
TValue
getLeftBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i1;
}


/**
.Function.getRightBoundary:
..cat:Miscellaneous
..summary:Get method for the right boundary.
..signature:leftBoundary(interval)
..param.interval:The interval to return the right boundary for.
...type:Class.IntervalAndCargo
..returns:The copy of the right boundary of the interval of type TValue.
..see:Function.rightBoundary
..see:Function.getLeftBoundary
..see:Function.leftBoundary
*/
template<typename TValue, typename TCargo>
TValue
getRightBoundary(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.i2;
}


/**
.Function.cargo:
..signature:cargo(me)
..param.me:
...type:Class.IntervalAndCargo
..see:Function.getCargo
*/
template<typename TValue, typename TCargo>
TCargo &
cargo(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.cargo;
}

/**
.Function.getCargo:
..param.me:
...type:Class.IntervalAndCargo
..see:Function.cargo
*/
template<typename TValue, typename TCargo>
TCargo
getCargo(IntervalAndCargo<TValue,TCargo> & interval)
{
SEQAN_CHECKPOINT
	return interval.cargo;
}


/////////////////// Metafunctions //////////////////////
    
///.Metafunction.Value.param.T.type:Class.IntervalAndCargo
template<typename TValue,typename TCargo>
struct Value<IntervalAndCargo<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalAndCargo
template<typename TValue,typename TCargo>
struct Cargo<IntervalAndCargo<TValue,TCargo> >
{
	typedef TCargo Type;
};


/////////////////////// Point and ID type ////////////////
/**
.Class.PointAndCargo:
..cat:Miscellaneous
..summary:Simple record class storing a point (one-value interval) and a cargo.
..signature:PointAndCargo<TValue, TCargo>
..param.TValue:
...default:int.
...metafunction:Metafunction.Value
..param.TCargo:
...default:int.
...metafunction:Metafunction.Value
..include:seqan/refinement.h
*/
template<typename TValue=int, typename TCargo=int>
class PointAndCargo {
public:
    /**
.Memvar.PointAndCargo#point:
..class:Class.PointAndCargo
..summary:The stored point of type TValue.
     */
	TValue point;

    /**
.Memvar.PointAndCargo#cargo:
..class:Class.PointAndCargo
..summary:The stored cargo of type TCargo.
     */
	TCargo cargo;

    /**
.Memfunc.PointAndCargo#PointAndCargo
..class:Class.PointAndCargo
..signature:PointAndCargo(point, cargo)
    */
	PointAndCargo() {
SEQAN_CHECKPOINT
	}

    /**
.Memfunc.PointAndCargo#PointAndCargo
..class:Class.PointAndCargo
..summary:Constructor.
..signature:PointAndCargo(point, cargo)
..param.point:
...summary:The point to store of type TValue.
..param.cargo:
...summary:The cargo to store of type TCargo.
    */
	PointAndCargo(TValue point, TCargo cargo):
		point(point), cargo(cargo)
	{
SEQAN_CHECKPOINT
	}
};


/**
.Function.leftBoundary:
..signature:leftBoundary(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TValue &
leftBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}


/**
.Function.rightBoundary:
..signature:rightBoundary(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TValue &
rightBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}


/**
.Function.getLeftBoundary:
..signature:getLeftBoundary(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TValue
getLeftBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}


/**
.Function.getRightBoundary:
..signature:getRightBoundary(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TValue
getRightBoundary(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.point;
}


/**
.Function.cargo:
..signature:cargo(point)
..param.point.type:Class.PointAndCargo
 */
template<typename TValue, typename TCargo>
TCargo &
cargo(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.cargo;
}


/**
.Function.cargo:
..signature:getCargo(point)
..param.point:
...type:Class.IntervalAndCargo
..see:Function.getCargo
*/
template<typename TValue, typename TCargo>
TCargo
getCargo(PointAndCargo<TValue,TCargo> & point)
{
SEQAN_CHECKPOINT
	return point.cargo;
}



////////////////// Metafunctions //////////////////
///.Metafunction.Value.param.T.type:Class.PointAndCargo
template<typename TValue,typename TCargo>
struct Value<PointAndCargo<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.PointAndCargo
template<typename TValue,typename TCargo>
struct Cargo<PointAndCargo<TValue,TCargo> >
{
	typedef TCargo Type;
};


//// Comparators
template <typename TPair>
bool _less_compI1_ITree(TPair const& p1, TPair const& p2){
SEQAN_CHECKPOINT
  return (leftBoundary(const_cast<TPair&>(p1)) < leftBoundary(const_cast<TPair&>(p2)));
}


template <typename TPair>
bool _greater_compI2_ITree(TPair const& p1, TPair const& p2){
SEQAN_CHECKPOINT
  return (rightBoundary(const_cast<TPair&>(p1)) > rightBoundary(const_cast<TPair&>(p2)));
}


///////////////////////////////////////////////////////////////////////////
/////////////////////////// IntervalTreeNode	///////////////////////////

/**
.Tag.IntervalTree Node Types
..summary:Tags to select the node type for @Class.IntervalTree@.
..cat:Miscellaneous

..tag.StorePointsOnly:The tree nodes store points.
*/
struct StorePointsOnly {};


///..tag.StoreIntervals:The tree nodes store intervals.
struct StoreIntervals {};


/**
.Class.IntervalTreeNode:
..cat:Miscellaneous
..summary:Element of @Class.IntervalTree@.
..signature:IntervalTreeNode<TInterval, TSpec>
..param.TInterval:The type of interval to store.
..param.TSpec:The type of interval to store.
...default:StorePointsOnly.
...metafunction:Metafunction.Spec
..include: seqan/refinement.h

.Memvar.IntervalTreeNode#center:
..class:Class.IntervalTreeNode
..summary:The center of the interval of type TValue.

.Memvar.IntervalTreeNode#list1
..class:Class.IntervalTreeNode
..summary:Sorted list of pointers to intervals, sorted in ascending according to their left boundary points.

.Memvar.IntervalTreeNode#list2
..class:Class.IntervalTreeNode
..summary:Sorted list of pointers to intervals, sorted in descending according to their right boundary points.
 */
template<typename TInterval, typename TSpec=StorePointsOnly>
class IntervalTreeNode;


/**
.Spec.Interval Tree Node
..cat:Miscellaneous
..general:Class.IntervalTreeNode
..summary:An Interval Tree Node that stores intervals explicitely in each node.
..signature:IntervalTreeNode<TInterval, StoreIntervals>
..param.TInterval:The interval type to store in the node.
 */
template<typename TInterval>
class IntervalTreeNode<TInterval, StoreIntervals> {
public:
	typedef typename Value<TInterval>::Type TValue;

	TValue center;
	String<TInterval> list1;
	String<TInterval> list2;
};


/**
.Spec.Points Only Tree Node
..cat:Miscellaneous
..general:Class.IntervalTreeNode
..summary:Spec for IntervalTreeNode that stores only the relevant point in each node meaning the endpoint of the interval in the list sorted by endpoints (list2) and only the beginpoint of the interval in the list sorted by beginpoints (list1).
..signature:IntervalTreeNode<TInterval, StorePointsOnly>
..param.TInterval:The interval type to store in the node.
 */
template<typename TInterval>
class IntervalTreeNode<TInterval, StorePointsOnly> {
public:
	typedef typename Cargo<TInterval>::Type TCargo;
	typedef typename Value<TInterval>::Type TValue;

	TValue center;
	String<PointAndCargo<TValue,TCargo> > list1;
	String<PointAndCargo<TValue,TCargo> > list2;

    /**
.Memfunc.IntervalTreeNode#IntervalTreeNode:
..class:Class.IntervalTreeNode
..summary:Default constructor.
..signature:IntervalTreeNode()
     */
    IntervalTreeNode()
    {
SEQAN_CHECKPOINT
    }

	IntervalTreeNode(IntervalTreeNode const & other):
		center(other.center),
		list1(other.list1),
		list2(other.list2)
	{
SEQAN_CHECKPOINT
	}
};


//internal set functions
template<typename TValue, typename TInterval>
void
_setIntervalTreeNode(IntervalTreeNode<TInterval,StoreIntervals> & knot,TValue center,TInterval & interval)
{
SEQAN_CHECKPOINT
	
	knot.center = center;
	appendValue(knot.list1,interval);
	appendValue(knot.list2,interval);

}

template<typename TInterval>
void
_appendIntervalTreeNodeLists(IntervalTreeNode<TInterval,StoreIntervals> & knot,TInterval & interval)
{
SEQAN_CHECKPOINT
	
	appendValue(knot.list1,interval);
	appendValue(knot.list2,interval);

}


//internal set functions
template<typename TValue, typename TInterval>
void
_setIntervalTreeNode(IntervalTreeNode<TInterval,StorePointsOnly> & knot,TValue center,TInterval & interval)
{
SEQAN_CHECKPOINT
	
	knot.center = center;
	appendValue(knot.list1,PointAndCargo<TValue,typename Cargo<TInterval>::Type>(leftBoundary(interval),cargo(interval)));
	appendValue(knot.list2,PointAndCargo<TValue,typename Cargo<TInterval>::Type>(rightBoundary(interval),cargo(interval)));
	

}


template<typename TInterval>
void
_appendIntervalTreeNodeLists(IntervalTreeNode<TInterval,StorePointsOnly> & knot,TInterval & interval)
{
SEQAN_CHECKPOINT
	
	appendValue(knot.list1,PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type>(leftBoundary(interval),cargo(interval)));
	appendValue(knot.list2,PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type>(rightBoundary(interval),cargo(interval)));
	

}

/////////////////// Metafunctions ///////////////////////
///.Metafunction.Value.param.T.type:Class.IntervalTreeNode
template<typename TInterval, typename TSpec>
struct Value<IntervalTreeNode<TInterval,TSpec> >
{
	typedef typename Value<TInterval>::Type Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalTreeNode
template<typename TInterval, typename TSpec>
struct Cargo<IntervalTreeNode<TInterval,TSpec> >
{
	typedef typename Cargo<TInterval>::Type Type;
};


/**
.Metafunction.ListType:
..cat:Miscellaneous
..signature:ListType<T>::ListType
..summary:Type of lists in tree nodes.
..param.T:The type to retrieve the list type for.
..returns:Returns the type of the the lists in @Class.IntervalTreeNode@ objects.
 */
template<typename T>
struct ListType;


///.Metafunction.ListType.param.T.type:Class.IntervalTreeNode
template<typename TInterval>
struct ListType<IntervalTreeNode<TInterval,StorePointsOnly> >
{
	typedef String<PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type> > Type;

};


///.Metafunction.ListType.param.T.type:Class.IntervalTreeNode
template<typename TInterval>
struct ListType<IntervalTreeNode<TInterval,StoreIntervals> >
{
	typedef String<IntervalAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type> > Type;

};

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // #ifndef SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_TYPES_H
