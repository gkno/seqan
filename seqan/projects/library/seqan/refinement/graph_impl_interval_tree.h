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

#ifndef SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_H
#define SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_H

// TODO(holtgrew): Remove crufty code?

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Interval Tree
//////////////////////////////////////////////////////////////////////////////

/**
.Class.IntervalTree:
..cat:Miscellaneous
..summary:A datastructure that efficiently stores intervals.
..signature:IntervalTree<TValue, TCargo>
..param.TValue:The value type.
..param.TCargo:The cargo/id type.
...default:int
...remarks:If the intervals are not associated with cargos/IDs, they will be numbered consecutively.
..include:seqan/refinement.h
*/
template<typename TValue=int, typename TCargo=unsigned int>
class IntervalTree
{
public:
	typedef Graph<Directed<void,WithoutEdgeId> > TGraph;
	typedef IntervalAndCargo<TValue,TCargo> TInterval;
	typedef IntervalTreeNode<TInterval> TNode;
	typedef String<TNode> TPropertyMap;

	TGraph g;
	TPropertyMap pm;
	size_t interval_counter;
	
	
	IntervalTree()
	{
SEQAN_CHECKPOINT
		interval_counter = 0;
	}
	
	/**
.Memfunc.IntervalTree#IntervalTree
..class:Class.IntervalTree
..summary:Constructor
..signature:IntervalTree(intervalBegins, intervalEnds, intervalCargos, len)
..param.intervalBegins:Iterator pointing to begin position of first interval.
..param.intervalEnds:Iterator pointing to end position of first interval.
..param.intervalCargos:Iterator pointing to cargos/ids for intervals.
..param.len:Number of intervals to store in tree.
     */
	template<typename TIterator,typename TCargoIterator>
	IntervalTree(TIterator interval_begins,
				 TIterator interval_ends, 
				 TCargoIterator interval_cargos, 
				 size_t len)	
	{
SEQAN_CHECKPOINT
		String<TInterval> intervals;
		resize(intervals,len);
		size_t i = 0;
		while(i<len)
		{
			intervals[i].i1 = value(interval_begins);
			++interval_begins;
			intervals[i].i2 = value(interval_ends);
			++interval_ends;
			intervals[i].cargo = value(interval_cargos);
			++interval_cargos;
			++i;
		}
		interval_counter = len;
		
		createIntervalTree(g,pm,intervals);
	}

	/**
..signature:IntervalTree(intervalBegins, intervalEnds, len)
     */
	template<typename TIterator>
	IntervalTree(TIterator interval_begins,
				 TIterator interval_ends,
				 size_t len)
	{
SEQAN_CHECKPOINT
		String<TInterval> intervals;
		resize(intervals,len);
		size_t i = 0;
		while(i<len)
		{
			intervals[i].i1 = value(interval_begins);
			++interval_begins;
			intervals[i].i2 = value(interval_ends);
			++interval_ends;
			intervals[i].cargo = i;
			++i;
		}
		interval_counter = len;
		createIntervalTree(g,pm,intervals);
	}
	
	/**
..signature:IntervalTree(String<TInterval> intervals)
     */
	IntervalTree(String<TInterval> intervals)	
	{
SEQAN_CHECKPOINT
		interval_counter = length(intervals);
		createIntervalTree(g,pm,intervals);
	}

	/**
..signature:IntervalTree(String<TInterval> intervals, Tag<TTagSpec> const tag)
     */
	template <typename TTagSpec>
	IntervalTree(String<TInterval> intervals, Tag<TTagSpec> const tag)
	{
SEQAN_CHECKPOINT
		interval_counter = length(intervals);
		createIntervalTree(g,pm,intervals,tag);
	}

    /**
..signature:IntervalTree(String<TInterval> intervals, TValue center)
    */
	IntervalTree(String<TInterval> intervals, TValue center)	
	{
SEQAN_CHECKPOINT
		interval_counter = length(intervals);
		createIntervalTree(g,pm,intervals,center);
	}
};



///////Specs for the way interval centers are determined
/**
.Tag.IntervalTree Centers
..cat:Miscellaneous
..summary:Tag to select a specific way to compute the center of an interval tree node.
..see:Class.IntervalTree
..include:seqan/refinement.h
 */


/**
..tag.ComputeCenter
...summary:For intervals that are more or less uniformly distributed in the value range, using the ComputeCenter tag may result in a more balanced tree compared to using the RandomCenter tag.
...signature:ComputeCenter
...remarks:center = minbegin + (maxend-minbegin)/2
 */
//template <typename TSpec = SpecPointAndCargo>
struct _TagComputeCenter;
typedef Tag<_TagComputeCenter> const ComputeCenter;


/**
..tag.RandomCenter
...summary:The RandomCenter tag guarantees that each node contains at least one interval, therefore the size of the tree is limited by the nummer of intervals. This may lead to an unbalanced tree, but is the most space-efficient and in practice the fastest method.
...signature:RandomCenter
...remarks:center = center of random interval
 */
//template <typename TSpec = SpecPointAndCargo>
struct _TagRandomCenter;
typedef Tag<_TagRandomCenter> const RandomCenter;

// fill the container interval_pointers with pointers to the corresponding objects in intervals.
// this is done to avoid copying and passing the whole IntervalAndCargo objects during interval tree construction
template<typename TIntervals, typename TIntervalPointers>
void
_makePointerInterval(TIntervals & intervals,TIntervalPointers & interval_pointers)
{
SEQAN_CHECKPOINT
	typedef typename Value<TIntervalPointers>::Type TIntervalPointer;
	typedef typename Iterator<TIntervalPointers, Rooted>::Type TIntervalPointerIterator;

	TIntervalPointer it;
	TIntervalPointerIterator iit = begin(interval_pointers);
	if(length(intervals)>0)
		for(it = &intervals[0]; it <= &intervals[length(intervals)-1]; ++it)
		{
			*iit = it;	
			++iit;
		}

}


/**
.Function.createIntervalTree
..summary:Create an interval tree.
..cat:Miscellaneous
..signature:createIntervalTree(TGraph &g, TPropertyMap &pm, TIntervals &intervals, Tag<TSpec> const tag)
..param.g:DirectedGraph to create interval tree in.
...type:Class.Graph
..param.pm:Property map to use for the created interval tree.
...type:Class.PropertyMap
..param.intervals:Container of intervals.
...type:Class.String
...remark:Should be a String of @Class.Interval@ or @Class.IntervalAndCargo@ objects.
..param.tag:Tag for tree construction method. @tag.RandomCenter@ or @tag.ComputeCenter@
..remark:center of root node is computed by _calcIntervalTreeRootCenter
..include:seqan/refinement.h
 */
template<typename TGraph, typename TPropertyMap, typename TIntervals, typename TSpec>
void 
createIntervalTree(TGraph & g, 
				   TPropertyMap & pm, 
				   TIntervals & intervals,
				   Tag<TSpec> const tag)
{
SEQAN_CHECKPOINT
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TIntervals>::Type TInterval;
	typedef typename Value<TInterval>::Type TValue;

	reserve(g.data_vertex,length(intervals));
	reserve(pm,length(intervals));
	
	TVertexDescriptor root = addVertex(g);
	resizeVertexMap(g,pm);
	
	TValue center =	_calcIntervalTreeRootCenter(intervals);
	
	std::sort(begin(intervals),end(intervals),_less_compI1_ITree<TInterval>);

	String<TInterval*> interval_pointers;
	resize(interval_pointers,length(intervals));
	_makePointerInterval(intervals,interval_pointers);

	_createIntervalTree(g,pm,interval_pointers,root,(TValue)0.0,center,length(intervals),tag);
	reserve(pm, length(pm), Exact());
	reserve(g.data_vertex, length(g.data_vertex), Exact());

}


/**
..signature:createIntervalTree(TGraph &g, TPropertyMap &pm, TIntervals &intervals)
..param.tag.default:Tag.IntervalTree Centers.tag.ComputeCenter
 */
// most user friendly interval tree construction for the moment...
// RandomCenter tag as default construction method
template<typename TGraph, typename TPropertyMap, typename TIntervals>
void 
createIntervalTree(TGraph & g, 
				   TPropertyMap & pm, 
				   TIntervals & intervals)
{
SEQAN_CHECKPOINT
	createIntervalTree(g,pm,intervals,RandomCenter());
}


/**
..signature:createIntervalTree(TGraph &g, TPropertyMap &pm, TIntervals &intervals, center, tag)
 */
template<typename TGraph, typename TPropertyMap, typename TIntervals, typename TSpec>
void 
createIntervalTree(TGraph & g,
 TPropertyMap & pm, 
				   TIntervals & intervals,
				   typename Value<typename Value<TIntervals>::Type>::Type center,
				   Tag<TSpec> const tag)
{
SEQAN_CHECKPOINT
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TIntervals>::Type TInterval;
	typedef typename Value<typename Value<TIntervals>::Type>::Type TValue;
	
	reserve(g.data_vertex,length(intervals));
	reserve(pm,length(intervals));
	
	TVertexDescriptor root = addVertex(g);
	resizeVertexMap(g,pm);
	
	TInterval a;
	typename Iterator<TIntervals, Standard>::Type begin_ = begin(intervals);
	typename Iterator<TIntervals, Standard>::Type end_ = end(intervals);
	std::sort(begin_, end_ ,_less_compI1_ITree<TInterval>);

	String<TInterval*> interval_pointers;
	resize(interval_pointers,length(intervals));

	_makePointerInterval(intervals,interval_pointers);

	if(length(intervals)==1)
		center = (rightBoundary(intervals[0])-leftBoundary(intervals[0]))/(TValue)2.0;

	_createIntervalTree(g,pm,interval_pointers,root,(TValue)0.0,center,length(intervals),tag);
		
	reserve(pm, length(pm), Exact());
	reserve(g.data_vertex, length(g.data_vertex), Exact());

}

/**
..signature:createIntervalTree(TGraph &g, TPropertyMap &pm, TIntervals &intervals, center)
 */
// RandomCenter tag as default construction method
template<typename TGraph, typename TPropertyMap, typename TIntervals>
void 
createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   TIntervals & intervals, 
				   typename Value<typename Value<TIntervals>::Type>::Type center)
{
SEQAN_CHECKPOINT
	createIntervalTree(g,pm,intervals,center,RandomCenter());
}



//////////////////////////////////////////////////////////////////////////////
//remembers minimum and maximum of point values in intervals and sets the center
//of each node to min+(max-min)/2
template<typename TGraph, typename TPropertyMap, typename TIntervalPointer, typename TValue>
void 
_createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   String<TIntervalPointer*> & intervals,
				   typename VertexDescriptor<TGraph>::Type & knot, 
				   TValue, 
				   TValue center, 
				   typename VertexDescriptor<TGraph>::Type len,
				   Tag<_TagComputeCenter> const tag)
{
SEQAN_CHECKPOINT
	//  Rekursionsanker
	if(len==1){
		_setIntervalTreeNode(value(pm,knot),center,*intervals[0]);
		return;
	}

	typedef typename Value<TPropertyMap>::Type TNode;
	typedef typename ListType<TNode>::Type TList;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TIntervalPointer*> TIntervalPointers;

	// one list of interval pointers for the intervals to the left of center
	TIntervalPointers S_left;
	// one list of interval pointers for the intervals to the right of center
	TIntervalPointers S_right;

	TValue min1 = supremumValue<TValue>();
	TValue min2 = supremumValue<TValue>();
	TValue max1 = infimumValue<TValue>();
	TValue max2 = infimumValue<TValue>();

	value(pm,knot).center = center;
	
 
	typedef typename Iterator<TIntervalPointers,Standard>::Type TIntervalIterator;
	TIntervalIterator it = begin(intervals,Standard());
	TIntervalIterator it_end = end(intervals,Standard());
	
	// walk through intervals
	while(it != it_end)
	{
		// interval belongs to the left list
		if((**it).i2<=center)
		{
			appendValue(S_left,*it, Generous());
			 //remember right most and left most point in left list
			if((**it).i2 > max1)
				max1 = (**it).i2;
			if((**it).i1 < min1)
				min1 = (**it).i1;
		}
		else
		{
			// interval belongs to the right list
			if((**it).i1>center)
			{
				appendValue(S_right,(*it), Generous());
				 //remember right most and left most point in right list
				if((**it).i2 > max2)
					max2 = (**it).i2;
				if ((**it).i1 < min2)
					min2 = (**it).i1;
			}
			else // interval belongs to this node
			{
				_appendIntervalTreeNodeLists(value(pm,knot),**it);
			}
		}
        ++it;
	}

//	std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
	std::sort(begin(value(pm,knot).list2),end(value(pm,knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);

	// build subtree to the left
	if(!empty(S_left))
	{
		TVertexDescriptor vd = addVertex(g);
		resize(pm, vd+1); 
		addEdge(g,knot,vd);
		_createIntervalTree(g,pm,S_left,vd,center,min1+(max1-min1)/2,length(S_left),tag);
	}
	// build subtree to the right
	if(!empty(S_right))
	{
		TVertexDescriptor vd = addVertex(g);
		resize(pm, vd+1); 
		addEdge(g,knot,vd);
		_createIntervalTree(g,pm,S_right,vd,center,min2+(max2-min2)/2,length(S_right),tag);
	}
}




//////////////////////////////////////////////////////////////////////////////
//createIntervalTree for all specs except CompCenter, the center value of each 
//node is determined by functions _calcIntervalTreeNodeCenterLeft and 
//_calcIntervalTreeNodeCenterRight
template<typename TGraph, typename TPropertyMap, typename TSpec, typename TInterval, typename TValue>
void 
_createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   String<TInterval*> & intervals, 
				   typename VertexDescriptor<TGraph>::Type & knot, 
				   TValue last_center, TValue center, 
				   typename VertexDescriptor<TGraph>::Type len,
				   Tag<TSpec> const tag)
{
SEQAN_CHECKPOINT
	// Rekursionsanker
	if(len==1){
		_setIntervalTreeNode(value(pm,knot),center,*value(intervals,0));
		return;
	}

	typedef typename Value<TPropertyMap>::Type TNode;
	typedef typename ListType<TNode>::Type TList;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TInterval*> TIntervalPointers;
	
	// one list of interval pointers for the intervals to the left of center
	TIntervalPointers S_left;
	// one list of interval pointers for the intervals to the right of center
	TIntervalPointers S_right;
		
	value(pm,knot).center = center;
	
	typedef typename Iterator<TIntervalPointers,Standard>::Type TIntervalIterator;
	TIntervalIterator it = begin(intervals,Standard());
	TIntervalIterator it_end = end(intervals,Standard());
	
	// walk through intervals
	while(it != it_end)
	{
		// interval belongs to the left list
		if((**it).i2<=center)
		{
			appendValue(S_left,*it, Generous());
		}
		else
		{	// interval belongs to the right list
			if((**it).i1>center)
			{
				appendValue(S_right,(*it), Generous());
			}
			else
			{
				// interval belongs to the current node
				_appendIntervalTreeNodeLists(value(pm,knot),**it);
			}
		}
        ++it;
	}

//	std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
	std::sort(begin(value(pm,knot).list2),end(value(pm,knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);

	// build subtree to the left
	if(!empty(S_left))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		addEdge(g,knot,vd);
		TValue next_center = _calcIntervalTreeNodeCenterLeft(S_left,last_center,center,tag);
		_createIntervalTree(g,pm,S_left,vd,center,next_center,length(S_left),tag);
	}
	// build subtree to the right
	if(!empty(S_right))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		addEdge(g,knot,vd);
		TValue next_center = _calcIntervalTreeNodeCenterRight(S_right,last_center,center,tag);
		_createIntervalTree(g,pm,S_right,vd,center,next_center,length(S_right),tag);
	}
}




//the RandomCenter spec way of chosing center values:
//pick a random interval from the list and take its center as the center value 
//for the left child node (during interval tree construction)
template<typename TIntervals, typename TValue>
TValue
_calcIntervalTreeNodeCenterLeft(TIntervals & intervals, TValue &, TValue &, Tag<_TagRandomCenter> const)
{
SEQAN_CHECKPOINT
	TValue rand_index = rand()%length(intervals);  
	return (rightBoundary(*value(intervals,rand_index))+leftBoundary(*value(intervals,rand_index)))/(TValue)2.0;
}

//the RandomCenter spec way of chosing center values:
//pick a random interval from the list and take its center as the center value 
//for the right child node (during interval tree construction)
template<typename TIntervals, typename TValue>
TValue
_calcIntervalTreeNodeCenterRight(TIntervals & intervals, TValue &, TValue &, Tag<_TagRandomCenter> const)
{
SEQAN_CHECKPOINT
	TValue rand_index = rand()%length(intervals);  
	return (rightBoundary(*value(intervals,rand_index))+leftBoundary(*value(intervals,rand_index)))/(TValue)2.0;
}


// if the center of the root is not given, it is placed in the "ComputeCenter way": in the middle of minValue and maxValue
// where minValue is the minimum left boundary and maxValue is the maximum right boundary of all intervals
template<typename TIntervals>
typename Value<typename Value<TIntervals>::Type>::Type
_calcIntervalTreeRootCenter(TIntervals & intervals)
{
SEQAN_CHECKPOINT
	
	typedef typename Value<typename Value<TIntervals>::Type>::Type TValue;
	typedef typename Iterator<TIntervals,Standard>::Type TIntervalIterator;

	TIntervalIterator it = begin(intervals);
	TIntervalIterator it_end = end(intervals);

	TValue min = supremumValue<TValue>();
	TValue max = infimumValue<TValue>();

	while(it != it_end)
	{
		if(leftBoundary(*it)<min) min = leftBoundary(*it);
		if(rightBoundary(*it)>max) max = rightBoundary(*it);
		++it;
	}
	
	return (min+(max-min)/(TValue)2.0);

}



/**
.Function.addInterval
..cat:Miscellaneous
..signature:addInterval(graph, propertyMap, interval)
..param.graph:The directed graph that contains the topography of the interval tree.
..param.propertyMap:The property map containing the node properties of the interval tree.
..param.interval:The interval to be added to the interval tree.
..summary:Adds an interval to an interval tree.
..include:seqan/refinement.h
*/
template<typename TGraph, typename TPropertyMap, typename TInterval>
void
addInterval(TGraph & g, TPropertyMap & pm, TInterval interval)
{
SEQAN_CHECKPOINT

	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Value<TInterval>::Type TValue;
	typedef typename ListType<TProperty>::Type TList;
	

	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if(act_prop.center < leftBoundary(interval))
		{
			if(atEnd(it)){
				TVertexDescriptor vd = addVertex(g);
				resizeVertexMap(g,pm);
				addEdge(g,act_knot,vd);
				_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/(TValue)2.0,interval);
				break;
			}
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)){
						TVertexDescriptor vd = addVertex(g);
						resizeVertexMap(g,pm);
						addEdge(g,act_knot,vd);
						_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/(TValue)2.0,interval);
						break;
					}
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(rightBoundary(interval) <= act_prop.center)
			{
				if(atEnd(it)){
					TVertexDescriptor vd = addVertex(g);
					resizeVertexMap(g,pm);
					addEdge(g,act_knot,vd);
					_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/2,interval);
					break;
				}
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)){
							TVertexDescriptor vd = addVertex(g);
							resizeVertexMap(g,pm);
							addEdge(g,act_knot,vd);
							_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/2,interval);
							break;
						}
					}
				}
				act_knot = targetVertex(it);
			}
			else{
				_appendIntervalTreeNodeLists(property(pm, act_knot),interval);
				std::sort(begin(property(pm,act_knot).list1),end(property(pm,act_knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
				std::sort(begin(property(pm,act_knot).list2),end(property(pm,act_knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);
				break;
			}
		}
	}

}


/**
..signature:addInterval(intervalTree, interval)
..param.intervalTree:The interval tree to add the interval to.
...type:Class.IntervalTree
 */
template<typename TValue, typename TCargo, typename TInterval>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TInterval interval)
{
SEQAN_CHECKPOINT

	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}


// TODO(holtgrewe): Is this begin/end in C++ style or is it first/last?
/**
..signature:addInterval(intervalTree, begin, end, cargo)
..param.begin:Begin position of interval of type TValue.
..param.end:End position of interval of type TValue.
..param.cargo:Cargo to attach to the interval.
 */
template<typename TValue, typename TCargo>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TValue begin, TValue end, TCargo cargo)
{
SEQAN_CHECKPOINT

	IntervalAndCargo<TValue,TCargo> interval;
	interval.i1 = begin;
	interval.i2 = end;
	interval.cargo = cargo;
	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}


/**
..signature:addInterval(intervalTree, begin, end)
 */
template<typename TValue, typename TCargo>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TValue begin, TValue end)
{
SEQAN_CHECKPOINT

	IntervalAndCargo<TValue,TCargo> interval;
	interval.i1 = begin;
	interval.i2 = end;
	interval.cargo = itree.interval_counter;
	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}


/**
.Function.findIntervals
..summary:Find all intervals that contain the query point or overlap with the query interval.
..signature:findIntervals(graph, propertyMap, query, result)
..param.query:A query point.
..include:seqan/refinement.h
*/
template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervals(TGraph & g, TPropertyMap & pm, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	resize(result,0);

	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if(act_prop.center < query)
		{
			unsigned int i = 0;
			while(i < length(act_prop.list2) && rightBoundary(value(act_prop.list2,i)) > query)
			{
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
				++i;	
			}
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(query < act_prop.center)
			{
				unsigned int i = 0;
				while(i < length(act_prop.list1) && leftBoundary(value(act_prop.list1,i)) <= query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
					++i;
				}
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{
				for(unsigned int i = 0; i < length(act_prop.list1); ++i)
                    appendValue(result, cargo(value(act_prop.list1,i)), Generous());
				break;
			}
		}
	}

}


/**
..signature:findIntervals(intervalTree, query, result)
..param.intervalTree:An interval tree
...type:Class.IntervalTree
*/
template<typename TValue,typename TCargo>
void
findIntervals(IntervalTree<TValue,TCargo> & it, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	findIntervals(it.g,it.pm,query,result);

}


/**
.Function.findIntervalsExcludeTouching
..cat:Miscellaneous
..summary::Find all intervals that contain the query point, exclude intervals that touch the query, i.e. where the query point equals the start or end point.
..signature:findIntervalsExcludeTouching(graph, propertyMap, query, result)
..param.graph:The directed graph that contains the topography of the interval tree.
..param.propertyMap:The property map containing the node properties of the interval tree
..param.query:The TValue to query here.
..param.result:The resulting string of cargos/ids of the intervals that contain the query point.
...type:Class.String
...remark:Should be a string of TCargo.
..include:seqan/refinement.h
 */
template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervalsExcludeTouching(TGraph & g, TPropertyMap & pm, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename Iterator<TGraph, OutEdgeIterator >::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	
	resize(result,0);

	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if( (TValue) act_prop.center < query)
		{
			int i = 0;
			while(i < (int) length(act_prop.list2) && (TValue) rightBoundary(value(act_prop.list2,i)) > query)
			{
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
				++i;	
			}
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(query < (TValue) act_prop.center)
			{
				int i = 0;
				while(i < (int) length(act_prop.list1) && (TValue) leftBoundary(value(act_prop.list1,i)) < query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
					++i;
				}
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{
				int i = 0;
				while(i < (int) length(act_prop.list1) && (TValue) leftBoundary(value(act_prop.list1,i)) < query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
					++i;
				}
				break;
			}
		}
	}

}


/**
..signature:findIntervalsExcludeTouching(intervalTree, query, result)
..param.intervalTree:An interval tree
...type:Class.IntervalTree
*/
template<typename TValue,typename TCargo>
void
findIntervalsExcludeTouching(IntervalTree<TValue,TCargo> & tree, TValue query, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	findIntervalsExcludeTouching(tree.g,tree.pm,query,result);

}




/**
.Function.findIntervals
..signature:findIntervals(intervalTree, query_begin, query_end, result)
..param.query_begin:The begin position of the query interval.
..param.query_end:The end position of the query interval.
*/
template<typename TValue,typename TCargo>
void
findIntervals(IntervalTree<TValue,TCargo> & tree, TValue query_begin, TValue query_end, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	findIntervals(tree.g,tree.pm,query_begin,query_end,result);

}



template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervals(TGraph & g, TPropertyMap & pm, TValue query_begin, TValue query_end, String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	resize(result,0);

	// start at root
	TVertexDescriptor act_knot = 0;
	findIntervals(g, pm, act_knot, query_begin, query_end, result);
}


template<typename TGraph, typename TPropertyMap, typename TValue,typename TCargo>
void
findIntervals(TGraph & g, 
			  TPropertyMap & pm, 
			  typename VertexDescriptor<TGraph>::Type & act_knot, 
			  TValue query_begin, 
			  TValue query_end, 
			  String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		//
		if(act_prop.center < query_begin) // query interval is to the right of node center
		{
			unsigned int i = 0;
			while(i < length(act_prop.list2) && rightBoundary(value(act_prop.list2,i)) > query_begin)
			{
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
				++i;	
			}
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(query_end < act_prop.center) // query interval is to the left of node center
			{
				unsigned int i = 0;
				while(i < length(act_prop.list1) && leftBoundary(value(act_prop.list1,i)) < query_end)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
					++i;
				}
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{//node center is contained in query interval
				for(unsigned int i = 0; i < length(act_prop.list1); ++i)
                    appendValue(result, cargo(value(act_prop.list1,i)), Generous());
				
				while(!atEnd(it))
				{
					TVertexDescriptor next_knot = targetVertex(it);
					findIntervals(g,pm, next_knot, query_begin, query_end, result);
					goNext(it);
				}
				break;

				//break; //dont break! continue in both subtrees!!
			}
		}
	}

}



/////////////////// Metafunctions ///////////////////////

///.Metafunction.Value.param.T.type:Class.IntervalTree
template<typename TValue, typename TCargo>
struct Value<IntervalTree<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalTree
template<typename TValue, typename TCargo>
struct Cargo<IntervalTree<TValue,TCargo> >
{
	typedef TCargo Type;
};

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  //#ifndef SEQAN_HEADER_GRAPH_IMPL_INTERVALTREE_H
