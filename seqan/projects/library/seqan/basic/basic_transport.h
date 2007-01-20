#ifndef SEQAN_HEADER_BASIC_TRANSPORT_H
#define SEQAN_HEADER_BASIC_TRANSPORT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//assign
//////////////////////////////////////////////////////////////////////////////

/**
.Function.assign:
..summary:Assigns one object to another object.
..cat:Content Manipulation
..signature:assign(target, source)
..signature:assign(target, source [, limit] [,resize_tag])
..param.target: Gets the content of $source$.
..param.source: Is copied to $target$.
..param.limit: The maximal length of $target$ after the operation. (optional)
...remarks:This arguments can be applied if $target$ is a container.
..param.resize_tag: Specifies the strategy that is applied if $target$ has not enough capacity to store the complete content. (optional)
...type:Tag.Overflow Strategy
...default:Specified by @Metafunction.DefaultOverflowImplicit@ of the $target$ type.
...remarks:This arguments can be applied if $target$ is a container.
..remarks:$assign(target, source)$ is semantically equivalent to $target = source$. 
*/

template <typename TTarget, typename TSource>
inline void
assign(TTarget & target,
	   TSource & source)
{
SEQAN_CHECKPOINT
	target = source;
}
template <typename TTarget, typename TSource>
inline void
assign(TTarget & target,
	   TSource const & source)
{
SEQAN_CHECKPOINT
	target = source;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TTargetSpec, typename TSource>
inline void 
assign(Proxy<TTargetSpec> & target,
	   TSource & source)
{
SEQAN_CHECKPOINT
	assignValue(iter(target), source);
}

template<typename TTargetSpec, typename TSource>
inline void 
assign(Proxy<TTargetSpec> & target,
	   TSource const & source)
{
SEQAN_CHECKPOINT
	assignValue(iter(target), source);
}

//////////////////////////////////////////////////////////////////////////////
// move
//////////////////////////////////////////////////////////////////////////////

/**
.Function.move:
..summary:Hands over content from one container to another container.
..cat:Content Manipulation
..signature:move(target, source)
..param.target:A container $source$ is moved to.
..param.source:A container that is moved to $target$.
..remarks:The function tries to hand over the contents of $source$ to $target$.
If this is possible, $source$ losts its content and will therefore be empty after this operation.
Otherwise, the function behaves like @Function.assign@ and $source$ is copied to $target$. 
..see:Function.assign
*/

template<typename TTarget, typename TSource>
inline void 
move(TTarget & target,
	 TSource & source)
{
SEQAN_CHECKPOINT
	assign(target, source);
}
template<typename TTarget, typename TSource>
inline void 
move(TTarget const & target,
	 TSource & source)
{
SEQAN_CHECKPOINT
	assign(target, source);
}
template<typename TTarget, typename TSource>
inline void 
move(TTarget & target,
	 TSource const & source)
{
SEQAN_CHECKPOINT
	assign(target, source);
}
template<typename TTarget, typename TSource>
inline void 
move(TTarget const & target,
	 TSource const & source)
{
SEQAN_CHECKPOINT
	assign(target, source);
}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...


