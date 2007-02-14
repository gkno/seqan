#ifndef SEQAN_HEADER_GRAPH_STACK_H
#define SEQAN_HEADER_GRAPH_STACK_H

namespace SEQAN_NAMESPACE_MAIN {

//////////////////////////////////////////////////////////////////////////////
// TAGS
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Block String:
..cat:Strings
..general:Class.String
..summary:String optimized for push_back, top, and pop (Stack behaviour).
..signature:String<TValue, BlockString<size> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.size:A positive integer that specifies the number of values in each allocated block.
...remarks: Size should be a power of 2, e.g., 1024.
*/
template<unsigned int SPACE = 1024>
struct BlockString;

//////////////////////////////////////////////////////////////////////////////


template<typename TValue, unsigned int SPACE>
class String<TValue, BlockString<SPACE> > 
{
	public:
		String<TValue*> blockBegins;
		TValue* lastValue;
    
//____________________________________________________________________________
      
    public:
		String() : lastValue(0) 
		{
			SEQAN_CHECKPOINT
		}

		template<typename TSource>
		String(TSource const& source) : lastValue(0) 
		{
			SEQAN_CHECKPOINT
			assign(*this,source);
		} 

		String(String const & source) : lastValue(0) 
		{
			SEQAN_CHECKPOINT
			assign(*this, source);
		}

		template<typename TSource>
		String & operator =(TSource const& source) 
		{
			SEQAN_CHECKPOINT
			// Clean-up old data
			while (lastValue != 0) {
				pop_back(*this);
			} 
			assign(*this, source);
			return *this;
		}

		String & operator =(String const& _other)	
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			// Clean-up old data
			while (lastValue != 0) {
				pop_back(*this);
			} 
			// Assign new data
			assign(*this, _other);
			return *this;
		}

		~String() 
		{
			while (lastValue != 0) {
				pop_back(*this);
			}
		}

//____________________________________________________________________________

	public:
		template<typename TPos>
		inline typename Reference<String>::Type 
			operator[] (TPos pos) 
		{
			SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template<typename TPos>
		inline typename Reference<String const>::Type 
			operator[] (TPos pos) const 
		{
			SEQAN_CHECKPOINT
			return value(*this, pos);
		}

};



//////////////////////////////////////////////////////////////////////////////
// BlockString metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Iterators
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.BlockString
template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
struct Iterator<String<TValue, BlockString<SPACE> >, TIteratorSpec> 
{
	typedef Iter<String<TValue, BlockString<SPACE> >, PositionIterator> Type;
};

template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
struct Iterator<String<TValue, BlockString<SPACE> > const, TIteratorSpec> 
{
	typedef Iter<String<TValue, BlockString<SPACE> > const, PositionIterator> Type;
};


///////////////////////////////////////////////////////////////
// BlockString interface
///////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// begin
//////////////////////////////////////////////////////////////////////////////
/**
.Function.begin:
..cat:Iteration
..cat:String
..summary:The begin of a container. 
..signature:begin(object [, tag])
..param.object:A container.
...type:Class.String
..param.tag:An @Tag.Iterator Spec.iterator spec@ tag that specifies the kind of the iterator returned. (optional)
...default:Given by @Metafunction.DefaultIteratorSpec@.
..returns:An iterator to the first item in $object$. 
$Iterator<T>::Type$ is the type of the iterator for $object$-type $T$.
..remarks.text:If the container does not contain any items at all, the function may return 0.
..see:Function.end
..see:Metafunction.Iterator
*/
template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
inline typename Iterator<String<TValue, BlockString<SPACE> >, Tag<TIteratorSpec> const>::Type 
begin(
	  String<TValue, BlockString<SPACE> >& me, 
	  Tag<TIteratorSpec> const) 
{
	SEQAN_CHECKPOINT
	return Iter<String<TValue, BlockString<SPACE> >, PositionIterator>(me,0);
}

template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
inline typename Iterator<String<TValue, BlockString<SPACE> > const, Tag<TIteratorSpec> const>::Type 
begin(
	  String<TValue, BlockString<SPACE> > const& me, 
	  Tag<TIteratorSpec> const) 
{
	SEQAN_CHECKPOINT
	return Iter<String<TValue, BlockString<SPACE> > const, PositionIterator>(me,0);
}


template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
inline typename Iterator<String<TValue, BlockString<SPACE> >, Tag<TIteratorSpec> const>::Type 
end(
	String<TValue, BlockString<SPACE> >& me, 
	Tag<TIteratorSpec> const)
{
	SEQAN_CHECKPOINT
    if (me.lastValue==0) return Iter<String<TValue, BlockString<SPACE> >, PositionIterator>(me,0);
    else {
		unsigned int endPos = (length(me.blockBegins)-1) * SPACE;
		TValue* endOfBlock = me.blockBegins[length(me.blockBegins)-1];
		while (endOfBlock != me.lastValue) {
			++endOfBlock;
			++endPos;
		}
		return Iter<String<TValue, BlockString<SPACE> >, PositionIterator>(me,endPos);
    }
}

template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
inline typename Iterator<String<TValue, BlockString<SPACE> > const, Tag<TIteratorSpec> const>::Type 
end(
	String<TValue, BlockString<SPACE> > const& me, 
	Tag<TIteratorSpec> const)
{
	SEQAN_CHECKPOINT
    if (me.lastValue==0) return Iter<String<TValue, BlockString<SPACE> > const, PositionIterator>(me,0);
    else {
		unsigned int endPos = (length(me.blockBegins)-1) * SPACE;
		TValue* endOfBlock = me.blockBegins[length(me.blockBegins)-1];
		while (endOfBlock != me.lastValue) {
			++endOfBlock;
			++endPos;
		}
		return Iter<String<TValue, BlockString<SPACE> > const, PositionIterator>(me,endPos);
    }
}

template<typename TValue, unsigned int SPACE, typename TSource>
inline void 
assign(
	   String<TValue, BlockString<SPACE> >& target, 
	   TSource const& source) 
{
    SEQAN_CHECKPOINT
    typedef typename Iterator<TSource const, Standard>::Type TIter;
    for(TIter it = begin(source, Standard()); !atEnd(it, source); goNext(it)) {
		push_back(target, *it);
    }
}


template<typename TValue, unsigned int SPACE, typename TPos>
inline typename Reference<String<TValue, BlockString<SPACE> > >::Type 
value(
	  String<TValue, BlockString<SPACE> >& stack, 
	  TPos const pos) 
{
	SEQAN_CHECKPOINT
    return *(stack.blockBegins[(pos / SPACE)] + (pos % SPACE));
}

template<typename TValue, unsigned int SPACE, typename TPos>
inline typename Reference<String<TValue, BlockString<SPACE> > >::Type 
value(
	  String<TValue, BlockString<SPACE> > const& stack, 
	  TPos const pos) 
{
	SEQAN_CHECKPOINT
    return *(stack.blockBegins[(pos / SPACE)] + (pos % SPACE));
}

template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
inline bool 
atEnd(
	  Iter<String<TValue, BlockString<SPACE> >, TIteratorSpec>& it, 
	  String<TValue, BlockString<SPACE> >& container) 
{
    SEQAN_CHECKPOINT
    typedef typename Iterator<String<TValue, BlockString<SPACE> >, TIteratorSpec>::Type TIter;
    TIter endIt = end(container, Standard());
    return (it == endIt);
}


template<typename TValue, unsigned int SPACE>
inline void 
clear(String<TValue, BlockString<SPACE> >& me)
{
	SEQAN_CHECKPOINT
	while (me.lastValue != 0) {
		pop_back(me);
	}
}

template<typename TValue, unsigned int SPACE, typename TSource, typename TExpand>
inline void 
append(String<TValue, BlockString<SPACE> >& me,
	   TSource const& source,
	   Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<TSource const, Standard>::Type TIter;
    for(TIter it = begin(source, Standard()); !atEnd(it, source); goNext(it)) {
		push_back(me, *it);
    }
}

template<typename TValue, unsigned int SPACE, typename TVal, typename TExpand>
inline void 
appendValue(String<TValue, BlockString<SPACE> >& me,
	   TVal const& source,
	   Tag<TExpand> const tag)
{
SEQAN_CHECKPOINT
	typedef typename Iterator<TVal const, Standard>::Type TIter;
    for(TIter it = begin(source, Standard()); !atEnd(it, source); goNext(it)) {
		push_back(me, *it);
    }
}

template<typename TValue, unsigned int SPACE>
inline void 
push_back(
		  String<TValue, BlockString<SPACE> >& stack, 
		  TValue const& element) 
{
	SEQAN_CHECKPOINT
    
    //Do we need to allocate a new block
    if((stack.lastValue == 0) || 
       (stack.lastValue == (stack.blockBegins[length(stack.blockBegins)-1]) + SPACE)) {
		TValue* mem=0;
		allocate(stack, mem, sizeof(TValue)*SPACE);
		arrayFill(mem, mem+sizeof(TValue)*SPACE, 0);
		appendValue(stack.blockBegins, mem);
		stack.lastValue=mem;
	}
      
    valueConstruct(stack.lastValue);
    *(stack.lastValue) = element;
    ++stack.lastValue;
}
	 

template<typename TValue, unsigned int SPACE>
inline TValue 
top(String<TValue, BlockString<SPACE> > const& stack) 
{
    SEQAN_CHECKPOINT
    if (stack.lastValue == 0) {
		static TValue dummy;
		return dummy;
    }
    return *(stack.lastValue - 1);
}

template<typename TValue, unsigned int SPACE>
inline void 
pop_back(String<TValue, BlockString<SPACE> >& stack) 
{
    SEQAN_CHECKPOINT
    if (stack.lastValue==0) return;
    --stack.lastValue;
    valueDestruct(stack.lastValue);
    
    if (length(stack.blockBegins)<=1) {
		if (stack.lastValue == stack.blockBegins[0]) {
      		deallocate(stack, stack.blockBegins[0], stack.blockBegins[0] + sizeof(TValue)*SPACE);
			erase(stack.blockBegins, length(stack.blockBegins)-1);
			stack.lastValue=0;
		}
    } else if(stack.lastValue == (stack.blockBegins[length(stack.blockBegins)-1])) {
		deallocate(stack, stack.blockBegins[length(stack.blockBegins)-1], stack.blockBegins[length(stack.blockBegins)-1] + sizeof(TValue)*SPACE);
		erase(stack.blockBegins, length(stack.blockBegins)-1);
		stack.lastValue = stack.blockBegins[length(stack.blockBegins)-1] + SPACE;
    }
    return;
}


template<typename TValue, unsigned int SPACE>
inline bool 
empty(String<TValue, BlockString<SPACE> > const& stack) 
{
    SEQAN_CHECKPOINT
    return (stack.lastValue==0);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
