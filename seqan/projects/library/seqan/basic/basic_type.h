#ifndef SEQAN_HEADER_BASIC_TYPE_H
#define SEQAN_HEADER_BASIC_TYPE_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Value:
..summary:Type of the items in the container. 
..signature:Value<T>::Type
..param.T:Type for which the value type is determined.
..returns.param.Type:Value type of $T$.
..remarks.text:The value type of a container $T$ is the type of the elements in $T$.
    For example, the value type of a sequence of $int$ is $int$.
..example.code:Value<String<char> >::Type c; //c has type char
*/

template <typename T, const int i = 0>
struct Value
{
	typedef T Type;
};
/*
template <typename T>
struct Value<T const>
{
	typedef T Type;
};
*/
//____________________________________________________________________________

/**
.Metafunction.GetValue:
..summary:Type for reading values. 
..signature:GetValue<T>::Type
..param.T:Type of container that holds a value.
..returns.param.Type:GetValue type of $T$.
..remarks.text:Depending on $T$, the $GetValue$-type can either be $Value<T>::Type &$ or $Value<T>::Type$.
..text:$GetValue$ is the return type of @Function.getValue@ that allows a (read-only) access to objects.
Do not confuse it with @Function.value@ that returns a @Metafunction.Reference.reference@ to the value.
..see:Metafunction.Value
..see:Function.getValue
*/
template <typename T>
struct GetValue
{
	typedef typename Value<T>::Type const & Type;
};
template <typename T>
struct GetValue<T const>:
	public GetValue<T>
{
};

//____________________________________________________________________________

/**
.Metafunction.Reference:
..summary:Reference type. 
..signature:Reference<T>::Type
..param.T:A Type.
..returns.param.Type:Either $T &$ or a proxy object @Class.Proxy@ for $T$.
..see:Metafunction.Value
..see:Metafunction.GetValue
*/
template <typename T>
struct Reference
{
	typedef typename Value<T>::Type & Type;
};
template <typename T>
struct Reference<T const>
{
	typedef typename Value<T>::Type const & Type;
};

//____________________________________________________________________________


/**
.Metafunction.Size:
..summary:Type of an object that is suitable to hold size information.
..signature:Size<T>::Type
..param.T:Type for which the size type is determined.
..returns.param.Type:Size type of $T$.
..remarks.text:In most cases this type is $size_t$.
*/
template <typename T>
struct Size
{
	typedef size_t Type;
};
template <typename T>
struct Size<T const>:
	Size<T>
{
};

//____________________________________________________________________________


/**
.Metafunction.Difference:
..summary:Type of an object that stores the difference between two iterators.
..signature:Difference<T>::Type
..param.T:Type for which the difference type is determined.
...type:Class.Iter
..returns.param.Type:Difference type of $T$.
..remarks.text:In most cases this type is $ptrdiff_t$.
..see:Metafunction.Size
*/
template <typename T>
struct Difference
{
	typedef ptrdiff_t Type;
};
template <typename T>
struct Difference<T const>:
	Difference<T>
{
};

//____________________________________________________________________________


/**
.Metafunction.Position:
..summary:Type of an object that represents a position in a container.
..signature:Position<T>::Type
..param.T:Type for which the position type is determined.
...type:Class.Iter
...type:Class.String
..returns.param.Type:Position type of $T$.
..see:Metafunction.Iterator
*/
template <typename T>
struct Position
{
	typedef typename Size<T>::Type Type;
};
template <typename T>
struct Position<T const>:
	Position<T>
{
};

//____________________________________________________________________________

/**
.Metafunction.Host:
..summary:Type of the object a given object depends on.
..signature:Host<T>::Type
..param.T:Type for which the host type is determined.
..returns.param.Type:Host type of $T$.
*/
template <typename T>
struct Host
{
	typedef T Type;
};

//____________________________________________________________________________

/**
.Metafunction.Spec:
..summary:The spec of a class. 
..signature:Spec<T>::Type
..param.T:Type for which the spec is determined.
..returns.param.Type:Spec of $T$.
..remarks:The spec of a SeqAn type is the class that is used in template subclassing 
 to specify the specialization. 
 For example, spec of $String<char, Alloc<> >$ is $Alloc<>$.
*/

template <typename T>
struct Spec;

//____________________________________________________________________________

/**
.Metafunction.Cargo:
..summary:Additional data of a class. 
..signature:Cargo<T>::Type
..param.T:Type for which the cargo is determined.
..returns.param.Type:Cargo of $T$.
..remarks:The definition of Cargo allows the addition of user specific data to existing data structures.
*/

template <typename T>
struct Cargo {
	typedef Nothing Type;
};
template <typename T>
struct Cargo<T const> {
	typedef typename Cargo<T>::Type const Type;
};

//____________________________________________________________________________

/**
.Metafunction.VertexDescriptor:
..summary:Type of an object that represents a vertex descriptor.
..signature:VertexDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs currently use ids as vertex descriptors.
..returns.param.Type:VertexDescriptor type.
..remarks.text:The vertex descriptor is a unique handle to a vertex in a graph.
It is used in various graph functions, e.g., to add edges, to create OutEdge Iterators or to remove a vertex.
It is also used to attach properties to vertices.
..example.code:VertexDescriptor<Graph<> >::Type vD; //vD is a vertex descriptor
*/

template <typename T>
struct VertexDescriptor {
	typedef void* Type;
};
template <typename T>
struct VertexDescriptor<T const>:
	public VertexDescriptor<T> {};


//____________________________________________________________________________

/**
.Internal._Parameter:
..cat:Metafunctions
..summary:Type for function parameters and return values.
..signature:_Parameter<T>::Type
..param.T:A type.
..returns.param.Type:The parameter type for arguments of type $T$.
...text:If $T$ is a pointer or array type, then $_Parameter<T>::Type$ is $T$, 
otherwise $_Parameter<T>::Type$ is $T &$.
*/
template <typename T>
struct _Parameter
{
	typedef T & Type;
};

template <typename T>
struct _Parameter<T *>
{
	typedef T * Type;
};
template <typename T, size_t I>
struct _Parameter<T [I]>
{
	typedef T * Type;
};


/**
.Internal._toParameter:
..cat:Functions
..summary:Transforms pointers to parameter types.
..signature:_toParameter<T>(pointer)
..param.pointer:A pointer.
..param.T:A Type.
...text:$object$ is transformed into the parameter type of $T$ that is given by @Internal._Parameter@.
...note:This type must be explicitely specified.
..returns:To $TParameter$ transformed $object$.
..see:Internal._Parameter
*/
template <typename T>
typename _Parameter<T>::Type
_toParameter(T * _object)
{
SEQAN_CHECKPOINT
	return * _object;
}
template <typename T>
typename _Parameter<T>::Type
_toParameter(T _object)
{
SEQAN_CHECKPOINT
	return _object;
}

//____________________________________________________________________________

/**
.Internal._ConstParameter:
..cat:Metafunctions
..summary:Type for constant function parameters and return values.
..signature:_ConstParameter<T>::Type
..param.T:A type.
..returns.param.Type:The const parameter type for arguments of type $T$.
...text:If $T$ is a pointer or array type, then $_Parameter<T>::Type$ is a pointer to a const array, 
otherwise $_Parameter<T>::Type$ is $T const &$.
..see:Internal._Parameter
*/
template <typename T>
struct _ConstParameter
{
	typedef T const & Type;
};
template <typename T>
struct _ConstParameter<T const>:
	public _ConstParameter<T> {};

template <typename T>
struct _ConstParameter<T *>
{
	typedef T const * Type;
};
template <typename T>
struct _ConstParameter<T const *>
{
	typedef T const * Type;
};

template <typename T, size_t I>
struct _ConstParameter<T [I]>
{
	typedef T const * Type;
};
template <typename T, size_t I>
struct _ConstParameter<T const [I]>
{
	typedef T const * Type;
};

//____________________________________________________________________________

/**
.Internal._Pointer:
..cat:Metafunctions
..summary:The associated pointer type.
..signature:_Pointer<T>::Type
..param.T:A type.
..returns.param.Type:A pointer type for $T$.
...text:if $T$ is already a pointer type, then $_Pointer<T>::Type$ is $T$,
otherwise $_Pointer<T>::Type$ is $T *$.
..see:Internal._Parameter
..see:Internal._toParameter
*/
template <typename T>
struct _Pointer
{
	typedef T * Type;
};

template <typename T>
struct _Pointer<T *>
{
	typedef T * Type;
};
template <typename T>
struct _Pointer<T * const>
{
	typedef T * const Type;
};

template <typename T, size_t I>
struct _Pointer<T [I]>
{
	typedef T * Type;
};

/**
.Internal._toPointer:
..cat:Functions
..summary:Transforms types into pointers.
..signature:_toPointer(object)
..param.object:An object.
..returns:$object$, transformed to a pointer. 
...text:The type of the returned pointer is given by @Internal._Pointer@.
..see:Internal._Pointer
*/
template <typename T>
typename _Pointer<T>::Type
_toPointer(T & _object)
{
SEQAN_CHECKPOINT
	return & _object;
}
template <typename T>
typename _Pointer<T const>::Type
_toPointer(T const & _object)
{
SEQAN_CHECKPOINT
	return & _object;
}

template <typename T>
typename _Pointer<T *>::Type
_toPointer(T * _object)
{
SEQAN_CHECKPOINT
	return _object;
}

//////////////////////////////////////////////////////////////////////////////

//Iterator: see basic_iterator.h

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
