#ifndef SEQAN_HEADER_GRAPH_PROPERTY_H
#define SEQAN_HEADER_GRAPH_PROPERTY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - External Property Manager
//////////////////////////////////////////////////////////////////////////////
struct EmptyMap_;
typedef Tag<EmptyMap_> const EmptyMap;


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec, typename TPropertyMap>
inline void
initVertexMap(Graph<TEdges, TSpec> const& g,
			  TPropertyMap& pm)
{
	SEQAN_CHECKPOINT
	resize(pm, getIdUpperBound(g.data_id_managerV));
}


template<typename TEdges, typename TSpec, typename TPropertyMap>
inline void
initEdgeMap(Graph<TEdges, TSpec> const& g,
			 TPropertyMap& pm)
{
	SEQAN_CHECKPOINT
	resize(pm, getIdUpperBound(g.data_id_managerE));
}


// Simple _getId function to get the id for a vertex descriptor which is the id!
template<typename TId>
inline TId
_getId(TId const id)
{
	SEQAN_CHECKPOINT
	return id;
}

template<typename TPropertyMap, typename TDescriptor, typename TValue>
inline void
assignProperty(TPropertyMap& pm,
			TDescriptor const d,
			TValue const val)
{
	SEQAN_CHECKPOINT
	assignValue(pm, _getId(d), val);
}

template<typename TPropertyMap, typename TDescriptor>
inline typename Value<TPropertyMap const>::Type
getProperty(TPropertyMap const& pm,
			TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return getValue(pm, _getId(d));
}

template<typename TPropertyMap, typename TDescriptor>
inline typename Value<TPropertyMap>::Type&
property(TPropertyMap& pm,
		TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return value(pm, _getId(d));
}

template<typename TPropertyMap, typename TDescriptor>
inline typename Value<TPropertyMap const>::Type&
property(TPropertyMap const& pm,
		TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return value(pm, _getId(d));
}




//////////////////////////////////////////////////////////////////////////////
// Graph - Internal Property Manager (only for edges!!!)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// First variant: Internal Property Manager using member Ids
//////////////////////////////////////////////////////////////////////////////
template<typename TContainer, unsigned int const MemberId = 0>
struct InternalMap 
{
};


//////////////////////////////////////////////////////////////////////////////
//	Internal Property Manager using member Ids - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 1> const> {
	typedef T1 const Type;
};

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 1> > {
	typedef T1 Type;
};

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 2> const> {
	typedef T2 const Type;
};

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 2> > {
	typedef T2 Type;
};

template<typename T>
struct Value<InternalMap<T, 0> const> {
	typedef T const Type;
};

template<typename T>
struct Value<InternalMap<T, 0> > {
	typedef T Type;
};


//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using member Ids - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec, typename TContainer, unsigned int const MemberId>
inline void
initEdgeMap(Graph<TEdges, TSpec> const& g,
			InternalMap<TContainer, MemberId>& pm)
{
	SEQAN_CHECKPOINT
}


template<typename T1, typename T2, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<Pair<T1, T2>, 1>& pm,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).i1 = val;
}

template<typename T1, typename T2, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<Pair<T1, T2>, 2>& pm,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).i2 = val;
}


template<typename T, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<T, 0>& pm,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	assignCargo(e, val);
}

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> >::Type&
property(InternalMap<Pair<T1, T2>, 2>& pm,
		TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i2;
}

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> const>::Type&
property(InternalMap<Pair<T1, T2>, 2> const& pm,
		TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i2;
}

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> >::Type&
property(InternalMap<Pair<T1, T2>, 1>& pm,
		TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i1;
}

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> const>::Type&
property(InternalMap<Pair<T1, T2>, 1> const& pm,
		TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i1;
}

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> >::Type&
property(InternalMap<T, 0>& pm,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return cargo(e);
}

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> const>::Type&
property(InternalMap<T, 0> const& pm,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return cargo(e);
}

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> const>::Type
getProperty(InternalMap<Pair<T1, T2>, 1> const& pm,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i1;
}

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> >::Type
getProperty(InternalMap<Pair<T1, T2>, 1>& pm,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i1;
}

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> const>::Type
getProperty(InternalMap<Pair<T1, T2>, 2> const& pm,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i2;
}

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> >::Type
getProperty(InternalMap<Pair<T1, T2>, 2>& pm,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i2;
}

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> const>::Type
getProperty(InternalMap<T, 0> const& pm,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return getCargo(e);
}

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> >::Type
getProperty(InternalMap<T, 0>& pm,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return getCargo(e);
}



//////////////////////////////////////////////////////////////////////////////
// Second variant: Internal Property Manager using pointer to members
//////////////////////////////////////////////////////////////////////////////
template <typename TPropmap, TPropmap const Instance> 
struct InternalPointerMap 
{
}; 


//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using pointer to members - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TClass, typename TValue, TValue TClass:: * TPMember>
struct Value<InternalPointerMap<TValue TClass::*, TPMember> const> {
	typedef TValue const Type;
};

template<typename TClass, typename TValue, TValue TClass:: * TPMember>
struct Value<InternalPointerMap<TValue TClass::*, TPMember> > {
	typedef TValue Type;
};



//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using pointer to members - Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec, typename TPropmap, TPropmap const Instance>
inline void
initEdgeMap(Graph<TEdges, TSpec> const& g,
			InternalPointerMap<TPropmap, Instance>& pm)
{
	SEQAN_CHECKPOINT
}


template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline void
assignProperty(InternalPointerMap<TValue TClass::*, TPMember>& pm,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).*TPMember = val;
}

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type&
property(InternalPointerMap<TValue TClass::*, TPMember>& pm,
		TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*TPMember;
}

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> const>::Type&
property(InternalPointerMap<TValue TClass::*, TPMember> const& pm,
		TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*TPMember;
}

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> const>::Type
getProperty(InternalPointerMap<TValue TClass::*, TPMember> const& pm,
			TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*TPMember;
}

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type
getProperty(InternalPointerMap<TValue TClass::*, TPMember>& pm,
			TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*TPMember;
}


//////////////////////////////////////////////////////////////////////////////
// Third variant: Raw pointer to member
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Raw pointer to member - Metafunctions
//////////////////////////////////////////////////////////////////////////////
template <typename TClass, typename TValue> 
struct Value<TValue TClass:: *> {
	typedef TValue Type;
};

template <typename TClass, typename TValue> 
struct Value<TValue TClass:: * const> {
	typedef TValue const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Raw pointer to member - Functions
//////////////////////////////////////////////////////////////////////////////
template <typename TEdges, typename TSpec, typename TClass, typename TValue> 
inline void
initEdgeMap(Graph<TEdges, TSpec> const& g,
			 TValue TClass:: * ptr_to_member)
{
}

template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline void 
assignProperty(TValue TClass:: * ptr_to_member, 
			TEdgeDescriptor const e, 
			TValue const val) 
{
	SEQAN_CHECKPOINT
	(cargo(e)).*ptr_to_member=val; 
}


template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline TValue& 
property(TValue TClass:: * const ptr_to_member, 
			TEdgeDescriptor const e) 
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*ptr_to_member; 
}

template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline TValue
getProperty(TValue TClass:: * const ptr_to_member, 
			TEdgeDescriptor const e) 
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*ptr_to_member; 
} 


//////////////////////////////////////////////////////////////////////////////
// Generic functions
//////////////////////////////////////////////////////////////////////////////
template<typename TEdges, typename TSpec, typename TPropertyMap, typename TProperties>
inline void
initVertexMap(Graph<TEdges, TSpec> const& g,
			  TPropertyMap& pm,
			  TProperties const& prop)
{
	SEQAN_CHECKPOINT
	initVertexMap(g,pm);
	unsigned int count = 0;
	for(unsigned int i=0;i<length(pm);++i) {
		if (idInUse(g.data_id_managerV,i)) {
			assignProperty(pm,i,prop[count]);
			++count;
		}
	}
}

template<typename TEdges, typename TSpec, typename TPropertyMap, typename TProperties>
inline void
initEdgeMap(Graph<TEdges, TSpec>& g,
			  TPropertyMap& pm,
			  TProperties const& prop)
{
	SEQAN_CHECKPOINT
	initEdgeMap(g,pm);
	typedef typename Id<Graph<TEdges, TSpec> >::Type TIdType;
	typedef typename VertexDescriptor<Graph<TEdges, TSpec> >::Type TVertexDescriptor;
	typedef typename EdgeType<Graph<TEdges, TSpec> >::Type TEdgeStump;
	for(TIdType id = getIdLowerBound(g.data_id_managerV);id<getIdUpperBound(g.data_id_managerV);++id) {
		if (!idInUse(g.data_id_managerV, id)) continue;
		TEdgeStump* current = getValue(g.data_vertex, id);
		while(current != (TEdgeStump*) 0) {
			assignProperty(pm,current,prop[_getId(current)]);
			current = current->data_next;
		}
	}
}
/*
template<typename TContainer, unsigned int const MemberId>
struct MemberType
{

};

template<typename T1, typename T2>
struct MemberType<Pair<T1, T2>, 1> {
	typedef T1 Pair<T1, T2>:: * Type;
};

template<typename T1, typename T2>
struct MemberType<Pair<T1, T2>, 2> {
	typedef T2 Pair<T1, T2>:: * Type;
};
*/



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
