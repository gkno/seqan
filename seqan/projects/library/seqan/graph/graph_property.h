#ifndef SEQAN_HEADER_GRAPH_PROPERTY_H
#define SEQAN_HEADER_GRAPH_PROPERTY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//	Graph - External property map
//////////////////////////////////////////////////////////////////////////////

/**
.Class.ExternalMap:
..cat:Graph
..summary:The external property map is assumed to be an instance of @Class.String@.
It is indexed via VertexDescriptors or EdgeDescriptors.
..signature:String<TValue, TSpec>
..param.TValue:The value type. That is the type of information stored in the property map.
...metafunction:Metafunction.Value
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Alloc<>$, see @Spec.Alloc String@.
..include:graph.h
*/


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.initVertexMap:
..cat:Graph
..summary:Initializes a vertex map. 
..signature:initVertexMap(g, pm)
..signature:initVertexMap(g, pm, prop)
..param.g:A Graph.
...type:Class.Graph
..param.pm:An External Property Map.
...type:Class.String
..param.prop:An optional array with properties that are to be assigned to the items in the property map.
...remarks:For every vertex descriptor there must be an entry in the array.
..returns:void
..see:Function.initEdgeMap
*/

///.Function.initVertexMap.param.pm.type:Class.ExternalMap

template<typename TSpec, typename TPropertyMap>
inline void
initVertexMap(Graph<TSpec> const& g,
			  TPropertyMap& pm)
{
	SEQAN_CHECKPOINT
	resize(pm, getIdUpperBound(g.data_id_managerV));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.initEdgeMap:
..cat:Graph
..summary:Initializes an edge map
..signature:initEdgeMap(g, pm)
..signature:initEdgeMap(g, pm, prop)
..param.g:A Graph.
...type:Class.Graph
..param.pm:An External Property Map.
...type:Class.String
..param.prop:An optional array with properties that are to be assigned to the items in the property map.
...remarks:For every edge id there must be an entry in the array.
..returns:void
..see:Function.initVertexMap
*/

///.Function.initEdgeMap.param.pm.type:Class.ExternalMap

template<typename TSpec, typename TPropertyMap>
inline void
initEdgeMap(Graph<TSpec> const& g,
			  TPropertyMap& pm)
{
	SEQAN_CHECKPOINT
	resize(pm, getIdUpperBound(g.data_id_managerE), Generous());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TPropertyMap>
inline void
initEdgeMap(Graph<Tree<TCargo, TSpec> > const& g,
			TPropertyMap& pm)
{
	SEQAN_CHECKPOINT
	resize(pm, getIdUpperBound(g.data_id_managerV), Generous());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TPropertyMap>
inline void
initEdgeMap(Graph<Tree<TCargo, TSpec> >& g,
			TPropertyMap& pm)
{
	SEQAN_CHECKPOINT
	resize(pm, getIdUpperBound(g.data_id_managerV), Generous());
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignProperty:
..cat:Graph
..summary:Assigns a property to an item in the property map.
..signature:assignProperty(pm, d, val)
..param.pm:An External Property Map.
...type:Class.String
..param.d:A vertex or edge descriptor.
...remarks:Identifies the item in the property map.
..param.val:The new value.
...remarks:Type of the new value must match the value type of the property map.
..returns:void
..see:Function.getProperty
..see:Function.property
*/

///.Function.assignProperty.param.pm.type:Class.ExternalMap

template<typename TPropertyMap, typename TDescriptor, typename TValue>
inline void
assignProperty(TPropertyMap& pm,
			TDescriptor const d,
			TValue const val)
{
	SEQAN_CHECKPOINT
	assignValue(pm, _getId(d), val);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.property:
..cat:Graph
..summary:Accesses the property of an item in the property map.
..signature:assignProperty(pm, d)
..param.pm:An External Property Map.
...type:Class.String
..param.d:A vertex or edge descriptor.
...remarks:Identifies the item in the property map.
..returns:Reference to the item in the property map.
..see:Function.getProperty
..see:Function.assignProperty
*/

///.Function.property.param.pm.type:Class.ExternalMap

template<typename TPropertyMap, typename TDescriptor>
inline typename Reference<TPropertyMap>::Type
property(TPropertyMap& pm,
		TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return value(pm, _getId(d));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TPropertyMap, typename TDescriptor>
inline typename Reference<TPropertyMap const>::Type
property(TPropertyMap const& pm,
		TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return value(pm, _getId(d));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getProperty:
..cat:Graph
..summary:Get method for an item's property.
..signature:getProperty(pm, d)
..param.pm:An External Property Map.
...type:Class.String
..param.d:A vertex or edge descriptor.
...remarks:Identifies the item in the property map.
..returns:Value of the item in the property map.
..see:Function.property
..see:Function.assignProperty
*/

///.Function.getProperty.param.pm.type:Class.ExternalMap

template<typename TPropertyMap, typename TDescriptor>
inline typename GetValue<TPropertyMap const>::Type
getProperty(TPropertyMap const& pm,
			TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return getValue(pm, _getId(d));
}






//////////////////////////////////////////////////////////////////////////////
// Graph - Internal Property Manager using member ids (only for edges!!!)
//////////////////////////////////////////////////////////////////////////////

/**
.Class.InternalMap:
..cat:Graph
..summary:An internal property map using member ids.
Internal property maps are used to access internal edge cargos.
..signature:InternalMap<TContainer, MemberId>
..param.TContainer:The cargo type.
...metafunction:Metafunction.Cargo
..param.MemberId:An unsigned int.
...remarks:Specifies the position of the member in the cargo.
Note: If zero it is assumed that the cargo is a simple type (e.g., int).
...default:$0$.
..include:graph.h
*/
template<typename TContainer, unsigned int const MemberId = 0>
struct InternalMap 
{
};


//////////////////////////////////////////////////////////////////////////////
//	Internal Property Manager using member Ids - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.InternalMap

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 1> const> {
	typedef T1 const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 1> > {
	typedef T1 Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 2> const> {
	typedef T2 const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 2> > {
	typedef T2 Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct Value<InternalMap<T, 0> const> {
	typedef T const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct Value<InternalMap<T, 0> > {
	typedef T Type;
};


//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using member Ids - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.initEdgeMap:
..cat:Graph
..summary:Initializes an edge map
..signature:initEdgeMap(g, pm)
..signature:initEdgeMap(g, pm, prop)
..param.g:A Graph.
...type:Class.Graph
..param.pm:An Internal Property Map.
...type:Class.InternalMap
..param.prop:An optional array with properties that are to be assigned to the items in the property map.
...remarks:For every edge id there must be an entry in the array.
..returns:void
*/

///.Function.initEdgeMap.param.pm.type:Class.InternalMap

template<typename TSpec, typename TContainer, unsigned int const MemberId>
inline void
initEdgeMap(Graph<TSpec> const& g,
			InternalMap<TContainer, MemberId>& pm)
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TContainer, unsigned int const MemberId>
inline void
initEdgeMap(Graph<TSpec>& g,
			InternalMap<TContainer, MemberId>& pm)
{
}

/**
.Function.assignProperty:
..cat:Graph
..summary:Assigns a property to an item in the property map.
..signature:assignProperty(pm, e, val)
..param.pm:An Internal Property Map.
...type:Class.InternalMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..param.val:The new value.
...remarks:Type of the new value must match the value type of the property map.
..returns:void
..see:Function.getProperty
..see:Function.property
*/

///.Function.assignProperty.param.pm.type:Class.InternalMap

template<typename T1, typename T2, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<Pair<T1, T2>, 1>& pm,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).i1 = val;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<Pair<T1, T2>, 2>& pm,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).i2 = val;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<T, 0>& pm,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	assignCargo(e, val);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.property:
..cat:Graph
..summary:Accesses the property to an item in the property map.
..signature:property(pm, e)
..param.pm:An Internal Property Map.
...type:Class.InternalMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..returns:Reference to the property.
..see:Function.getProperty
..see:Function.assignProperty
*/

///.Function.property.param.pm.type:Class.InternalMap

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> >::Type&
property(InternalMap<Pair<T1, T2>, 2>& pm,
		TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> const>::Type&
property(InternalMap<Pair<T1, T2>, 2> const& pm,
		TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> >::Type&
property(InternalMap<Pair<T1, T2>, 1>& pm,
		TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> const>::Type&
property(InternalMap<Pair<T1, T2>, 1> const& pm,
		TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> >::Type&
property(InternalMap<T, 0>& pm,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return cargo(e);
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> const>::Type&
property(InternalMap<T, 0> const& pm,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return cargo(e);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getProperty:
..cat:Graph
..summary:Get method for an item's property.
..signature:getProperty(pm, e)
..param.pm:An Internal Property Map.
...type:Class.InternalMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..returns:The value of that item's property.
..see:Function.property
..see:Function.assignProperty
*/

///.Function.property.param.pm.type:Class.InternalMap

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> const>::Type
getProperty(InternalMap<Pair<T1, T2>, 1> const& pm,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> >::Type
getProperty(InternalMap<Pair<T1, T2>, 1>& pm,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> const>::Type
getProperty(InternalMap<Pair<T1, T2>, 2> const& pm,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> >::Type
getProperty(InternalMap<Pair<T1, T2>, 2>& pm,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> const>::Type
getProperty(InternalMap<T, 0> const& pm,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return getCargo(e);
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> >::Type
getProperty(InternalMap<T, 0>& pm,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return getCargo(e);
}


//////////////////////////////////////////////////////////////////////////////
// Graph - Internal Property Manager using pointer to members (only for edges!!!)
//////////////////////////////////////////////////////////////////////////////

/**
.Class.InternalPointerMap:
..cat:Graph
..summary:An internal property map using pointer to members.
Internal property maps are used to access internal edge cargos.
..signature:InternalPointerMap<TPropmap, Instance>
..param.TPropmap:A pointer to member type.
..param.Instance:A pointer to a member of type TPropmap.
..include:graph.h
*/
template <typename TPropmap, TPropmap const Instance> 
struct InternalPointerMap 
{
}; 

//////////////////////////////////////////////////////////////////////////////
//	Internal Property Manager using member Ids - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.InternalPointerMap

template<typename TClass, typename TValue, TValue TClass:: * TPMember>
struct Value<InternalPointerMap<TValue TClass::*, TPMember> const> {
	typedef TValue const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TClass, typename TValue, TValue TClass:: * TPMember>
struct Value<InternalPointerMap<TValue TClass::*, TPMember> > {
	typedef TValue Type;
};


//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using pointer to members - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.initEdgeMap:
..cat:Graph
..summary:Initializes an edge map
..signature:initEdgeMap(g, pm)
..signature:initEdgeMap(g, pm, prop)
..param.g:A Graph.
...type:Class.Graph
..param.pm:An Internal Property Map.
...type:Class.InternalPointerMap
..param.prop:An optional array with properties that are to be assigned to the items in the property map.
...remarks:For every edge id there must be an entry in the array.
..returns:void
*/

///.Function.initEdgeMap.param.pm.type:Class.InternalPointerMap

template<typename TSpec, typename TPropmap, TPropmap const Instance>
inline void
initEdgeMap(Graph<TSpec>& g,
			InternalPointerMap<TPropmap, Instance>& pm)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPropmap, TPropmap const Instance>
inline void
initEdgeMap(Graph<TSpec> const& g,
			InternalPointerMap<TPropmap, Instance>& pm)
{
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignProperty:
..cat:Graph
..summary:Assigns a property to an item in the property map.
..signature:assignProperty(pm, e, val)
..param.pm:An Internal Property Map.
...type:Class.InternalPointerMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..param.val:The new value.
...remarks:Type of the new value must match the value type of the property map.
..returns:void
..see:Function.getProperty
..see:Function.property
*/

///.Function.assignProperty.param.pm.type:Class.InternalPointerMap

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline void
assignProperty(InternalPointerMap<TValue TClass::*, TPMember>& pm,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).*TPMember = val;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.property:
..cat:Graph
..summary:Accesses the property to an item in the property map.
..signature:property(pm, e)
..param.pm:An Internal Property Map.
...type:Class.InternalPointerMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..returns:Reference to the property.
..see:Function.getProperty
..see:Function.assignProperty
*/

///.Function.property.param.pm.type:Class.InternalPointerMap

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type&
property(InternalPointerMap<TValue TClass::*, TPMember>& pm,
		TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*TPMember;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> const>::Type&
property(InternalPointerMap<TValue TClass::*, TPMember> const& pm,
		TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*TPMember;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getProperty:
..cat:Graph
..summary:Get method for an item's property.
..signature:getProperty(pm, e)
..param.pm:An Internal Property Map.
...type:Class.InternalPointerMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..returns:The value of that item's property.
..see:Function.property
..see:Function.assignProperty
*/

///.Function.getProperty.param.pm.type:Class.InternalPointerMap

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> const>::Type
getProperty(InternalPointerMap<TValue TClass::*, TPMember> const& pm,
			TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*TPMember;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type
getProperty(InternalPointerMap<TValue TClass::*, TPMember>& pm,
			TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*TPMember;
}



//////////////////////////////////////////////////////////////////////////////
// Graph - Internal Property Manager using raw pointer to member (only for edges!!!)
//////////////////////////////////////////////////////////////////////////////

/**
.Class.InternalRawMap:
..cat:Graph
..summary:An internal property map using raw pointer to members.
Internal property maps are used to access internal edge cargos.
..include:graph.h
*/

//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using raw pointer to member - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.InternalRawMap

template <typename TClass, typename TValue> 
struct Value<TValue TClass:: *> {
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TClass, typename TValue> 
struct Value<TValue TClass:: * const> {
	typedef TValue const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Raw pointer to member - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.initEdgeMap:
..cat:Graph
..summary:Initializes an edge map
..signature:initEdgeMap(g, pm)
..signature:initEdgeMap(g, pm, prop)
..param.g:A Graph.
...type:Class.Graph
..param.pm:An Internal Property Map.
...type:Class.InternalRawMap
..param.prop:An optional array with properties that are to be assigned to the items in the property map.
...remarks:For every edge id there must be an entry in the array.
..returns:void
*/

///.Function.initEdgeMap.param.pm.type:Class.InternalRawMap

template <typename TSpec, typename TClass, typename TValue> 
inline void
initEdgeMap(Graph<TSpec> const& g,
			 TValue TClass:: * ptr_to_member)
{
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TClass, typename TValue> 
inline void
initEdgeMap(Graph<TSpec>& g,
			 TValue TClass:: * ptr_to_member)
{
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignProperty:
..cat:Graph
..summary:Assigns a property to an item in the property map.
..signature:assignProperty(pm, e, val)
..param.pm:An Internal Property Map.
...type:Class.InternalRawMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..param.val:The new value.
...remarks:Type of the new value must match the value type of the property map.
..returns:void
..see:Function.getProperty
..see:Function.property
*/

///.Function.assignProperty.param.pm.type:Class.InternalRawMap

template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline void 
assignProperty(TValue TClass:: * ptr_to_member, 
			TEdgeDescriptor const e, 
			TValue const val) 
{
	SEQAN_CHECKPOINT
	(cargo(e)).*ptr_to_member=val; 
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.property:
..cat:Graph
..summary:Accesses the property to an item in the property map.
..signature:property(pm, e)
..param.pm:An Internal Property Map.
...type:Class.InternalRawMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..returns:Reference to the property.
..see:Function.getProperty
..see:Function.assignProperty
*/

///.Function.property.param.pm.type:Class.InternalRawMap


template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline TValue& 
property(TValue TClass:: * const ptr_to_member, 
			TEdgeDescriptor const e) 
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*ptr_to_member; 
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getProperty:
..cat:Graph
..summary:Get method for an item's property.
..signature:getProperty(pm, e)
..param.pm:An Internal Property Map.
...type:Class.InternalRawMap
..param.e:An edge descriptor.
...remarks:Identifies the item in the property map.
..returns:The value of that item's property.
..see:Function.property
..see:Function.assignProperty
*/

///.Function.getProperty.param.pm.type:Class.InternalRawMap

template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline TValue
getProperty(TValue TClass:: * const ptr_to_member, 
			TEdgeDescriptor const e) 
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*ptr_to_member; 
} 



//////////////////////////////////////////////////////////////////////////////
// Init functions for all Maps - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPropertyMap, typename TProperties>
inline void
initVertexMap(Graph<TSpec> const& g,
			  TPropertyMap& pm,
			  TProperties const& prop)
{
	SEQAN_CHECKPOINT
	initVertexMap(g,pm);
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, VertexIterator<> >::Type TVertexIterator;
	TVertexIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(pm,*it,prop[_getId(*it)]);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPropertyMap, typename TProperties>
inline void
initEdgeMap(Graph<TSpec> const& g,
			TPropertyMap& pm,
			TProperties const& prop)
{
	SEQAN_CHECKPOINT
	initEdgeMap(g,pm);
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, EdgeIterator<> >::Type TEdgeIterator;
	TEdgeIterator it(g);
	for(;!atEnd(it);goNext(it)) {
		assignProperty(pm,*it,prop[_getId(*it)]);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
