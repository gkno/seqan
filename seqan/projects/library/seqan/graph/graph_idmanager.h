#ifndef SEQAN_HEADER_GRAPH_IDMANAGER_H
#define SEQAN_HEADER_GRAPH_IDMANAGER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Default IdManager
//////////////////////////////////////////////////////////////////////////////
template<typename TIdType, typename TSpec>
class IdManager 
{
	public:
		String<TIdType> data_freeIds;  
		String<bool> data_in_use;   //1 = in use, 0 = not in use

//____________________________________________________________________________	
	public:
		IdManager()
		{
			SEQAN_CHECKPOINT
		}

		~IdManager() 
		{
			SEQAN_CHECKPOINT
		}

		IdManager(IdManager const & _other)
		{
			SEQAN_CHECKPOINT
			data_freeIds = _other.data_freeIds;
			data_in_use = _other.data_in_use;
		}

		IdManager const& 
		operator = (IdManager const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_freeIds = _other.data_freeIds;
			data_in_use = _other.data_in_use;
			return *this;
		}

//____________________________________________________________________________
};
	

//////////////////////////////////////////////////////////////////////////////
// IdManager Specific Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TIdType, typename TSpec> 
struct Value<IdManager<TIdType, TSpec> > 
{
	typedef TIdType Type;
};

template<typename TIdType, typename TSpec> 
struct Value<IdManager<TIdType, TSpec> const> 
{
	typedef TIdType Type;
};


//////////////////////////////////////////////////////////////////////////////
//	IdManager INTERFACE
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TIdType, typename TSpec>
inline typename Value<IdManager<TIdType, TSpec> >::Type 
obtainId(IdManager<TIdType, TSpec>& idm) 
{
	SEQAN_CHECKPOINT

	TIdType id;
	if (!empty(idm.data_freeIds)) {
		id = getValue(idm.data_freeIds, length(idm.data_freeIds) - 1);
		_setLength(idm.data_freeIds, length(idm.data_freeIds) - 1);
		assignValue(idm.data_in_use, id, true);
	} else {
		id = length(idm.data_in_use);
		resize(idm.data_in_use, id + 1, Generous());
		assignValue(idm.data_in_use, id, true);
	}
	return id;
}

template<typename TIdType, typename TSpec, typename TId>
inline void 
releaseId(IdManager<TIdType, TSpec>& idm, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(idInUse(idm,id) == true)
	if (id == (TId) length(idm.data_in_use) - 1) {
		_setLength(idm.data_in_use, length(idm.data_in_use) - 1);
	} else {
		assignValue(idm.data_in_use, id, false);
		appendValue(idm.data_freeIds, id);
	}
}

template<typename TIdType, typename TSpec>
inline void 
releaseAll(IdManager<TIdType, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	clear(idm.data_freeIds);
	clear(idm.data_in_use);
}

template<typename TIdType, typename TSpec>
inline typename Size<IdManager<TIdType, TSpec> >::Type 
getIdUpperBound(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return length(idm.data_in_use);
}

template<typename TIdType, typename TSpec>
inline typename Size<IdManager<TIdType, TSpec> >::Type 
getIdLowerBound(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	typedef typename Size<IdManager<TIdType, TSpec> >::Type TSize;
	for(TSize it = 0; it < length(idm.data_in_use); ++it) {
		if (getValue(idm.data_in_use, it)) return it;
	}
	return 0;
}

template<typename TIdType, typename TSpec>
inline typename Size<IdManager<TIdType, TSpec> >::Type 
idCount(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return (length(idm.data_in_use) - length(idm.data_freeIds));
}

template<typename TIdType, typename TSpec, typename TId>
inline bool 
idInUse(IdManager<TIdType, TSpec> const& idm, 
		TId const id)
{
	SEQAN_CHECKPOINT
	SEQAN_ASSERT(id < (TId) length(idm.data_in_use))
	return getValue(idm.data_in_use, id);
}


//////////////////////////////////////////////////////////////////////////////
// Dummy IdManager
//////////////////////////////////////////////////////////////////////////////
template<typename TSpec>
class IdManager<void, TSpec> 
{
	public:
		typedef typename Id<IdManager>::Type TIdType;
		TIdType data_idCount;

//____________________________________________________________________________	
	public:
		IdManager() : data_idCount(0) 
		{
			SEQAN_CHECKPOINT
		}

		~IdManager() 
		{
			SEQAN_CHECKPOINT
		}

		IdManager(IdManager const & _other) : data_idCount(_other.data_idCount) 
		{
			SEQAN_CHECKPOINT
		}

		IdManager const& 
		operator = (IdManager const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_idCount = _other.data_idCount;
			return *this;
		}

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Dummy IdManager Specific Metafunctions
//////////////////////////////////////////////////////////////////////////////
template<typename TSpec> 
struct Value<IdManager<void, TSpec> > {
	typedef unsigned int Type;
};

template<typename TSpec> 
struct Value<IdManager<void, TSpec> const> {
	typedef unsigned int Type;
};


//////////////////////////////////////////////////////////////////////////////
//	Dummy IdManager INTERFACE
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////
template<typename TSpec>
inline typename Value<IdManager<void, TSpec> >::Type 
obtainId(IdManager<void, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	++idm.data_idCount;
	return 0;
}

template <typename TSpec, typename TId>
inline void 
releaseId(IdManager<void, TSpec>& idm, 
		  TId const id) 
{
	SEQAN_CHECKPOINT
	if (idm.data_idCount > 0) --idm.data_idCount;
}

template<typename TSpec>
inline void 
releaseAll(IdManager<void, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	idm.data_idCount = 0;
}

template<typename TSpec>
inline typename Size<IdManager<void, TSpec> >::Type 
getIdUpperBound(IdManager<void, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	// Must be data_idCount in order to resize property maps!!!
	// Don't change to 0
	return idm.data_idCount;
}

template<typename TSpec>
inline typename Size<IdManager<void, TSpec> >::Type 
getIdLowerBound(IdManager<void, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return 0;
}

template <typename TSpec>
inline typename Size<IdManager<void, TSpec> >::Type 
idCount(IdManager<void, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return idm.data_idCount;
}

template <typename TSpec, typename TId>
inline bool 
idInUse(IdManager<void, TSpec> const& idm, 
		TId const id) 
{
	SEQAN_CHECKPOINT
	return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
