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
		std::deque< TIdType > data_freeIds;
		TIdType data_endId;

//____________________________________________________________________________	
	public:
		IdManager() : data_endId(0) 
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
			data_endId = _other.data_endId;
		}

		IdManager const& 
		operator = (IdManager const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_freeIds = _other.data_freeIds;
			data_endId = _other.data_endId;
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
	if (!idm.data_freeIds.empty()) {
		id = idm.data_freeIds.front();
		idm.data_freeIds.pop_front();
	} else {
		id = idm.data_endId;
		++idm.data_endId;
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

	if (id == idm.data_endId - 1) {
		--idm.data_endId;
		// Trim the data_freeIds if possible
		while ( (idm.data_endId != 0) && 
				(!idm.data_freeIds.empty()) &&
				(idm.data_freeIds.back() == idm.data_endId - 1)) 
		{
				idm.data_freeIds.pop_back();
				--idm.data_endId;
		}
	} else {
		if (idm.data_freeIds.empty()) idm.data_freeIds.push_back(id);
		else if (id < idm.data_freeIds.front()) idm.data_freeIds.push_front(id);
		else if (id > idm.data_freeIds.back()) idm.data_freeIds.push_back(id);
		else {
			typedef typename std::deque<TIdType>::iterator TDequeIterator;
			for(TDequeIterator pos = idm.data_freeIds.begin(); pos != idm.data_freeIds.end(); ++pos) {
				if (id < *pos) {
					idm.data_freeIds.insert(pos, id);
					break;
				}
			}
		}
	}
}

template<typename TIdType, typename TSpec>
inline void 
releaseAll(IdManager<TIdType, TSpec>& idm) 
{
	SEQAN_CHECKPOINT
	idm.data_freeIds.clear();
	idm.data_endId = 0;
}

template<typename TIdType, typename TSpec>
inline typename Size<IdManager<TIdType, TSpec> >::Type 
getIdUpperBound(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return idm.data_endId;
}

template<typename TIdType, typename TSpec>
inline typename Size<IdManager<TIdType, TSpec> >::Type 
getIdLowerBound(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	typedef typename Size<IdManager<TIdType, TSpec> >::Type TSize;
	TSize count = 0;
	for(TSize it = 0; it < idm.data_freeIds.size(); ++it) {
		if (idm.data_freeIds[it] != count) return count;
		++count;
	}
	return count;
}

template<typename TIdType, typename TSpec>
inline typename Size<IdManager<TIdType, TSpec> >::Type 
idCount(IdManager<TIdType, TSpec> const& idm)
{
	SEQAN_CHECKPOINT
	return (idm.data_endId - idm.data_freeIds.size());
}

template<typename TIdType, typename TSpec, typename TId>
inline bool 
idInUse(IdManager<TIdType, TSpec> const& idm, 
		TId const id)
{
	SEQAN_CHECKPOINT
	if (id >= idm.data_endId) return false;
	for(unsigned int it = 0; ((it < idm.data_freeIds.size()) && (idm.data_freeIds[it]<=id)); ++it) {
		if (idm.data_freeIds[it] == id) return false;
	}
	return true;
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
