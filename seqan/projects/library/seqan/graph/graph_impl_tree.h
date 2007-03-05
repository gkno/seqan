#ifndef SEQAN_HEADER_GRAPH_IMPL_EDGELISTU_H
#define SEQAN_HEADER_GRAPH_IMPL_EDGELISTU_H

namespace SEQAN_NAMESPACE_MAIN
{
////////////////
// Tree that stores the tree edges in a list
////////////////
template<typename TCargo, typename TEdgeSpec, typename TSpec>
class Graph<EdgeListT<TCargo, TEdgeSpec>, TSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdgeStumpT;	
		typedef Allocator<SinglePool<sizeof(TEdgeStumpT)> > TAllocator;
		
		String<TEdgeStumpT*> data_vertex;			// Pointers to EdgeStumpT lists
		IdManager<TIdType> data_id_managerV;
		TAllocator data_allocator;

//____________________________________________________________________________


		Graph() {
			SEQAN_CHECKPOINT
		}


		~Graph() {
			SEQAN_CHECKPOINT
			//clear(*this);
		}

		Graph(Graph const & _other) :
			data_allocator(_other.data_allocator)
		{
			SEQAN_CHECKPOINT
			//_copyGraph(_other, *this);		
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			/*
			clear(*this);
			data_allocator = _other.data_allocator;
			_copyGraph(_other, *this);
			*/
			return *this;
		}
};


//////////////////////////////////////////////////////////////////////////////
// EdgeListT specific graph functions
//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
