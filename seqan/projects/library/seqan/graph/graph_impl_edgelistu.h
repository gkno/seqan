#ifndef SEQAN_HEADER_GRAPH_IMPL_EDGELISTU_H
#define SEQAN_HEADER_GRAPH_IMPL_EDGELISTU_H

namespace SEQAN_NAMESPACE_MAIN
{
////////////////
// Undirected graph which stores the edges in a list
////////////////
template<typename TCargo, typename TEdgeSpec, typename TSpec>
class Graph<EdgeListU<TCargo, TEdgeSpec>, TSpec> 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename EdgeType<Graph>::Type TEdgeStump;	
		typedef typename IdHandler<TEdgeStump, TIdType>::Type TEdgeIdManager;
		typedef Allocator<SinglePool<sizeof(TEdgeStump)> > TAllocator;
		
		String<TEdgeStump*> data_vertex;			// Pointers to EdgeStump lists
		IdManager<TIdType> data_id_managerV;
		TEdgeIdManager data_id_managerE;		
		TAllocator data_allocator;

//____________________________________________________________________________


		Graph() {
			SEQAN_CHECKPOINT
		}


		template<typename TEdgeArray, typename TSize>
		Graph(TEdgeArray edges, TSize size) {
			SEQAN_CHECKPOINT
			//_copyGraph(edges,size,*this);
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
// EdgeListU specific graph functions
//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
