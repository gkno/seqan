#ifndef SEQAN_HEADER_GRAPH_IMPL_ALIGN_H
#define SEQAN_HEADER_GRAPH_IMPL_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Alignment
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TSize>
class SegmentInfo {
public:
	TId data_seq_id;
	TSize data_begin;
	TSize data_length;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec>
class Graph<Alignment<TAlphabet, TCargo, TSpec> > 
{
	public:
		typedef typename Id<Graph>::Type TIdType;
		typedef typename Size<Graph>::Type TSize;
		typedef String<TAlphabet> TSequence;

		// Alignment graph
		Graph<Undirected<TCargo, TSpec> >* align;

		// Sequences
		IdManager<TIdType> data_id_managerS;
		String<TSequence*> data_seq;

		// Alignment specific members
		String<SegmentInfo<TIdType, TSize> > data_nodeMap;

//____________________________________________________________________________


		Graph() {
			SEQAN_CHECKPOINT
		}


		~Graph() {
			SEQAN_CHECKPOINT
		}

		Graph(Graph const & _other) 
		{
			SEQAN_CHECKPOINT
		}
	
		Graph const& operator = (Graph const & _other) {
			SEQAN_CHECKPOINT
			return *this;
		}
};

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
