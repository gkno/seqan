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

	SegmentInfo() :
		data_seq_id(0),
		data_begin(0),
		data_length(0)
	{
	}

	SegmentInfo(TId id, TSize beg, TSize len) :
		data_seq_id(id),
		data_begin(beg),
		data_length(len)
	{
	}

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
		Graph<Undirected<TCargo, TSpec> > data_align;

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

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec> 
inline unsigned int
addSequence(Graph<Alignment<TAlphabet, TCargo, TSpec> >& g,
			String<TAlphabet>& str) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TAlphabet, TCargo, TSpec> > TGraph;
	typedef String<TAlphabet> TSequence;
	typedef typename Id<TGraph>::Type TIdType;
	TIdType sId = obtainId(g.data_id_managerS);
	if (sId == length(g.data_seq)) {
		appendValue(g.data_seq, (TSequence*) &str); 
	} else {
		value(g.data_seq, sId) = (TSequence*) &str;
	}
	return sId;
}

template<typename TAlphabet, typename TCargo, typename TSpec, typename TSegment> 
inline typename VertexDescriptor<Graph<Alignment<TAlphabet, TCargo, TSpec> > >::Type 
addVertex(Graph<Alignment<TAlphabet, TCargo, TSpec> >& g,
		  TSegment& seg) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor vd = addVertex(g.data_align);
	if (length(g.data_nodeMap) <= vd) resize(g.data_nodeMap, vd + 1, Generous());
	assignProperty(g.data_nodeMap, vd, seg);
	return vd;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
