#ifndef SEQAN_HEADER_GRAPH_ALIGN_CHAINING_H
#define SEQAN_HEADER_GRAPH_ALIGN_CHAINING_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Progressive Chaining
//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Weight;


//////////////////////////////////////////////////////////////////////////////

struct _TagChaining;
typedef Tag<_TagChaining> const Chaining;


//PseudoFragment datastructure for chaining
template<typename TGraph, typename TSpec = Default>
class PseudoFragment;

template<typename TGraph>
class PseudoFragment<TGraph,Default>
{

public:
  
  TGraph const * data_g;
  String<int> data_vertices; //mischung aus vertexdescriptoren und negativen positionen
  double data_weight;
  unsigned int data_pos;
  
unsigned int _dim;
  
  PseudoFragment(){}

  PseudoFragment(unsigned int dim)
  {
	  resize(data_vertices,dim);
	  for(unsigned int i = 0; i < dim; ++i)
		data_vertices[i] = -1;		
  }

  PseudoFragment(TGraph & ali_graph, String<int> vertices,
  	   unsigned int dim, double weight, unsigned int pos)
  {
    data_g = &ali_graph;
    resize(data_vertices,dim);
    for(unsigned int i = 0; i < dim; ++i)
      data_vertices[i] = vertices[i];		
    data_weight = weight;	
	data_pos = pos;
  }

  PseudoFragment(TGraph const & ali_graph, String<int> vertices,
  	   unsigned int dim, double weight, unsigned int pos)
  {
    data_g = &ali_graph;
    resize(data_vertices,dim);
    for(unsigned int i = 0; i < dim; ++i)
      data_vertices[i] = vertices[i];		
    data_weight = weight;	
	data_pos = pos;
  }


	PseudoFragment & operator=( const PseudoFragment & old )
	{
		if( this == &old) 
			return *this;
		data_pos = old.data_pos;
		data_g = old.data_g;
		data_weight =  old.data_weight;
		resize(data_vertices,dimension(old));
		data_vertices = old.data_vertices;	
		return *this;
	}

	PseudoFragment( const PseudoFragment & old )
	{
		data_pos = old.data_pos;
		data_g = old.data_g;
		data_weight =  old.data_weight;
		resize(data_vertices,dimension(old));
		data_vertices = old.data_vertices;	
	}


	~PseudoFragment()
	{
		clear(data_vertices);
		data_pos = 0;
		data_weight = 0.0;
	}



};




template< typename TGraph > 
inline size_t
dimension( PseudoFragment<TGraph,Default> & me )
{
	return length(me.data_vertices);
}

template< typename TGraph > 
inline size_t
dimension( const PseudoFragment<TGraph,Default> & me )
{
	return length(me.data_vertices);
}

template< typename TGraph, typename TSize > 
inline int 
leftPosition( PseudoFragment<TGraph,Default>  & me, 
				TSize dim )
{
	if(me.data_vertices[dim] >= 0)
        return fragmentBegin(*me.data_g,me.data_vertices[dim])+1;
	else
		return -(me.data_vertices[dim]);

}

template< typename TGraph, typename TSize > 
inline int 
leftPosition( const PseudoFragment<TGraph,Default> & me, 
				TSize dim )
{
	if(me.data_vertices[dim] >= 0)
        return fragmentBegin(*me.data_g,me.data_vertices[dim])+1;
	else
		return -(me.data_vertices[dim]);
}

template< typename TGraph, typename TSize >
inline int 
rightPosition( PseudoFragment<TGraph,Default>  & me, 
				TSize dim )
{
	if(me.data_vertices[dim] >= 0)
        return fragmentBegin(*me.data_g,me.data_vertices[dim]) + fragmentLength(*me.data_g,me.data_vertices[dim])/*-1 ??*/;
	else
		return -(me.data_vertices[dim])-1;
}

template< typename TGraph, typename TSize > 
inline int 
rightPosition( const PseudoFragment<TGraph,Default>  & me,
				TSize dim )
{
	if(me.data_vertices[dim] >= 0)
        return fragmentBegin(*me.data_g,me.data_vertices[dim]) + fragmentLength(*me.data_g,me.data_vertices[dim])/*-1 ??*/;
	else
		return -(me.data_vertices[dim])-1;
}

template< typename TGraph > 
inline double
weight( PseudoFragment<TGraph,Default> & me )
{
	return me.data_weight;
}

template< typename TGraph > 
inline double
weight( const PseudoFragment<TGraph,Default> & me )
{
	return me.data_weight;
}


template< typename TGraph, typename TWeight >
inline typename Weight< PseudoFragment< TGraph,Default > >::Type 
setWeight(PseudoFragment<TGraph,Default> & me,
			TWeight weight )
{
	return me.data_weight = weight;
}

template< typename TGraph, typename TSize, typename TPosition > 
inline void 
_setLeftPosition( PseudoFragment<TGraph,Default> & me,
					TSize dim,
					TPosition val )
{
	SEQAN_TASSERT( dim >= 0 && dim < dimension(me));
	if(me.data_vertices[dim]<0)
		me.data_vertices[dim] = -(val+1);
}

template< typename TGraph, typename TSize, typename TPosition > 
inline void 
_setRightPosition( PseudoFragment<TGraph,Default> & me,
					TSize dim, 
					TPosition val )
{
	SEQAN_TASSERT( dim >= 0 && dim < dimension(me));
	if(me.data_vertices[dim]<0)
		me.data_vertices[dim] = -(val+1);
}




template< typename TGraph >
struct Size< PseudoFragment<TGraph,Default> >
{
	typedef size_t Type;
};


template< typename TGraph >
struct Weight< PseudoFragment<TGraph,Default> >
{
	typedef double Type;
};



template< typename TGraph >
struct Key< PseudoFragment<TGraph,Default> >
{
	typedef int Type;
};





template<typename TStringSet, typename TCargo, typename TSpec, typename TSequenceIn, typename TSequence>
inline TCargo 
_alignStringsAccordingToGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							  TSequenceIn const& seq1,
							  TSequenceIn const& seq2,
							  TSequence& alignSeq,
							  Chaining)
{
	SEQAN_CHECKPOINT
	// Clear output parameter
	clear(alignSeq);

	return chainSubsequence(g, seq1, seq2, alignSeq);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence>
inline TCargo 
_recursiveProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							   TGuideTree& tree,
							   TVertexDescriptor const root,
							   TSequence& alignSeq,
							   Chaining tag)
{
	SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;
	typedef int TPseudoVertex;

	TStringSet& str = stringSet(g);
	TCargo score = 0;

	if(isLeaf(tree, root)) {
		TId seqId = positionToId(str, root);
		TSize i = 0;
		TSize lenRoot = length(str[root]);
		while(i<lenRoot) {
			TVertexDescriptor nextVertex = findVertex(g, seqId, i);
			appendValue(alignSeq, String<TPseudoVertex>(nextVertex));
			i += fragmentLength(g, nextVertex);
		}
	} else {
		// Align the two children (Binary tree)
		typedef String<String<TPseudoVertex> > TSegmentString;
		TSegmentString seq1;
		TSegmentString seq2;
		TAdjacencyIterator adjIt(tree, root);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq1, tag);
		goNext(adjIt);
		_recursiveProgressiveAlignment(g,tree, *adjIt, seq2,tag);
		score = _alignStringsAccordingToGraph(g,seq1,seq2,alignSeq, tag);

	}
	return score;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TScore, typename TOutGraph>
inline typename Value<TScore>::Type
iterativeProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							  TGuideTree& tree,
							  TScore const& score_type,
							  TOutGraph& gOut,
							  Chaining tag)			 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	
	typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
	typedef int TPseudoVertex;

	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	// Initialization
	TSize nseq = length(stringSet(g));

	// Perform initial progressive alignment
	typedef String<TPseudoVertex> TVertexString;
	typedef String<TVertexString> TSegmentString;
	TSegmentString alignSeq;

	_recursiveProgressiveAlignment(g,tree,getRoot(tree),alignSeq,tag);
	TScoreValue maxSumScore = sumOfPairsScore(g, alignSeq, score_type,tag);


	// Build cutting order of edges
	String<TVertexDescriptor, Block<> > finalOrder;
	std::deque<TVertexDescriptor> toDoList;
	toDoList.push_back(getRoot(tree));
	while(!toDoList.empty()) {
		TVertexDescriptor v = toDoList[0];
		toDoList.pop_front();
		push_back(finalOrder, v);
		if (!isLeaf(tree, v)) {
			TAdjacencyIterator adjIt(tree, v);
			toDoList.push_back(*adjIt);
			goNext(adjIt);
			toDoList.push_back(*adjIt);
		}
	}

	// Iterative alignment of profiles
	TSize len = length(finalOrder) - 1;  // Don't touch the root again
	TSize iterationsWithoutImprovement = 0;
	for(TSize edge_count=0; ((edge_count<len) && (iterationsWithoutImprovement<5));++edge_count) {
		TVertexDescriptor subtree_root = top(finalOrder);
		pop(finalOrder);

		// Collect all vertex descriptors that belong to the subtree
		std::set<TVertexDescriptor> subtree;
		toDoList.clear();
		toDoList.push_back(subtree_root);
		TSize numSeqs2 = 0;
		while(!toDoList.empty()) {
			TVertexDescriptor v = toDoList[0];
			toDoList.pop_front();
			if (!isLeaf(tree, v)) {
				TAdjacencyIterator adjIt(tree, v);
				toDoList.push_back(*adjIt);
				goNext(adjIt);
				toDoList.push_back(*adjIt);
			} else {
				// Insert the sequence id into the subtree set
				subtree.insert(positionToId(stringSet(g), v));
				++numSeqs2;
			}
		}

		// Build position -> id map
		TSize numSeqs1 = nseq - numSeqs2;
		TSize alignSeqLen = length(alignSeq);
		std::map<unsigned int,TId> profileMap;
		typename Iterator<TSegmentString>::Type alignSeq_it = begin(alignSeq);
		typename Iterator<TSegmentString>::Type alignSeq_end = end(alignSeq);
		TSize vertexSetLen = length(value(alignSeq_it));
		for(TSize j=0; j<vertexSetLen; ++j) {
			TPseudoVertex v = getValue(value(alignSeq_it), j);
			while(alignSeq_it != alignSeq_end && v < 0)
			{
				++alignSeq_it;
				v = getValue(value(alignSeq_it), j);
			}
			if (subtree.find(sequenceId(g,v)) == subtree.end()) profileMap[j] = 0;
			else profileMap[j] = 1;
		}
	
		// Build the 2 profile strings
		TSegmentString seq1;
		TSegmentString seq2;
		for(TSize i = 0; i<alignSeqLen;++i) {
			TVertexString set1; TVertexString set2;
			TSize count1 = 0; TSize count2 = 0;
			bool gibts1 = false, gibts2 = false;
			TVertexString& alignSeq_i = alignSeq[i];
			TSize vertexSetLen = length(alignSeq_i);
			for(TSize j=0; j<vertexSetLen; ++j) {
				resize(set1, numSeqs1);
				resize(set2, numSeqs2);
				TPseudoVertex v = getValue(alignSeq_i, j);
				if (v < 0)
				{
					typename std::map<unsigned int, TId>::iterator pmIt = profileMap.find(j);
					if(pmIt->second == 0)
						set1[count1++] = v;
					else if(pmIt->second == 1)
						set2[count2++] = v;
					else std::cout << "schlimm\n";
					continue;
				}
				else if (subtree.find(sequenceId(g,v)) == subtree.end()){ gibts1 = true; set1[count1++] = v;}
				else {gibts2 = true; set2[count2++] = v;}
			}
			if (gibts1) appendValue(seq1, set1);
			if (gibts2) appendValue(seq2, set2);
		}


		// Align profile strings
		TSegmentString localAlignSeq;
		_alignStringsAccordingToGraph(g,seq1,seq2,localAlignSeq,tag);
		TScoreValue localSumScore = sumOfPairsScore(g, localAlignSeq, score_type,tag);

		// New maximum?
		if (localSumScore > maxSumScore) {
			iterationsWithoutImprovement = 0;
			maxSumScore = localSumScore;
			alignSeq = localAlignSeq;
		} else {
			++iterationsWithoutImprovement;
		}
	}

	// Create the alignment graph
	TSize alignSeqLen = length(alignSeq);
	for(TSize i = 0; i<alignSeqLen;++i) {
		TVertexString& alignSeq_i = alignSeq[i];
		TSize len_i = length(alignSeq_i);
		for(TSize j=0; j<len_i; ++j) {
			TPseudoVertex v = getValue(alignSeq_i, j);
			if (v < 0) continue;
			SEQAN_TASSERT(fragmentBegin(g,v) < length(getValueById(stringSet(g), sequenceId(g,v))))
			SEQAN_TASSERT(fragmentLength(g,v) > 0)
			SEQAN_TASSERT(fragmentBegin(g,v) + fragmentLength(g,v) <= length(getValueById(stringSet(g), sequenceId(g,v))))
			TVertexDescriptor l = addVertex(gOut, sequenceId(g, v), fragmentBegin(g,v), fragmentLength(g,v));
			//std::cout << l << label(gOut, l) << ',';
			TSize count = 1;
			for(TSize k = j; k>0; --k) {
				if (getValue(alignSeq_i, k - 1) >= 0) {
					SEQAN_TASSERT(fragmentLength(gOut,l) == fragmentLength(gOut,l - count))
					addEdge(gOut, l - count, l);
					++count;
				}
			}
		}
		//std::cout << std::endl;
	}
	return maxSumScore;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TString>
inline TCargo
chainSubsequence(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
						  TString const& str1, 
						  TString const& str2,
						  TString& align) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	//pseudo vertex descriptors
	//negative values replace nilVertices
	//in that case -value-1 indicates the position on the sequence
    typedef int TPseudoVertex;
	typedef typename Iterator<TString const>::Type TStringIter;
	typedef typename Iterator<TString>::Type TSIter;
	typedef typename Value<TString>::Type TVertexSet;
	typedef typename Iterator<TVertexSet const>::Type TVertexSetIter;
	typedef typename Iterator<TVertexSet>::Type TIter;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	
	// Size of the sequences
	// Note for profile alignments every member of the sequence is a String!!! of pseudo vertex descriptors
	TSize m = length(str1);  // How many sets of vertex descriptors in seq1
	TSize n = length(str2);  // How many sets of vertex descriptors in seq1
	TSize seqsInStr1 = length(str1[0]);	 // #Vertex descriptors per node
	TSize seqsInStr2 = length(str2[0]);	
	TSize dim = seqsInStr1 + seqsInStr2;
	
	// Fill the vertex to position map for str1
	// For each vertex descriptor remember the position in the sequence
	typedef std::map<TPseudoVertex, TSize> TVertexToPosMap;
	typedef typename TVertexToPosMap::const_iterator TVertexToPosMapIter;
	TVertexToPosMap map;
	TStringIter itStrEnd1 = end(str1);
	for(TStringIter itStr1 = begin(str1);itStr1 != itStrEnd1;++itStr1) {
		TVertexSetIter itVEnd = end(getValue(itStr1));
		for(TVertexSetIter itV = begin(getValue(itStr1));itV != itVEnd;++itV) {
			if (*itV >= 0) map.insert(std::make_pair(*itV, position(itStr1)));
		}
	}

	// all edges of the full graph are created
	// Initially every edge receives weight=0
	String<double> weights;
	// For profile alignments, take the average weight
	double divider = (double) seqsInStr1 * (double) seqsInStr2;
	fill(weights, n*m, 0);

	// Walk through str2 and fill in the weights of the actual edges
	TStringIter itStrEnd2 = end(str2);
	for(TStringIter itStr2 = begin(str2);itStr2 != itStrEnd2;++itStr2) {
		TVertexSetIter itVEnd = end(getValue(itStr2));
		for(TVertexSetIter itV = begin(getValue(itStr2));itV != itVEnd;++itV) {
			if (*itV >= 0) {
				TOutEdgeIterator itOut(g, *itV);
				for(;!atEnd(itOut); ++itOut) {
					// Target vertex must be in the map
					TVertexToPosMapIter pPos = map.find(targetVertex(itOut));
					if (pPos != map.end()) {
						// Calculate the edge index
						TSize index = pPos->second + (m * position(itStr2) );
						weights[index] += (double) cargo(*itOut) / divider;	
					}
				}
			}
		}
	}
	map.clear();

	//we only need those edges that have a positive score
	typedef PseudoFragment<TGraph> TFragment;
	typedef String<TFragment> TFragString;
	typedef typename Iterator<TFragString,Standard>::Type TFragStringIter;
	String<PseudoFragment<TGraph> > pseudo_fragments;
	for(unsigned int i = 0; i < length(weights); ++i)
	{
		if(weights[i] > 0)
		{
			TSize index1 = i%m;
			TSize index2 = i/m; //floor
			String<TPseudoVertex> vertices;
			append(vertices,str1[index1]);
			append(vertices,str2[index2]);
			appendValue(pseudo_fragments,TFragment(g,vertices,dim,weights[i],i));
		}
	}
	clear(weights);

	

	String<PseudoFragment<TGraph> > chained_fragments;
	reserve(chained_fragments,length(pseudo_fragments)+2);

	double score = compute_chain(pseudo_fragments,chained_fragments,Score<double, Zero>(),SemiDeferred());
	


	// Create the alignment sequence
	TSize numMatches = length(chained_fragments)-2;
	TSize alignLength = numMatches + (n - numMatches) + (m - numMatches);
	clear(align);
	resize(align, alignLength);
	TSIter pointerAlign = begin(align);
	TSIter pointerAlignEnd = end(align);
	TStringIter pointerStr1 = begin(str1);
	TStringIter pointerStr2 = begin(str2);
	TFragStringIter chainIt = begin(chained_fragments,Standard()) + 1;
	TFragStringIter chainEnd = end(chained_fragments,Standard()) -1;

	String<TPseudoVertex> pseudo_nils1;
	fill(pseudo_nils1, seqsInStr1, -1);
	String<TPseudoVertex> pseudo_nils2;
	fill(pseudo_nils2, seqsInStr2, -1);

	int counter = 0;
	int counterAlign = 0;

	TSize ii = 0; 
	TSize jj = 0;
	while(pointerAlign != pointerAlignEnd && chainIt != chainEnd) {
		counter++;
	    counterAlign++;
		TSize index1 = value(chainIt).data_pos%m;
		TSize index2 = value(chainIt).data_pos/m;

		while(ii != index1)
		{
			value(pointerAlign) = str1[ii];
			append(*pointerAlign,pseudo_nils2);
			++ii;++pointerAlign;
			counterAlign++;
		}
		while(jj != index2)
		{
			value(pointerAlign) = pseudo_nils1;
			append(*pointerAlign,str2[jj]);
            ++jj;++pointerAlign;
			counterAlign++;
		}
		append(*pointerAlign,value(chainIt).data_vertices);
		++ii;++jj;
		for(unsigned int dim_i = 0; dim_i < seqsInStr1; ++dim_i)
			pseudo_nils1[dim_i] = -(rightPosition(*chainIt,dim_i)+1);
		for(unsigned int dim_j = 0; dim_j < seqsInStr2; ++dim_j)
			pseudo_nils2[dim_j] = -(rightPosition(*chainIt,seqsInStr1+dim_j)+1);
		++pointerAlign;
		++chainIt;
	}
	if(pointerAlign != pointerAlignEnd)
	{
		TSize index1 = m; 
		TSize index2 = n;
		while(ii != index1)
		{
			value(pointerAlign) = str1[ii];
			append(*pointerAlign,pseudo_nils2);
			++ii;++pointerAlign;
			counterAlign++;
		}
		while(jj != index2)
		{
			value(pointerAlign) = pseudo_nils1;
			append(*pointerAlign,str2[jj]);
            ++jj;++pointerAlign;
			counterAlign++;
		}
	}
	SEQAN_TASSERT(chainIt == chainEnd)
	SEQAN_TASSERT(position(pointerAlign) == length(align))

	return (TCargo)score;
}



template<typename TStringSet, typename TCargo, typename TSpec, typename TSegmentString, typename TScore> 
inline typename Value<TScore>::Type
sumOfPairsScore(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				TSegmentString const& alignSeq,
				TScore const& score_type,
				Chaining)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef int TVertexDescriptor;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	// Initialization
	TScoreValue gap = scoreGapExtend(score_type);
	TScoreValue gapOpen = scoreGapOpen(score_type);
	TScoreValue total = 0;

	// Sum of pair scores
	TSize alignSeqLen = length(alignSeq);
	for(TSize i = 0; i<alignSeqLen;++i) {
		TSize count = 0;
		TSize vertexSetLen = length(alignSeq[i]);
		TSize fragLen = 0;
		for(TSize j=0; j<vertexSetLen;++j) {
			TVertexDescriptor v1 = getValue(alignSeq[i], j);
			if ((fragLen == 0) && (v1 >= 0)) fragLen = fragmentLength(g, v1);
			for(TSize k=j+1; k<vertexSetLen;++k) {
				TVertexDescriptor v2 = getValue(alignSeq[i], k);
				//negative values indicate positions on the sequences
				if ((v1 < 0) ||
					(v2 < 0))
				{
					// Count number of pairs where one is a vertex and the other is just a position 
					if (!((v1 < 0) &&
						(v2 < 0))) ++count;
					continue;
				}
				typedef typename Iterator<TInfix,Standard>::Type TInfixIter;
				TInfix inf1 = label(g,v1);
				TInfix inf2 = label(g,v2);
				TInfixIter sIt1 = begin(inf1,Standard());
				TInfixIter sIt2 = begin(inf2,Standard());
				TInfixIter sItEnd1 = end(inf2,Standard());
				TInfixIter sItEnd2 = end(inf2,Standard());
				while((sIt1 != sItEnd1) || (sIt2 != sItEnd2)) {
					total += score(const_cast<TScore&>(score_type), *sIt1, *sIt2);
					goNext(sIt1); goNext(sIt2);
				}
			}
		}
		SEQAN_TASSERT(fragLen > 0)
		// How many sequences are left? --> Aligned with gaps
		total += (( (TScoreValue) (count) ) * (gapOpen + ( (TScoreValue) fragLen - 1) * gap));
	}
	return total;
}



//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...




