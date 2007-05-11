#ifndef SEQAN_HEADER_GRAPH_TEMP_REFINENEW_H
#define SEQAN_HEADER_GRAPH_TEMP_REFINENEW_H

namespace SEQAN_NAMESPACE_MAIN
{


 
///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functions for Align<TSource,TSpec>
//project onto other sequence 
template<typename TSource,typename TSpec,typename TValue,typename TMap>
bool
_getOtherSequenceAndProject(Align<TSource,TSpec> & segment, 
						   TMap & seq_map, TValue seq_i, TValue node_i, TValue & seq_j, TValue & node_j)
{
SEQAN_CHECKPOINT
	TValue ali_seq_0 = seq_map[id(source(row(segment,0)))];
	if(seq_i == ali_seq_0)
	{
		seq_j = seq_map[id(source(row(segment,1)))];
		node_j = toSourcePosition(row(segment,1),toViewPosition(row(segment,0),node_i));
		return true;
	}
	else
	{
		seq_j = ali_seq_0;
		node_j = toSourcePosition(row(segment,0),toViewPosition(row(segment,1),node_i));
		return false;
	}

}


//unspektakuläre funktion, die die int ID zurückgibt (braucht man damit es für alle alignment typen geht)
template<typename TSource,typename TSpec, typename TValue, typename TSeqMap>					
int 
_getSeqMapId(TSeqMap & seq_map,
			Align<TSource,TSpec> & segment,
			TValue seq_i)
{
SEQAN_CHECKPOINT
	return seq_map[id(source(row(segment,seq_i)))];
}



//given seq and segment, get the sequenceId (seq_i) and its begin and end
//if seq = 0 get first sequence (that takes part in the segment match)
//if seq = 1 get second sequence
template<typename TAliSource,typename TAliSpec, typename TValue>
void
getSeqBeginAndEnd(Align<TAliSource,TAliSpec> & segment,
				  std::map<const void * ,int> & seq_map, 
				  TValue & seq_i, 
				  TValue & begin_i, 
				  TValue & end_i,
				  TValue seq)
{
	seq_i = seq_map[id(source(row(segment,seq)))];
	begin_i = sourceBeginPosition(row(segment,seq));
	end_i = sourceEndPosition(row(segment,seq));
}




///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functios for Align Graphs
//project onto other sequence for Graph<Alignment>
template<typename TAlignment,typename TValue, typename TMap>
bool
_getOtherSequenceAndProject(Graph<TAlignment> & segment, TMap & seq_map, TValue seq_i, TValue pos_i, TValue & seq_j, TValue & pos_j)
{
SEQAN_CHECKPOINT

	typedef Graph<TAlignment> TGraph;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertex;
	
	TVertex vd1 = findVertex(segment,seq_i,pos_i);
	TEdgeIterator it(segment,vd1);
	TVertex vd2 = targetVertex(it);
	seq_j = sequenceId(segment,vd2);
	pos_j = fragmentBegin(segment,vd2) + (pos_i - fragmentBegin(segment,vd1));

	if(vd1 == 0)
		return true;
	else
		return false;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functions for Fragments
//project onto other sequence for Graph<Alignment>
template<typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec,typename TValue, typename TMap>
bool
_getOtherSequenceAndProject(Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment,
						   TMap & seq_map,
						   TValue seq_i,
						   TValue pos_i,
						   TValue & seq_j,
						   TValue & pos_j)
{
SEQAN_CHECKPOINT

	pos_j = getProjectedPosition(segment,seq_i, pos_i);

	if(seq_i == sequenceId(segment,0))
	{
		seq_j = sequenceId(segment,1);
		return true;
	}
	else
	{
		seq_j = sequenceId(segment,0);
		return false;
	}

}

//////////////////////////////////////////////////////////////////////////////////////////
//Functions for all Alignment types except Align<TSource,TSpec>
//unspektakuläre funktion, die die int ID zurückgibt (braucht man damit es für alle alignment typen geht)
template<typename TAlignment, typename TValue,typename TVertexDescriptor, typename TSeqMap>					
int
_getSeqMapId(TSeqMap & seq_map,
			TAlignment & segment,
			TVertexDescriptor seq_i)
{
SEQAN_CHECKPOINT
	return sequenceId(segment,seq_i);
}


//given seq and segment, get the sequenceId (seq_i) and its begin and end
//if seq = 0 get first sequence (that takes part in the segment match)
//if seq = 1 get second sequence
template<typename TValue,typename TAlign>
void
getSeqBeginAndEnd(TAlign & segment,
				  std::map<const void * ,int> & seq_map, 
				  TValue & seq_i, 
				  TValue & begin_i, 
				  TValue & end_i,
				  TValue seq)
{
SEQAN_CHECKPOINT
	seq_i = sequenceId(segment,seq);
	begin_i = fragmentBegin(segment,seq_i);
	end_i = begin_i + fragmentLength(segment,seq_i);
}








///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Recursive Refinement
//refine position node_i on sequence seq_i
template<typename TValue, typename TAlignmentString, typename TGraph, typename TPropertyMap,typename TSeqMap>
void
_refine(TValue node_i, 
	 TValue seq_i, 
	 TSeqMap & seq_map,
	 TAlignmentString & alis, 
	 String<TGraph> & gs, 
	 String<TPropertyMap> & pms, 
     String<std::set<TValue> > & all_nodes)
{
SEQAN_CHECKPOINT
	typedef typename Cargo<typename Value<TPropertyMap>::Type>::Type TAlignmentPointer;
	typedef typename Iterator<String<TAlignmentPointer> >::Type TSegmentIterator;

	//find all segment matches that contain the current position (node_i)
	String<TAlignmentPointer> relevant_segments;
	findIntervalsExcludeTouching(gs[seq_i],pms[seq_i],node_i,relevant_segments);
	
	TSegmentIterator segment_it = begin(relevant_segments);
	TSegmentIterator segment_end = end(relevant_segments);

	//foreach of those segments
	while(segment_it != segment_end)
	{

		//get the sequence that node_i needs to be projected onto (seq_j)
		//and get the projected position (pos_j)
		TValue seq_j, node_j;
		_getOtherSequenceAndProject(alis[*segment_it],seq_map,seq_i,node_i,seq_j,node_j);

		typename std::set<TValue>::iterator iter;
		iter = all_nodes[seq_j].find(node_j);
		
		//if node does not exist yet ---> insert and continue cutting
		if(iter == all_nodes[seq_j].end())
		{
			//TODO Abbruch: if(!stop(all_nodes,seq_check,node_j,seq_j,TStop()))
			all_nodes[seq_j].insert(node_j);
			_refine(node_j,seq_j,seq_map,alis,gs,pms,all_nodes);
			//TODO: else //verschmelzen, abschneiden und übergehen, erst später... 	
			//do nothing or resolve problems  
		}
	
		++segment_it;
	}




}

template<typename TValue>
inline bool
cutIsOk(String<std::set<TValue> > & all_nodes,
		TValue seq_i,
		TValue pos_i,
		typename std::set<TValue>::iterator iter,
		TValue min_len)
{
SEQAN_CHECKPOINT
	
	if(iter != all_nodes[seq_i].end())
		return false;
//	if(min_len == 0)
//		return true;

	typename std::set<TValue>::iterator tmp_iter = all_nodes[seq_i].upperBound(pos_i);
	if(tmp_iter != all_nodes[seq_i].end())
		if((*tmp_iter - pos_i) < min_len)
			return false;
	if(tmp_iter != all_nodes[seq_i].begin())
	{
		--tmp_iter;
		if((pos_i - *tmp_iter) < min_len)
			return false;
	}

	return true;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Recursive Refinement
//refine position node_i on sequence seq_i
template<typename TValue, typename TAlignmentString, typename TGraph, typename TPropertyMap,typename TSeqMap>
void
_refine(TValue node_i, 
	 TValue seq_i, 
	 TSeqMap & seq_map,
	 TAlignmentString & alis, 
	 String<TGraph> & gs, 
	 String<TPropertyMap> & pms, 
     String<std::set<TValue> > & all_nodes, 
	 TValue min_len)
{
SEQAN_CHECKPOINT
	typedef typename Cargo<typename Value<TPropertyMap>::Type>::Type TAlignmentPointer;
	typedef typename Iterator<String<TAlignmentPointer> >::Type TSegmentIterator;

	//find all segment matches that contain the current position (node_i)
	String<TAlignmentPointer> relevant_segments;
	findIntervalsExcludeTouching(gs[seq_i],pms[seq_i],node_i,relevant_segments);
	
	TSegmentIterator segment_it = begin(relevant_segments);
	TSegmentIterator segment_end = end(relevant_segments);

	//foreach of those segments
	while(segment_it != segment_end)
	{

		//get the sequence that node_i needs to be projected onto (seq_j)
		//and get the projected position (pos_j)
		TValue seq_j, node_j;
		_getOtherSequenceAndProject(alis[*segment_it],seq_map,seq_i,node_i,seq_j,node_j);

		typename std::set<TValue>::iterator iter;
		iter = all_nodes[seq_j].find(node_j);
		
		//if node does not exist yet ---> insert and continue cutting
		if(cutIsOk(all_nodes,seq_j,pos_i,iter,min_len))
		{
			all_nodes[seq_j].insert(node_j);
			_refine(node_j,seq_j,seq_map,alis,gs,pms,all_nodes,min_len);
			//TODO: else //verschmelzen, abschneiden und übergehen, erst später... 	
			//do nothing or resolve problems  
		}
	
		++segment_it;
	}




}
	

///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Construct interval trees 
////////////////////////////////////////////////////////////////////////////////////////////////////


//construct intervals from allignments for each sequence (other Alignment types)
template<typename TInterval, typename TAlignmentString, typename TSeqMap>
void
_buildIntervalsForAllSequences(TAlignmentString & alis, 
							   String<String<TInterval> > & intervals, 
							   TSeqMap & seq_map)

{
SEQAN_CHECKPOINT
	
	typedef typename Value<TInterval>::Type TValue;
	typedef typename Cargo<TInterval>::Type TCargo;

	typedef typename Iterator<TAlignmentString,Standard>::Type TAliIterator;
	TAliIterator ali_it = begin(alis,Standard());
	TAliIterator ali_end = end(alis,Standard());

	TValue ali_counter = 0;
	//foreach alignment
	while(ali_it != ali_end)
	{
		TValue seq_i,begin_,end_;
	
		//get the first sequence (and its begin and end) that takes part in the alignment (seq_i)
		getSeqBeginAndEnd(*ali_it,seq_map,seq_i,begin_,end_,0);
		//and append the interval (ali_begin, ali_end) with cargo ali* to the list of intervals of seq_i
		appendValue(intervals[seq_i],IntervalAndCargo<TValue,TCargo>(begin_,end_,ali_counter)); 
	
		//get the second sequence (and its begin and end) that takes part in the alignment (seq_i)
		getSeqBeginAndEnd(*ali_it,seq_map,seq_i,begin_,end_,1);
		//and again append the interval (ali_begin, ali_end) with cargo ali* to the list of intervals of seq_i
		appendValue(intervals[seq_i],IntervalAndCargo<TValue,TCargo>(begin_,end_,ali_counter)); 
	
		++ali_counter;
		++ali_it;
	}


}

//get all intervals from the alignments and construct an interval tree for each sequence
template<typename TGraph, typename TPropertyMap, typename TAlignmentString, typename TSpec, typename TSequence, typename TSetSpec, typename TValue, typename TSeqMap>
void
createTreesForAllSequences(String<TGraph> & gs, 
						   String<TPropertyMap> & pms, 
						   TAlignmentString & alis, 
						   StringSet<TSequence,TSetSpec> & seqs,
                           TSeqMap & seq_map,
						   TValue numSequences,
						   Tag<TSpec> const tag)
{
SEQAN_CHECKPOINT

	typedef typename Value<TAlignmentString>::Type TAlignment;
//	typedef TAlignment* TCargo;
	typedef TValue TCargo;
	typedef IntervalAndCargo<int,TCargo> TInterval;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	
	std::cout <<"create interval trees...";
	clock_t start, finish1;
	double duration;
	start = clock();

	//one tree for each sequence
	resize(gs,numSequences);
	resize(pms,numSequences);
	
	//and one string of intervals for each sequence
	String<String<TInterval> > intervals;
	resize(intervals,numSequences);

	//fill intervals
	_buildIntervalsForAllSequences(alis,intervals,seq_map);
	
	TValue i = 0;
	
	while(i < numSequences)
	{
		std::cout << (numSequences-i) <<" more ("<<length(intervals[i])<<" intervals)... ";
		//vllt zum speicher sparen: numSequences mal alle alis durchgehen
		//und jedes mal nur buildIntervalsForJustOneSequence(); 
		TValue center = length(seqs[i])/2; // center raus, hat hier nix zu suchen

		//create interval tree!
		createIntervalTree(gs[i],pms[i],intervals[i],center,tag);
		
		//intervals for sequence i are not needed anymore
		clear(intervals[i]);
		++i;
	}

	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	std::cout << "\ntook " << duration << " seconds.\n";



}



////////////////////////////////////////////////////////////////////////////////////////
// 50000 getScore Functions
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////
//for Align<TAliSource,TAliSpec>

//get score for alignment of length len starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true)
template<typename TScore,typename TStringSet,typename TAliSource,typename TAliSpec,typename TValue>
typename Value<TScore>::Type
getScore(TScore & score_type,
		 TStringSet & seqs,
		 Align<TAliSource,TAliSpec> & segment,
		 bool i_am_first, 
		 TValue pos_i,
		 TValue pos_j,
		 TValue len)
{
SEQAN_CHECKPOINT
	typedef Align<TAliSource,TAliSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TIterator;	

	TIterator row0_it, row1_it;
	TValue len1;
	if(i_am_first)
	{
		row0_it = iter(row(segment,0),toViewPosition(row(segment,0),pos_i));
		row1_it = iter(row(segment,1),toViewPosition(row(segment,1),pos_j));
		len = toViewPosition(row(segment,0),pos_i + len) - toViewPosition(row(segment,0),pos_i);
	}
	else
	{
		row0_it = iter(row(segment,0),toViewPosition(row(segment,0),pos_j));
		row1_it = iter(row(segment,1),toViewPosition(row(segment,1),pos_i));
		len = toViewPosition(row(segment,1),pos_i + len) - toViewPosition(row(segment,1),pos_i);
	}

	int i = 0;
	typename Value<TScore>::Type ret_score = 0;
	

	while(i < len)
	{
		ret_score += score(score_type,getValue(row0_it),getValue(row1_it));
		++i;
		++row0_it;
		++row1_it; 
	}

	return ret_score;
}				
					
//get score for alignment starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true), if len1!=len2 then the refinement
//process was stopped (the cut is not exact)
template<typename TScore,typename TStringSet, typename TAliSource,typename TAliSpec,typename TValue>
typename Value<TScore>::Type
getScore(TScore & score_type,
		 TStringSet & seqs, 
		 Align<TAliSource,TAliSpec> & segment,
		 bool i_am_first, 
		 TValue pos_i,
		 TValue pos_j,
		 TValue len1,
		 TValue len2)
{
SEQAN_CHECKPOINT
	typedef Align<TAliSource,TAliSpec> TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TIterator;	

	TIterator row0_it, row1_it;
	TValue len;
	if(i_am_first)
	{
		row0_it = iter(row(segment,0),toViewPosition(row(segment,0),pos_i));
		row1_it = iter(row(segment,1),toViewPosition(row(segment,1),pos_j));
		len1 = toViewPosition(row(segment,0),pos_i + len1) - toViewPosition(row(segment,0),pos_i);
		len2 = toViewPosition(row(segment,1),pos_j + len2) - toViewPosition(row(segment,1),pos_j);
		len = (len1 < len2) ? len1 : len2;
	}
	else
	{
		row0_it = iter(row(segment,0),toViewPosition(row(segment,0),pos_j));
		row1_it = iter(row(segment,1),toViewPosition(row(segment,1),pos_i));
		len1 = toViewPosition(row(segment,1),pos_i + len1) - toViewPosition(row(segment,1),pos_i);
		len2 = toViewPosition(row(segment,0),pos_j + len2) - toViewPosition(row(segment,0),pos_j);
		len = (len1 < len2) ? len1 : len2;
	}

	int i = 0;
	typename Value<TScore>::Type ret_score = 0;
	
	//calculate score for aligned region
	while(i < len)
	{
		ret_score += score(score_type,getValue(row0_it),getValue(row1_it));
		++i;
		++row0_it;
		++row1_it;
	}
	//fill up with gaps if one sequence is longer than the other
	len = (len1 > len2) ? len1 : len2;
	ret_score += (len - i) * scoreGapExtend(score_type);
	
	return ret_score;
}				

//////////////////////////
//for Graph<TAlign>

//vorsichtig! noch nicht richtig, bis jetzt nur ungapped exact matches...
template<typename TScore,typename TStringSet,typename TAlignment,typename TValue>
typename Value<TScore>::Type
getScore(TScore & score_type,
		 TStringSet & seqs,
		 Graph<TAlignment> & segment,
		 bool i_am_first, 
		 TValue pos_i,
		 TValue pos_j,
		 TValue len)
{
SEQAN_CHECKPOINT

	//typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,0);
	//typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,1);

	int i = 0;
	typename Value<TScore>::Type ret_score = 0;

//	while(i < len)
//	{
//		ret_score += score(score_type,label0[i],label1[i]);
////		ret_score += score(score_type,label0[i],label1[i]);
//		++i;
//	}
	ret_score = scoreMatch(score_type);
	ret_score *= len;

	return ret_score;
}				

//////////////////////////
//for Fragment

//get score for alignment starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true), if len1!=len2 then the refinement
//process was stopped (the cut is not exact)
template<typename TScore,typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec,typename TValue>
typename Value<TScore>::Type
getScore(TScore & score_type, 
		 TStringSet & seqs,
		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment, 
		 bool i_am_first, 
		 TValue pos_i, 
		 TValue pos_j,
		 TValue len1, 
		 TValue len2)
{
SEQAN_CHECKPOINT

	typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,seqs,sequenceId(segment,0));
	typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,seqs,sequenceId(segment,1));

	int i = 0;
	typename Value<TScore>::Type ret_score = 0;
	TValue len = (len1 < len2) ? len1 : len2;

	while(i < len)
	{
		ret_score += score(score_type,label0[i],label1[i]);
		++i;
	}
	len = (len1 > len2) ? len1 : len2;
	ret_score += (len - i) * scoreGapExtend(score_type);


	return ret_score;
}				



/////////////////////////////////////////////////////////////////////////////////////////////
//die nächsten beiden funktionen: für Fragmente und Score vom Typ Simple
//für den fall dass es keine mismatches innerhalb der segmente gibt und Score vom typ Simple ist
//TODO: müsste für einen bestimmten TFragSpec sein (Exact oder noMismatches)


//get score for alignment starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true), if len1!=len2 then the refinement
//process was stopped (the cut is not exact)
template<typename TScoreValue,typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec>
TScoreValue
getScore(Score<TScoreValue, Simple> & score_type,
		 TStringSet & seqs, 
		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment, 
		 bool i_am_first, 
		 TFragPos pos_i, 
		 TFragPos pos_j, 
		 TFragSize len1, 
		 TFragSize len2)
{
SEQAN_CHECKPOINT

	typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,seqs,sequenceId(segment,0));
	typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,seqs,sequenceId(segment,1));

	TScoreValue ret_score = 0;
	TFragSize len;
	if (len1 < len2) len = len1;
	else len = len2;

	if(len1 <= len2)
	{
		ret_score += len1 * scoreMatch(score_type);
		ret_score += (len2 - len1) * scoreGapExtend(score_type);
	}
	else{
		ret_score += len2 * scoreMatch(score_type);
		ret_score += (len1 - len2) * scoreGapExtend(score_type);
	}

	return ret_score;
}				

//get score for alignment of length len starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true)
template<typename TScoreValue,typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec>
TScoreValue
getScore(Score<TScoreValue, Simple> & score_type,
		 TStringSet & seqs,
		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment,
		 bool i_am_first, 
		 TFragPos pos_i,
		 TFragPos pos_j,
		 TFragSize len)
{
SEQAN_CHECKPOINT

	TScoreValue ret_score = 0;

	ret_score = scoreMatch(score_type);
	ret_score *= len;

	return ret_score;
}				



////////////////////////////////////////////////////////////////////////////////////////
//build refined alignment graph
////////////////////////////////////////////////////////////////////////////////////////
//nodes are numbered ascendingly:
//seq1   0  1  2  3  4 
//seq2   5  6  7  8  9 10
//seq3  11 12 13 14 15 
template<typename TValue,typename TAlignmentString,typename TScore,typename TSequence, typename TSetSpec,typename TIntervalTreeGraph,typename TPropertyMap, typename TAliGraph,typename TSeqMap>
void
_makeAlignmentGraphFromRefinedSegmentsSimple(String<std::set<TValue> > & all_nodes,
				   TAlignmentString & alis,
				   TScore & score_type,
				   StringSet<TSequence, TSetSpec> & seqs,
				   TSeqMap & seq_map,
				   String<TIntervalTreeGraph> & gs, 
				   String<TPropertyMap> & pms, 
				   TAliGraph & ali_g)
{
SEQAN_CHECKPOINT
	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Iterator<TAlignmentString>::Type TAliIterator;
	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
	typedef typename std::set<TValue>::iterator TSetIterator;

	std::cout << "making refined alignment graph...";
	clock_t start, finish1;
	double duration;
	start = clock();

	//make nodes
	//for each sequence
	for(int seq_i = 0; seq_i < (int) length(seqs); ++seq_i)
	{
		TSetIterator it = all_nodes[seq_i].begin();
		TSetIterator end_it = all_nodes[seq_i].end();
		TSetIterator next_it = it;
		if(next_it != end_it)
			++next_it;
		
		//a new node for each interval
		while(next_it != end_it)
		{
			TValue pos_i = *it;
			addVertex(ali_g, seq_i, pos_i, *next_it - pos_i);
			++it;
			++next_it;
		}

		all_nodes[seq_i].clear();
	}

	//make edges
	TAliIterator ali_it = begin(alis);
	TAliIterator ali_end = end(alis);


	//for each segment/fragement/alignment
	while(ali_it != ali_end)
	{
		//get sequence, begin position and end position
		TValue seq,begin_pos,end_pos;
		getSeqBeginAndEnd(*ali_it,seq_map,seq,begin_pos,end_pos,(TValue)0);
		
		//get the node represents the current interval (begin_pos until next_cut_pos or end_pos)
		TVertexDescriptor act_knot = findVertex(ali_g,seq,begin_pos);
		TValue act_pos = begin_pos;
	
		//for each interval that lies within the current segment/fragement/alignment
		while(act_pos < end_pos)
		{

			//get other sequence and projected position
			TValue seq_j,pos_j;
			bool i_am_first = _getOtherSequenceAndProject(*ali_it,seq_map,seq,act_pos,seq_j,pos_j);
			
			//find node that contains the projected position (pos_j)
			TVertexDescriptor vd = findVertex(ali_g,seq_j,pos_j);
			
			//if exact cut exists 
			if(fragmentBegin(ali_g,vd)==pos_j)
			{
				TValue score = getScore(score_type,seqs,*ali_it,i_am_first,act_pos,pos_j,fragmentLength(ali_g,act_knot));//,fragmentLength(ali_g,vd));
				addEdge(ali_g,act_knot,vd,score);
			}
			else //if the refinement was stopped here
			{
				std::cout << "probleme...\n";
			}

			//prepare for next interval
			act_pos += fragmentLength(ali_g,act_knot);
			act_knot = findVertex(ali_g,seq,act_pos);
		
		}
		++ali_it;
	}

	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	std::cout << "\ntook " << duration << " seconds.\n";


}



////////////////////////////////////////////////////////////////////////////////////////
//The big matchRefinement function that does everything: build interval trees, do the 
//refinement and construct a refined alignment graph
////////////////////////////////////////////////////////////////////////////////////////

template<typename TAlignmentString, typename TOutGraph, typename TSequence, typename TSetSpec, typename TScore>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				TScore & score_type,
				TOutGraph & ali_graph)
{
SEQAN_CHECKPOINT

	////////////////////////////////////////////////////////////////
	//typedefs

	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Iterator<TAlignmentString>::Type TAliIterator;
	typedef typename Size<TAlign>::Type TValue;
//	typedef TAlign* TCargo;
	typedef TValue TCargo;
	typedef IntervalAndCargo<int,TCargo> TInterval;
	typedef Graph<Directed<void,WithoutEdgeId> > TGraph;
	typedef IntervalTreeNode<TInterval> TNode;
	typedef String<TNode> TPropertyMap;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TCargo> TList;
	
	////////////////////////////////////////////////////////////////

	TValue numSequences = length(seq);
	//weird ID --> good ID map
	std::map<const void * ,int> seq_map;
	for(int i = 0; i < (int) numSequences; ++i)
		seq_map[id(seq[i])] = i;


	////////////////////////////////////////////////////////////////
	//build interval trees
	String<TGraph> gs;
	String<TPropertyMap> pms;

	createTreesForAllSequences(gs, pms, alis, seq, seq_map, numSequences, ComputeCenter());


	////////////////////////////////////////////////////////////////
	//do refinement
	std::cout <<"refining..."<<std::flush;
	clock_t start, finish1;
	double duration;
	start = clock();
	
	//all_nodes = set of all cut positions
	String<std::set<TValue> > all_nodes;
	resize(all_nodes,numSequences);

	//call function _refine for each startknoten
	TAliIterator ali_it = begin(alis);
	TAliIterator ali_end = end(alis);

	//for each segment/fragement/alignment
	while(ali_it != ali_end)
	{
		//for each of the two sequences
		for(TValue i = 0; i < 2; ++i)
		{
			TValue seq_i,begin_i,end_i;
			getSeqBeginAndEnd(*ali_it,seq_map,seq_i,begin_i,end_i,i);
	
			//refine begin
			if(all_nodes[seq_i].find(begin_i) == all_nodes[seq_i].end())
			{
				all_nodes[seq_i].insert(begin_i);
				_refine(begin_i, seq_i, seq_map, alis, gs,pms,all_nodes);//TStop());
			}
			//and end position
			if(all_nodes[seq_i].find(end_i) == all_nodes[seq_i].end())
			{
				all_nodes[seq_i].insert(end_i);
				_refine(end_i, seq_i, seq_map, alis, gs,pms,all_nodes);//TStop());
			}
		}	
		++ali_it;
	}

	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	std::cout << "\ntook " << duration << " seconds.\n";



	//for(int seq_i = 0; seq_i < length(seq); ++seq_i)
	//{
	//	std::set<TValue>::iterator it = all_nodes[seq_i].begin();
	//	std::set<TValue>::iterator end_it = all_nodes[seq_i].end();
	//
	//	while(it != end_it)
	//	{
	//		std::cout << *it << ",";
	//		++it;
	//	}
	//	std::cout << "\n";
	//}


	////////////////////////////////////////////////////////////////
	//build refined alignment graph
	_makeAlignmentGraphFromRefinedSegmentsSimple(all_nodes,alis,score_type,seq,seq_map,gs,pms,ali_graph);




}







////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// UNFINISHED PART //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////////////
//build refined alignment graph
////////////////////////////////////////////////////////////////////////////////////////
//nodes are numbered ascendingly:
//seq1   0  1  2  3  4 
//seq2   5  6  7  8  9 10
//seq3  11 12 13 14 15 
template<typename TValue,typename TAlignmentString,typename TScore,typename TSequence, typename TSetSpec,typename TIntervalTreeGraph,typename TPropertyMap, typename TAliGraph,typename TSeqMap>
void
_makeAlignmentGraphFromRefinedSegmentsUnclean(String<std::set<TValue> > & all_nodes,
				   TAlignmentString & alis,
				   TScore & score_type,
				   StringSet<TSequence, TSetSpec> & seqs,
				   TSeqMap & seq_map,
				   String<TIntervalTreeGraph> & gs, 
				   String<TPropertyMap> & pms, 
				   TAliGraph & ali_g)
{
SEQAN_CHECKPOINT
	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Iterator<TAlignmentString>::Type TAliIterator;
	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
	typedef typename std::set<TValue>::iterator TSetIterator;

	std::cout << "making refined alignment graph...";
	clock_t start, finish1;
	double duration;
	start = clock();

	//make nodes
	//for each sequence
	for(int seq_i = 0; seq_i < (int) length(seqs); ++seq_i)
	{
		TSetIterator it = all_nodes[seq_i].begin();
		TSetIterator end_it = all_nodes[seq_i].end();
		TSetIterator next_it = it;
		if(next_it != end_it)
			++next_it;
		
		//a new node for each interval
		while(next_it != end_it)
		{
			TValue pos_i = *it;
			addVertex(ali_g, seq_i, pos_i, *next_it - pos_i);
			++it;
			++next_it;
		}

		all_nodes[seq_i].clear();
	}

	//make edges
	TAliIterator ali_it = begin(alis);
	TAliIterator ali_end = end(alis);


	//for each segment/fragement/alignment
	while(ali_it != ali_end)
	{
		//get sequence, begin position and end position
		TValue seq,begin_pos,end_pos;
		getSeqBeginAndEnd(*ali_it,seq_map,seq,begin_pos,end_pos,0);
		
		//get the node that represents the current interval (begin_pos until next_cut_pos or end_pos)
		TVertexDescriptor act_knot = findVertex(ali_g,seq,begin_pos);
		if(begin_pos == fragmentBegin(ali_g,act_knot))


		TValue act_pos = begin_pos;
	
		//for each interval that lies within the current segment/fragement/alignment
		while(act_pos < end_pos)
		{
			//get other sequence and projected position
			TValue seq_j,pos_j;
			bool i_am_first = _getOtherSequenceAndProject(*ali_it,seq_map,seq,act_pos,seq_j,pos_j);
			
			//find node that contains the projected position (pos_j)
			TVertexDescriptor vd = findVertex(ali_g,seq_j,pos_j);
			
			//if exact cut exists 
			if(fragmentBegin(ali_g,vd)==pos_j)
			{
				TValue score = getScore(score_type,seqs,*ali_it,i_am_first,act_pos,pos_j,fragmentLength(ali_g,act_knot));//,fragmentLength(ali_g,vd));
				addEdge(ali_g,act_knot,vd,score);
			}
			else //if the refinement was stopped here
			{
				std::cout << "probleme...\n";
			}

			//prepare for next interval
			act_pos += fragmentLength(ali_g,act_knot);
			act_knot = findVertex(ali_g,seq,act_pos);
		
		}
		++ali_it;
	}

	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	std::cout << "\ntook " << duration << " seconds.\n";


}





template<typename TAlignmentString, typename TOutGraph, typename TSequence, typename TSetSpec, typename TScore>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				TScore & score_type,
				TOutGraph & ali_graph,
				unsigned int min_fragment_len)
{
SEQAN_CHECKPOINT

	////////////////////////////////////////////////////////////////
	//typedefs

	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Iterator<TAlignmentString>::Type TAliIterator;
	typedef typename Size<TAlign>::Type TValue;
//	typedef TAlign* TCargo;
	typedef TValue TCargo;
	typedef IntervalAndCargo<int,TCargo> TInterval;
	typedef Graph<Directed<void,WithoutEdgeId> > TGraph;
	typedef IntervalTreeNode<TInterval> TNode;
	typedef String<TNode> TPropertyMap;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TCargo> TList;
	typedef typename std::set<TValue>::iterator TSetIterator;

	////////////////////////////////////////////////////////////////

	TValue numSequences = length(seq);
	//weird ID --> good ID map
	std::map<const void * ,int> seq_map;
	for(int i = 0; i < (int) numSequences; ++i)
		seq_map[id(seq[i])] = i;


	////////////////////////////////////////////////////////////////
	//build interval trees
	String<TGraph> gs;
	String<TPropertyMap> pms;

	createTreesForAllSequences(gs, pms, alis, seq, seq_map, numSequences, ComputeCenter());


	////////////////////////////////////////////////////////////////
	//do refinement
	std::cout <<"refining..."<<std::flush;
	clock_t start, finish1;
	double duration;
	start = clock();
	
	//all_nodes = set of all cut positions
	String<std::set<TValue> > all_nodes;
	resize(all_nodes,numSequences);

	//call function _refine for each startknoten
	TAliIterator ali_it = begin(alis);
	TAliIterator ali_end = end(alis);

	//for each segment/fragement/alignment
	while(ali_it != ali_end)
	{
		//for each of the two sequences
		for(TValue i = 0; i < 2; ++i)
		{
			TValue seq_i,begin_i,end_i;
			getSeqBeginAndEnd(*ali_it,seq_map,seq_i,begin_i,end_i,i);
	
			//refine begin
			TSetIterator iter = all_nodes[seq_i].find(begin_i);		
			if(cutIsOk(all_nodes,seq_i,begin_i,iter,min_fragment_len))
			{
				all_nodes[seq_i].insert(begin_i);
				_refine(begin_i, seq_i, seq_map, alis, gs, pms, all_nodes, min_fragment_len);//TStop());
			}
			//and end position
			iter = all_nodes[seq_i].find(end_i);		
			if(cutIsOk(all_nodes,seq_i,end_i,iter,min_fragment_len))
			{
				all_nodes[seq_i].insert(end_i);
				_refine(end_i, seq_i, seq_map, alis, gs, pms, all_nodes, min_fragment_len);//TStop());
			}
		}	
		++ali_it;
	}

	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	std::cout << "\ntook " << duration << " seconds.\n";



	//for(int seq_i = 0; seq_i < length(seq); ++seq_i)
	//{
	//	std::set<TValue>::iterator it = all_nodes[seq_i].begin();
	//	std::set<TValue>::iterator end_it = all_nodes[seq_i].end();
	//
	//	while(it != end_it)
	//	{
	//		std::cout << *it << ",";
	//		++it;
	//	}
	//	std::cout << "\n";
	//}


	////////////////////////////////////////////////////////////////
	//build refined alignment graph
	_makeAlignmentGraphFromRefinedSegmentsUnclean(all_nodes,alis,score_type,seq,seq_map,gs,pms,ali_graph);




}
	


}
#endif //#ifndef SEQAN_HEADER_...
