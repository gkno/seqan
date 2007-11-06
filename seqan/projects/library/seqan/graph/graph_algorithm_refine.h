 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_TEMP_REFINE_H
#define SEQAN_HEADER_GRAPH_TEMP_REFINE_H

namespace SEQAN_NAMESPACE_MAIN
{


 
///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functions for Align<TSource,TSpec>
//project onto other sequence 
template<typename TSource,typename TSpec,typename TStringSet,typename TValue,typename TMap>
bool
_getOtherSequenceAndProject(Align<TSource,TSpec> & segment, 
						    TStringSet & seqs,
							TMap & seq_map, 
						   TValue seq_i, 
						   TValue node_i, 
						   TValue & seq_j, 
						   TValue & node_j)
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
template<typename TAlignment,typename TStringSet,typename TValue, typename TMap>
bool
_getOtherSequenceAndProject(Graph<TAlignment> & segment, 
							TStringSet & seqs,
                            TMap & seq_map, 
							TValue seq_i, 
							TValue pos_i, 
							TValue & seq_j, 
							TValue & pos_j)
{
SEQAN_CHECKPOINT

	seq_i = positionToId(seqs, seq_i);
	pos_j = getProjectedPosition(segment,seq_i, pos_i);
	
	if(seq_i == sequenceId(segment,0))
	{
		seq_j = sequenceId(segment,1);
		seq_j = idToPosition(seqs, seq_j);
		return true;
	}
	else
	{
		seq_j = sequenceId(segment,0);
		seq_j = idToPosition(seqs, seq_j);
		return false;
	}

	//typedef Graph<TAlignment> TGraph;
	//typedef typename Iterator<TGraph, OutEdgeIterator>::Type TEdgeIterator;
	//typedef typename VertexDescriptor<TGraph>::Type TVertex;
	//
	//seq_i = sequenceId(segment,seq_i);
	//TVertex vd1 = findVertex(segment,seq_i,pos_i);
	//TEdgeIterator it(segment,vd1);
	//TVertex vd2 = targetVertex(it);
	//seq_j = sequenceId(segment,vd2);
	//pos_j = fragmentBegin(segment,vd2) + (pos_i - fragmentBegin(segment,vd1));

	//if(vd1 == 0)
	//	return true;
	//else
	//	return false;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functions for Fragments
//project onto other sequence for Graph<Alignment>
template<typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec,typename TValue, typename TMap>
bool
_getOtherSequenceAndProject(Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment,
							TStringSet & seqs,
						   TMap &,
						   TValue seq_i,
						   TValue pos_i,
						   TValue & seq_j,
						   TValue & pos_j)
{
SEQAN_CHECKPOINT
	seq_i = positionToId(seqs, seq_i);
	getProjectedPosition(segment,seq_i, pos_i, seq_j, pos_j);
	
	if(seq_i == sequenceId(segment,0))
	{
		seq_j = sequenceId(segment,1);
		seq_j = idToPosition(seqs, seq_j);
		return true;
	}
	else
	{
		seq_j = sequenceId(segment,0);
		seq_j = idToPosition(seqs, seq_j);
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
				  std::map<const void * ,int> &, 
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
template<typename TValue, typename TAlignmentString, typename TStringSet, typename TGraph, typename TPropertyMap,typename TSeqMap>
void
_refine(TValue node_i, 
	 TValue seq_i, 
	 TStringSet & seqs,
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
		_getOtherSequenceAndProject(alis[*segment_it],seqs,seq_map,seq_i,node_i,seq_j,node_j);
		
		typename std::set<TValue>::iterator iter;
		iter = all_nodes[seq_j].find(node_j);
		
		//if node does not exist yet ---> insert and continue cutting
		if(iter == all_nodes[seq_j].end())
		{
			//TODO Abbruch: if(!stop(all_nodes,seq_check,node_j,seq_j,TStop()))
			all_nodes[seq_j].insert(node_j);
			_refine(node_j,seq_j,seqs,seq_map,alis,gs,pms,all_nodes);
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
	
	//cut already exists
	if(iter != all_nodes[seq_i].end())
		return false;
//	if(min_len == 0)
//		return true;

	typename std::set<TValue>::iterator tmp_iter = all_nodes[seq_i].upper_bound(pos_i);
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
template<typename TValue, typename TAlignmentString, typename TStringSet,typename TGraph, typename TPropertyMap,typename TSeqMap>
void
_refine(TValue node_i, 
	 TValue seq_i, 
	 TStringSet & seqs,
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
		_getOtherSequenceAndProject(alis[*segment_it],seqs,seq_map,seq_i,node_i,seq_j,node_j);
		
		typename std::set<TValue>::iterator iter;
		iter = all_nodes[seq_j].find(node_j);
		
		//if node does not exist yet ---> insert and continue cutting
		if(cutIsOk(all_nodes,seq_j,node_j,iter,min_len))
		{
			all_nodes[seq_j].insert(node_j);
			_refine(node_j,seq_j,seqs,seq_map,alis,gs,pms,all_nodes,min_len);
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
template<typename TInterval, typename TStringSet, typename TAlignmentString, typename TSeqMap>
void
_buildIntervalsForAllSequences(TAlignmentString & alis, 
							   String<String<TInterval> > & intervals, 
	   						   TStringSet & seqs,
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
		seq_i = idToPosition(seqs, seq_i);

		//and append the interval (ali_begin, ali_end) with cargo ali* to the list of intervals of seq_i
		appendValue(intervals[seq_i],IntervalAndCargo<TValue,TCargo>(begin_,end_,ali_counter)); 
	
		//get the second sequence (and its begin and end) that takes part in the alignment (seq_i)
		getSeqBeginAndEnd(*ali_it,seq_map,seq_i,begin_,end_,1);
		seq_i = idToPosition(seqs, seq_i);
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
	
	//std::cout <<"create interval trees...";
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
	_buildIntervalsForAllSequences(alis,intervals,seqs,seq_map);
	
	TValue i = 0;
	
	while(i < numSequences)
	{
		//std::cout << (numSequences-i) <<" more ("<<length(intervals[i])<<" intervals)... ";
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
	//std::cout << "\ntook " << duration << " seconds.\n";



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

	typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,stringSet(segment)[0]);
	typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,stringSet(segment)[1]);
	//typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,0);
	//typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,1);

	int i = 0;
	typename Value<TScore>::Type ret_score = 0;

	while(i < len)
	{
		ret_score += score(score_type,label0[i],label1[i]);
		++i;
	}
	//ret_score = scoreMatch(score_type);
	//ret_score *= len;

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
template<typename TScoreValue, typename TSpec, typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec>
TScoreValue
getScore(Score<TScoreValue, TSpec> & score_type,
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
template<typename TScoreValue, typename TSpec, typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec>
TScoreValue
getScore(Score<TScoreValue, TSpec> & score_type,
		 TStringSet & seqs,
		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment,
		 bool, 
		 TFragPos pos_i,
		 TFragPos pos_j,
		 TFragSize len)
{
SEQAN_CHECKPOINT
	typedef typename Infix<typename Value<TStringSet>::Type>::Type TSegmentLabel;
	TSegmentLabel label0 = label(segment,seqs, sequenceId(segment, 0));
	TSegmentLabel label1 = label(segment,seqs, sequenceId(segment, 1));

	typename Iterator<TSegmentLabel>::Type label_it0 = begin(label0) + (pos_i - fragmentBegin(segment,sequenceId(segment,0)));
	typename Iterator<TSegmentLabel>::Type label_it1 = begin(label1) + (pos_j - fragmentBegin(segment,sequenceId(segment,1)));

	int i = 0;
	TScoreValue ret_score = 0;

	while(i < (int) len)
	{
		ret_score += score(score_type,*label_it0,*label_it1);
		++label_it0;
		++label_it1;
		++i;
	}

	return ret_score;

}				


struct FakeScore;

template<typename TValue>
class Score<TValue,FakeScore>
{
public:
	Score() {
	}

	Score(Score const & other) {
	}

	Score & operator = (Score const & other) {
		return *this;
	}
};

template <typename TValue>
inline TValue
scoreGapExtend(Score<TValue,FakeScore> & me)
{
	return 1;
}


template <typename TValue>
inline TValue const
scoreGapExtend(Score<TValue,FakeScore> const & me)
{
	return 1;
}

template <typename TValue>
inline TValue
scoreGapOpen(Score<TValue,FakeScore> & me)
{
	return 1;
}
template <typename TValue>
inline TValue const
scoreGapOpen(Score<TValue,FakeScore> const & me)
{
	return 1;
}

template <typename TValue, typename T>
inline TValue
score(Score<TValue,FakeScore> & me,
	  T const & left,
	  T const & right)
{
	return 1;
}

template <typename TValue>
struct Value< Score<TValue, FakeScore> >
{
	typedef TValue Type;
};


//fake score function 
template<typename TScoreValue,typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec>
TScoreValue
getScore(Score<TScoreValue,FakeScore> &,
		 TStringSet &,
		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> &,
		 bool, 
		 TFragPos,
		 TFragPos,
		 TFragSize)
{
SEQAN_CHECKPOINT
	return 1;
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
				   String<TIntervalTreeGraph> &, 
				   String<TPropertyMap> &, 
				   TAliGraph & ali_g)
{
SEQAN_CHECKPOINT
	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Iterator<TAlignmentString>::Type TAliIterator;
	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAliGraph>::Type TEdgeDescriptor;
	typedef typename std::set<TValue>::iterator TSetIterator;
	typedef typename Cargo<TAliGraph>::Type TCargo;

	//std::cout << "making refined alignment graph...";
	clock_t start, finish1;
	double duration;
	start = clock();

	
	//make nodes
	//for each sequence
	for(int seq_i = 0; seq_i < (int) length(seqs); ++seq_i)
	{
		TValue seq_id = positionToId(stringSet(ali_g), seq_i);
		TSetIterator it = all_nodes[seq_i].begin();
		TSetIterator end_it = all_nodes[seq_i].end();
		TSetIterator next_it = it;
		if(next_it != end_it)
			++next_it;
		
		//first unaligned node
		if(it != end_it && *it != 0)
			addVertex(ali_g, seq_id, 0, *it);

		//a new node for each interval
		while(next_it != end_it)
		{
			TValue pos_i = *it;
			addVertex(ali_g, seq_id, pos_i, *next_it - pos_i); 
			++it;
			++next_it;
		}

		//last unaligned node
		if(it !=end_it && *it<length(seqs[seq_i])-1)
			addVertex(ali_g, seq_id, *it, (length(seqs[seq_i])-1) - *it);


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
		seq = idToPosition(seqs,seq);		

		TValue act_pos = begin_pos;
	
		//for each interval that lies within the current segment/fragement/alignment
		while(act_pos < end_pos)
		{

			//get other sequence and projected position
			TValue seq_j,pos_j;
			bool i_am_first = _getOtherSequenceAndProject(*ali_it,seqs,seq_map,seq,act_pos,seq_j,pos_j);
			seq_j = positionToId(seqs, seq_j);
			//find node that contains the projected position (pos_j)
			TVertexDescriptor vd = findVertex(ali_g, seq_j, pos_j);
			
			SEQAN_TASSERT(fragmentBegin(ali_g,vd)==pos_j)
			typename Value<TScore>::Type score = getScore(score_type,seqs,*ali_it,i_am_first,act_pos,pos_j,fragmentLength(ali_g,act_knot));//,fragmentLength(ali_g,vd));
			//this needs to be generalized (makes sense for positive scores only)
			if(score > 0)
			{
				if (findEdge(ali_g, act_knot, vd) == 0) addEdge(ali_g,act_knot,vd,score);
				else {
					TEdgeDescriptor ed = findEdge(ali_g, act_knot, vd);
					if((TCargo) score > getCargo(ed))
						assignCargo(ed, score);
					// ToDo: Adapt score of the edge
				}
			}
			//prepare for next interval
			act_pos += fragmentLength(ali_g,act_knot);
			act_knot = findVertex(ali_g,positionToId(seqs, seq),act_pos);
		
		}
		++ali_it;
	}

	//std::cout << "check\n";
	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";


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
	//std::cout <<"refining..."<<std::flush;
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
			seq_i = idToPosition(seq,seq_i);
			
			//refine begin
			if(all_nodes[seq_i].find(begin_i) == all_nodes[seq_i].end())
			{
				all_nodes[seq_i].insert(begin_i);
				_refine(begin_i, seq_i, seq, seq_map, alis, gs,pms,all_nodes);//TStop());
			}
			//and end position
			if(all_nodes[seq_i].find(end_i) == all_nodes[seq_i].end())
			{
				all_nodes[seq_i].insert(end_i);
				_refine(end_i, seq_i, seq, seq_map, alis, gs,pms,all_nodes);//TStop());
			}
		}	
		++ali_it;
	}

	finish1 = clock();
	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";



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



template<typename TAlignmentString, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				TOutGraph & ali_graph)
{
SEQAN_CHECKPOINT

	Score<int,FakeScore > score_type;
	matchRefinement(alis,seq,score_type,ali_graph);

}





////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// UNFINISHED PART //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//template<typename TAliGraph, typename TVertexDescriptor,typename TValue>
//TValue  
//_getClosestRefinedNeighbor(TAliGraph & ali_g,
//						   TVertexDescriptor & vd,
//						   TValue seq,
//						   TValue pos)
//{
//SEQAN_CHECKPOINT
//
//	if(pos-fragmentBegin(ali_g,vd) < fragmentBegin(ali_g,vd)+fragmentLength(ali_g,vd)-pos)
//		return fragmentBegin(ali_g,vd);
//	else
//		return fragmentBegin(ali_g,vd) + fragmentLength(ali_g,vd);
//
//}
//
//template<typename TAliGraph,typename TValue>
//void
//_getCutEndPos(TAliGraph & ali_g, 
//			  typename VertexDescriptor<TAliGraph>::Type & end_knot,
//			  TValue seq,
//			  TValue end_pos,
//			  TValue & cut_end_pos)
//{
//SEQAN_CHECKPOINT
//
//	end_knot = findVertex(ali_g,seq,end_pos-1);//end_pos1 is the first position of the next node
//	if(end_pos == fragmentBegin(ali_g,end_knot) + fragmentBegin(ali_g,end_knot))
//		cut_end_pos = end_pos;
//	else
//	{
//		cut_end_pos = _getClosestRefinedNeighbor(ali_g,end_knot,seq,end_pos);
//		end_knot =  findVertex(ali_g,seq,cut_end_pos-1);
//		SEQAN_TASSERT(cut_end_pos == fragmentBegin(ali_g,end_knot)+fragmentLength(ali_g,end_knot))
//	}
//		
//
//}
//
//template<typename TAliGraph,typename TValue>
//void
//_getCutBeginPos(TAliGraph & ali_g, 
//			  typename VertexDescriptor<TAliGraph>::Type & act_knot,
//			  TValue seq,
//			  TValue act_pos,
//			  TValue & cut_act_pos)
//{
//SEQAN_CHECKPOINT
//	
//	act_knot = findVertex(ali_g,seq,act_pos);
//	//if completely refined
//	if(act_pos == fragmentBegin(ali_g,act_knot))
//		cut_act_pos = act_pos;
//	else //if incompletely refined
//	{
//		cut_act_pos = _getClosestRefinedNeighbor(ali_g,act_knot,seq,act_pos);
//		act_knot =  findVertex(ali_g,seq,cut_act_pos);
//		SEQAN_TASSERT(cut_act_pos == fragmentBegin(ali_g,act_knot))
//	}
//
//
//}
//
//
////template<typename TScore, typename TStringSet, typename TValue>
////typename Value<TScore>::Type 
////_getLeftRestScore(TScore & score_type,
////						  TStringSet & seqs,
////						  TValue seq1,
////						  TValue seq2,
////						  TValue act_pos1,
////						  TValue cut_act_pos1,
////						  TValue act_pos2,
////						  TValue cut_act_pos2)
////{
////SEQAN_CHECKPOINT
////					  
////	typedef typename Iterator<typename Value<TStringSet>::Type>::Type TStringIterator;
////
////	if(cut_act_pos1 == act_pos1 && cut_act_pos2 == act_pos2)
////		return 0;
////	
////	typename Value<TScore>::Type retscore = 0;
////
////	if(cut_act_pos1 <= act_pos1 && cut_act_pos2 <= act_pos2)
////	{
////		TStringIterator it1 = seqs[seq1][--act_pos1];
////		TStringIterator end_it1 = seqs[seq1][cut_act_pos1];
////		TStringIterator it2 = seqs[seq2][--act_pos2];
////		TStringIterator end_it2 = seqs[seq2][cut_act_pos2];
////
////		while(it1>=end_it1 && it2>=end_it2)
////			retscore += score(score_type,*it1--,*it2--);
////		retscore += (it1-end_it1)*scoreGapExtend(score_type);
////		retscore += (it2-end_it2)*scoreGapExtend(score_type);
////
////	}
////}
//
//////////////////////////////////////////////////////////////////////////////////////////
////build refined alignment graph
//////////////////////////////////////////////////////////////////////////////////////////
////nodes are numbered ascendingly:
////seq1   0  1  2  3  4 
////seq2   5  6  7  8  9 10
////seq3  11 12 13 14 15 
//template<typename TValue,typename TAlignmentString,typename TScore,typename TSequence, typename TSetSpec,typename TIntervalTreeGraph,typename TPropertyMap, typename TAliGraph,typename TSeqMap>
//void
//_makeAlignmentGraphFromRefinedSegmentsUnclean(String<std::set<TValue> > & all_nodes,
//				   TAlignmentString & alis,
//				   TScore & score_type,
//				   StringSet<TSequence, TSetSpec> & seqs,
//				   TSeqMap & seq_map,
//				   String<TIntervalTreeGraph> & gs, 
//				   String<TPropertyMap> & pms, 
//				   TAliGraph & ali_g)
//{
//SEQAN_CHECKPOINT
//	typedef typename Value<TAlignmentString>::Type TAlign;
//	typedef typename Iterator<TAlignmentString>::Type TAliIterator;
//	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
//	typedef typename std::set<TValue>::iterator TSetIterator;
//
//	std::cout << "making refined alignment graph...";
//	clock_t start, finish1;
//	double duration;
//	start = clock();
//
//	//make nodes
//	//for each sequence
//	for(int seq_i = 0; seq_i < (int) length(seqs); ++seq_i)
//	{
//		TSetIterator it = all_nodes[seq_i].begin();
//		TSetIterator end_it = all_nodes[seq_i].end();
//		TSetIterator next_it = it;
//		if(next_it != end_it)
//			++next_it;
//		
//		//a new node for each interval
//		while(next_it != end_it)
//		{
//			TValue pos_i = *it;
//			addVertex(ali_g, seq_i, pos_i, *next_it - pos_i);
//			++it;
//			++next_it;
//		}
//
//		all_nodes[seq_i].clear();
//	}
//
//	//make edges
//	TAliIterator ali_it = begin(alis);
//	TAliIterator ali_end = end(alis);
//
//
//	//for each segment/fragment/alignment
//	while(ali_it != ali_end)
//	{
//		//get first sequence that takes part in the alignment + boundaries of the ali
//		TValue seq1,begin_pos1,end_pos1;
//		getSeqBeginAndEnd(*ali_it,seq_map,seq1,begin_pos1,end_pos1,(TValue)0);
//
//		//get the last node that is within the current ali
//		TVertexDescriptor end_knot1;
//		TValue cut_end_pos1;
//		_getCutEndPos(ali_g,end_knot1,seq1,end_pos1,cut_end_pos1);
//	
//		//get the node that represents the current interval (begin_pos until next_cut_pos or end_pos)
//		TVertexDescriptor act_knot1;
//		TValue cut_act_pos1,act_pos1;
//		act_pos1 = begin_pos1;
//		_getCutBeginPos(ali_g,act_knot1,seq1,act_pos1,cut_act_pos1);
//
//		TValue act_end_pos1 = cut_act_pos1 + fragmentLength(ali_g,act_knot1);
//
//		//walk through cuts on the first sequence
////		while (act_end_pos1 <= cut_end_pos1)
//		while (true)
//		{
//			//get other sequence and projected position
//			TValue seq2,act_pos2;
//			_getOtherSequenceAndProject(*ali_it,seq_map,seq1,act_pos1,seq2,act_pos2);
//		
//			//get node that corresponds to that position
//			TVertexDescriptor act_knot2;
//			TValue cut_act_pos2;
//			_getCutBeginPos(ali_g,act_knot2,seq2,act_pos2,cut_act_pos2);
//
//			//corresponding end on seq2 (there might be more than one node on seq2 that corresponds
//			//to the same interval (=node) on seq1)
//			TValue act_end_pos2;
//			_getOtherSequenceAndProject(*ali_it,seq_map,seq1,act_end_pos1-1,seq2,act_end_pos2);
//			++act_end_pos2;
//			TVertexDescriptor act_end_knot2;
//			TValue cut_act_end_pos2;
//			_getCutEndPos(ali_g,act_end_knot2,seq2,act_end_pos2,cut_act_end_pos2);
//			
//			if(cut_act_pos2 == cut_act_end_pos2)
//				break;
//			while(true)
//			{
//				//should at the moment return score for:
//				//
//				//seq1 = ....cr...rc....
//				//            ||||||
//				//seq2 = ...c.r...rc....
//				//bzw
//				//seq1 = ..cr.....x....   man will aber nur    ..cr......x....
//				//          |||||||-							 ---||||||  
//				//seq2 = ...r.c...rc... 					   ...r.c...rc....
//				typename Value<TScore>::Type score = 0;
//				score = getScore(score_type,seqs,*ali_it,i_am_first,act_pos1,act_pos2,act_end_pos1-act_pos1,cut_act_end_pos2);
//				//add score for
//				//
//				//seq1 = ...-cr....x....
//				//          ||
//				//seq2 = ...c.r...rc....
////					score += getLeftRestScore(score_type,seqs,seq1,seq2,act_pos1,cut_act_pos1,act_pos2,cut_act_pos2);
//
//				if(score > 0)
//					if(findEdge(ali_g,act_knot1,act_knot2)==0)
//						addEdge(ali_g,act_knot1,act_knot2,score);
//				
//				if(act_knot2==act_end_knot2)
//					break;
//
//				act_pos2 = cut_act_pos2 + fragmentLength(ali_g,act_knot2);
//				_getCutBeginPos(ali_g,act_knot2,seq2,act_pos2,cut_act_pos2);
//			}
//
//			if(act_knot1 == end_knot1)
//				break;
//
//			act_pos1 = act_end_pos1;
//			act_knot1 = findVertex(ali_g,seq1,act_pos1);
//			cut_act_pos1 = act_pos1;
//			act_end_pos1 = cut_act_pos1 + fragmentLength(ali_g,act_knot1);
//		}
//		++ali_it;
//	}
//
//	finish1 = clock();
//	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
//	std::cout << "\ntook " << duration << " seconds.\n";
//
//
//}
//
//
//
//
//
//template<typename TAlignmentString, typename TOutGraph, typename TSequence, typename TSetSpec, typename TScore>
//void
//matchRefinement(TAlignmentString & alis,
//				StringSet<TSequence, TSetSpec> & seq, 
//				TScore & score_type,
//				TOutGraph & ali_graph,
//				unsigned int min_fragment_len)
//{
//SEQAN_CHECKPOINT
//
//	////////////////////////////////////////////////////////////////
//	//typedefs
//
//	typedef typename Value<TAlignmentString>::Type TAlign;
//	typedef typename Iterator<TAlignmentString>::Type TAliIterator;
//	typedef typename Size<TAlign>::Type TValue;
////	typedef TAlign* TCargo;
//	typedef TValue TCargo;
//	typedef IntervalAndCargo<int,TCargo> TInterval;
//	typedef Graph<Directed<void,WithoutEdgeId> > TGraph;
//	typedef IntervalTreeNode<TInterval> TNode;
//	typedef String<TNode> TPropertyMap;
//	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
//	typedef String<TCargo> TList;
//	typedef typename std::set<TValue>::iterator TSetIterator;
//
//	////////////////////////////////////////////////////////////////
//
//	TValue numSequences = length(seq);
//	//weird ID --> good ID map
//	std::map<const void * ,int> seq_map;
//	for(int i = 0; i < (int) numSequences; ++i)
//		seq_map[id(seq[i])] = i;
//
//
//	////////////////////////////////////////////////////////////////
//	//build interval trees
//	String<TGraph> gs;
//	String<TPropertyMap> pms;
//
//	createTreesForAllSequences(gs, pms, alis, seq, seq_map, numSequences, ComputeCenter());
//
//
//	////////////////////////////////////////////////////////////////
//	//do refinement
//	std::cout <<"refining..."<<std::flush;
//	clock_t start, finish1;
//	double duration;
//	start = clock();
//	
//	//all_nodes = set of all cut positions
//	String<std::set<TValue> > all_nodes;
//	resize(all_nodes,numSequences);
//
//	//call function _refine for each startknoten
//	TAliIterator ali_it = begin(alis);
//	TAliIterator ali_end = end(alis);
//
//	//for each segment/fragement/alignment
//	while(ali_it != ali_end)
//	{
//		//for each of the two sequences
//		for(TValue i = 0; i < 2; ++i)
//		{
//			TValue seq_i,begin_i,end_i;
//			getSeqBeginAndEnd(*ali_it,seq_map,seq_i,begin_i,end_i,i);
//			seq_i = idToPosition(seq,seq_i);
//
//			//refine begin
//			TSetIterator iter = all_nodes[seq_i].find(begin_i);		
//			if(cutIsOk(all_nodes,seq_i,begin_i,iter,min_fragment_len))
//			{
//				all_nodes[seq_i].insert(begin_i);
//				_refine(begin_i, seq_i, seq, seq_map, alis, gs, pms, all_nodes, min_fragment_len);//TStop());
//			}
//			//and end position
//			iter = all_nodes[seq_i].find(end_i);		
//			if(cutIsOk(all_nodes,seq_i,end_i,iter,min_fragment_len))
//			{
//				all_nodes[seq_i].insert(end_i);
//				_refine(end_i, seq_i, seq, seq_map, alis, gs, pms, all_nodes, min_fragment_len);//TStop());
//			}
//		}	
//		++ali_it;
//	}
//
//	finish1 = clock();
//	duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
//	std::cout << "\ntook " << duration << " seconds.\n";
//	
//	int insgesamt = 0;
//	for(TValue ii = 0; ii < length(seq); ++ii)
//		insgesamt += all_nodes[ii].size();
//	std::cout << "Number of cuts: " << insgesamt << "\n";
//
//	//for(int seq_i = 0; seq_i < length(seq); ++seq_i)
//	//{
//	//	typename std::set<TValue>::iterator it = all_nodes[seq_i].begin();
//	//	typename std::set<TValue>::iterator end_it = all_nodes[seq_i].end();
//	//
//	//	while(it != end_it)
//	//	{
//	//		std::cout << *it << ",";
//	//		++it;
//	//	}
//	//	std::cout << "\n";
//	//}
//
//
//	////////////////////////////////////////////////////////////////
//	//build refined alignment graph
//	_makeAlignmentGraphFromRefinedSegmentsUnclean(all_nodes,alis,score_type,seq,seq_map,gs,pms,ali_graph);
//
//
//
//
//}
	


}
#endif //#ifndef SEQAN_HEADER_...
