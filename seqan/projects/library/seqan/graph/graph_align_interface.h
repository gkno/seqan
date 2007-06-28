#ifndef SEQAN_HEADER_GRAPH_ALIGN_INTERFACE_H
#define SEQAN_HEADER_GRAPH_ALIGN_INTERFACE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Function.globalAlignment:
..summary:Computes the best global alignment of the two sequences.
..cat:Alignments
..signature:
globalAlignment(strSet, score, tag)
globalAlignment(graph, score, tag)
globalAlignment(file, strSet, score, tag)
..param.strSet:A string set with 2 sequences.
...type:Class.StringSet
...remarks: If an alignment graph is used that graph must contain a string set with two sequences
..param.graph:The alignment graph having 2 sequences.
...type:Class.Graph Alignment
..param.file:A file stream or std::cout to write a textual alignment
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.tag:A tag indicating the alignment algorithm to use
...remarks:NeedlemanWunsch, Gotoh, Hirschberg, or MyersBitVector.
...remarks:If MyersBitVector is used leave out score because the scoring is fixed in this algorithm.
..returns:The maximum score of the best global alignment.
*/
template<typename TFile, typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
TScoreValue
globalAlignment(TFile& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(file,str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
TScoreValue
globalAlignment(TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TTag>
TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _globalAlignment(g,stringSet(g),sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet>
unsigned int
globalAlignment(TStringSet const& str,
				MyersBitVector)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(str,Score<unsigned int>(0,1,1,0),MyersBitVector());
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.localAlignment:
..summary:Computes the best local alignment of two sequences.
..cat:Alignments
..signature:
localAlignment(strSet, score, tag)
localAlignment(graph, score, tag)
localAlignment(file, strSet, score, tag)
..param.strSet:A string set with 2 sequences.
...type:Class.StringSet
...remarks: If an alignment graph is used that graph must contain a string set with two sequences
..param.graph:The alignment graph having 2 sequences.
...type:Class.Graph Alignment
..param.file:A file stream or std::cout to write a textual alignment
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.tag:A tag indicating the alignment algorithm to use
...remarks:SmithWaterman
..returns:The maximum score of the best local alignment.
*/
template<typename TFile, typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
TScoreValue
localAlignment(TFile& file,
			   TStringSet const& str,
			   Score<TScoreValue, TSpec> const& sc,
			   TTag)
{
	SEQAN_CHECKPOINT
	return _localAlignment(file,str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
TScoreValue
localAlignment(TStringSet const& str,
			   Score<TScoreValue, TSpec> const& sc,
			   TTag)
{
	SEQAN_CHECKPOINT
	return _localAlignment(str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TTag>
TScoreValue
localAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
			   Score<TScoreValue, TSpec2> const& sc,
			   TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _localAlignment(g,stringSet(g),sc,TTag());
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
