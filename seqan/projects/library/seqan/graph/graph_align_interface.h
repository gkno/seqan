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
globalAlignment(g, score, tag)
globalAlignment(file, str, score, tag)
..param.g:The alignment graph having 2 sequences.
...type:Class.Graph Alignment
..param.str:A string set with 2 sequences.
...type:Class.StringSet
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.tag:A tag indicating the alignment algorithm to use
...remarks:Either NeedlemanWunsch or Gotoh.
..returns:The score value of the best scoring global alignment.
*/
template<typename TFile, typename TStringSet, typename TScoreValue, typename TTag>
TScoreValue
globalAlignment(TFile& file,
				TStringSet const& str,
				Score<TScoreValue, Simple> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(file,str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TTag>
TScoreValue
globalAlignment(TStringSet const& str,
				Score<TScoreValue, Simple> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TTag>
TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, Simple> const& sc,
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

template<typename TStringSet, typename TCargo, typename TSpec, typename TTag>
unsigned int
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
		TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _globalAlignment(g,stringSet(g),Score<unsigned int>(0,1,1,0),Hirschberg());
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
