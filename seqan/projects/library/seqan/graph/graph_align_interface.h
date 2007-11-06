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
 globalAlignment(g, score_type, Gotoh());
globalAlignment( [strSet | graph | align, strSet], score [, align_config], tag)
..param.strSet:A string set with 2 sequences.
...type:Class.StringSet
..param.graph:An alignment graph containing 2 sequences.
...type:Class.Graph Alignment
..param.align:Different alignment data structures, e.g., a file stream or std::cout for a textual alignment or a FragmentString for all the matches.
..param.score:The score values to be used for computing the alignment.
...type:Class.Score
..param.align_config:Alignment configuration options.
...type:Class.AlignConfig
...remarks:The class AlignConfig has four boolean parameters, i.e., TTop, TLeft, TRight, and TBottom.
If TTop is true the first row of the DP Matrix is initialized with 0's. If TLeft is true the first
column is initialized with 0's. If TRight is true, the maximum is search in the last column. If TBottom
is true, the maximum is searched in the last row. All options can be combined in all possible ways.
The Hirschberg algorithm currently doesn't support this kind of configuration.
..param.tag:A tag indicating the alignment algorithm to use
...remarks:NeedlemanWunsch, Gotoh, or Hirschberg.
..returns:The maximum score of the best global alignment.
*/
template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
inline TScoreValue
globalAlignment(TAlign& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return globalAlignment(file,str,sc, AlignConfig<>(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TAlignConfig, typename TTag>
inline TScoreValue
globalAlignment(TAlign& file,
				TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TAlignConfig const,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(file,str,sc, TAlignConfig(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
inline TScoreValue
globalAlignment(TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	return globalAlignment(str,sc,AlignConfig<>(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TAlignConfig, typename TTag>
inline TScoreValue
globalAlignment(TStringSet const& str,
				Score<TScoreValue, TSpec> const& sc,
				TAlignConfig const,
				TTag)
{
	SEQAN_CHECKPOINT
	return _globalAlignment(str,sc, TAlignConfig(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TTag>
inline TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return globalAlignment(g,stringSet(g),sc, AlignConfig<>(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TAlignConfig, typename TTag>
inline TScoreValue
globalAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
				Score<TScoreValue, TSpec2> const& sc,
				TAlignConfig const,
				TTag)
{
	SEQAN_CHECKPOINT
	clearVertices(g);
	return _globalAlignment(g,stringSet(g),sc, TAlignConfig(), TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet>
inline unsigned int
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
template<typename TAlign, typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
inline TScoreValue
localAlignment(TAlign& file,
			   TStringSet const& str,
			   Score<TScoreValue, TSpec> const& sc,
			   TTag)
{
	return _localAlignment(file,str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TScoreValue, typename TSpec, typename TTag>
inline TScoreValue
localAlignment(TStringSet const& str,
			   Score<TScoreValue, TSpec> const& sc,
			   TTag)
{
	SEQAN_CHECKPOINT
	return _localAlignment(str,sc,TTag());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TTag>
inline TScoreValue
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
