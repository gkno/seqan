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

#ifndef SEQAN_HEADER_GRAPH_IMPL_ORACLE_H
#define SEQAN_HEADER_GRAPH_IMPL_ORACLE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Oracle
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Oracle:
..cat:Graph
..general:Class.Graph
..summary:A factor oracle.
..remarks:A factor oracle is a special automaton and thus, it is not implemented in its own class.
It solely provides create functions where based upon a string an oracle is created.
..signature:Graph<Automaton<TAlphabet, TCargo, TSpec> > 
..param.TAlphabet:The alphabet type that is used for the transition labels.
...metafunction:Metafunction.Alphabet
...remarks:Use @Metafunction.Alphabet@ to get the type of the labels in an automaton.
...default:$char$
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use @Metafunction.Cargo@ to get the cargo type of an undirected graph.
...default:$void$
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
...remarks:Use WithoutEdgeId here to omit edge ids.
Note: If edges do not store ids external property maps do not work.
...default:$Default$, see @Tag.Default@.
..include:graph.h
*/

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TPropertyMap, typename TChar>
inline void
addLetterToOracle(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
				  TPropertyMap& supplyState,
				  TChar const c)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TVertexDescriptor newState = addVertex(g);
	TVertexDescriptor pred = newState - 1;
	addEdge(g, pred, newState, c);
	TVertexDescriptor k = getProperty(supplyState, pred);
	while ((k!=nilVal) &&
			(getTarget(&g.data_vertex[k].data_edge[(TSize) TAlphabet(c)])==nilVal))
	{
		addEdge(g,k,newState,c);
		k = getProperty(supplyState, k);
	}
	TVertexDescriptor s;
	if (k==nilVal) s=0;
	else s = getTarget(&g.data_vertex[k].data_edge[(TSize) TAlphabet(c)]);
	assignProperty(supplyState, newState, s);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.createOracle:
..cat:Spec.Oracle
..summary:Creates a factor oracle.
..signature:createOracle(g, text)
..param.g:Out-parameter: An automaton.
...type:Spec.Oracle
..param.text:In-parameter: A string.
...type:Class.String
..returns:void
..see:Function.createOracleOnReverse
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TText>
inline void
createOracle(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
			 TText const text)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TText>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TSize len = length(text);
	String<TVertexDescriptor> supplyState;
	resize(supplyState, len+1);
	TVertexDescriptor v1 = addVertex(g);
	assignProperty(supplyState, v1, nilVal);
	for(TSize i = 0; i<len; ++i) {
		addLetterToOracle(g, supplyState, getValue(text,i));
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.createOracleOnReverse:
..cat:Spec.Oracle
..summary:Creates a factor oracle for the reversed string.
..signature:createOracleOnReverse(g, text)
..param.g:Out-parameter: An automaton.
...type:Spec.Oracle
..param.text:In-parameter: A string.
...type:Class.String
..returns:void
..see:Function.createOracle
*/
template<typename TAlphabet, typename TCargo, typename TSpec, typename TText>
inline void
createOracleOnReverse(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
					  TText const text)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Size<TText>::Type TSize;
	TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
	TSize len = length(text);
	String<TVertexDescriptor> supplyState;
	resize(supplyState, len+1);
	TVertexDescriptor v1 = addVertex(g);
	assignProperty(supplyState, v1, nilVal);
	for(TSize i = len-1; i>0; --i) {
		addLetterToOracle(g, supplyState, getValue(text,i));
	}
	addLetterToOracle(g, supplyState, getValue(text,0));
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
