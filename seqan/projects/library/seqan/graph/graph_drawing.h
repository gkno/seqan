#ifndef SEQAN_HEADER_GRAPH_DRAWING_H
#define SEQAN_HEADER_GRAPH_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Drawing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TNodeMap>
void _createNodeNames(Graph<TSpec> const& g,
					  TNodeMap& nodeMap)
{
	SEQAN_CHECKPOINT
    typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	initVertexMap(g, nodeMap);
	char strV[BitsPerValue<TVertexDescriptor>::VALUE];

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		sprintf(strV, "\"%d\"", *it);
		assignProperty(nodeMap, *it, String<char>(strV));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeMap>
void _createEdgeNames(Graph<Directed<TCargo, TSpec> > const& g,
					  TEdgeMap& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	initEdgeMap(g, edgeMap);
	char strE[20];

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		sprintf(strE, "(%d,%d)", sourceVertex(itEd), targetVertex(itEd));
		assignProperty(edgeMap, *itEd, String<char>(strE));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeMap>
void _createEdgeNames(Graph<Tree<TCargo, TSpec> > const& g,
					  TEdgeMap& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	initEdgeMap(g, edgeMap);
	char strE[20];

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		sprintf(strE, "(%d,%d)", sourceVertex(itEd), targetVertex(itEd));
		assignProperty(edgeMap, *itEd, String<char>(strE));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeMap>
void _createEdgeNames(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
					  TEdgeMap& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	initEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		assignProperty(edgeMap, *itEd, label(itEd));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeMap>
void _createEdgeNames(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > const& g,
					  TEdgeMap& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > TGraph;
	initEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		String<TAlphabet> labelTmp = getCargo(*itEd);
		String<char> str;
		resize(str,length(labelTmp)+1);
		value(str,0) = label(itEd);
		typename Iterator<String<TAlphabet> >::Type it = begin(labelTmp);
		for(;!atEnd(it);++it) {
			char c = convert<char>(getValue(it));
			value(str,position(it) + 1) = c;
		}
		assignProperty(edgeMap, *itEd, str);
	}
}


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec, typename TNodeMap, typename TEdgeMap>
void write(TFile & file, 
	   Graph<TSpec> const& g,
	   TNodeMap const& nodeMap,
	   TEdgeMap const& edgeMap,
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	_streamWrite(file, "digraph G {\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Graph Attributes */\n");
	_streamWrite(file, "graph [rankdir = LR];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Node Attributes */\n");
	_streamWrite(file, "node [shape = ellipse, fillcolor = lightgrey, style = filled, fontname = \"Times-Italic\"];\n");
	_streamPut(file, '\n');
	_streamWrite(file, "/* Edge Attributes */\n");
	_streamWrite(file, "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n");
	_streamPut(file, '\n');

	_streamWrite(file, "/* Nodes */\n");
	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		_streamWrite(file, getProperty(nodeMap, *it));
		_streamPut(file, ';');
		_streamPut(file, '\n');
	}
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		TVertexDescriptor sc = sourceVertex(itEd);
		TVertexDescriptor tr = targetVertex(itEd);
		_streamWrite(file, getProperty(nodeMap, sc));
		_streamWrite(file, " -> ");
		_streamWrite(file, getProperty(nodeMap, tr));
		_streamWrite(file, " [label = \"");
		_streamWrite(file, getProperty(edgeMap, *itEd));
		_streamWrite(file, "\"];\n");
	}
	_streamPut(file, '\n');

	_streamWrite(file, "}\n");
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec>
void write(TFile & file, 
	   Graph<TSpec> const& g, 
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	String<String<char> > nodeMap;
	_createNodeNames(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeNames(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
