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

template <typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
void 
_markRootVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
	if (isRoot(g,v)) {
		append(str, ", shape = box");
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
void 
_markRootVertex(Graph<Directed<TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
void 
_markRootVertex(Graph<Undirected<TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
void 
_markRootVertex(Graph<Tree<TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
	if (isRoot(g,v)) {
		append(str, ", shape = box");
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TNodeAttributes>
void _createNodeAttributes(Graph<TSpec> const& g,
						   TNodeAttributes& nodeMap)
{
	SEQAN_CHECKPOINT
    typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	resizeVertexMap(g, nodeMap);
	char strV[BitsPerValue<TVertexDescriptor>::VALUE];

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		String<char> tmp("label = \"");
		sprintf(strV, "%d", *it);
		append(tmp, strV);
		append(tmp, "\"");
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}


//////////////////////////////////////////////////////////////////////////////
template<typename TSpec, typename TEdgeAttributes>
void _createEmptyEdgeAttributes(Graph<TSpec> const& g,
								TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		assignProperty(edgeMap, *itEd, String<char>(""));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
void _createEdgeAttributes(Graph<Directed<TCargo, TSpec> > const& g,
						   TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
void _createEdgeAttributes(Graph<Undirected<TCargo, TSpec> > const& g,
						   TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
void _createEdgeAttributes(Graph<Tree<TCargo, TSpec> > const& g,
						   TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeAttributes>
void _createEdgeAttributes(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
						   TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		String<char> tmp("label = \"");
		append(tmp, label(itEd));
		append(tmp, "\"");
		assignProperty(edgeMap, *itEd, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeAttributes>
void _createEdgeAttributes(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > const& g,
						   TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > TGraph;
	resizeEdgeMap(g, edgeMap);

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
		String<char> tmp("label = \"");
		append(tmp, str);
		append(tmp, "\"");
		assignProperty(edgeMap, *itEd, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
void _writeGraphType(TFile & file,
					 Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
					 DotDrawing)
{
	_streamWrite(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
void _writeGraphType(TFile & file,
					 Graph<Directed<TCargo, TSpec> > const& g,
					 DotDrawing)
{
	_streamWrite(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
void _writeGraphType(TFile & file,
					 Graph<Undirected<TCargo, TSpec> > const& g,
					 DotDrawing)
{
	_streamWrite(file, "graph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
void _writeGraphType(TFile & file,
					 Graph<Tree<TCargo, TSpec> > const& g,
					 DotDrawing)
{
	_streamWrite(file, "graph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
void _writeEdgeType(TFile & file,
					Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
					DotDrawing)
{
	_streamWrite(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
void _writeEdgeType(TFile & file,
					Graph<Directed<TCargo, TSpec> > const& g,
					DotDrawing)
{
	_streamWrite(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
void _writeEdgeType(TFile & file,
					Graph<Undirected<TCargo, TSpec> > const& g,
					DotDrawing)
{
	_streamWrite(file, " -- ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
void _writeEdgeType(TFile & file,
					Graph<Tree<TCargo, TSpec> > const& g,
					DotDrawing)
{
	_streamWrite(file, " -- ");
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TSpec, typename TNodeAttributes, typename TEdgeAttributes>
void write(TFile & file, 
		   Graph<TSpec> const& g,
		   TNodeAttributes const& nodeMap,
		   TEdgeAttributes const& edgeMap,
		   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	_writeGraphType(file,g,DotDrawing());
	_streamWrite(file, " G {\n");
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
		_streamPutInt(file, *it);
		_streamWrite(file, " [");
		_streamWrite(file, getProperty(nodeMap, *it));
		_streamWrite(file, "];\n");
	}
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		TVertexDescriptor sc = sourceVertex(itEd);
		TVertexDescriptor tr = targetVertex(itEd);
		_streamPutInt(file, sc);
		_writeEdgeType(file, g, DotDrawing());
		_streamPutInt(file, tr);
		_streamWrite(file, " [");
		_streamWrite(file, getProperty(edgeMap, *itEd));
		_streamWrite(file, "];\n");
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
	_createNodeAttributes(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
