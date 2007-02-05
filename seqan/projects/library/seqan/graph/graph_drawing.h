#ifndef SEQAN_HEADER_GRAPH_DRAWING_H
#define SEQAN_HEADER_GRAPH_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

struct TagDotDrawing_;
typedef Tag<TagDotDrawing_> const DotDrawing;

template <typename TFile, typename TEdges, typename TSpec, typename TNodeMap, typename TEdgeMap>
void write(TFile & file, 
	   Graph<TEdges, TSpec> const& g,
	   TNodeMap const& nodeMap,
	   TEdgeMap const& edgeMap,
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	typedef Graph<TEdges, TSpec> TGraph;
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
	typedef typename Iterator<TGraph, VertexIterator<> >::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		_streamWrite(file, getProperty(nodeMap, *it));
		_streamPut(file, ';');
		_streamPut(file, '\n');
	}
	_streamPut(file, '\n');

	_streamWrite(file, "/* Edges */\n");
	typedef typename Iterator<TGraph, EdgeIterator<> >::Type TConstEdIter;
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

template <typename TEdges, typename TSpec, typename TNodeMap>
void _createNodeNames(Graph<TEdges, TSpec> const& g,
		      TNodeMap& nodeMap)
{
        typedef Graph<TEdges, TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	initVertexMap(g, nodeMap);
	char strV[BitsPerValue<TVertexDescriptor>::VALUE];

	typedef typename Iterator<TGraph, VertexIterator<> >::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		sprintf(strV, "%d", *it);
		assignProperty(nodeMap, *it, String<char>(strV));
	}
}

template <typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeMap>
void _createEdgeNames(Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> const& g,
		      TEdgeMap& edgeMap)
{
	typedef Graph<EdgeList<TCargo, TEdgeSpec>, TSpec> TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	initEdgeMap(g, edgeMap);
	char strE[3+2*(BitsPerValue<TEdgeDescriptor>::VALUE)];

	typedef typename Iterator<TGraph, EdgeIterator<> >::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		sprintf(strE, "(%d,%d)", sourceVertex(itEd), targetVertex(itEd));
		assignProperty(edgeMap, *itEd, String<char>(strE));
	}
}

template <typename TAlphabet, typename TCargo, typename TEdgeSpec, typename TSpec, typename TEdgeMap>
void _createEdgeNames(Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> const& g,
		      TEdgeMap& edgeMap)
{
	typedef Graph<Automaton<TAlphabet, TCargo, TEdgeSpec>, TSpec> TGraph;
	initEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator<> >::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		assignProperty(edgeMap, *itEd, (*itEd).i2);
	}
}

template <typename TAlphabet, typename TCargo, typename TSpec, typename TGraphSpec, typename TEdgeMap>
void _createEdgeNames(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> const& g,
		      TEdgeMap& edgeMap)
{
	typedef Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> >, TGraphSpec> TGraph;
	initEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator<> >::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		String<TAlphabet> label = getProperty(g.data_edge_label, *itEd);
		String<char> str;
		resize(str,length(label)+1);
		value(str,0) = (*itEd).i2;
		typename Iterator<String<TAlphabet> >::Type it = begin(label);
		for(;!atEnd(it);++it) {
			char c = getValue(it);
			value(str,position(it) + 1) = c;
		}
		assignProperty(edgeMap, *itEd, str);
	}
}

template <typename TFile, typename TEdges, typename TSpec>
void write(TFile & file, 
	   Graph<TEdges, TSpec> const& g, 
	   DotDrawing) 
{
	SEQAN_CHECKPOINT
	String<String<char> > nodeMap;
	_createNodeNames(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeNames(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}



/*
struct TagPsTricks_;
typedef Tag<TagPsTricks_> const PsTricks;

template<typename TType> 
inline TType 
getUIntRand(TType inclusive_min, 
			TType inclusive_max) 
{
  return inclusive_min+(TType)((inclusive_max-inclusive_min+1)*((double)rand()/(RAND_MAX+1.0)));
}

inline double 
getUniform() {
  return ( (double) rand() / (RAND_MAX+1.0) );
}

template <typename TFile>
inline void 
_makePsTricksHeader(TFile & file) {
  _streamWrite(file, "\\documentclass{article}\n");
  _streamWrite(file, "\\usepackage{pst-all}\n");
  _streamWrite(file, "\\usepackage{pst-poly}\n");
  _streamWrite(file, "\\usepackage{multido}\n");
  _streamWrite(file, "\\usepackage{pstricks}\n");
  _streamWrite(file, "\\pagestyle{empty}\n");
  _streamWrite(file, "\\parindent=0pt\n");
  _streamWrite(file, "\\begin{document}\n");
  _streamWrite(file, "\\psset{xunit=0.1cm,yunit=0.1cm,runit=0.1cm}\n");
}

template <typename TFile>
inline void 
_makePsTricksFooter(TFile & file) {
  _streamWrite(file, "\\end{document}\n");
}

template <typename TFile, typename TVertexDescriptor, typename TPos>
inline void 
_makeNode(TFile & file,
		  TVertexDescriptor const v,
		  TPos x,
		  TPos y) 
{
  _streamWrite(file, "\\cnodeput(");
  _streamPutInt(file, x);
  _streamPut(file, ',');
  _streamPutInt(file, y);
  _streamWrite(file, "){");
  _streamPutInt(file, v);
  _streamWrite(file, "}{");
  _streamPutInt(file, v);
  _streamWrite(file, "}\n");
}

template <typename TFile, typename TVertexDescriptor>
inline void 
_makeEdge(TFile & file,
		  TVertexDescriptor const source,
		  TVertexDescriptor const sink) 
{
  _streamWrite(file, "\\ncarc{->}{");
  _streamPutInt(file, source);
  _streamWrite(file, "}{");
  _streamPutInt(file, sink);
  _streamWrite(file, "}\n");
}


template <typename TFile, typename TEdges, typename TSpec>
void _write_impl(TFile & file, 
				 Graph<TEdges, TSpec> const& g, 
				 PsTricks) 
{
	SEQAN_CHECKPOINT
	typedef unsigned int TCoord;
	typedef Size<TCoord>::Type TSize;
	srand((unsigned int)time(NULL));

	TCoord x = 0;
	TCoord y = 0;
	TSize width = 100;
	TSize height = 100;

	_makePsTricksHeader(file);
	_streamWrite(file, "\\begin{pspicture}(");
	_streamPutInt(file, x);
	_streamPut(file, ',');
	_streamPutInt(file, y);
	_streamWrite(file, ")(");
	_streamPutInt(file, width);
	_streamPut(file, ',');
	_streamPutInt(file, height);
	_streamWrite(file, ")\n");
	
	typedef typename VertexDescriptor<Graph<TEdges, TSpec> >::Type TVertexDescriptor;
	typedef typename Iterator<Graph<TEdges, TSpec>, EdgeIterator<> >::Type TEdgeIterator;
	TEdgeIterator it(g);
	std::set<TVertexDescriptor> vertexBag;
	for(;!atEnd(it);goNext(it)) {
		TVertexDescriptor source = sourceVertex(g, getValue(it));
		TVertexDescriptor sink = targetVertex(g, getValue(it));
		unsigned int posX;
		unsigned int posY;
		if (vertexBag.find(source) == vertexBag.end()) {
			posX = getUIntRand(x+1,width-1);
			posY = getUIntRand(y+1,height-1);
			_makeNode(file,source,posX,posY);
			vertexBag.insert(source);
		}
		if (vertexBag.find(sink) == vertexBag.end()) {
			posX = getUIntRand(x+1,width-1);
			posY = getUIntRand(y+1,height-1);
			_makeNode(file,sink,posX,posY);
			vertexBag.insert(sink);
		}
		_makeEdge(file,source,sink);
	}
	_streamWrite(file, "\\end{pspicture}\n");
	_makePsTricksFooter(file);
}
//____________________________________________________________________________

template <typename TFile, typename TEdges, typename TSpec>
void write(TFile & file, 
		   Graph<TEdges, TSpec> const& g, 
		   PsTricks) 
{
	SEQAN_CHECKPOINT
	_write_impl(file, g, PsTricks());
}

*/


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
