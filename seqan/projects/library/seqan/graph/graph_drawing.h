#ifndef SEQAN_HEADER_GRAPH_DRAWING_H
#define SEQAN_HEADER_GRAPH_DRAWING_H

namespace SEQAN_NAMESPACE_MAIN
{

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




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
