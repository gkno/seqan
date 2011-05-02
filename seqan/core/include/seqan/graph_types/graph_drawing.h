// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_GRAPH_DRAWING_H
#define SEQAN_HEADER_GRAPH_DRAWING_H

#include <seqan/align.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Drawing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// WRITING
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
	SEQAN_CHECKPOINT
	if (isRoot(g,v)) {
		append(str, ", shape = doublecircle");
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Directed<TCargo, TSpec> > const&,
				TVertexDescriptor const&,
				TAttributes&)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Undirected<TCargo, TSpec> > const&,
				TVertexDescriptor const&,
				TAttributes&)
{
	SEQAN_CHECKPOINT
}


//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TVertexDescriptor, typename TAttributes>
inline void 
_markRootVertex(Graph<Tree<TCargo, TSpec> > const& g,
				TVertexDescriptor const& v,
				TAttributes& str)
{
	SEQAN_CHECKPOINT
	if (isRoot(g,v)) {
		append(str, ", shape = doublecircle");
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TPosition, typename TNodeMap>
inline void
_createTrieNodeAttributes(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
						  String<String<TPosition> > pos,
						  TNodeMap& nodeMap)
{
	SEQAN_CHECKPOINT
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	resizeVertexMap(g, nodeMap);
	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		String<char> tmp;
		std::stringstream s;
		s << *it;
		String<TPosition> endPositions = getProperty(pos,*it);
		if (!empty(endPositions)) {
			s <<  " {";
			append(tmp, "shape = box, ");
			typename Iterator<String<TPosition>, Rooted>::Type itP = begin(endPositions);
			typename Iterator<String<TPosition>, Rooted>::Type beginP = itP;
			for(;!atEnd(itP);goNext(itP)) {
				if (beginP != itP) s << ", ";
				s << *itP;
			}
			s << "}";
		}
		
		append(tmp, "label = \"");
		append(tmp, s.str().c_str());
		append(tmp, "\"");
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TNodeAttributes>
inline void
_createNodeAttributes(Graph<TSpec> const& g,
					  TNodeAttributes& nodeMap)
{
	SEQAN_CHECKPOINT
    typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	resizeVertexMap(g, nodeMap);

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		std::ostringstream outs; 
		outs << "label = \"";
		outs << *it;
		outs << "\"";
		String<char> tmp;
		append(tmp, outs.str().c_str());
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TNodeAttributes, typename TNameMap>
inline void
_createNodeAttributes(Graph<TSpec> const& g,
					  TNodeAttributes& nodeMap,
					  TNameMap const& nameMap)
{
    typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	resizeVertexMap(g, nodeMap);

	typedef typename Iterator<TGraph, VertexIterator>::Type TConstIter;
	TConstIter it(g);
	for(;!atEnd(it);++it) {
		std::ostringstream outs; 
		outs << "label = \"";
        outs << getProperty(nameMap,*it);
        outs << "\"";
		String<char> tmp;
		append(tmp, outs.str().c_str());
		_markRootVertex(g, *it, tmp);
		assignProperty(nodeMap, *it, tmp);
	}
}

//////////////////////////////////////////////////////////////////////////////
template<typename TSpec, typename TEdgeAttributes>
inline void
_createEmptyEdgeAttributes(Graph<TSpec> const& g,
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
inline void
_createEdgeAttributes(Graph<Directed<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Undirected<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Tree<void, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	SEQAN_CHECKPOINT
	_createEmptyEdgeAttributes(g,edgeMap);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void 
_createEdgeAttributes(Graph<Tree<TCargo, TSpec> > const& g,
					  TEdgeAttributes& edgeMap)
{
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	resizeEdgeMap(g, edgeMap);

	typedef typename Iterator<TGraph, EdgeIterator>::Type TConstEdIter;
	TConstEdIter itEd(g);
	for(;!atEnd(itEd);++itEd) {
		std::ostringstream outs; 
		outs << "label = \"";
		outs << (TCargo) getCargo(*itEd);
		outs << "\"";
		append(property(edgeMap, *itEd), outs.str().c_str());		
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TCargo, typename TSpec, typename TEdgeAttributes>
inline void
_createEdgeAttributes(Graph<Automaton<TAlphabet, TCargo, TSpec> > const& g,
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
inline void
_createEdgeAttributes(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > > const& g,
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
		typename Iterator<String<TAlphabet>, Rooted>::Type it = begin(labelTmp);
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

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Directed<TCargo, TSpec> > const&,
				  DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Undirected<TCargo, TSpec> > const&,
				  DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Tree<TCargo, TSpec> > const&,
				  DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeGraphFooter(TFile &,
				  Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
				  DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
				DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
	_streamWrite(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Directed<TCargo, TSpec> > const&,
				DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
	_streamWrite(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Undirected<TCargo, TSpec> > const&,
				DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
	_streamWrite(file, "graph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeGraphType(TFile & file,
				Graph<Tree<TCargo, TSpec> > const&,
				DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
	_streamWrite(file, "digraph");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TAlphabet, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Automaton<TAlphabet, TCargo, TSpec> > const&,
			   DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
	_streamWrite(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Directed<TCargo, TSpec> > const&,
			   DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
	_streamWrite(file, " -> ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Undirected<TCargo, TSpec> > const&,
			   DotDrawing)
{
//IOREV
	SEQAN_CHECKPOINT
	_streamWrite(file, " -- ");
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TCargo, typename TSpec>
inline void
_writeEdgeType(TFile & file,
			   Graph<Tree<TCargo, TSpec> > const&,
			   DotDrawing)
{
   //IOREV
	SEQAN_CHECKPOINT
	_streamWrite(file, " -> ");
}



//////////////////////////////////////////////////////////////////////////////
// internal functions for cgviz output
//////////////////////////////////////////////////////////////////////////////

    
//////////////////////////////////////////////////////////////////////////////  
// compute node numbers in cgviz file
    
unsigned _seqData(unsigned i) { return 3*i; }
    
unsigned _pairData(unsigned i, unsigned j, unsigned k) { 
    // initialize with the next free node after segement nodes
    unsigned ret= 3*k; 
        
    for(unsigned l=0; l<i; ++l)
        for(unsigned m=l+1; m<k; ++m)
            ret += 3;
        
    return (ret+(j-i-1)*3);
        
}
    
template <typename TFile>
inline void
_writeSegmentGlyph(TFile & file,
                   unsigned nodenum,
                   CGViz)
{
    //IOREV
    SEQAN_CHECKPOINT
    _streamWrite(file, "{GLYPH SegmentGlyph_node_");
    _streamPutInt(file,nodenum); 
    _streamWrite(file, "\n");
    _streamWrite(file, "  color=0,255,0,128\n");   
    _streamWrite(file, "  drawerName=Blocks\n");   
    _streamWrite(file, "  drawerProperties=[centered=false closed=true filled=false upsideDown=false vertical=false width=1]\n");   
    _streamWrite(file, "  highlightColor=255,0,0,255\n");   
    _streamWrite(file, "  highlightSelected=true\n");   
    _streamWrite(file, "  lineWidth=4\n");   
    _streamWrite(file, "  swapAxes=false\n");   
    _streamWrite(file, "  visibility=true\n"); 
    _streamWrite(file, "}\n");
}
    
//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void
_writeSegmentPane(TFile & file,
                   unsigned nodenum,
                   unsigned num,
                   unsigned length,
                   CGViz)
{
    //IOREV
    SEQAN_CHECKPOINT
    _streamWrite(file, "{PANE SegmentPane_node_");
    _streamPutInt(file,nodenum); 
    _streamWrite(file, "\n");
    _streamWrite(file, "  antiAliasing=false\n");
    _streamWrite(file, "  axesLogarithmic=false,false\n");
    _streamWrite(file, "  axisLabels=Seq");
    _streamPutInt(file,num); 
    _streamWrite(file, ",\n");
    _streamWrite(file, "  axisLabelsOpposite=,\n");
    _streamWrite(file, "  axisStarts=0.0,0.0\n");
    _streamWrite(file, "  axisStops=");
    _streamPutFloat(file,length); 
    _streamWrite(file, ",");
    _streamPutFloat(file,length); 
    _streamWrite(file, "\n");
    _streamWrite(file, "  axisTicks=true,false\n");
    _streamWrite(file, "  axisTicksOpposite=false,false\n");
    _streamWrite(file, "  color=255,255,255,255\n");
    _streamWrite(file, "  coordinatesParallel=false\n");
    _streamWrite(file, "  lockHeight=false\n");
    _streamWrite(file, "  lockWidth=false\n");
    _streamWrite(file, "  mayZoomScroll=true\n");
    _streamWrite(file, "  pixelHeight=30\n");
    _streamWrite(file, "  pixelWidth=800\n");
    _streamWrite(file, "  shapeName=DefaultPaneShape\n");
    _streamWrite(file, "  shapeProperties=[]\n");
    _streamWrite(file, "  swapXY=false\n");
    _streamWrite(file, "  uvLocked=false\n");
    _streamWrite(file, "  visibility=true\n");
    _streamWrite(file, "}\n");
}


    
//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void
_writeSegmentMatchGlyph(TFile & file, unsigned nodenum,  unsigned r, unsigned b, unsigned g, CGViz)
{
        //IOREV
        SEQAN_CHECKPOINT
        _streamWrite(file, "{GLYPH SegmentMatchGlyph_node_");
        _streamPutInt(file,nodenum); 
        _streamWrite(file, "\n");
        _streamWrite(file, "  color=");
        _streamPutInt(file,r); 
        _streamWrite(file, ",");
        _streamPutInt(file,b); 
        _streamWrite(file, ",");
        _streamPutInt(file,g); 
        _streamWrite(file, ",153\n");
        _streamWrite(file, "  drawerName=Matches\n");   
        _streamWrite(file, "  drawerProperties=[gap=5 wrapAround=true]\n");   
        _streamWrite(file, "  highlightColor=255,0,0,204\n");   
        _streamWrite(file, "  highlightSelected=true\n");   
        _streamWrite(file, "  lineWidth=2\n");   
        _streamWrite(file, "  swapAxes=false\n");   
        _streamWrite(file, "  visibility=true\n"); 
        _streamWrite(file, "}\n");
        
}
    
//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void
_writeSegmentMatchPane(TFile & file, unsigned nodenum, unsigned i, unsigned j, CGViz)
    {
        //IOREV
        SEQAN_CHECKPOINT
        _streamWrite(file, "{PANE SegmentMatchPane_node_");
        _streamPutInt(file,nodenum); 
        _streamWrite(file, "\n");
        _streamWrite(file, "  antiAliasing=false\n");
        _streamWrite(file, "  axesLogarithmic=false,false\n");
        _streamWrite(file, "  axisLabels=,\n");
        _streamWrite(file, "  axisLabelsOpposite=,\n");
        _streamWrite(file, "  axisStarts=1.0,1.0\n");
        _streamWrite(file, "  axisStops=700.0,700.0\n");    
        _streamWrite(file, "  axisTicks=false,false\n");
        _streamWrite(file, "  axisTicksOpposite=false,false\n");
        _streamWrite(file, "  color=255,255,255,0\n");
        _streamWrite(file, "  coordinatesParallel=true\n");
        _streamWrite(file, "  lockHeight=false\n");
        _streamWrite(file, "  lockWidth=false\n");
        _streamWrite(file, "  mayZoomScroll=false\n");
        // we hve to accommodate different heights    
        _streamWrite(file, "  pixelHeight=");
        _streamPutInt(file,(j-i-1)*200+170); 
        _streamWrite(file, "\n");
        _streamWrite(file, "  pixelWidth=800\n");
        _streamWrite(file, "  shapeName=DefaultPaneShape\n");
        _streamWrite(file, "  shapeProperties=[]\n");
        _streamWrite(file, "  swapXY=false\n");
        _streamWrite(file, "  uvLocked=false\n");
        _streamWrite(file, "  visibility=true\n");
        _streamWrite(file, "}\n");
  
}
    
//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void
_writeWindow(TFile & file, unsigned k, CGViz)
    {
        //IOREV
        SEQAN_CHECKPOINT
        _streamWrite(file, "{WINDOW AlignmentGraphWindow");
        _streamWrite(file, "\n");
        _streamWrite(file, "  color=255,255,255,255\n");
        _streamWrite(file, "  contentHeight=");
        _streamPutInt(file,200+k*200); 
        _streamWrite(file, "\n");
        _streamWrite(file, "  contentWidth=900\n");
        _streamWrite(file, "  sourceFile=dummy\n");
        _streamWrite(file, "  title=AlignmentGraph\n");
        _streamWrite(file, "  zoomPanesWithWindow=true\n");
        _streamWrite(file, "}\n");
}
    
//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void
_writeFeeder(TFile & file, int i, CGViz)
{

    //IOREV
    SEQAN_CHECKPOINT
    _streamWrite(file, "{FEEDER Feeder_for_seq_");
    _streamPutInt(file,i);
    _streamWrite(file, " ");
    _streamPutInt(file,_seqData(i));
    _streamWrite(file, " ");
    _streamPutInt(file,_seqData(i)+1);
    _streamWrite(file, "\n");
    _streamWrite(file, "  filter=\n");
    _streamWrite(file, "}\n");
}

//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void
_writeThreader(TFile & file, int i, CGViz)
{
        //IOREV
        SEQAN_CHECKPOINT
        _streamWrite(file, "{THREADER THREADER_for_seq_");
        _streamPutInt(file,i);
        _streamWrite(file, " ");
        _streamPutInt(file,_seqData(i)+1);
        _streamWrite(file, " ");
        _streamPutInt(file,_seqData(i)+2);
        _streamWrite(file, "\n");
        _streamWrite(file, "  depth=10\n");
        _streamWrite(file, "  filter=\n");
        _streamWrite(file, "  trackNumber=0\n");
        _streamWrite(file, "  trackSpan=1\n");
        _streamWrite(file, "  trackWeight=1.0\n");
        _streamWrite(file, "}\n");
        
}
    
    
//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void _writeAnchor(TFile & file, int i, int k, CGViz)
{
        //IOREV
        SEQAN_CHECKPOINT
        _streamWrite(file, "{ANCHOR ANCHOR_for_seq_");
        _streamPutInt(file,i);
        _streamWrite(file, " ");
        _streamPutInt(file,_seqData(i)+2);
        _streamWrite(file, " ");
        _streamPutInt(file,_pairData(k-2,k-1,k)+3);   
        _streamWrite(file, "\n");
        _streamWrite(file, "  depth=9\n");
        _streamWrite(file, "  xPos=30\n");
        _streamWrite(file, "  yPos=");
        _streamPutInt(file,(i+1)*200);
        _streamWrite(file, "\n");
        _streamWrite(file, "}\n");
        
}


//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void
_writeFeeder(TFile & file, int i, int j, int k, CGViz)
{
        
    //IOREV
    SEQAN_CHECKPOINT
    _streamWrite(file, "{FEEDER Feeder_for_seq_pair");
    _streamPutInt(file,i);
    _streamWrite(file, "_");
    _streamPutInt(file,j);
    _streamWrite(file, " ");
    _streamPutInt(file,_pairData(i,j,k));
    _streamWrite(file, " ");
    _streamPutInt(file,_pairData(i,j,k)+1);
    _streamWrite(file, "\n");
    _streamWrite(file, "  filter=\n");
    _streamWrite(file, "}\n");
}
    
//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void
_writeThreader(TFile & file, int i, int j, int k, CGViz)
{
    //IOREV
    SEQAN_CHECKPOINT
    _streamWrite(file, "{THREADER THREADER_for_seq_pair");
    _streamPutInt(file,i);
    _streamWrite(file, "_");
    _streamPutInt(file,j);
    _streamWrite(file, " ");
    _streamPutInt(file,_pairData(i,j,k)+1);
    _streamWrite(file, " ");
    _streamPutInt(file,_pairData(i,j,k)+2);
    _streamWrite(file, "\n");
    _streamWrite(file, "  depth=10\n");
    _streamWrite(file, "  filter=\n");
    _streamWrite(file, "  trackNumber=0\n");
    _streamWrite(file, "  trackSpan=1\n");
    _streamWrite(file, "  trackWeight=1.0\n");
    _streamWrite(file, "}\n");
        
}
    
    
//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
    inline void _writeAnchor(TFile & file, int i, int j, int k, CGViz)
    {
    //IOREV
    SEQAN_CHECKPOINT
    _streamWrite(file, "{ANCHOR ANCHOR_for_seq_pair");
    _streamPutInt(file,i);
    _streamWrite(file, "_");
    _streamPutInt(file,j);
    _streamWrite(file, " ");
    _streamPutInt(file,_pairData(i,j,k)+2);
    _streamWrite(file, " ");
    _streamPutInt(file,_pairData(k-2,k-1,k)+3);   
    _streamWrite(file, "\n");
    _streamWrite(file, "  depth=9\n");
    _streamWrite(file, "  xPos=30\n");
    _streamWrite(file, "  yPos=");
    _streamPutInt(file,(i+1)*200+30);
    _streamWrite(file, "\n");
    _streamWrite(file, "}\n");
        
}
    

//////////////////////////////////////////////////////////////////////////////  
template <typename TFile>
inline void _writePaneLocks(TFile & file, int i, int j, int k, CGViz)
{
    //IOREV
    SEQAN_CHECKPOINT
    _streamWrite(file, "{PANELOCK PANELOCK_1_for_seq_pair");
    _streamPutInt(file,i);
    _streamWrite(file, "_");
    _streamPutInt(file,j);
    _streamWrite(file, " ");
    _streamPutInt(file,_seqData(i)+2);
    _streamWrite(file, " ");
    _streamPutInt(file,_pairData(i,j,k)+2);   
    _streamWrite(file, "\n");
    _streamWrite(file, "  bidirectional=false\n");
    _streamWrite(file, "  lockCenterU=true\n");
    _streamWrite(file, "  lockCenterV=false\n");
    _streamWrite(file, "  lockScaleU=true\n");
    _streamWrite(file, "  lockScaleV=false\n");
    _streamWrite(file, "  swapUV=false\n");
    _streamWrite(file, "}\n");
        
    _streamWrite(file, "{PANELOCK PANELOCK_2_for_seq_pair");
    _streamPutInt(file,i);
    _streamWrite(file, "_");
    _streamPutInt(file,j);
    _streamWrite(file, " ");
    _streamPutInt(file,_seqData(j)+2);
    _streamWrite(file, " ");
    _streamPutInt(file,_pairData(i,j,k)+2);   
    _streamWrite(file, "\n");
    _streamWrite(file, "  bidirectional=false\n");
    _streamWrite(file, "  lockCenterU=false\n");
    _streamWrite(file, "  lockCenterV=true\n");
    _streamWrite(file, "  lockScaleU=false\n");
    _streamWrite(file, "  lockScaleV=true\n");
    _streamWrite(file, "  swapUV=true\n");
    _streamWrite(file, "}\n");
}
    
    

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.write:
..signature:write(file, graph, nodeMap, edgeMap, tag)
..param.graph:The graph to write out.
...type:Class.Graph
..param.nodeMap:A mapping from vertex descriptor to vertex label.
..param.edgeMap:A mapping from edge descriptor to edge label.
..param.tag:A tag to select the output format.
...type:Tag.DotDrawing
..include:seqan/graph_types.h
 */
template <typename TFile, typename TSpec, typename TNodeAttributes, typename TEdgeAttributes>
void 
write(TFile & file, 
	  Graph<TSpec> const& g,
	  TNodeAttributes const& nodeMap,
	  TEdgeAttributes const& edgeMap,
	  DotDrawing) 
{
//IOREV _doc_ _batchreading_
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
	_streamWrite(file, "node [shape = circle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n");
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

	_writeGraphFooter(file,g,DotDrawing());

	_streamWrite(file, "}\n");
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.write:
..signature:write(file, graph, nodeMap, tag)
..include:seqan/graph_types.h
 */
template <typename TFile, typename TSpec, typename TNodeAttributes>
inline void
write(TFile & file,
	  Graph<TSpec> const& g, 
	  TNodeAttributes const& nodeMap,
	  DotDrawing) 
{
//IOREV _doc_ _batchreading_
	SEQAN_CHECKPOINT
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.write:
..signature:write(file, graph, tag)
..include:seqan/graph_types.h
 */
template <typename TFile, typename TSpec>
inline void
write(TFile & file,
	  Graph<TSpec> const& g, 
	  DotDrawing) 
{
//IOREV _doc_ _batchreading_
	SEQAN_CHECKPOINT
	String<String<char> > nodeMap;
	_createNodeAttributes(g,nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(g,edgeMap);
	write(file,g,nodeMap,edgeMap,DotDrawing());
}

    
    
    
    
    
/**
.Function.write:
..signature:write(file, graph, nodeMap, edgeMap, tag)
..param.graph:The graph to write out.
...type:Class.Graph
..param.nodeMap:A mapping from vertex descriptor to vertex label.
..param.edgeMap:A mapping from edge descriptor to edge label.
..param.tag:A tag to select the output format.
...type:Tag.DotDrawing
..include:seqan/graph_types.h
*/
template <typename TFile,  typename TAlignmentGraph, typename TNodeAttributes, typename TEdgeAttributes>
void 
write(TFile & file, 
      TAlignmentGraph & g,
      TNodeAttributes const& nodeMap,
      TEdgeAttributes const& edgeMap,
      CGViz) 
{
    //IOREV _doc_ _batchreading_
    SEQAN_CHECKPOINT
 
    /* 
    First write the DATA nodes:
     * for all sequences all nodes are consecutively read and lines in a CGviz DATA node are filled.   
     * the DATA nodes are output during the mapping. Make sure that it in right to left order and in order of the sequences. There is a
     * DATA node for each sequence. 
     
     * secondly the edges are traversed for each pair of sequences. A DATA node is generated for each pair of sequences. Write accessor functions
     that return the node number for sequences and pairs of sequences like:
     seq_data(int k)               
     pair_data(int i, int j)     
     the numbers are always increased by 3 since each data node has a glyph and panel
        
     Finally, a window node is inserted, then tehe edges  and Panelock nodes between sequence panels and pair panels
     */
    

    typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TAlignmentGraph>::Type TEdgeDescriptor;
    typedef typename Cargo<TAlignmentGraph>::Type TCargo;
    typedef typename Iterator<TAlignmentGraph,OutEdgeIterator>::Type TOutEdgeIterator;
     
    typedef typename StringSetType<TAlignmentGraph>::Type TStringSet;
    typedef typename Id<TAlignmentGraph>::Type TId;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename Position<typename Value<TStringSet>::Type>::Type TPosition;
     
    typedef Pair<TVertexDescriptor, TVertexDescriptor> TVertexPair;
    typedef Pair<TId, TPosition> TVertexInfo; // seqId, beginPos
    typedef Pair<TVertexInfo> TEdgeInfo; // source, target
     
    
    // iterate over all sequences and in the sequences form left to right over the nodes
     for (TSize i = 0; i < length(stringSet(g)); ++i) {
         // sequence ID
         TId id_i = positionToId(stringSet(g), i);
         //        std::cout << "id = " << id_i << std::endl;
         
         TPosition pos = beginPosition(value(stringSet(g), i));
         TVertexDescriptor vertex = findVertex(g, id_i, pos);
         
         // first interval
         // **** OUTPUT
         _streamWrite(file, "{DATA segments_for_seq");
         _streamPutInt(file, id_i);
         _streamWrite(file, "_node_");
         _streamPutInt(file, _seqData(i));
         _streamWrite(file, "\n  [__GLOBAL__]\n");

  
         std::cout << "string " << id_i << " interval " << pos << " " <<  pos + fragmentLength(g, vertex) << std::endl;
         _streamWrite(file, "[segment annotation] : ");
         _streamPutInt(file, pos);
         _streamWrite(file, " ");
         _streamPutInt(file, pos+fragmentLength(g, vertex));
         _streamWrite(file, "\n");  

         pos += fragmentLength(g, vertex);
         while (pos < endPosition(value(stringSet(g), i))) {
             vertex = findVertex(g, id_i, pos);
             // **** OUTPUT
             std::cout << "string " << id_i << " interval " << pos << " " <<  pos + fragmentLength(g, vertex) << std::endl;
             _streamWrite(file, "[segment annotation] : ");
             _streamPutInt(file, pos);
             _streamWrite(file, " ");
             _streamPutInt(file, pos+fragmentLength(g, vertex));
             _streamWrite(file, "\n");  

             pos += fragmentLength(g, vertex);
        }
         std::cout << "End Data node " << std::endl;
         _streamWrite(file, "}\n"); 
         // write glyph and panel for segments
         _writeSegmentGlyph(file,_seqData(i)+1,CGViz());
         _writeSegmentPane(file,_seqData(i)+2,i,pos,CGViz());   
     }
    
    
    // traverse edges
    for (TSize i = 0; i < length(stringSet(g))-1; ++i) {
        for(TSize j=i+1; j < length(stringSet(g)); j++) {
            std::cout << "Traversing sequence pair " << i << "," << j << std::endl;

            TId id_i = positionToId(stringSet(g), i);
            TPosition pos_i = beginPosition(value(stringSet(g), i));
            TVertexDescriptor vertex = findVertex(g, id_i, pos_i);

            // we look at the pair i,j  j>i
            // sequence ID
            TId id_j = positionToId(stringSet(g), j);
            
            _streamWrite(file, "{DATA segment_matches_for_pair_");
            _streamPutInt(file, id_i);
            _streamWrite(file, "_");
            _streamPutInt(file, id_j);
            _streamWrite(file, "_node_");
            _streamPutInt(file, _pairData(i,j,length(stringSet(g))));
            _streamWrite(file, "\n  [__GLOBAL__]\n");
            
            //       std::cout << "id_i, id_j = " << id_i << "," << id_j << std::endl;
            //      std::cout << "pos, id, degree " << pos_i << "," << vertex << "," << outDegree(g,vertex) << std::endl;
            TOutEdgeIterator outEdgeIt(g, vertex);

            while( !atEnd(outEdgeIt) ){
                TVertexDescriptor target = targetVertex(outEdgeIt);
                if( sequenceId(g,target) == id_j ){
                    TPosition targetPos = fragmentBegin(g, target);
                    TPosition targetEnd = fragmentBegin(g, target)+fragmentLength(g, target);
                    /// ****** OUTPUT 
                    std::cout << "segment match " << pos_i << " " << targetPos << " " << pos_i + fragmentLength(g, vertex) << " " << targetEnd << std::endl;

                    _streamWrite(file, "[segment match annotation] : ");
                    _streamPutInt(file, pos_i);
                    _streamWrite(file, " ");
                    _streamPutInt(file, targetPos);
                    _streamWrite(file, " ");
                    _streamPutInt(file, pos_i + fragmentLength(g, vertex));
                    _streamWrite(file, " ");
                    _streamPutInt(file, targetEnd);
                    _streamWrite(file, "\n");  
                }
                goNext(outEdgeIt);
                   
            }

            pos_i += fragmentLength(g, vertex);

            while (pos_i < endPosition(value(stringSet(g), i))) {
                vertex = findVertex(g, id_i, pos_i);
                // we look at the pair i,j  j>i
                // sequence ID
                TId id_j = positionToId(stringSet(g), j);
                
                //         std::cout << "id_i, id_j = " << id_i << "," << id_j << std::endl;
                //         std::cout << "pos, id, degree " << pos_i << "," << vertex << "," << outDegree(g,vertex) << std::endl;
                TOutEdgeIterator outEdgeIt(g, vertex);
                while( !atEnd(outEdgeIt) ){
                    TVertexDescriptor target = targetVertex(outEdgeIt);
                    if( sequenceId(g,target) == id_j ){
                        TPosition targetPos = fragmentBegin(g, target);
                        TPosition targetEnd = fragmentBegin(g, target)+fragmentLength(g, target);
                        // **** OUTPUT
                        std::cout << "segment match " << pos_i << " " << targetPos << " " << pos_i + fragmentLength(g, vertex) << " " << targetEnd << std::endl;
                        _streamWrite(file, "[segment match annotation] : ");
                        _streamPutInt(file, pos_i);
                        _streamWrite(file, " ");
                        _streamPutInt(file, targetPos);
                        _streamWrite(file, " ");
                        _streamPutInt(file, pos_i + fragmentLength(g, vertex));
                        _streamWrite(file, " ");
                        _streamPutInt(file, targetEnd);
                        _streamWrite(file, "\n");  
                                                
                    }   
                    goNext(outEdgeIt);
                }
                pos_i += fragmentLength(g, vertex);
            }
            std::cout << "End Data node " << std::endl;
            _streamWrite(file, "}\n"); 
            // output glyph and panel for matches
            _writeSegmentMatchGlyph(file,_pairData(i,j,length(stringSet(g)))+1,(i*80)%255,(j*40)%255,((i+j)*60)%255,CGViz());
            _writeSegmentMatchPane(file,_pairData(i,j,length(stringSet(g)))+2,i,j,CGViz());               
            
        }   

    }
    
    // write the window 
    _writeWindow(file,length(stringSet(g)),CGViz());
    
    // finally write the edges
    // first the feeders, threaders and anchors for the data nodes
    for (TSize i = 0; i < length(stringSet(g)); ++i) {
        _writeFeeder(file,i, CGViz());
        _writeThreader(file,i, CGViz());      
        _writeAnchor(file,i,length(stringSet(g)), CGViz());
    }
    // then the feeders, threaders and anchors for the match nodes
    for (TSize i = 0; i < length(stringSet(g))-1; ++i) {
        for(TSize j=i+1; j < length(stringSet(g)); j++) {
            _writeFeeder(file,i,j, length(stringSet(g)),CGViz());
            _writeThreader(file,i,j,length(stringSet(g)), CGViz());      
            _writeAnchor(file,i,j,length(stringSet(g)), CGViz());
        }  
    }
    // finally the Panelocks
    for (TSize i = 0; i < length(stringSet(g))-1; ++i) 
        for(TSize j=i+1; j < length(stringSet(g)); j++)
            _writePaneLocks(file,i,j, length(stringSet(g)),CGViz());
    
}


//////////////////////////////////////////////////////////////////////////////
    
/**
.Function.write:
 ..signature:write(file, graph, tag)
..include:seqan/graph_types.h
*/
template <typename TFile, typename TAlignmentGraph>
inline void
write(TFile & file,
      TAlignmentGraph & g, 
      CGViz) 
{
    //IOREV _doc_ _batchreading_
    SEQAN_CHECKPOINT
    String<String<char> > nodeMap;
    _createNodeAttributes(g,nodeMap);
    String<String<char> > edgeMap;
    _createEdgeAttributes(g,edgeMap);
    write(file,g,nodeMap,edgeMap,CGViz());
}


//////////////////////////////////////////////////////////////////////////////
// READING
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_addNode(Graph<TSpec>& g,
		 TStatement& node_id,
		 TStatement& attr_list,
		 TNodeAttributes& nodeMap,
		 TEdgeAttributes&,			  
		 TNodeIdMap& nodeIdMap)
{
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	if (nodeIdMap.find(node_id) == nodeIdMap.end()) {
		TVertexDescriptor _id = addVertex(g);
		nodeIdMap.insert(std::make_pair(node_id, _id));
		resizeVertexMap(g, nodeMap);
		assignProperty(nodeMap, _id, attr_list);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Directed<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Directed<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Undirected<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Undirected<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Tree<TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	TEdgeDescriptor e = addEdge(g, sourceV, targetV);
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TString>
inline typename Alphabet<Graph<Automaton<TAlphabet, TCargo, TSpec> > >::Type
_getInternalLabel(Graph<Automaton<TAlphabet, TCargo, TSpec> >&,
				  TString& str)
{
	return str[0];
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TString>
inline String<TAlphabet>
_getInternalLabel(Graph<Automaton<TAlphabet, TCargo, WordGraph<TSpec> > >&,
				  TString& str)
{
	return str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TVertexDescriptor, typename TNodeAttributes, typename TEdgeAttributes, typename TStatement>
inline void
_addEdge(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
		 TVertexDescriptor sourceV,
		 TVertexDescriptor targetV,
		 TNodeAttributes&,
		 TEdgeAttributes& edgeMap,
		 TStatement& attr_list)
{
	typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	// We need the label
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	typedef typename Position<TIter>::Type TPos;
	
	String<TValue> label;
	TIter it = begin(attr_list);
	bool found = false;
	for(;!atEnd(it);goNext(it)) {
		TPos pos = position(it);
		if (*it == ',') {
			found = false;
		} else if (found) {
			append(label, *it);
		} else if ((pos + 5 < length(attr_list)) &&
			(infix(attr_list, it, it + 5) == "label")) 
		{
				found = true;
				it += 5;
		}
	}
	TEdgeDescriptor e = addEdge(g, sourceV, targetV, _getInternalLabel(g, label));
	resizeEdgeMap(g, edgeMap);
	assignProperty(edgeMap, e, attr_list);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_addEdge(Graph<TSpec>& g,
		 TStatement& left_node_id,
		 TStatement& right_node_id,
		 TStatement& attr_list,
		 TNodeAttributes& nodeMap,
		 TEdgeAttributes& edgeMap,
		 TNodeIdMap& nodeIdMap)
{
	typedef Graph<TSpec> TGraph;
	typedef typename Value<TStatement>::Type TValue;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef std::map<String<TValue>, TVertexDescriptor> TMap;

	TVertexDescriptor sourceV;
	TVertexDescriptor targetV;

	typename TMap::iterator pos;
	pos = nodeIdMap.find(left_node_id);
	if (pos == nodeIdMap.end()) return;
	else sourceV = pos->second;

	pos = nodeIdMap.find(right_node_id);
	if (pos == nodeIdMap.end()) return;
	else targetV = pos->second;

	_addEdge(g, sourceV, targetV, nodeMap, edgeMap, attr_list);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_processNodeStatement(Graph<TSpec>& g,
					  TStatement& stmt,
					  TNodeAttributes& nodeMap,
					  TEdgeAttributes& edgeMap,
					  TNodeIdMap& nodeIdMap) 
{
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	
	String<TValue> node_id;
	String<TValue> attr_list;  // Multiple attribute lists are ignored
	bool inAttr = false;
	TIter it = begin(stmt);
	for(;!atEnd(it);goNext(it)) {
		if (*it == '[') {
			inAttr = true;
			continue;
		} else if (*it == ']') {
			// Finished
			break;
		} else if ((*it == ' ') ||
			(*it == '"')) {
			continue;
		}
		if (inAttr) {
			append(attr_list, *it);
		} else {
			append(node_id, *it);
		}
	}
	_addNode(g, node_id, attr_list, nodeMap, edgeMap, nodeIdMap);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TPosition, typename TNodeIdMap>
inline void
_processEdgeStatement(Graph<TSpec>& g,
					  TStatement& stmt,
					  TNodeAttributes& nodeMap,
					  TEdgeAttributes& edgeMap,
					  TPosition pos,
					  TNodeIdMap& nodeIdMap) 
{
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;
	
	String<TValue> left_node_id;
	String<TValue> right_node_id;
	String<TValue> attr_list;  // Multiple attribute lists are ignored
	bool inAttr = false;
	TIter it = begin(stmt);
	unsigned int localPos = 0;
	for(;!atEnd(it);goNext(it), ++localPos) {
		if (*it == '[') {
			inAttr = true;
			continue;
		} else if (*it == ']') {
			// Finished
			break;
		} else if ((*it == ' ') ||
			(*it == '"')) {
			continue;
		}
		if (inAttr) {
			append(attr_list, *it);
		} else if (localPos < pos) {
			append(left_node_id, *it);
		} else if (localPos > pos+1) {
			append(right_node_id, *it);
		}
	}
	//std::cout << left_node_id << "," << right_node_id << "," << std::endl;
	_addEdge(g, left_node_id, right_node_id, attr_list, nodeMap, edgeMap, nodeIdMap);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TStatement, typename TNodeAttributes, typename TEdgeAttributes, typename TNodeIdMap>
inline void
_processStatement(Graph<TSpec>& g,
				  TStatement& stmt,
				  TNodeAttributes& nodeMap,
				  TEdgeAttributes& edgeMap,
				  TNodeIdMap& nodeIdMap) 
{
	// Clear everything up to the last line
	typedef typename Value<TStatement>::Type TValue;
	typedef typename Iterator<TStatement>::Type TIter;

	// Exclude header and empty lines
	TIter it = begin(stmt);
	String<TValue> _id;
	for(;!atEnd(it);goNext(it)) {
	  if ((*it != '\t') && (*it != ' ') && (*it != '\n') && (*it != '\r')) {
	    append(_id, *it);
	  } else {
	    // Exclude any graph, subgraph, node and edge processing attributes
	    if ((_id == "graph") || (_id == "node") || (_id == "edge") || (_id == "subgraph") || (length(_id)<1)) {
	      clear(stmt);
	      return;
	    } else break; 
	  }
	}

	// Process Edges
	it = begin(stmt);
	clear(_id);
	_id = "00";
	unsigned int pos = 0;
	for(;!atEnd(it);goNext(it), ++pos) {
	  _id[pos % 2] = *it;
	  if ((_id == "--") || (_id == "->")) {
	    //std::cout << stmt << std::endl;
	    _processEdgeStatement(g, stmt, nodeMap, edgeMap, pos - 1, nodeIdMap);
	    clear(stmt);
	    return;
	  }
	}

	// Process nodes
	//std::cout << stmt << std::endl;
	_processNodeStatement(g, stmt, nodeMap, edgeMap, nodeIdMap);
	clear(stmt);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TNodeAttributes, typename TEdgeAttributes>
void read(TFile & file,
		  Graph<TSpec>& g,
		  TNodeAttributes& nodeMap,
		  TEdgeAttributes& edgeMap,
		  DotDrawing) 
{
//IOREV _batchreading_ uses custom EOL code
	typedef Graph<TSpec> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef std::map<String<TValue>, TVertexDescriptor> TMap;
	TMap nodeIdMap;

	TValue c;
	String<TValue> stmt;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		
		if (c == ';') _processStatement(g,stmt, nodeMap, edgeMap, nodeIdMap);
		else if ((c == '\n') ||
				(c == '\r')) {
					clear(stmt);
		}
		else append(stmt,c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec>
void read(TFile & file,
		  Graph<TSpec>& g,
		  DotDrawing) 
{
//IOREV _batchreading_ 
	String<String<char> > nodeMap;
	String<String<char> > edgeMap;
	read(file,g,nodeMap,edgeMap,DotDrawing());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
