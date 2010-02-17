// The lines FRAGMENT(fragment-name) mark the begin of fragments for being
// included in the documentation.  You can ignore them.
// FRAGMENT(includes)
#include <iostream>
#include <seqan/graph_types.h>
using namespace seqan;
// FRAGMENT(main-typedefs)
int main ()
{
	typedef unsigned int TCargo;
	typedef Graph<Undirected<TCargo> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
// FRAGMENT(main-graph-construction)
	TGraph g;

	TVertexDescriptor vertBerlin = addVertex(g);
	TVertexDescriptor vertHamburg = addVertex(g);
	TVertexDescriptor vertHannover = addVertex(g);
	TVertexDescriptor vertMainz = addVertex(g);
	TVertexDescriptor vertMuenchen = addVertex(g);

	addEdge(g, vertBerlin, vertHamburg, 289);
	addEdge(g, vertBerlin, vertHannover, 286);
	addEdge(g, vertBerlin, vertMainz, 573);
	addEdge(g, vertBerlin, vertMuenchen, 586);
	addEdge(g, vertHannover, vertMuenchen, 572);
	addEdge(g, vertHamburg, vertMainz, 521);
// FRAGMENT(main-graph-io)
	FILE* strmWrite = fopen("graph.dot", "w");
	write(strmWrite, g, DotDrawing());
	fclose(strmWrite);

	return 0;
}
