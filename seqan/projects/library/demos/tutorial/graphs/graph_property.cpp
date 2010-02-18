#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;

int main ()
{
	typedef unsigned int TCargo;
	typedef Graph<Undirected<TCargo> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

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
// FRAGMENT(definition-property-map)
	typedef String<char> TCityName;
	typedef String<TCityName> TProperties;

	TProperties cityNames;
	resizeVertexMap(g, cityNames);
// FRAGMENT(enter-properties)
	assignProperty(cityNames, vertBerlin, "Berlin");
	assignProperty(cityNames, vertHamburg, "Hamburg");
	assignProperty(cityNames, vertMuenchen, "Munich");
	assignProperty(cityNames, vertMainz, "Mainz");
	assignProperty(cityNames, vertHannover, "Hannover");
// FRAGMENT(iterate-and-output)
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(g);

	for(;!atEnd(itV);goNext(itV)) {
		::std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << ::std::endl;
	}

	return 0;
}
