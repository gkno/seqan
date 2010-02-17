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


	typedef String<char> TCityName;
	typedef String<TCityName> TProperties;
	TProperties cityNames;
	resizeVertexMap(g, cityNames);
	assignProperty(cityNames, vertBerlin, "Berlin");
	assignProperty(cityNames, vertHamburg, "Hamburg");
	assignProperty(cityNames, vertMuenchen, "Munich");
	assignProperty(cityNames, vertMainz, "Mainz");
	assignProperty(cityNames, vertHannover, "Hannover");

	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(g);
	for(;!atEnd(itV);goNext(itV)) {
		::std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << ::std::endl;
	}
// part 2
	typedef Size<TGraph>::Type TSize;
	InternalMap<TCargo> cargoMap;
	resizeEdgeMap(g,cargoMap);
	String<TSize> predMap;
	String<TSize> distMap;
// part 3
	dijkstra(g,vertHannover,cargoMap,predMap,distMap);
// part 4
	TVertexIterator itV2(g);
	while(!atEnd(itV2)) {
		::std::cout << "Shortest path from " << property(cityNames, vertHannover) << " to " << property(cityNames, value(itV2)) << ": ";
		::std::cout << property(distMap, value(itV2)) << ::std::endl;
		goNext(itV2);
	}

	return 0;
}
