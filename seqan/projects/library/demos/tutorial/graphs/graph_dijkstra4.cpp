	TVertexIterator itV2(g);
	while(!atEnd(itV2)) {
		::std::cout << "Shortest path from " << property(cityNames, vertHannover) << " to " << property(cityNames, value(itV2)) << ": ";
		::std::cout << property(distMap, value(itV2)) << ::std::endl;
		goNext(itV2);
	}

	return 0;
}
