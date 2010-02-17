	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator itV(g);

	for(;!atEnd(itV);goNext(itV)) {
		::std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << ::std::endl;
	}

	return 0;
}
