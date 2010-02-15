	typedef Iterator<TCounterString>::Type TCounterIter;
	TCounterIter countIt = begin(counter);
	TCounterIter countItEnd = end(counter);
	TSize pos = 0;
	for(;countIt != countItEnd; ++countIt, ++pos) {
		::std::cout << AminoAcid(pos) << ':' << value(countIt) << ::std::endl;
	}


	return 0;
}

