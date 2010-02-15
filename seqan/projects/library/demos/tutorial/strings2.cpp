	typedef Iterator<TAminoAcidString>::Type TIter;

	TIter it = begin(sourceSeq);
	TIter itEnd = end(sourceSeq);
	for(;it != itEnd; goNext(it)) {
		if (value(it) == 'R') value(it) = 'A';
		::std::cout << value(it) << ',';
	}
	::std::cout << ::std::endl;
