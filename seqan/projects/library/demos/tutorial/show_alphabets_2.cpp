template<typename TAlphabet> void showAllLetterOfMyAlphabet(TAlphabet const&) {
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	for(TSize i = 0; i < alphSize; ++i) {
		::std::cout << TAlphabet(i) << ',';
	}
	::std::cout << ::std::endl;
}
