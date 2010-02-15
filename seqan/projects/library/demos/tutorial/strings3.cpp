	typedef Size<TAminoAcidString>::Type TSize;
	typedef String<TSize> TCounterString;
	TCounterString counter;
	TSize alphSize = ValueSize<AminoAcid>::VALUE;
	fill(counter, alphSize, 0);
	it = begin(sourceSeq);
	for(;it != itEnd; goNext(it)) {
		value(counter, ordValue(value(it))) += 1;
	}
