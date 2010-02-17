	typedef Size<TGraph>::Type TSize;
	InternalMap<TCargo> cargoMap;
	resizeEdgeMap(g,cargoMap);
	String<TSize> predMap;
	String<TSize> distMap;
