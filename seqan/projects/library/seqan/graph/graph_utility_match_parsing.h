#ifndef SEQAN_HEADER_GRAPH_UTILITY_MATCH_PARSING_H
#define SEQAN_HEADER_GRAPH_UTILITY_MATCH_PARSING_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TFragment, typename TSpec, typename TSize>
void 
read(TFile & file,
	 String<TFragment, TSpec>& matches,
	 TSize const minMatchSize,
	 AtacMatches) 
{
	SEQAN_CHECKPOINT
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	TSize count = length(matches);
	while (!_streamEOF(file)) {
		if (c == '!') {
			// Ignore comments
			_parse_skipLine(file, c);
		} else if (c == '/') {
			// Ignore comments
			_parse_skipLine(file, c);
		} else {
			_parse_readWord(file, c);
			_parse_skipWhitespace(file,c);
			_parse_readWord(file, c);
			_parse_skipWhitespace(file,c);
			_parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			// Ignore the single dot
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			// Get the identifier
			char identifier1 = c;
			// Ignore the colon
			c = _streamGet(file);
			c = _streamGet(file);
			int seqId1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			int pos1 = _parse_readNumber(file, c);
			int len1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			// Ignore reverse matches
			if (c == '-') {
				_parse_skipLine(file, c);
				continue;
				//exit(-1);
			}
			if (len1 < minMatchSize) {
				_parse_skipLine(file, c);
				continue;
			}
			_parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			// Get the identifier
			char identifier2 = c;
			// Ignore the colon
			c = _streamGet(file);
			c = _streamGet(file);
			int seqId2 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			int pos2 = _parse_readNumber(file, c);
			int len2 = _parse_readNumber(file, c);
			SEQAN_ASSERT(len1 == len2)
			// Ignore reverse matches
			_parse_skipWhitespace(file,c);
			if (c == '-') {
				_parse_skipLine(file, c);
				continue;
				//exit(-1);
			}
			_parse_skipLine(file, c);
			std::cout << identifier1 << ":" << seqId1 << "," << pos1 << "," << len1 << "," << identifier2 << ":" << seqId2 << "," << pos2 << "," << len2 << std::endl;
			if (identifier1 == 'W') seqId1 += 24;
			else if (identifier1 == 'B') seqId1 += 48;
			if (identifier2 == 'W') seqId2 += 24;
			else if (identifier2 == 'B') seqId2 += 48;
			push_back(matches, TFragment(seqId1,pos1,seqId2,pos2,len1));
			//TStringSet strSet;
			//assignValueById(strSet, initialSet, seqId1);
			//assignValueById(strSet, initialSet, seqId2);
			//assignStringSet(matches[count], strSet);
			//addEdge(matches[count], addVertex(matches[count], seqId1, pos1, len1), addVertex(matches[count], seqId2, pos2, len2));
			++count;
			SEQAN_TASSERT(len1 == len2)
		}
	}
	std::cout << std::endl;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
