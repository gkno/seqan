#ifndef SEQAN_HEADER_GRAPH_UTILITY_MATCH_PARSING_H
#define SEQAN_HEADER_GRAPH_UTILITY_MATCH_PARSING_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template<typename TFile>
void 
read(TFile & file,
	 AtacMatches) 
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;


	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	
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
			char identifier = c;
			// Ignore the colon
			c = _streamGet(file);
			int seqId = _parse_readNumber(file, c);
			std::cout << c << ":" << seqId << std::endl;
			_parse_skipLine(file, c);
		}
	}
		
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
