 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_utility_parsing.h 1911 2008-05-02 09:28:04Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_UTILITY_PARSING_H
#define SEQAN_HEADER_GRAPH_UTILITY_PARSING_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// General parsing funtions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline void 
_parse_skipLine(TFile& file, TChar& c)
{
	if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
		c = _streamGet(file);
		return;
	}
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
	}
	c = _streamGet(file);
}
//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TChar>
inline void 
_parse_skipWhitespace(TFile& file, TChar& c)
{
	if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) return;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) break;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_parse_isDigit(TChar const c)
{
	//return (((unsigned) c >=  48) && ((unsigned) c <=  57));
	return ((c == '0') || (c == '1') || (c == '2') || (c == '3') || (c == '4') || 
		    (c == '5') || (c == '6') || (c == '7') || (c == '8') || (c == '9'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_parse_isLetter(TChar const c)
{
	//return ((((unsigned) c >=  97) && ((unsigned) c <=  122)) || (((unsigned) c >=  65) && ((unsigned) c <=  90)));
	return ((c == 'a') || (c == 'b') || (c == 'c') || (c == 'd') || (c == 'e') || 
			(c == 'f') || (c == 'g') || (c == 'h') || (c == 'i') || (c == 'j') ||
			(c == 'k') || (c == 'l') || (c == 'm') || (c == 'n') || (c == 'o') || 
			(c == 'p') || (c == 'q') || (c == 'r') || (c == 's') || (c == 't') ||
			(c == 'u') || (c == 'v') || (c == 'w') || (c == 'x') || (c == 'y') || 
			(c == 'z') || (c == 'A') || (c == 'B') || (c == 'C') || (c == 'D') ||
			(c == 'E') || (c == 'F') || (c == 'G') || (c == 'H') || (c == 'I') || 
			(c == 'J') || (c == 'K') || (c == 'L') || (c == 'M') || (c == 'N') ||
			(c == 'O') || (c == 'P') || (c == 'Q') || (c == 'R') || (c == 'S') || 
			(c == 'T') || (c == 'U') || (c == 'V') || (c == 'W') || (c == 'X') ||
			(c == 'Y') || (c == 'Z'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TChar>
inline bool
_parse_isAlphanumericChar(TChar const c)
{
	return ((_parse_isDigit(c)) || (_parse_isLetter(c)) || (c == '_') || (c == '.') || (c == '-') || (c == '|'));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline int
_parse_readNumber(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isDigit(c)) break;
		append(str, c);
	}
 	return atoi(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline double
_parse_readDouble(TFile & file, TChar& c)
{
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isDigit(c) && (c != '.')) break;
		append(str, c);
	}
 	return atof(toCString(str));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline String<char>
_parse_readIdentifier(TFile & file, TChar& c)
{
	// Read identifier
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isAlphanumericChar(c)) break;
		append(str, c);
	}
	return str;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TChar>
inline void
_parse_readIdentifier(TFile & file, TString& str, TChar& c)
{
	// Read identifier
	append(str, c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isAlphanumericChar(c)) break;
		append(str, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar>
inline String<char>
_parse_readWord(TFile & file, TChar& c)
{
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(str, c);
	}
	return str;
}


//read filename (read line and trim trailing whitespaces)
template<typename TFile, typename TChar>
inline String<char>
_parse_readFilepath(TFile& file, TChar& c)
{
	String<char> str(c);
	if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) {
		c = _streamGet(file);
		return str;
	}
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
		append(str, c);
	}
	typename Iterator<String<char>,Rooted >::Type str_it = end(str);	
	while(str_it != begin(str)) {
		--str_it;
		if(*str_it != ' ' && *str_it != '\t'){
		++str_it;
		break;
		}
	}
	resize(str,position(str_it));
	return str;
}









//////////////////////////////////////////////////////////////////////////////
// Command line parsing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template<typename TKey, typename TValue>
class ConfigOptions {
public:
	std::map<TKey, TValue> option;
	TValue help;

	ConfigOptions() {}
};

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue>
struct Iterator<ConfigOptions<TKey, TValue> >
{	
	typedef std::map<TKey, TValue> TMap;
	typedef typename TMap::iterator Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue>
struct Iterator<ConfigOptions<TKey, TValue> const>
{	
	typedef std::map<TKey, TValue> const TMap;
	typedef typename TMap::const_iterator Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TContainer, typename TSize>
inline void
assignKeys(ConfigOptions<TKey, TValue>& cfgOpt,
		   TContainer& params,
		   TSize numKeys) 
{
	for(TSize positionIt = 0;positionIt < numKeys;++positionIt) cfgOpt.option.insert(std::make_pair(value(params, positionIt), TValue()));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TContainer>
inline void
assignKeys(ConfigOptions<TKey, TValue>& cfgOpt,
		   TContainer& params) 
{
	assignKeys(cfgOpt, params, length(params));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TKey1, typename TValue1>
inline void
assign(ConfigOptions<TKey, TValue>& cfgOpt,
	   TKey1& key1,
	   TValue1& value1)
{
	cfgOpt.option[key1] = value1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TKey1>
inline TValue&
value(ConfigOptions<TKey, TValue>& cfgOpt,
	  TKey1 const& key1)
{
	if (cfgOpt.option.find(key1) == cfgOpt.option.end()) {
		static TValue tmp;
		return tmp;
	}
	else return (cfgOpt.option.find(key1))->second;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TKey1>
inline TValue&
value(ConfigOptions<TKey, TValue> const& cfgOpt,
	  TKey1 const& key1)
{
	if (cfgOpt.option.find(key1) == cfgOpt.option.end()) {
		static TValue tmp;
		return tmp;
	}
	else return ((const_cast<ConfigOptions<TKey, TValue>&>(cfgOpt)).option.find(key1))->second;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TKey1>
inline TValue
getValue(ConfigOptions<TKey, TValue>& cfgOpt,
		 TKey1& key1)
{
	return value(cfgOpt, key1);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue>
inline typename Iterator<ConfigOptions<TKey, TValue> >::Type
begin(ConfigOptions<TKey, TValue>& cfgOpt)
{
	return cfgOpt.option.begin();
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue>
inline typename Iterator<ConfigOptions<TKey, TValue> >::Type
end(ConfigOptions<TKey, TValue>& cfgOpt)
{
	return cfgOpt.option.end();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TKey, typename TValue, typename TIDString>
inline void
write(TFile & target,
	  ConfigOptions<TKey, TValue>& cfgOpt,
	  TIDString const &,
	  Raw)
{
	typedef typename Iterator<ConfigOptions<TKey, TValue> >::Type TIter;
	TIter it = begin(cfgOpt);
	TIter itEnd = end(cfgOpt);
	for(;it!=itEnd;++it) {
		_streamWrite(target, (*it).first);
		_streamPut(target, '=');
		_streamWrite(target, (*it).second);
		_streamPut(target, '\n');
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TKey, typename TValue>
inline TStream &
operator << (TStream & target, 
			 ConfigOptions<TKey, TValue>& source)
{
	write(target, source);
	return target;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue,typename TValue1>
inline void
assignHelp(ConfigOptions<TKey, TValue>& cfgOpt,
	   	   TValue1& value1)
{
	cfgOpt.help = value1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue>
inline TValue&
valueHelp(ConfigOptions<TKey, TValue> const& cfgOpt)
{
	return (const_cast<ConfigOptions<TKey, TValue>&>(cfgOpt)).help;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue>
inline TValue&
valueHelp(ConfigOptions<TKey, TValue>& cfgOpt)
{
	return cfgOpt.help;
}


//////////////////////////////////////////////////////////////////////////////


template<typename TKey, typename TValue>
inline bool
parseCmdLine(int argc, const char *argv[], ConfigOptions<TKey, TValue>& cfgOpt) {
	TKey lastKey = "";
	for(int i = 1; i < argc; ++i) {
		if (argv[i][0] == '-') {
			TKey key = &argv[i][1];
			if (cfgOpt.option.find(key) != end(cfgOpt)) {
				if (i + 1 == argc) {
					std::cerr << valueHelp(cfgOpt) << std::endl;
					return false;
				}
				++i;
				cfgOpt.option[key] = argv[i];
				lastKey = key;
				continue;
			} else {
				std::cerr << valueHelp(cfgOpt) << std::endl;
				return false;
			}
		} else {
			if (length(lastKey)) {
				TValue tmp = cfgOpt.option[lastKey];
				append(tmp, &argv[i][0]);
				cfgOpt.option[lastKey] = tmp;
			} else {
				std::cerr << valueHelp(cfgOpt) << std::endl;
				return false;
			}
		}
	}
	return true;
}


template <typename T, typename TString>
inline T
_stringToNumber(TString& str)
{
	T val = 0;
	//reserve(str, length(str) + 1);
	std::stringstream ssStream1(toCString(str));
	ssStream1 >> val; 
	return val;
}




//////////////////////////////////////////////////////////////////////////////
// File reading and timing functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TPath, typename TStringSet, typename TNames>
inline unsigned int
_loadSequences(TPath const& in_path, 
					  TStringSet& origStrSet,
					  TNames& names)
{
	typedef typename Size<TStringSet>::Type TSize;
	
	// Count sequences and read names
	TSize seqCount = 0;
	std::ifstream file;
	std::stringstream input;
	input << in_path;
	file.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
	if (!file.is_open()) return 0;
	while (!_streamEOF(file)) {
		String<char> id;
		readID(file, id, Fasta());
		appendValue(names, id);
		goNext(file, Fasta());
		++seqCount;
	}

	// Load sequences
	file.clear();
	file.seekg(0, std::ios_base::beg);
	resize(origStrSet, seqCount);
	TSize count = 0;
	for(TSize i = 0; (i < seqCount) && !_streamEOF(file); ++i) 	{
		read(file, origStrSet[i], Fasta());
		count += length(origStrSet[i]);
	}
    file.close();

	return count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TText>
inline void
_alignTiming(std::clock_t& startTime,
			 TText const& text)
{
	std::clock_t endTime=clock();
	double time=((float)(endTime-startTime)/CLOCKS_PER_SEC);
	startTime = endTime;
	std::cout << text << time << " sec" << std::endl;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
