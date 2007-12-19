#ifndef SEQAN_HEADER_BLAST_PARSING_H
#define SEQAN_HEADER_BLAST_PARSING_H

//SEQAN_NO_DDDOC: do not generate documentation for this file

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Internal Blast Parsing Functions
// remark: uses a lot of functions from graph_utility_parsing.h
//////////////////////////////////////////////////////////////////////////////


	//die beiden raus hier, gehören in streamAlgorithms.h

		template <typename TStream>
		inline void
		_streamPutFloatBlast(TStream & target,
					float number, 
					char const * format_string)
		{
		SEQAN_CHECKPOINT
			char str[32];
			sprintf(str, format_string, number);
			_streamWrite(target, str);
		}

		template <typename TStream>
		inline void
		_streamPutFloatBlast(TStream & target,
					float number)
		{
		SEQAN_CHECKPOINT
			_streamPutFloatBlast(target, number, "%f");
		}



//template<typename TFile, typename TChar>
//inline void 
//_parse_skipLine(TFile& file, TChar& c)
//{
//	if (c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) return;
//	while (!_streamEOF(file)) {
//		c = _streamGet(file);
//		if (c == '\n'|| (c == '\r' && _streamPeek(file) != '\n')) break;
//	}
//	c = _streamGet(file);
//}
//


/////////////////////////////////////////////////////////////////////////////////
// read expect values into double 
// steht am ende hinter der zahl
template<typename TFile, typename TChar>
inline double
_parse_readEValue(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
	// Read number
	String<char> str(c);
	bool e = false;
	double val1;
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if(!e && c == 'e'){
			e = true;
			val1 = atof(toCString(str));
			c = _streamGet(file);
			resize(str,0);
		}
		if (!_parse_isDigit(c) && c != '.' && c != '-') break;
		append(str, c);
	}
	if(e)
	{
		return val1 * pow((double)10.0,(double)atof(toCString(str)));
	}	
 	else 
		return (double)atof(toCString(str));
}

/////////////////////////////////////////////////////////////////////////////////
// read alignment string (letters and gaps)
// steht am ende dahinter 
template<typename TFile, typename TChar>
inline String<char>
_parse_readAlignmentString(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
	// Read word
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c) && !(c=='-')) break;
		append(str, c);
	}
	return str;
}


/////////////////////////////////////////////////////////////////////////////////
// read floating point number
// steht am ende dahinter 
template<typename TFile, typename TChar>
inline float
_parse_readFloat(TFile & file, TChar& c)
{
SEQAN_CHECKPOINT
	// Read number
	String<char> str(c);
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (c != '.' && c != ',' && !_parse_isDigit(c)) break;
		append(str, c);
	}
 	return atof(toCString(str));
}


// steht am ende dahinter!!!
template<typename TFile, typename TChar, typename TSize>
inline String<char>
_parse_readWord(TFile & file, TChar& c, TSize max_len)
{
	// Read word
	String<char> str(c);
	--max_len;
	TSize i = 0;
	while (!_streamEOF(file) ) {
		c = _streamGet(file);
		if (!_parse_isLetter(c) || i >= max_len) break;
		append(str, c);
		++i;
	}
	return str;
}


//// steht am ende dahinter!!!
//template<typename TFile, typename TChar>
//inline String<char>
//_parse_readIdLine(TFile & file, TChar& c)
//{
//	// Read word
//	String<char> str(c);
//	
//	typename Position<TFile>::Type start_pos, end_pos;
//	start_pos = _streamTellG(file);
//	String<char> search = "Length";
//	_parse_untilBeginLine(file,c,search,6);
//	end_pos = _streamTellG(file);
//	_streamSeekG(file,start_pos);
//	while (!_streamEOF(file) && _streamTellG(file)!=end_pos) {
//		c = _streamGet(file);
//		if (c == '\n' || c == '\r')	append(str, ' ');
//		else append(str, c);
//	}
//	return str;
//}
//

/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with character x (skip whitespaces)
// zeigt am ende darauf!!!
template<typename TFile, typename TChar>
inline bool
_parse_untilBeginLine(TFile & file, TChar& c, TChar x)
{
SEQAN_CHECKPOINT
	_parse_skipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file) && c != x){
		_parse_skipLine(file, c);
		_parse_skipWhitespace(file,c);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with word
//zeigt am ende dahinter!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parse_untilBeginLine(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	_parse_skipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parse_readWord(file,c,len))
				break;
		_parse_skipLine(file, c);
		_parse_skipWhitespace(file,c);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with word (parse no more than num_lines lines)
//zeigt am ende dahinter!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parse_untilBeginLine(TFile & file, TChar& c, String<TChar> & word, TSize len, TSize num_lines)
{
SEQAN_CHECKPOINT
	_parse_skipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	TSize i = 0;
	bool found = false;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parse_readWord(file,c,len))
			{
				found = true;
				break;
			}
		if(i >= num_lines)
			break;
		++i;
		_parse_skipLine(file, c);
		_parse_skipWhitespace(file,c);
	}
	if(!_streamEOF(file) && found) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until line begins with one of the characters in string x (skip whitespaces)
//zeigt am ende darauf!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parse_untilBeginLineOneOf(TFile & file, TChar& c, String<TChar> & x, TSize len)
{
SEQAN_CHECKPOINT
	_parse_skipWhitespace(file,c);
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	bool found = false;
	while (!_streamEOF(file)){
		for(int i = 0; i < len; ++i)
			if(c == x[i]) 
			{
				found = true;
				break;
			}
		if(found) break;
		_parse_skipLine(file, c);
		_parse_skipWhitespace(file,c);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until c == x
//zeigt am ende darauf!
template<typename TFile, typename TChar>
inline bool
_parse_until(TFile & file, TChar& c, TChar x)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file) && c != x){
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}



/////////////////////////////////////////////////////////////////////////////////
//parse until word
//zeigt am ende dahinter!
template<typename TFile, typename TChar, typename TSize>
inline bool
_parse_until(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
			if(word == _parse_readWord(file,c,len))
				break;
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse until c == x or new line
//zeigt am ende darauf!
template<typename TFile, typename TChar>
inline bool
_parse_lineUntil(TFile & file, TChar& c, TChar x)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file) && c != x){
		if (c == '\n' || c == '\r')
		{
			_streamSeekG(file,pos);
			c = c_before;
			return false;
		}
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}


/////////////////////////////////////////////////////////////////////////////////
//parse this line until word
//zeigt am ende hinter wort if true, oder auf ende der zeile
template<typename TFile, typename TChar, typename TSize>
inline bool
_parse_lineUntil(TFile & file, TChar& c, String<TChar> & word, TSize len)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	while (!_streamEOF(file)){
		if(c == word[0])
		{	if(word == _parse_readWord(file,c,len))
				break;
		}
		else if (c == '\n' || c == '\r')
			{
				_streamSeekG(file,pos);
				c = c_before;
				return false;
			}
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}



/////////////////////////////////////////////////////////////////////////////////
//parse until c is any of the characters in x
//zeigt am ende darauf
template<typename TFile, typename TChar, typename TSize>
inline bool
_parse_untilOneOf(TFile & file, TChar& c, String<TChar> x, TSize len)
{
SEQAN_CHECKPOINT
	typename Position<TFile>::Type pos = _streamTellG(file);
	TChar c_before = c;
	bool found = false;
	while (!_streamEOF(file)){
		for(int i = 0; i < len; ++i)
			if(c == x[i]) 
			{
				found = true;
				break;
			}
		if(found) break;
		c = _streamGet(file);
	}
	if(!_streamEOF(file)) return true;
	_streamSeekG(file,pos);
	c = c_before;
	return false;
}

////////////////////////////////////////////////////////////////////////////
// parses query and database name (file should be pointing to the beginning of the file)
template<typename TFile, typename TChar>
typename Position<TFile>::Type
_parse_readQueryAndDBName(TFile & file,
						  TChar & c,
						  String<char> & query_name,
						  String<char> & db_name)
{
SEQAN_CHECKPOINT
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;

	TChar c_before = c;
	TPosition act_pos = _streamTellG(file);
	TPosition query_pos,db_pos;

	//String<char> delim = "DQ";
	//_parse_untilBeginLineOneOf(file,c,delim,2);


	String<char> query = "Query";
	if(_parse_untilBeginLine(file,c,query,6))
		query_pos = _streamTellG(file);
	else
		return act_pos;
	_streamSeekG(file,act_pos);
	c = c_before;
	String<char> database = "Database";
	if(_parse_untilBeginLine(file,c,database,8))
		db_pos = _streamTellG(file);
	else
		return act_pos;
	
	
	//getQueryName
	_streamSeekG(file,query_pos);
	_parse_skipWhitespace(file,c);
	c = _streamGet(file);
	_parse_skipWhitespace(file,c);
	query_name = _parse_readWord(file, c);
	while (!_streamEOF(file) && c != '\n' && c != '\r')
		query_name += _parse_readWord(file, c);
	
	//getDBName
	_streamSeekG(file,db_pos);
	c = _streamGet(file);
	_parse_skipWhitespace(file,c);
	db_name = _parse_readWord(file, c);
	while (!_streamEOF(file) && c != '\n' && c != '\r')
		db_name += _parse_readWord(file, c);
	_parse_skipWhitespace(file,c);
		
	c = _streamGet(file);

	return _streamTellG(file); 
	
}








}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
