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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_FILE_EMBL_H
#define SEQAN_HEADER_FILE_EMBL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Embl
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Embl:EMBL format for sequences.
*/
struct TagEmbl_;
typedef Tag<TagEmbl_> const Embl;



//////////////////////////////////////////////////////////////////////////////
// FileReader Iterator
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Embl, TFile2, TSpec> > & it, 
bool skip_meta = true)
{
SEQAN_CHECKPOINT
	if (skip_meta)
	{
		while (true)
		{
			if (_streamEOF(host(it)))
			{
				it.data_eof = true;
				return;
			}
			if (it.data_char == '/')
			{//end of record
				_stream_skipLine(host(it), it.data_char);
				it.data_eof = true;
				return;
			}
			if (it.data_char == ' ')
			{
				break;
			}
			//skip meta line
			_stream_skipLine(host(it), it.data_char);
		}
	}

	//find first character
	while (true)
	{
		if (_streamEOF(host(it)))
		{
			it.data_eof = true;
			return;
		}
		if ((it.data_char != ' ') && ((it.data_char < '0') || (it.data_char > '9')))
		{
			if ((it.data_char != '\n') && (it.data_char != '\r')) break;

			it.data_char = _streamGet(host(it));
			if (it.data_char == '/')
			{//end of record
				_stream_skipLine(host(it), it.data_char);
				it.data_eof = true;
				return;
			}
		}
		else
		{
			it.data_char = _streamGet(host(it));
		}
	}

//	it.data_file_pos = _streamTellG(host(it));
	it.data_file_pos -= 1;
	it.data_eof = _streamEOF(host(it));
}


template <typename TFile, typename TFile2, typename TSpec>
inline void
goNext(Iter<TFile, FileReader<Embl, TFile2, TSpec> > & it)
{
SEQAN_CHECKPOINT
	do
	{
		it.data_char = _streamGet(host(it));
		if (_streamEOF(host(it)))
		{
			it.data_eof = true;
			return;
		}
		it.data_file_pos += 1;

		if ((it.data_char == '\n') || (it.data_char == '\r'))
		{//linebreak detected: find begin of next line
			do
			{
				it.data_char = _streamGet(host(it));
				if (_streamEOF(host(it)))
				{
					it.data_eof = true;
					return;
				}
				it.data_file_pos += 1;
			} while ((it.data_char == '\n') || (it.data_char == '\r'));

			if (it.data_char == '/')
			{//end of record
				_stream_skipLine(host(it), it.data_char);
				_streamUnget(host(it));
				it.data_eof = true;
				return;
			}
		}
	} while ((it.data_char == ' ') || ((it.data_char >= '0') && (it.data_char <= '9')));
}

//////////////////////////////////////////////////////////////////////////////
// File Format Access Function
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
read(TFile & file,
	 TData & data,
	 Embl)
{
SEQAN_CHECKPOINT
	Iter<TFile, FileReader<Embl> > it(file);

	clear(data);
	while (!atEnd(it))
	{
		appendValue(data, getValue(it));
		goNext(it);
	}
}

template <typename TFile, typename TData, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Embl)
{
SEQAN_CHECKPOINT
	typename Size<TData>::Type siz = length(data);
	Iter<TFile, FileReader<Embl> > it(file);

	clear(data);
	while (!atEnd(it) && (siz < limit))
	{
		appendValue(data, getValue(it));
		goNext(it);
	}
	while (!atEnd(it))
	{
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TMeta>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 Embl)
{
SEQAN_CHECKPOINT
	typedef typename Value<TMeta>::Type TValue;

	clear(meta);
	if (_streamEOF(file))
	{
		return;
	}

	TValue c = _streamGet(file);

	while (!_streamEOF(file))
	{
		if (c == ' ')
		{//end of meta data
			_streamUnget(file);
			return;
		}
		if (c == '/')
		{//end of record
			_stream_skipLine(file, c);
			_streamUnget(file);
			return;
		}

		_stream_appendLine(file, meta, c);
		appendValue(meta, '\n');
	}
}

//////////////////////////////////////////////////////////////////////////////



/**
.Function.readLineType:
..cat:File
..summary:Reads the information belonging to the two-character line code specified.
..signature:readLineType(file,data,key,Embl);
..param.file:The input file or string.
...remarks:This function works on an open file stream or on the string data obtained from calling Function.readMeta
..param.data:The target container that will be filled.
..param.key:The two-character code specifying the file entry to be read, e.g. "AC" for the acession number line or "DE" for the description line. 
..see:Function.readMeta
..see:Function.readFeature
*/
template<typename TFile, typename TData, typename TKey>
inline void
readLineType(TFile & file,
			 TData & data,
			 TKey key,
			 Embl)

{
SEQAN_CHECKPOINT

	//this function is meant to be used for two letter codes only 
	SEQAN_TASSERT(length(key)==2); 

	typedef typename Value<TFile>::Type TValue;
	typedef typename Position<TFile>::Type TPosition;

	clear(data);
	if(_streamEOF(file))
		return;
	
	TPosition pos = _streamTellG(file);
	TValue c = _streamGet(file);
	while (!_streamEOF(file))
	{
		if(c == '/')
		{
			_streamSeekG(file,pos);
			return;
		}
		if(c == key[0])
		{
			c = _streamGet(file);
			if(c == key[1])
			{
				for(unsigned int i = 0; i < 4; ++i)
					c = _streamGet(file);
				_stream_appendLine(file, data, c);
				while(!_streamEOF(file) && _stream_readWord(file,c) == key)
				{
					appendValue(data, '\n');
					for(unsigned int i = 0; i < 3; ++i)
						c = _streamGet(file);
					_stream_appendLine(file, data, c);
				}
				_streamSeekG(file,pos);
				return;
			}
		}
		_stream_skipLine(file, c);
	}

	_streamSeekG(file,pos);
	
}



template<typename TData, typename TValue, typename TSpec, typename TKey>
inline void
readLineType(String<TValue,TSpec> & meta,
			 TData & data,
			 TKey key,
			 Embl)

{
SEQAN_CHECKPOINT

	//this function is meant to be used for two letter codes only 
	SEQAN_TASSERT(length(key)==2); 

	typedef typename Iterator<String<TValue,TSpec>,Standard>::Type TIterator;
	typedef typename Position<String<TValue,TSpec> >::Type TPosition;

	clear(data);
	if(empty(meta))
		return;
	
	TIterator it = begin(meta,Standard());
	TIterator end_it = end(meta,Standard());

	while (it != end_it)
	{
		if(*it == '/')
			return;

		if(*it == key[0])
		{
			++it;
			if(*it == key[1])
			{
				it+=4;
				_string_appendLine(meta, data, it);
				while(it!=end_it && *it==key[0] && *(++it)==key[1])
				{
					appendValue(data, '\n');
					it+=4;
					_string_appendLine(meta, data, it);
				}
				return;
			}
		}
		_string_skipLine(meta, it);
	}

	
}


/**
.Function.readFeature:
..cat:File
..summary:Finds the first feature specified by 'key' starting from position 'start' in the feature table (the feature table can be
obtained by calling readLineType with the two-character code "FT").
..signature:readFeature(ft_string,start,data,key,Embl);
..param.ft_string:The feature table.
..param.start:Position in feature table where search begins.
..param.data:The target container that will be filled.
..param.key:The key word specifying the feature to be read, e.g. "mRNA" or "CDS".
..return:The position behind the feature if found, 0 otherwise.
..see:Function.readMeta
..see:Function.readFeature
*/
//read parts of feature table (those that belong to key)
template<typename TData, typename TKey, typename TString>
inline typename Position<TString>::Type
readFeature(TString & str,
			typename Position<TString>::Type start_pos,
			TData & data,
			TKey key,
			Embl)

{
SEQAN_CHECKPOINT

	typedef typename Iterator<TString,Standard>::Type TIterator;
	typedef typename Position<TString>::Type TPosition;

	clear(data);
	if(empty(str) || start_pos >= length(str))
		return 0;
	
	TIterator it = iter(str,start_pos,Standard());
	TIterator end_it = end(str,Standard());

	while (it != end_it)
	{
		if(*it == key[0])
		{
			++it;
			bool found = true;
			for(unsigned int i = 1; i < length(key); ++i)
			{
				if(key[i] != *it)
				{
					found = false;
					break;
				}
				++it;
			}
			if(found)
			{
				_string_skipWhitespace(str,it);
				_string_appendLine(str, data, it);
				while(it!=end_it && *it == ' ')
				{
					appendValue(data, '\n');
					_string_skipWhitespace(str,it);
					_string_appendLine(str, data, it);
				}
				return position(it,str);
			}
		}
		_string_skipLine(str, it);
	}

	return 0;
	
}
	


//////////////////////////////////////////////////////////////////////////////
template <typename TFile>
inline void
goNext(TFile & file,
	   Embl)
{
SEQAN_CHECKPOINT
	typedef typename Value<TFile>::Type TValue;

	if (_streamEOF(file))
	{
		return;
	}

	while (!_streamEOF(file))
	{
		TValue c = _streamGet(file);
		if (c == '/')
		{//end of record
			_stream_skipLine(file, c);
			_streamUnget(file);
			return;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
/*
template <typename TFile>
inline void
length(TFile & file,
	   Embl)
{
SEQAN_CHECKPOINT
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
write(TFile & file,
	  TData & data,
	  Embl)
{
SEQAN_CHECKPOINT
	enum
	{
		BLOCK_SIZE = 10,
		BLOCKS_PER_LINE = 6,
		FIRST_INDENT = 5
	};
	char const * NUM_BLOCK_FORMAT = "%9d";

	typedef typename Size<TData>::Type TSize;
	typedef typename Iterator<TData, Standard>::Type TIterator;

	TSize count = 0;
	int block_count = 0;
	int char_in_block_count = 0;
	TIterator it = begin(data, Standard());
	TIterator it_end = end(data, Standard());

	while (it != it_end)
	{
		//write indent
		for (int j = 0; j < FIRST_INDENT; ++j) _streamPut(file, ' ');

		//write rest of line
		while (true)
		{
			//write next character
			if (it != it_end)
			{
				_streamPut(file, *it);
				++it;
				++count;
			}
			else
			{//no more chars left: fill up with ' '
				_streamPut(file, ' ');
			}
			++char_in_block_count;

			if (char_in_block_count >= BLOCK_SIZE)
			{//end of block
				char_in_block_count = 0;
				_streamPut(file, ' ');
				++block_count;

				if (block_count >= BLOCKS_PER_LINE)
				{//enough blocks: write end of line
					block_count = 0;
					_streamPutInt(file, count, NUM_BLOCK_FORMAT);
					_streamPut(file, '\n');
					break; //end of line
				}
			}
		}
	}

	write(file, "//\n");
}

template <typename TFile, typename TData, typename TMeta>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  Embl)
{
SEQAN_CHECKPOINT
	write(file, meta);
	write(file, data, Embl());
}

//////////////////////////////////////////////////////////////////////////////

/*leere Vorlage
template <typename TFile, typename TData>
inline void
read(TFile & file,
	 TData & data,
	 Embl)
{
SEQAN_CHECKPOINT
}

template <typename TFile, typename TData, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Embl)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TMeta>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 Embl)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
inline void
goNext(TFile & file,
	   Embl)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
inline void
length(TFile & file,
	   Embl)
{
SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData>
inline void
write(TFile & file,
	  TData & data,
	  Embl)
{
SEQAN_CHECKPOINT
}
template <typename TFile, typename TData, typename TMeta>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  Embl)
{
SEQAN_CHECKPOINT
}
*/

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
