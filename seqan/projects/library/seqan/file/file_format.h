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

#ifndef SEQAN_HEADER_FILE_FORMAT_H
#define SEQAN_HEADER_FILE_FORMAT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format:
..summary:A file format.
*/


//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//Base Class for all FileFormat classes
//////////////////////////////////////////////////////////////////////////////

/**
.Class.FileFormat:
..cat:Input/Output
..summary:Object that stores a file format.
..signature:FileFormat<File, Data [, Format [, Meta] ]>
..see:Tag.File Format
*/

template <
	typename TFile, 
	typename TData,
	typename TMeta,
	typename TFormat = void >
struct FileFormat:
	public FileFormat<TFile, TData, TMeta, void>
{
public:
	typedef typename Size<TData>::Type TSize;

	FileFormat() {}
	FileFormat(FileFormat const &) {}
	~FileFormat() {}
	FileFormat const & operator =(FileFormat const &) {}

	inline void * 
	formatID_() const
	{
SEQAN_CHECKPOINT
		return _ClassIdentifier<TFormat>::getID();
	}

	virtual void
	read_(TFile & file, TData & data) const
	{
SEQAN_CHECKPOINT
		read(file, data, TFormat());
	}
	virtual void
	read_(TFile & file, TData & data, TSize limit) const
	{
SEQAN_CHECKPOINT
		read(file, data, limit, TFormat());
	}

	virtual void
	readMeta_(TFile & file, TMeta & meta) const
	{
SEQAN_CHECKPOINT
		readMeta(file, meta, TFormat());
	}

	virtual void
	goNext_(TFile & file) const
	{
SEQAN_CHECKPOINT
		goNext(file, TFormat());
	}

	virtual TSize
	length_(TFile & file) const
	{
SEQAN_CHECKPOINT
		length(file, TFormat());
	}

	virtual void
	write_(TFile & file, TData & data) const
	{
SEQAN_CHECKPOINT
		write(file, data, TFormat());
	}
	virtual void
	write_(TFile & file, TData & data, TMeta & meta) const
	{
SEQAN_CHECKPOINT
		write(file, data, meta, TFormat());
	}
};

//____________________________________________________________________________

//base class for all file format classes 

template <typename TFile, typename TData, typename TMeta>
struct FileFormat<TFile, TData, TMeta, void>
{
public:
	typedef typename Size<TData>::Type TSize;

	FileFormat() {}
	FileFormat(FileFormat const &) {}
	~FileFormat() {};
	FileFormat const & operator =(FileFormat const &) {}

	virtual void *
	formatID_() const = 0;

	virtual void
	read_(TFile & file, TData & data) const = 0;
	virtual void
	read_(TFile & file, TData & data, TSize limit) const = 0;

	virtual void
	readMeta_(TFile & file, TMeta & meta) const = 0;

	virtual void
	goNext_(TFile & file) const = 0;

	virtual TSize
	length_(TFile & file) const = 0;

	virtual void
	write_(TFile & file, TData & data) const = 0;
	virtual void
	write_(TFile & file, TData & data, TMeta & meta) const = 0;

};

//////////////////////////////////////////////////////////////////////////////
// Wrapper for functions to virtuals
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void *
formatID(FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	return file_format.formatID_();
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.read:
..cat:Input/Output
..summary:Loads a record from file.
..signature:read(file, data [, meta], format)
..signature:read(file, data [, meta], tag)
..param.file:An input file.
..param.data:A container that gets the data read from $file$.
..param.meta:A container that gets meta data from $file$. (optional)
..param.format:A file format object.
...type:Class.FileFormat.File Format object
..param.tag:A file format tag.
...type:Tag.File Format.File Format tag
..remarks:The result of this operation is stored in $data$.
..remarks:The function leaves $file$ at the position for reading the next record.
..see:Function.assign
*/
template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
read(TFile & file,
	 TData & data,
	 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.read_(file, data);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.read_(file, data, limit);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.readMeta:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.readMeta_(file, meta);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
goNext(TFile & file,
	   FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.goNext_(file);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.length:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
length(TFile & file,
	   FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.length_(file);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.write:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
write(TFile & file,
	  TData & data,
	  FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.write_(file, data);
}
template <typename TFile, typename TData, typename TMeta, typename TFormat>
inline void
write(TFile & file,
	  TData & data,
	  TMeta & meta,
	  FileFormat<TFile, TData, TMeta, TFormat> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.write_(file, data, meta);
}




//////////////////////////////////////////////////////////////////////////////
// Comparison of two FileFormat objects
//////////////////////////////////////////////////////////////////////////////

template <typename TFileLeft, typename TDataLeft, typename TMetaLeft, typename TFormatLeft, typename TFileRight, typename TDataRight, typename TMetaRight, typename TFormatRight>
inline bool
operator == (FileFormat<TFileLeft, TDataLeft, TMetaLeft, TFormatLeft> const & left, 
			 FileFormat<TFileRight, TDataRight, TMetaRight, TFormatRight> const & right)
{
SEQAN_CHECKPOINT
	return formatID(left) == formatID(right);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator == (FileFormat<TFile, TData, TMeta, TFormat> const & left, 
			 Tag<TFormat2> const)
{
SEQAN_CHECKPOINT
	return formatID(left) == _ClassIdentifier<Tag<TFormat2> const>::getID();
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator == (Tag<TFormat2> const,
			 FileFormat<TFile, TData, TMeta, TFormat> const & right)
{
SEQAN_CHECKPOINT
	return _ClassIdentifier<Tag<TFormat2> const>::getID() == formatID(right);
}

//____________________________________________________________________________

template <typename TFileLeft, typename TDataLeft, typename TMetaLeft, typename TFormatLeft, typename TFileRight, typename TDataRight, typename TMetaRight, typename TFormatRight>
inline bool
operator != (FileFormat<TFileLeft, TDataLeft, TMetaLeft, TFormatLeft> const & left, 
			 FileFormat<TFileRight, TDataRight, TMetaRight, TFormatRight> const & right)
{
SEQAN_CHECKPOINT
	return formatID(left) != formatID(right);
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator != (FileFormat<TFile, TData, TMeta, TFormat> const & left, 
			 Tag<TFormat2> const)
{
SEQAN_CHECKPOINT
	return formatID(left) != _ClassIdentifier<Tag<TFormat2> const>::getID();
}

template <typename TFile, typename TData, typename TMeta, typename TFormat, typename TFormat2>
inline bool
operator != (Tag<TFormat2> const,
			 FileFormat<TFile, TData, TMeta, TFormat> const & right)
{
SEQAN_CHECKPOINT
	return _ClassIdentifier<Tag<TFormat2> const>::getID() != formatID(right);
}

//////////////////////////////////////////////////////////////////////////////
// allgemeine Funktionen fuer Streams
//////////////////////////////////////////////////////////////////////////////
//TODO??? Das muss in eine extra Datei

/*DISABLED
.Function.write:
..summary:Writes to stream.
..cat:Input/Output
..signature:write(stream, source)
..signature:write(stream, begin, end)
..param.stream: A stream object.
...type:Adaption."std::iostream"
..param.source: Container that is written to $stream$.
..param.begin: Iterator to the first character of the range.
..param.end: Iterator behind the last character of the range.
..param.limit: The maximal number of charactes written to $stream$. (optional)
..remarks:The content of $source$ is written 'as-is' to $stream$.
??? TODO: save as file format
*/
/*
template <typename TStream, typename TIterator>
inline void
write(TStream & target,
	  TIterator begin_,
	  TIterator end_)
{
	while (begin_ != end_)
	{
		_streamPut(target, convert<char>(*begin_));
		++begin_;
	}
}

//____________________________________________________________________________

template <typename TStream, typename TSource>
inline void
write(TStream & target,
	  TSource const & source)
{
	write(target, begin(source), end(source));
}
//TODO???: Spezialisierungen zum blockweise schreiben bei contiguous strings von char
//Anmerkungen: write wird nach dem zweiten Argument (source) spezialisiert!

//____________________________________________________________________________

template <typename TStream, typename TSource>
inline void
write(TStream & target,
	  TSource const & source,
	  typename Size<TSource>::Type limit_)
{
	if (length(source) > limit_)
	{
		write(target, begin(source), begin(source) + limit_);
	}
	else
	{
		write(target, begin(source), end(source));
	}
}

*/
//////////////////////////////////////////////////////////////////////////////

// Helper function for scanning a stream
// c = next character, pass it to the next call of the function

template <typename TFile, typename TString, typename TChar>
inline void
_stream_appendLine(TFile & file,
				   TString & str,
				   TChar & c)
{
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		appendValue(str, c);

		c = _streamGet(file);
	}
}
//____________________________________________________________________________

template <typename TFile, typename TChar>
inline void
_stream_countLine(TFile & file,
				  TChar & c)

{
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		c = _streamGet(file);
	}
}

//____________________________________________________________________________

template <typename TFile, typename TChar>
inline typename Size<TFile>::Type
_stream_skipLine(TFile & file,
				 TChar & c)

{
	typename Size<TFile>::Type count = 0;
	while (true)
	{
		if (_streamEOF(file)) break;

		if (c == '\r')
		{
			c = _streamGet(file);
			if (c == '\n') 
			{
				c = _streamGet(file);
			}
			break;
		}
		if (c == '\n')
		{
			c = _streamGet(file);
			break;
		}

		++count;

		c = _streamGet(file);
	}

	return count;
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
