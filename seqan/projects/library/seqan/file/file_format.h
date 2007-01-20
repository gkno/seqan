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

/**
.Metafunction.Meta:
..summary:Metadata storage class.
..signature:Meta<Value, Format>::Type
..param.Value:The value type of the file.
..param.Format:A file format tag.
...value:Tag.File Format
..returns.param.Type:Datastructure for storing the metadata of a record in files of format $Format$.
..see:Tag.File Format
..see:Class.FileFormat
..see:Class.Metadata
*/

template <typename TValue, typename TFormat>
struct Meta
{
	typedef Metadata<TValue, TFormat> Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.RecordID:
..summary:ID type of file records.
..signature:RecordID<Value, Format>::Type
..param.Value:The value type of the file.
..param.Format:A file format tag.
...value:Tag.File Format
..returns.param.Type:Datastructure for storing the ID of record in files of format $Format$.
..see:Tag.File Format
..see:Class.FileFormat
*/

template <typename TValue, typename TFormat>
struct RecordID
{
	typedef String<TValue> Type;
};


//////////////////////////////////////////////////////////////////////////////
//Base Class for all FileFormat classes
//////////////////////////////////////////////////////////////////////////////

/**
.Class.FileFormat:
..cat:Input/Output
..summary:Object that stores a file format.
..signature:FileFormat<File, Data [, Format]>
..see:Tag.File Format
*/

template <
	typename TFile, 
	typename TData = String<typename Value<TFile>::Type>,
	typename TFormat = void,
	typename TMeta = Metadata<typename Value<TFile>::Type> >
struct FileFormat:
	public FileFormat<TFile, TData, void>
{
public:
	typedef String<typename Value<TFile>::Type> TString;
	typedef typename Size<TData>::Type TSize;
//	typedef typename Meta<TFile, TFormat>::Type TMeta;

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
	readID_(TFile & file, TString & id) const
	{
SEQAN_CHECKPOINT
		readID(file, id, TFormat());
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

	virtual void
	write_(TFile & file, TData & data, TString & id) const
	{
SEQAN_CHECKPOINT
		write(file, data, id, TFormat());
	}
	virtual void
	write_(TFile & file, TData & data, TString & id, TMeta & meta) const
	{
SEQAN_CHECKPOINT
		write(file, data, id, meta, TFormat());
	}
};

//____________________________________________________________________________

//base class for all file format classes 

template <typename TFile, typename TData, typename TMeta>
struct FileFormat<TFile, TData, void, TMeta>
{
public:
	typedef String<typename Value<TFile>::Type> TString;
	typedef typename Size<TData>::Type TSize;

	FileFormat() {}
	FileFormat(FileFormat const &) {}
	~FileFormat() {}
	FileFormat const & operator =(FileFormat const &) {}

	virtual void *
	formatID_() const = 0;

	virtual void
	read_(TFile & file, TData & data) const = 0;
	virtual void
	read_(TFile & file, TData & data, TSize limit) const = 0;

	virtual void
	readID_(TFile & file, TString & id) const = 0;

	virtual void
	readMeta_(TFile & file, TMeta & meta) const = 0;

	virtual void
	goNext_(TFile & file) const = 0;

	virtual void
	write_(TFile & file, TData & data, TString & id) const = 0;
	virtual void
	write_(TFile & file, TData & data, TString & id, TMeta & meta) const = 0;

};

//////////////////////////////////////////////////////////////////////////////
// Wrapper for functions to virtuals
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData, typename TFormat, typename TMeta>
inline void *
formatID(FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
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
template <typename TFile, typename TData, typename TFormat, typename TMeta>
inline void
read(TFile & file,
	 TData & data,
	 FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.read_(file, data);
}

template <typename TFile, typename TData, typename TFormat, typename TMeta, typename TSize>
inline void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.read_(file, data, limit);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.readID:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TFormat, typename TMeta, typename TString>
inline void
readID(TFile & file,
	   TString & id,
	   FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	String<typename Value<TFile>::Type> str;
	file_format.readID_(file, str);
	assign(id, str);
}

template <typename TFile, typename TData, typename TFormat, typename TMeta>
inline void
readID(TFile & file,
	   String<typename Value<TFile>::Type> & id,
	   FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.readID_(file, id);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.readMeta:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TFormat, typename TMeta>
inline void
readMeta(TFile & file,
		 TMeta & meta,
		 FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.readMeta_(file, meta);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TFormat, typename TMeta>
inline void
goNext(TFile & file,
	   FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	file_format.goNext_(file);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.write:
..cat:Input/Output
*/

template <typename TFile, typename TData, typename TFormat, typename TMeta>
inline void
write(TFile & file,
	  TData & data,
	  FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	typedef String<typename Value<TFile>::Type> TString_2;
	TString_2 id_2;
	file_format.write_(file, data, id_2);
}
template <typename TFile, typename TData, typename TFormat, typename TMeta, typename TString>
inline void
write(TFile & file,
	  TData & data,
	  TString & id,
	  FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	typedef String<typename Value<TFile>::Type> TString_2;
	TString_2 id_2(id);
	file_format.write_(file, data, id_2);
}
template <typename TFile, typename TData, typename TFormat, typename TMeta, typename TString>
inline void
write(TFile & file,
	  TData & data,
	  TString & id,
	  TMeta & meta,
	  FileFormat<TFile, TData, TFormat, TMeta> const & file_format)
{
SEQAN_CHECKPOINT
	typedef String<typename Value<TFile>::Type> TString_2;
	TString_2 id_2(id);
	file_format.write_(file, data, id_2, meta);
}




//////////////////////////////////////////////////////////////////////////////
// Comparison of two FileFormat objects
//////////////////////////////////////////////////////////////////////////////

template <typename TFileLeft, typename TDataLeft, typename TMetaLeft, typename TFormatLeft, typename TFileRight, typename TDataRight, typename TFormatRight, typename TMetaRight>
inline bool
operator == (FileFormat<TFileLeft, TDataLeft, TFormatLeft, TMetaLeft> const & left, 
			 FileFormat<TFileRight, TDataRight, TFormatRight, TMetaRight> const & right)
{
SEQAN_CHECKPOINT
	return formatID(left) == formatID(right);
}

template <typename TFile, typename TData, typename TFormat, typename TMeta, typename TFormat2>
inline bool
operator == (FileFormat<TFile, TData, TFormat, TMeta> const & left, 
			 Tag<TFormat2> const)
{
SEQAN_CHECKPOINT
	return formatID(left) == _ClassIdentifier<Tag<TFormat2> const>::getID();
}

template <typename TFile, typename TData, typename TFormat, typename TMeta, typename TFormat2>
inline bool
operator == (Tag<TFormat2> const,
			 FileFormat<TFile, TData, TFormat, TMeta> const & right)
{
SEQAN_CHECKPOINT
	return _ClassIdentifier<Tag<TFormat2> const>::getID() == formatID(right);
}

//____________________________________________________________________________

template <typename TFileLeft, typename TDataLeft, typename TMetaLeft, typename TFormatLeft, typename TFileRight, typename TDataRight, typename TFormatRight, typename TMetaRight>
inline bool
operator != (FileFormat<TFileLeft, TDataLeft, TFormatLeft, TMetaLeft> const & left, 
			 FileFormat<TFileRight, TDataRight, TFormatRight, TMetaRight> const & right)
{
SEQAN_CHECKPOINT
	return formatID(left) != formatID(right);
}

template <typename TFile, typename TData, typename TFormat, typename TMeta, typename TFormat2>
inline bool
operator != (FileFormat<TFile, TData, TFormat, TMeta> const & left, 
			 Tag<TFormat2> const)
{
SEQAN_CHECKPOINT
	return formatID(left) != _ClassIdentifier<Tag<TFormat2> const>::getID();
}

template <typename TFile, typename TData, typename TFormat, typename TMeta, typename TFormat2>
inline bool
operator != (Tag<TFormat2> const,
			 FileFormat<TFile, TData, TFormat, TMeta> const & right)
{
SEQAN_CHECKPOINT
	return _ClassIdentifier<Tag<TFormat2> const>::getID() != formatID(right);
}

//////////////////////////////////////////////////////////////////////////////
// allgemeine Funktionen fuer Streams
//////////////////////////////////////////////////////////////////////////////
//TODO??? Das muss in eine extra Datei

/**
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

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
