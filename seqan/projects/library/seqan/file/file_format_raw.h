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

#ifndef SEQAN_HEADER_FILE_RAW_H
#define SEQAN_HEADER_FILE_RAW_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Raw
//////////////////////////////////////////////////////////////////////////////
/**
.Tag.File Format.tag.Raw:
	The file contains data in a raw format.
..remark:It is supposed that the file contains one single piece of data, 
that is the file cannot store multiple records.
..include:seqan/file.h
*/

struct TagRaw_;
typedef Tag<TagRaw_> const Raw;



//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData, typename TTag>
struct _Read_Raw;

//____________________________________________________________________________

template <typename TFile, typename TData>
struct _Read_Raw<TFile, TData, True>
{
	template <typename TSize>
	inline static void
	read_(TFile & file,
		TData & data,
		TSize _limit)
	{
SEQAN_CHECKPOINT
		SEQAN_ASSERT(!_streamEOF(file))

		typename Size<TData>::Type limit = _limit;

		//determine length
		typename Position<TFile>::Type begin_pos = _streamTellG(file);
		typename Size<TData>::Type count;
		typename Position<TData>::Type pos;

		for (count = 0; !_streamEOF(file) && count != limit; ++count)
			_streamGet(file);

		//reserve space
		resize(data, count);

		if (!count) return;

		if (length(data) < count)
			count = length(data);

		//read sequence
		_streamSeekG(file, begin_pos);

		for (pos = 0; pos < count; ++pos)
			assignValue(data, pos, _streamGet(file));
	}

//____________________________________________________________________________

	inline static void
	read_(TFile & file,	
		TData & data)
	{
SEQAN_CHECKPOINT
		typedef typename Size<TData>::Type TSize;
		read(file, data, supremumValue<TSize>());
	}

};

//____________________________________________________________________________

template <typename TFile, typename TData>
struct _Read_Raw<TFile, TData, False>
{
	static void
	read_(TFile & file,
		TData & data)
	{
SEQAN_CHECKPOINT

		clear(data);
		if (!_streamEOF(file))
		{
SEQAN_CHECKPOINT
			_ChunkCollector<TData> chunk_collector(data);
			assign(chunk_collector, file);
			append(data, chunk_collector);
		}
	}

//____________________________________________________________________________

	template <typename TSize>
	static void
	read_(TFile & file,
		TData & data,
		TSize limit)
	{
SEQAN_CHECKPOINT

		clear(data);
		if (!_streamEOF(file))
		{
SEQAN_CHECKPOINT
			_ChunkCollector<TData> chunk_collector(data);
			assign(chunk_collector, file, limit);
			append(data, chunk_collector, limit);
		}
	}
};

//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData>
void
read(TFile & file,
	 TData & data,
	 Raw)
{
SEQAN_CHECKPOINT
	_Read_Raw<TFile, TData, typename _IsTellSeekStream<TFile>::Type>::read_(file, data);
}

//____________________________________________________________________________

template <typename TFile, typename TData, typename TSize>
void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Raw)
{
SEQAN_CHECKPOINT
	_Read_Raw<TFile, TData, typename _IsTellSeekStream<TFile>::Type>::read_(file, data, limit);
}


//////////////////////////////////////////////////////////////////////////////
// readID
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TString>
void
readID(TFile & /*file*/,
	   TString & id,
	   Raw)
{
SEQAN_CHECKPOINT
	clear(id);
}

//////////////////////////////////////////////////////////////////////////////
// readMeta
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TMeta>
void
readMeta(TFile & /*file*/,
		 TMeta & meta,
		 Raw)
{
SEQAN_CHECKPOINT
	clear(meta);
}


//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
void
goNext(TFile & file,
	   Raw)
{
SEQAN_CHECKPOINT
  (void) file;  // When compiled without assertions.
	SEQAN_ASSERT_TRUE(!_streamEOF(file));

//??? TODO: set file to eof
}



//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData>
void
write(TFile & file,
	  TData const & data,
	  Raw)
{
SEQAN_CHECKPOINT
	_streamWrite(file, data);
}

//____________________________________________________________________________

template <typename TFile, typename TData, typename TString>
void
write(TFile & file,
	  TData const & data,
	  TString const &,
	  Raw)
{
SEQAN_CHECKPOINT
	_streamWrite(file, data);
}

//____________________________________________________________________________

template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData const & data,
	  TString const &,
	  TMeta const &,
	  Raw)
{
SEQAN_CHECKPOINT
	_streamWrite(file, data);
}

//////////////////////////////////////////////////////////////////////////////
// default functions
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TData>
void
read(TFile & file,
	 TData & data)
{
	read(file, data, Raw());
}

template <typename TFile, typename TData, typename TSize>
void
read(TFile & file,
	 TData & data,
	 TSize limit)
{
SEQAN_CHECKPOINT
	read(file, data, limit, Raw());
}

//____________________________________________________________________________

template <typename TFile, typename TData>
void
write(TFile & file,
	  TData & data)
{
SEQAN_CHECKPOINT
	write(file, data, "", Raw());
}
template <typename TFile, typename TData>
void
write(TFile & file,
	  TData const & data)
{
SEQAN_CHECKPOINT
	write(file, data, "", Raw());
}

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
