#ifndef SEQAN_HEADER_FILE_EMBL_H
#define SEQAN_HEADER_FILE_EMBL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - EMBL
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.value.EMBL:EMBL format for sequences.
*/
struct TagEmbl_;
typedef Tag<TagEmbl_> const Embl;



//////////////////////////////////////////////////////////////////////////////
// FileReader Iterator
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TFile2, typename TSpec>
inline void
goBegin(Iter<TFile, FileReader<Embl, TFile2, TSpec> > & it)
{
SEQAN_CHECKPOINT
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

	it.data_file_pos = _streamTellG(host(it)) - 1;
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
		++it.data_file_pos;

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
				++it.data_file_pos;
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
