#ifndef SEQAN_HEADER_FILE_FASTA_H
#define SEQAN_HEADER_FILE_FASTA_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Fasta
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.value.Fasta:
	FASTA file format for sequences.
*/
struct TagFasta_;
typedef Tag<TagFasta_> const Fasta;

//////////////////////////////////////////////////////////////////////////////

// File Reader Iterator

template <typename TFormat>
struct FileReader;


template <typename TFile>
class Iter<TFile, FileReader<Fasta> >
{
public:
	typedef typename Value<TFile>::Type TValue;

	TFile * data_host;
	TValue data_char;
//	TFilePosition data_begin_pos;

	Iter(TFile & file_):
		data_host(& file_),
		data_char(0)
	{
		if (_streamEOF(data_host)) return;
		_stream_skipLine(data_host, data_char);
	}
	Iter(Iter const & other_):
		data_host(other_.data_host),
		data_char(other_.data_char)
	{
	}
	~Iter() 
	{
	}

	Iter const &
	operator = (Iter const & other_)
	{
		data_host = other_.data_host;
		data_char = other_.data_char;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
struct Value< Iter<TFile, FileReader<Fasta> > >:
	Value<TFile>
{
};

template <typename TFile>
struct GetValue< Iter<TFile, FileReader<Fasta> > >
{
	typedef typename Value< Iter<TFile, FileReader<Fasta> > >::Type Type;
};

template <typename TFile>
struct Reference< Iter<TFile, FileReader<Fasta> > >
{
	typedef typename Value< Iter<TFile, FileReader<Fasta> > >::Type & Type;
};


//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
inline typename Reference<Iter<TFile, FileReader<Fasta> > >::Type
value(Iter<TFile, FileReader<Fasta> > & it)
{
	return it.data_char;
}

template <typename TFile>
inline typename GetValue<Iter<TFile, FileReader<Fasta> > >::Type
getValue(Iter<TFile, FileReader<Fasta> > & it)
{
	return it.data_char;
}

template <typename TFile>
inline void
goNext(Iter<TFile, FileReader<Fasta> > & it)
{
	do
	{
		if (_streamEOF(it.data_file)) return;
		it.data_char = _streamGet(it.data_file);
	} while ((it.data_char == '\n') || (it.data_char == '\r'));
}

template <typename TFile>
inline bool
atEnd(Iter<TFile, FileReader<Fasta> > & it)
{
	return _streamEOF(it.data_file) || (it.data_char == '>');
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////
//count_valid: zaehlt die nicht-Zeilenumbrueche (input/output)
//count_all: zaehlt alle Zeichen incl. Zeilenumbrueche (input/output)
//returns: zuletzt gelesenes Zeichen = das erste hinter dem Zeilenumbruch bzw. eof
//the last read char is not counted!
//count_valid and count_all are not resetted but counted up
template <typename TFile, typename TSize>
inline typename Value<TFile>::Type
_fasta_scan_line(TFile & file,
				 TSize & count_valid,
				 TSize & count_all)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	TSize count = 0;

	while (true)
	{
		typename Value<TFile>::Type c = _streamGet(file);

		if (_streamEOF(file))
		{
			count_valid += count;
			count_all += count;
			return c;
		}

		if ((c == '\n') || (c == '\r'))
		{
			do
			{
				++count_all;
				c = _streamGet(file);
			} while ((c == '\n') || (c == '\r'));

			count_valid += count;
			count_all += count;
			return c;
		}

		if (c != '\r')
		{
			++count;
		}
	}
}


/////////////////////////////////////////////////////////////////////////
template <typename TFile, typename TSize>
inline void
_read_n_chars_from_file(TFile & file, TSize count)
{
SEQAN_CHECKPOINT
	for (TSize i = 0; i < count; ++i)
	{
		_streamGet(file);
	}
}


//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData, typename TSize>
void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Fasta)
{
SEQAN_CHECKPOINT

	SEQAN_ASSERT(!_streamEOF(file))

	//determine begin position
	typename Value<TFile>::Type c_first = _streamGet(file);
	SEQAN_ASSERT(!_streamEOF(file))

	typename Position<TFile>::Type begin_pos = _streamTellG(file);
	typename Size<TData>::Type count_valid = 1; //"valid" characters read (without line breaks)
	typename Size<TData>::Type count_all = 1;	//all characters read (with line breaks)

	if (c_first == '>')
	{//there is an id line: skip it
		c_first = _fasta_scan_line(file, count_valid, count_all);
	}

	if (c_first == '>') 
	{//another id line = empty entry
		clear(data);
		_streamSeekG(file, begin_pos);
		_read_n_chars_from_file(file, count_all);
		return;
	}

	begin_pos = _streamTellG(file);

	count_valid = 1;
	count_all = 1;
	typename Value<TFile>::Type c;

	//determine length
	while (true)
	{
		c = _fasta_scan_line(file, count_valid, count_all);
		if (_streamEOF(file)) 
		{//end of file: stop searching
			break;
		}
		if (c == '>')
		{//next entry found: stop seaching
			break;
		}
		if ((c != '\n') && (c != '\r'))
		{
			++count_valid; //count c
		}
		++count_all;
	}

	//reserve space
	typename Size<TData>::Type count = count_valid;
	if (count > limit)
	{
		count = limit;
	}
	resize(data, count);
	if (length(data) < count)
	{
		count = length(data);
	}

	//read sequence
	_streamSeekG(file, begin_pos);

	typename Position<TData>::Type pos = 0;
	c = c_first;
	while (true)
	{
		if ((c != '\n') && (c != '\r'))
		{
			data[pos] = c;
			++pos;
		}
		if (pos >= count) break;

		c =  _streamGet(file);
		--count_all;
	}

	//move file ptr to next entry
	_read_n_chars_from_file(file, count_all - 1);
}

//____________________________________________________________________________

template <typename TFile, typename TData>
void
read(TFile & file,
	 TData & data,
	 Fasta tag)
{
SEQAN_CHECKPOINT
	typedef typename Size<TData>::Type TSize;
	read(file, data, supremumValue<TSize>(), tag);
}


//////////////////////////////////////////////////////////////////////////////
// readID
//////////////////////////////////////////////////////////////////////////////
 
//the ID is the complete first line (without the leading '>'-sign)

template <typename TFile, typename TString>
void
readID(TFile & file,
	   TString & id,
	   Fasta)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	typename Position<TFile>::Type start_pos = _streamTellG(file);

	typename Value<TFile>::Type c = _streamGet(file);
	if (c != '>')
	{
		clear(id);
	}
	else
	{
		typename Size<TString>::Type count_valid = 0;
		typename Size<TString>::Type count_all = 0;
		_fasta_scan_line(file, count_valid, count_all);

		if (! count_valid)
		{
			clear(id);
		}
		else
		{
			resize(id, count_valid);
			if (length(id) < count_valid)
			{
				count_valid = length(id);
			}

			_streamSeekG(file, start_pos);
			c = _streamGet(file); //pop the '>' character
			for (typename Position<TString>::Type pos = 0; count_valid; --count_valid)
			{
				id[pos] = _streamGet(file);
				++pos;
			}
		}
	}
	_streamSeekG(file, start_pos);
}

//////////////////////////////////////////////////////////////////////////////
// readMeta
//////////////////////////////////////////////////////////////////////////////

//Fasta file records have no meta data

template <typename TFile, typename TMeta>
void
readMeta(TFile & file,
		 TMeta & meta,
		 Fasta)
{
SEQAN_CHECKPOINT
//	clear(meta);
}


//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
void
goNext(TFile & file,
	   Fasta)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	bool found_data = false;
	while (true)
	{
		typename Value<TFile>::Type c = _streamGet(file);

		if (_streamEOF(file)) return;

		if (c == '\n' || c == '\r')
		{
			do {
				c = _streamGet(file);
				if (_streamEOF(file)) return;
			} while (c == '\n' || c == '\r');

			if (c != '>')
			{
				found_data = true;
			}
			else if (found_data)
			{
				_streamUnget(file);
				return;
			}
		}
	}
}



//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////


template <typename TFile, typename TString, typename TData>
void
_write_impl(TFile & file,
			TData & data,
			TString & id,
			Fasta)
{
SEQAN_CHECKPOINT
	_streamPut(file, '>');
	_streamWrite(file, id);
	_streamPut(file, '\n');

	typename Iterator<TData, Standard>::Type it = begin(data, Standard());
	typename Iterator<TData, Standard>::Type it_end = end(data, Standard());

	int i = 0;

	for (; it < it_end; ++it)
	{
		if (i == 60)
		{
			_streamPut(file, '\n');
			i = 0;
		}
		++i;

		_streamPut(file, *it);
	}
	_streamPut(file, '\n');
}

//____________________________________________________________________________

template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData & data,
	  Fasta)
{
SEQAN_CHECKPOINT
	_write_impl(file, data, "", Fasta());
}

//____________________________________________________________________________

template <typename TFile, typename TString, typename TData>
void
write(TFile & file,
	  TData & data,
	  TString & id,
	  Fasta)
{
SEQAN_CHECKPOINT
	_write_impl(file, data, id, Fasta());
}


//VisualC++ const array bug workaround
template <typename TFile, typename TString, typename TDataValue>
void
write(TFile & file,
	  TDataValue * data,
	  TString & id,
	  Fasta)
{
SEQAN_CHECKPOINT
	_write_impl(file, data, id, Fasta());

}

//____________________________________________________________________________

template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData & data,
	  TString & id,
	  TMeta &,
	  Fasta)
{
SEQAN_CHECKPOINT
	_write_impl(file, data, id, Fasta());
}


//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
