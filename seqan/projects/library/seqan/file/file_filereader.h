#ifndef SEQAN_HEADER_FILE_FILEREADER_H
#define SEQAN_HEADER_FILE_FILEREADER_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// FileReader String
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
class String<TValue, FileReader<TFormat, TFile, TSpec> >
{
public:
	enum
	{
		BLOCK_SIZE = 0x1000
	};

	typedef typename Position<TFile>::Type TFilePosition;

	typedef typename Size<TFile>::Type TFileSize;
	typedef String<TFileSize> TABL;

	typedef typename Position<TABL>::Type TABLPosition;

	typedef String<TValue> TBuf;

	Holder<TFile, Tristate2> data_file;
	TFilePosition data_file_begin;		//file pointer to begin of data in file
	TABL data_abl;						//accumulated block lengths
	TABLPosition data_active_block;		//number of active block
	TFileSize data_active_block_begin;	//begin position of active block
	TFileSize data_active_block_end;	//end position of active block
	TBuf data_buf;						//data of active block
	bool data_scanned;					//true if the complete string was scanned


	//TODO
	//String()...

	String(TFile & fl_)
		: data_scanned(false)
	{
		setValue(data_file, fl_);
		_FileReaderString_construct(*this);
	}
	~String()
	{
		if (!dependent(data_file))
		{
			_streamClose(value(data_file));
			clear(data_file);
		}
	}
};
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
inline TFile &
_dataFile(String<TValue, FileReader<TFormat, TFile, TSpec> > & me)
{
	return value(me.data_file);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TPosition>
inline void
_FileReaderString_loadblock(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
							TPosition blocknum)
{
	typedef String<TValue, FileReader<TFormat, TFile, TSpec> > TString;

	if (blocknum > length(me.data_abl))
	{
		for (TPosition bp = length(me.data_abl); !me.data_scanned && (bp <= blocknum); ++bp)
		{
			_FileReaderString_loadblock(me, bp);
		}
	}
	else
	{
		typedef Iter<TFile, FileReader<TFormat> > TFileReaderIt;
		_streamSeekG(_dataFile(me), me.data_file_begin + blocknum * TString::BLOCK_SIZE);
		TPosition end_filepos = me.data_file_begin + (blocknum + 1) * TString::BLOCK_SIZE;
		TFileReaderIt fit(_dataFile(me), false);

		clear(me.data_buf);

		for (; !atEnd(fit) && (fit.data_file_pos < end_filepos); goNext(fit))
		{
			appendValue(me.data_buf, value(fit));
		}

		if (blocknum == length(me.data_abl))
		{
			if (blocknum > 0)
			{
				appendValue(me.data_abl, me.data_abl[blocknum-1] + length(me.data_buf));
			}
			else
			{
				appendValue(me.data_abl, length(me.data_buf));
			}

			if (atEnd(fit))
			{
				me.data_scanned = true;
			}
		}

		me.data_active_block = blocknum;
		me.data_active_block_begin = (blocknum) ? me.data_abl[blocknum - 1] : 0;
		me.data_active_block_end = me.data_abl[blocknum];
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TPosition>
inline unsigned int
_FileReaderString_findblock(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
							TPosition pos)
{
	typedef typename Size<TFile>::Type TFileSize;

	while (!me.data_scanned && (me.data_abl[length(me.data_abl) - 1] <= pos))
	{
		_FileReaderString_loadblock(me, length(me.data_abl));
	}

	SEQAN_ASSERT2(me.data_abl[length(me.data_abl) - 1] > pos, "range error")

	return ::std::lower_bound(begin(me.data_abl, Standard()), end(me.data_abl, Standard()), (TFileSize) pos) - begin(me.data_abl, Standard());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
inline void
_FileReaderString_construct(String<TValue, FileReader<TFormat, TFile, TSpec> > & me)
{
	//find begin of data in file
	typedef Iter<TFile, FileReader<TFormat> > TFileReaderIt;
	TFileReaderIt fit(_dataFile(me));
	me.data_file_begin = fit.data_file_pos;

	_FileReaderString_loadblock(me, 0);
}

//////////////////////////////////////////////////////////////////////////////

//tests whether block_number will exist when the file was scanned completely

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TUint>
inline bool
_FileReaderString_isValidBlock(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
							   TUint block_number)
{
	typedef typename Size<TFile>::Type TFileSize;

	while (!me.data_scanned && (length(me.data_abl) <= block_number))
	{
		_FileReaderString_loadblock(me, length(me.data_abl));
	}
	return (length(me.data_abl) > block_number);

}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec, typename TPos>
inline TValue
value(String<TValue, FileReader<TFormat, TFile, TSpec> > & me,
	  TPos pos)
{
	if ((me.data_active_block_begin > pos) || (me.data_active_block_end <= pos))
	{//change block
		_FileReaderString_loadblock(me, _FileReaderString_findblock(me, pos));
	}
	return me.data_buf[pos - me.data_active_block_begin];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat, typename TFile, typename TSpec>
inline typename Size< String<TValue, FileReader<TFormat, TFile, TSpec> > >::Type
length(String<TValue, FileReader<TFormat, TFile, TSpec> > & me)
{
	if (!me.data_scanned)
	{//scan the whole sequence
		typedef typename Position<TFile>::Type TPosition;
		_FileReaderString_loadblock(me, supremumValue<TPosition>());
	}

	return me.data_abl[length(me.data_abl) - 1];
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Iterator for FileReader String
// (note: do not confuse with FileReader Iterator, see file_filereaderiterator.h)
//////////////////////////////////////////////////////////////////////////////

struct FileReaderIterator;

template <typename TContainer>
class Iter<TContainer, FileReaderIterator>
{
public:
	typedef typename TContainer::TABLPosition TABLPosition;

	typedef typename TContainer::TBuf TBuf;
	typedef typename Position<TBuf>::Type TBufPosition;
	typedef typename Size<TBuf>::Type TBufSize;

	struct TContainer * data_container;
	struct TABLPosition data_abl_pos;	//number of block
	struct TBufPosition data_buf_pos;	//number of char in block
	struct TBufSize data_buf_len;		//length of block
	bool data_atEnd;					//true if iterator is atEnd
};

//////////////////////////////////////////////////////////////////////////////

template <typename TContainer>
inline void
goNext(Iter<TContainer, FileReaderIterator> & it)
{
	++it.data_buf_pos;
	if (it.data_buf_pos >= it.data_buf_len)
	{
		if (!it.data_atEnd)
		{
			it.data_buf_pos = 0;
			++it.data_abl_pos;

			TContainer & cont = *(it.data_container);
			it.data_atEnd = _FileReaderString_isValidBlock(cont, it.data_abl_pos);
			if (it.data_atEnd)
			{//at end
				it.data_buf_len = 0
			}
			else
			{//not at end
				if (cont.data_active_block != it.data_abl_pos)
				{
					_FileReaderString_loadblock(cont, it.data_abl_pos);
				}
				it.data_buf_len = length(it.data_buf);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
