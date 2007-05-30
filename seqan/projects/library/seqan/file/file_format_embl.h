#ifndef SEQAN_HEADER_FILE_EMBL_H
#define SEQAN_HEADER_FILE_EMBL_H

#include <map>

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

// typedefs 


//////////////////////////////////////////////////////////////////////////////
// Internal functions
//////////////////////////////////////////////////////////////////////////////

// reads one line until \n. similar to fasta_scan_line

/**
.Internal._embl_scan_line:
..summary: Reads one line of a $file$ and increases $count$.
...remarks: Do not moves back to the beginning of the line.
..signature:_embl_scan_line(file, count)
..param.file:The input file.
..param.count:The counter which is increased by the number of signs in the given line. 
..remarks.text:Same as $_genbank_scan_line$.
..see:Internal._genbank_scan_line.
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/

template <typename TFile, typename TSize>
void
_embl_scan_line(	TFile & file,
			TSize & count)
{
SEQAN_CHECKPOINT
//	SEQAN_ASSERT(!_streamEOF(file))
	if (_streamEOF(file)) return;

	int pos = 0;
	while (true)
	{
		typename Value<TFile>::Type c = _streamGet(file);

		if (_streamEOF(file)) return;
		if (c == '\n') return;
		if (c != '\r')
		{
			++count;
		}
	}	
}

/**
.Internal._embl_scan_tag:
..summary: Reads characters from a stream until the next '\n', 't' or " " and increases 
$count$ (Internal use only).
...remarks: Do not moves back to the beginning of the line.
..signature:_embl_scan_tag(file, count)
..param.file: The input file.
..param.count: The counter which is increased by the number of characters read. 
..see:Tag.File Format.value.EMBL
*/

template <typename TFile, typename TSize>
void
_embl_scan_tag(	TFile & file,
			TSize & count)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))
	
	typename Value<TFile>::Type c = _streamGet(file);
	while ((c!='\n')&&(c!=' ')&&(!_streamEOF(file))&&(c!='\t'))
	{
		if (c!='\r') ++count;
		c = _streamGet(file);
	}
}


//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////

/**
.Function.read:
..summary: Reads the sequence data of a $EMBL$ file and stores it in $data$. 
..signature:read(file, data [, limit], Embl)
..param.file:The input file.
..param.data:The data object in which the sequence data is parsed.
..param.meta:The meta data container.
..param.Embl: Embl Tag
...Type:Tag.File Format.EMBL
...remarks.text: Can be guessed by $guessFileFormat$.
..see:Function.readID
..see:Function.readMeta
..see:Function.write
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/

// call of 'read' with 'limit' // bleibt noch hängen, wenn eof??
template <typename TFile, typename TData, typename TSize>
void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Embl)
{
SEQAN_CHECKPOINT

	SEQAN_ASSERT(!_streamEOF(file))

	typename Position<TFile>::Type start_pos = _streamTellG(file);	// default data start position
	typename Position<TFile>::Type end_pos;							// data end position
	typename Size<TData>::Type count = 0;							// default TData size
	typename Value<TFile>::Type c;									// holds one character from stream
	String<char,Array<2> > origin_tag = "  ";						// holds six characters
	int origin_found = 0;

	// find ORIGIN tag
	while ((origin_found==0)&&(!_streamEOF(file)))
	{
		
		// reads the first six characters from line
		for (int i=0; i<length(origin_tag);++i)
		{
			origin_tag[i] = _streamGet(file);
		}

		// if ORIGIN tag found, set flag = 1
		if (origin_tag=="SQ")
		{
			// set flag = true
			origin_found=1;
			start_pos = _streamTellG(file);
		}
		else // go 'til end of line
		{
			c = _streamGet(file);
			while((c!='\n')&&(!_streamEOF(file)))
			{
				c = _streamGet(file);
			}
		} 
	} // while
	
	clear(data);
	
	// reads over 'ORIGIN'
	_embl_scan_line(file, count); 
	start_pos = _streamTellG(file);
	
	count = 0;
	
	// determine size of sequence
	int exit_flag = 0;
	while ((exit_flag==0)&&(!_streamEOF(file)))
	{
		c = _streamGet(file);

		int c_int = c;
		if (((c_int>64)&&(c_int<91))||((c_int>96)&&(c_int<123)))
		{
			++count;
		}

		if (c=='/') 
		{	
			c = _streamGet(file);
			if (c=='/') 
			exit_flag=1;
			end_pos = _streamTellG(file);
		}
	}

	// if limit set
	if (count > limit)
	{
		count = limit;
	}
	
	//reserve space
	resize(data, count);
	if (length(data) < count)
	{
		count = length(data);
	}
	
	//read sequence, go to starting position
	_streamSeekG(file, start_pos);

	// read the sequence into data
	typename Position<TData>::Type pos = 0;
	while((_streamTellG(file)!=end_pos)&&(!_streamEOF(file)))
	{
		c = _streamGet(file);
		int c_int = c;
		if (((c_int>64)&&(c_int<91))||((c_int>96)&&(c_int<123)))
		{
			data[pos] = c;
			++pos;
		}
	}
	//read sequence, go to starting position
	_streamSeekG(file, end_pos);
	
}

// call of 'read' without 'limit'
template <typename TFile, typename TData>
void
read(TFile & file,
	 TData & data,
	 Embl tag)
{
SEQAN_CHECKPOINT
	typedef typename Size<TData>::Type TSize;
	read(file, data, supremumValue<TSize>(), tag);
}


//////////////////////////////////////////////////////////////////////////////
// readID
//////////////////////////////////////////////////////////////////////////////

/**
.Function.readID:
..summary: Reads the id of a $EMBL$ sequence record and stores it in $id$. 
..signature:readID(file, id, Embl)
..param.file:The input file.
..param.id:The object in which the sequence $id$ is stored.
..param.Embl: Embl Tag
...Type:Tag.File Format.EMBL
...remarks.text: Can be guessed by $guessFileFormat$. 
..remarks.text: The ID in an EMBL sequence record is the first entry in the line beginning with 'AC', 
not with 'ID'.
..see:Function.read
..see:Function.readMeta
..see:Function.write
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/

template <typename TFile, typename TString>
void
readID(TFile & file,
	   TString & id,
	   Embl)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	typename Position<TFile>::Type begin_pos = _streamTellG(file);	// contains start position of gb entry
	typename Position<TFile>::Type id_pos = _streamTellG(file);	    // contains start position of id
	typename Size<TString>::Type count = 0;							// length of id
	typename Value<TFile>::Type c;									// holds one character from stream
	String<char,Array<2> > locus_tag = "  ";					    // holds the AC tag
	
	clear(id); 

	bool exit_flag = false;
	while (!exit_flag)
	{
		//  read first two chars from file
		for (int i =0 ;i<length(locus_tag);++i)
		{
			locus_tag[i] = _streamGet(file);
		}

	//  check if its ==ACESSION
	if (locus_tag=="AC")
	{
		// -1 before id, 0 in id, 1 id ends
		int escape_flag = -1;
		while(escape_flag<1)
		{
			// read one character from file
			c = _streamGet(file);

			// id found if char stream is != tab or space
			if ((escape_flag<0)&&((c!=' ')&&(c!='\t')))
			{
				escape_flag++;
				_streamSeek2G(file, -1);
				id_pos = _streamTellG(file); // set begin pos of id
				_streamSeek2G(file, +1);
				count = 0;
			}

			// increase count if id found
			if (escape_flag==0)
			{
				if (c!='\n')
				{
					count++;
				}
				else 
				{
					escape_flag=1;
				} // if 
			} // if

		} //while 
		
		// allocate memory for id
		resize(id,count);

		// read id
		_streamSeekG(file, id_pos);
		for (typename Position<TString>::Type pos = 0; pos<count; pos++)
		{
			id[pos] = _streamGet(file);
		}
		exit_flag = true;

	} // if not AC empty id and exit
	else 
	{
		if (!_streamEOF(file))
		{
			_embl_scan_line(file, count);
		}
		else
		{
			exit_flag = true;
		}
	}
	
	} //while
	// move stream to begin of file
	_streamSeekG(file, begin_pos);
} // readID


//////////////////////////////////////////////////////////////////////////////
// readMeta
//////////////////////////////////////////////////////////////////////////////

/**
.Function.readMeta:
..summary: Reads the id of a $EMBL$ file and stores it in $meta$. 
..signature:readMeta(file, data [, meta], Embl)
..param.file:The input file.
..param.id:The object in which the sequence $id$ is stored.
..param.meta:The meta data container.
..param.Embl: Embl Tag
...Type:Tag.File Format.EMBL
...remarks.text: Can be guessed by $guessFileFormat$. 
..remarks.text: Metadata is everything but the $id$ and the $sequence$.
..see:Function.read
..see:Function.readID
..see:Function.write
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/

template <typename TFile, typename TMeta>
void
readMeta(TFile & file,
		 TMeta & meta,
		 Embl)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	typename Position<TFile>::Type start_pos	= _streamTellG(file); // store starting position of stream
	typename Position<TFile>::Type act_pos		= _streamTellG(file); // temp: position of stream
	typename Value<TFile>::Type c;

	typename Size<TMeta>::Type count = 0;
	typename Position<String<char, Alloc<void> > >::Type pos = 0;

	// temp for the line tag
	String<char, Array<2> >		line_tag	= "  ";
	String<char, Alloc<void> >	tag			= "  ";
	String<char, Alloc<void> >	old_tag		= "  ";
	String<char, Alloc<void> >	temp_tag	= "  ";
	String<char, Alloc<void> >	feat_tag	= "  ";
	String<char> content;
	
	int line_type = -1;
	bool first_appearance = true;
	bool empty_feature = false;
	int feature_gap_length = 0;
	bool got_feat_name = false;
	int fh_count = 0;

	clear(meta);


	while(!_streamEOF(file))
	{
		empty_feature = false;

		// get line tag
		line_tag[0] = _streamGet(file);
		line_tag[1] = _streamGet(file);
		
		// check for empty feature
		c = _streamGet(file);
		if (c=='\n') empty_feature = true;
		_streamSeek2G(file, -1);

		// exit if sequence reached
		if (line_tag == "SQ") 
		{
			_streamSeek2G(file,-2);
			typename Position<TFile>::Type sequ_start = _streamTellG(file);
			// add meta tag for sequence
			
			_streamSeek2G(file,+2);
			count = 0;
			do
			{
				c = _streamGet(file);
				if (c=='\n') 
				{ 
					count = 1;
				}
			}
			while (c==' ');
			_streamSeek2G(file,-1);
			
			if (count==0)
			{
				act_pos = _streamTellG(file);
				_embl_scan_line(file, count);
				resize (content, count+1);
				_streamSeekG(file, act_pos);
				for (pos=0;pos<count;pos++)
				{
					content[pos] = _streamGet(file); 
				}
				content[length(content)-1] = '\n';

				meta.insert(MetadataEntry(line_tag, content));
			}
			else
			{
				resize(content,0);
				meta.insert(MetadataEntry(line_tag, content));
				clear(content);
				resize(content,0);
			}

			_streamSeekG(file, sequ_start);

			return;
		}

		if (line_tag=="FH")
		{
			fh_count++;
			if (fh_count==2)
			{
				meta.insert( MetadataEntry(line_tag, content));
			}
		
		}
		
		if ((line_tag == "FT")) // "FT"-block
		{
			// if first time, number of spaces before a feature key appears
			if (first_appearance)
			{
				act_pos = _streamTellG(file);
				
				// examine how many spaces are between 'FT' and the feature key 
				do
				{
					c = _streamGet(file);
					feature_gap_length++;
				}
				while ((c==' ')||(c=='\t'));
				_streamSeekG(file, act_pos);

				first_appearance = false;
				feature_gap_length--;
			}

			// skip spaces in FT-block
			for (int i=0;i<feature_gap_length;i++)
			{
				c = _streamGet(file);
			}

			c = _streamGet(file);
			_streamSeek2G(file,-1);
			
			if (c==' ') // add to old feature
			{
				// no new tag is read
				resize(old_tag,length(tag));
				old_tag = tag;

				int old_content_length = length(content);
				
				// skip leading spaces
				do
				{
					c = _streamGet(file);
					if (c=='\n') empty_feature=true;
				}
				while (c==' ');
				_streamSeek2G(file,-1);

				if (!empty_feature)
				{
					act_pos = _streamTellG(file);
					count = 0;
					_embl_scan_line(file, count);
					resize (content, count+old_content_length+1);

					_streamSeekG(file, act_pos);
					for (pos=old_content_length;pos<length(content)-1;pos++)
					{
						content[pos] = _streamGet(file); 
					}
					content[length(content)-1] = '\n';
				}

			}
			else // new feature
			{
				// add completed meta information onto multimap
				if (fh_count!=3) 
				{
					meta.insert (MetadataEntry(tag,content));
				}
				else // bugfix
				{
					fh_count++;
				}

				// examnie length of new tag
				act_pos = _streamTellG(file);
				
				count = 0;
				_embl_scan_tag(file, count);
				resize(tag,count+1);
				clear(tag);

				_streamSeekG(file, act_pos);

				// read new tag
				tag += '|';
				for (pos=1;pos<count+1;pos++)
				{
					c = _streamGet(file);
					tag += c;
				}

				resize(old_tag, length(tag));
				old_tag = tag;

				// read over spaces between feature key and content
				do
				{
					c = _streamGet(file);
					if (c=='\n') empty_feature=true;
				}
				while (c==' ');
				_streamSeek2G(file,-1);
				
				if (!empty_feature)
				{
					// examine length of feature content
					count = 0;
					act_pos = _streamTellG(file);
					_embl_scan_line(file,count);

					resize(content,count+1);
					clear(content);

					_streamSeekG(file, act_pos);
					for (pos=0;pos<count;pos++)
					{
						c = _streamGet(file); 
						content += c;
					}
					content[length(content)-1] = '\n';
				}
			}
		}
		else // normal lines
		{
			resize(temp_tag, length(tag));
			temp_tag = tag;

			resize(tag, 2);
			tag = line_tag;
			
			if  (fh_count==2) 
			{
				//std::cout << fh_count << '\n';
				fh_count++;
			}
		/*	std::cout << fh_count << '\n'; */

			
			if (tag!="XX") // add to old feature
			{
				int old_content_length = length(content);
				
				if (got_feat_name==false)
				{
					resize(old_tag,length(tag));
					old_tag = tag;
					got_feat_name = true;

					// skip leading spaces
					do
					{
						c = _streamGet(file);
						if (c=='\n') empty_feature=true;
					}
					while (c==' ');

					_streamSeek2G(file,-1);

				}
				else
				{
					_streamSeek2G(file, -2);
				}

				if (!empty_feature)
				{
					count = 0;
					act_pos = _streamTellG(file);
					_embl_scan_line(file,count);
				
					resize(content,(count+old_content_length+1));
				
					_streamSeekG(file, act_pos);
					for (pos=old_content_length;pos<length(content)-1;pos++)
					{
						content[pos] = _streamGet(file); 
					}
					content[length(content)-1] = '\n';
				}
				
			}
			else // add new feature
			{ 
		
				if ((temp_tag!="XX")) // don't add empty lines
				{
					 meta.insert (MetadataEntry(old_tag,content));
				}
				resize (tag, length(old_tag));
				old_tag = tag;
				got_feat_name = false;
			
				do
				{
					c = _streamGet(file);
					if (c=='\n') empty_feature=true;
				}
				while ((c==' ')||(c=='\t'));
				_streamSeek2G(file,-1);
				
				clear(content);
				if (!empty_feature)
				{
					// examine length of feature content
					count = 0;
					act_pos = _streamTellG(file);
					_embl_scan_line(file,count);

					resize(content,count+1);
					_streamSeekG(file, act_pos);
					for (pos=0;pos<count;pos++)
					{
						content[pos] = _streamGet(file); 
					}
					content[length(content)-1] = '\n';
				}			
			}	
		}

	 c = _streamGet(file); // overread linebreak
	} //while
}// readMeta


//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

/**
.Function.goNext:
..summary: Moves to the next sequence record in multi $EMBL$ file. 
..signature:readID(file, Embl)
..param.file:The input file.
..param.Embl: Embl Tag
...Type:Tag.File Format.EMBL
...remarks.text: Can be guessed by $guessFileFormat$. 
..remarks.text: $EOF$ if no record is found, 
not with 'ID'.
..see:Function.read
..see:Function.readMeta
..see:Function.write
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/

template <typename TFile>
void
goNext(TFile & file,
	   Embl)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))
		
	typename Position<TFile>::Type end_pos = _streamTellG(file);
	String<char,Array<2> > locus_tag = "  ";	
	typename Value<TFile>::Type c;

	while (true)
	{
		c = _streamGet(file);
		
		if (_streamEOF(file))
		{
			end_pos = _streamTellG(file);
			_streamSeekG(file, end_pos);
			return;
		} // file ends
		
		if (c=='I')
		{
			locus_tag[0] = c;
			locus_tag[1] = _streamGet(file);
		
			if (locus_tag == "ID")
			{
				_streamSeek2G(file, -2);
				end_pos = _streamTellG(file);
				_streamSeekG(file, end_pos);
				return;
			}
		}	
	}
	 // go to eof or next dataset
}

//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////

/**
.Internal._embl_write_impl:
..summary: Implementation of the $write$ function for $EMBL$ sequence records. 
..signature:_embl_write_impl(file, data [, id] [, meta], Embl)
..param.file:The input file.
..param.data:The sequence data.
..param.id:The sequence id.
..param.meta:The meta data container.
..param.Embl: Embl Tag.
..see:Function.write
*/

template <typename TFile, typename TString, typename TData >
void
_embl_write_impl(TFile & file,
			TData & data,
			TString & id,
			MetadataMap & meta,
			Embl)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))
	
	// constants
	int const row_size	 = 60;  // length of one row
	int const block_size = 10;  // size of one block 
	int const tag_gap	 = 3;   // # of spaces between the ID-tag and the id, has to be >= 1
	int const seq_space_count	= 5; 
	int const feat_space_count	= 21;
	int const feat_lead_count	= 5;

	// write ID
	_streamWrite(file, "ID");
	for (int i=0;i<tag_gap;i++)
	{
		_streamPut(file,' ');
	}
	MetadataMap::iterator act_entry = meta.find("ID");

	_streamWrite(file, act_entry->second);
	_streamWrite(file, "XX\n");

	// write AC
	_streamWrite(file, "AC");
	for (int i=0;i<tag_gap;i++)
	{
		_streamPut(file,' ');
	}
	_streamWrite(file, id);	
	_streamPut(file, '\n');
	_streamWrite(file, "XX\n");

	// write meta data
	MetadataMap::iterator meta_iter;
	String<char> feat_tag;

	for ( meta_iter = meta.begin(); meta_iter != meta.end(); meta_iter++ )
	{
		resize(feat_tag, length(meta_iter->first));
		clear(feat_tag);
		
		bool feat = false;
		if (meta_iter->first[0]=='|') 
		{ 
			feat = true; 
		}

		if ((feat==false)&&(meta_iter->first!="ID")&&(meta_iter->first!="AC")&&(meta_iter->first!="SQ"))
		{
			_streamWrite(file, meta_iter->first);
		
			// write spaces
			for (int i=0;i<tag_gap;i++)
			{
				_streamPut(file,' ');
			}
			_streamWrite(file, meta_iter->second);
			_streamWrite(file, "XX\n");
		}
	}
	
	typename Value<String<char> >::Type c = ' ';
	for ( meta_iter = meta.begin(); meta_iter != meta.end(); meta_iter++ )
	{
		if (meta_iter->first[0]=='|') 
		{
			_streamWrite(file, "FT");
			for (int i=0;i<tag_gap;i++)
			{
				_streamPut(file, ' ');
			}
			
			for (int j=1;j<length(meta_iter->first);j++)
			{
				_streamPut(file, meta_iter->first[j]);
			}
			
			for (int j=0;j<(feat_space_count-feat_lead_count-(length(meta_iter->first)-1));j++)
			{
				_streamPut(file, ' ');
			}

			for (int i=0;i<length(meta_iter->second);i++)
			{
				c = meta_iter->second[i];
				_streamPut(file, c); 
			
				if (c=='\n')
				{
					if (i!=(length(meta_iter->second)-1))	
					{
						_streamWrite(file, "FT");
						for (int j=0;j<feat_space_count-2;j++)
						{
							_streamPut(file, ' ');
						}
					}
				}
			}	
		}
	}
	_streamWrite(file, "XX\n");

	// write sequence data
	_streamWrite(file, "SQ");
	for (int i=0;i<tag_gap;i++)
	{
		_streamPut(file,' ');
	}
	act_entry = meta.find("SQ");
	_streamWrite(file, act_entry->second);

	typename Iterator<TData>::Type it		= begin(data);
	typename Iterator<TData>::Type it_end	= end(data);

	int j = 0; // coloumn counter
	int r = 1; // row counter

	for (; it < it_end; ++it)
	{	
		// space, if block completed
		if ((j % block_size) == 0)
		{
			if (j!=0) 
			{
				_streamPut(file, ' '); 
			}
		}
		
		// if eol, print line counter and line break
		if ((j % row_size) == 0)		
		{	
			int current_row_nr = row_size*r;
			
			if (j!=0) 
			{
				_streamPut(file, '\n');
			}
			
			// leading space of new line
			for (int i=0;i<(seq_space_count);i++)
			{
				_streamPut(file,' ');
			}
			++r;
		}
		++j;

		_streamPut(file, *it);
	}
	_streamPut(file, '\n');

	_streamWrite(file, "//");
}

//____________________________________________________________________________


/**
.Function.write:
..summary: Creates and writes a sequence file in $EMBL$ format. 
..signature:write(file, data [, id] [, meta], Embl)
..param.file:The input file.
..param.data:The sequence data.
..param.id:The sequence id.
..param.meta:The meta data container.
..param.Embl: Embl Tag
...Type:Tag.File Format.EMBL
...remarks.text: Can be guessed by $guessFileFormat$.
..see:Function.readID
..see:Function.readMeta
..see:Function.write
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/


template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData & data,
	  Embl)
{
SEQAN_CHECKPOINT
	MetaMultMap dummy_map;
	_embl_write_impl(file, data, "", dummy_map, Embl());
}

template <typename TFile, typename TString, typename TData>
void
write(TFile & file,
	  TData & data,
	  TString & id,
	  Embl)
{
SEQAN_CHECKPOINT
	MetadataMap dummy_map;
	_embl_write_impl(file, data, id, dummy_map, Embl());
}

template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData & data,
	  TString & id,
	  TMeta & meta,
	  Embl)
{
SEQAN_CHECKPOINT
	_embl_write_impl(file, data, id, meta, Embl());
}

//VisualC++ const array bug workaround
template <typename TFile, typename TString, typename TDataValue>
void
write(TFile & file,
	  TDataValue * data,
	  Embl)
{
SEQAN_CHECKPOINT
	MetadataMap dummy_map;
	_embl_write_impl(file, data, "", dummy_map, Embl());

}


//VisualC++ const array bug workaround
template <typename TFile, typename TString, typename TDataValue>
void
write(TFile & file,
	  TDataValue * data,
	  TString & id,
	  Embl)
{
SEQAN_CHECKPOINT
	MetadataMap dummy_map;
	_embl_write_impl(file, data, id, dummy_map, Embl());

}

//VisualC++ const array bug workaround
template <typename TFile, typename TString, typename TDataValue, typename TMeta>
void
write(TFile & file,
	  TDataValue * data,
	  TString & id,
	  TMeta & meta,
	  Embl)
{
SEQAN_CHECKPOINT
	_embl_write_impl(file, data, id, meta, Embl());

}



//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
