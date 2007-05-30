#ifndef SEQAN_HEADER_FILE_GENBANK_H
#define SEQAN_HEADER_FILE_GENBANK_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - Genbank
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.value.Genbank:
	Genbank format for sequences from the Genbank database.
*/
struct TagGenbank_;
typedef Tag<TagGenbank_> const Genbank;



//////////////////////////////////////////////////////////////////////////////
// Internal functions ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Internal._genbank_scan_line:
..summary: Reads one line of a $file$ and increases $count$ (Internal use only).
...remarks: Do not moves back to the beginning of the line.
..signature:_genbank_scan_line(file, count)
..param.file: The input file.
..param.count: The counter which is increased by the number of signs in the given line. 
..remarks.text: Same as $_embl_scan_line$.
..see:Internal._embl_scan_line.
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/
template <typename TFile, typename TSize>
void
_genbank_scan_line(	TFile & file,
			TSize & count)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	typename Value<TFile>::Type c = _streamGet(file);
	
	while ((c!='\n')&&(!_streamEOF(file)))
	{
		if (c!='\r') ++count;
		c = _streamGet(file);
	}
}

/**
.Internal._genbank_scan_tag:
..summary: Reads characters from a stream until the next '\n', 't' or " " and increases 
$count$ (Internal use only).
...remarks: Do not moves back to the beginning of the line.
..signature:_genbank_scan_tag(file, count)
..param.file: The input file.
..param.count: The counter which is increased by the number of characters read. 
..see:Tag.File Format.value.Genbank
*/

template <typename TFile, typename TSize>
void
_genbank_scan_tag(	TFile & file,
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
} // _genbank_scan_tag


/**
.Internal._genbank_det_line:
..summary: Determines thy type of one line in a $Genbank$ record.
...remarks: Moves back to the beginning of the line after exaamining the line.
..signature:_genbank_det_line(file)
..param.file: The input file.
..see:Tag.File Format.value.Genbank
*/

template <typename TFile>
int // 0 = new tag, 1 = old tag, 2 = sequence begins , 3 = start of FEATURES block
_genbank_det_line(TFile & file)
{	// must beginning of new line
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	String<char,Array<8> > current_tag = "        ";
	typename Position<TFile>::Type line_start_pos = _streamTellG(file);
	typename Value<TFile>::Type c = _streamGet(file);
	
	if (c == 'F')
	{
		current_tag[0] = c;
		current_tag[1] = _streamGet(file);
		current_tag[2] = _streamGet(file);
		current_tag[3] = _streamGet(file);
		current_tag[4] = _streamGet(file);
		current_tag[5] = _streamGet(file);
		current_tag[6] = _streamGet(file);
		current_tag[7] = _streamGet(file);
		if (current_tag == "FEATURES") 
		{
			_streamSeekG(file, line_start_pos);
			return 3;
		}
		else
		{
			_streamSeekG(file, line_start_pos);
			return 0;
		}
	}
	
	// if sequence begins return 2
	if (c == 'O')
	{
		current_tag[0] = c;
		current_tag[1] = _streamGet(file);
		current_tag[2] = _streamGet(file);
		current_tag[3] = _streamGet(file);
		current_tag[4] = _streamGet(file);
		current_tag[5] = _streamGet(file);
		current_tag[6] = ' ';
		current_tag[7] = ' ';
		if (current_tag == "ORIGIN  ") 
		{
			_streamSeekG(file, line_start_pos);
			return 2;
		}
		else
		{
			_streamSeekG(file, line_start_pos);
			return 0;
		}
	}
	
	// if space or tab return 1, else 0
	if ((c == ' ')||(c == '\t'))
	{
		_streamSeekG(file, line_start_pos);
		return 1;  // add to old feature
	}
	else 
	{
		_streamSeekG(file, line_start_pos);
		return 0; // start new feature, move back to beginning of line
	}
} // _genbank_det_line


//////////////////////////////////////////////////////////////////////////////
// read
//////////////////////////////////////////////////////////////////////////////

/**
.Function.read:
..summary: Reads the sequence data of a $_Genbank$ file and stores it in $data$. 
..signature:read(file, data [, limit], Genbank)
..param.file:The input file.
..param.data:The data object in which the sequence data is parsed.
..param.meta:The meta data container.
..param.Embl: Embl Tag
...Type:Tag.File Format.Genbank
...remarks.text: Can be guessed by $guessFileFormat$.
..see:Function.readID
..see:Function.readMeta
..see:Function.write
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/

template <typename TFile, typename TData, typename TSize>
void
read(TFile & file,
	 TData & data,
	 TSize limit,
	 Genbank)
{
SEQAN_CHECKPOINT

	SEQAN_ASSERT(!_streamEOF(file))

	typename Position<TFile>::Type start_pos = _streamTellG(file);	// default data start position
	typename Position<TFile>::Type end_pos;							// data end position
	typename Size<TData>::Type count = 0;							// default TData size
	typename Value<TFile>::Type c;									// holds one character from stream
	String<char,Array<6> > origin_tag = "      ";					// holds six characters
	bool origin_found = false;

	// find ORIGIN tag
	while ((!origin_found)&&(!_streamEOF(file)))
	{
		// reads the first six characters from line
		for (int i=0; i<length(origin_tag);++i)
		{
			origin_tag[i] = _streamGet(file);
		}

		// if ORIGIN tag found, set flag = true
		if (origin_tag=="ORIGIN")
		{
			origin_found = true;
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
	_genbank_scan_line(file, count); 
	start_pos = _streamTellG(file);
	
	count = 0;
	
	// determine size of sequence
	bool exit_flag = false;
	while ((!exit_flag)&&(!_streamEOF(file)))
	{
		c = _streamGet(file);
		int c_int = c;
		if (((c_int>64)&&(c_int<91))||((c_int>96)&&(c_int<123)))
		{
			++count;
		}
		
		// check for end of file
		if (c=='/') 
		{	
			c = _streamGet(file);
			if (c=='/') 
			exit_flag = true;
			_streamSeek2G(file, -2);
	
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
	 Genbank tag)
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
..summary: Reads the id of a $Genbank$ sequence record and stores it in $id$. 
..signature:readID(file, id, Genbank)
..param.file:The input file.
..param.id:The object in which the sequence $id$ is stored.
..param.Genbank: Embl Tag
...Type:Tag.File Format.Genbank
...remarks.text: Can be guessed by $guessFileFormat$. 
..remarks.text: The ID in an EMBL sequence record is the first entry in the line beginning with 'AC', 
not with 'ID'. Moves back to the beginning of file.
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
	   Genbank)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	typename Position<TFile>::Type begin_pos = _streamTellG(file);	// contains start position of gb entry
	typename Position<TFile>::Type id_pos = _streamTellG(file);	    // contains start position of id
	typename Size<TString>::Type count = 0;							// length of id
	typename Value<TFile>::Type c;									// holds one character from stream
	String<char,Array<9> > locus_tag = "         ";					// holds the ACCESSION tag
	
	bool exit_flag = false;
	while (!exit_flag)
	{
		//  read first nine chars from file
		for (int i =0 ;i<length(locus_tag);++i)
		{	
			locus_tag[i] = _streamGet(file);	
		}

	//  check if its ==ACESSION
	if (locus_tag=="ACCESSION")
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
				
				if ((c!=' ')&&(c!='\n'))
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

	} // if not LOCUS empty id and exit
	else 
	{
		if (!_streamEOF(file))
		{
			_genbank_scan_line(file, count);
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
..summary: Reads the meta data of a $Genbank$ file and stores it in $meta$. 
..signature:readMeta(file, data [, meta], Genbank)
..param.file:The input file.
..param.meta:The meta data container.
..param.Embl: Genbank Tag
...Type:Tag.File Format.Genbank
...remarks.text: Can be guessed by $guessFileFormat$. 
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
		 Genbank)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	typename Position<TFile>::Type act_pos = _streamTellG(file); // store starting position of stream
	typename Value<TFile>::Type c = '\n';						// forces to check the first line 
	typename Size<TMeta>::Type count = 0;
	typename Position<String<char, Alloc<char> > >::Type pos = 0;
	
	clear(meta); // empties old meta object 

	String<char> tag;
	String<char> content;

	unsigned int line_type = 0;
	bool features_found = false;
	short features_space_count = 0;

	// loop while not eof
	while (!_streamEOF(file))
	{
		// if eol reached, check what kind  is the next line of
		if (c == '\n') 
		{
			// skip leading spaces of subfeatures if in FEATURES block
			line_type = _genbank_det_line(file);
			
			// determines line type if already in tabbed features block
			if ((features_found==true)&&(line_type!=2)) 
			{
				for (short i=0;i<features_space_count;++i)
				{
					c = _streamGet(file);
				}
				line_type = _genbank_det_line(file);
			}
			if (line_type!=1) // add old tag to multimap 
			{
				meta.insert(MetadataEntry(tag, content));
			}

			// check for FEATURES block
			if (line_type == 3) 
			{
				features_found = true;

				// store 'FEATURE' feature
				_genbank_scan_tag(file, count); // overread 'FEATURES'
				
				c = _streamGet(file);
				while((c==' ')||(c=='\t'))
				{
					c = _streamGet(file);
				}
				_streamSeek2G(file, -1);

				count = 0;
				act_pos = _streamTellG(file);
				_genbank_scan_line(file, count);
				clear(content);
				resize(content, count+1); 
				
				_streamSeekG(file, act_pos);
				for (pos=0;pos<count;pos++)
				{
					content[pos] = _streamGet(file); 
				}
				content[length(content)-1] = '\n';
				
				meta.insert(MetadataEntry("FEATURES",content));
				
				// examine count of leading spaces of the sub-features of FEATURES
				_genbank_scan_line(file,count);  // overread 'FEATURE' line 
				act_pos = _streamTellG(file);
				
				count = 0;
				c = _streamGet(file);
				while (c==' ')
				{
					c = _streamGet(file);
					++count;
				}
				
				features_space_count = count;
				_streamSeekG(file, act_pos); // move back
				
				for (short i=0;i<features_space_count;++i)
				{
					c = _streamGet(file);
				}

			}

			// check if meta dat ends and sequence begins
			if (line_type == 2) 
			{
				return; // exit if sequence start found
			}

			// if new tag
			if ((line_type == 0)||(line_type == 3))
			{
				// read the new tag
				count = 0;
				act_pos = _streamTellG(file);
				_genbank_scan_tag(file, count);
				clear(tag);

				int start_pos = 0;
				int end_pos = count;

				// mark tag, if in FEATURES block
				if (features_found) 
				{
					resize (tag,count+1);
					start_pos = 1;
					end_pos++;
					tag[0] = '|'; // internal marker for feature, that lies within the feature block	
				}
				else
				{
					resize(tag,count);
				}

				_streamSeekG(file,act_pos);
				for (pos=start_pos;pos<end_pos;pos++)
				{
					tag[pos] = _streamGet(file);
				}

				// read over space and tabs
				c = _streamGet(file);
				while((c==' ')||(c=='\t'))
				{
					c = _streamGet(file);
				}
				_streamSeek2G(file, -1);
				
				// read the new feature
				count = 0;
				act_pos = _streamTellG(file);
				_genbank_scan_line(file, count);
				clear(content);
				resize(content, count+1); 
				
				_streamSeekG(file, act_pos);
				for (pos=0;pos<count;pos++)
				{
					content[pos] = _streamGet(file); 
				}
				content[length(content)-1] = '\n';
				
			}
			if (line_type == 1)
			{
				// read over leading space and tabs
				c = _streamGet(file);
				while((c==' ')||(c=='\t'))
				{
					c = _streamGet(file);
				}
				_streamSeek2G(file, -1);
				
				count = 0;
				act_pos = _streamTellG(file);
				_genbank_scan_line(file, count);

				int old_length = length(content) ; // must be int or specific type?? 
				resize(content, (length(content)+count)+1); 

				_streamSeekG(file, act_pos);
				for (pos=old_length;pos<length(content)-1;pos++)
				{
					content[pos] = _streamGet(file); 
				}
				content[length(content)-1]= '\n';
			}
		}// if 
		c = _streamGet(file);
	}
} // readMeta


/**
.Function.goNext:
..summary: Moves to the next sequence record in multi $Genbank$ file. 
..signature:goNext(file, Genbank)
..param.file:The input file.
..param.Genbank: Genbank Tag
...Type:Tag.File Format.Genbank
...remarks.text: Can be guessed by $guessFileFormat$. 
..remarks.text: $EOF$ if no record is found.
..see:Function.read
..see:Function.readID
..see:Function.readMeta
..see:Function.write
..see:Tag.File Format.value.Genbank
..see:Tag.File Format.value.EMBL
*/

template <typename TFile>
void
goNext(TFile & file,
	   Genbank)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))
	typename Position<TFile>::Type end_pos = _streamTellG(file);
	typename Value<TFile>::Type c;	

	String<char,Array<5> > locus_tag = "     ";	

	while (true)
	{
		c = _streamGet(file);
		
		// either return if eof
		if (_streamEOF(file))
		{
			end_pos = _streamTellG(file); // set stream to eof if no next dataset found
			_streamSeekG(file, end_pos);
			return;
		} // file ends
		
		// or if next data set found
		if (c=='L')
		{
			locus_tag[0] = c;
			
			for (int pos=1;pos<5;++pos)
			{
				locus_tag[pos] = _streamGet(file);
			}

			if (locus_tag == "LOCUS")
			{
				_streamSeek2G(file, -5);  // set stream to next dataset if found
				end_pos = _streamTellG(file);
				_streamSeekG(file, end_pos);
				return;
			}
		}	
	}
	
} // goNext


//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////

/**
.Internal._genbank_write_impl:
..summary: Implementation of the $write$ function for $Genbank$ sequence records. 
..signature:_embl_write_impl(file, data [, id] [, meta], Genbank)
..param.file:The input file.
..param.data:The sequence data.
..param.id:The sequence id.
..param.meta:The meta data container.
..param.Genbank: Genbank Tag.
..see:Function.write
*/

template <typename TFile, typename TString, typename TData, typename TMeta>
void
_genbank_write_impl(TFile & file,
			TData & data,
			TString & id,
			TMeta & meta,
			Genbank)
{
SEQAN_CHECKPOINT
	SEQAN_ASSERT(!_streamEOF(file))

	// constants due to genbank standard
	int const row_size			= 60;  // length of one row
	int const block_size		= 10;  // size of one block 
	int const tag_gap			= 3;   // # of spaces between the ID-tag and the id, has to be >= 1
	int const seq_space_count	= 5; 
	int const feat_space_count	= 21;
	int const feat_lead_count	= 5;

	//write LOCUS
	_streamWrite(file, "LOCUS");
	for (int i=0;i<tag_gap;i++)
	{
		_streamPut(file,' ');
	}

	MetadataMap::iterator act_entry = meta.find("LOCUS");
	String<char> content;
//	resize(content, length(act_entry->second));
//	clear(content);
//	content = act_entry->second;
	_streamWrite(file,content);

	//write DEFINITION 
	_streamWrite(file, "DEFINITION");
	for (int i=0;i<tag_gap;i++)
	{
		_streamPut(file,' ');
	}
	
	act_entry = meta.find("DEFINITION");
	_streamWrite(file,act_entry->second);

	
	// write ID (ACCESSION)
	_streamWrite(file, "ACCESSION");
	for (int i=0;i<tag_gap;i++)
	{
		_streamPut(file,' ');
	}
	_streamWrite(file, id);
	_streamPut(file, '\n');	

	//writing meta data 

	String<char> feat_tag;				// holds the feature tag = key
	String<char> meta_content;			// holds the content of the meta data
	MetadataMap::iterator meta_iter;	// generic iterator used for meta-data
	typename Value<String<char> >::Type c = ' ';

	for (meta_iter = meta.begin(); meta_iter != meta.end(); meta_iter++ )
	{
		// write meta key
		bool feat = false;	// set if FEATURES meta data found
		clear(feat_tag);
		resize(feat_tag, length(meta_iter->first));
		feat_tag	= meta_iter->first;
		
		if (length(feat_tag)>0) // check for empty tag
		{		
			if (feat_tag[0]=='|') 
			{ 
				feat = true; 
			}

			if ((meta_iter->first!="LOCUS")&&(meta_iter->first!="DEFINITION")&&(meta_iter->first!="ACCESSION")&&(meta_iter->first!="FEATURES")&&(feat==false))
			{
				_streamWrite(file, meta_iter->first);
			
				int leading_key_spaces;
				if (length(meta_iter->first)>9)
				{
					leading_key_spaces = tag_gap;
				}
				else 
				{
					leading_key_spaces = block_size - length(meta_iter->first);
				}
		
				for (int i=0;i<leading_key_spaces;i++)
				{
					_streamPut(file,' ');
				}
		
				// writing meta content

				clear(meta_content);
				resize(meta_content, length(meta_iter->second));
				meta_content = meta_iter->second;
			
				for (int i=0;i<length(meta_content);i++)
				{
					c = meta_content[i];
			
					_streamPut(file, c); 
			
					if (c=='\n')
					{
						if (i!=(length(meta_content)-1))
						{
							for (int j=0;j<block_size;j++)
							{
								_streamPut(file, ' ');
							}
						}
					}
	
				}
			}
		} // if length(feat_tag)>0
	}

	// writing FEATURES meta data
	_streamWrite(file, "FEATURES");
	for (int i=0;i<(feat_space_count-8);i++)  // 8==lenght("FEATURES")
	{
		_streamPut(file,' ');
	}
	act_entry = meta.find("FEATURES");
	_streamWrite(file,act_entry->second);

	for (meta_iter = meta.begin(); meta_iter != meta.end(); meta_iter++ )
	{
		bool feat = false;	// set if FEATURES meta data found
		/*clear(feat_tag);
		resize(feat_tag, length(meta_iter->first));
		feat_tag	= meta_iter->first;*/

		if (length(meta_iter->first)>0)  // check for empty tag
		{
			if (meta_iter->first[0]=='|')
			{
				feat = true;
			}

			if (feat==true)
			{
				// write meta key
				for (int i=0;i<feat_lead_count;i++) 
				{
					_streamPut(file, ' ');
				}
				for (int i=1;i<length(meta_iter->first);i++)
				{
					_streamPut(file,meta_iter->first[i]);
				}
				
				for (int j=0;j<(feat_space_count-(length(meta_iter->first)-1)-feat_lead_count);j++)
				{	
					_streamPut(file, ' ');
				}
				
				// write meta data content
				/*clear(meta_content);
				resize(meta_content, length(meta_iter->second));
				meta_content = meta_iter->second;*/

				for (int i=0;i<length(meta_iter->second);i++)
				{
					c = meta_iter->second[i];
					_streamPut(file, c); 
			
					if (c=='\n')
					{
						if (i!=(length(meta_iter->second)-1))
						{
							for (int j=0;j<feat_space_count;j++)
							{
								_streamPut(file, ' ');
							}
						}
					}
				}

			}
		}

	}

	// write sequence data
	_streamWrite(file, "ORIGIN");
	
	typename Iterator<TData>::Type it		= begin(data);
	typename Iterator<TData>::Type it_end	= end(data);

	int j = 0; // coloumn counter
	int r = 1; // row counter

	for (; it < it_end; ++it)
	{	
		// space, if block completed
		
		if ((j % block_size) == 0)
		{
			_streamPut(file, ' '); 
		}
		
		// if eol, print line counter and line break
		if ((j % row_size) == 0)		
		{
			int current_row_nr = row_size*r;
			_streamPut(file, '\n');
			
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


template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData & data,
	  Genbank)
{
SEQAN_CHECKPOINT
	MetadataMap dummy_map;
	_genbank_write_impl(file, data, "",dummy_map Genbank());
}

template <typename TFile, typename TString, typename TData>
void
write(TFile & file,
	  TData & data,
	  TString & id,
	  Genbank)
{
SEQAN_CHECKPOINT
	MetadataMap dummy_map;
	_genbank_write_impl(file, data, id, dummy_map, Genbank());
}



template <typename TFile, typename TString, typename TData, typename TMeta>
void
write(TFile & file,
	  TData & data,
	  TString & id,
	  TMeta & meta,
	  Genbank)
{
SEQAN_CHECKPOINT
	_genbank_write_impl(file, data, id, meta, Genbank());
}

//VisualC++ const array bug workaround
template <typename TFile, typename TString, typename TDataValue, typename TMeta>
void
write(TFile & file,
	  TDataValue * data,
	  TString & id,
	  TMeta & meta,
	  Genbank)
{
SEQAN_CHECKPOINT
	_genbank_write_impl(file, data, id, meta, Genbank());
}

//VisualC++ const array bug workaround
template <typename TFile, typename TString, typename TDataValue>
void
write(TFile & file,
	  TDataValue * data,
	  TString & id,
	  Genbank)
{
SEQAN_CHECKPOINT
	MetadataMap dummy_map;
	_genbank_write_impl(file, data, id, dummy_map, Genbank());
}

//VisualC++ const array bug workaround
template <typename TFile, typename TString, typename TDataValue>
void
write(TFile & file,
	  TDataValue * data,
	  Genbank)
{
SEQAN_CHECKPOINT
	MetadataMap dummy_map;
	_genbank_write_impl(file, data, "" , dummy_map, Genbank());
}
//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
