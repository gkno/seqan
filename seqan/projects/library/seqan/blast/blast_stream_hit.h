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
#ifndef SEQAN_HEADER_BLAST_STREAM_HIT_H
#define SEQAN_HEADER_BLAST_STREAM_HIT_H


namespace SEQAN_NAMESPACE_MAIN
{


////////////////////////////////////////////////////////////////////////////////////////////
//  Blast Hit storing only one hsp at a time
////////////////////////////////////////////////////////////////////////////////////////////

template<typename TBlastHsp, typename TFile>
class BlastHit<TBlastHsp, StreamReport<TFile> > 
{
	public:
		typedef typename Position<TFile>::Type TPosition_;

		String<char> name;
		unsigned int length; //length of whole sequence  
		TBlastHsp act_hsp;
		TPosition_ begin_pos, first_hsp_pos;
		
		BlastReport<TBlastHsp,StreamReport<TFile> >* data_host;

	
		BlastHit()
		{
		}

		~BlastHit()
		{
		}

};




//parse BlastHit
template<typename TFile, typename TChar, typename TBlastSpec>
inline typename Position<TFile>::Type
_parseBlastHit(TFile & file,
			TChar & c, 
			BlastHit<TBlastSpec,StreamReport<TFile> > & hit)
{
	typedef typename Position<TFile>::Type TPosition;
	typedef BlastHit<TBlastSpec,StreamReport<TFile> > TBlastHit;
	typedef typename Hsp<TBlastHit>::Type TBlastHsp;

	String<char> pword;
	int pint;
	TPosition start_pos,act_pos;
	act_pos = _streamTellG(file);

	if(_parseUntilBeginLine(file,c,'>'))
	{
		hit.begin_pos = _streamTellG(file);
		c = _streamGet(file);
		pword = _parseReadWord(file, c);
		while (!_streamEOF(file) && c != '\n' && c != '\r')
			pword += _parseReadWord(file, c);
		if(pword[length(pword)-1] == ' ')
			resize(pword,length(pword)-1);
		hit.name = pword;
		_parseSkipWhitespace(file,c);
		String<char> search = "Length";
		if(_parseUntilBeginLine(file,c,search,6))
		{
			_parseSkipWhitespace(file,c);
			if(c == '=')
				c = _streamGet(file);
			_parseSkipWhitespace(file,c);
			pint = _parseReadNumber(file, c);
			hit.length = pint;
		}
//		TPosition temp = _streamTellG(file);
		//foreach Hsp
		//if(_parseUntilBeginLine(file,c,'S') && _parseReadWord(file,c)=="Score")
		search = "Score";
		if(_parseUntilBeginLine(file,c,search,5))
		{
			hit.first_hsp_pos = _streamTellG(file);
			//if(_parseUntilBeginLine(file,c,'>'))
			//	return _streamTellG(file);
		}

		
	}//end hit
	_streamSeekG(file,act_pos);
	c = '>';
	return act_pos;
}







}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
